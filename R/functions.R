#' Returns a summary list for a quantile regression fit.
#'
#' This function is slightly adapted based on the summary.rq function in the quantreg package
mysummary_rq= function (object, se = NULL, covariance = FALSE, hs = TRUE,
                        U = NULL, gamma = 0.7, ...) {
  if (object$method == "lasso")
    stop("no inference for lasso'd rq fitting: try rqss (if brave, or credulous)")
  if (object$method == "conquer")
    se = "conquer"
  mt <- terms(object)
  m <- model.frame(object)
  y <- model.response(m)
  dots <- list(...)
  method <- object$method
  if (object$method == "sfn") {
    x <- object$model$x
    vnames <- names(object$coef)
    ctrl <- object$control
  }
  else {
    x <- model.matrix(mt, m, contrasts = object$contrasts)
    vnames <- dimnames(x)[[2]]
  }
  wt <- as.vector(model.weights(object$model))
  tau <- object$tau
  eps <- .Machine$double.eps^(1/2)
  coef <- coefficients(object)
  if (is.matrix(coef))
    coef <- coef[, 1]
  resid <- object$residuals
  n <- length(y)
  p <- length(coef)
  rdf <- n - p
  if (!is.null(wt)) {
    resid <- resid * wt
    x <- x * wt
    y <- y * wt
  }
  if (is.null(se)) {
    if (n < 1001 & covariance == FALSE)
      se <- "rank"
    else se <- "nid"
  }
  if (se == "rank") {
    f <- rq.fit.br(x, y, tau = tau, ci = TRUE, ...)
  }
  if (se == "iid") {
    xxinv <- diag(p)
    xxinv <- backsolve(qr(x)$qr[1:p, 1:p, drop = FALSE],
                       xxinv)
    xxinv <- xxinv %*% t(xxinv)
    pz <- sum(abs(resid) < eps)
    h <- max(p + 1, ceiling(n * bandwidth.rq(tau, n, hs = hs)))
    ir <- (pz + 1):(h + pz + 1)
    ord.resid <- sort(resid[order(abs(resid))][ir])
    xt <- ir/(n - p)
    sparsity <- rq(ord.resid ~ xt)$coef[2]
    cov <- sparsity^2 * xxinv * tau * (1 - tau)
    scale <- 1/sparsity
    serr <- sqrt(diag(cov))
  }
  else if (se == "nid") {
    h <- bandwidth.rq(tau, n, hs = hs)
    while ((tau - h < 0) || (tau + h > 1)) h <- h/2
    bhi <- rq.fit(x, y, tau = tau + h, method = method)$coef
    blo <- rq.fit(x, y, tau = tau - h, method = method)$coef
    dyhat <- x %*% (bhi - blo)
    if (any(dyhat <= 0))
      warning(paste(sum(dyhat <= 0), "non-positive fis"))
    f <- pmax(0, (2 * h)/(dyhat - eps))
    fxxinv <- diag(p)
    if (method == "sfn") {
      D <- t(x) %*% (f * x)
      D <- chol(0.5 * (D + t(D)), nsubmax = ctrl$nsubmax,
                nnzlmax = ctrl$nnzlmax, tmpmax = ctrl$tmpmax)
      fxxinv <- backsolve(D, fxxinv)
    }
    else {
      fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p,
                                             drop = FALSE], fxxinv)
      fxxinv <- fxxinv %*% t(fxxinv)
    }
    xx <- t(x) %*% x
    cov <- tau * (1 - tau) * fxxinv %*% xx %*% fxxinv
    scale <- mean(f)
    serr <- sqrt(diag(cov))
  }
  else if (se == "ker") {
    h <- bandwidth.rq(tau, n, hs = hs)
    while ((tau - h < 0) || (tau + h > 1)) h <- h/2
    uhat <- c(y - x %*% coef)
    h <- (qnorm(tau + h) - qnorm(tau - h)) * min(sqrt(var(uhat)),
                                                 (quantile(uhat, 0.75) - quantile(uhat, 0.25))/1.34)
    f <- dnorm(uhat/h)/h
    fxxinv <- diag(p)
    fxxinv <- backsolve(qr(sqrt(f) * x)$qr[1:p, 1:p, drop = FALSE],
                        fxxinv)
    fxxinv <- fxxinv %*% t(fxxinv)
    cov <- tau * (1 - tau) * fxxinv %*% crossprod(x) %*%
      fxxinv
    scale <- mean(f)
    serr <- sqrt(diag(cov))
  }
  else if (se == "boot") {
    if ("cluster" %in% names(dots)) {
      bargs <- modifyList(list(x = x, y = y, tau = tau),
                          dots)
      if (length(object$na.action)) {
        cluster <- dots$cluster[-object$na.action]
        bargs <- modifyList(bargs, list(cluster = cluster))
      }
      if (class(bargs$x)[1] == "matrix.csr")
        bargs <- modifyList(bargs, list(control = ctrl))
      B <- do.call(boot.rq, bargs)
    }
    else B <- boot.rq(x, y, tau, coef = coef, ...)
    cov <- cov(B$B)
    serr <- sqrt(diag(cov))
  }
  else if (se == "BLB") {
    n <- length(y)
    b <- ceiling(n^gamma)
    S <- n%/%b
    V <- matrix(sample(1:n, b * S), b, S)
    Z <- matrix(0, NCOL(x), S)
    for (i in 1:S) {
      v <- V[, i]
      B <- boot.rq(x[v, ], y[v], tau, bsmethod = "BLB",
                   blbn = n, ...)
      Z[, i] <- sqrt(diag(cov(B$B)))
    }
    cov <- cov(B$B)
    serr <- apply(Z, 1, mean)
  }
  else if (se == "extreme") {
    tau0 <- tau
    if (tau > 0.5) {
      y <- -y
      tau <- 1 - tau
    }
    if (length(dots$mofn))
      mofn = dots$mofn
    else mofn = floor(n/5)
    if (length(dots$mofn))
      kex = dots$kex
    else kex = 20
    if (length(dots$alpha))
      alpha = dots$alpha
    else alpha = 0.1
    if (length(dots$R))
      R = dots$R
    else R = 200
    m <- (tau * n + kex)/(tau * n)
    taub <- min(tau * n/mofn, tau + (0.5 - tau)/3)
    xbar <- apply(x, 2, mean)
    b0 <- rq.fit(x, y, tau, method = method)$coef
    bm <- rq.fit(x, y, tau = m * tau, method = method)$coef
    An <- (m - 1) * tau * sqrt(n/(tau * (1 - tau)))/c(crossprod(xbar,
                                                                bm - b0))
    bt <- rq.fit(x, y, tau = taub, method = method)$coef
    s <- matrix(sample(1:n, mofn * R, replace = T), mofn,
                R)
    mbe <- (taub * mofn + kex)/(taub * mofn)
    bmbeb <- rq.fit(x, y, tau = mbe * taub, method = method)$coef
    B0 <- boot.rq.pxy(x, y, s, taub, bt, method = method)
    Bm <- boot.rq.pxy(x, y, s, tau = mbe * taub, bmbeb,
                      method = method)
    B <- (mbe - 1) * taub * sqrt(mofn/(taub * (1 - taub))) *
      (B0 - b0)/c((Bm - B0) %*% xbar)
    if (tau0 <= 0.5) {
      bbc <- b0 - apply(B, 2, quantile, 0.5, na.rm = TRUE)/An
      ciL <- b0 - apply(B, 2, quantile, 1 - alpha/2, na.rm = TRUE)/An
      ciU <- b0 - apply(B, 2, quantile, alpha/2, na.rm = TRUE)/An
    }
    else {
      bbc <- -(b0 - apply(B, 2, quantile, 0.5, na.rm = TRUE)/An)
      ciL <- -(b0 - apply(B, 2, quantile, alpha/2, na.rm = TRUE)/An)
      ciU <- -(b0 - apply(B, 2, quantile, 1 - alpha/2,
                          na.rm = TRUE)/An)
    }
    B <- R - sum(is.na(B[, 1]))
    coef <- cbind(b0, bbc, ciL, ciU)
    if (tau0 > 0.5) {
      coef <- -coef
      tau <- tau0
    }
    dimnames(coef) = list(dimnames(x)[[2]], c("coef", "BCcoef",
                                              "ciL", "ciU"))
  }
  else if (se == "conquer") {
    if (length(dots$R))
      R = dots$R
    else R = 200
    Z <- conquer(x[, -1], y, tau, ci = TRUE, B = R)
    coef <- cbind(Z$coef, Z$perCI)
    cnames <- c("coefficients", "lower bd", "upper bd")
    dimnames(coef) <- list(vnames, cnames)
    resid <- y - x %*% Z$coef
  }
  if (se == "rank") {
    coef <- f$coef
  }
  else if (!(se %in% c("conquer", "extreme"))) {
    coef <- array(coef, c(p, 4))
    dimnames(coef) <- list(vnames, c("Value", "Std. Error",
                                     "t value", "Pr(>|t|)"))
    coef[, 2] <- serr
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- if (rdf > 0)
      2 * (1 - pt(abs(coef[, 3]), rdf))
    else NA
  }
  object <- object[c("call", "terms")]
  if (covariance == TRUE) {
    if (se != "rank")
      object$cov <- cov
    if (se == "iid")
      object$scale <- scale
    if (se %in% c("nid", "ker")) {
      object$Hinv <- fxxinv
      object$J <- crossprod(x)
      object$scale <- scale
    }
    else if (se == "boot") {
      object$B <- B$B
      object$U <- B$U
    }
  }
  object$coefficients <- coef
  object$residuals <- resid
  object$rdf <- rdf
  object$tau <- tau
  object$cov <- cov
  class(object) <- "summary.rq"
  object
}



#' Generate a simulated dataset.
#' @param n sample size
#' @return a simulated data set with K=3 time periods
data.gen=function(n=5000) {
  # data generation
  X1=cbind(rnorm(n,mean=0,sd=1),rnorm(n,mean=0,sd=1))
  colnames(X1)=c("X11","X12")
  A0=rep(0,n)
  pi1=1/(1+exp(-( -1 -A0+2*(X1[,1]>0)+1*(X1[,2]>0))))
  A1= as.numeric( runif(n) < pi1 )
  U2= A1
  X2=cbind(rnorm(n,mean=0,sd=1)+U2,rnorm(n,mean=0,sd=1)+2*U2)
  colnames(X2)=c("X21","X22")
  pi2=1/(1+exp(-( -1.5-1.5*A1+2*(X2[,1]>0)+1*(X2[,2]>0))))
  A2= as.numeric( runif(n) < pi2 )
  U3=A2
  X3=cbind(rnorm(n,mean=0,sd=1)+U3,rnorm(n,mean=0,sd=1)+2*U3)
  colnames(X3)=c("X31","X32")
  pi3=1/(1+exp(-( -2-1.5*A2+2*(X3[,1]>0)+1*(X3[,2]>0))))
  A3= as.numeric( runif(n) < pi3 )
  sd.e=2+2*A3
  epsilonY=rnorm(n,sd=sd.e)
  deltaX1=2*X1[,1]+2*X1[,2]
  deltaX2=2*X2[,1]+2*X2[,2]
  deltaX3=2*X3[,1]+2*X3[,2]
  Y=10-10*(A1+A2+A3)+deltaX1+deltaX2+deltaX3+epsilonY
  mydata=as.data.frame(cbind(X1,X2,X3,A1,A2,A3,Y))
  mydata
}


#' Naive estimator
naive.f=function(data,q,MSQM.formula) {
  mod=rq(formula=MSQM.formula,tau=q,data=data,method="fn")
  theta.est = mod$coefficients
  myrq=mysummary_rq(mod,se="ker")
  theta.var = myrq$coefficients[,2]
  res=list(Point=theta.est,SE=theta.var,
           SE.TE = sqrt(abs((t(c(0,rep(1,length(theta.est)-1))) %*% myrq$cov %*% c(0,rep(1,length(theta.est)-1)))[1,1])),
           tau=sd(mod$residuals))
  res
}



#' IPW estimator
IPW.f=function(data,q,PS.formula,S.formula,MSQM.formula,Outcome.name,start.value,tau_n) {
  n.steps=length(PS.formula)
  n=dim(data)[1]
  alpha = alpha_sep = c()

  PS.model=S.model=list()
  for (j in (1:n.steps)) {
    ## Logistic model
    PS.model[[j]] = glm(PS.formula[[j]],family=binomial(link="logit"),data=data)
    S.model[[j]] = glm(S.formula[[j]],family=binomial(link="logit"),data=data)
    ## parameter
    alpha = c(alpha,PS.model[[j]]$coefficients,S.model[[j]]$coefficients)
    alpha_sep = c(alpha_sep,length(PS.model[[j]]$coefficients),length(S.model[[j]]$coefficients))
    ## sandwich variance terms for PS models
    if (j==1) {
      Ua = cbind(estfun(PS.model[[j]]),estfun(S.model[[j]]));
      Ia=magic::adiag(bread(PS.model[[j]]),bread(S.model[[j]]))
    } else {
      Ua = cbind(Ua,estfun(PS.model[[j]]),estfun(S.model[[j]]));
      Ia=magic::adiag(Ia,bread(PS.model[[j]]),bread(S.model[[j]]))
    }
  }
  alpha_sep=c(0,cumsum(alpha_sep))
  Xmat = model.matrix(MSQM.formula,data=data)
  n.theta=dim(Xmat)[2]
  # type 1 estimation, type 2 meat, type 3 bread
  est.eq=function(theta,alpha=1,type=1) {
    if (type==3) {
      alpha = theta[-(1:dim(Xmat)[2])];
      theta = theta[1:dim(Xmat)[2]]
    }
    # calculate the weight
    IPW.mat = SW.mat = matrix(NA,ncol=n.steps,nrow=n)
    for (j in (1:n.steps)) {
      alpha0 = alpha[(alpha_sep[2*j-1]+1):alpha_sep[2*j]]
      alphas0 = alpha[(alpha_sep[2*j]+1):alpha_sep[2*j+1]]
      # PS model
      p.logis=plogis(c(model.matrix(PS.model[[j]]) %*% alpha0))
      IPW.mat[,j]=(p.logis^(PS.model[[j]]$y))*((1-p.logis)^(1-PS.model[[j]]$y))
      # S model
      p.logis=plogis(c(model.matrix(S.model[[j]]) %*% alphas0))
      SW.mat[,j]=(p.logis^(S.model[[j]]$y))*((1-p.logis)^(1-S.model[[j]]$y))
    }
    sw = apply(SW.mat,1,prod)/apply(IPW.mat,1,prod)
    # calculate the estimating equation
    expit=function(x,m,gamma=tau_n) {plogis(m-x,scale=gamma)}
    y=as.numeric(Xmat %*% theta)
    Pval=expit(x=data[,Outcome.name],m=y)
    out=matrix(NA,ncol=n.theta,nrow=n)
    for (m in (1:n.theta)) {
      out[,m] = Xmat[,m]*sw*(Pval-q)
    }
    if (type==1) {return(apply(out,2,sum))}
    if (type==2) {return(out)}
    if (type==3) {return(apply(out,2,sum))}
  }
  ### point estimation
  theta.est=nleqslv::nleqslv(x=start.value,fn=est.eq,alpha=alpha,type=1)$x
  ### sandwich variance
  ## U_theta
  Ut = est.eq(theta=theta.est,alpha=alpha,type=2)
  ## I_theta and S_theta_alpha
  I_esteq = numDeriv::jacobian(est.eq, x=c(theta.est,alpha), method="simple",
                               method.args=list(eps=1e-4),
                               type=3)/n
  I_theta = I_esteq[,1:dim(Xmat)[2]]
  S_theta_alpha = I_esteq[,-c(1:dim(Xmat)[2])]
  ## calculate the sandwich variance
  # meat
  V_ipw = crossprod(Ut + t(S_theta_alpha %*% Ia %*% t(Ua)))/n
  theta.vcov = (solve(I_theta) %*% V_ipw %*% t(solve(I_theta)))/n

  # extract alpha and alphas
  alpha0=list(PS.model[[1]]$coefficients)
  alphas0=list(S.model[[1]]$coefficients)
  for (j in (2:n.steps)) {
    alpha0[[j]] = PS.model[[j]]$coefficients
    alphas0[[j]] = S.model[[j]]$coefficients
  }
  res=list(Point=theta.est,
           SE=sqrt(abs(diag(theta.vcov))),
           SE.TE = sqrt(abs((t(c(0,rep(1,n.steps))) %*% theta.vcov %*% c(0,rep(1,n.steps)))[1,1])),
           alpha=alpha0,
           alphas=alphas0,
           Ua=Ua,
           Ia=Ia,
           SE.mat = theta.vcov)
  res
}

#' ICR estimator (point estimation)
ICR.f.point=function(data,q,
                     Outcome.formula,
                     Var.formula,
                     MSQM.formula,
                     A.names,Outcome.name,start.value,icr_coef=0) {
  n.steps=length(A.names)
  n=dim(data)[1]
  coef.list=sigma2.list=Outcome.m.list=Var.m.list=list()
  Outcome.m.list[[1]] = lm(Outcome.formula[[n.steps]],data=data)
  coef.list[[1]]=coef(Outcome.m.list[[1]]);
  EY=predict(Outcome.m.list[[1]],newdata=data)
  data$Sigma2 = (Outcome.m.list[[1]]$residuals)^2
  Var.m.list[[1]] = glm(as.formula(paste("Sigma2 ~",as.character(Var.formula[[n.steps]])[2])),family = quasipoisson(link = "log"),data=data)
  sigma2.list[[1]]=Var.m.list[[1]]$coefficients
  for (j in (1:n.steps)) {
    # generate the extended dataset
    Alist=list()
    for (i in (1:j)) { Alist[[i]]=c(0,1) }
    Alist=expand.grid(Alist);
    Aname.now=A.names[(n.steps-j+1):n.steps]
    colnames(Alist)=Aname.now
    K=2^j
    data_e=data[rep(seq_len(nrow(data)), times = K),]
    for (k in (1:K)) {
      m=(n*(k-1)+1):(n*k)
      for (h in (1:j)) { data_e[m,Aname.now[h]]=Alist[k,h] }
    }
    if (j<n.steps) {
      EY=predict(Outcome.m.list[[j]],newdata=data_e)
      #EY2=sigma2.list[[j]]+EY^2
      data_e[,Outcome.name]=EY
      Outcome.m.list[[j+1]] = lm(Outcome.formula[[n.steps-(j+1)+1]],data=data_e)
      coef.list[[j+1]] = coef(Outcome.m.list[[j+1]])
      my=predict(Outcome.m.list[[j+1]],newdata=data_e)
      data_e$Sigma2 = EY^2 + predict(Var.m.list[[j]],newdata=data_e,type="response") - 2*EY*my +my^2
      Var.m.list[[j+1]] = glm(as.formula(paste("Sigma2 ~",as.character(Var.formula[[n.steps-(j+1)+1]])[2])),
                              family = quasipoisson(link = "log"),data=data_e)
      sigma2.list[[j+1]] = Var.m.list[[j+1]]$coefficients
    } else {
      EY=predict(Outcome.m.list[[j]],newdata=data_e)
      Sigma2=predict(Var.m.list[[j]],newdata=data_e,type="response")
      Xmat = model.matrix(MSQM.formula,data=data_e)
      n.theta=dim(Xmat)[2]
      est.eq=function(theta) {
        y=as.numeric(Xmat %*% theta)
        Pval=pnorm(y,mean=EY,sd=sqrt(Sigma2))
        out=rep(NA,n.theta)
        for (m in (1:n.theta)) {
          out[m] = sum(Xmat[,m]*(Pval-q))
        }
        out
      }
      theta.est=nleqslv::nleqslv(x=start.value,fn=est.eq)$x
    }
  }
  return(list(Point=theta.est,delta=coef.list,eta=sigma2.list))
}

#' ICR estimator (standard error)
ICR.f.SE = function(data,q,
                    Outcome.formula,
                    Var.formula,
                    MSQM.formula,
                    A.names,Outcome.name,delta,eta,theta,tau_n) {
  n.steps=length(A.names)
  n=dim(data)[1]
  beta_sep=cumsum(c(0,c(rbind(sapply(delta,length),sapply(eta,length)))))
  beta = c(unlist(do.call(c, Map(list, delta, eta))))
  # create the extended datasets and Outcome.m.list and Var.m.list
  data_e.list =list()
  data_e.list[[1]] = data
  for (j in (1:n.steps)) {
    Alist=list()
    for (i in (1:j)) { Alist[[i]]=c(0,1) }
    Alist=expand.grid(Alist);
    Aname.now=A.names[(n.steps-j+1):n.steps]
    colnames(Alist)=Aname.now
    K=2^j
    data_e=data[rep(seq_len(nrow(data)), times = K),]
    for (k in (1:K)) {
      m=(n*(k-1)+1):(n*k)
      for (h in (1:j)) { data_e[m,Aname.now[h]]=Alist[k,h] }
    }
    data_e.list[[j+1]] = data_e
  }
  ######## estimating equation of beta
  # type 2 meat, type 3 bread
  est.eq.beta = function(beta,type=2) {
    #out = matrix(NA,ncol=length(beta),nrow=n)
    for (j in (1:n.steps)) {
      if (j==1) {
        delta.now=beta[(beta_sep[2*j-1]+1):(beta_sep[2*j])];eta.now=beta[(beta_sep[2*j]+1):(beta_sep[2*j+1])]
        # estimating equation for delta
        X =model.matrix(Outcome.formula[[n.steps-j+1]],data=data_e.list[[j]])
        Y = data_e.list[[j]][,Outcome.name]
        resi= c(Y - c(X %*% delta.now))
        eq.delta.now = X * resi
        # estimating equation for eta
        data_e.list[[j]]$Sigma2 = resi^2
        Xv = model.matrix(as.formula(paste("Sigma2 ~",as.character(Var.formula[[n.steps]])[2])),data=data_e.list[[j]])
        Yv = data_e.list[[j]][,"Sigma2"]
        resiv = c(Yv - exp(Xv %*% eta.now))
        eq.eta.now = Xv * resiv
        out = cbind(eq.delta.now,eq.eta.now)
      } else {
        X = model.matrix(Outcome.formula[[n.steps-j+2]],data=data_e.list[[j]])
        EY = c(X %*% delta.now)
        Xv = model.matrix(as.formula(Var.formula[[n.steps-j+2]]),data=data_e.list[[j]])
        EYv = exp(Xv %*% eta.now)
        data_e.list[[j]][,Outcome.name] = EY
        delta.now=beta[(beta_sep[2*j-1]+1):(beta_sep[2*j])];eta.now=beta[(beta_sep[2*j]+1):(beta_sep[2*j+1])]
        # estimating equation for delta
        X =model.matrix(Outcome.formula[[n.steps-j+1]],data=data_e.list[[j]])
        Y = data_e.list[[j]][,Outcome.name]
        resi= c(Y - c(X %*% delta.now))
        eq.delta.now = t(X * resi)
        dim(eq.delta.now)=c(dim(eq.delta.now)[1],n,2^(j-1))
        eq.delta.now = aperm(eq.delta.now,c(2,1,3))
        eq.delta.now = rowSums(eq.delta.now, dims = 2)
        # estimating equation for eta
        my= c(X %*% delta.now)
        data_e.list[[j]]$Sigma2 = EY^2 + EYv - 2*EY*my +my^2
        Xv = model.matrix(as.formula(paste("Sigma2 ~",as.character(Var.formula[[n.steps-j+1]])[2])),data=data_e.list[[j]])
        Yv = data_e.list[[j]][,"Sigma2"]
        resiv = c(Yv - exp(Xv %*% eta.now))
        eq.eta.now = t(Xv * resiv)
        dim(eq.eta.now)=c(dim(eq.eta.now)[1],n,2^(j-1))
        eq.eta.now = aperm(eq.eta.now,c(2,1,3))
        eq.eta.now = rowSums(eq.eta.now, dims = 2)
        out = cbind(out,eq.delta.now,eq.eta.now)
      }
    }
    if (type==2) {return(out)}
    if (type==3) {return(apply(out,2,sum))}
  }
  ## U_beta
  Ub = est.eq.beta(beta,type=2)
  ## I_theta and S_theta_alpha
  I_beta = numDeriv::jacobian(est.eq.beta, x=beta, method="simple", type=3)/n
  I_beta1 = solve(I_beta,tol=10^(-23))[(beta_sep[2*n.steps-1]+1):beta_sep[length(beta_sep)],]
  ####### estimating equation of theta
  param_sep= c(0,cumsum(c(length(theta),length(delta[[n.steps]]),length(eta[[n.steps]]))))
  param = c(theta,delta[[n.steps]],eta[[n.steps]])
  # estimating equation of theta
  # type 2 meat, type 3 bread
  Xmat = model.matrix(MSQM.formula,data=data_e.list[[n.steps+1]])
  est.eq=function(theta,type=2) {
    delta.now = theta[(param_sep[2]+1):param_sep[3]]
    eta.now = theta[(param_sep[3]+1):param_sep[4]]
    theta = theta[(param_sep[1]+1):param_sep[2]]
    X = model.matrix(Outcome.formula[[1]],data=data_e.list[[n.steps+1]])
    EY = c(X %*% delta.now)
    Xv = model.matrix(as.formula(Var.formula[[1]]),data=data_e.list[[n.steps+1]])
    EYv = exp(Xv %*% eta.now)
    y=as.numeric(Xmat %*% theta)
    Pval=pnorm(y,mean=EY,sd=sqrt(EYv))
    eq.theta.now = t(Xmat * (Pval-q))
    dim(eq.theta.now)=c(dim(eq.theta.now)[1],n,2^n.steps)
    eq.theta.now = aperm(eq.theta.now,c(2,1,3))
    out = rowSums(eq.theta.now, dims = 2)
    if (type==2) {return(out)}
    if (type==3) {return(apply(out,2,sum))}
  }
  ### sandwich variance
  ## U_theta
  Ut = est.eq(theta=param,type=2)
  ## I_theta and S_theta_beta
  I_esteq = numDeriv::jacobian(est.eq, x=param, method="simple", type=3)/n
  I_theta = I_esteq[,(param_sep[1]+1):param_sep[2]]
  S_theta_beta = I_esteq[,-c((param_sep[1]+1):param_sep[2])]
  ## calculate the sandwich variance
  # meat
  V_icr = crossprod(Ut - t((S_theta_beta) %*% I_beta1 %*% t(Ub)))/n
  theta.vcov = (solve(I_theta) %*% V_icr %*% t(solve(I_theta)))/n
  return(list(SE=sqrt(abs(diag(theta.vcov))),
              SE.TE = sqrt(abs((t(c(0,rep(1,n.steps))) %*% theta.vcov %*% c(0,rep(1,n.steps)))[1,1])),
              Ub=Ub,
              I_beta=I_beta,
              SE.mat = theta.vcov))
}

#' DR estimator (point estimation)
DR.f.point=function(data,q,
                    Outcome.formula,
                    PS.formula,S.formula,
                    Var.formula=Var.formula,MSQM.formula,
                    A.names,Outcome.name,start.value,tau_n) {
  n.steps=length(PS.formula)
  n=dim(data)[1]
  # Obtain PS estimates
  PS.m.list=S.m.list=list()
  for (j in (1:n.steps)) {
    PS.m.list[[j]] = glm(PS.formula[[j]],family=binomial(link="logit"),data=data)
    S.m.list[[j]] = glm(S.formula[[j]],family=binomial(link="logit"),data=data)
  }
  # Obtain Outcome model estimates
  coef.list=sigma2.list=Outcome.m.list=Var.m.list=list()
  Outcome.m.list[[1]] = lm(Outcome.formula[[n.steps]],data=data)
  coef.list[[1]]=coef(Outcome.m.list[[1]]);
  EY=predict(Outcome.m.list[[1]],newdata=data)
  data$Sigma2 = (Outcome.m.list[[1]]$residuals)^2
  Var.m.list[[1]] = glm(as.formula(paste("Sigma2 ~",as.character(Var.formula[[n.steps]])[2])),family=quasipoisson(link="log"),data=data)
  sigma2.list[[1]]=Var.m.list[[1]]$coefficients
  data_e.list=list();data_e.list[[1]]=data
  for (j in (1:(n.steps))) {
    # generate the extended dataset
    Alist=list()
    for (i in (1:j)) { Alist[[i]]=c(0,1) }
    Alist=expand.grid(Alist);
    Aname.now=A.names[(n.steps-j+1):n.steps]
    colnames(Alist)=Aname.now
    K=2^j
    data_e=data[rep(seq_len(nrow(data)), times = K),]
    for (k in (1:K)) {
      m=(n*(k-1)+1):(n*k)
      for (h in (1:j)) { data_e[m,Aname.now[h]]=Alist[k,h] }
    }
    data_e.list[[j+1]]=data_e
    if (j<n.steps) {
      EY=predict(Outcome.m.list[[j]],newdata=data_e)
      #EY2=sigma2.list[[j]]+EY^2
      data_e[,Outcome.name]=EY
      Outcome.m.list[[j+1]] = lm(Outcome.formula[[n.steps-(j+1)+1]],data=data_e)
      coef.list[[j+1]] = coef(Outcome.m.list[[j+1]])
      my=predict(Outcome.m.list[[j+1]],newdata=data_e)
      data_e$Sigma2 = EY^2 + predict(Var.m.list[[j]],newdata=data_e,type="response") - 2*EY*my +my^2
      Var.m.list[[j+1]] = glm(as.formula(paste("Sigma2 ~",as.character(Var.formula[[n.steps-(j+1)+1]])[2])),
                              family=quasipoisson(link="log"),data=data_e)
      sigma2.list[[j+1]]=Var.m.list[[j+1]]$coefficients
    }
  }
  # Calculate the IPW and conditional mean for each step
  sw.list=s.list=list();mean.list.r=mean.list.l=list()
  var.list.r=var.list.l=list()
  for (j in (1:n.steps)) {
    # calculate the IPW
    S.mat=matrix(NA,ncol=n.steps,nrow=dim(data_e.list[[j]]))
    IPW.mat = matrix(NA,ncol=n.steps-j+1,nrow=dim(data_e.list[[j]]))
    for (m in (1:(n.steps-j+1))) {
      ps.now = predict(PS.m.list[[m]],newdata=data_e.list[[j]],type="response")
      IPW.mat[,m]=(ps.now^(data_e.list[[j]][,A.names[m]]))*((1-ps.now)^(1-data_e.list[[j]][,A.names[m]]))
    }
    for (m in (1:(n.steps))) {
      s.now = predict(S.m.list[[m]],newdata=data_e.list[[j]],type="response")
      S.mat[,m]=(s.now^(data_e.list[[j]][,A.names[m]]))*((1-s.now)^(1-data_e.list[[j]][,A.names[m]]))
    }
    s.list[[j]] = apply(S.mat,1,prod)
    sw.list[[j]]=1/apply(IPW.mat,1,prod)
    # calculate the conditional mean and variance
    mean.list.r[[j]]=predict(Outcome.m.list[[j]],newdata=data_e.list[[j]])
    mean.list.l[[j]]=predict(Outcome.m.list[[j]],newdata=data_e.list[[j+1]])
    var.list.r[[j]]=predict(Var.m.list[[j]],newdata=data_e.list[[j]],type="response")
    var.list.l[[j]]=predict(Var.m.list[[j]],newdata=data_e.list[[j+1]],type="response")
  }
  S.mat=matrix(NA,ncol=n.steps,nrow=dim(data_e.list[[n.steps+1]]))
  for (m in (1:(n.steps))) {
    s.now = predict(S.m.list[[m]],newdata=data_e.list[[n.steps+1]],type="response")
    S.mat[,m]=(s.now^(data_e.list[[n.steps+1]][,A.names[m]]))*((1-s.now)^(1-data_e.list[[n.steps+1]][,A.names[m]]))
  }
  s.list[[n.steps+1]] = apply(S.mat,1,prod)
  # Construct the estimating equation
  Xmat.list=list()
  for (j in (1:(n.steps+1))) {
    Xmat.list[[j]]=model.matrix(MSQM.formula,data=data_e.list[[j]])
  }
  n.theta=dim(Xmat.list[[1]])[2]
  est.eq=function(theta) {
    expit=function(x,m,gamma=tau_n) {plogis(m-x,scale=gamma)}
    y=as.numeric(Xmat.list[[1]] %*% theta)
    Pval=expit(x=data[,Outcome.name],m=y)
    Fval=pnorm(y,mean=mean.list.r[[1]],sd=sqrt(var.list.r[[1]]))
    out=rep(NA,n.theta)
    for (m in (1:n.theta)) {
      out[m] = sum(Xmat.list[[1]][,m]*s.list[[1]]*sw.list[[1]]*(Pval-Fval))
    }
    for (j in (1:(n.steps-1))) {
      y=as.numeric(Xmat.list[[j+1]] %*% theta)
      Pval=pnorm(y,mean=mean.list.l[[j]],sd=sqrt(var.list.l[[j]]))
      Fval=pnorm(y,mean=mean.list.r[[j+1]],sd=sqrt(var.list.r[[j+1]]))
      for (m in (1:n.theta)) {
        out[m] = out[m] + sum(Xmat.list[[j+1]][,m]*s.list[[j+1]]*sw.list[[j+1]]*(Pval-Fval))
      }
    }
    y=as.numeric(Xmat.list[[n.steps+1]] %*% theta)
    Pval=pnorm(y,mean=mean.list.l[[n.steps]],sd=sqrt(var.list.l[[n.steps]]))
    for (m in (1:n.theta)) {
      out[m] = out[m]+sum(Xmat.list[[n.steps+1]][,m]*s.list[[j+2]]*(Pval-q))
    }
    out
  }
  theta.est=nleqslv::nleqslv(x=start.value,fn=est.eq)$x
  return(list(Point=theta.est))
}


#' DR estimator (standard error)
DR.f.SE=function(data,q,
                 Outcome.formula,
                 PS.formula,S.formula,
                 Var.formula=Var.formula,MSQM.formula,
                 A.names,Outcome.name,tau_n,
                 theta,alpha,alphas,delta,eta,
                 Ua,Ia,Ub,Ib) {
  n.steps=length(PS.formula)
  n=dim(data)[1]
  # generate the extended dataset
  data_e.list=list();data_e.list[[1]]=data
  for (j in (1:(n.steps))) {
    Alist=list()
    for (i in (1:j)) { Alist[[i]]=c(0,1) }
    Alist=expand.grid(Alist);
    Aname.now=A.names[(n.steps-j+1):n.steps]
    colnames(Alist)=Aname.now
    K=2^j
    data_e=data[rep(seq_len(nrow(data)), times = K),]
    for (k in (1:K)) {
      m=(n*(k-1)+1):(n*k)
      for (h in (1:j)) { data_e[m,Aname.now[h]]=Alist[k,h] }
    }
    data_e.list[[j+1]]=data_e
  }
  # all parameters for the MSQM, PS model and Outcome regression models
  param_sep=cumsum(c(0,length(theta),c(rbind(sapply(alpha,length),sapply(alphas,length),
                                             sapply(delta,length),sapply(eta,length)))))
  param = c(theta,unlist(do.call(c, Map(list, alpha,alphas,delta, eta))))
  # estimating equation of theta
  # type 2 meat, type 3 bread
  est.eq=function(theta,type=2) {
    # restore the parameters
    alpha=alphas=delta=eta=list()
    for (j in (1:n.steps)) {
      alpha[[j]] = theta[(param_sep[4*(j-1)+2]+1):param_sep[4*(j-1)+3]]
      alphas[[j]] = theta[(param_sep[4*(j-1)+3]+1):param_sep[4*(j-1)+4]]
      delta[[j]] = theta[(param_sep[4*(j-1)+4]+1):param_sep[4*(j-1)+5]]
      eta[[j]] = theta[(param_sep[4*(j-1)+5]+1):param_sep[4*(j-1)+6]]
    }
    delta = rev(delta);eta=rev(eta)
    theta = theta[1:param_sep[2]]
    # Calculate the IPW and conditional mean for each step
    sw.list=s.list=list();mean.list.r=mean.list.l=list()
    var.list.r=var.list.l=list()
    for (j in (1:n.steps)) {
      # calculate the IPW
      S.mat=matrix(NA,ncol=n.steps,nrow=dim(data_e.list[[j]]))
      IPW.mat = matrix(NA,ncol=n.steps-j+1,nrow=dim(data_e.list[[j]]))
      for (m in (1:(n.steps-j+1))) {
        ps.now = plogis(c(model.matrix(PS.formula[[m]],data=data_e.list[[j]]) %*% alpha[[m]]))
        IPW.mat[,m]=(ps.now^(data_e.list[[j]][,A.names[m]]))*((1-ps.now)^(1-data_e.list[[j]][,A.names[m]]))
      }
      for (m in (1:(n.steps))) {
        s.now = plogis(c(model.matrix(S.formula[[m]],data=data_e.list[[j]]) %*% alphas[[m]]))
        S.mat[,m]=(s.now^(data_e.list[[j]][,A.names[m]]))*((1-s.now)^(1-data_e.list[[j]][,A.names[m]]))
      }
      s.list[[j]] = apply(S.mat,1,prod)
      sw.list[[j]]=1/apply(IPW.mat,1,prod)
      # calculate the conditional mean and variance
      mean.list.r[[j]]= c(model.matrix(Outcome.formula[[n.steps-j+1]],data=data_e.list[[j]]) %*% delta[[n.steps-j+1]])
      mean.list.l[[j]]= c(model.matrix(Outcome.formula[[n.steps-j+1]],data=data_e.list[[j+1]]) %*% delta[[n.steps-j+1]])
      var.list.r[[j]] = exp(model.matrix(Var.formula[[n.steps-j+1]],data=data_e.list[[j]]) %*% eta[[n.steps-j+1]])
      var.list.l[[j]]=  exp(model.matrix(Var.formula[[n.steps-j+1]],data=data_e.list[[j+1]]) %*% eta[[n.steps-j+1]])
    }
    S.mat=matrix(NA,ncol=n.steps,nrow=dim(data_e.list[[n.steps+1]]))
    for (m in (1:(n.steps))) {
      s.now = plogis(c(model.matrix(S.formula[[m]],data=data_e.list[[n.steps+1]]) %*% alphas[[m]]))
      S.mat[,m]=(s.now^(data_e.list[[n.steps+1]][,A.names[m]]))*((1-s.now)^(1-data_e.list[[n.steps+1]][,A.names[m]]))
    }
    s.list[[n.steps+1]] = apply(S.mat,1,prod)
    # Construct the estimating equation
    Xmat.list=list()
    for (j in (1:(n.steps+1))) {
      Xmat.list[[j]]=model.matrix(MSQM.formula,data=data_e.list[[j]])
    }
    n.theta=param_sep[2]
    # construct the estimating equation
    expit=function(x,m,gamma=tau_n) {plogis(m-x,scale=gamma)}
    y=as.numeric(Xmat.list[[1]] %*% theta)
    Pval=expit(x=data[,Outcome.name],m=y)
    Fval=pnorm(y,mean=mean.list.r[[1]],sd=sqrt(var.list.r[[1]]))
    est.theta.eq= Xmat.list[[1]] * s.list[[1]]*sw.list[[1]]* (Pval-Fval)
    for (j in (1:(n.steps-1))) {
      y=as.numeric(Xmat.list[[j+1]] %*% theta)
      Pval=pnorm(y,mean=mean.list.l[[j]],sd=sqrt(var.list.l[[j]]))
      Fval=pnorm(y,mean=mean.list.r[[j+1]],sd=sqrt(var.list.r[[j+1]]))
      eq.theta.now = t(Xmat.list[[j+1]] * s.list[[j+1]]*sw.list[[j+1]]* (Pval-Fval))
      dim(eq.theta.now)=c(dim(eq.theta.now)[1],n,2^j)
      eq.theta.now = aperm(eq.theta.now,c(2,1,3))
      est.theta.eq = est.theta.eq + rowSums(eq.theta.now, dims = 2)
    }
    y=as.numeric(Xmat.list[[n.steps+1]] %*% theta)
    Pval=pnorm(y,mean=mean.list.l[[n.steps]],sd=sqrt(var.list.l[[n.steps]]))
    eq.theta.now = t(Xmat.list[[n.steps+1]] * s.list[[n.steps+1]]* (Pval-q))
    dim(eq.theta.now)=c(dim(eq.theta.now)[1],n,2^(n.steps))
    eq.theta.now = aperm(eq.theta.now,c(2,1,3))
    est.theta.eq = est.theta.eq + rowSums(eq.theta.now, dims = 2)
    if (type==2) {return(est.theta.eq)}
    if (type==3) {return(apply(est.theta.eq,2,sum))}
  }
  ### sandwich variance
  ## U_theta
  Ut = est.eq(theta=param,type=2)
  ## I_theta and S_theta_beta
  I_esteq = numDeriv::jacobian(est.eq, x=param, method="simple", type=3)/n
  I_theta = I_esteq[,(param_sep[1]+1):param_sep[2]]
  alpha_loc=beta_loc=c()
  for (j in (1:n.steps)) {
    alpha_loc = c(alpha_loc,(param_sep[4*(j-1)+1+1]+1):param_sep[4*(j-1)+1+3])
    beta_loc = c(beta_loc,(param_sep[4*(j-1)+1+3]+1):param_sep[4*(j-1)+1+5])
  }
  S_theta_alpha = I_esteq[,alpha_loc]
  S_theta_beta  = I_esteq[,beta_loc]
  ## calculate the sandwich variance
  # meat
  V_dr = crossprod(Ut - t((S_theta_beta) %*% solve(Ib,tol=10^(-23)) %*% t(Ub)) + t(S_theta_alpha %*% Ia %*% t(Ua))   )/n
  theta.vcov = (solve(I_theta,tol=10^(-23)) %*% V_dr %*% t(solve(I_theta,tol=10^(-23))))/n
  return(list(SE=sqrt(abs(diag(theta.vcov))),
              SE.TE = sqrt(abs((t(c(0,rep(1,n.steps))) %*% theta.vcov %*% c(0,rep(1,n.steps)))[1,1])),
              SE.mat = theta.vcov
  ))
}


#' Estimate the Marginal Structural Quantile Model (MSQM)
#'
#' Estimate and inference on MSQM based on four approaches (IPW, ICR, and DR)
#'
#' @param data a data frame
#' @param q desired lower quantile (ranged between 0 and 1)
#' @param PS.formula a list of formulas for the propensity score models in each time period
#' @param S.formula a list of formulas for the stabilized weights in each time period
#' @param Outcome.formula a list of formulas for mean structure of the outcome model in each time period
#' @param Var.formula a list of formulas for the variance structure of the outcome model in each time period
#' @param MSQM.formula Specification of the MSQM model
#' @param A.names the column name of treatment from period 1 to period K
#' @param Outcome.names the column name of the outcome (measured at the end of period K)
#' @return Point, standard error, and confidence interval estimates based on the inverse
#' probability weighting (IPW) approach, the iterative conditional regression (ICR) approach, and the
#' doubly robust (DR) approach.
#' @references Cheng, C., Hu, L., & Li, F. (2022). Doubly robust estimation and sensitivity analysis for marginal structural quantile models. Accepted by Biometrics (arXiv preprint arXiv:2210.04100).
#' @examples
#' #### Step 1: generate data
#' set.seed(12345)
#' data=data.gen(n=500)
#'
#' head(round(data,2))
#' #     X11   X12   X21  X22   X31   X32 A1 A2 A3      Y
#' # 1 -1.21 -0.44  0.23 1.65 -0.03  0.87  1  1  0 -10.84
#' # 2  1.79  1.23  1.30 0.97  2.35  5.36  1  1  0  11.59
#' # 3  0.81  1.29  1.05 1.09  0.86  2.30  1  1  1   0.43
#' # 4 -0.03 -1.84  2.63 2.77 -0.08  1.87  1  1  0  -4.04
#' # 5 -3.88 -0.57  1.37 1.05 -0.48  0.09  1  0  0  -3.59
#' # 6 -1.07  1.43 -0.67 0.76  0.57 -1.08  1  0  1 -17.87
#' # Note: (X11, X12): baseline covariates; A1: treatment in period 1; (X21, X22): covariates at beginning of period 2;
#' # A2: treatment in period 2; (X31, X32): covariates at beginning of period 2; A3: treatment in period 3;
#' # Y: outcome measured after A3
#'
#' #### Step 2: Specification
#' MSQM.formula=as.formula(Y~A1+A2+A3) # the MSQM model
#' q=0.5 # the desired quantile of interest (median here)
#' A.names=c("A1","A2","A3") # treatment variables
#' Outcome.name="Y" # outcome variable
#'
#' # specification of the stabilized weights (k=1,2,3 sequentially)
#' S.formula=list(as.formula(A1~1),
#'                as.formula(A2~A1),
#'                as.formula(A3~A2))
#' # specification of the principal score model (k=1,2,3 sequentially)
#' PS.formula=list(as.formula(A1~X11+X12),
#'                 as.formula(A2~A1+X21+X22),
#'                 as.formula(A3~A2+X31+X32))
#' # specification of outcome mean (k=1,2,3 sequentially)
#' Outcome.formula=list(as.formula(Y~A1+A2+A3+X11+X12),
#'                      as.formula(Y~A1+A2+A3+X11+X12+X21+X22),
#'                      as.formula(Y~A1+A2+A3+X11+X12+X21+X22+X31+X32))
#' # specification of outcome variance (k=1,2,3 sequentially)
#' Var.formula =list(as.formula(~A1+A2+A3),
#'                   as.formula(~A1+A2+A3),
#'                   as.formula(~A1+A2+A3))
#'
#' #### Step 3: fit the MSQM
#' obj = msqm.fit(data,q,PS.formula,S.formula,Outcome.formula,Var.formula,
#'          MSQM.formula,A.names,Outcome.name)
#' # point estimates, SE, and 95% CI of the MSQM coefficients
#' obj$MSQM_coef
#' # the variance-covariance matrix of the MSQM coefficients
#' obj$MSQM_vcov
msqm.fit=function(data,q,
                    PS.formula,S.formula,
                    Outcome.formula,Var.formula,
                    MSQM.formula,A.names,Outcome.name) {

  # naive estimator
  naive.e = naive.f(data,q,MSQM.formula)
  tau_n = naive.e$tau*dim(data)[1]^(-0.26) #naive.e$tau*
  start.value=naive.e$Point
  # IPW estimator
  ipw.e=IPW.f(data,q,PS.formula,S.formula,MSQM.formula,Outcome.name,start.value,tau_n)
  start.value=ipw.e$Point
  # ICR estimator
  icr.e = ICR.f.point(data,q,Outcome.formula,Var.formula,MSQM.formula,A.names,
                      Outcome.name,start.value,icr_coef=0)
  delta = icr.e$delta
  eta = icr.e$eta
  theta=icr.e$Point
  theta.icr.se=ICR.f.SE(data,q, Outcome.formula,Var.formula,MSQM.formula,A.names,Outcome.name,delta,eta,theta)
  # DR estimator
  dr.e = DR.f.point(data,q,Outcome.formula,PS.formula,S.formula,Var.formula=Var.formula,MSQM.formula,
                    A.names,Outcome.name,start.value,tau_n)
  theta=dr.e$Point
  alpha=ipw.e$alpha
  alphas = ipw.e$alphas
  delta = icr.e$delta
  eta = icr.e$eta
  Ua = ipw.e$Ua
  Ia = ipw.e$Ia
  Ub = theta.icr.se$Ub
  Ib = theta.icr.se$I_beta
  theta.dr.se=  DR.f.SE(data,q,Outcome.formula,PS.formula,S.formula,Var.formula=Var.formula,MSQM.formula,
                        A.names,Outcome.name,tau_n,theta,alpha,alphas,delta,eta,Ua,Ia,Ub,Ib)


  IPW.res=rbind(ipw.e$Point,ipw.e$SE,ipw.e$Point-1.96*ipw.e$SE,ipw.e$Point+1.96*ipw.e$SE)
  ICR.res=rbind(icr.e$Point,theta.icr.se$SE,icr.e$Point-1.96*theta.icr.se$SE,icr.e$Point+1.96*theta.icr.se$SE)
  DR.res =rbind(dr.e$Point,theta.dr.se$SE,dr.e$Point-1.96*theta.dr.se$SE,dr.e$Point+1.96*theta.dr.se$SE)

  rownames(IPW.res)=rownames(ICR.res)=rownames(DR.res)=c("point","SE","CI.lower","CI.upper")
  MSQM_coef = list(IPW=IPW.res,ICR=ICR.res,DR=DR.res)
  MSQM_vcov = list(IPW=ipw.e$SE.mat,ICR=theta.icr.se$SE.mat,DR=theta.dr.se$SE.mat)
  return(list(MSQM_coef=MSQM_coef,MSQM_vcov=MSQM_vcov))
}
