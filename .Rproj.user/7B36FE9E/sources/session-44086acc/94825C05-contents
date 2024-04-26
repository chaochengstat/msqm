####################################################################
#
# Functions for Case 1
#
####################################################################


#' Mediation analysis function for continuous outcome and mediator
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @return A list containing NIE, NDE and MP point and interval estimates.
mediate_contY_contM=function(data,
                             outcome="Y",
                             mediator="M",
                             exposure="X",
                             covariateY=c("X1","X2"),
                             covariateM=c("X1","X2"),x0=0,x1=1) {


  #data=dat1;outcome=Outcome;mediator=Mediator;covariateY=covariateM=CovariateY;exposure=Exposure

  data = as.data.frame(data)
  ## (1) formula for outcome model and mediator model
  if (is.null(covariateY)) {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
  } else {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
  }
  if (is.null(covariateM)) {
    formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
  } else {
    formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
  }
  ## (2) estimate the outcome model and mediator model
  model_Y=summary(lm(formula_Y,data=data))
  model_M=summary(lm(formula_M,data=data))

  beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
  gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled

  ## (3) covariance matrix of (beta,gamma)
  nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
  S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
  S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
  colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

  ## (4) approximate method
  ##### NIE TE and PM expressions
  if (is.null(covariateY)==0) {
    names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
    beta_c=paste("beta_",covariateY,sep="")
  }
  if (is.null(covariateM)==0) {
    names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
    gamma_c=paste("gamma_",covariateM,sep="")
  }

  NIEa_fun = function() {
    output = "beta2*gamma1*(x1-x0)"
    return(output)
  }
  variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
  NIEa_D=deriv(parse(text=NIEa_fun()),variable)
  gamma0=gamma[1];gamma1=gamma[2];
  if(is.null(covariateM)==0) {
    for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
  }
  beta0=beta[1];beta1=beta[2];beta2=beta[3]
  if(is.null(covariateY)==0) {
    for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
  }

  TEa_fun = function() {
    output = "(beta2*gamma1+beta1)*(x1-x0)"
    return(output)
  }
  TEa_D=deriv(parse(text=TEa_fun()),variable)

  PMa_fun = function() {
    .UP = "beta2*gamma1"
    .BOT = "beta2*gamma1+beta1"
    output=paste("(",.UP,")/(",.BOT,")")
    return(output)
  }
  PMa_D=deriv(parse(text=PMa_fun()),variable)


  NIEa_D = eval(NIEa_D)
  NIEa_p = NIEa_D[1]
  lambda= t(attr(NIEa_D,"gradient"))
  V_NIEa = as.vector(t(lambda) %*% S %*% lambda)

  TEa_D = eval(TEa_D)
  TEa_p = TEa_D[1]
  lambda= t(attr(TEa_D,"gradient"))
  V_TEa = as.vector(t(lambda) %*% S %*% lambda)

  PMa_D = eval(PMa_D)
  PMa_p = PMa_D[1]
  lambda= t(attr(PMa_D,"gradient"))
  V_PMa = as.vector(t(lambda) %*% S %*% lambda)


  point_est = c(NIEa_p,TEa_p,PMa_p);
  names(point_est)=c("NIE","TE","PM")
  var_est = c(V_NIEa,V_TEa,V_PMa);
  names(var_est)=c("NIE","TE","PM")
  sd_est = sqrt(var_est)
  names(sd_est)=c("NIE","TE","PM")
  ci_est = rbind(point_est-1.96*sd_est,point_est+1.96*sd_est)
  rownames(ci_est) = c("Lower boundary","Upper boundary")

  return(list(point_est=point_est,var_est=var_est,sd_est=sd_est,ci_est=ci_est))
}


#' Using bootstrap to calculate the confidence intervals in the scenario of continuous outcome and continuous mediator.
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param R The number of replications.
#' @return The 95% confidence interval by the bootstrap approach
Mediate_contY_contM_bootci=function(data,
                                    outcome="Y",
                                    mediator="M",
                                    exposure="X",
                                    covariateY=c("X1","X2"),
                                    covariateM=c("X1","X2"),
                                    x0=0,x1=1,R=1000) {

  data = as.data.frame(data)
  get_par_boot=function(data=data,indices) {
    data=data[indices,]
    ## (1) formula for outcome model and mediator model
    if (is.null(covariateY)) {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
    } else {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
    }
    if (is.null(covariateM)) {
      formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
    } else {
      formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
    }
    ## (2) estimate the outcome model and mediator model
    model_Y=summary(lm(formula_Y,data=data))
    model_M=summary(lm(formula_M,data=data))

    beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
    gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled
    ## (3) covariance matrix of (beta,gamma)
    nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
    S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
    S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
    colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

    ## (4) approximate method
    ##### NIE TE and PM expressions
    if (is.null(covariateY)==0) {
      names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
      beta_c=paste("beta_",covariateY,sep="")
    }
    if (is.null(covariateM)==0) {
      names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
      gamma_c=paste("gamma_",covariateM,sep="")
    }

    NIEa_fun = function() {
      output = "beta2*gamma1*(x1-x0)"
      return(output)
    }
    variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
    gamma0=gamma[1];gamma1=gamma[2];
    if(is.null(covariateM)==0) {
      for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
    }
    beta0=beta[1];beta1=beta[2];beta2=beta[3]
    if(is.null(covariateY)==0) {
      for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
    }

    TEa_fun = function() {
      output = "(beta2*gamma1+beta1)*(x1-x0)"
      return(output)
    }


    PMa_fun = function() {
      .UP = "beta2*gamma1*(x1-x0)"
      .BOT = "(beta2*gamma1+beta1)*(x1-x0)"
      output=paste("(",.UP,")/(",.BOT,")")
      return(output)
    }

    NIEa_p = eval(parse(text=NIEa_fun()))
    TEa_p = eval(parse(text=TEa_fun()))
    PMa_p = eval(parse(text=PMa_fun()))


    point_est = c(NIEa_p,TEa_p,PMa_p);
    names(point_est)=c("NIE","TE","PM")

    return(point_est)
  }


  boot.par=boot::boot(data=data, statistic=get_par_boot, R=R)
  boot.parciNIEa <- boot::boot.ci(boot.par, index=1, type=c("perc"))
  boot.parciTEa <- boot::boot.ci(boot.par, index=2, type=c("perc"))
  boot.parciPMa <- boot::boot.ci(boot.par, index=3, type=c("perc"))


  ci_est_prec = c(boot.parciNIEa$percent[4:5],
                  boot.parciTEa$percent[4:5],
                  boot.parciPMa$percent[4:5])
  names(ci_est_prec)=c(paste(rep("CI_",6),rep(c("NIE","TE","PM"),each=2),rep(c("_Low","_High"),times=3),sep=""))
  return(ci_est_prec)
}


####################################################################
#
# Functions for Case 2
#
####################################################################

#' Mediation analysis function for continuous outcome and binary mediator
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param cY conditional levels of covariateY
#' @param cM conditional levels of covariateM
#' @return A list containing NIE, NDE and MP point and interval estimates.
mediate_contY_binaM=function(data,
                             outcome="Y",
                             mediator="M",
                             exposure="X",
                             covariateY=c("X1","X2"),
                             covariateM=c("X1","X2"),
                             x0=0,x1=1,cY=c(0,0),cM=c(0,0)) {


  #data=dat1;outcome=Outcome;mediator=Mediator;covariateY=covariateM=CovariateY;exposure=Exposure

  data = as.data.frame(data)
  ## (1) formula for outcome model and mediator model
  if (is.null(covariateY)) {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
  } else {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
  }
  if (is.null(covariateM)) {
    formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
  } else {
    formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
  }
  ## (2) estimate the outcome model and mediator model
  model_Y=summary(lm(formula_Y,data=data))
  model_M=summary(glm(formula_M,family=binomial(link="logit"),data=data))

  beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
  gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled

  ## (3) covariance matrix of (beta,gamma)
  nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
  S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
  S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
  colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

  ## (4) approximate method
  ##### NIE TE and PM expressions
  if (is.null(covariateY)==0) {
    names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
    beta_c=paste("beta_",covariateY,sep="")
  }
  if (is.null(covariateM)==0) {
    names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
    gamma_c=paste("gamma_",covariateM,sep="")
  }
  .A = paste("exp(gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .B= paste("exp(gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")

  NIEa_fun = function() {
    output=paste("beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
    return(output)
  }
  variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
  NIEa_D=deriv(parse(text=NIEa_fun()),variable)
  gamma0=gamma[1];gamma1=gamma[2];
  if(is.null(covariateM)==0) {
    for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
  }
  beta0=beta[1];beta1=beta[2];beta2=beta[3]
  if(is.null(covariateY)==0) {
    for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
  }

  TEa_fun = function() {
    output = paste("beta1*",x1-x0,"+","beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
    return(output)
  }
  TEa_D=deriv(parse(text=TEa_fun()),variable)

  PMa_fun = function() {
    .UP = paste("beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
    .BOT = paste("beta1*",x1-x0,"+","beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
    output=paste("(",.UP,")/(",.BOT,")")
    return(output)
  }
  PMa_D=deriv(parse(text=PMa_fun()),variable)


  NIEa_D = eval(NIEa_D)
  NIEa_p = NIEa_D[1]
  lambda= t(attr(NIEa_D,"gradient"))
  V_NIEa = as.vector(t(lambda) %*% S %*% lambda)

  TEa_D = eval(TEa_D)
  TEa_p = TEa_D[1]
  lambda= t(attr(TEa_D,"gradient"))
  V_TEa = as.vector(t(lambda) %*% S %*% lambda)

  PMa_D = eval(PMa_D)
  PMa_p = PMa_D[1]
  lambda= t(attr(PMa_D,"gradient"))
  V_PMa = as.vector(t(lambda) %*% S %*% lambda)

  point_est = c(NIEa_p,TEa_p,PMa_p);
  names(point_est)=c("NIE","TE","PM")
  var_est = c(V_NIEa,V_TEa,V_PMa);
  names(var_est)=c("NIE","TE","PM")
  sd_est = sqrt(var_est)
  names(sd_est)=c("NIE","TE","PM")
  ci_est = rbind(point_est-1.96*sd_est,point_est+1.96*sd_est)
  rownames(ci_est) = c("Lower boundary","Upper boundary")

  return(list(point_est=point_est,var_est=var_est,sd_est=sd_est,ci_est=ci_est))
}


#' Using bootstrap to calculate the confidence intervals in the scenario of continuous outcome and binary mediator.
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param cY conditional levels of covariateY
#' @param cM conditional levels of covariateM
#' @param R The number of replications.
#' @return The 95% confidence interval by the bootstrap approach
Mediate_contY_binaM_bootci=function(data,
                                    outcome="Y",
                                    mediator="M",
                                    exposure="X",
                                    covariateY=c("X1","X2"),
                                    covariateM=c("X1","X2"),
                                    x0=0,x1=1,cY=c(0,0),cM=c(0,0),R=1000) {
  #data=mydata;
  #outcome="Y";
  #mediator="M";
  #exposure="X";
  #covariateY=c("X1","X2","X3","X4","X5","X6","X7","X8");
  #covariateM=c("X1","X2","X3","X4","X5","X6","X7");
  #x0=0;x1=5;cY=rep(0,8);cM=rep(0,7)

  data = as.data.frame(data)
  get_par_boot=function(data=data,indices) {
    data=data[indices,]
    ## (1) formula for outcome model and mediator model
    if (is.null(covariateY)) {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
    } else {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
    }
    if (is.null(covariateM)) {
      formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
    } else {
      formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
    }
    ## (2) estimate the outcome model and mediator model
    model_Y=summary(lm(formula_Y,data=data))
    model_M=summary(glm(formula_M,family=binomial(link="logit"),data=data))

    beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
    gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled
    ## (3) covariance matrix of (beta,gamma)
    nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
    S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
    S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
    colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

    ## (4) approximate method
    ##### NIE TE and PM expressions
    if (is.null(covariateY)==0) {
      names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
      beta_c=paste("beta_",covariateY,sep="")
    }
    if (is.null(covariateM)==0) {
      names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
      gamma_c=paste("gamma_",covariateM,sep="")
    }
    .A = paste("exp(gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .B= paste("exp(gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")


    NIEa_fun = function() {
      output = paste("beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
      return(output)
    }
    variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
    gamma0=gamma[1];gamma1=gamma[2];
    if(is.null(covariateM)==0) {
      for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
    }
    beta0=beta[1];beta1=beta[2];beta2=beta[3]
    if(is.null(covariateY)==0) {
      for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
    }

    TEa_fun = function() {
      output = paste("beta1*",x1-x0,"+","beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
      return(output)
    }


    PMa_fun = function() {
      .UP = paste("beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
      .BOT = paste("beta1*",x1-x0,"+","beta2*((",.A,")/(1+",.A,")-(",.B,")/(1+",.B,"))")
      output=paste("(",.UP,")/(",.BOT,")")
      return(output)
    }

    NIEa_p = eval(parse(text=NIEa_fun()))
    TEa_p = eval(parse(text=TEa_fun()))
    PMa_p = eval(parse(text=PMa_fun()))


    point_est = c(NIEa_p,TEa_p,PMa_p);
    names(point_est)=c("NIE","TE","PM")

    return(point_est)
  }

  boot.par=boot::boot(data=data, statistic=get_par_boot, R=R)
  boot.parciNIEa <- boot::boot.ci(boot.par, index=1, type=c("perc"))
  boot.parciTEa <- boot::boot.ci(boot.par, index=2, type=c("perc"))
  boot.parciPMa <- boot::boot.ci(boot.par, index=3, type=c("perc"))


  ci_est_prec = c(boot.parciNIEa$percent[4:5],
                  boot.parciTEa$percent[4:5],
                  boot.parciPMa$percent[4:5])
  names(ci_est_prec)=c(paste(rep("CI_",6),rep(c("NIE","TE","PM"),each=2),rep(c("_Low","_High"),times=3),sep=""))
  return(ci_est_prec)
}


####################################################################
#
# Functions for Case 3
#
####################################################################



#' Mediation analysis function for binary outcome and continuous mediator
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param cY conditional levels of covariateY
#' @param cM conditional levels of covariateM
#' @return A list containing NIE, NDE and MP point and interval estimates.
mediate_binaY_contM=function(data,
                             outcome="Y",
                             mediator="M",
                             exposure="X",
                             covariateY=c("X1","X2"),
                             covariateM=c("X1","X2"),
                             x0=0,x1=1,cY=c(0,0),cM=c(0,0)) {

  data = as.data.frame(data)

  HermiteCoefs=function (order) {
    x <- 1
    if (order > 0)
      for (n in 1:order) x <- c(0, 2 * x) - c(((0:(n - 1)) *
                                                 x)[-1L], 0, 0)
      return(x)
  }

  gauss.hermite=function (f, mu = 0, sd = 1, ..., order = 5) {
    stopifnot(is.function(f))
    stopifnot(length(mu) == 1)
    stopifnot(length(sd) == 1)
    Hn <- HermiteCoefs(order)
    Hn1 <- HermiteCoefs(order - 1)
    x <- sort(Re(polyroot(Hn)))
    Hn1x <- matrix(Hn1, nrow = 1) %*% t(outer(x, 0:(order - 1),
                                              "^"))
    w <- 2^(order - 1) * factorial(order) * sqrt(pi)/(order *
                                                        Hn1x)^2
    ww <- w/sqrt(pi)
    xx <- mu + sd * sqrt(2) * x
    ans <- 0
    for (i in seq_along(x)) ans <- ans + ww[i] * f(xx[i], ...)
    return(ans)
  }

  mygrad=function (f, x0,heps = 1e-5, ...) {
    if (!is.numeric(x0))
      stop("Argument 'x0' must be a numeric value.")
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    p =length(f(x0))
    n <- length(x0)
    hh <- rep(0, n)
    gr <- matrix(0,nrow=n,ncol=p)
    for (i in 1:n) {
      hh[i] <- heps
      gr[i,] <- (f(x0 + hh) - f(x0 - hh))/(2 * heps)
      hh[i] <- 0
    }
    return(gr)
  }

  NIE_unbiased = function(theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                          x0s=x0s,x1=x1,cY=cY,cM=cM) {
    gamma0 = theta[1]
    gamma1 = theta[2]
    gamma_c = theta[loc_gamma_c]
    beta0  = theta[loc_beta_0[1]]
    beta1  = theta[loc_beta_0[2]]
    beta2  = theta[loc_beta_0[3]]
    beta_c = theta[loc_beta_c]
    sigma2 = theta[length(theta)]

    if (is.null(loc_beta_c)) {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
    } else {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
    }

    if (is.null(loc_gamma_c)) {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1,sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
    } else {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
    }

    output = log((p11s)/(1-p11s)) - log((p10s)/(1-p10s))
    return(output)
  }



  TE_unbiased = function(theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                         x0s=x0s,x1=x1,cY=cY,cM=cM) {
    gamma0 = theta[1]
    gamma1 = theta[2]
    gamma_c = theta[loc_gamma_c]
    beta0  = theta[loc_beta_0[1]]
    beta1  = theta[loc_beta_0[2]]
    beta2  = theta[loc_beta_0[3]]
    beta_c = theta[loc_beta_c]
    sigma2 = theta[length(theta)]

    if (is.null(loc_beta_c)) {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x)/(1+exp(beta0+beta1*x0s+beta2*x))}
    } else {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY)))}
    }

    if (is.null(loc_gamma_c)) {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1,sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
    } else {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
    }

    output = log((p11s)/(1-p11s)) - log((p00s)/(1-p00s))
    return(output)
  }

  PM_unbiased=function(theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                       x0s=x0s,x1=x1,cY=cY,cM=cM) {
    gamma0 = theta[1]
    gamma1 = theta[2]
    gamma_c = theta[loc_gamma_c]
    beta0  = theta[loc_beta_0[1]]
    beta1  = theta[loc_beta_0[2]]
    beta2  = theta[loc_beta_0[3]]
    beta_c = theta[loc_beta_c]
    sigma2 = theta[length(theta)]

    if (is.null(loc_beta_c)) {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x)/(1+exp(beta0+beta1*x0s+beta2*x))}
    } else {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY)))}
    }

    if (is.null(loc_gamma_c)) {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1,sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
    } else {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
    }

    te  = log((p11s)/(1-p11s)) - log((p00s)/(1-p00s))
    nie = log((p11s)/(1-p11s)) - log((p10s)/(1-p10s))
    return(nie/te)
  }

  ## (1) formula for outcome model and mediator model
  if (is.null(covariateY)) {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
  } else {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
  }
  if (is.null(covariateM)) {
    formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
  } else {
    formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
  }
  ## (2) estimate the outcome model and mediator model
  model_Y=summary(glm(formula_Y,family=binomial(link="logit"),data=data))
  model_M=summary(lm(formula_M,data=data))

  beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
  gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled

  ## (3) covariance matrix of (beta,gamma)
  nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
  S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
  S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
  colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))
  ## (4) approximate method
  ### NIE
  NIE=beta[3]*gamma[2]
  V_NIE = gamma[2]^2 * cov_beta[3,3] + beta[3]^2 * cov_gamma[2,2]
  ### TE
  TE = beta[2]+beta[3]*gamma[2]
  lambda = matrix(c(0,beta[3],rep(0,length(covariateM)),0,1,gamma[2],rep(0,length(covariateY))),ncol=1)
  V_TE=as.vector(t(lambda) %*% S %*% lambda)
  ### PM
  lambda = matrix(c(0,beta[2]*beta[3]/(beta[2]+beta[3]*gamma[2])^2,rep(0,length(covariateM)),0,-beta[3]*gamma[2]/(beta[2]+beta[3]*gamma[2])^2,beta[2]*gamma[2]/(beta[2]+beta[3]*gamma[2])^2,rep(0,length(covariateY))),ncol=1)
  PM = beta[3]*gamma[2]/(beta[2]+beta[3]*gamma[2])
  V_PM=as.vector(t(lambda) %*% S %*% lambda)
  ## (5) exact method
  ### preliminary parameters
  sigma2=model_M$sigma
  var_sigma2=2*model_M$sigma^2*(1/model_M$df[2])
  theta=c(gamma,beta,sigma2)
  S=matrix(0,ncol=nbeta+ngamma+1,nrow=nbeta+ngamma+1)
  S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
  S[(nbeta+ngamma+1),(nbeta+ngamma+1)]=var_sigma2
  colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""),"sigma2_M")
  loc_gamma_0 = 1:2
  loc_gamma_c = 3:(3+length(covariateM)-1)
  loc_beta_0 = (ngamma+1):(ngamma+3)
  loc_beta_c  = (ngamma+4):(ngamma+4+length(covariateY)-1)
  if (is.null(covariateM)) {loc_gamma_c=NULL}
  if (is.null(covariateY)) {loc_beta_c=NULL}
  ### NIE
  x0s=x0
  NIE_nonrare_p = NIE_unbiased(theta=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                               x0s=x0s,x1=x1,cY=cY,cM=cM)
  lambda= mygrad(NIE_unbiased,x0=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                 x0s=x0s,x1=x1,cY=cY,cM=cM)
  V_NIE_nonrare = as.vector(t(lambda) %*% S %*% lambda)

  TE_nonrare_p = TE_unbiased(theta=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                             x0s=x0s,x1=x1,cY=cY,cM=cM)
  lambda= mygrad(TE_unbiased,x0=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                 x0s=x0s,x1=x1,cY=cY,cM=cM)
  V_TE_nonrare = as.vector(t(lambda) %*% S %*% lambda)

  PM_nonrare_p = PM_unbiased(theta=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                             x0s=x0s,x1=x1,cY=cY,cM=cM)
  lambda= mygrad(PM_unbiased,x0=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                 x0s=x0s,x1=x1,cY=cY,cM=cM)
  V_PM_nonrare = as.vector(t(lambda) %*% S %*% lambda)

  point_est = c(NIE,TE,PM,NIE_nonrare_p,TE_nonrare_p,PM_nonrare_p);
  names(point_est)=c("NIE","TE","PM","NIE_nonrare","TE_nonrare","PM_nonrare")
  var_est = c(V_NIE,V_TE,V_PM,V_NIE_nonrare,V_TE_nonrare,V_PM_nonrare);
  names(var_est)=c("NIE","TE","PM","NIE_nonrare","TE_nonrare","PM_nonrare")
  sd_est = sqrt(var_est)
  names(sd_est)=c("NIE","TE","PM","NIE_nonrare","TE_nonrare","PM_nonrare")
  ci_est = rbind(point_est-1.96*sd_est,point_est+1.96*sd_est)
  rownames(ci_est) = c("Lower boundary","Upper boundary")

  return(list(point_est=point_est,var_est=var_est,sd_est=sd_est,ci_est=ci_est))
}

#' Using bootstrap to calculate the confidence intervals in the scenario of binary outcome and continuous mediator.
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param cY conditional levels of covariateY
#' @param cM conditional levels of covariateM
#' @param R The number of replications.
#' @return The 95% confidence interval by the bootstrap approach
Mediate_binaY_contM_bootci=function(data,
                                    outcome="Y",
                                    mediator="M",
                                    exposure="X",
                                    covariateY=c("X1","X2"),
                                    covariateM=c("X1","X2"),
                                    x0=0,x1=1,cY=c(0,0),cM=c(0,0),R=1000) {
  HermiteCoefs=function (order) {
    x <- 1
    if (order > 0)
      for (n in 1:order) x <- c(0, 2 * x) - c(((0:(n - 1)) *
                                                 x)[-1L], 0, 0)
      return(x)
  }

  gauss.hermite=function (f, mu = 0, sd = 1, ..., order = 5) {
    stopifnot(is.function(f))
    stopifnot(length(mu) == 1)
    stopifnot(length(sd) == 1)
    Hn <- HermiteCoefs(order)
    Hn1 <- HermiteCoefs(order - 1)
    x <- sort(Re(polyroot(Hn)))
    Hn1x <- matrix(Hn1, nrow = 1) %*% t(outer(x, 0:(order - 1),
                                              "^"))
    w <- 2^(order - 1) * factorial(order) * sqrt(pi)/(order *
                                                        Hn1x)^2
    ww <- w/sqrt(pi)
    xx <- mu + sd * sqrt(2) * x
    ans <- 0
    for (i in seq_along(x)) ans <- ans + ww[i] * f(xx[i], ...)
    return(ans)
  }

  mygrad=function (f, x0,heps = 1e-5, ...) {
    if (!is.numeric(x0))
      stop("Argument 'x0' must be a numeric value.")
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    p =length(f(x0))
    n <- length(x0)
    hh <- rep(0, n)
    gr <- matrix(0,nrow=n,ncol=p)
    for (i in 1:n) {
      hh[i] <- heps
      gr[i,] <- (f(x0 + hh) - f(x0 - hh))/(2 * heps)
      hh[i] <- 0
    }
    return(gr)
  }

  NIE_unbiased = function(theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                          x0s=x0s,x1=x1,cY=cY,cM=cM) {
    gamma0 = theta[1]
    gamma1 = theta[2]
    gamma_c = theta[loc_gamma_c]
    beta0  = theta[loc_beta_0[1]]
    beta1  = theta[loc_beta_0[2]]
    beta2  = theta[loc_beta_0[3]]
    beta_c = theta[loc_beta_c]
    sigma2 = theta[length(theta)]

    if (is.null(loc_beta_c)) {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
    } else {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
    }

    if (is.null(loc_gamma_c)) {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1,sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
    } else {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
    }

    output = log((p11s)/(1-p11s)) - log((p10s)/(1-p10s))
    return(output)
  }

  #NIE_unbiased_D = mygrad(NIE_unbiased,x0=theta)

  TE_unbiased = function(theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                         x0s=x0s,x1=x1,cY=cY,cM=cM) {
    gamma0 = theta[1]
    gamma1 = theta[2]
    gamma_c = theta[loc_gamma_c]
    beta0  = theta[loc_beta_0[1]]
    beta1  = theta[loc_beta_0[2]]
    beta2  = theta[loc_beta_0[3]]
    beta_c = theta[loc_beta_c]
    sigma2 = theta[length(theta)]

    if (is.null(loc_beta_c)) {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x)/(1+exp(beta0+beta1*x0s+beta2*x))}
    } else {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY)))}
    }

    if (is.null(loc_gamma_c)) {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1,sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=20)
    } else {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
    }

    output = log((p11s)/(1-p11s)) - log((p00s)/(1-p00s))
    return(output)
  }
  #TE_unbiased_D=mygrad(TE_unbiased,x0=theta)

  PM_unbiased=function(theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                       x0s=x0s,x1=x1,cY=cY,cM=cM) {
    gamma0 = theta[1]
    gamma1 = theta[2]
    gamma_c = theta[loc_gamma_c]
    beta0  = theta[loc_beta_0[1]]
    beta1  = theta[loc_beta_0[2]]
    beta2  = theta[loc_beta_0[3]]
    beta_c = theta[loc_beta_c]
    sigma2 = theta[length(theta)]

    if (is.null(loc_beta_c)) {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x)/(1+exp(beta0+beta1*x1+beta2*x))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x)/(1+exp(beta0+beta1*x0s+beta2*x))}
    } else {
      f11 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f10 = function(x) {exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x1+beta2*x+sum(beta_c*cY)))}
      f00 = function(x) {exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY))/(1+exp(beta0+beta1*x0s+beta2*x+sum(beta_c*cY)))}
    }

    if (is.null(loc_gamma_c)) {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1,sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s,sd=sqrt(sigma2),order=40)
    } else {
      p11s= gauss.hermite(f=f11,mu=gamma0+gamma1*x1+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p10s= gauss.hermite(f=f10,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
      p00s= gauss.hermite(f=f00,mu=gamma0+gamma1*x0s+sum(gamma_c*cM),sd=sqrt(sigma2),order=40)
    }

    te  = log((p11s)/(1-p11s)) - log((p00s)/(1-p00s))
    nie = log((p11s)/(1-p11s)) - log((p10s)/(1-p10s))
    return(nie/te)
  }

  data = as.data.frame(data)
  get_par_boot=function(data=data,indices) {
    data=data[indices,]
    ## (1) formula for outcome model and mediator model
    if (is.null(covariateY)) {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
    } else {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
    }
    if (is.null(covariateM)) {
      formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
    } else {
      formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
    }
    ## (2) estimate the outcome model and mediator model
    model_Y=summary(glm(formula_Y,family=binomial(link="logit"),data=data))
    model_M=summary(lm(formula_M,data=data))

    beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
    gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled

    ## (3) covariance matrix of (beta,gamma)
    nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
    S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
    S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
    colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))
    ## (4) approximate method
    ### NIE
    NIE=beta[3]*gamma[2]
    V_NIE = gamma[2]^2 * cov_beta[3,3] + beta[3]^2 * cov_gamma[2,2]
    ### TE
    TE = beta[2]+beta[3]*gamma[2]
    lambda = matrix(c(0,beta[3],rep(0,length(covariateM)),0,1,gamma[2],rep(0,length(covariateY))),ncol=1)
    V_TE=as.vector(t(lambda) %*% S %*% lambda)
    ### PM
    lambda = matrix(c(0,beta[2]*beta[3]/(beta[2]+beta[3]*gamma[2])^2,rep(0,length(covariateM)),0,-beta[3]*gamma[2]/(beta[2]+beta[3]*gamma[2])^2,beta[2]*gamma[2]/(beta[2]+beta[3]*gamma[2])^2,rep(0,length(covariateY))),ncol=1)
    PM = beta[3]*gamma[2]/(beta[2]+beta[3]*gamma[2])
    V_PM=as.vector(t(lambda) %*% S %*% lambda)
    ## (5) exact method
    ### preliminary parameters
    sigma2=model_M$sigma
    var_sigma2=2*model_M$sigma^2*(1/model_M$df[2])
    theta=c(gamma,beta,sigma2)
    S=matrix(0,ncol=nbeta+ngamma+1,nrow=nbeta+ngamma+1)
    S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
    S[(nbeta+ngamma+1),(nbeta+ngamma+1)]=var_sigma2
    colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""),"sigma2_M")
    loc_gamma_0 = 1:2
    loc_gamma_c = 3:(3+length(covariateM)-1)
    loc_beta_0 = (ngamma+1):(ngamma+3)
    loc_beta_c  = (ngamma+4):(ngamma+4+length(covariateY)-1)
    if (is.null(covariateM)) {loc_gamma_c=NULL}
    if (is.null(covariateY)) {loc_beta_c=NULL}
    ### NIE
    ### NIE
    x0s=x0
    NIE_nonrare_p = NIE_unbiased(theta=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                                 x0s=x0s,x1=x1,cY=cY,cM=cM)
    lambda= mygrad(NIE_unbiased,x0=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                   x0s=x0s,x1=x1,cY=cY,cM=cM)
    V_NIE_nonrare = as.vector(t(lambda) %*% S %*% lambda)

    TE_nonrare_p = TE_unbiased(theta=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                               x0s=x0s,x1=x1,cY=cY,cM=cM)
    lambda= mygrad(TE_unbiased,x0=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                   x0s=x0s,x1=x1,cY=cY,cM=cM)
    V_TE_nonrare = as.vector(t(lambda) %*% S %*% lambda)

    PM_nonrare_p = PM_unbiased(theta=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                               x0s=x0s,x1=x1,cY=cY,cM=cM)
    lambda= mygrad(PM_unbiased,x0=theta,loc_gamma_0=loc_gamma_0,loc_gamma_c=loc_gamma_c,loc_beta_0=loc_beta_0,loc_beta_c=loc_beta_c,
                   x0s=x0s,x1=x1,cY=cY,cM=cM)
    V_PM_nonrare = as.vector(t(lambda) %*% S %*% lambda)

    point_est = c(NIE,TE,PM,NIE_nonrare_p,TE_nonrare_p,PM_nonrare_p);
    names(point_est)=c("NIE","TE","PM","NIE_nonrare","TE_nonrare","PM_nonrare")

    return(point_est)
  }

  boot.par=boot::boot(data=data, statistic=get_par_boot, R=R)
  boot.parciNIE <- boot::boot.ci(boot.par, index=1, type=c("perc"))
  boot.parciTE <- boot::boot.ci(boot.par, index=2, type=c("perc"))
  boot.parciPM <- boot::boot.ci(boot.par, index=3, type=c("perc"))
  boot.parciNIE_nonrare <- boot::boot.ci(boot.par, index=4, type=c("perc"))
  boot.parciTE_nonrare <- boot::boot.ci(boot.par, index=5, type=c("perc"))
  boot.parciPM_nonrare <- boot::boot.ci(boot.par, index=6, type=c("perc"))


  ci_est_prec = c(boot.parciNIE$percent[4:5],
                  boot.parciTE$percent[4:5],
                  boot.parciPM$percent[4:5],
                  boot.parciNIE_nonrare$percent[4:5],
                  boot.parciTE_nonrare$percent[4:5],
                  boot.parciPM_nonrare$percent[4:5])
  names(ci_est_prec)=c(paste(rep("CI_",6),rep(c("NIE","TE","PM"),each=2),rep(c("_Low","_High"),times=3),sep=""),
                       paste(rep("CI_",6),rep(c("NIE_nonrare","TE_nonrare","PM_nonrare"),each=2),rep(c("_Low","_High"),times=3),sep=""))
  return(ci_est_prec)
}

####################################################################
#
# Functions for Case 4
#
####################################################################

#' Mediation analysis function for binary outcome and binary mediator
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param cY conditional levels of covariateY
#' @param cM conditional levels of covariateM
#' @return A list containing NIE, NDE and MP point and interval estimates.
mediate_binaY_binaM=function(data,
                             outcome="Y",
                             mediator="M",
                             exposure="X",
                             covariateY=c("X1","X2"),
                             covariateM=c("X1","X2"),
                             x0=0,x1=1,cY=c(0,0),cM=c(0,0)) {


  #data=dat1;outcome=Outcome;mediator=Mediator;covariateY=covariateM=CovariateY;exposure=Exposure

  data = as.data.frame(data)
  ## (1) formula for outcome model and mediator model
  if (is.null(covariateY)) {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
  } else {
    formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
  }
  if (is.null(covariateM)) {
    formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
  } else {
    formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
  }
  ## (2) estimate the outcome model and mediator model
  model_Y=summary(glm(formula_Y,family=binomial(link="logit"),data=data))
  model_M=summary(glm(formula_M,family=binomial(link="logit"),data=data))

  beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
  gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled

  ## (3) covariance matrix of (beta,gamma)
  nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
  S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
  S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
  colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

  ## (4) approximate method
  ##### NIE TE and PM expressions
  if (is.null(covariateY)==0) {
    names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
    beta_c=paste("beta_",covariateY,sep="")
  }
  if (is.null(covariateM)==0) {
    names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
    gamma_c=paste("gamma_",covariateM,sep="")
  }
  .A = paste("exp(gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .B= paste("exp(beta2+gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .C=paste("exp(gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .D= paste("exp(beta2+gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")

  NIEa_fun = function() {
    output = paste("log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
    return(output)
  }
  variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
  NIEa_D=deriv(parse(text=NIEa_fun()),variable)
  gamma0=gamma[1];gamma1=gamma[2];
  if(is.null(covariateM)==0) {
    for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
  }
  beta0=beta[1];beta1=beta[2];beta2=beta[3]
  if(is.null(covariateY)==0) {
    for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
  }

  TEa_fun = function() {
    output = paste("beta1*",(x1-x0),"+","log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
    return(output)
  }
  TEa_D=deriv(parse(text=TEa_fun()),variable)

  PMa_fun = function() {
    .UP = paste("log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
    .BOT = paste("beta1*",(x1-x0),"+","log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
    output=paste("(",.UP,")/(",.BOT,")")
    return(output)
  }
  PMa_D=deriv(parse(text=PMa_fun()),variable)


  NIEa_D = eval(NIEa_D)
  NIEa_p = NIEa_D[1]
  lambda= t(attr(NIEa_D,"gradient"))
  V_NIEa = as.vector(t(lambda) %*% S %*% lambda)

  TEa_D = eval(TEa_D)
  TEa_p = TEa_D[1]
  lambda= t(attr(TEa_D,"gradient"))
  V_TEa = as.vector(t(lambda) %*% S %*% lambda)

  PMa_D = eval(PMa_D)
  PMa_p = PMa_D[1]
  lambda= t(attr(PMa_D,"gradient"))
  V_PMa = as.vector(t(lambda) %*% S %*% lambda)

  ## (5) exact method
  ##### NIE TE and PM expressions
  .A = paste("exp(gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .B = paste("exp(beta2+gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .C =paste("exp(gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .D = paste("exp(beta2+gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
    paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
    ,")")
  .E = paste("exp(beta0+beta1*",x1,if(is.null(covariateY)==0) {
    paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
    ,")")
  .F = paste("exp(beta0+beta2+beta1*",x1,if(is.null(covariateY)==0) {
    paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
    ,")")
  .G = paste("exp(beta0+beta1*",x0,if(is.null(covariateY)==0) {
    paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
    ,")")
  .H = paste("exp(beta0+beta2+beta1*",x0,if(is.null(covariateY)==0) {
    paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
    ,")")

  NIE_fun = function() {
    .A1 = paste("(1+",.A,"+",.E,"*",.A,"+",.F,")")
    .A2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
    .B1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
    .B2 = paste("(1+", .F,"+",.D,"*(1+",.E,"))")
    output = paste("log(",.A1,"/",.A2,")+","log(",.B1,"/",.B2,")")
    #output = paste("log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
    return(output)
  }
  NIE_D=deriv(parse(text=NIE_fun()),variable)

  TE_fun = function() {
    .C1 = paste("(1+",.A,"+",.G,"*",.A,"+",.H,")")
    .C2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
    .D1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
    .D2 = paste("(1+", .H,"+",.D,"*(1+",.G,"))")
    output = paste("beta1*",(x1-x0),"+log(",.C1,"/",.C2,")+","log(",.D1,"/",.D2,")")
    return(output)
  }
  TE_D=deriv(parse(text=TE_fun()),variable)


  PM_fun = function() {
    .A1 = paste("(1+",.A,"+",.E,"*",.A,"+",.F,")")
    .A2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
    .B1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
    .B2 = paste("(1+", .F,"+",.D,"*(1+",.E,"))")
    .C1 = paste("(1+",.A,"+",.G,"*",.A,"+",.H,")")
    .C2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
    .D1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
    .D2 = paste("(1+", .H,"+",.D,"*(1+",.G,"))")
    output1 = paste("(log(",.A1,"/",.A2,")+","log(",.B1,"/",.B2,"))")
    output2 = paste("(beta1*",(x1-x0),"+log(",.C1,"/",.C2,")+","log(",.D1,"/",.D2,"))")
    return(paste(output1,"/",output2))
  }
  PM_D=deriv(parse(text=PM_fun()),variable)

  NIE_D = eval(NIE_D)
  NIE_p = NIE_D[1]
  lambda= t(attr(NIE_D,"gradient"))
  V_NIE = as.vector(t(lambda) %*% S %*% lambda)

  TE_D = eval(TE_D)
  TE_p = TE_D[1]
  lambda= t(attr(TE_D,"gradient"))
  V_TE = as.vector(t(lambda) %*% S %*% lambda)

  PM_D = eval(PM_D)
  PM_p = PM_D[1]
  lambda= t(attr(PM_D,"gradient"))
  V_PM = as.vector(t(lambda) %*% S %*% lambda)

  point_est = c(NIEa_p,TEa_p,PMa_p,NIE_p,TE_p,PM_p);
  names(point_est)=c("NIEa","TEa","PMa","NIE","TE","PM")
  var_est = c(V_NIEa,V_TEa,V_PMa,V_NIE,V_TE,V_PM);
  names(var_est)=c("NIEa","TEa","PMa","NIE","TE","PM")
  sd_est = sqrt(var_est)
  names(sd_est)=c("NIEa","TEa","PMa","NIE","TE","PM")
  ci_est = rbind(point_est-1.96*sd_est,point_est+1.96*sd_est)
  rownames(ci_est) = c("Lower boundary","Upper boundary")

  return(list(point_est=point_est,var_est=var_est,sd_est=sd_est,ci_est=ci_est))
}

#' Using bootstrap to calculate the confidence intervals in the scenario of binary outcome and binary mediator.
#' @param data A dataset.
#' @param outcome The outcome variable.
#' @param mediator The mediator variable.
#' @param exposure The exposure variable.
#' @param covariateY A vector of confounders in the outcome regression.
#' @param covariateM A vector of confounders in the mediator regression.
#' @param x0 The baseline exposure level.
#' @param x1 The new exposure level.
#' @param cY conditional levels of covariateY
#' @param cM conditional levels of covariateM
#' @param R The number of replications.
#' @return The 95% confidence interval by the bootstrap approach
Mediate_binaY_binaM_bootci=function(data,
                                    outcome="Y",
                                    mediator="M",
                                    exposure="X",
                                    covariateY=c("X1","X2"),
                                    covariateM=c("X1","X2"),
                                    x0=0,x1=1,cY=c(0,0),cM=c(0,0),R=1000) {
  #data=mydata;
  #outcome="Y";
  #mediator="M";
  #exposure="X";
  #covariateY=c("X1","X2","X3","X4","X5","X6","X7","X8");
  #covariateM=c("X1","X2","X3","X4","X5","X6","X7");
  #x0=0;x1=5;cY=rep(0,8);cM=rep(0,7)

  data = as.data.frame(data)
  get_par_boot=function(data=data,indices) {
    data=data[indices,]
    ## (1) formula for outcome model and mediator model
    if (is.null(covariateY)) {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,sep=""))
    } else {
      formula_Y=as.formula(paste(outcome,"~",exposure,"+",mediator,"+",paste(covariateY,collapse="+"),sep=""))
    }
    if (is.null(covariateM)) {
      formula_M=as.formula(paste(mediator,"~",exposure,sep=""))
    } else {
      formula_M=as.formula(paste(mediator,"~",exposure,"+",paste(covariateM,collapse="+"),sep=""))
    }
    ## (2) estimate the outcome model and mediator model
    model_Y=summary(glm(formula_Y,family=binomial(link="logit"),data=data))
    model_M=summary(glm(formula_M,family=binomial(link="logit"),data=data))

    beta=model_Y$coef[,1];cov_beta=model_Y$cov.unscaled
    gamma=model_M$coef[,1];cov_gamma=model_M$cov.unscaled
    ## (3) covariance matrix of (beta,gamma)
    nbeta=dim(cov_beta)[1];ngamma=dim(cov_gamma)[1]
    S=matrix(0,ncol=nbeta+ngamma,nrow=nbeta+ngamma)
    S[1:ngamma,1:ngamma]=cov_gamma; S[(ngamma+1):(nbeta+ngamma),(ngamma+1):(nbeta+ngamma)]=cov_beta
    colnames(S)=rownames(S)=c(paste(names(gamma),"_gamma",sep=""),paste(names(beta),"_beta",sep=""))

    ## (4) approximate method
    ##### NIE TE and PM expressions
    if (is.null(covariateY)==0) {
      names(cY) = paste(names(beta),"_betafix",sep="")[-c(1:3)]
      beta_c=paste("beta_",covariateY,sep="")
    }
    if (is.null(covariateM)==0) {
      names(cM) = paste(names(gamma),"_gammafix",sep="")[-c(1:2)]
      gamma_c=paste("gamma_",covariateM,sep="")
    }
    .A = paste("exp(gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .B= paste("exp(beta2+gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .C=paste("exp(gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .D= paste("exp(beta2+gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")

    NIEa_fun = function() {
      output = paste("log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
      return(output)
    }
    variable=c("gamma0","gamma1",if(is.null(covariateM)==0) {gamma_c},"beta0","beta1","beta2",if(is.null(covariateY)==0) {beta_c})
    gamma0=gamma[1];gamma1=gamma[2];
    if(is.null(covariateM)==0) {
      for (i in (1:length(covariateM))) {assign(gamma_c[i],gamma[2+i])}
    }
    beta0=beta[1];beta1=beta[2];beta2=beta[3]
    if(is.null(covariateY)==0) {
      for (i in (1:length(covariateY))) {assign(beta_c[i],beta[3+i])}
    }

    TEa_fun = function() {
      output = paste("beta1*",(x1-x0),"+","log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
      return(output)
    }


    PMa_fun = function() {
      .UP = paste("log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
      .BOT = paste("beta1*",(x1-x0),"+","log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
      output=paste("(",.UP,")/(",.BOT,")")
      return(output)
    }



    NIEa_p = eval(parse(text=NIEa_fun()))
    TEa_p = eval(parse(text=TEa_fun()))
    PMa_p = eval(parse(text=PMa_fun()))
    ## (5) exact method
    ##### NIE TE and PM expressions
    .A = paste("exp(gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .B = paste("exp(beta2+gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .C =paste("exp(gamma0+gamma1*",x1,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .D = paste("exp(beta2+gamma0+gamma1*",x0,if(is.null(covariateM)==0) {
      paste("+",paste(paste(gamma_c,"*",cM),collapse = "+"))}
      ,")")
    .E = paste("exp(beta0+beta1*",x1,if(is.null(covariateY)==0) {
      paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
      ,")")
    .F = paste("exp(beta0+beta2+beta1*",x1,if(is.null(covariateY)==0) {
      paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
      ,")")
    .G = paste("exp(beta0+beta1*",x0,if(is.null(covariateY)==0) {
      paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
      ,")")
    .H = paste("exp(beta0+beta2+beta1*",x0,if(is.null(covariateY)==0) {
      paste("+",paste(paste(beta_c,"*",cY),collapse = "+"))}
      ,")")

    NIE_fun = function() {
      .A1 = paste("(1+",.A,"+",.E,"*",.A,"+",.F,")")
      .A2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
      .B1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
      .B2 = paste("(1+", .F,"+",.D,"*(1+",.E,"))")
      output = paste("log(",.A1,"/",.A2,")+","log(",.B1,"/",.B2,")")
      #output = paste("log(","(1+",.A,")*","(1+",.B,")/((","1+",.C,")*(","1+",.D,"))",")")
      return(output)
    }

    TE_fun = function() {
      .C1 = paste("(1+",.A,"+",.G,"*",.A,"+",.H,")")
      .C2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
      .D1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
      .D2 = paste("(1+", .H,"+",.D,"*(1+",.G,"))")
      output = paste("beta1*",(x1-x0),"+log(",.C1,"/",.C2,")+","log(",.D1,"/",.D2,")")
      return(output)
    }


    PM_fun = function() {
      .A1 = paste("(1+",.A,"+",.E,"*",.A,"+",.F,")")
      .A2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
      .B1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
      .B2 = paste("(1+", .F,"+",.D,"*(1+",.E,"))")
      .C1 = paste("(1+",.A,"+",.G,"*",.A,"+",.H,")")
      .C2 = paste("(1+",.C,"+",.E,"*",.C,"+",.F,")")
      .D1 = paste("(1+", .F,"+",.B,"*(1+",.E,"))")
      .D2 = paste("(1+", .H,"+",.D,"*(1+",.G,"))")
      output1 = paste("(log(",.A1,"/",.A2,")+","log(",.B1,"/",.B2,"))")
      output2 = paste("(beta1*",(x1-x0),"+log(",.C1,"/",.C2,")+","log(",.D1,"/",.D2,"))")
      return(paste(output1,"/",output2))
    }

    NIE_p = eval(parse(text=NIE_fun()))
    TE_p = eval(parse(text=TE_fun()))
    PM_p = eval(parse(text=PM_fun()))


    point_est = c(NIEa_p,TEa_p,PMa_p,NIE_p,TE_p,PM_p);
    names(point_est)=c("NIEa","TEa","PMa","NIE","TE","PM")

    return(point_est)
  }

  boot.par=boot::boot(data=data, statistic=get_par_boot, R=R)
  boot.parciNIEa <- boot::boot.ci(boot.par, index=1, type=c("perc"))
  boot.parciTEa <- boot::boot.ci(boot.par, index=2, type=c("perc"))
  boot.parciPMa <- boot::boot.ci(boot.par, index=3, type=c("perc"))
  boot.parciNIE<- boot::boot.ci(boot.par, index=4, type=c("perc"))
  boot.parciTE <- boot::boot.ci(boot.par, index=5, type=c("perc"))
  boot.parciPM <- boot::boot.ci(boot.par, index=6, type=c("perc"))


  ci_est_prec = c(boot.parciNIEa$percent[4:5],
                  boot.parciTEa$percent[4:5],
                  boot.parciPMa$percent[4:5],
                  boot.parciNIE$percent[4:5],
                  boot.parciTE$percent[4:5],
                  boot.parciPM$percent[4:5])
  names(ci_est_prec)=c(paste(rep("CI_",6),rep(c("NIEa","TEa","PMa"),each=2),rep(c("_Low","_High"),times=3),sep=""),
                       paste(rep("CI_",6),rep(c("NIE","TE","PM"),each=2),rep(c("_Low","_High"),times=3),sep=""))
  return(ci_est_prec)
}

###########################################################
#
# Main Function
#
###########################################################



mediate=function(data,
                 outcome="Y1",
                 mediator="Mc",
                 exposure="X",
                 binary.outcome=0,
                 binary.mediator=0,
                 covariate.outcome=c("C1","C2"),
                 covariate.mediator=c("C1","C2"),
                 x0=0,
                 x1=1,
                 c.outcome=c(0,0),
                 c.mediator=c(0,0),
                 boot=0,
                 R=2000) {


  data=as.data.frame(data)
  covariateY=covariate.outcome
  covariateM=covariate.mediator
  cY<-c.outcome;cM<-c.mediator

  if (binary.outcome==0 & binary.mediator==0) {
    delta_res=mediate_contY_contM(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                  covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1)
    if (boot==1) {
      boot_res=Mediate_contY_contM_bootci(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                          covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,R=R)
    }
  }

  if (binary.outcome==0 & binary.mediator==1) {
    delta_res=mediate_contY_binaM(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                  covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,cY=cY,cM=cM)
    if (boot==1) {
      boot_res=Mediate_contY_binaM_bootci(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                          covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,R=R,cY=cY,cM=cM)
    }
  }

  if (binary.outcome==1 & binary.mediator==0) {
    delta_res=mediate_binaY_contM(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                  covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,cY=cY,cM=cM)
    if (boot==1) {
      boot_res=Mediate_binaY_contM_bootci(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                          covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,R=R,cY=cY,cM=cM)
    }
  }

  if (binary.outcome==1 & binary.mediator==1) {
    delta_res=mediate_binaY_binaM(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                  covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,cY=cY,cM=cM)
    if (boot==1) {
      boot_res=Mediate_binaY_binaM_bootci(data=data,outcome=outcome,mediator=mediator,exposure=exposure,
                                          covariateY=covariateY,covariateM=covariateM,x0=x0,x1=x1,R=R,cY=cY,cM=cM)
    }
  }

  if (binary.outcome==1) {
    point=delta_res[[1]]
    serror=delta_res[[3]]
    ci_delta = delta_res[[4]]
    res=as.data.frame(rbind(point,serror,ci_delta))
    colnames(res)=c("Approximate NIE","Approximate TE","Approximate MP",
                    "Exact NIE","Exact TE","Exact MP")
    rownames(res)=c("point estimate","S.E by Delta Method","CI Lower by Delta Method",
                    "CI Upper by Delta Method")
    if (boot==1) {
      ci_boot=as.data.frame(rbind(boot_res[c(1,3,5,7,9,11)],boot_res[c(2,4,6,8,10,12)]))
      colnames(ci_boot)=c("Approximate NIE","Approximate TE","Approximate MP",
                          "Exact NIE","Exact TE","Exact MP")
      rownames(ci_boot)=c("CI Lower by Bootstrap Method",
                          "CI Upper by Bootstrap Method")
      res=rbind(res,ci_boot)
    }
  }

  if (binary.outcome==0) {
    point=delta_res[[1]]
    serror=delta_res[[3]]
    ci_delta = delta_res[[4]]
    res=as.data.frame(rbind(point,serror,ci_delta))
    colnames(res)=c("NIE","TE","MP")
    rownames(res)=c("point estimate","S.E by Delta Method","CI Lower by Delta Method",
                    "CI Upper by Delta Method")
    if (boot==1) {
      ci_boot=as.data.frame(rbind(boot_res[c(1,3,5)],boot_res[c(2,4,6)]))
      colnames(ci_boot)=c("NIE","TE","MP")
      rownames(ci_boot)=c("CI Lower by Bootstrap Method",
                          "CI Upper by Bootstrap Method")
      res=rbind(res,ci_boot)
    }
  }
  res=list(res=res,class="mediate")
  attr(res, "class") <- "mediate"
  #class(res)<-"mediateP"
  print.mediate(res)

  invisible(res)
}


#' Print mediation results
#' @param x An object of class `mediate`.
#' @param ... other parameters in `print`.
print.mediate=function(x, ...) {
  res=format(x$res, digits=3)
  isboot=ifelse(dim(res)[1]==4,0,1)
  iscontY=ifelse(dim(res)[2]==3,1,0)
  if (iscontY==1) {num.row=3} else {num.row=6}
  if (isboot==1) {num.col=3} else {num.col=2}

  out=as.data.frame(matrix(0,ncol=num.col,nrow=num.row))
  if (isboot==1) {
    colnames(out)=c("Point (S.E.)"," 95% CI by Delta Approach"," 95% CI by Bootstrap")
  } else {
    colnames(out)=c("Point (S.E.)"," 95% CI by Delta Approach")
  }
  if (iscontY==1) {
    rownames(out)=c("NIE:  ","TE:   ","MP:   ")
  } else {
    rownames(out)=c("NIE: Approximate  ","TE:  Approximate  ","MP:  Approximate  ",
                    "NIE:       Exact  ","TE:        Exact  ","MP:        Exact  ")
  }
  for (i in (1:num.row)) {
    out[i,1]=paste(res[1,i]," (",res[2,i],")",sep="")
    out[i,2]=paste("(",res[3,i],",",res[4,i],")",sep="")
    if (isboot==1) {
      out[i,3]=paste("(",res[5,i],",",res[6,i],")",sep="")
    }
  }
  cat(paste("Mediation Analysis Results\n"))
  print.data.frame(out)
}
