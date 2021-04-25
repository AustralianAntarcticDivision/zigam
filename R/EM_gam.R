## version 0.1 respresents source file EM_gam_source_8-March-2016.R



#' Zero-inflated Poisson GAM
#'
#' Fit a zero-inflated Poisson Generalized Additive Model
#' using the EM Algorithm
#' @param lambda.formula formula for the count model
#' @param pi.formula formula for the binary model
#' @param data a data frame or list containing the model
#'  response variable and covariates required by the formula.
#' @param lambda initial lambda vector
#' @param pi initial pi vector
#' @param gamma.pi binary model gamma
#' @param gamma.lambda count model gamma
#' @param select include model selection penalty
#' @param method method for selecting smoothing parameters
#' @param min.em minimum number of EM iterations
#' @param max.em maximum number of EM iterations
#' @param tol tolerance (default=1.0E-2)
#' @export
zipgam <- function(lambda.formula,pi.formula,data,
                    lambda=NULL,pi=NULL,
                    gamma.pi=1,gamma.lambda=1,select=FALSE,
                    method="GCV.Cp",min.em=5,max.em=50,tol=1.0E-2) {
  ## Log density
  dzip.log <- function(x, lambda, pi) {
    logp <- log(pi) + dpois(x, lambda, log=TRUE)
    logp[x==0] <- log(exp(logp[x==0]) + (1-pi[x==0]))
    logp
  }

  cl <- match.call()
  ## Extract the response y
  mf <- model.frame(update(lambda.formula, .~1),data=data)
  y <- model.response(mf)
  N <- length(y)
  ## Set initial pi, lambda
  if(is.null(lambda)) lambda <- mean(y)
  if(is.null(pi)) pi <- mean(y>0)
  ## Response for pi component is the weights
  pi.formula <- update(pi.formula, w ~ .)
  environment(pi.formula) <- environment()
  environment(lambda.formula) <- environment()
  logL <- double(max.em)
  ## Evaluate initial weights
  w <- ifelse(y==0,pi*dpois(0,lambda)/(1-pi+pi*dpois(0,lambda)),1)
  ## Setup models for fitting
  G.pi <- suppressWarnings(gam(pi.formula,family=binomial(),
                               select=select,fit=FALSE,data=data))
  G.lambda <- suppressWarnings(gam(lambda.formula,weights=w,family=poisson(),
                                   select=select,fit=FALSE,data=data))
  for(k in 1:max.em) {
    ## Update models for current iteration
    G.pi$y <- ifelse(y==0,pi*dpois(0,lambda)/(1-pi+pi*dpois(0,lambda)),1)
    fit.pi <- suppressWarnings(gam(G=G.pi,gamma=gamma.pi,method=method))
    pi <- predict(fit.pi,type="response")
    G.lambda$w <- ifelse(y==0,pi*dpois(0,lambda)/(1-pi+pi*dpois(0,lambda)),1)
    fit.lambda <- suppressWarnings(gam(G=G.lambda,gamma=gamma.lambda,method=method))
    lambda <- predict(fit.lambda,type="response")
    ## Evaluate likelihood
    logL[k] <- sum(dzip.log(y,lambda,pi))
    #print(logL[1:k])
    if(k>min.em && abs(logL[k]-logL[k-1]) < tol) {
      logL <- logL[1:k]
      break
    }
  }
  ## Calculate degrees of freedom and aic
  df <- attr(logLik(fit.pi),"df")+attr(logLik(fit.lambda),"df")
  aic <- 2*(df-logL[length(logL)])
  ## Return results
  fit <- list(fit.lambda=fit.lambda,
              fit.pi=fit.pi,
              lambda=lambda,
              pi=pi,
              w=ifelse(y==0,pi*dpois(0,lambda)/(1-pi+pi*dpois(0,lambda)),1),
              aic=aic,
              logL=logL,
              call=cl)
  class(fit) <- "zipgam"
  fit
}

#' Simulate response from the fitted model
#'
#' Simulate response from a model of class zipgam.
#' @param object an object representing a fitted model of class zipgam.
#' @param nsim number of response vectors to simulate. Defaults to 1.
#' @param seed an object specifying if and how the random number
#' generator should be initialized ('seeded').
#' @param ... additional optional arguments.
#' @export
simulate.zipgam <- function(object,nsim=1,seed=NULL,...) {
  suppressWarnings(simulate(object$fit.lambda,nsim=nsim,seed=seed,...)*
                     simulate(object$fit.pi,nsim=nsim,seed=seed,...))
}

#' Predict from the fitted model
#'
#' Obtain predictions from a fitted model of class zipgam.
#' @param object an object representing a fitted model of class zipgam.
#' @param newdata data to predict from.
#' @param type the type of prediction required.
#' @param ... additional optional arguments.
#' @export
predict.zipgam <- function(object,newdata,type=c("response","link"),...) {
  type <- match.arg(type)
  switch(type,
         response={
           z <- predict(object$fit.pi,newdata,type="response")
           y <- predict(object$fit.lambda,newdata,type="response")
           z*y},
         pi.response=predict(object$fit.pi,newdata,type="response"),
         lambda.response=predict(object$fit.lambda,newdata,type="response"),
         pi.link=predict(object$fit.pi,newdata,type="link"),
         lambda.link=predict(object$fit.lambda,newdata,type="link"))
}

#' Zero-inflated Negative Binomial GAM
#'
#' Fit a zero-inflated Negative Binomial Generalized Additive
#' Model using the EM Algorithm
#' @param mu.formula formula for the count model
#' @param pi.formula formula for the binary model
#' @param data a data frame or list containing the model
#'  response variable and covariates required by the formula.
#' @param mu initial mu vector
#' @param pi intial pi vector
#' @param theta initial theta value
#' @param gamma.mu count model gamma
#' @param gamma.pi binary model gamma
#' @param select include model selection penalty
#' @param method method for selecting smoothing parameters
#' @param min.em minimum number of EM iterations
#' @param max.em maximum number of EM iterations
#' @param tol tolerance (default=1.0E-2)
#' @export
zinbgam <- function(mu.formula,pi.formula,data,
                     mu=NULL,pi=NULL,theta=1,
                     gamma.pi=1,gamma.mu=1,select=FALSE,
                     method="GCV.Cp",min.em=5,max.em=50,tol=1.0E-2) {
  ## Log density
  dzinb.log <- function(x,mu,pi,shape) {
    logp <- log(pi)+dnbinom(x,size=shape,mu=mu,log=T)
    logp[x==0] <- log(exp(logp[x==0])+(1-pi[x==0]))
    logp
  }

  cl <- match.call()
  ## Extract the response y
  mf <- model.frame(update(mu.formula,.~1),data=data)
  y <- model.response(mf)
  N <- length(y)
  ## Set initial pi, mu
  if(is.null(mu)) mu <- mean(y)
  if(is.null(pi)) pi <- mean(y>0)
  ## Response for pi component is the weights
  pi.formula <- update(pi.formula, w ~ .)
  environment(pi.formula) <- environment()
  environment(mu.formula) <- environment()
  logL <- double(max.em)
  ## Evaluate initial weights
  w <- ifelse(y==0,pi*dnbinom(0,size=theta,mu=mu)/(1-pi+pi*dnbinom(0,size=theta,mu=mu)),1)
  ## Setup models for fitting
  G.pi <- suppressWarnings(gam(pi.formula,family=binomial(),
                               select=select,fit=FALSE,data=data))
  G.mu <- suppressWarnings(gam(mu.formula,weights=w,family=nb(),
                               select=select,fit=FALSE,data=data))
  for(k in 1:max.em) {
    ## Update models for current iteration
    G.pi$y <- ifelse(y==0,pi*dnbinom(0,size=theta,mu=mu)/(1-pi+pi*dnbinom(0,size=theta,mu=mu)),1)
    fit.pi <- suppressWarnings(gam(G=G.pi,gamma=gamma.pi,method=method))
    pi <- predict(fit.pi,type="response")
    G.mu$w <- ifelse(y==0,pi*dnbinom(0,size=theta,mu=mu)/(1-pi+pi*dnbinom(0,size=theta,mu=mu)),1)
    fit.mu <- suppressWarnings(gam(G=G.mu,gamma=gamma.mu,method=method))
    mu <- predict(fit.mu,type="response")
    theta <- fit.mu$family$getTheta(TRUE)
    ## Evaluate likelihood
    logL[k] <- sum(dzinb.log(y,mu,pi,theta))
    if(k>min.em && abs(logL[k]-logL[k-1]) < tol) {
      logL <- logL[1:k]
      break
    }
  }
  ## Calculate degrees of freedom and aic
  df <- attr(logLik(fit.pi),"df")+attr(logLik(fit.mu),"df")
  aic <- 2*(df-logL[length(logL)])
  ## Return results
  fit <- list(fit.mu=fit.mu,
              fit.pi=fit.pi,
              mu=mu,
              pi=pi,
              w=ifelse(y==0,pi*dnbinom(0,size=theta,mu=mu)/(1-pi+pi*dnbinom(0,size=theta,mu=mu)),1),
              aic=aic,
              logL=logL,
              theta=theta,
              call=cl)
  class(fit) <- "zinbgam"
  fit
}


## no simulate method for negative binomial family
# Simulate response from the fitted model
#
# Simulate response from a model of class zinbgam.
# @param object an object representing a fitted model of class zinbgam.
# @param nsim number of response vectors to simulate. Defaults to 1.
# @param seed an object specifying if and how the random number
# generator should be initialized (‘seeded’).
# @param ... additional optional arguments.
# @importFrom stats simulate.lm
# @export
# simulate.zinbgam <- function(object,nsim=1,seed=NULL,...) {
#   suppressWarnings(simulate(object$fit.mu,nsim=nsim,seed=seed,...)*
#                      simulate(object$fit.pi,nsim=nsim,seed=seed,...))
# }

#' Predict from the fitted model
#'
#' Obtain predictions from a fitted model of class zinbgam.
#' @param object an object representing a fitted model of class zinbgam.
#' @param newdata data to predict from.
#' @param type the type of prediction required.
#' @param ... additional optional arguments.
#' @export
predict.zinbgam <- function(object,newdata,type=c("response","link"),...) {
  type <- match.arg(type)
  switch(type,
         response={
           z <- predict(object$fit.pi,newdata,type="response")
           y <- predict(object$fit.mu,newdata,type="response")
           z*y},
         pi.response=predict(object$fit.pi,newdata,type="response"),
         mu.response=predict(object$fit.mu,newdata,type="response"),
         pi.link=predict(object$fit.pi,newdata,type="link"),
         mu.link=predict(object$fit.mu,newdata,type="link"))
}
