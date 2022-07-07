## ------------------------------------------------------------------------
set.seed(12032016)

## ------------------------------------------------------------------------
n <- 5000
x1 <- runif(n)
x2 <- runif(n)
ilogit <- function(x) 1/(1+exp(-x))
p <- ilogit(sin(2*pi*x1)*sin(2*pi*x2))
lambda <- exp(10*x1*(1-x1))
z <- rbinom(n,1,p)
y <- rpois(n,lambda*z)
d <- data.frame(x1=x1,x2=x2,y=y)

## ------------------------------------------------------------------------
library(mgcv)
library(zigam)

## ------------------------------------------------------------------------
fit <- zipgam(y ~ s(x1,x2),
              ~ s(x1,x2),
              data=d)

## ---- fig.height=9, fig.width=9------------------------------------------
plot(fit$fit.pi,select=1,scheme=2)

## ---- fig.height=9, fig.width=9------------------------------------------
plot(fit$fit.lambda,select=1,scheme=2)

## ------------------------------------------------------------------------
n <- 200
x1 <- runif(n)
x2 <- runif(n)
ilogit <- function(x) 1/(1+exp(-x))
p <- ilogit(1+sin(2*pi*x2))
lambda <- exp(10*x1*(1-x1))
z <- rbinom(n,1,p)
y <- rpois(n,lambda*z)
d <- data.frame(x1=x1,x2=x2,y=y)

## ------------------------------------------------------------------------
fit <- zipgam(y ~ s(x1),
              ~ s(x2),
              data=d)

## ------------------------------------------------------------------------
x1.pred <- seq(0,1,length.out=40)
d.pred <- data.frame(x1=x1.pred,x2=rep(0,40))

## ------------------------------------------------------------------------
d.boot <- d
boot <- sapply(simulate(fit,20),function(s) {
  ## Update response
  d.boot[,"y"] <- s
  ## Refit model and generate predictions
  fit.boot <- update(fit,data=d.boot)
  predict(fit.boot,d.pred,type="response")
})

## ---- fig.height=9, fig.width=9------------------------------------------
pr.ci <- t(apply(boot,1,quantile,prob=c(0.025,0.5,0.975)))
matplot(x1.pred,pr.ci,type="l",lty=1,
        col=c("dodgerblue1","dodgerblue3","dodgerblue1"))

