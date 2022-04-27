library("fdrtool")
library(isotone)
library(stats)
library(data.table)
library(xtable)
library(mvtnorm)
library("msm")
library(dplyr)

source("main functions.R")
load("Ainv.rda")

# define true regression functions
#define f1
f1 <- function(x){
  return (x^2 + x/5)
}

#define f2
f2 <- function(x){
  xx <- exp(4 * (x - 0.5))
  return (xx / (1 + xx))
}
##############################################

# set data generation parameters
set.seed(100)
n = 500; f=f1; t0=0.5; sigma=0.1;
d0 <- f(t0) # d0=0.35 for f=f1, t0=0.5

# generate data X and Y; X is on (0,1)
x <- seq(1/n, 1, 1/n)
y <- f(x) + rnorm(n, 0, sigma) 

# fit projection-posterior function
fit1stage.rq = reg.quantile.projection(x, y, d0 = d0, credibility.level=0.95, nPostSamp=1000)

# uncorrected credible interval
fit1stage.rq[[1]]

# recalibrated credible interval 
fit1stage.rq[[2]]

# true value of the parameter
t0

###################################################################

###################  2-stage example run ###########################

# set data generation parameters
set.seed(100)
n = 500; f=f1; t0=0.5; sigma=0.1;
d0 <- f(t0) # d0=0.35 for f=f1, t0=0.5

n1 = floor(n/2)

# generate first-stage data X and Y; X is on (0,1)
x <- seq(1/n1, 1, 1/n1)
y <- f(x) + rnorm(n1, 0, sigma) 

# fit projection-posterior function
fit1.rq = reg.quantile.projection(x, y, d0 = d0, credibility.level=0.99, nPostSamp=1000)
interval.1stage = fit1.rq[[1]]

# generate data on the obtained sampling interval
n2 = n-n1
firstStageUpper = min(interval.1stage[2],1)
firstStageLower = max(interval.1stage[1],0)
mu.tilde = (firstStageUpper+firstStageLower)/2

x <- runif(n2, firstStageLower, firstStageUpper)
y <- f(x) + rnorm(n2, 0, sigma)

fit2.rq = reg.quantile.2stage(x, y, d0 = d0, mu.tilde = mu.tilde, 
                              sigmaHatSq = fit1.rq[[3]], nPostSamp=1000)

# 95% 2-stage credible interval for mu
c(quantile(fit2.rq, 0.025), quantile(fit2.rq, 0.975))

# 2-stage point-estimator of mu
mean(fit2.rq)

# true value of mu
t0