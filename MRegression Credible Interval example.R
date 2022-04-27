# title: "Monotone regression: point-wise credible interval"
# author: "Moumita Chakraborty"

library("msm")
library("fdrtool")
library(dplyr)
library(isotone)

source("main functions.R")
#load the A-inverse function for the recalibration values
load("Ainv.rda")


# true regression functions 
#define f1
f1 <- function(x){
  return (x^2 + x/5)
}

#define f2
f2 <- function(x){
  xx <- exp(4 * (x - 0.5))
  return (xx / (1 + xx))
}
####################################################################

# Following is an example to construct a credible interval with the 
# specified $100(1-\alpha)\%$ projection-posterior credibility, and a 
# corrected credible interval with $100(1-\alpha)\%$ asymptotic coverage. 

set.seed(1010)
n <- 100
sigma <- 0.1
x0 <- 0.5
f <- f1

x <- runif(n)
y <- f(x) + rnorm(n, 0, sigma)
cred0 <- projection.credible.interval(x, y, x0 = x0, credibility.level=0.95,
                                      plot.function=FALSE, nPostSamp=1000)

# 100(1-\alpha)\%$ projection-posterior credible interval:
cred0[[1]] # unadjusted credible interval

# corrected projection-posterior credible interval with asymptotic coverage 
#  $100(1-\alpha)\%$ :
cred0[[2]] # adjusted credible interval

# the true value of $f(x_0)$ :
f(x0)
