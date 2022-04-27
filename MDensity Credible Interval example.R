
# title: "Monotone density point-wise credible interval"
# author: "Moumita Chakraborty"


# load relevant packages
library(fdrtool)
library(stats)
library(data.table)
library(MCMCpack) # for rdirichlet
library(hbmem)
library(dplyr)
library(isotone)

# load source functions
source("main functions.R")
# load the A-inverse function for the recalibration values
load("Ainv.rda")

# Example on (0,1)
# generate data

set.seed(101)
x0 <- 0.4
d00 <- dbeta(x0, 1, 3)
n <- 100;
X <- rbeta(n, 1,3)

# run function
cred0 <- density.credInt.bounded(X, x0=x0, credibility.level=0.95,
                                 nPostSamp = 1000)

# The $100(1-\alpha)\%$ projection-posterior credible interval is
cred0[[1]] # unadjusted credible interval

# The corrected projection-posterior credible interval with asymptotic coverage  $100(1-\alpha)\%$ is
cred0[[2]] # adjusted credible interval

# The $95\%$ confidence interval based on the sieve-MLE method is
cred0[[3]] # confidence interval using sieve-MLE

#The true value of $g(x_0)$ is
dbeta(x0, 1, 3)


# Example (unbounded)
# generate data
set.seed(101)
n=1000
X <- rexp(n)
x0 <- 1.5

# run function
credible.interval.unbounded(X, x0=x0, credibility.level=0.95, nPostSamp = 1000, plot.function = TRUE)

# The true value of $g(x_0)$ is
dexp(x0)

