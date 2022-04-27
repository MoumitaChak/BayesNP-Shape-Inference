
library("bnpmr")
library(stats)
library(isotone)

source("main functions.R")

#define the true functions
f1 <- function(x){
  return (x*0)
}
f2 <- function(x){x^2 + x/5}
f3 <- function(x){ifelse(x<=0.6, 0, 1)}
f5 <- function(x){x^7/(x^7+0.4^7)}
f4 <- function(x){
  (pbeta(x, 1, 1)+pbeta(x, 200, 80)+pbeta(x, 80, 200))/3
}
#fmu3 <- function(x){ifelse(x<=0.4, x/0.4, 1)}


# set data
set.seed(101)
grd <- seq(0, 1, 0.01)
n=100; sigma=0.1;
f=f2
x <- c(1:n)/(n+1)
y <- f(x) + rnorm(n, 0, sigma)

# estimate f using our method
iso.obj <- pr.posterior.samples(x, y, Jsize="> n power 1/3", nPostSamp = 1000)
iso.pred <- predict.iso(grd, iso.obj[[1]], iso.obj[[2]])


iso050 <- apply(iso.pred, 2, mean) # projection-posterior mean
# projection-posterior quantiles
iso005 <- apply(iso.pred, 2, quantile, prob = 0.025)
iso095 <- apply(iso.pred, 2, quantile, prob = 0.975)


# estimate f using isotonic regression
fit.nb <- gpava(x, y, solver=weighted.median)

# estimate f using bnpmr
res <- bnpmr(y, x, burnIn = 5000, niter = 20000)
aa <- pred.bnpmr(grd, res)

# L1-distanve of projection-posterior f from monotone class
dist.cal(iso050, grd, f)

# L1-distanve of bnpmr f from monotone class
#out050 <- apply(aa, 2, mean)
#dist.cal(out050, grd, f)

# L1-distanve of isotonic regression f from monotone class
#dist.cal(stepfun(c(1:n)/n, c(fit.nb$x[1], fit.nb$x))(grd), grd, f)

######### plot results  ##################
plot(x,y, col="lightblue", main="f=f2", pch=19)
lines(grd, iso005,lty=2, lwd=2)
lines(grd, iso050, lwd=2)
lines(grd, iso095, lty=2, lwd=2)
curve(f, add=TRUE, col=2, lwd=2)
legend(0, 1.1, legend=c("True f", "Projection-posterior estimate", "95% credible region"),
       col=c("red", "black", "black"), lty=c(1, 1, 2), cex=0.4)
############################################

##############  compare run times ###################

# `bnpmr' method
#start.time = proc.time()
#res <- bnpmr(y, x, burnIn = 5000, niter = 20000)
#print(proc.time() - start.time)

# our method
start.time = proc.time()
iso.obj <- pr.posterior.samples(x, y, Jsize="> n power 1/3", nPostSamp = 1000)
print(proc.time() - start.time)
######################################################

