# title: "Monotone regression testing"
# author: "Moumita Chakraborty"

library(stats)
library(isotone)

source("main functions.R")

#define the true functions

#monotone functions
f1 <- function(x){
  return (x*0)
}
f2 <- function(x){ifelse(x<=0.6, 0, 0.2)}
f3 <- function(x){ifelse(x<=0.4, x/0.4, 1)}
f4 <- function(x){
  (pbeta(x, 1, 1)+pbeta(x, 200, 80)+pbeta(x, 80, 200))/3
}
f5 <- function(x){x^7/(10*(x^7)+0.7^7)}


#non-monotone functions
nm1 <- function(x){ifelse(x<0.5, 4*(x-0.5)^3 + 0.1*(x-0.5) - 0.25*exp(-250*(x-0.25)^2),
                          0.1*(x-0.5) - 0.25*exp(-250*(x-0.25)^2))}
nm2 <- function(x){ -0.1*x}
nm3 <- function(x){-0.1*exp(-50*(x-0.5)^2)}
nm4 <- function(x){0.1*cos(6*pi*x)}
nm5 <- function(x){0.2*x + nm3(x)}
nm6 <- function(x){ 0.2*x + nm4(x)}
nm7 <- function(x){x+ 0.415*exp(-50*(x)^2)}
nm8 <- function(x){x+1 - 0.45*exp(-50*(x-0.5)^2)}
f7 <- function(x){stepfun(c(0.45, 0.55), c(0, 0.2, 0))(x)}
f8 <- function(x){sin(x*3)}


# Example run on monotone functions

nRepeat <- 500  # number of replications
sigma=0.1

set.seed(100) # 1000,1010

fvec = c(f1)#, f2, f4)#, f5)
nvec <- c(50) #, 100, 200, 500)
for(l in 1:length(fvec)){
  f = fvec[[l]]
  print(paste("f=m", sep="", l))
  for(l.n in 1:length(nvec)){
    n = nvec[l.n]
    #cut.0 <- M0 / (n^(1/3))
    print(paste("n=", sep="", n))
    
    testVec = sal.vec = NULL
    for(i in 1:nRepeat){
      x = seq(0,1,length = n+1)[-(n+1)]
      y = f(x) + rnorm(n,sd = sigma)
      testVec <- c(testVec, testL1(x, y, nPostSamp=1000))
    }
    
    print(sum(testVec)/nRepeat)
  }
}


# Example run on non-monotone functions

set.seed(100)

fvec = c(nm2)#, nm3,nm4, nm6, nm7, nm8)
#fvec = c(f1, f2, f4, f5, f6, f9, f7, f8)
nvec <- c(200)#, 500, 600, 700)
for(l in 1:length(fvec)){
  f = fvec[[l]]
  print(paste("f=m", sep="", l))
  for(l.n in 1:length(nvec)){
    n = nvec[l.n]
    #cut.0 <- M0 / (n^(1/3))
    print(paste("n=", sep="", n))
    
    testVec = sal.vec = NULL
    for(i in 1:nRepeat){
      x = seq(0,1,length = n+1)[-(n+1)]
      y = f(x) + rnorm(n,sd = sigma)
      testVec <- c(testVec, testL1(x, y, nPostSamp=1000))
    }
    
    print(sum(testVec)/nRepeat)
  }
}


############## Calculating the time elapsed ################
set.seed(100)
sigma=0.1
n=500
f=f1
x = seq(0,1,length = n+1)[-(n+1)]
y = f(x) + rnorm(n,sd = sigma)

#L_1 projection method
start.time = proc.time()
testL1(x, y, nPostSamp=2500)

#time elapsed L_1 projection method
print(proc.time() - start.time)


################  Comparison of run time with Salomond's method ###########
library(dplyr)
library(MASS)
library(invgamma)
library(pbapply)
source("JBS code.R")

sal.run = function(f,n,sd = 0.1,verbose = F,autoM0 = F,prior,lambda,M0,mu){
  X = seq(0,1,length = n+1)[-(n+1)]
  y = f(X) + rnorm(n,sd = sd)
  return(test(X,y,prior,alpha = 100, beta = .5,mu = mu,lambda = lambda,Kmax= n-1,
              Nk = 2500,M0 = M0,verbose = verbose,autoM0 = autoM0) ) 
}

set.seed(100)

# calculating time
sigma=0.1
n=500
f=f1
x = seq(0,1,length = n+1)[-(n+1)]
y = f(x) + rnorm(n,sd = sigma)

start.time = proc.time()
sal.test(x, -y,prior = "geom",alpha = 100, beta = 0.5,
                                 lambda = .4,M0 = 1,mu = 1,autoM0 = T,
                                 verbose = F, Kmax= n-1,Nk = 2500)

#time elapsed in Salomond's method
print(proc.time() - start.time)

