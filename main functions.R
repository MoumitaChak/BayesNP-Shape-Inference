# isotone.f function: to evaluate the isotonic regression estimator
# the vector y should be sorted wrt x's
#return values: r gives the isotonized step-function, t is its value at t1
isotone.f <- function(y, w, t1){
  J <- length(w)
  xval <- c(1:J) / J
  
  index.w <- which(w==0)
  w0 <- w
  
  y <- y[w>0]
  xval <- xval[w>0]
  w <- w[w>0]
  
  fit <- gpava(xval, y, weights = w, solver = weighted.mean,
               ties = "primary", p = NA)
  
  r <- stepfun(xval, c(fit[[1]][1], fit[[1]]))
  return(r)
}

# isotone2 function: to evaluate the isotonic regression estimator
# the vector y should be sorted wrt x's
#return values: r gives the isotonized step-function, t is its value at t1
isotone2 <- function(y, w, t1){
  J <- length(w)
  xval <- c(1:J) / J
  
  index.w <- which(w==0)
  w0 <- w
  
  y <- y[w>0]
  xval <- xval[w>0]
  w <- w[w>0]
  
  fit <- gpava(xval, y, weights = w, solver = weighted.mean,
               ties = "primary", p = NA)
  
  r <- stepfun(xval, c(fit[[1]][1], fit[[1]]))
  t <- r(t1) #value of the isotonized step function at given t1
  return(list(r, t))
}

##########  extract Nj  ##########################
# obtain Nj, yBar and sigmahat
extract.nj.sigmahat <- function(x, y, J, zi, tausq){
  n <- length(x)
  ybar <- c(1:J); xbar <- c(1:J) / J; nj <- c(1:J);
  sigmaHatSq <- NULL;
  
  for(j in 1:J){
    xInIj <- x > (j-1)/J & x <= j/J #indices present in j-th interval.
    nj[j] <- sum(xInIj) #number of Xj's in j-th interval.
    
    if(nj[j] == 0){
      ybar[j] <- 0
      next # go to next j if ther's no data in this interval
    }
    else{
      yInIj <-  y[xInIj]
      ybar[j] <- mean(yInIj) # define ybar
      sigmaHatSqIj <- t(yInIj-zi[j])%*%solve(diag(nj[j])+
                                               matrix(tausq[j],nj[j],nj[j]))%*%(yInIj-zi[j])
      sigmaHatSq <- c(sigmaHatSq, sigmaHatSqIj)
    }
  }
  sigmaHatSq <- sum(sigmaHatSq)/n
  sigmaHat <- sqrt(sigmaHatSq)
  return(list(nj, ybar, sigmaHat))
}
##########################################################################

# function to obtain a projection-posterior credible interval for a given x_0
# for monotone regression function f in (0,1)

# The function `projection.credible.interval' takes $x, y, x_0$ and the credibility level 
# as inputs, and optionslly generates a plot of the data along with a credible region. 
# The input $x$ must be a vector in $[0,1]$.
projection.credible.interval <- function(x, y, x0 = 0.5, credibility.level=0.95,
                                         plot.function=TRUE, nPostSamp=1000){
  n <- length(x)
  
  #sort y according to x
  xy.matrix.sorted <- arrange(as.data.frame(cbind(x, y)), x)
  x <- xy.matrix.sorted[, 1]
  y <- xy.matrix.sorted[, 2]
  
  # set J 
  J <- floor((n ^ (1/3))*log(n))
  
  #prior hyperparameters
  zi <- rep(0,J); tausq <- rep(100, J)
  
  # calculate Nj, YjBar and sigmahat
  obj1 <- extract.nj.sigmahat(x, y, J, zi, tausq)
  nj <- obj1[[1]]
  ybar <- obj1[[2]]
  sigmaHat <- obj1[[3]]
  
  alpha.level <-  (1 - credibility.level) / 2
  corrected.cred <- A.inv(credibility.level)
  correct.level <- (1-corrected.cred) / 2
  
  # posterior mean and variance
  postMean <- (nj*ybar + zi/tausq)/(nj + 1/tausq)
  postSd <- sigmaHat/sqrt(nj + 1/tausq)
  
  if(plot.function==TRUE){
    # plot the data 
    plot(x, y, xlab="x", ylab="y")
    xx <- seq(0.001, 0.999, 0.001)
    function.sample <- NULL
    sample.fx0 <- NULL
    for(i in 1: nPostSamp){
      theta <- rnorm(J, postMean, postSd)
      theta1<- isotone2(theta, nj, x0) 
      sample.fx0 <- c(sample.fx0, theta1[[2]])
      
      yy1 <- theta1[[1]](xx)
      function.sample <- cbind(function.sample, yy1)
      #lines(xx, yy1, xlab="x", ylab="fstar(x)", col="green")
    }
    lines(xx, apply(function.sample,1, FUN=function(x){quantile(x, 1-alpha.level)}))
    lines(xx, apply(function.sample,1, FUN=function(x){quantile(x, alpha.level)}))
    lines(xx, rowMeans(function.sample), lty="dashed")
    legend(0.02, max(y), legend=c(paste(100*credibility.level, "% projection-posterior credible limits"),"Projection-posterior mean"), lty=1:2, cex=0.7)
  }
  
  if(plot.function != TRUE){
    sample.fx0 <- NULL
    # draw samples from projection-posterior
    for(i in 1: nPostSamp){
      theta <- rnorm(J, postMean, postSd)
      theta1<- isotone2(theta, nj, x0) 
      sample.fx0 <- c(sample.fx0, theta1[[2]])
    }
  }
  
  lower.limit <- quantile(sample.fx0, alpha.level)
  upper.limit <- quantile(sample.fx0, 1 - alpha.level)
  interval.obs <- cbind(lower.limit, upper.limit)
  
  lower.limit <- quantile(sample.fx0, correct.level)
  upper.limit <- quantile(sample.fx0, 1 - correct.level)
  interval.correct <- cbind(lower.limit, upper.limit)
  
  return(list(interval.obs, interval.correct))
  
}


#######################################################

#################  monotone density on (0,1)  ##################

#The function `credInt.bounded()' takes $X$ in $(0, 1)$ as input, and gives the unadjusted
#and adjusted projection-posterior credible intervals at a desired target coverage level
# (default is 0.95)

grenEstimateBayes.bdd <- function(y, w, t1){
  J <- length(w)
  w <- w / sum(w) # normalizing the weights
  xval=cumsum(w)
  fit <- gpava(xval, -y)
  r <- stepfun(xval, -c(fit$x, fit$x[J]))
  t <- r(t1) #value of the isotonized step function at given t1
  return(list(r, t))
}


density.credInt.bounded <- function(X, x0=0.4, credibility.level=0.95, sieve.interval.95=T, nPostSamp=1000){
  n <- length(X)
  J <- floor((n^(1/3))*log(n))
  xjKnots <- c(1 : J) * 1 / J
  Nj <- NULL
  for (j in 1:J){
    Nj <- c(Nj, sum(X <= xjKnots[j] & X > (xjKnots[j] - 1/J)))
  }
  alpha <- rep(1, J)
  alpha.level <- (1-credibility.level) / 2
  corrected.cred <- A.inv(credibility.level)
  correct.level <- (1-corrected.cred) / 2
  
  d.x <- density(X)
  d.x.fun <- stepfun(d.x$x, c(d.x$y[1],d.x$y))
  d.x.der <- stepfun(d.x$x[-1], c(0, d.x$y[-1]-d.x$y[-length(d.x$y)])/
                       c(1, d.x$x[-1]-d.x$x[-length(d.x$x)]))
  a <- sqrt(d.x.fun(x0))
  b <- abs(d.x.der(x0))/2
  C0 <- 2*b*(a/b)^(2/3)
  
  gx0.sample <- NULL
  for(i in 1: nPostSamp){
    theta <- rdirichlet(1, alpha+Nj)
    theta.function <- grenEstimateBayes.bdd(c(theta), rep(1, J), x0)
    gx0.sample <- c(gx0.sample, J*theta.function[[2]])
  }
  
  lower.limit <- quantile(gx0.sample, alpha.level)
  upper.limit <- quantile(gx0.sample, 1 - alpha.level)
  #cvrg <- ifelse((upper.limit-d0)*(d0-lower.limit) > 0, 1, 0)
  interval.obs <- cbind(lower.limit, upper.limit)
  
  lower.limit <- quantile(gx0.sample, correct.level)
  upper.limit <- quantile(gx0.sample, 1 - correct.level)
  #cvrg <- ifelse((upper.limit-d0)*(d0-lower.limit) > 0, 1, 0)
  interval.correct <- cbind(lower.limit, upper.limit)
  
  interval.sieve = NULL
  if(sieve.interval.95 == T){
    chern.quant=0.998 # from quantiles of the Chernoff distribution
    sieve.est <- J*grenEstimateBayes.bdd(Nj/n, rep(1, J), x0)[[2]]
    lower.limit <- sieve.est - n^(-1/3)*chern.quant*C0
    upper.limit <- sieve.est + n^(-1/3)*chern.quant*C0
    #cvrg <- ifelse((upper.limit-d0)*(d0-lower.limit) > 0, 1, 0)
    interval.sieve <- cbind(lower.limit, upper.limit)
  }
  
  return(list(interval.obs, interval.correct, interval.sieve))
  #return(list(c(lower.limit, upper.limit), credibility.level))
} 


##################################################################

################ monotone density on the positive real line  #############
# The function `credible.interval.unbounded()' takes the data $X>0$ as input, and gives
# the unadjusted and adjusted projection-posterior credible intervals at a desired 
# target coverage level (default is 0.95)

grenEstimateBayes.unbd <- function(y, w, grid.points){
  J <- length(w)
  w <- w / sum(w) # normalizing the weights
  y1 <- cumsum(y*w) # cumulative sum vector of y
  ll <- gcmlcm(cumsum(w), y1, type = "lcm") # least concave majorant of 
  # the cumulative sum vector of y
  if(length(ll$x.knots)>1){
    y11 <- c(ll$slope.knots[1], ll$slope.knots, ll$slope.knots[length(ll$slope.knots)])
    stepfunction <- stepfun(ll$x.knots, y11)
  }
  if(length(ll$x.knots)<=1){
    stepfunction <- stepfun(c(0), c(ll$slope.knots, ll$slope.knots), right=TRUE)
  }
  fw <- stepfunction(cumsum(w))
  r <- stepfun(grid.points, c(fw[1], fw))
  #t <- r(cumsum(w)[min(which(cumsum(w)>=t1))]) #value of the isotonized step function at given t1
  return(r)
}

credible.interval.unbounded <- function(X, x0=1.5, credibility.level=0.95,
                                        nPostSamp=1000, plot.function=FALSE){
  a <- 2
  n <- length(X)
  Kn <- floor(log(n, base=a)/2)
  J <- floor((n^(1/3))*(log(n)))
  J1 <- J*(Kn+1)
  Tn <- a^Kn
  
  xjKnots.dots <- c(0, a^(c(0, 1:Kn)))
  xjKnots <- NULL#seq(1/J, 1, len=J)
  for(i in 1: (length(xjKnots.dots)-1)){
    xjKnots <- c(xjKnots, seq(xjKnots.dots[i]+1/J, xjKnots.dots[i+1], len=J))
  }
  
  Nj <- sum(X <= xjKnots[1] & X > 0)
  for (j in 2:J1){
    Nj <- c(Nj, sum(X <= xjKnots[j] & X > xjKnots[j-1]))
  }
  alpha <- rep(1, J1)
  wt <- rep(0, J1)
  wt[1:J] <- J
  for(i in 1:Kn){
    wt[(i*J + 1) : ((i+1)*J)] <- J / (a^i - a^(i-1))
  }
  postMean <- alpha+Nj
  alpha.level <- (1-credibility.level) / 2
  
  if(plot.function==TRUE){
    xx <- seq(0.01, Tn+0.5, 0.01)
    xval <- xx
    L0 <- rep(0, length(xval))
    L0[between(xval, 0, 1)] <- 1
    L0[xval==1] <- 1
    for(i in 1:Kn){
      L0[between(xval, a^(i-1), a^i)] <- a^i - a^(i-1)
      L0[xval==a^i] <- a^i - a^(i-1)
    }
    J.by.L0 <- J / L0
    
    hist(X, probability = TRUE)
    function.sample <- NULL
    gx0.sample <- NULL
    L0.x0 <- ifelse(x0>1, a^(ceiling(log(x0, base=a))) - a^(floor(log(x0, base=a))), 1)
    J.by.L0.x0 <- J/L0.x0 
    for(i in 1:nPostSamp){
      theta <- rdirichlet(1, postMean)
      theta1.function <- grenEstimateBayes.unbd(theta, wt, xjKnots)
      theta1 <- J.by.L0*theta1.function(xx)
      #lines(xx, yy11, xlab="x", ylab="f(x)", col="grey")
      function.sample <- cbind(function.sample, theta1)
      gx0.sample <- c(gx0.sample, J.by.L0.x0 * theta1.function(x0))
    }
    lines(xx, apply(function.sample,1, FUN=function(x){quantile(x, alpha.level)}), col="blue", lwd=2)
    lines(xx, apply(function.sample,1, FUN=function(x){quantile(x, 1-alpha.level)}), col="blue", lwd=2)
    lines(xx, rowMeans(function.sample), col="blue", lty="dashed")
    legend(4, 0.6, legend=c(paste(100*credibility.level, "% projection-posterior credible limits"),"Projection-posterior mean"), lty=1:2, cex=0.7)
  }
  
  if(plot.function==FALSE){
    gx0.sample <- NULL
    L0.x0 <- ifelse(x0>1, a^(ceiling(log(x0, base=a))) - a^(floor(log(x0, base=a))), 1)
    J.by.L0.x0 <- J/L0.x0 
    for(i in 1:nPostSamp){
      theta <- rdirichlet(1, postMean)
      theta1.function <- grenEstimateBayes.unbd(theta, wt, xjKnots)
      gx0.sample <- c(gx0.sample, J.by.L0.x0 * theta1.function(x0))
    }
  }
  
  lower.limit <- quantile(gx0.sample, alpha.level)
  upper.limit <- quantile(gx0.sample, 1 - alpha.level)
  print(paste(credibility.level*100,"% point-wise posterior-projection credible interval= (",
              sep=" ", lower.limit,",", upper.limit, ")"))
  #return(list(c(lower.limit, upper.limit), credibility.level))
}

######################################################################

######## Estimation of the full monotone regression function on (0,1) #########

# returns posterior samples of the step heights of f
pr.posterior.samples <- function(x, y, Jsize="n power 1/3", nPostSamp=1000){
  if(length(x)!=length(y)) stop("X and Y lengths are different")
  n=length(x)
  #set J
  J <- floor(n ^ (1/3))
  if(Jsize!="n power 1/3")
    J <- floor((n ^ (1/3))*log(n))
  
  # prior hyperparameters
  zi <- rep(0,J); tausq <- rep(100, J)
  
  xbar <- c(1:J) / J; nj <- c(1:J)
  xInIj <- list()
  ybar <- c(1:J)
  sigmaHatSq <- NULL
  for(j in 1:J){
    xInIj[[j]] <- x > (j-1)/J & x <= j/J #indices present in j-th interval.
    nj[j] <- sum(xInIj[[j]]) #number of Xj's in j-th interval.
    yInIj <- y[xInIj[[j]]]
    ybar[j] <- mean(yInIj)
    sigmaHatSqIj <- t(yInIj-zi[j])%*%solve(diag(nj[j])+
                                             matrix(tausq[j],nj[j],nj[j]))%*%(yInIj-zi[j])
    sigmaHatSq <- c(sigmaHatSq, sigmaHatSqIj)
  }
  sigmaHatSq <- sum(sigmaHatSq)/n
  sigmaHat <- sqrt(sigmaHatSq)
  njByn <- nj/n
  distSamp <- NULL
  
  #obtain posterior mean and variance.
  postMean <- (nj*ybar + zi/tausq)/(nj + 1/tausq)
  postSd <- sigmaHat/sqrt(nj + 1/tausq) 
  postSamp <- NULL # initializing vector of posterior samples.
  
  for(k in 1:nPostSamp){
    #loop for creating posterior samples
    theta <- rnorm(J, postMean, postSd) #generate posterior sample.
    fit.iso <- gpava(z=xbar, y=theta, solver=weighted.median, weights = nj)
    fstar.pieces <- fit.iso$x
    postSamp <- rbind(postSamp, fstar.pieces) # add fstar to the list of posterior samples
  }
  return(list(postSamp, xbar))
}

# values of the step function on a given point
predict.iso <- function(x.values, yval.mat, x.bar){
  if(ncol(yval.mat)!=length(x.bar)) stop("Dimensions do not match")
  fstar.values <- NULL
  for(i in 1:nrow(yval.mat)){
    fstar <- stepfun(x.bar, c(yval.mat[i, 1], yval.mat[i, ]))(x.values)
    fstar.values <- rbind(fstar.values, fstar)
  }
  return(fstar.values)
}

# L1 distance between the two step functions
dist.cal <- function(f.vec, x.vec, ff){
  k <- length(x.vec)
  if(length(f.vec)!=k) stop("Lengths do not match")
  return(sum(abs(f.vec - ff(x.vec)))/k)
}

############################################################################

############ Hypothesis testing for monotonicity of regression function #############

testL1 <- function(x, y, nPostSamp=1000){
  if(length(x)!=length(y)) stop("X and Y lengths are different")
  n=length(x)
  #set J
  J <- floor(n ^ (1/3))
  #J <- floor((n ^ (1/3))*log(n))
  
  # prior hyperparameters
  zi <- rep(0,J); tausq <- rep(100, J)
  
  xbar <- c(1:J) / J; nj <- c(1:J)
  xInIj <- list()
  ybar <- c(1:J)
  sigmaHatSq <- NULL
  for(j in 1:J){
    xInIj[[j]] <- x > (j-1)/J & x <= j/J #indices present in j-th interval.
    nj[j] <- sum(xInIj[[j]]) #number of Xj's in j-th interval.
    yInIj <- y[xInIj[[j]]]
    ybar[j] <- mean(yInIj)
    sigmaHatSqIj <- t(yInIj-zi[j])%*%solve(diag(nj[j])+
                                             matrix(tausq[j],nj[j],nj[j]))%*%(yInIj-zi[j])
    sigmaHatSq <- c(sigmaHatSq, sigmaHatSqIj)
  }
  sigmaHatSq <- sum(sigmaHatSq)/n
  sigmaHat <- sqrt(sigmaHatSq)
  njByn <- nj/n
  distSamp <- NULL
  
  #obtain posterior mean and variance.
  postMean <- (nj*ybar + zi/tausq)/(nj + 1/tausq)
  postSd <- sigmaHat/sqrt(nj + 1/tausq) 
  
  for(k in 1:nPostSamp){
    #loop for creating posterior samples
    theta <- rnorm(J, postMean, postSd) #generate posterior sample.
    fit.iso <- gpava(z=xbar, y=theta, solver=weighted.median, weights = nj)
    #fstar <- sum(abs(fit.iso$y-fit.iso$x))
    fstar <- sum(abs(fit.iso$y-fit.iso$x) * njByn)
    distSamp <- c(distSamp, fstar) # add fstar to the list of distances
  }
  #M0 <- 0.09283178
  tau.cut <- 0.8*sigmaHat*log(n)^0.1/(n^(1/3)) 
  #tau.cut <- 0.7*sigmaHat*n^(1/15)/(n^(1/3)) 
  return(sum(distSamp > tau.cut)/nPostSamp >=0.5)
}

###############################################################################

############### single-stage projection-posterior for regression quantile ###########

reg.quantile.projection = function(x, y, d0 = 0.35, credibility.level=0.95,nPostSamp=1000){ 
  
  n <- length(x)
  
  #sort y according to x
  xy.matrix.sorted <- arrange(as.data.frame(cbind(x, y)), x)
  x <- xy.matrix.sorted[, 1]
  y <- xy.matrix.sorted[, 2]
  
  # set J 
  J <- floor((n ^ (1/3))*log(n))
  
  #prior hyperparameters
  zi <- rep(0,J); tausq <- rep(100, J)
  
  # calculate Nj, YjBar and sigmahat
  obj1 <- extract.nj.sigmahat(x, y, J, zi, tausq)
  nj <- obj1[[1]]
  ybar <- obj1[[2]]
  sigmaHat <- obj1[[3]]
  
  alpha.level <-  (1 - credibility.level) / 2
  corrected.cred <- A.inv(credibility.level)
  correct.level <- (1-corrected.cred) / 2
  
  # posterior mean and variance
  postMean <- (nj*ybar + zi/tausq)/(nj + 1/tausq)
  postSd <- sigmaHat/sqrt(nj + 1/tausq)
  
  postSamp = NULL 
  for(i in 1:nPostSamp){ 
    theta <- rnorm(J, postMean, postSd) #generate posterior sample.
    fstar.0 <- isotone.f(theta, nj, d0)(seq(1/J,1,1/J))
    fstar <- ifelse(sum(fstar.0>=d0)>0, min(which(fstar.0>=d0))/J - 1/J, 1)# fstar is projected posterior value at t0.
    #fstar.unrestricted = ifelse(sum(theta<=d0)>0, max(which(theta<=d0))/J, 0)
    postSamp <- c(postSamp, fstar) # add fstar to the list of posterior samples.
  }
  lower.limit <- quantile(postSamp, alpha.level)
  upper.limit <- quantile(postSamp, 1 - alpha.level)
  interval.obs <- cbind(lower.limit, upper.limit)
  
  lower.limit <- quantile(postSamp, correct.level)
  upper.limit <- quantile(postSamp, 1 - correct.level)
  interval.correct <- cbind(lower.limit, upper.limit)
  
  return(list(interval.obs, interval.correct, sigmaHatSquare=sigmaHat^2))
}#end of 1stage projection-posterior function

###########################################################

################# regression quantile two-stage posterior samples #############

# function to obtain posterior samples of regression quantile from the second-stage data
reg.quantile.2stage = function(x, y, d0 = d0, mu.tilde, sigmaHatSq, nPostSamp=1000){
  
  n2 = length(x)
  X <- cbind(rep(1, n2), x-mu.tilde)
  
  mu0=c(0, 0); Lambda0=diag(100,2); 
  Lambda0[2,2] = ifelse(max(x) - mu.tilde > 0.001,
                        (max(x) - mu.tilde)^(-2), Lambda0[2,2])
  Lambda0[1,1]=1
  
  LambdaNInv <- solve(t(X)%*%X+solve(Lambda0))
  muN <- LambdaNInv%*%(t(X)%*%y + solve(Lambda0)%*%mu0)
  postSamp1 <- NULL
  for(i in 1:nPostSamp){
    beta <- rmvnorm(1,muN, sigmaHatSq*LambdaNInv)
    v <- (d0 - beta[1])/beta[2] + mu.tilde
    v1 <- ifelse(v<=1, v, 1)
    v1 <- ifelse(v1>=0, v1, 0)
    postSamp1 <- c(postSamp1, v1)
  }
  return(postSamp1)
}

