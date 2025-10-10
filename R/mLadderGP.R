
library(Rcpp)
library(RcppArmadillo)
library(globpso)

sourceCpp(file.path('src/cppFunc.cpp'))


################################################################################
### Multiplicative Ladder GP
################################################################################
# y: vector
# x: matrix
# z: integer vector
mLadderFit <- function(yList, xList, zType = c("o", "n"),
                       contiParRange = 10^c(-3, .5), 
                       varParRange = 10^c(-3, .5),
                       nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE) {
 
  cputime <- system.time({
    xDims <- sapply(1:length(xList), function(k) ncol(xList[[k]]))
    xzDim <- min(xDims)
    zs <- xDims/xzDim
    x <- matrix(0, nrow = 0, ncol = max(xDims))
    y <- z <- dimCheck <- c()
    for (i in 1:length(yList)) {
      n <- length(yList[[i]])
      x <- rbind(x, cbind(xList[[i]], matrix(-1, n, ncol(x) - ncol(xList[[i]]))))
      y <- c(y, yList[[i]])
      z <- c(z, rep(xDims[i]/xzDim, n))
      dimCheck[i] <- xDims[i] %% xzDim
    }
    stopifnot(all(dimCheck == 0))
    
    nContiPar <- ncol(x) + 1
    nVarPar <- max(z) + 1
    low_bound <- c(rep(min(contiParRange), nContiPar),
                   rep(min(varParRange), nVarPar))
    upp_bound <- c(rep(max(contiParRange), nContiPar),
                   rep(max(varParRange), nVarPar))
    
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")

    res <- globpso(objFunc = ahgpObjCpp, lower = low_bound, upper = upp_bound,
                   PSO_INFO = alg_setting, verbose = optVerbose,
                   y = y, x = x, z = z, xzDim = xzDim, nugget = nugget)
    #res$val
    mdl <- ahgpModel(param = res$par, y = y, x = x, z = z, xzDim = xzDim, nugget = nugget)
    mdl$data <- list(y = y, x = x, z = z, xzDim = xzDim)
  })[3]
  mdl$cputime <- cputime
  cat(sprintf("AHGP FIT CPU time: %.2f seconds.\n", cputime))
  return(mdl)
}


mLadderPred <- function(gpMdl, x0List, ei_alpha = 0.5, min_y = NULL) {
  
  cputime <- system.time({
    
    xDims <- sapply(1:length(x0List), function(k) ncol(x0List[[k]]))
    xzDim <- min(xDims)
    zs <- xDims/xzDim
    x <- matrix(0, nrow = 0, ncol = max(xDims))
    y <- z <- dimCheck <- c()
    for (i in 1:length(yList)) {
      n <- length(yList[[i]])
      x <- rbind(x, cbind(x0List[[i]], matrix(-1, n, ncol(x) - ncol(x0List[[i]]))))
      y <- c(y, yList[[i]])
      z <- c(z, rep(xDims[i]/xzDim, n))
      dimCheck[i] <- xDims[i] %% xzDim
    }
    stopifnot(all(dimCheck == 0))
    
    
    if (is.null(min_y)) { min_y <- min(gpMdl$data$y) }
    
    z0 <- ncol(x0)/gpMdl$data$xzDim
    
    pred <- ahgpPred(x0, z0, gpMdl$data$y, gpMdl$data$x, gpMdl$data$z, gpMdl$data$xzDim,
                     gpMdl$vecParams, gpMdl$invPsi, gpMdl$mu, ei_alpha, min_y)
  })[3]
  return(pred)
} 



mLadderMaxEi <- function(gp, ei_alpha = 0.5, min_y = NULL, nSwarm = 64, maxIter = 200, optVerbose = TRUE) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  cputime <- system.time({
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")
    #
    low_bound <- rep(0, gp$data$xDim)
    upp_bound <- rep(1, gp$data$xDim)
    #
    res <- globpso(objFunc = ahgpMaxEiObj, lower = low_bound, upper = upp_bound,
                   PSO_INFO = alg_setting, verbose = optVerbose,
                   gp = gp, ei_alpha = ei_alpha, min_y = min_y)
    
    rx <- matrix(res$par, 1, gp$data$xDim)
    rdata <- list(x = rx)
  })[3]
  cat(sprintf("AHGP MAXEI CPU time: %.2f seconds.\n", cputime))
  return(list(eiVal = exp(-res$val),
              newpoint = rdata,
              cputime = cputime))
}

MladderMaxEiObj <- function(xx, gp, ei_alpha = 0.5, min_y = NULL) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  pred <- ahgpPred(xx, gp$data$y, gp$data$x, 
                   gp$alpha, gp$invPsi, gp$mu, gp$sigma, ei_alpha, min_y)
  
  return( -log(pred$ei[1,1]) )
}


