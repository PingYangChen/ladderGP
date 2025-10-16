
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
                       categParRange = c(0.15, 0.5), 
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
    
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")
    
    if (zType == "o") {
      nContiPar <- ncol(x) + 1
      low_bound <- c(rep(min(contiParRange), nContiPar))
      upp_bound <- c(rep(max(contiParRange), nContiPar))
      
      res <- globpso(objFunc = lgpOBObjCpp, lower = low_bound, upper = upp_bound,
                     PSO_INFO = alg_setting, verbose = optVerbose,
                     y = y, x = x, z = z, xzDim = xzDim, nugget = nugget)
      #res$val
      mdl <- lgpOBModel(param = res$par, y = y, x = x, z = z, xzDim = xzDim, nugget = nugget)
      
    } else if (zType == "n") {
      nContiPar <- ncol(x)
      nCategPar <- (0.5*max(z)*(max(z) - 1))
      low_bound <- c(rep(min(contiParRange), nContiPar),
                     rep(min(categParRange), nCategPar))
      upp_bound <- c(rep(max(contiParRange), nContiPar),
                     rep(max(categParRange), nCategPar))
      
      res <- globpso(objFunc = lgpNBObjCpp, lower = low_bound, upper = upp_bound,
                     PSO_INFO = alg_setting, verbose = optVerbose,
                     y = y, x = x, z = z, xzDim = xzDim, nugget = nugget)
      #res$val
      mdl <- lgpNBModel(param = res$par, y = y, x = x, z = z, xzDim = xzDim, nugget = nugget)
    } else {
      stop(sprintf("ERROR: zType should be 'n' or 'o'\n"))
    }
    mdl$zType <- zType
    mdl$data <- list(y = y, x = x, z = z, xzDim = xzDim)
  })[3]
  mdl$cputime <- cputime
  cat(sprintf("mLadderGP(%s) FIT CPU time: %.2f seconds.\n", zType, cputime))
  return(mdl)
}


mLadderPred <- function(gpMdl, x0List, ei_alpha = 0.5, min_y = NULL) {
  
  cputime <- system.time({
    
    xDims <- sapply(1:length(x0List), function(k) ncol(x0List[[k]]))
    xzDim <- min(xDims)
    zs <- xDims/xzDim
    x0 <- matrix(0, nrow = 0, ncol = max(xDims))
    z0 <- dimCheck <- c()
    for (i in 1:length(yList)) {
      n <- nrow(x0List[[i]])
      x0 <- rbind(x0, cbind(x0List[[i]], matrix(-1, n, ncol(x0) - ncol(x0List[[i]]))))
      z0 <- c(z0, rep(xDims[i]/xzDim, n))
      dimCheck[i] <- xDims[i] %% xzDim
    }
    stopifnot(all(dimCheck == 0))
    
    if (is.null(min_y)) { min_y <- min(gpMdl$data$y) }
    
    zType <- gpMdl$zType
    if (zType == "o") {
      pred <- lgpOBPred(x0, z0, gpMdl$data$y, gpMdl$data$x, gpMdl$data$z, gpMdl$data$xzDim,
                        gpMdl$vecParams, gpMdl$invPsi, gpMdl$mu, gpMdl$sigma2, ei_alpha, min_y)
    } else if (zType == "n") {
      pred <- lgpNBPred(x0, z0, gpMdl$data$y, gpMdl$data$x, gpMdl$data$z, gpMdl$data$xzDim,
                        gpMdl$vecParams, gpMdl$invPsi, gpMdl$mu, gpMdl$sigma2, ei_alpha, min_y)
    }
   
  })[3]
  return(pred)
} 

