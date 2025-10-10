

library(Rcpp)
library(RcppArmadillo)
library(globpso)

sourceCpp(file.path('src/cppFunc.cpp'))

################################################################################
### GP
################################################################################
# y: vector
# x: matrix
gpFit <- function(y, x, 
                  contiParRange = 10^c(-3, .5), 
                  nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE) {
  
  cputime <- system.time({
    xDim <- ncol(x)
    nContiPar <- xDim
    
    low_bound <- rep(min(contiParRange), nContiPar)
    upp_bound <- rep(max(contiParRange), nContiPar)
    
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")
    
    res <- globpso(objFunc = gpObjCpp, lower = low_bound, upper = upp_bound,
                   PSO_INFO = alg_setting, verbose = optVerbose,
                   y = y, x = x, nugget = nugget)
    #res$val
    mdl <- gpModel(param = res$par, y = y, x = x, nugget = nugget)
    mdl$data <- list(y = y, x = x, xDim = xDim)
  })[3]
  mdl$cputime <- cputime
  cat(sprintf("GP FIT CPU time: %.2f seconds.\n", cputime))
  return(mdl)
}

gpPredict <- function(gp, x0, ei_alpha = 0.5, min_y = NULL) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  pred <- gpPred(x0, gp$data$y, gp$data$x, 
                 gp$alpha, gp$invPsi, gp$mu, gp$sigma, ei_alpha, min_y)
  return(pred)
} 


################################################################################
### GP
################################################################################
gpMaxEi <- function(gp, ei_alpha = 0.5, min_y = NULL, nSwarm = 64, maxIter = 200, optVerbose = TRUE) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  cputime <- system.time({
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")
    #
    low_bound <- rep(0, gp$data$xDim)
    upp_bound <- rep(1, gp$data$xDim)
    #
    res <- globpso(objFunc = gpMaxEiObj, lower = low_bound, upper = upp_bound,
                   PSO_INFO = alg_setting, verbose = optVerbose,
                   gp = gp, ei_alpha = ei_alpha, min_y = min_y)
    
    rx <- matrix(res$par, 1, gp$data$xDim)
    rdata <- list(x = rx)
  })[3]
  cat(sprintf("BNGP MAXEI CPU time: %.2f seconds.\n", cputime))
  return(list(eiVal = exp(-res$val),
              newpoint = rdata,
              cputime = cputime))
}

gpMaxEiObj <- function(xx, gp, ei_alpha = 0.5, min_y = NULL) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  pred <- gpPred(xx, gp$data$y, gp$data$x, 
                 gp$alpha, gp$invPsi, gp$mu, gp$sigma, ei_alpha, min_y)
  
  return( -log(pred$ei[1,1]) )
}


