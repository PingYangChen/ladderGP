

mLadderMaxEi <- function(gp, ei_alpha = 0.5, min_y = NULL, nSwarm = 64, maxIter = 200, optVerbose = TRUE) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  cputime <- system.time({
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")
    #
    low_bound <- rep(0, gp$data$xDim)
    upp_bound <- rep(1, gp$data$xDim)
    #
    res <- globpso(objFunc = mladderMaxEiObj, lower = low_bound, upper = upp_bound,
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

mladderMaxEiObj <- function(xx, gp, ei_alpha = 0.5, min_y = NULL) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  
  zType <- gp$zType
  if (zType == "o") {
    pred <- lgpOBPred(x0, z0, gp$data$y, gp$data$x, gp$data$z, gp$data$xzDim,
                      gp$vecParams, gp$invPsi, gp$mu, ei_alpha, min_y)
  } else if (zType == "n") {
    pred <- lgpNBPred(x0, z0, gp$data$y, gp$data$x, gp$data$z, gp$data$xzDim,
                      gp$vecParams, gp$invPsi, gp$mu, ei_alpha, min_y)
  }
  return( -log(pred$ei[1,1]) )
}



aLadderMaxEi <- function(gp, ei_alpha = 0.5, min_y = NULL, nSwarm = 64, maxIter = 200, optVerbose = TRUE) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  cputime <- system.time({
    alg_setting <- getPSOInfo(nSwarm = nSwarm, maxIter = maxIter, psoType = "quantum")
    #
    low_bound <- rep(0, gp$data$xDim)
    upp_bound <- rep(1, gp$data$xDim)
    #
    res <- globpso(objFunc = aLadderMaxEiObj, lower = low_bound, upper = upp_bound,
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

aLadderMaxEiObj <- function(xx, gp, ei_alpha = 0.5, min_y = NULL) {
  
  if (is.null(min_y)) { min_y <- min(gp$data$y) }
  pred <- aIntPred(xx, gp$data$y, gp$data$x, 
                   gp$alpha, gp$invPsi, gp$mu, gp$sigma, ei_alpha, min_y)
  
  return( -log(pred$ei[1,1]) )
}


