source('R/models.R')
source('R/testfunction.R')

#library(lhs)
library(SFDesign)

p_scen <- c(3, 6, 9)
n_scen <- c(10, 10, 10)

init_d <- lapply(1:length(p_scen), function(k) {
  maxproLHD(n_scen[k], p_scen[k])$design
})

init_y <- lapply(1:length(p_scen), function(k) {
  tmp <- numeric(nrow(init_d[[k]]))
  for (i in 1:nrow(init_d[[k]])) {
    tmp[i] <- Rastrigin(init_d[[k]][i,])
  }
  tmp
})


xList <- init_d
yList <- init_y
cbind(c(yList[[1]], yList[[2]], yList[[3]]),
      rbind(cbind(xList[[1]], matrix(-1, 10, 6)),
            cbind(xList[[2]], matrix(-1, 10, 3)), xList[[3]]), rep(1:3, each = 10)
)

cbind(c(init_y[[1]], init_y[[2]], init_y[[3]]),
rbind(cbind(init_d[[1]], matrix(NA, 10, 6)),
      cbind(init_d[[2]], matrix(NA, 10, 3)), init_d[[3]]), rep(1:3, each = 10)
)
exp_x <- init_x
exp_y <- init_y




gpm <- gpFit(exp_y, exp_x, contiParRange = 10^c(-3, .5), 
             nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE)


gpm
