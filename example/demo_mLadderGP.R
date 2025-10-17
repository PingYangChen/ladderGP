
source('R/mladderGP.R')
source('R/testfunction.R')

library(SFDesign)

p_data <- c(1, 2, 3)
n_train <- c(10, 10, 10)
n_test <- c(3, 3, 3)

xList <- lapply(1:length(p_data), function(k) {
  maxproLHD(n_train[k], p_data[k])$design
})

yList <- lapply(1:length(p_data), function(k) {
  tmp <- numeric(nrow(xList[[k]]))
  for (i in 1:nrow(xList[[k]])) {
    tmp[i] <- Rastrigin(xList[[k]][i,])
  }
  tmp
})

x0List <- lapply(1:length(p_data), function(k) {
  maxproLHD(n_test[k], p_data[k])$design
})
y0List <- lapply(1:length(p_data), function(k) {
  tmp <- numeric(nrow(x0List[[k]]))
  for (i in 1:nrow(x0List[[k]])) {
    tmp[i] <- Rastrigin(x0List[[k]][i,])
  }
  tmp
})



mLadderMdl_o <- mLadderFit(yList, xList, zType = "o", 
                           contiParRange = 10^c(-3, .5), categParRange = c(0.15, 0.5), 
                           nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE)

mLadderPred(mLadderMdl_o, x0List)$pred


mLadderMdl_n <- mLadderFit(yList, xList, zType = "n", 
                           contiParRange = 10^c(-3, .5), categParRange = c(0.15, 0.5), 
                           nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE)

mLadderPred(mLadderMdl_n, x0List)$pred

