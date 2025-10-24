
source('R/aladderGP.R')
source('R/testfunction.R')

library(SFDesign)

p_scen <- c(3, 6, 9)
n_scen <- c(10, 10, 10)

xList <- lapply(1:length(p_scen), function(k) {
  maxproLHD(n_scen[k], p_scen[k])$design
})

yList <- lapply(1:length(p_scen), function(k) {
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


aLadderMdl <- aLadderFit(yList, xList, 
                         contiParRange = 10^c(-3, .5), varParRange = 10^c(-3, .5), 
                         nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE)

aLadderMdl

names(aLadderMdl)


aLadderPred(aLadderMdl, x0List)$pred

