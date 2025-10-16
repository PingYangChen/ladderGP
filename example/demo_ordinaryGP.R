source('R/ordinaryGP.R')
source('R/testfunction.R')

#library(lhs)
library(SFDesign)

init_n <- 10
init_x <- maxproLHD(init_n, 2)$design

init_y <- numeric(init_n)
for (i in 1:init_n) {
  init_y[i] <- BarninHoo(init_x[i,])
}

exp_x <- init_x
exp_y <- init_y

gpm <- gpFit(exp_y, exp_x, contiParRange = 10^c(-3, .5), 
             nSwarm = 64, maxIter = 200, nugget = 0., optVerbose = TRUE)


gpm
