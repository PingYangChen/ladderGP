### 引入 Demo 須要用的 Space Filling Design 套件
library(SFDesign)
### 引入測試用函數集
source('R/testfunction.R')
### 引入 GP 模組
source('R/nvladderGP.R')

### 設定各階段資料維度
p_data <- c(2, 4, 6)     # 各階段實驗的因子數
n_train <- c(10, 10, 10) # 訓練資料各階段實驗的樣本數

### 使用 maxproLHD 產生測試資料的解釋變數
xList <- lapply(1:length(p_data), function(k) {
  maxproLHD(n_train[k], p_data[k])$design
})

### 以 Rastrigin 函數產生訓練資料各階段模擬資料的反應變數 
yList <- lapply(1:length(p_data), function(k) {
  tmp <- numeric(nrow(xList[[k]]))
  for (i in 1:nrow(xList[[k]])) {
    tmp[i] <- Rastrigin(xList[[k]][i,])
  }
  tmp
})

### 測試資料各階段實驗的樣本數
n_test <- c(3, 3, 3) 

### 使用 maxproLHD 產生測試資料的解釋變數
x0List <- lapply(1:length(p_data), function(k) {
  maxproLHD(n_test[k], p_data[k])$design
})

### 以 Rastrigin 函數產生測試資料各階段模擬資料的反應變數 
y0List <- lapply(1:length(p_data), function(k) {
  tmp <- numeric(nrow(x0List[[k]]))
  for (i in 1:nrow(x0List[[k]])) {
    tmp[i] <- Rastrigin(x0List[[k]][i,])
  }
  tmp
})

### 建立不包含描述兩階段數相異實驗資料間的關聯參數之 GP 模型 
nvLadderMdl <- nvLadderFit(yList, xList, contiParLogRange = c(-6.5, 1.5),
                           nSwarm = 64, maxIter = 200, psoType = "basic", nugget = 1e-6, optVerbose = FALSE)

### 以 nvLadderMdl GP 模型對測試資料進行預測並檢視 RMSE
nvPred <- nvLadderPred(nvLadderMdl, x0List, y0listTrue = y0List)
cat(sprintf("RMSE(nvLadderMdl) = %.4f\n", sqrt(sum((nvPred$pred - nvPred$y_true)^2))))

nvLadderMdl$nugget
nvLadderMdl$negloglik
det(nvLadderMdl$psi)
det(nvLadderMdl$invPsi)


