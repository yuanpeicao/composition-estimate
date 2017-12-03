set.seed(150)
## real data analysis
setwd("~/comp-est-R")
source("tune_proximal_gradient.R")
source("basic_functions.R")
require(foreach)
require(doParallel)

# load count data
d <- as.matrix(read.table("/Users/yuanpeicao/comp-est-R/dataset/BMI_and_Counts.txt", sep = " ", header = TRUE))
w <- d[, 3:89]
w_b <- w > 0
freq <- colSums(w_b)
sumCount <- colSums(w)

# # filter out rarest species
# w_select <- w
# w_select <- w[, -which(freq < 10)]
w_select <- w[, -which(freq < 10 | freq == 98)]
write.csv(w_select, file = '/Users/yuanpeicao/comp-est-R/real-data-result/w_select.csv')

p <- ncol(w_select)
X_select <- w_select/(rowSums(w_select)%*%matrix(1,1,p))

## proximal gradient 
X_pg_res <- autoTuneProxGradient(W = w_select, n_grid = 10)
X_pg_res_hat <- X_pg_res$X_hat
write.csv(X_pg_res_hat, file = '/Users/yuanpeicao/comp-est-R/real-data-result/Xhat.csv')

## ZR
X_zr_res <- autoTuneZr(w_select)
X_zr_res_hat <- X_zr_res$X_hat
write.csv(X_zr_res_hat, file = '/Users/yuanpeicao/comp-est-R/real-data-result/Xhat_zr.csv')

## SVT
X_svt_res <- autoTuneSvt(W = w_select, r = 10)
X_svt_res_hat <- X_svt_res$X_hat
write.csv(X_svt_res_hat, file = '/Users/yuanpeicao/comp-est-R/real-data-result/Xhat_svt.csv')