require(ggplot2)
require(reshape)
require(gridExtra)
require(grid)
require(cowplot)
require(latex2exp)
setwd("~/comp-est-R")
source("naive_methods.R")

####################################################################################################################
# Load count matrix
####################################################################################################################
wSelect <- as.matrix(read.csv('/Users/yuanpeicao/comp-est-R/real-data-result/w_select.csv'))
wSelect <- wSelect[, 2:ncol(wSelect)]

####################################################################################################################
# Load proximal gradient or cvx results
####################################################################################################################
## from cvx result
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.1_c_0.25_p42.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.5_c_0.1_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.1_c_0.5_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.5_c_0.5_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_316_c_1e_100_p45.dat', header = FALSE))
# xHatRare <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_09_c_1e_5_p42.dat', header = FALSE))
xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_05_c_1e_5_p45.dat', header = FALSE))

## from proximal gradient result
# xHat <- as.matrix(read.csv('/Users/yuanpeicao/comp-est-R/real-data-result/Xhat.csv'))
# xHat <- xHat[, 2:ncol(xHat)]

## re-correct the composition for the species whose original counts are 0s
# replace composition estimiation from zero count by max(composition estimation, 0.5 * minimum row minimum estimated composition from non-zero count index)
xHat_tmp <- xHat
xHat_tmp[wSelect == 0] <- 100
rowMins <- matrix(apply(xHat_tmp, 1, min), nrow(xHat_tmp), 1)
replaceMatrix <- rowMins %*% matrix(1,1,ncol(xHat_tmp))
xHat[wSelect == 0] <- pmax(replaceMatrix[wSelect == 0] * 0.5, xHat[wSelect == 0])
xHat <- xHat/(rowSums(xHat)%*%matrix(1,1,ncol(xHat)))

####################################################################################################################
# Load zero replacement result
####################################################################################################################
xHat_zr <- zr(W = wSelect, alpha = 0.5)

####################################################################################################################
# Load SVT result
####################################################################################################################
xHat_svt <- svt(W = wSelect, r = 10, alpha = 0.5)

####################################################################################################################
# filter out common species
####################################################################################################################
w_b <- wSelect > 0
freq <- colSums(w_b)

xHatRare <- xHat[, -which(freq == nrow(wSelect))]
wSelectRare <- wSelect[, -which(freq == nrow(wSelect))]

xHatZrRare <- xHat_zr[, -which(freq == nrow(wSelect))]
xHatSvtRare <- xHat_svt[, -which(freq == nrow(wSelect))]

####################################################################################################################
# Load the species name
####################################################################################################################
speciesNames <- colnames(wSelectRare)
taxa <- matrix(unlist(strsplit(speciesNames, "\\.")), 6)
phyla <- taxa[2, ]
genera <- taxa[6, ]

####################################################################################################################
# Split the compsition into rare and common
####################################################################################################################
xZeroIndex <- (wSelectRare == 0)

n <- nrow(xHatRare)
p <- ncol(xHatRare)

####################################################################################################################
# calculate the sequence depth
####################################################################################################################
depth <- log(rowSums(wSelect))
depth_mat <- depth %*% matrix(1,1,ncol(xHatRare))

####################################################################################################################
# record the composition estimtator from zero count
####################################################################################################################
seq_depth <- c()
xHat_v <- c()
xHat_svt_v <- c()
xHat_zr_v <- c()

for (i in 1:n){
  for (j in 1:p){
    if (wSelectRare[i,j] == 0){
      seq_depth <- c(seq_depth, depth_mat[i,j])
      xHat_v <- c(xHat_v, xHatRare[i,j])
      xHat_svt_v <- c(xHat_svt_v, xHatSvtRare[i,j])
      xHat_zr_v <- c(xHat_zr_v, xHatZrRare[i,j])
    }
  }
}

####################################################################################################################
# plot proximal gradient result
####################################################################################################################
df <- data.frame(seq_depth = seq_depth, xHat_v = xHat_v, d = densCols(seq_depth, xHat_v, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

p <- ggplot(df) + geom_hex(aes(seq_depth, xHat_v), bins = 100) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))
print(p)

####################################################################################################################
# plot SVT result
####################################################################################################################
df_svt <- data.frame(seq_depth = seq_depth, xHat_svt_v = xHat_svt_v, d = densCols(seq_depth, xHat_svt_v, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

p_svt <- ggplot(df_svt) + geom_hex(aes(seq_depth, xHat_svt_v), bins = 100) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))
print(p_svt)

####################################################################################################################
# plot zero replacement result
####################################################################################################################
df_zr <- data.frame(seq_depth = seq_depth, xHat_zr_v = xHat_zr_v, d = densCols(seq_depth, xHat_svt_v, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))

p_zr <- ggplot(df_zr) + geom_hex(aes(seq_depth, xHat_zr_v), bins = 100) +
  scale_fill_gradientn("", colours = rev(rainbow(10, end = 4/6)))
print(p_zr)