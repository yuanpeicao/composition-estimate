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
# plot proximal gradient result
####################################################################################################################
dat <- data.frame(Genera = c(rep(genera,each = n)),zero_index = factor(c(xZeroIndex)), Composition = c(log(xHatRare)))
dat$zero_index <- factor(dat$zero_index, labels = c(expression("italic(M)^c","italic(M)")))

## boxplot
setEPS(width = 18, height = 10)
postscript(file="boxplotSR.eps")
ggplot(dat, aes(x = Genera, y = Composition), mar = c(0.5,0.5,0.5,0.5)) +
  geom_boxplot(aes(fill = zero_index), coef = 10) +
  ylab(expression(log(widehat(X)))) + xlab("Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1), legend.title = element_blank(), legend.text = element_text(size = 30)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"),
                    labels = c(expression(Omega),expression(Omega^c))) +
  theme(axis.title=element_text(size = 30), legend.text.align = 0, strip.text = element_text(size = 30), panel.grid.major = element_blank(),
        panel.background = element_blank(), panel.grid.minor = element_blank())
dev.off()

####################################################################################################################
# plot zero replacement result
####################################################################################################################
dat_zr <- data.frame(Genera = c(rep(genera,each = n)),zero_index = factor(c(xZeroIndex)), Composition = c(log(xHatZrRare)))
dat_zr$zero_index <- factor(dat$zero_index, labels = c(expression("italic(M)^c","italic(M)")))

## boxplot
setEPS(width = 18, height = 10)
postscript(file="boxplotSR_zr.eps")
ggplot(dat_zr, aes(x = Genera, y = Composition), mar = c(0.5,0.5,0.5,0.5)) + 
  geom_boxplot(aes(fill = zero_index), coef = 10) +
  ylab(expression(log(widehat(X)))) + xlab("Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1), legend.title = element_blank(), legend.text = element_text(size = 30)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"),
                    labels = c(expression(Omega),expression(Omega^c))) +
  theme(axis.title=element_text(size = 30), legend.text.align = 0, strip.text = element_text(size = 30), panel.grid.major = element_blank(),
        panel.background = element_blank(), panel.grid.minor = element_blank())
dev.off()

####################################################################################################################
# SVT result
####################################################################################################################
dat_svt <- data.frame(Genera = c(rep(genera,each = n)),zero_index = factor(c(xZeroIndex)), Composition = c(log(xHatSvtRare)))
dat_svt$zero_index <- factor(dat$zero_index, labels = c(expression("italic(M)^c","italic(M)")))

## boxplot
setEPS(width = 18, height = 10)
postscript(file="boxplotSR_svt.eps")
ggplot(dat_svt, aes(x = Genera, y = Composition), mar = c(0.5,0.5,0.5,0.5)) + 
  geom_boxplot(aes(fill = zero_index), coef = 10) +
  ylab(expression(log(widehat(X)))) + xlab("Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1), legend.title = element_blank(), legend.text = element_text(size = 30)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"),
                    labels = c(expression(Omega),expression(Omega^c))) +
  theme(axis.title=element_text(size = 30), legend.text.align = 0, strip.text = element_text(size = 30), panel.grid.major = element_blank(),
        panel.background = element_blank(), panel.grid.minor = element_blank())
dev.off()