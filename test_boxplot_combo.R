require(ggplot2)
require(reshape)
require(gridExtra)
require(grid)
require(cowplot)
require(latex2exp)
setwd("~/comp-est-R")

# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.1_c_0.25_p42.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.5_c_0.1_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.1_c_0.5_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/lambda_0.5_c_0.5_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_316_c_1e_100_p45.dat', header = FALSE))
# xHatRare <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_09_c_1e_5_p42.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_01_c_0_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_01_c_05zr_p45.dat', header = FALSE))
xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_1_c_05zr_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_05_c_05zr_p45.dat', header = FALSE))
# xHat <- as.matrix(read.table('/Users/yuanpeicao/comp-est-R/dataset/l_0_075_c_05zr_p45.dat', header = FALSE))

# ## replace zero with the row minimum * 0.5
# zero_threshold <- 1e-4
# xHat_tmp <- xHat
# xHat_tmp[xHat_tmp < zero_threshold] <-1e+10
# xHat_min <- apply(xHat_tmp,1,min)
# xHat_min_matrix <- xHat_min %*% matrix(c(rep(1,ncol(xHat))),1,ncol(xHat)) * 0.5
# xHat[xHat < zero_threshold] <- xHat_min_matrix[xHat < zero_threshold]
# xHat <- xHat/(rowSums(xHat)%*%matrix(1,1,ncol(xHat)))

## Load the composition estimate information
# xHat <- as.matrix(read.csv('/Users/yuanpeicao/comp-est-R/real-data-result/Xhat.csv'))
# xHat <- xHat[, 2:ncol(xHat)]

wSelect <- as.matrix(read.csv('/Users/yuanpeicao/comp-est-R/real-data-result/w_select.csv'))
wSelect <- wSelect[, 2:ncol(wSelect)]

## re-correct the composition for the species whose original counts are 0s
xHat_tmp <- xHat
xHat_tmp[wSelect == 0] <- 100
rowMins <- matrix(apply(xHat_tmp, 1, min), nrow(xHat_tmp), 1)
replaceMatrix <- rowMins %*% matrix(1,1,ncol(xHat_tmp))
xHat[wSelect == 0] <- pmax(replaceMatrix[wSelect == 0] * 0.5, xHat[wSelect == 0])
xHat <- xHat/(rowSums(xHat)%*%matrix(1,1,ncol(xHat)))

## filter out common species
w_b <- wSelect > 0
freq <- colSums(w_b)

xHatRare <- xHat[, -which(freq == nrow(wSelect))]
wSelectRare <- wSelect[, -which(freq == nrow(wSelect))]

# xHatRare <- xHat
# wSelectRare <- wSelect

## Load the species name
speciesNames <- colnames(wSelectRare)
taxa <- matrix(unlist(strsplit(speciesNames, "\\.")), 6)
phyla <- taxa[2, ]
genera <- taxa[6, ]

## Split the compsition into rare and common
xZeroIndex <- (wSelectRare == 0)

n <- nrow(xHatRare)
p <- ncol(xHatRare)

## proximal gradient result
dat <- data.frame(Genera = c(rep(genera,each = n)),zero_index = factor(c(xZeroIndex)), Composition = c(log(xHatRare)))
dat$zero_index <- factor(dat$zero_index, labels = c(expression("italic(M)^c","italic(M)")))

## boxplot
# setEPS(width = 18, height = 10)
# postscript(file="boxplotSR_test.eps")
ggplot(dat, aes(x = Genera, y = Composition), mar = c(0.5,0.5,0.5,0.5)) + 
  geom_boxplot(aes(fill = zero_index), coef = 10) +
  ylab(expression(log(widehat(X)))) + xlab("Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1), legend.title = element_blank(), legend.text = element_text(size = 30)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"),
                    labels = c(expression(Omega),expression(Omega^c))) +
  theme(axis.title=element_text(size = 30), legend.text.align = 0, strip.text = element_text(size = 30), panel.grid.major = element_blank(),
        panel.background = element_blank(), panel.grid.minor = element_blank())
# dev.off()