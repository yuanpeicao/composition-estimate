setwd("~/comp-est-R")

source("tune_proximal_gradient.R")
source("naive_methods.R")
source("real_data_analysis_plot.R")
require(ggplot2)

## prepare the estimated composition
d <- as.matrix(read.table("/Users/yuanpeicao/comp-est-R/dataset/BMI_and_Counts.txt", sep = " ", header = TRUE))
w <- d[, 3:89]

w_b <- w > 0
freq <- colSums(w_b)
sumCount <- colSums(w)

w_common <- w[, which(freq == 98)]
w_common_sum <- rowSums(w_common)
w_common_avg <- matrix(1, nrow(w_common), 1) %*% colSums(w_common)/nrow(w_common)
w_common[w_common == 0] <- w_common_avg[w_common == 0]
x_common <- w_common/(rowSums(w_common)%*%matrix(1,1,ncol(w_common)))

w_rare <- w[,-which(freq < 10 | sumCount < 100 | freq == 98)]
w_rare_sum <- rowSums(w_rare)

percentage_common <- w_common_sum / (w_common_sum + w_rare_sum)
percentage_rare <- w_rare_sum / (w_common_sum + w_rare_sum)

## estimate composition using selected parameter
p <- ncol(w_rare)
X_rare <- w_rare/(rowSums(w_rare)%*%matrix(1,1,p))
lambda <- 1.5
c <- 1
alphaX <- min(min(X_rare[X_rare>0])*c*p, 1)
X_pg_res <- proxGradient(W = w_rare, lambda = lambda, alphaX = alphaX,
                         betaX = p, L = 1e-4, gamma = 2, iterMax = 1e+4,
                         epsilon = 1e-5)

## combine two estimation together
x_common_hat <- percentage_common %*%matrix(1,1,ncol(x_common))
x_rare_hat <- percentage_rare %*%matrix(1,1,ncol(X_pg_res$XHat))
xHat <- cbind(x_common_hat, x_rare_hat)

w_select <- cbind(w[, which(freq == 98)], w[, -which(freq < 10 | sumCount < 100 | freq == 98)])

#### shannon scatter plot
X_select <- w_select/(rowSums(w_select)%*%matrix(1,1,ncol(w_select)))
comp_shannon_scatterplot(XHat1 = X_select, XHat2 = xHat,
                         NameEst1 = "zero_replacement", NameEst2 = "proximal_gradient")

comp_simpson_scatterplot(XHat1 = X_select, XHat2 = xHat,
                         NameEst1 = "zero_replacement", NameEst2 = "proximal_gradient")



## Load the species name
speciesNames <- colnames(w_select)
taxa <- matrix(unlist(strsplit(speciesNames, "\\.")), 6)
phyla <- taxa[2, ]
genera <- taxa[6, ]

## Split the compsition into zero and non-zero
xZeroIndex <- (w_select == 0)

n <- nrow(xHat)
p <- ncol(xHat)

dat <- data.frame(Genera = c(rep(genera,each = n)),zero_index = factor(c(xZeroIndex)), Composition = c(xHat))
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

## zero-replacement


setEPS(width = 18, height = 10)
postscript(file="boxplotHat.eps")
ggplot(log_dat, aes(x = Genera, y = Composition), mar = c(0.5,0.5,0.5,0.5)) + ylim(-10, 0) +
  geom_boxplot(aes(fill = zero_index), coef = 10) +
  ylab(expression(log(widehat(X)))) + xlab("Genera") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1), legend.title = element_blank(), legend.text = element_text(size = 30)) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9"),
                    labels = c(expression(Omega),expression(Omega^c))) +
  theme(axis.title=element_text(size = 30), legend.text.align = 0, strip.text = element_text(size = 30), panel.grid.major = element_blank(),
        panel.background = element_blank(), panel.grid.minor = element_blank())
dev.off()