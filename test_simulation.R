# test simulation 
setwd("~/comp-est-R")
source("simulation.R")

n <- 100
p <- 50
para <- 20
gamma <- 3

N <- floor(gamma*n*p)

X0 <- generateLowRankComp(n, p, para)
X0_v <- as.vector(X0)

R0 <- generateR(n)

N_row <- floor(N*R0)

W <- t(sapply(seq(n), function(i){rmultinom(1, N_row[i], X0[i, ])}))

Xhat_pg_res <- autoTuneProxGradient(W = W, n_grid = 10)

xhat <- Xhat_pg_res$X_hat
xhat_v <- as.vector(xhat)

plot(X0_v, xhat_v)

X0_zr <- W
X0_zr[X0_zr == 0] <- 0.5
X0_zr <- X0_zr/(rowSums(X0_zr)%*%matrix(1,1,ncol(X0_zr)))
X0_zr_v <- as.vector(X0_zr)

plot(X0_v, X0_zr_v)
