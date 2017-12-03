set.seed(150)
## real data analysis
setwd("~/comp-est-R")
source("tune_proximal_gradient.R")
source("naive_methods.R")
source("basic_functions.R")
require(foreach)
require(doParallel)

# load count data
d <- as.matrix(read.table("/Users/yuanpeicao/comp-est-R/dataset/BMI_and_Counts.txt", sep = " ", header = TRUE))
w <- d[, 3:89]
w_b <- w > 0
freq <- colSums(w_b)
sumCount <- colSums(w)

# filter out most rare species
w_select <- w[, -which(freq < 10)]
n <- nrow(w_select)
p <- ncol(w_select)
X_select <- w_select/(rowSums(w_select)%*%matrix(1,1,p))

# test parameter
lambda_opt <- 0.3
c_opt <- 1e-3
alphaX_opt <- min(min(X_select[X_select>0])*c_opt*p, 1)
X_pg_res <- proxGradient(W = w_select, lambda = lambda_opt, alphaX = alphaX_opt,
                         betaX = p, L = 1e-4, gamma = 5, iterMax = 1e+4,
                         epsilon = 1e-5)
X_pg_res_hat <- X_pg_res$XHat

write.csv(X_pg_res_hat, file = '/Users/yuanpeicao/comp-est-R/real-data-result/Xhat.csv')
write.csv(w_select, file = '/Users/yuanpeicao/comp-est-R/real-data-result/w_select.csv')

# loss function: composition estimation from zero with original composition from non-zero (column average comparision) (remove common species)
# original composition from non-zero
x_test <- X_select
freq2 <- colSums(w_select > 0)
x_test <- x_test[, -which(freq2 == nrow(w_select))]
avg_x_test <- colSums(x_test) / colSums(x_test > 0)

# composition estimation from zero index
x_train <- X_pg_res_hat
x_train <- x_train[, -which(freq2 == nrow(w_select))]
x_train[x_test > 0] <- 0
avg_x_train <- colSums(x_train) / colSums(x_train > 0)

# RMSE
mse <- mean((avg_x_test - avg_x_train)^2)
print(mse)
# MAE
mae <- mean(abs(avg_x_test - avg_x_train))
print(mae)