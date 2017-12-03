set.seed(150)

# set.seed(200)

## real data analysis
setwd("~/comp-est-R")
source("tune_proximal_gradient.R")
source("naive_methods.R")
require(foreach)
require(doParallel)

# load count data
d <- as.matrix(read.table("/Users/yuanpeicao/comp-est-R/dataset/BMI_and_Counts.txt", sep = " ", header = TRUE))
w <- d[, 3:89]

# filter out rare species
w_b <- w > 0
freq <- colSums(w_b)
sumCount <- colSums(w)
w_select <- w[, -which(freq < 10 | sumCount < 100)]
n <- nrow(w_select)
p <- ncol(w_select)

# x_res <- autoTuneProxGradient(W = w_select, n_grid = 3)

## Tune parameter
# Split the dataset into training and test set
trainTestData <- split_train_test(w_select, fracTest = 1/4)
w_temp <- trainTestData$W_train
ind_test_row <- trainTestData$ind_test_row
ind_test_col_mat <- trainTestData$ind_test_col_mat
X <- w_temp/(rowSums(w_temp)%*%matrix(1,1,p))

## prepare for the test composition
X_select <- w_select/(rowSums(w_select)%*%matrix(1,1,p))
X_test <- matrix(0, n, p)
for (i in 1:length(ind_test_col_mat)){
  ind_test_col <- unlist(ind_test_col_mat[i], use.names=FALSE)
  X_test[ind_test_row[i], ind_test_col] <- X_select[ind_test_row[i], ind_test_col]
}

## tuning parameter (single)
lambda <- 50
c <- 1.5
alphaX <-  min(min(X_select[X_select>0])*c*p, 1)

## proximal gradient method
X_pg_temp_res <- proxGradient(W = w_temp, lambda = lambda, alphaX = alphaX,
                              betaX = p, L = 1e-4, gamma = 2, iterMax = 1e+4,
                              epsilon = 1e-5)
X_hat_temp <- X_pg_temp_res$XHat
rowSumComp <- rowSums(X_hat_temp)
if (any(abs(rowSumComp-1)>1e-10)){
  stop("Bug! Row sum is not 1!")
}
kl <- KL(X_test[ind_test_row,], X_hat_temp[ind_test_row,])
frob_e <- norm(X_test[ind_test_row,] - X_hat_temp[ind_test_row,], 'F')

print(kl)
print(frob_e)

# ## tuning parameter (multiple)
# # set tuning parameter range
# n_grid <- 30
# index <- 1:n_grid
# 
# lambda_min <- 70
# lambda_max <- 99
# lambda_range <- lambda_min + (lambda_max-lambda_min)*((index-1)/(n_grid-1))
# 
# c_min <- 0.2
# c_max <- 10
# c_range <- c_min + (c_max-c_min)*((index-1)/(n_grid-1))
# 
# c_range <- c(rep(1,n_grid))
# 
# # # Parallel computing for all replications
# # threads <- detectCores() - 2
# # cl <- makeCluster(threads)
# # registerDoParallel(cl)
# 
# # foreach(k = 1:(n_grid * n_grid - 1)) %dopar% {
# #   source("tune_proximal_gradient.R")
# #   source("naive_methods.R")
# # 
# #   i <- floor(k / n_grid) + 1
# #   j <- k %% n_grid + 1
# 
# for (i in 1:n_grid){
#   for (j in 1:n_grid){
#       start <- Sys.time()
#       lambda <- lambda_range[i]
#       c <- c_range[j]
#       alphaX <-  min(min(X_select[X_select>0])*c*p, 1)
# 
#       ## proximal gradient method
#       X_pg_temp_res <- proxGradient(W = w_temp, lambda = lambda, alphaX = alphaX,
#                                     betaX = p, L = 1e-4, gamma = 2, iterMax = 1e+4,
#                                     epsilon = 1e-5)
#       X_hat_temp <- X_pg_temp_res$XHat
#       rowSumComp <- rowSums(X_hat_temp)
#       if (any(abs(rowSumComp-1)>1e-10)){
#         stop("Bug! Row sum is not 1!")
#       }
#       kl <- KL(X_test[ind_test_row,], X_hat_temp[ind_test_row,])
#       frob_e <- norm(X_test[ind_test_row,] - X_hat_temp[ind_test_row,], 'F')
# 
#       end <- Sys.time()
#       run_time <- (end - start)
# 
#       # # path for saving the result
#       # pathRealRes <- paste("~/comp-est-R/realDataAnalysis/lambda", toString(lambda), "-alphaC", toString(c), ".RData", sep = "")
#       # save(kl, frob_e, run_time, file = pathRealRes)
# 
#       print("lambda = ")
#       print(lambda)
#       print("kl = ")
#       print(kl)
#       print("frob error = ")
#       print(frob_e)
#       end <- Sys.time()
#       print("run time = ")
#       print(end - start)
#   }
# }
# # stopCluster(cl)