setwd("~/comp-est-R")

require(MASS)
require(foreach)
require(doParallel)
source("proximal_gradient.R")
source("tune_proximal_gradient.R")
source("naive_methods.R")
#----------------------------------------------------------------------------------------
#  Generate replications given X0 and R0
#  Input:
#                X0 ------ n x p true composition matrix (row/column is sample/variable)
#                R0 ------ n x 1 row probabiltiy vector
#              nRep ------ number of replications
#             gamma ------ total count fraction
#       settingName ------ "lowRank" or "appRank"
#       containZero ------ contain zero if containZero = 1; else, X0 does not contain zero 
#  Output:
#             X_hat ------ estimated zero-replacement composition
#-----------------------------------------------------------------------------------------
generateRep <- function(X0, R0, nRep = 50, gamma = 1, settingName = "lowRank", containZero = 0){
  set.seed(1)
  n <- nrow(X0)
  p <- ncol(X0)
  N <- floor(gamma*n*p)
  
  if (containZero == 0){
    for (i in 1:nRep){
      N_row <- floor(N*R0)
      W <- t(sapply(seq(n), function(i){rmultinom(1, N_row[i], X0[i, ])}))
      pathName <- paste("~/comp-est/simulation/nonzero/p", toString(p), "/", settingName, "N", toString(gamma), "/replicate/rep", toString(i), ".RData", sep = "")
      save(W, file = pathName)
      }
  }else if (containZero == 1){
    for (i in 1:nRep){
      N_row <- floor(N*R0)
      W <- t(sapply(seq(n), function(i){rmultinom(1, N_row[i], X0[i, ])}))
      pathName <- paste("~/comp-est/simulation/zero/p", toString(p), "/", settingName, "N", toString(gamma), "/replicate/rep", toString(i), ".RData", sep = "")
      save(W, file = pathName)
    }
  }
}

#----------------------------------------------------------------------------------------
#  Generate row probability vector R0 from unifrom [a, b]
#  Input:
#                 n ------ number of rows
#                 a ------ low bound
#                 b ------ upper bound
#                
#  Output:
#                R0 ------ n x 1 row probabiltiy vector
#-----------------------------------------------------------------------------------------
generateR <- function(n, a = 1, b = 10){
  R0 <- runif(n, min = a, max = b)
  R0 <- R0/sum(R0)
  return(R0 = R0)
}

#----------------------------------------------------------------------------------------
#  Generate n x p low rank composition matrix X0
#  Input:
#                 n ------ number of rows
#                 p ------ number of columns
#                 r ------ rank of matrix
#                
#  Output:
#                X0 ------ n x p composition matrix
#-----------------------------------------------------------------------------------------
generateLowRankComp <- function(n, p, r){
  while (TRUE){
    U <- abs(matrix(rnorm(n*r, mean = 0, sd = 1), n, r))
    V <- matrix((runif(p*p) < 0.3) + 0, p, p)
    V <- V - diag(diag(V)) + diag(p)
    V <- V[,1:r] + matrix(rnorm(p*r, mean = 0, sd = 1e-3), p, r)
    Z <- U %*% t(V)
    if (sum(Z <= 0) == 0){
      X0 <- Z/(rowSums(Z)%*%matrix(1,1,p))
      return(X0 = X0)
    }
  }
}

#----------------------------------------------------------------------------------------
#  Generate n x p low rank composition matrix X0 containing zero
#  Input:
#                 n ------ number of rows
#                 p ------ number of columns
#                 r ------ rank of matrix
#          zeroProp ------ the proportion of zeros
#                
#  Output:
#                X0 ------ n x p composition matrix
#-----------------------------------------------------------------------------------------
generateLowRankCompWZero <- function(n, p, r, zeroProp = 0.1){
  X0 = generateLowRankComp(n, p, r)
  while (TRUE){
    for (i in 1:p){
      zeroRowInd <- sample(1:n, floor(n*zeroProp))
      X0[zeroRowInd, i] = 0
    }
    # check whether whole row is zero
    if (sum(rowSums(X0) == 0) == 0){
      X0 <- X0/(rowSums(X0)%*%matrix(1,1,p))
      return(X0 = X0)
    }
  }
}

#----------------------------------------------------------------------------------------
#  Generate n x p approximately low rank composition matrix X0
#  Input:
#                 n ------ number of rows
#                 p ------ number of columns
#                 a ------ approximately low rank parameter
#                
#  Output:
#                X0 ------ n x p composition matrix
#-----------------------------------------------------------------------------------------
generateAppLowRankComp <- function(n, p, a){
  while (TRUE){
    U <- abs(matrix(rnorm(n*n, mean = 0, sd = 1), n, n))
    V <- matrix((runif(p*p) < 0.3) + 0, p, p)
    V <- V - diag(diag(V)) + diag(p) + matrix(rnorm(p*p, mean = 0, sd = 1e-3), p, p)
    D <- matrix(0, n, p)
    diag(D) <- (1:min(n,p))^(-a)
    Z <- U %*% D %*% t(V)
    if (sum(Z <= 0) == 0){
      X0 <- Z/(rowSums(Z)%*%matrix(1,1,p))
      return(X0 = X0)
    }
  }
}

#----------------------------------------------------------------------------------------
#  Generate n x p approximately low rank composition matrix X0 containing zero
#  Input:
#                 n ------ number of rows
#                 p ------ number of columns
#                 a ------ approximately low rank parameter
#          zeroProp ------ the proportion of zeros
#                
#  Output:
#                X0 ------ n x p composition matrix
#-----------------------------------------------------------------------------------------
generateAppLowRankCompWZero <- function(n, p, a, zeroProp = 0.1){
  X0 = generateAppLowRankComp(n, p, a)
  while (TRUE){
    for (i in 1:p){
      zeroRowInd <- sample(1:n, floor(n*zeroProp))
      X0[zeroRowInd, i] = 0
    }
    # check whether whole row is zero
    if (sum(rowSums(X0) == 0) == 0){
      X0 <- X0/(rowSums(X0)%*%matrix(1,1,p))
      return(X0 = X0)
    }
  }
}

#----------------------------------------------------------------------------------------
#  Generate replications main function
#  Input:
#                 n ------ number of rows
#                 p ------ number of columns
#           lowRank ------ 1 if low rank matrix; 0 if approxmiately low rank matrix
#              para ------ the rank if lowRank == 1; apprxomiatly low rank parameter if lowRank == 0
#       containZero ------ contain zero if containZero = 1; else, X0 does not contain zero 
#              nRep ------ number of replications
#      
#             gamma ------ total count fraction 
#-----------------------------------------------------------------------------------------
mainGenerate <- function(n, p, lowRank, para, containZero = 0, nRep = 10, gamma = 1){
  R0 <- generateR(n)
  if ((lowRank == 1) && (containZero == 0)){
    X0 <- generateLowRankComp(n, p, para)
    pathTrueComp <- paste("~/comp-est/simulation/nonzero/p", toString(p), "/lowRankN", toString(gamma), "/trueComp.RData", sep = "")
    save(X0, file = pathTrueComp)
    generateRep(X0, R0, nRep, gamma, settingName = "lowRank", containZero = 0)
  }else if((lowRank == 1) && (containZero == 1)){
    X0 <- generateLowRankComp(n, p, para)
    pathTrueComp <- paste("~/comp-est/simulation/zero/p", toString(p), "/lowRankN", toString(gamma), "/trueComp.RData", sep = "")
    save(X0, file = pathTrueComp)
    generateRep(X0, R0, nRep, gamma, settingName = "lowRank", containZero = 1)
  }else if((lowRank == 0) && (containZero == 0)){
    X0 <- generateAppLowRankComp(n, p, para)
    pathTrueComp <- paste("~/comp-est/simulation/nonzero/p", toString(p), "/appRankN", toString(gamma), "/trueComp.RData", sep = "")
    save(X0, file = pathTrueComp)
    generateRep(X0, R0, nRep, gamma, settingName = "appRank", containZero = 0)
  }else if((lowRank == 0) && (containZero == 1)){
    X0 <- generateAppLowRankComp(n, p, para)
    pathTrueComp <- paste("~/comp-est/simulation/zero/p", toString(p), "/appRankN", toString(gamma), "/trueComp.RData", sep = "")
    save(X0, file = pathTrueComp)
    generateRep(X0, R0, nRep, gamma, settingName = "appRank", containZero = 1)
  }
}

#----------------------------------------------------------------------------------------
#  main function for composition estimation
#  Input:
#              nRep ------ number of replications
#                 p ------ number of columns
#       containZero ------ contain zero if containZero = 1; else, X0 does not contain zero 
#           lowRank ------ 1 if low rank matrix; 0 if approxmiately low rank matrix
#                 r ------ the true rank in the low-rank setting
#             gamma ------ total count fraction 
#-----------------------------------------------------------------------------------------
mainEst <- function(nRep = 50, p = 50, containZero = 0, lowRank = 1, r = 0, gamma = 1){
  
  start <- Sys.time()
  
  if ((lowRank == 1) && (containZero == 0)){
    pathCompCount = paste("~/comp-est/simulation/nonzero/p", toString(p), "/lowRankN", toString(gamma), "/", sep = "")
  }else if((lowRank == 1) && (containZero == 1)){
    pathCompCount = paste("~/comp-est/simulation/zero/p", toString(p), "/lowRankN", toString(gamma), "/", sep = "")
  }else if((lowRank == 0) && (containZero == 0)){
    pathCompCount = paste("~/comp-est/simulation/nonzero/p", toString(p), "/appRankN", toString(gamma), "/", sep = "")
  }else if((lowRank == 0) && (containZero == 1)){
    pathCompCount = paste("~/comp-est/simulation/zero/p", toString(p), "/appRankN", toString(gamma), "/", sep = "")
  }
  
  pathTrueComp <- paste(pathCompCount, "trueComp.RData", sep = "")
  load(pathTrueComp)
  
  # diversity index
  d0 <- diversity(X0)
  n <- nrow(X0)
  
  # Parallel computing for all replications
  threads <- detectCores() - 1
  cl <- makeCluster(threads)
  registerDoParallel(cl)
  
  foreach(i = 1:nRep) %dopar% {
    source("proximal_gradient.R")
    source("tune_proximal_gradient.R")
    source("naive_methods.R")
    
    pathCount <- paste(pathCompCount, "replicate/rep", toString(i), ".RData", sep = "")
    load(pathCount)
    
    # Proximal Gradient estimator
    Xhat_pg_res <- autoTuneProxGradient(W = W, n_grid = 3)
    Xhat_pg <- Xhat_pg_res$X_hat
    
    # Zero-replacement
    Xhat_zr05 <- zr(W = W, alpha = 0.5) # replaced by 0.5
    Xhat_zrCv_res <- autoTuneZr(W = W, fracTest = 1/4)
    Xhat_zrCv <- Xhat_zrCv_res$X_hat
    alpha_zr_oracle <- tune_zr_oracle(W = W, X0 = X0, n_grid = 1e+2)
    Xhat_zrOralce <- alpha_zr_oracle$X_hat # replaced by oracle parameter
    
    ## Frobeius norm loss
    # Proximal Gradient estimator
    frob_e_pg <- norm(Xhat_pg - X0, 'F')
    # Zero-replacement estimator
    frob_e_zr05 <- norm(Xhat_zr05 - X0, 'F')
    frob_e_zrCv <- norm(Xhat_zrCv - X0, 'F')
    frob_e_zrOracle <- norm(Xhat_zrOralce - X0, 'F')
    
    ## KL divergence
    # Proximal Gradient estimator
    kl_e_pg <- KL(X0, Xhat_pg)
    # Zero-replacement estimator
    kl_e_zr05 <- KL(X0, Xhat_zr05)
    kl_e_zrCv <- KL(X0, Xhat_zrCv)
    kl_e_zrOracle <- KL(X0, Xhat_zrOralce)
    
    ## error for diversity index
    # Proximal Gradient estimator
    d_pg <- diversity(Xhat_pg)
    # Zero-replacement estimator
    d_zr05 <- diversity(Xhat_zr05)
    d_zrCv <- diversity(Xhat_zrCv)
    d_zrOracle <- diversity(Xhat_zrOralce)
    
    ## error for Shannon index
    # Proximal Gradient estimator
    sse_sh_pg <- (norm(t(d_pg$sh - d0$sh), 'F'))^2/n
    # Zero-replacement estimator
    sse_sh_zr05 <- (norm(t(d_zr05$sh - d0$sh), 'F'))^2/n
    sse_sh_zrCv <- (norm(t(d_zrCv$sh - d0$sh), 'F'))^2/n
    sse_sh_zrOracle <- (norm(t(d_zrOracle$sh - d0$sh), 'F'))^2/n
    
    ## error for Simpson index
    # Proximal Gradient estimator
    sse_sp_pg <- (norm(t(d_pg$sp - d0$sp), 'F'))^2/n
    # Zero-replacement estimator
    sse_sp_zr05 <- (norm(t(d_zr05$sp - d0$sp), 'F'))^2/n
    sse_sp_zrCv <- (norm(t(d_zrCv$sp - d0$sp), 'F'))^2/n
    sse_sp_zrOracle <- (norm(t(d_zrOracle$sp - d0$sp), 'F'))^2/n
    
    ## error for Bray Curtis index
    # Proximal Gradient estimator
    sse_bc_pg <- (norm(t(d_pg$bc - d0$bc), 'F'))^2/(n*(n-1))
    # Zero-replacement estimator
    sse_bc_zr05 <- (norm(t(d_zr05$bc - d0$bc), 'F'))^2/(n*(n-1))
    sse_bc_zrCv <- (norm(t(d_zrCv$bc - d0$bc), 'F'))^2/(n*(n-1))
    sse_bc_zrOracle <- (norm(t(d_zrOracle$bc - d0$bc), 'F'))^2/(n*(n-1))
    
    # SVT
    if (r > 0){
      Xhat_svt05 <- svt(W = W, r = r, alpha = 0.5) # replaced by 0.5
      Xhat_svtCv_res <- autoTuneSvt(W = W, r = r, fracTest = 1/4)
      Xhat_svtCv <- Xhat_svtCv_res$X_hat
      alpha_svt_oracle <- tune_svt_oracle(W = W, r = r, X0 = X0, n_grid = 1e+2)
      Xhat_svtOralce <- alpha_svt_oracle$X_hat # replaced by oracle parameter
      
      ## Frobeius norm loss
      frob_e_svt05 <- norm(Xhat_svt05 - X0, 'F')
      frob_e_svtCv <- norm(Xhat_svtCv - X0, 'F')
      frob_e_svtOracle <- norm(Xhat_svtOralce - X0, 'F')
      
      ## KL divergence
      kl_e_svt05 <- KL(X0, Xhat_svt05)
      kl_e_svtCv <- KL(X0, Xhat_svtCv)
      kl_e_svtOracle <- KL(X0, Xhat_svtOralce)
      
      ## error for diversity index
      d_svt05 <- diversity(Xhat_svt05)
      d_svtCv <- diversity(Xhat_svtCv)
      d_svtOracle <- diversity(Xhat_svtOralce)
      
      ## error for Shannon index
      sse_sh_svt05 <- (norm(t(d_svt05$sh - d0$sh), 'F'))^2/n
      sse_sh_svtCv <- (norm(t(d_svtCv$sh - d0$sh), 'F'))^2/n
      sse_sh_svtOracle <- (norm(t(d_svtOracle$sh - d0$sh), 'F'))^2/n
      
      ## error for Simpson index
      sse_sp_svt05 <- (norm(t(d_svt05$sp - d0$sp), 'F'))^2/n
      sse_sp_svtCv <- (norm(t(d_svtCv$sp - d0$sp), 'F'))^2/n
      sse_sp_svtOracle <- (norm(t(d_svtOracle$sp - d0$sp), 'F'))^2/n
      
      ## error for Bray Curtis index
      sse_bc_svt05 <- (norm(t(d_svt05$bc - d0$bc), 'F'))^2/(n*(n-1))
      sse_bc_svtCv <- (norm(t(d_svtCv$bc - d0$bc), 'F'))^2/(n*(n-1))
      sse_bc_svtOracle <- (norm(t(d_svtOracle$bc - d0$bc), 'F'))^2/(n*(n-1))
    }
    
    
    # path for saving the result
    pathSimRes <- paste(pathCompCount, "result/rep", toString(i), ".RData", sep = "")
    if (r == 0){
      save(frob_e_pg, frob_e_zr05, frob_e_zrCv, frob_e_zrOracle,
           kl_e_pg, kl_e_zr05, kl_e_zrCv, kl_e_zrOracle,
           sse_sh_pg, sse_sh_zr05, sse_sh_zrCv, sse_sh_zrOracle,
           sse_sp_pg, sse_sp_zr05, sse_sp_zrCv, sse_sp_zrOracle,
           sse_bc_pg, sse_bc_zr05, sse_bc_zrCv, sse_bc_zrOracle,
           file = pathSimRes)
    }else if (r > 0){
      save(frob_e_pg, frob_e_zr05, frob_e_zrCv, frob_e_zrOracle, frob_e_svt05, frob_e_svtCv, frob_e_svtOracle,
           kl_e_pg, kl_e_zr05, kl_e_zrCv, kl_e_zrOracle, kl_e_svt05, kl_e_svtCv, kl_e_svtOracle,
           sse_sh_pg, sse_sh_zr05, sse_sh_zrCv, sse_sh_zrOracle, sse_sh_svt05, sse_sh_svtCv, sse_sh_svtOracle,
           sse_sp_pg, sse_sp_zr05, sse_sp_zrCv, sse_sp_zrOracle, sse_sp_svt05, sse_sp_svtCv, sse_sp_svtOracle,
           sse_bc_pg, sse_bc_zr05, sse_bc_zrCv, sse_bc_zrOracle, sse_bc_svt05, sse_bc_svtCv, sse_bc_svtOracle,
           file = pathSimRes)
    }
  }
  
  stopCluster(cl)
  end <- Sys.time()
  print("Main function ends")
  print(end-start)
}

#----------------------------------------------------------------------------------------
#  Record the summary of composition estimation
#  Input:
#              nRep ------ number of replications
#                 p ------ number of columns
#       containZero ------ contain zero if containZero = 1; else, X0 does not contain zero 
#           lowRank ------ 1 if low rank matrix; 0 if approxmiately low rank matrix
#                 r ------ the true rank in the low-rank setting
#             gamma ------ total count fraction 
#-----------------------------------------------------------------------------------------
mainEstSummary <- function(nRep = 50, p = 50, containZero = 0, lowRank = 1, r = 0, gamma = 1){
  
  start <- Sys.time()

  if ((lowRank == 1) && (containZero == 0)){
    pathCompCount = paste("~/comp-est/simulation/nonzero/p", toString(p), "/lowRankN", toString(gamma), "/", sep = "")
  }else if((lowRank == 1) && (containZero == 1)){
    pathCompCount = paste("~/comp-est/simulation/zero/p", toString(p), "/lowRankN", toString(gamma), "/", sep = "")
  }else if((lowRank == 0) && (containZero == 0)){
    pathCompCount = paste("~/comp-est/simulation/nonzero/p", toString(p), "/appRankN", toString(gamma), "/", sep = "")
  }else if((lowRank == 0) && (containZero == 1)){
    pathCompCount = paste("~/comp-est/simulation/zero/p", toString(p), "/appRankN", toString(gamma), "/", sep = "")
  }
  
  folderCompCount = paste(pathCompCount, "result", sep = "")
  
  frob_e_pg_v <- rep(0, nRep)
  frob_e_zr05_v <- rep(0, nRep)
  frob_e_zrCv_v <- rep(0, nRep)
  frob_e_zrOracle_v <- rep(0, nRep)
  
  kl_e_pg_v <- rep(0, nRep)
  kl_e_zr05_v <- rep(0, nRep)
  kl_e_zrCv_v <- rep(0, nRep)
  kl_e_zrOracle_v <- rep(0, nRep)
  
  sse_sh_pg_v <- rep(0, nRep)
  sse_sh_zr05_v <- rep(0, nRep)
  sse_sh_zrCv_v <- rep(0, nRep)
  sse_sh_zrOracle_v <- rep(0, nRep)
  
  sse_sp_pg_v <- rep(0, nRep)
  sse_sp_zr05_v <- rep(0, nRep)
  sse_sp_zrCv_v <- rep(0, nRep)
  sse_sp_zrOracle_v <- rep(0, nRep)
  
  sse_bc_pg_v <- rep(0, nRep)
  sse_bc_zr05_v <- rep(0, nRep)
  sse_bc_zrCv_v <- rep(0, nRep)
  sse_bc_zrOracle_v <- rep(0, nRep)
  
  if(r > 0){
    frob_e_svt05_v <- rep(0, nRep)
    frob_e_svtCv_v <- rep(0, nRep)
    frob_e_svtOracle_v <- rep(0, nRep)
    
    kl_e_svt05_v <- rep(0, nRep)
    kl_e_svtCv_v <- rep(0, nRep)
    kl_e_svtOracle_v <- rep(0, nRep)
    
    sse_sh_svt05_v <- rep(0, nRep)
    sse_sh_svtCv_v <- rep(0, nRep)
    sse_sh_svtOracle_v <- rep(0, nRep)
    
    sse_sp_svt05_v <- rep(0, nRep)
    sse_sp_svtCv_v <- rep(0, nRep)
    sse_sp_svtOracle_v <- rep(0, nRep)
    
    sse_bc_svt05_v <- rep(0, nRep)
    sse_bc_svtCv_v <- rep(0, nRep)
    sse_bc_svtOracle_v <- rep(0, nRep)
  }
  
  
  for (i in 1:nRep){
    pathCompEstRes <- paste(folderCompCount, "/rep", toString(i), ".RData", sep = "")
    load(pathCompEstRes)
    
    frob_e_pg_v[i] <- frob_e_pg
    frob_e_zr05_v[i] <- frob_e_zr05
    frob_e_zrCv_v[i] <- frob_e_zrCv
    frob_e_zrOracle_v[i] <- frob_e_zrOracle
    
    kl_e_pg_v[i] <- kl_e_pg
    kl_e_zr05_v[i] <- kl_e_zr05
    kl_e_zrCv_v[i] <- kl_e_zrCv
    kl_e_zrOracle_v[i] <- kl_e_zrOracle
    
    sse_sh_pg_v[i] <- sse_sh_pg
    sse_sh_zr05_v[i] <- sse_sh_zr05
    sse_sh_zrCv_v[i] <- sse_sh_zrCv
    sse_sh_zrOracle_v[i] <- sse_sh_zrOracle
    
    sse_sp_pg_v[i] <- sse_sp_pg
    sse_sp_zr05_v[i] <- sse_sp_zr05
    sse_sp_zrCv_v[i] <- sse_sp_zrCv
    sse_sp_zrOracle_v[i] <- sse_sp_zrOracle
    
    sse_bc_pg_v[i] <- sse_bc_pg
    sse_bc_zr05_v[i] <- sse_bc_zr05
    sse_bc_zrCv_v[i] <- sse_bc_zrCv
    sse_bc_zrOracle_v[i] <- sse_bc_zrOracle
    
    if(r > 0){
      frob_e_svt05_v[i] <- frob_e_svt05
      frob_e_svtCv_v[i] <- frob_e_svtCv
      frob_e_svtOracle_v[i] <- frob_e_svtOracle
      
      kl_e_svt05_v[i] <- kl_e_svt05
      kl_e_svtCv_v[i] <- kl_e_svtCv
      kl_e_svtOracle_v[i] <- kl_e_svtOracle
      
      sse_sh_svt05_v[i] <- sse_sh_svt05
      sse_sh_svtCv_v[i] <- sse_sh_svtCv
      sse_sh_svtOracle_v[i] <- sse_sh_svtOracle
      
      sse_sp_svt05_v[i] <- sse_sp_svt05
      sse_sp_svtCv_v[i] <- sse_sp_svtCv
      sse_sp_svtOracle_v[i] <- sse_sp_svtOracle
      
      sse_bc_svt05_v[i] <- sse_bc_svt05
      sse_bc_svtCv_v[i] <- sse_bc_svtCv
      sse_bc_svtOracle_v[i] <- sse_bc_svtOracle
    }
  }
  
  frob_e_pg_mean <- mean(frob_e_pg_v)
  frob_e_zr05_mean <- mean(frob_e_zr05_v)
  frob_e_zrCv_mean <- mean(frob_e_zrCv_v)
  frob_e_zrOracle_mean <- mean(frob_e_zrOracle_v)
  
  kl_e_pg_mean <- mean(kl_e_pg_v)
  kl_e_zr05_mean <- mean(kl_e_zr05_v)
  kl_e_zrCv_mean <- mean(kl_e_zrCv_v)
  kl_e_zrOracle_mean <- mean(kl_e_zrOracle_v)
  
  sse_sh_pg_mean <- mean(sse_sh_pg_v)
  sse_sh_zr05_mean <- mean(sse_sh_zr05_v)
  sse_sh_zrCv_mean <- mean(sse_sh_zrCv_v)
  sse_sh_zrOracle_mean <- mean(sse_sh_zrOracle_v)
  
  sse_sp_pg_mean <- mean(sse_sp_pg_v)
  sse_sp_zr05_mean <- mean(sse_sp_zr05_v)
  sse_sp_zrCv_mean <- mean(sse_sp_zrCv_v)
  sse_sp_zrOracle_mean <- mean(sse_sp_zrOracle_v)
  
  sse_bc_pg_mean <- mean(sse_bc_pg_v)
  sse_bc_zr05_mean <- mean(sse_bc_zr05_v)
  sse_bc_zrCv_mean <- mean(sse_bc_zrCv_v)
  sse_bc_zrOracle_mean <- mean(sse_bc_zrOracle_v)
  
  frob_e_pg_se <- sd(frob_e_pg_v)/sqrt(nRep)
  frob_e_zr05_se <- sd(frob_e_zr05_v)/sqrt(nRep)
  frob_e_zrCv_se <- sd(frob_e_zrCv_v)/sqrt(nRep)
  frob_e_zrOracle_se <- sd(frob_e_zrOracle_v)/sqrt(nRep)
  
  kl_e_pg_se <- sd(kl_e_pg_v)/sqrt(nRep)
  kl_e_zr05_se <- sd(kl_e_zr05_v)/sqrt(nRep)
  kl_e_zrCv_se <- sd(kl_e_zrCv_v)/sqrt(nRep)
  kl_e_zrOracle_se <- sd(kl_e_zrOracle_v)/sqrt(nRep)
  
  sse_sh_pg_se <- sd(sse_sh_pg_v)/sqrt(nRep)
  sse_sh_zr05_se <- sd(sse_sh_zr05_v)/sqrt(nRep)
  sse_sh_zrCv_se <- sd(sse_sh_zrCv_v)/sqrt(nRep)
  sse_sh_zrOracle_se <- sd(sse_sh_zrOracle_v)/sqrt(nRep)
  
  sse_sp_pg_se <- sd(sse_sp_pg_v)/sqrt(nRep)
  sse_sp_zr05_se <- sd(sse_sp_zr05_v)/sqrt(nRep)
  sse_sp_zrCv_se <- sd(sse_sp_zrCv_v)/sqrt(nRep)
  sse_sp_zrOracle_se <- sd(sse_sp_zrOracle_v)/sqrt(nRep)
  
  sse_bc_pg_se <- sd(sse_bc_pg_v)/sqrt(nRep)
  sse_bc_zr05_se <- sd(sse_bc_zr05_v)/sqrt(nRep)
  sse_bc_zrCv_se <- sd(sse_bc_zrCv_v)/sqrt(nRep)
  sse_bc_zrOracle_se <- sd(sse_bc_zrOracle_v)/sqrt(nRep)
  
  if(r > 0){
    frob_e_svt05_mean <- mean(frob_e_svt05_v)
    frob_e_svtCv_mean <- mean(frob_e_svtCv_v)
    frob_e_svtOracle_mean <- mean(frob_e_svtOracle_v)
    
    kl_e_svt05_mean <- mean(kl_e_svt05_v)
    kl_e_svtCv_mean <- mean(kl_e_svtCv_v)
    kl_e_svtOracle_mean <- mean(kl_e_svtOracle_v)
    
    sse_sh_svt05_mean <- mean(sse_sh_svt05_v)
    sse_sh_svtCv_mean <- mean(sse_sh_svtCv_v)
    sse_sh_svtOracle_mean <- mean(sse_sh_svtOracle_v)
    
    sse_sp_svt05_mean <- mean(sse_sp_svt05_v)
    sse_sp_svtCv_mean <- mean(sse_sp_svtCv_v)
    sse_sp_svtOracle_mean <- mean(sse_sp_svtOracle_v)
    
    sse_bc_svt05_mean <- mean(sse_bc_svt05_v)
    sse_bc_svtCv_mean <- mean(sse_bc_svtCv_v)
    sse_bc_svtOracle_mean <- mean(sse_bc_svtOracle_v)
    
    frob_e_svt05_se <- sd(frob_e_svt05_v)/sqrt(nRep)
    frob_e_svtCv_se <- sd(frob_e_svtCv_v)/sqrt(nRep)
    frob_e_svtOracle_se <- sd(frob_e_svtOracle_v)/sqrt(nRep)
    
    kl_e_svt05_se <- sd(kl_e_svt05_v)/sqrt(nRep)
    kl_e_svtCv_se <- sd(kl_e_svtCv_v)/sqrt(nRep)
    kl_e_svtOracle_se <- sd(kl_e_svtOracle_v)/sqrt(nRep)
    
    sse_sh_svt05_se <- sd(sse_sh_svt05_v)/sqrt(nRep)
    sse_sh_svtCv_se <- sd(sse_sh_svtCv_v)/sqrt(nRep)
    sse_sh_svtOracle_se <- sd(sse_sh_svtOracle_v)/sqrt(nRep)
    
    sse_sp_svt05_se <- sd(sse_sp_svt05_v)/sqrt(nRep)
    sse_sp_svtCv_se <- sd(sse_sp_svtCv_v)/sqrt(nRep)
    sse_sp_svtOracle_se <- sd(sse_sp_svtOracle_v)/sqrt(nRep)
    
    sse_bc_svt05_se <- sd(sse_bc_svt05_v)/sqrt(nRep)
    sse_bc_svtCv_se <- sd(sse_bc_svtCv_v)/sqrt(nRep)
    sse_bc_svtOracle_se <- sd(sse_bc_svtOracle_v)/sqrt(nRep)
  }
  
  if (r > 0){
    summaryDf <- data.frame(cbind(frob_e_pg_mean = frob_e_pg_mean, frob_e_pg_se = frob_e_pg_se, 
                                  frob_e_zr05_mean = frob_e_zr05_mean, frob_e_zr05_se = frob_e_zr05_se, 
                                  frob_e_zrCv_mean = frob_e_zrCv_mean, frob_e_zrCv_se = frob_e_zrCv_se, 
                                frob_e_zrOracle_mean = frob_e_zrOracle_mean, frob_e_zrOracle_se = frob_e_zrOracle_se, 
                                frob_e_svt05_mean = frob_e_svt05_mean, frob_e_svt05_se =frob_e_svt05_se, 
                                frob_e_svtCv_mean = frob_e_svtCv_mean, frob_e_svtCv_se = frob_e_svtCv_se, 
                                frob_e_svtOracle_mean = frob_e_svtOracle_mean, frob_e_svtOracle_se = frob_e_svtOracle_se, 
                                kl_e_pg_mean = kl_e_pg_mean, kl_e_pg_se = kl_e_pg_se, 
                                kl_e_zr05_mean = kl_e_zr05_mean, kl_e_zr05_se = kl_e_zr05_se, 
                                kl_e_zrCv_mean = kl_e_zrCv_mean, kl_e_zrCv_se = kl_e_zrCv_se, 
                                kl_e_zrOracle_mean = kl_e_zrOracle_mean, kl_e_zrOracle_se = kl_e_zrOracle_se, 
                                kl_e_svt05_mean = kl_e_svt05_mean, kl_e_svt05_se = kl_e_svt05_se, 
                                kl_e_svtCv_mean = kl_e_svtCv_mean, kl_e_svtCv_se = kl_e_svtCv_se, 
                                kl_e_svtOracle_mean = kl_e_svtOracle_mean, kl_e_svtOracle_se = kl_e_svtOracle_se, 
                                sse_sh_pg_mean = sse_sh_pg_mean, sse_sh_pg_se = sse_sh_pg_se, 
                                sse_sh_zr05_mean = sse_sh_zr05_mean, sse_sh_zr05_se = sse_sh_zr05_se, 
                                sse_sh_zrCv_mean = sse_sh_zrCv_mean, sse_sh_zrCv_se = sse_sh_zrCv_se, 
                                sse_sh_zrOracle_mean = sse_sh_zrOracle_mean, sse_sh_zrOracle_se = sse_sh_zrOracle_se, 
                                sse_sh_svt05_mean = sse_sh_svt05_mean, sse_sh_svt05_se = sse_sh_svt05_se, 
                                sse_sh_svtCv_mean = sse_sh_svtCv_mean, sse_sh_svtCv_se = sse_sh_svtCv_se,
                                sse_sh_svtOracle_mean = sse_sh_svtOracle_mean, sse_sh_svtOracle_se = sse_sh_svtOracle_se, 
                                sse_sp_pg_mean = sse_sp_pg_mean, sse_sp_pg_se = sse_sp_pg_se, 
                                sse_sp_zr05_mean = sse_sp_zr05_mean, sse_sp_zr05_se = sse_sp_zr05_se, 
                                sse_sp_zrCv_mean = sse_sp_zrCv_mean, sse_sp_zrCv_se = sse_sp_zrCv_se, 
                                sse_sp_zrOracle_mean = sse_sp_zrOracle_mean, sse_sp_zrOracle_se = sse_sp_zrOracle_se, 
                                sse_sp_svt05_mean = sse_sp_svt05_mean, sse_sp_svt05_se = sse_sp_svt05_se, 
                                sse_sp_svtCv_mean = sse_sp_svtCv_mean, sse_sp_svtCv_se = sse_sp_svtCv_se, 
                                sse_sp_svtOracle_mean = sse_sp_svtOracle_mean, sse_sp_svtOracle_se = sse_sp_svtOracle_se,
                                sse_bc_pg_mean = sse_bc_pg_mean, sse_bc_pg_se = sse_bc_pg_se, 
                                sse_bc_zr05_mean = sse_bc_zr05_mean, sse_bc_zr05_se = sse_bc_zr05_se, 
                                sse_bc_zrCv_mean = sse_bc_zrCv_mean, sse_bc_zrCv_se = sse_bc_zrCv_se, 
                                sse_bc_zrOracle_mean = sse_bc_zrOracle_mean, sse_bc_zrOracle_se = sse_bc_zrOracle_se, 
                                sse_bc_svt05_mean = sse_bc_svt05_mean, sse_bc_svt05_se = sse_bc_svt05_se, 
                                sse_bc_svtCv_mean = sse_bc_svtCv_mean, sse_bc_svtCv_se = sse_bc_svtCv_se,
                                sse_bc_svtOracle_mean = sse_bc_svtOracle_mean, sse_bc_svtOracle_se = sse_bc_svtOracle_se))
  }else if (r == 0){
    summaryDf <- data.frame(cbind(frob_e_pg_mean = frob_e_pg_mean, frob_e_pg_se = frob_e_pg_se, 
                                  frob_e_zr05_mean = frob_e_zr05_mean, frob_e_zr05_se = frob_e_zr05_se, 
                                  frob_e_zrCv_mean = frob_e_zrCv_mean, frob_e_zrCv_se = frob_e_zrCv_se, 
                                  frob_e_zrOracle_mean = frob_e_zrOracle_mean, frob_e_zrOracle_se = frob_e_zrOracle_se, 
                                  kl_e_pg_mean = kl_e_pg_mean, kl_e_pg_se = kl_e_pg_se, 
                                  kl_e_zr05_mean = kl_e_zr05_mean, kl_e_zr05_se = kl_e_zr05_se, 
                                  kl_e_zrCv_mean = kl_e_zrCv_mean, kl_e_zrCv_se = kl_e_zrCv_se, 
                                  kl_e_zrOracle_mean = kl_e_zrOracle_mean, kl_e_zrOracle_se = kl_e_zrOracle_se, 
                                  sse_sh_pg_mean = sse_sh_pg_mean, sse_sh_pg_se = sse_sh_pg_se, 
                                  sse_sh_zr05_mean = sse_sh_zr05_mean, sse_sh_zr05_se = sse_sh_zr05_se, 
                                  sse_sh_zrCv_mean = sse_sh_zrCv_mean, sse_sh_zrCv_se = sse_sh_zrCv_se, 
                                  sse_sh_zrOracle_mean = sse_sh_zrOracle_mean, sse_sh_zrOracle_se = sse_sh_zrOracle_se, 
                                  sse_sp_pg_mean = sse_sp_pg_mean, sse_sp_pg_se = sse_sp_pg_se, 
                                  sse_sp_zr05_mean = sse_sp_zr05_mean, sse_sp_zr05_se = sse_sp_zr05_se, 
                                  sse_sp_zrCv_mean = sse_sp_zrCv_mean, sse_sp_zrCv_se = sse_sp_zrCv_se, 
                                  sse_sp_zrOracle_mean = sse_sp_zrOracle_mean, sse_sp_zrOracle_se = sse_sp_zrOracle_se,
                                  sse_bc_pg_mean = sse_bc_pg_mean, sse_bc_pg_se = sse_bc_pg_se, 
                                  sse_bc_zr05_mean = sse_bc_zr05_mean, sse_bc_zr05_se = sse_bc_zr05_se, 
                                  sse_bc_zrCv_mean = sse_bc_zrCv_mean, sse_bc_zrCv_se = sse_bc_zrCv_se, 
                                  sse_bc_zrOracle_mean = sse_bc_zrOracle_mean, sse_bc_zrOracle_se = sse_bc_zrOracle_se))
  }
  
  pathSimRes <- paste(pathCompCount, "/summary.csv", sep = "")
  write.csv(summaryDf, file = pathSimRes)
  
  end <- Sys.time()
  print("Main function ends")
  print(end-start)
  return(summaryDf = summaryDf)
}

#----------------------------------------------------------------------------------------
#  Combine the summary of composition estimation
#  Input:
#              nRep ------ number of replications
#                 p ------ number of columns
#       containZero ------ contain zero if containZero = 1; else, X0 does not contain zero 
#           lowRank ------ 1 if low rank matrix; 0 if approxmiately low rank matrix
#          gammaMax ------ the maximum of total count fraction 
#-----------------------------------------------------------------------------------------
CombSummary <- function(nRep = 50, p = 50, containZero = 0, lowRank = 1, gammaMax = 6){
  
  if ((lowRank == 1) && (containZero == 0)){
    pathCompCount = paste("~/comp-est/simulation/nonzero/p", toString(p),  "/lowRankN", sep = "")
  }else if((lowRank == 1) && (containZero == 1)){
    pathCompCount = paste("~/comp-est/simulation/zero/p", toString(p), "/lowRankN", sep = "")
  }else if((lowRank == 0) && (containZero == 0)){
    pathCompCount = paste("~/comp-est/simulation/nonzero/p", toString(p), "/appRankN", sep = "")
  }else if((lowRank == 0) && (containZero == 1)){
    pathCompCount = paste("~/comp-est/simulation/zero/p", toString(p), "/appRankN", sep = "")
  }
  
  fileCompCountRes = paste(pathCompCount, "1/summary.csv", sep = "")
  dfSummary <- read.csv(fileCompCountRes)
  
  for (i in 2:gammaMax){
    fileCompCountRes = paste(pathCompCount, toString(i), "/summary.csv", sep = "")
    dfSummaryUpdate <- read.csv(fileCompCountRes)
    dfSummary <- rbind(dfSummary, dfSummaryUpdate)
  }
  
  if ((lowRank == 1) && (containZero == 0)){
    pathSimRes = paste("~/comp-est/simulation/nonzero/p", toString(p),  "/lowRankSummary.csv", sep = "")
  }else if((lowRank == 1) && (containZero == 1)){
    pathSimRes = paste("~/comp-est/simulation/zero/p", toString(p), "/lowRankSummary.csv", sep = "")
  }else if((lowRank == 0) && (containZero == 0)){
    pathSimRes = paste("~/comp-est/simulation/nonzero/p", toString(p), "/appRankSummary.csv", sep = "")
  }else if((lowRank == 0) && (containZero == 1)){
    pathSimRes = paste("~/comp-est/simulation/zero/p", toString(p), "/appRankSummary.csv", sep = "")
  }
  
  write.csv(dfSummary, file = pathSimRes)
  return(dfSummary = dfSummary)
}

#----------------------------------------------------------------------------------------
#  main function to run the simulation testing
#  Input:
#              nRep ------ number of replications
#                 n ------ number of rows
#                 p ------ number of columns
#       containZero ------ contain zero if containZero = 1; else, X0 does not contain zero 
#           lowRank ------ 1 if low rank matrix; 0 if approxmiately low rank matrix
#                 r ------ the true rank in the low-rank setting or the parameter in approximately low-rank setting
#          gammaMax ------ the maximum of total count fraction   
#-----------------------------------------------------------------------------------------
mainRun <- function(nRep = 10, n = 100, p = 50, containZero = 0, lowRank = 1, r = 20, gammaMax = 6){
  
  for (i in (1: gammaMax)){
    # Generate the replication
    mainGenerate(n, p, lowRank, para = r, containZero = containZero, nRep = nRep, gamma = i)
    
    # Composition estimation
    if (lowRank == 1){
      mainEst(nRep, p, containZero, lowRank, r, gamma = i) 
    }else{
      mainEst(nRep, p, containZero, lowRank, r = 0, gamma = i)  
    }
    
    # Record the result
    if (lowRank == 1){
      mainEstSummary(nRep, p, containZero, lowRank, r, gamma = i)
    }else{
      mainEstSummary(nRep, p, containZero, lowRank, r = 0, gamma = i)
      }
  }
  
  # Summarize the result 
  CombSummary(nRep, p, containZero, lowRank, gammaMax)
}

#----------------------------------------------------------------------------------------
# main function begins

# mainRun(nRep = 5, n = 50, p = 50, containZero = 0, lowRank = 1, r = 20, gammaMax = 1)
  
  