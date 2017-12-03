## scatter plot and box plot for real data analysis
require(ggplot2)
require(reshape)

#----------------------------------------------------------------------------------------
#  Diversity index
#  Input:
#                 X ------ n x p composition matrix (row/column is sample/variable)
#  Output:
#                sh ------ shannon index (n x 1 vector)
#                sp ------ simpson index (n x 1 vector)
#                bc ------ bray-curtis index (n(n-1)/2 x 1 vector)
#-----------------------------------------------------------------------------------------
diversity <- function(X){
  # shannon index
  sh_ele <- -X * log(X)
  sh_ele[is.nan(sh_ele)] <- 0
  sh <- rowSums(sh_ele)
  
  # simpson index
  sp <- rowSums(X * X)
  
  # bray-curtis index
  n <- nrow(X)
  bc <- matrix(vegdist(X), n*(n-1)/2, 1)
  return(list(sh = sh, sp = sp, bc = bc))
}

# boxplot for composition estimation
comp_boxplot <- function(W, XHat, Output){
  n <- nrow(XHat)
  p <- ncol(XHat)
  
  ## Load the species name
  speciesNames <- colnames(W)
  taxa <- matrix(unlist(strsplit(speciesNames, "\\.")), 6)
  phyla <- taxa[2, ]
  genera <- taxa[6, ]
  
  ## Split the compsition into zero and non-zero
  xZeroIndex <- (W == 0)
  
  ## Prepare the dataframe
  dat <- data.frame(Genera = c(rep(genera,each = n)),zero_index = factor(c(xZeroIndex)), Composition = c(log(XHat)))
  dat$zero_index <- factor(dat$zero_index, labels = c(expression("italic(M)^c","italic(M)")))
  
  ## boxplot
  setEPS(width = 18, height = 10)
  postscript(file = Output)
  print(ggplot(dat, aes(x = Genera, y = Composition), mar = c(0.5,0.5,0.5,0.5)) + 
    geom_boxplot(aes(fill = zero_index), coef = 10) +
    ylab(expression(log(widehat(X)))) + xlab("Genera") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15, angle = 60, hjust = 1), legend.title = element_blank(), legend.text = element_text(size = 30)) +
    scale_fill_manual(values = c("#D55E00", "#56B4E9"),
                      labels = c(expression(Omega),expression(Omega^c))) +
    theme(axis.title = element_text(size = 30), legend.text.align = 0, strip.text = element_text(size = 30), panel.grid.major = element_blank(),
          panel.background = element_blank(), panel.grid.minor = element_blank()))
  dev.off()
}

# scatter plot for sequecning depth of the estimated composition
comp_seq_depth_scatterplot <- function(W, XHat, Output){
  
  # get sequencing depth
  sd <- log(rowSums(W))
  sd_v <- as.vector(sd)
  
  # get the zero index
  xZeroIndex <- (W == 0)
  xZeroIndex_v <- as.vector(xZeroIndex) 
  
  # combine three matrices together to get the composition estimation for each of sequecning depth
  XHat_v <- as.vector(XHat) 
  df_all <- cbind(sd_v, xZeroIndex_v)
  df_all <- cbind(df_all, XHat_v)
  df_zero_ind_comp <- df_all[df_all[,2] == 1,]
  
  plot(df_zero_ind_comp[,1], df_zero_ind_comp[,3])
}

# scatter plot for Shannon index comparison of the estimated composition
comp_shannon_scatterplot <- function(XHat1, XHat2, NameEst1, NameEst2){
  
  # calculate the Shannon's index
  d1 <- diversity(XHat1)
  d2 <- diversity(XHat2)
  
  sh1 <- d1$sh
  sh2 <- d2$sh
  
  # combine two Shannon index together
  x_name <- NameEst1
  y_name <- NameEst2
  df_shannon <- melt(data.frame(sh1, sh2))
  colnames(df_shannon) <- c(x_name, y_name)
  
  plot(sh1, sh2)
}

# scatter plot for simpson index comparison of the estimated composition
comp_simpson_scatterplot <- function(XHat1, XHat2, NameEst1, NameEst2){
  
  # calculate the Simpson's index
  d1 <- diversity(XHat1)
  d2 <- diversity(XHat2)
  
  sp1 <- d1$sp
  sp2 <- d2$sp
  
  # combine two Simpson index together
  x_name <- NameEst1
  y_name <- NameEst2
  df_simpson <- melt(data.frame(sp1, sp2))
  colnames(df_simpson) <- c(x_name, y_name)
  
  plot(sp1, sp2)
}