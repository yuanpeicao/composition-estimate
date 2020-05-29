# Multi-sample Estimation of Bacterial Composition Matrix in Metagenomics Data

Metagenomics sequencing is routinely applied to quantify bacterial abundances in microbiome studies, where the bacterial composition is estimated based on the sequencing
read counts. Due to limited sequencing depth and DNA dropouts, many rare bacterial taxa might not be captured in the final sequencing reads, which results in many
zero counts. Naive composition estimation using count normalization leads to many zero
proportions, which tend to result in inaccurate estimates of bacterial abundance and diversity. In this repo, we propose a multi-sample approach to estimation of bacterial abundances
in order to borrow information across samples and across species. Empirical results from
real data sets suggest that the composition matrix over multiple samples is approximately
low rank, which motivates a regularized maximum likelihood estimation with a nuclear
norm penalty. An efficient optimization algorithm using the generalized accelerated proximal gradient and Euclidean projection onto simplex space is developed. 

Reference: https://academic.oup.com/biomet/article-abstract/107/1/75/5663561
