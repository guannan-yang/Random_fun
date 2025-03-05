# Permutation test for DIABLO model with a certain configuration
## This function is to be used in conjunction with the block.splsda function from the mixOmics package
## Empirical p-values are computed for each variable in each omic and for each component in the DIABLO model. The empirical p-values are calculated based on the null distribution of the loadings obtained by permuted Y, and adjusted for multiple testing using the false discovery rate (FDR) method.
## n_per should larger than 1/(target.fdr/max(ncol(X[[i]]))) to ensure the stability of the p-values; if not, the adjusted p-values may be unreliable, the raw p-values can be used instead.
## Parallel computing is used.

diablo.perm.test <- function(X, Y, ncomp, keepX, design, ..., n_per = 1000, seed = 1, target.fdr = 0.05) {

  # Fit the original DIABLO model
  diablo.orig <- block.splsda(X, Y, ncomp = ncomp, keepX = keepX, design = design, ...)
  
  # Run permutation tests in parallel
  perm.diablo <- foreach(i = 1:n_per, .packages = c("mixOmics"), .options.RNG = seed) %dorng% {
    perm.y <- sample(Y)
    diablo <- block.splsda(X, perm.y, ncomp = ncomp, keepX = keepX, design = design, ...)
    return(diablo)
  }
  
  # Extract loading matrices for each omic
  omics.names <- names(X)
  n_omics <- length(omics.names)
  
  perm.loadings <- vector("list", n_omics)
  for (i in 1:n_omics) {
    perm.loadings[[i]] <- array(NA, dim = c(n_per, ncol(X[[i]]), ncomp),
                                dimnames = list(NULL, colnames(X[[i]]), paste0("comp", 1:ncomp)))
    
    for (j in 1:n_per) {
      perm.loadings[[i]][j, , ] <- perm.diablo[[j]]$loadings[[i]][, 1:ncomp]
    }
  }
  
  # Compute empirical p-values
  pvalues <- vector("list", n_omics)
  adj.p <- vector("list", n_omics)
  
  for (i in 1:n_omics) {
    pvalues[[i]] <- array(NA, dim = c(ncol(X[[i]]), ncomp),
                          dimnames = list(colnames(X[[i]]), paste0("comp", 1:ncomp)))
    adj.p[[i]] <- array(NA, dim = c(ncol(X[[i]]), ncomp),
                        dimnames = list(colnames(X[[i]]), paste0("comp", 1:ncomp)))
    for (k in 1:ncomp) {
      for (j in 1:ncol(X[[i]])) {
        pvalues[[i]][j, k] <- sum(abs(perm.loadings[[i]][, j, k]) >= abs(diablo.orig$loadings[[i]][j, k])) / n_per
        adj.p[[i]][j, k] <- 0
      }
      # adjust p-values for multiple testing -- ncol(X[[i]]) tests for each component
      adj.p[[i]][, k] <- p.adjust(pvalues[[i]][, k], method = "fdr")
    }
    
  }
  
  # Select significant variables based on target FDR
  selected.vars <- vector("list", n_omics)
  for (i in 1:n_omics) {
    selected.vars[[i]] <- array(FALSE, dim = c(ncol(X[[i]]), ncomp),
                                dimnames = list(colnames(X[[i]]), paste0("comp", 1:ncomp)))
    
    for (j in 1:ncol(X[[i]])) {
      for (k in 1:ncomp) {
        if (pvalues[[i]][j, k] <= target.fdr) {
          selected.vars[[i]][j, k] <- TRUE
        }
      }
    }
  }
  
  # Extract names of selected variables and their corresponding p-values
  selected.var.names <- vector("list", n_omics)
  for (i in 1:n_omics) {
    selected.var.names[[i]] <- vector("list", ncomp)

    for (j in 1:ncomp) {

      selected.var.names[[i]][[j]] <- data.frame(name = colnames(X[[i]])[selected.vars[[i]][, j]], raw.p = pvalues[[i]][, j][selected.vars[[i]][, j]], adj.p = adj.p[[i]][, j][selected.vars[[i]][, j]])
    }
    names(selected.var.names[[i]]) <- paste0("comp", 1:ncomp)
  }
  
  # Extract the selected in the original DIABLO model and their corresponding p-values
  original.selected.var.names <- vector("list", n_omics)
  for (i in 1:n_omics) {
    original.selected.var.names[[i]] <- vector("list", ncomp)
    
    for (j in 1:ncomp) {
      original.selected.var.names[[i]][[j]] <- data.frame(name = selectVar(diablo.orig, comp = j)[[i]]$name, 
                                                          loading = selectVar(diablo.orig, comp = j)[[i]]$value,
                                                          raw.p = pvalues[[i]][, j][selectVar(diablo.orig, comp = j)[[i]]$name],
                                                          adj.p = adj.p[[i]][, j][selectVar(diablo.orig, comp = j)[[i]]$name])
    }
    names(original.selected.var.names[[i]]) <- paste0("comp", 1:ncomp)
  }
  
  # name the results
  names(pvalues) <- names(selected.var.names) <- names(original.selected.var.names) <- omics.names
  
  return(list(pvalues = pvalues, selected.var = selected.var.names, original.selected.var = original.selected.var.names))
}
