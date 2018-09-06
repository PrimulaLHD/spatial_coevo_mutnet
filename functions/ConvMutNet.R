#-----------------------------------------------------------------------------------------------------#

ConvMutNet = function(n_sp, n_row, n_col, z, method, alpha) {
  # Calculates trait convergence in a mutualistic network.
  #
  # Args:
  #   n_sp: species richness
  #   n_row: number of species in the rows in the bipartite adjacency matrix
  #   n_col: number of species in the columns in the bipartite adjacency matrix
  #   z: vector of species trait values, first row species then column species
  #   method: a string indicating the metric to be used: "guimaraes", "nuismer", or "exponential"
  #   alpha: constant used for the exponential metric  
  #
  # Returns:
  #   A vector containing two values of trait convergence using the chosen metric. The value for the
  #   row species is on the first position and the value for the column species is on the second position.
  #   Metric details:
  #   guimaraes: -log of mean absolute difference of traits of species in the same group
  #              (Guimarães et al., 2011, Ecol. Lett.) 
  #   nuismer: variance of trait values of species in the same group
  #            (similar to Nuismer et al., 2013, Evolution)
  #   exponential: exponential expression based on trait differences (exp(-alpha*(z[i]-z[j])^2)) 
  
  z_row = z[1:n_row] # z values of row species 
  z_col = z[(n_row+1):n_sp] # z values of column species
  if(length(z_row) > 1) {
    pair_diff_row = abs(combn(z_row, 2, FUN = diff)) # all possible pairwise differences between row species
  } else {
    pair_diff_row = NA
  }
  if (length(z_col) > 1) {
    pair_diff_col = abs(combn(z_col, 2, FUN = diff)) # all possible pairwise differences between column species
  } else {
    pair_diff_col = NA
  }
  # Guimarães et al. (2011) metric
  if (method == "guimaraes") {
    conv_guima = c(-log(mean(pair_diff_row)), -log(mean(pair_diff_col)))
    return(conv_guima)
  }
  # Nuismer et al. (2013) metric
  if (method == "nuismer") {
    conv_nuismer = c(var(z_row), var(z_col))
    return(conv_nuismer)
  }
  # Exponential metric
  if (method == "exponential") {
    conv_exp = c(mean(exp(-alpha*(pair_diff_row)^2)), mean(exp(-alpha*(pair_diff_col)^2)))
    return(conv_exp)
  }
}

#-----------------------------------------------------------------------------------------------------#