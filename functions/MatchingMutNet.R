#-----------------------------------------------------------------------------------------------------#

MatchingMutNet = function(n_sp, n_row, n_col, f, z, method, alpha) {
  # Calculates species trait matching in a mutualistic network.
  #
  # Args:
  #   n_sp: species richness
  #   n_row: number of species in the rows in the bipartite adjacency matrix
  #   n_col: number of species in the columns in the bipartite adjacency matrix
  #   f: square adjacency matrix (n_sp by n_sp)
  #   z: vector of species trait values, first row species then column species
  #   method: a string indicating the metric to be used: "guimaraes", "nuismer", or "exponential"
  #   alpha: constant used for the exponential metric  
  #           
  # Returns:
  #   A list containing the mean value of trait matching at the first position and the pairwise
  #   values at the second position for "guimaraes" and "exponential".
  #   Metric details:
  #   guimaraes: -log of mean absolute difference of traits of interacting species 
  #              (Guimarães et al., 2011, Ecol. Lett.) 
  #   nuismer: correlation between traits of interacting species 
  #            (similar to Nuismer et al., 2013, Evolution)
  #   exponential: exponential expression based on trait differences (exp(-alpha*(z[i]-z[j])^2))
  
  f_traits = f * z # matrix f with trait values
  t_vec_traits = c(t(f_traits)) # vectorizing the transpose matrix
  row_sp_traits = t_vec_traits[1:(n_sp*n_row)] # row sp trait values
  row_sp_traits = row_sp_traits[row_sp_traits > 0]
  vec_traits = c(f_traits)
  col_sp_traits = vec_traits[1:(n_sp*n_row)] # col sp trait values
  col_sp_traits = col_sp_traits[col_sp_traits > 0]
  trait_diff = row_sp_traits - col_sp_traits # pairwise trait differences 
  
  row_names = paste("R", 1:n_row, sep = "") # codes for rows
  col_names = paste("C", 1:n_col, sep = "") # codes for columns
  names_rep = rep(c(row_names, col_names), n_sp) # getting codes for links
  names_node_1 = names_rep[which(f_traits != 0)]
  k = apply(f, 1, sum) # node degrees
  names_node_2 = rep(c(row_names, col_names), k)
  n_links = sum(f)/2 # total number of links
  results_df = data.frame(column = names_node_1[1:n_links], # putting results in data.frame
                          row = names_node_2[1:n_links],
                          matching = rep(NA, n_links))
  # Guimarães et al. (2011) metric
  if (method == "guimaraes") {
    results_df$matching = - log(abs(trait_diff))
    mean_match_guima = mean(results_df$matching)
    return(list(mean_match_guima, results_df))
  }
  # Nuismer et al. (2013) metric
  if (method == "nuismer") {
    match_nuismer = cor(row_sp_traits, col_sp_traits)
    return(list(match_nuismer))
  }
  # Exponential metric
  if (method == "exponential") {
    results_df$matching = exp(-alpha*(trait_diff)^2)
    mean_match_exp = mean(results_df$matching)
    return(list(mean_match_exp, results_df))
  }
}

#-----------------------------------------------------------------------------------------------------#