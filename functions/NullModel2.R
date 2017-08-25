#-----------------------------------------------------------------------------------------------------#

NullModel2 = function(mat, mat_name, write = TRUE, R = 100) {
  # Generates theoretical networks according to null model 2 (Bascompte et al 2003). Modified from 
  # M. M. Pires script.
  #
  # Args:
  #   mat: a bipartite adjacency matrix
  #   mat_name: character string containing the matrix name
  #   write: should the matrices be saved in the current working directory?
  #   R: number of matrices to be generated
  #
  # Returns:
  #   R theoretical matrices generated according to null model 2.    

  m = dim(mat)[1] # number of rows
  n = dim(mat)[2] # number of columns
  rmarg = apply(mat, 1, sum) # degrees of row species
  cmarg = apply(mat, 2, sum) # degrees of column species
  
  # Generating a matrix with probablities for each cell
  prob = matrix(0, m, n)
  for (i in 1:m)
    for (j in 1:n)
      prob[i, j] = ((rmarg[i]/n) + (cmarg[j]/m))/2
  
  # Filling theoretical matrices
  mat_t = list() # To store theoretical matrices
  s = 1
  while (s <= R) {
    rand = matrix(runif(m*n), m, n) # Matrix with numbers from a uniform distribution between 0 and 1
    aux = matrix(0, nrow = m, ncol = n) # Just to avoid indexing problems
    aux[which(rand < prob)] = 1 # Filling empty matrix
    
    # To avoid zeroed rows or columns 
  	rm_aux = apply(aux, 1, sum) # Row degrees
  	if (any(rm_aux == 0))	{
  	  zeroed_rows = which(rm_aux == 0) # Indexes of zeroed rows
      cols = sample(1:n, sum(rm_aux == 0), replace = TRUE) # Randomly selecting columns
  	  for (i in 1:length(zeroed_rows))
  	    aux[zeroed_rows[i], cols[i]] = 1 # Filling in the sampled positions
  	}
  	cm_aux = apply(aux, 2, sum) # Column degrees
  	if (any(cm_aux == 0))	{
  	  zeroed_cols = which(cm_aux == 0) # Indexes of zeroed columns
      rows = sample(1:m, sum(cm_aux == 0), replace = TRUE) # Randomly selecting rows
  	  for (i in 1:length(zeroed_cols))
  	    aux[rows[i], zeroed_cols[i]] = 1 # Filling in the sampled positions
  	}
  	
  	# Storing matrices
    if (sum(aux) == sum(mat)) { # Runs againf if the resulting matrix has a different number of interactions
      mat_t[[s]] = aux	# Store the matrix 
      s = s + 1
    }
  }
  
  	# Saving matrices
  	if (write == TRUE) {
      curr_dir = getwd()
      new_dir = paste(mat_name, "Null", sep = "_")
      full_dir = paste(curr_dir, new_dir, sep = "/")
      dir.create(full_dir)
      for(i in 1:R) {
        name = paste(new_dir, "_", i, ".txt", sep = "")
        write.table(mat_t[[i]], file = paste(full_dir, name, sep="/"), row.names = FALSE, col.names = FALSE,
                    append = TRUE, quote = FALSE)  
      }
  	}
    
    return(mat_t)
}

#-----------------------------------------------------------------------------------------------------#