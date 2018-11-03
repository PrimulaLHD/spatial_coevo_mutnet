#-----------------------------------------------------------------------------------------------------#

CoevoMutNet2SitesTurnover = function(n_sp_A, n_sp_B, sp_labels_A, sp_labels_B, f_A, f_B, g_A, g_B, alpha, phi_A, phi_B, theta_A, theta_B,
                                     init_A, init_B, m_A, m_B, epsilon, t_max) {
  # Simulates the coevolutionary dynamics on two different mutualistic network of species interactions 
  # connected by gene flow
  #
  # Args:
  #   n_sp_A: total number of species at site A
  #   n_sp_B: total number of species at site B
  #   sp_labels_A: labels of species at site A
  #   sp_labels_B: labels of species at site B
  #   f_A: square adjacency matrix representing the mutualistic network at site A (n_sp_A x n_sp_A)
  #   f_B: square adjacency matrix representing the mutualistic network at site B (n_sp_B x n_sp_B)
  #   g_A: vector with gene flow values from site A to site B (length: n_sp_A)
  #   g_B: vector with gene flow values from site B to site A (length: n_sp_B)
  #   alpha: parameter alpha, sensitivity of selection to trait matching
  #   phi_A: vector of phi values for site A (length: n_sp_A)
  #   phi_B: vector of phi values for site B (length: n_sp_B)
  #   theta_A: vector of environmental optimum values for site A (length: n_sp_A)
  #   theta_B: vector of environmental optimum values for site B (length: n_sp_B)
  #   init_A: vector of initial trait values for site A (length: n_sp_A)
  #   init_B: vector of initial trait values for site B (length: n_sp_B)
  #   m_A: vector of proportion of selection due to mutualism at site A (length: n_sp_A)
  #   m_B: vector of proportion of selection due to mutualism at site B (length: n_sp_B)
  #   epsilon: value to determine when equilibrium is reached
  #   t_max: maximum number of timesteps allowed
  #   
  # Obs: 
  #   If the interaction data comes from a bipartite matrix, all vectors need to 
  #   have first the row species attributes and then the column species attributes 
  #   (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
  #
  # Returns:
  #   A list containing two matrices. The first matrix contains, in each row t, the trait values (z)
  #   of all species in site A at time t. The second matrix contains, in each row, the trait values 
  #   of all species in site B at time t.
  
  # matrix to store z values
  z_mat_A = matrix(NA, nrow = t_max, ncol = n_sp_A)
  z_mat_B = matrix(NA, nrow = t_max, ncol = n_sp_B) 
  # initial trait values 
  z_mat_A[1, ] = init_A 
  z_mat_B[1, ] = init_B 
  
  # simulation runs for a maximum of t_max timesteps
  for (t in 1:(t_max - 1)) {
    # current z values
    z_A = z_mat_A[t, ]
    z_B = z_mat_B[t, ]
    # matrix with all trait differences
    z_dif_A = t(f_A*z_A) - f_A*z_A
    z_dif_B = t(f_B*z_B) - f_B*z_B
    # calculating matrix q
    q_A = f_A*(exp(-alpha * (z_dif_A^2))) 
    q_B = f_B*(exp(-alpha * (z_dif_B^2)))
    # normalizing the matrix
    q_n_A = q_A / apply(q_A, 1, sum) 
    q_n_B = q_B / apply(q_B, 1, sum)
    # multiplying each row i of matrix q by m[i]
    q_m_A = q_n_A * m_A 
    q_m_B = q_n_B * m_B 
    # calculating selection differentials
    sel_dif_A = q_m_A * z_dif_A 
    sel_dif_B = q_m_B * z_dif_B
    # response to selection related to mutualism
    r_mut_A = phi_A * apply(sel_dif_A, 1, sum) 
    r_mut_B = phi_B * apply(sel_dif_B, 1, sum)
    # response to selection related to the environment
    r_env_A = phi_A * (1 - m_A) * (theta_A - z_A) 
    r_env_B = phi_B * (1 - m_B) * (theta_B - z_B) 
    # trait values for next generation
    z_A_next = (z_A + r_mut_A + r_env_A) 
    z_B_next = (z_B + r_mut_B + r_env_B)
    # labels of species present at both sites
    sp_labels_A_B = intersect(sp_labels_A, sp_labels_B) 
    if (length(sp_labels_A_B) > 0) {
      # weighting trait changes by gene flow of common species
      z_A_next[sp_labels_A_B] = (1 - g_A)[sp_labels_A_B] * z_A_next[sp_labels_A_B] + 
        g_B[sp_labels_A_B] * z_B_next[sp_labels_A_B]
      z_B_next[sp_labels_A_B] = (1 - g_B)[sp_labels_A_B] * z_B_next[sp_labels_A_B] +
        g_A[sp_labels_A_B] * z_A_next[sp_labels_A_B]
    }
    # adding new z values to the matrix
    z_mat_A[t+1, ] = z_A_next 
    z_mat_B[t+1, ] = z_B_next
    # computing the mean difference between old and new z values
    dif_A = mean(abs(z_A - z_mat_A[t+1, ]))
    dif_B = mean(abs(z_B - z_mat_B[t+1, ]))
    if ((dif_A < epsilon) & (dif_B < epsilon))
      break
  }
  
  z_list = list(z_mat_A[1:(t+1), ], z_mat_B[1:(t+1), ])
  return(z_list)
}

#-----------------------------------------------------------------------------------------------------#