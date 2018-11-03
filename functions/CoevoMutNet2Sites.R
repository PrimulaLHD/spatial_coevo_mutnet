#-----------------------------------------------------------------------------------------------------#

CoevoMutNet2Sites = function(n_sp, f, g, alpha, phi_A, phi_B, theta_A, theta_B, init_A, init_B, m_A, m_B,
                             epsilon, t_max) {
  # Simulates the coevolutionary dynamics on a mutualistic network of species interactions at 2 sites
  # connected by gene flow
  #
  # Args:
  #   n_sp: total number of species in the mutualistic network
  #   f: square adjacency matrix representing the mutualistic network
  #   g: vector of gene flow values (symmetric gene flow between site A and B) for all species
  #   alpha: parameter alpha, sensitivity of selection to trait matching
  #   phi_A: vector of phi values for site A
  #   phi_B: vector of phi values for site B
  #   theta_A: vector of environmental optimum values for site A
  #   theta_B: vector of environmental optimum values for site B
  #   init_A: vector of initial trait values for site A
  #   init_B: vector of initial trait values for site B
  #   m_A: vector of proportion of selection due to mutualism at site A
  #   m_B: vector of proportion of selection due to mutualism at site B
  #   epsilon: value to determine when equilibrium is reached
  #   t_max: maximum number of timesteps allowed
  #   
  # Obs: 
  #   All vectors need to have first the row species attributes and then the column species
  #   attributes (e.g. c(row_sp[1], ..., row_sp[nrow], col_sp[1], ..., col_sp[ncol]))
  #
  # Returns:
  #   A list containing two matrices. The first matrix contains, in each row t, the trait values (z)
  #   of all species in site A at time t. The second matrix contains, in each row, the trait values 
  #   of all species in site B at time t.
  
  z_mat_A = matrix(NA, nrow = t_max, ncol = n_sp) # matrix to store z values
  z_mat_B = matrix(NA, nrow = t_max, ncol = n_sp) 
  z_mat_A[1, ] = init_A # initial trait values 
  z_mat_B[1, ] = init_B 
  
  for (t in 1:(t_max - 1)) { # simulation runs for a maximum of t_max timesteps
    z_A = z_mat_A[t, ] # current z values
    z_B = z_mat_B[t, ]
    z_dif_A = t(f*z_A) - f*z_A # matrix with all trait differences
    z_dif_B = t(f*z_B) - f*z_B
    q_A = f*(exp(-alpha * (z_dif_A^2))) # calculating matrix q
    q_B = f*(exp(-alpha * (z_dif_B^2)))
    q_n_A = q_A / apply(q_A, 1, sum) # normalizing the matrix
    q_n_B = q_B / apply(q_B, 1, sum)
    q_m_A = q_n_A * m_A # multiplying each row i of matrix q by m[i]
    q_m_B = q_n_B * m_B 
    sel_dif_A = q_m_A * z_dif_A # calculating selection differentials
    sel_dif_B = q_m_B * z_dif_B
    r_mut_A = phi_A * apply(sel_dif_A, 1, sum) # response to selection related to mutualism
    r_env_A = phi_A * (1 - m_A) * (theta_A - z_A) # response to selection related to the environment
    r_mut_B = phi_B * apply(sel_dif_B, 1, sum)
    r_env_B = phi_B * (1 - m_B) * (theta_B - z_B) 
    z_mat_A[t+1, ] = (1 - g) * (z_A + r_mut_A + r_env_A) + g * (z_B + r_mut_B + r_env_B) # updating z values
    z_mat_B[t+1, ] = (1 - g) * (z_B + r_mut_B + r_env_B) + g * (z_A + r_mut_A + r_env_A) 
    dif_A = mean(abs(z_A - z_mat_A[t+1, ])) # computing the mean difference between old and new z values
    dif_B = mean(abs(z_B - z_mat_B[t+1, ]))
    if ((dif_A < epsilon) & (dif_B < epsilon))
      break
  }
  
  z_list = list(z_mat_A[1:(t+1), ], z_mat_B[1:(t+1), ])
  return(z_list)
}

#-----------------------------------------------------------------------------------------------------#