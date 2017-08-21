#-----------------------------------------------------------------------------------------------------#

CoevoMutNet2Sites_array <- function(n_sp, f, g, h, alpha, theta_A, theta_B, init_A, init_B, m_A, m_B,
                                    epsilon, t_max) {
  # Simulates the coevolutionary dynamics on a mutualistic network of species interactions at 2 sites
  # connected by gene flow
  #
  # Args:
  #   n_sp: total number of species in the mutualistic network
  #   f: square adjacency matrix representing the mutualistic network
  #   g: vector of gene flow values (symmetric gene flow between site A and B) for all species
  #   h: vector of heritability values
  #   alpha: parameter alpha, sensitivity of selection to trait matching
  #   init_A: vector of initial trait values for site A
  #   init_B: vector of initial trait values for site B
  #   theta_A: vector of environmental optimum values for site A
  #   theta_B: vector of environmental optimum values for site B
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
  #   An array containing two sheets. The first sheet is a matrix that contains, in each row t, the trait values (z)
  #   of all species in site A at time t. The second sheet is a matrix that contains, in each row, the trait values 
  #   of all species in site B at time t.

    Z <- array(0, c(t_max, n_sp, 2))
    dimnames(Z) [[3]] <- c('A', 'B')
    
    ## initial trait values 
    Z [1, , 'A'] <- init_A
    Z [1, , 'B'] <- init_B

    ## simulation runs for a maximum of t_max timesteps
    for (t in 1:(t_max - 1))
    {
        ## current z values
        z_A <- Z[t, , 'A'] 
        z_B <- Z[t, , 'B']

        ## matrix with all trait differences
        z_dif_A <- t(f*z_A) - f*z_A 
        z_dif_B <- t(f*z_B) - f*z_B

        ## calculating matrix q
        q_A <- f*(exp(-alpha * (z_dif_A^2))) 
        q_B <- f*(exp(-alpha * (z_dif_B^2)))

        ## normalizing the matrix
        q_n_A <- q_A / apply(q_A, 1, sum) 
        q_n_B <- q_B / apply(q_B, 1, sum)

        ## multiplying each row i of matrix q by m[i]
        q_m_A <- q_n_A * m_A 
        q_m_B <- q_n_B * m_B

        ## calculating selection differentials
        sel_dif_A <- q_m_A * z_dif_A 
        sel_dif_B <- q_m_B * z_dif_B

        ## response to selection related to mutualism
        r_mut_A <- h * apply(sel_dif_A, 1, sum) 
        r_mut_B <- h * apply(sel_dif_B, 1, sum)

        ## response to selection related to the environment
        r_env_A <- h * (1 - m_A) * (theta_A - z_A) 
        r_env_B <- h * (1 - m_B) * (theta_B - z_B)

        ## updating z values
        Z[t+1, , 'A'] <- (1 - g) * (z_A + r_mut_A + r_env_A) + g * (z_B + r_mut_B + r_env_B) 
        Z[t+1, , 'B'] <- (1 - g) * (z_B + r_mut_B + r_env_B) + g * (z_A + r_mut_A + r_env_A) 

        ## computing the mean difference between old and new z values
        dif_A <- mean(abs(Z[t+1, , 'A'] - Z[t, , 'A'])) 
        dif_B <- mean(abs(Z[t+1, , 'B'] - Z[t, , 'B']))

        ## condition for stopping iterations
        if ((dif_A < epsilon) & (dif_B < epsilon))
            break
  }
    return(Z)
}

#-----------------------------------------------------------------------------------------------------#
