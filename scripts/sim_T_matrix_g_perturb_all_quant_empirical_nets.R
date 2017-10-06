#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Simulations of pertubations (disruption of gene flow) on the T-matrix for 29 empirical
#   quantitative mutualistic networks.
#
# Returns:
#   Saves a csv spreadsheet with the results.

# trait matching function
source("functions/MatchingMutNet.R")

# network files and names
net_files = dir("data/empirical_networks/weighted/")
net_names = substr(net_files, start = 1, stop = nchar(net_files) - 4) 

# defining folder to store results
folder = "~/LUCAS/spatial_coevo_mutnet_results/data/simulations_geneflow_removal"

# define m_A and m_B sequences
mA_seq = c(0.1, 0.3, 0.3, 0.5, 0.5, 0.5, 0.7, 0.7, 0.7, 0.7,
           0.9, 0.9, 0.9, 0.9, 0.9)
mB_seq = c(0.1, 0.1, 0.3, 0.1, 0.3, 0.5, 0.1, 0.3, 0.5, 0.7,
           0.1, 0.3, 0.5, 0.7, 0.9)
# define gene flow value
g = 0.3
# for the file name
g_char = gsub(".", "", as.character(g), fixed = TRUE)
# define perturbation sequence
pert_g_seq = seq(0, 1, by = 0.05)
# defining theta distribution
theta_A_min = 0
theta_A_max = 10
theta_B_min = 10
theta_B_max = 20
# number of theta vectors 
n_theta = 30
# number of replicates of perturbation simulations
n_rep_pert = 30

for(i in 1:length(net_files)) {
  
  print(net_names[i])
  
  # read empirical matrix
  mat = as.matrix(read.table(paste("data/empirical_networks/weighted/", net_files[i], sep = ""), 
                             sep = " ", header = FALSE))
  n_row = nrow(mat)
  n_col = ncol(mat)
  n_sp = n_row + n_col
  
  results_df = data.frame(network = rep(net_names[i], times = n_rep_pert*n_theta*2*length(mA_seq)*length(pert_g_seq)),
                          rep_pert = rep(rep(1:n_rep_pert, each = length(pert_g_seq)*n_theta*2), length(mA_seq)),
                          rep_theta = rep(rep((1:n_theta), each = 2), times = n_rep_pert*length(mA_seq)*length(pert_g_seq)),
                          site = rep(rep(c("A", "B"), times = n_theta), times = n_rep_pert*length(mA_seq)*length(pert_g_seq)),
                          pert = rep(rep(pert_g_seq, each = n_theta*2), times = n_rep_pert*length(mA_seq)),
                          m_A = rep(rep(mA_seq, each = length(pert_g_seq)), each = n_rep_pert*n_theta*2),
                          m_B = rep(rep(mB_seq, each = length(pert_g_seq)), each = n_rep_pert*n_theta*2),
                          matching = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2),
                          mean_eigen = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2),
                          sd_eigen = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2),
                          eigen1 = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2),
                          eigen2 = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2),
                          det_T = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2),
                          mean_col_sd = rep(NA, times = length(mA_seq)*length(pert_g_seq)*n_rep_pert*n_theta*2))

  # Q and F matrices
  Q = rbind(cbind(matrix(0, n_row, n_row), mat),
            cbind(t(mat), matrix(0, n_col, n_col)))
  F_matrix = Q
  F_matrix[F_matrix != 0] = 1
  Q_norm = Q/apply(Q, 1, sum)
  Q_2sites = rbind(cbind(Q_norm, matrix(0, n_sp, n_sp)),
                   cbind(matrix(0, n_sp, n_sp), Q_norm))

  # identity matrix
  I = diag(rep(1, 2*n_sp))
  
  for (j in 1:length(mA_seq)) {
    
    # matrix M
    m_A = rep(mA_seq[j], n_sp)
    m_B = rep(mB_seq[j], n_sp)
    M = diag(c(m_A, m_B))
    
    print(mA_seq[j])
    print(mB_seq[j])
    
    for (k in 1:n_rep_pert) {
      
      # matrix G
      g_mult = rep(g, n_sp)
      G = rbind(cbind(diag(1 - g_mult), diag(g_mult)),
                cbind(diag(g_mult), diag(1 - g_mult)))
      
      # sampling species perturbation order
      sp_order = sample(1:n_sp, replace = FALSE)
      
      for (l in 1:length(pert_g_seq)) {
        
        # perturbing matrix G
        sp_index = round(pert_g_seq[l]*n_sp)
        if (sp_index != 0) {
          pert_sp = sp_order[1:sp_index]
          G[c(pert_sp, pert_sp + n_sp), c(pert_sp, pert_sp + n_sp)] = diag(1, length(pert_sp)*2) 
        }
        
        # build T matrix
        T_matrix = solve(solve(G) -  M %*% Q_2sites) %*% (I - M)
        # T matrix eigenvalues
        eigens = Re(eigen(T_matrix)$values)
        # mean of eigenvalues
        mean_eigens = mean(eigens)
        # sd of eigenvalues
        sd_eigens = sd(eigens)
        # 1st and 2nd eigenvalues
        eigen1 = eigens[1]
        eigen2 = eigens[2]
        # determinant
        det_T = det(T_matrix)
        # mean column sd
        mean_col_sd = mean(apply(T_matrix, 2, sd))
        
        for (t in 1:n_theta) {
          
          # theta vector to calculate trait matching
          theta = c(runif(n_sp, min = theta_A_min, max = theta_A_max),
                    runif(n_sp, min = theta_B_min, max = theta_B_max))
          
          # equilibrium trait values
          z = T_matrix %*% theta
          
          # calculate trait matching
          match_A = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix, 
                                   z = z[1:n_sp], method = "exponential", alpha = 0.2)
          match_B = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix, 
                                   z = z[(n_sp+1):(2*n_sp)], method = "exponential", alpha = 0.2)
          
          # mean matching
          mean_mat_A = match_A[[1]]
          mean_mat_B = match_B[[1]]
          
          # adding results to data frame
          results_df$matching[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                                results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                                results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(mean_mat_A, mean_mat_B)
          
          results_df$mean_eigen[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                                  results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                                  results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(mean_eigens, mean_eigens)
          
          results_df$sd_eigen[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                                results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                                results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(sd_eigens, sd_eigens)
          
          results_df$eigen1[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                              results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                              results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(eigen1, eigen1)
          
          results_df$eigen2[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                              results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                              results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(eigen2, eigen2)
          
          results_df$det_T[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                             results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                             results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(det_T, det_T)
          
          results_df$mean_col_sd[results_df$network ==  net_names[i] & results_df$m_A ==  mA_seq[j] & 
                                   results_df$m_B ==  mB_seq[j] & results_df$rep_pert == k & 
                                   results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = c(mean_col_sd, mean_col_sd)
           
        }
      }
    }
  }
  
  # saving results
  write.csv(results_df, file = paste(folder, "/", net_names[i], "_geneflow_removal_from_g", 
                                     g_char, ".csv", sep = ""), row.names = FALSE)
  
}

#-----------------------------------------------------------------------------------------------------#