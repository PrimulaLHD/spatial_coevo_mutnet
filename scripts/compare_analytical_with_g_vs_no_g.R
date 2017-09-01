#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Computes difference in several metrics using analytical equilibrium values with and without gene flow.
#
# Returns:
#   Saves a csv spreadsheet with the results.

source("functions/MatchingMutNet.R")

# network files and names
net_files = dir("data/empirical_networks/weighted/")
net_names = substr(net_files, start = 1, stop = nchar(net_files) - 4) 
 
# defining m values   
m_seq = seq(0.1, 0.9, by = 0.2)
m_df = data.frame(m_A = m_seq, m_B = m_seq)
m_comb = as.data.frame(t(combn(m_seq, 2)))
names(m_comb) = c("m_A", "m_B")
m_df = rbind(m_df, m_comb)
# defining g values   
g_seq = seq(0, 0.4, by = 0.02)
# defining theta distribution
theta_A_min = 0
theta_A_max = 10
theta_B_min = 10
theta_B_max = 20
# defining number of analyses per combination of m and g
n = 10

# creating data.frame to store results
summary_df = data.frame(network = rep(net_names, each = n*nrow(m_df)*length(g_seq)),
                        mutualism = rep(substr(net_names, start = 1, stop = 2), 
                                        each = n*nrow(m_df)*length(g_seq)),
                        g = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        m_A = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        m_B = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        dif_eigen_1 = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        dif_eigen_2 = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        cos_angle_A = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        cos_angle_B = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        dif_mean_match_A = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)),
                        dif_mean_match_B = rep(NA, length(net_names)*n*nrow(m_df)*length(g_seq)))

# to count the current line of the data frame
counter = 1

for (i in 1:length(net_files)) {
  
  print(net_names[i])
  
  # matrix F and Q
  mat = as.matrix(read.table(paste("data/empirical_networks/weighted/", net_files[i], sep = ""), 
                             sep = " ", header = FALSE))
  n_row = nrow(mat)
  n_col = ncol(mat)
  n_sp = n_row + n_col
  Q = rbind(cbind(matrix(0, n_row, n_row), mat),
            cbind(t(mat), matrix(0, n_col, n_col)))
  F_matrix = Q
  F_matrix[F_matrix != 0] = 1
  Q_norm = Q/apply(Q, 1, sum)
  Q_2sites = rbind(cbind(Q_norm, matrix(0, n_sp, n_sp)),
                   cbind(matrix(0, n_sp, n_sp), Q_norm))
  
  for (j in 1:nrow(m_df)) {
    for (k in 1:length(g_seq)) {
      
      # matrix M
      m_A = rep(m_df$m_A[j], n_sp)
      m_B = rep(m_df$m_B[j], n_sp)
      M = diag(c(m_A, m_B))
      # matrix G
      g = rep(g_seq[k], n_sp)
      G = rbind(cbind(diag(1 - g), diag(g)),
                cbind(diag(g), diag(1 - g)))
      # identity matrix
      I = diag(rep(1, 2*n_sp))
      
      for (l in 1:n) {
        
        summary_df$m_A[counter] = m_df$m_A[j]
        summary_df$m_B[counter] = m_df$m_B[j]
        summary_df$g[counter] = g_seq[k]
  
        # sampling thetas
        theta_A = runif(n_sp, min = theta_A_min, max = theta_A_max) 
        theta_B = runif(n_sp, min = theta_B_min, max = theta_B_max)
        theta = c(theta_A, theta_B)
      
        # T matrix: no gene flow
        T_matrix = solve(I -  M %*% Q_2sites)
        # spectra of T_matrix: no gene flow
        eigen = Re(eigen(T_matrix)$values)
        # 1st and 2nd eigenvalues
        eigen_1 = eigen[1]
        eigen_2 = eigen[2]
        # analytical solution: no gene flow
        z = as.vector(T_matrix %*% (I - M) %*% theta)
        
        # T matrix: with gene flow
        T_matrix_g = solve(solve(G) -  M %*% Q_2sites)
        # spectra of T_matrix: no gene flow
        eigen_g = Re(eigen(T_matrix_g)$values)
        # 1st and 2nd eigenvalues
        eigen_g_1 = eigen_g[1]
        eigen_g_2 = eigen_g[2]
        # analytical solution with gene flow
        z_g = as.vector(T_matrix_g %*% (I - M) %*% theta)
        
        # difference between eigenvalues
        summary_df$dif_eigen_1[counter] = eigen_g_1 - eigen_1
        summary_df$dif_eigen_2[counter] = eigen_g_2 - eigen_2
        
        # vector z at site A
        z_A = as.vector(z[1:n_sp])
        z_g_A = as.vector(z_g[1:n_sp])
        # vector z at site B
        z_B = as.vector(z[(n_sp+1):(2*n_sp)])
        z_g_B = as.vector(z_g[(n_sp+1):(2*n_sp)])
        
        # cosine of the angle between z vectors at site A
        cos_angle_A = (t(z_A) %*% z_g_A)/(sqrt(sum(z_A^2))*sqrt(sum(z_g_A^2)))
        summary_df$cos_angle_A[counter] = cos_angle_A
        
        # cosine of the angle between z vectors at site B
        cos_angle_B = (t(z_B) %*% z_g_B)/(sqrt(sum(z_B^2))*sqrt(sum(z_g_B^2)))
        summary_df$cos_angle_B[counter] = cos_angle_B
        
        # difference between trait matching of both vectors at site A
        match_z_A = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix,
                                 z = z_A, method = "exponential", alpha = 0.2)
        mean_match_z_A = match_z_A[[1]]
        match_z_g_A = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix,
                                     z = z_g_A, method = "exponential", alpha = 0.2)
        mean_match_z_g_A = match_z_g_A[[1]]
        dif_mean_match_A = mean_match_z_g_A - mean_match_z_A
        summary_df$dif_mean_match_A[counter] = dif_mean_match_A
        
        # difference between trait matching of both vectors at site B
        match_z_B = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix,
                                   z = z_B, method = "exponential", alpha = 0.2)
        mean_match_z_B = match_z_B[[1]]
        match_z_g_B = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix,
                                     z = z_g_B, method = "exponential", alpha = 0.2)
        mean_match_z_g_B = match_z_g_B[[1]]
        dif_mean_match_B = mean_match_z_g_B - mean_match_z_B
        summary_df$dif_mean_match_B[counter] = dif_mean_match_B
        
        counter = counter + 1
      
      }
    }
  }
}

# saving results
write.csv(summary_df, row.names = FALSE, 
          file = "~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/results_comparing_analytical_with_g_vs_no_g.csv")

#-----------------------------------------------------------------------------------------------------#