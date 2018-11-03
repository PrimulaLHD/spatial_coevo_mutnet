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
folder = "output/data/simulations_empirical_networks/simulations_geneflow_removal"

# define m_A and m_B sequences
m_df = expand.grid(seq(0.1, 0.9, 0.2), seq(0.1, 0.9, 0.2))
m_df = m_df[m_df[ , 1] >= m_df[ , 2], ]
names(m_df) = c("mA", "mB")

# define gene flow values
g_high = 0.3
g_low = 0
# for the file name
g_high_char = gsub(".", "", as.character(g_high), fixed = TRUE)
g_low_char = gsub(".", "", as.character(g_low), fixed = TRUE)
# define perturbation sequence
pert_g_seq = seq(0, 1, by = 0.05)
# defining theta distribution
theta_A_min = 0
theta_A_max = 10
theta_B_min = 10
theta_B_max = 20
# number of theta vectors 
n_theta = 10
# number of replicates of perturbation simulations
n_rep_pert = 10

for(i in 1:length(net_files)) {
  # read empirical matrix
  mat = as.matrix(read.table(paste("data/empirical_networks/weighted/", net_files[i], sep = ""), 
                             sep = " ", header = FALSE))
  n_row = nrow(mat)
  n_col = ncol(mat)
  n_sp = n_row + n_col
  
  results_df = data.frame(network = rep(net_names[i], times = n_rep_pert*n_theta*nrow(m_df)*length(pert_g_seq)),
                          rep_pert = rep(rep(1:n_rep_pert, each = length(pert_g_seq)*n_theta), nrow(m_df)),
                          rep_theta = rep(1:n_theta, times = n_rep_pert*nrow(m_df)*length(pert_g_seq)),
                          pert = rep(rep(pert_g_seq, each = n_theta), times = n_rep_pert*nrow(m_df)),
                          m_A = rep(rep(m_df$mA, each = length(pert_g_seq)), each = n_rep_pert*n_theta),
                          m_B = rep(rep(m_df$mB, each = length(pert_g_seq)), each = n_rep_pert*n_theta),
                          matching_A = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          matching_B = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          dim_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          sum_eigen_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          mean_eigen_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          sd_eigen_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          eigen1_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          eigen2_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          det_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          mean_row_corr_T = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta),
                          sd_g = rep(NA, times = nrow(m_df)*length(pert_g_seq)*n_rep_pert*n_theta))

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
  
  for (j in 1:nrow(m_df)) {
    
    # matrix M
    m_A = rep(m_df$mA[j], n_sp)
    m_B = rep(m_df$mB[j], n_sp)
    M = diag(c(m_A, m_B))
    
    for (k in 1:n_rep_pert) {
      
      # sampling species perturbation order
      sp_order = sample(1:n_sp, replace = FALSE)
      
      for (l in 1:length(pert_g_seq)) {
        # number of species that will loose gene flow
        sp_number = round(pert_g_seq[l]*n_sp)
        # creating gene flow vector 
        g_mult = rep(g_high, n_sp)
        
        if (sp_number != 0) {
          # indexes of species that will loose gene flow
          pert_sp = sp_order[1:sp_number]
          # modifying gene flow vector
          g_mult[pert_sp] = g_low
        }
        # creating G matrix
        G = rbind(cbind(diag(1 - g_mult), diag(g_mult)),
                  cbind(diag(g_mult), diag(1 - g_mult)))
        
        # standard deviation of gene flow values
        sd_g = sd(g_mult)
        
        # T matrix
        T_matrix = solve(solve(G) -  M %*% Q_2sites) %*% (I - M)
        # dimension T
        dim_T = nrow(T_matrix)
        # T matrix eigenvalues
        eigen_T = Re(eigen(T_matrix)$values)
        # sum of eigenvalues
        sum_eigen_T = sum(eigen_T)
        # mean of eigenvalues
        mean_eigen_T = mean(eigen_T)
        # sd of eigenvalues
        sd_eigen_T = sd(eigen_T)
        # 1st and 2nd eigenvalues
        eigen1_T = eigen_T[1]
        eigen2_T = eigen_T[2]
        # determinant
        det_T = prod(eigen_T)
        
        # function to obtain vector norm
        norm = function(x) sqrt(sum(x^2))
        # normalize T matrix
        T_norm = T_matrix/as.vector(apply(T_matrix, 1, norm))
        # mean correlation among rows
        corr_mat = T_norm %*% t(T_norm)
        mean_row_corr_T = mean(corr_mat[lower.tri(corr_mat)])
        
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
          results_df$matching_A[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                  results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                  results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = mean_mat_A
          
          results_df$matching_B[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                  results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                  results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = mean_mat_B
          
          results_df$dim_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                             results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                             results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = dim_T
          
          results_df$sum_eigen_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                   results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                   results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = sum_eigen_T
          
          results_df$mean_eigen_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                    results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                    results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = mean_eigen_T
          
          results_df$sd_eigen_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                  results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                  results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = sd_eigen_T
          
          results_df$eigen1_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = eigen1_T
          
          results_df$eigen2_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = eigen2_T
          
          results_df$det_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                             results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                             results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = det_T
          
          results_df$mean_row_corr_T[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                                     results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                                     results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = mean_row_corr_T
          
          results_df$sd_g[results_df$network ==  net_names[i] & results_df$m_A ==  m_df$mA[j] & 
                            results_df$m_B ==  m_df$mB[j] & results_df$rep_pert == k & 
                            results_df$pert == pert_g_seq[l] & results_df$rep_theta == t] = sd_g
           
        }
      }
    }
  }
  
  # saving results
  write.csv(results_df, file = paste(folder, "/", net_names[i], "_geneflow_removal_g", 
                                     g_high_char, "-", g_low_char, ".csv", sep = ""), row.names = FALSE)
  
}

#-----------------------------------------------------------------------------------------------------#