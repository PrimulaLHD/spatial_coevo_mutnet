#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Computes difference in several metrics using analytical equilibrium values with and without gene flow.
#
# Returns:
#   Each calculated metric.

source("functions/MatchingMutNet.R")

# defining a matrix
mat_file = "PP_Olesen_etal_2002_Aigrettes_wgt.txt"

# matrix F and Q
mat = as.matrix(read.table(paste("data/empirical_networks/weighted/", mat_file, sep = ""), 
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

# matrix M
m_A = rep(0.7, n_sp)
m_B = rep(0.3, n_sp)
M = diag(c(m_A, m_B))

# matrix G
g = rep(0.2, n_sp)
G = rbind(cbind(diag(1 - g), diag(g)),
          cbind(diag(g), diag(1 - g)))

# identity matrix
I = diag(rep(1, 2*n_sp))

# sampling thetas
theta_A = runif(n_sp, min = 0, max = 10) 
theta_B = runif(n_sp, min = 10, max = 20)
theta = c(theta_A, theta_B)

# T matrix: no gene flow
T_matrix = solve(I -  M %*% Q_2sites)
# spectra of T_matrix: no gene flow
eigen = Re(eigen(T_matrix)$values)
# 1st and 2nd eigenvalues
eigen[1]
eigen[2]
# analytical solution: no gene flow
z = T_matrix %*% (I - M) %*% theta

# T matrix: with gene flow
T_matrix_g = solve(solve(G) -  M %*% Q_2sites)
# spectra of T_matrix: no gene flow
eigen_g = Re(eigen(T_matrix_g)$values)
# 1st and 2nd eigenvalues
eigen_g[1]
eigen_g[2]
# analytical solution with gene flow
z_g = T_matrix_g %*% (I - M) %*% theta

# comparing z vectors
(cos_angle = (t(z) %*% z_g)/(sqrt(sum(z^2))*sqrt(sum(z_g^2))))

# comparing trait matching of both vectors at site A
z_A = as.vector(z[1:n_sp])
z_g_A = as.vector(z_g[1:n_sp])

match_z = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix,
                         z = z_A, method = "exponential", alpha = 0.2)
mean_match_z = match_z[[1]]

match_z_g = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = F_matrix,
                           z = z_g_A, method = "exponential", alpha = 0.2)
mean_match_z_g = match_z_g[[1]]

(dif_mean_match = mean_match_z_g - mean_match_z)

#-----------------------------------------------------------------------------------------------------#