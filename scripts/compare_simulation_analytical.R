#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Compares the analytical expression for the equilibrium point of the coevolutionary model with the
#   simulation results for a random matrix.
#
# Returns:
#   A plot of the two equilibrium vectors and a correlation value.

source("functions/CoevoMutNet2Sites.R")
source("functions/RemoveNodesNoLinks.R")
source("functions/CompMutNet.R")

# creating adjacency matrix
n_row = 50
n_col = 50
n_sp = n_row + n_col
mat = matrix(data = sample(c(0, 1), size = n_row*n_col, replace = TRUE, prob = c(0.7, 0.3)),
             nrow = n_row, ncol = n_col)
mat = RemoveNodesNoLinks(mat)
n_row = nrow(mat)
n_col = ncol(mat)
n_sp = n_row + n_col
f = rbind(cbind(matrix(0, n_row, n_row), mat),
          cbind(t(mat), matrix(0, n_col, n_col)))

# defining simulation parameters
alpha = 0.2
# herdability
h_mean = 0.2
h_sd = 0.01
h = rnorm(n_sp, h_mean, h_sd) 
while (any(h < 0 | h > 1)) 
  h = rnorm(n_sp, h_mean, h_sd)
# mutualistic selection 
m_A_mean = 0.8
m_A_sd = 0.01
m_B_mean = 0.8
m_B_sd = 0.01
m_A = rnorm(n_sp, m_A_mean, m_A_sd)
while (any(m_A < 0 | m_A > 1)) 
  m_A = rnorm(n_sp, m_A_mean, m_A_sd)
m_B = rnorm(n_sp, m_B_mean, m_B_sd)
while (any(m_B < 0 | m_B > 1)) 
  m_B = rnorm(n_sp, m_B_mean, m_B_sd)
# gene flow
g_mean = 0.2
g_sd = 0.001
g = rnorm(n_sp, g_mean, g_sd)
while (any(g < 0 | g > 1)) 
  g = rnorm(n_sp, g_mean, g_sd)
# environmental optimum
theta_A_min = 0
theta_A_max = 10
theta_B_min = 10
theta_B_max = 20
theta_A = runif(n_sp, min = theta_A_min, max = theta_A_max) 
theta_B = runif(n_sp, min = theta_B_min, max = theta_B_max)
# initial conditions
init_A = runif(n_sp, min = theta_A_min, max = theta_A_max) 
init_B = runif(n_sp, min = theta_B_min, max = theta_B_max)
# tolerance to stop simulation
epsilon = 0.000001
# maximum number of generations
t_max = 10000

# running simulation 
z_list = CoevoMutNet2Sites(n_sp = n_sp, f = f, g = g, h = h, alpha = alpha, 
                           theta_A = theta_A, theta_B = theta_B,
                           init_A = init_A, init_B = init_B,
                           m_A = m_A, m_B = m_B, epsilon = epsilon, t_max = t_max)

# matrix Q at equilibrium
z_A = z_list[[1]][nrow(z_list[[1]]), ]
z_B = z_list[[2]][nrow(z_list[[2]]), ]
z_dif_A = t(f*z_A) - f*z_A
z_dif_B = t(f*z_B) - f*z_B
q_A = f*(exp(-alpha * (z_dif_A^2)))
q_B = f*(exp(-alpha * (z_dif_B^2)))
q_n_A = q_A / apply(q_A, 1, sum)
q_n_B = q_B / apply(q_B, 1, sum)
Q = rbind(cbind(q_n_A, matrix(0, n_sp, n_sp)),
          cbind(matrix(0, n_sp, n_sp), q_n_B))
# vector z
z = c(z_A, z_B)
# vector theta
theta = c(theta_A, theta_B)
# matrix M
M = diag(c(m_A, m_B))
# gene flow matrix
G = rbind(cbind(diag(1 - g), diag(g)),
          cbind(diag(g), diag(1 - g)))
# matrix Phi
Phi = diag(rep(h, times = 2))
# identity matrix
I = diag(rep(1, 2*n_sp))

# calculting solution using the analytical expression
z_analytic = solve(solve(G) - I + Phi %*% (I - M %*% Q)) %*% Phi %*% (I - M) %*% theta

# comparing the simulation results with analytical solution
plot(z, c(z_analytic))
abline(a = 0, b = 1, col = "red")
cor(z, c(z_analytic))

#-----------------------------------------------------------------------------------------------------#