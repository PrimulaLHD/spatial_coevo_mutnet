#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Runs n_sim simulations per empirical network for several different mean gene flow values (g) using
#   a given choice of parameters.
#
# Returns:
#   Saves two csv files (one for each site) containing the environmental optimum values in the first row
#   and all trait values through time in the subsequent rows.

# loading functions
source("functions/CoevoMutNet2Sites.R")

# defining number of simulations per network
n_sim = 50
# defining alpha
alpha = 0.2
# defining phi distribution
phi_mean = 0.5
phi_sd = 0.01
# defining mutualism selection distribution
m_A_mean = 0.7
m_A_sd = 0.01
m_B_mean = 0.1
m_B_sd = 0.01
# defining environmental optimum range
theta_A_min = 0
theta_A_max = 10
theta_B_min = 10
theta_B_max = 20
# defining gene flow range and distribution
g_min = 0
g_max = 0.3
g_dif = 0.025
g_means = seq(g_min, g_max, by = g_dif)
g_sd = 0.001
# for the file name
alpha_char = gsub(".", "", as.character(alpha), fixed = TRUE)
phi_char = gsub(".", "", as.character(phi_mean), fixed = TRUE)
m_A_char = gsub(".", "", as.character(m_A_mean), fixed = TRUE)
m_B_char = gsub(".", "", as.character(m_B_mean), fixed = TRUE)
theta_A_min_char = gsub(".", "", as.character(theta_A_min), fixed = TRUE)
theta_A_max_char = gsub(".", "", as.character(theta_A_max), fixed = TRUE)
theta_B_min_char = gsub(".", "", as.character(theta_B_min), fixed = TRUE)
theta_B_max_char = gsub(".", "", as.character(theta_B_max), fixed = TRUE)
g_min_char = gsub(".", "", as.character(g_min), fixed = TRUE)
g_max_char = gsub(".", "", as.character(g_max), fixed = TRUE)
# defining tolerace to stop simulation
epsilon = 0.000001
# defining maximum number of generations
t_max = 10000

# defining folder to store results
folder = "~/simulations/"

# creating folder to store results
dir.create(path = paste(folder, "all_networks", "_mA", m_A_char, "_mB", m_B_char, "_g",
                        g_min_char, "-", g_max_char, sep = ""))

# all network files names
net_files = list.files("data/empirical_networks/binary/")
# removing the extension .txt
net_names = substr(net_files, start = 1, stop = nchar(net_files) - 4) 

for (i in 1:length(net_files)) {
  # creating directory to store results from this network
  dir.create(path = paste(folder, "all_networks", "_mA", m_A_char, "_mB", m_B_char, "_g", 
                          g_min_char, "-", g_max_char, "/", net_names[i], sep = ""))
  # reading network
  mat = as.matrix(read.table(paste("data/empirical_networks/binary/", net_files[i], sep = ""), 
                             sep = " ", header = FALSE))
  # defining number of rows, columns and species
  n_row = nrow(mat)
  n_col = ncol(mat)
  n_sp = n_row + n_col
  # building the square adjacency matrix f
  f = rbind(cbind(matrix(0, n_row, n_row), mat),
            cbind(t(mat), matrix(0, n_col, n_col)))
  
  # simulations
  for (j in 1:length(g_means)) {
    for (k in 1:n_sim) {
      # sampling phi values
      phi_A = rnorm(n_sp, phi_mean, phi_sd) 
      while (any(phi_A < 0)) 
        phi_A = rnorm(n_sp, phi_mean, phi_sd) 
      phi_B = rnorm(n_sp, phi_mean, phi_sd) 
      while (any(phi_B < 0)) 
        phi_B = rnorm(n_sp, phi_mean, phi_sd) 
      # sampling mutualism selection values
      m_A = rnorm(n_sp, m_A_mean, m_A_sd)
      while (any(m_A < 0 | m_A > 1)) 
        m_A = rnorm(n_sp, m_A_mean, m_A_sd)
      m_B = rnorm(n_sp, m_B_mean, m_B_sd)
      while (any(m_B < 0 | m_B > 1)) 
        m_B = rnorm(n_sp, m_B_mean, m_B_sd)
      # sampling thetas
      theta_A = runif(n_sp, min = theta_A_min, max = theta_A_max) 
      theta_B = runif(n_sp, min = theta_B_min, max = theta_B_max)
      # sampling initial conditions
      init_A = runif(n_sp, min = theta_A_min, max = theta_A_max) 
      init_B = runif(n_sp, min = theta_B_min, max = theta_B_max)
      # sampling gene flow
      if (g_means[j] == 0) {
        # gene flow is the same for every species if it is 0
        g = rep(0, n_sp)
      } else {
        g = rnorm(n_sp, g_means[j], g_sd)
        while (any(g < 0 | g > 1)) 
          g = rnorm(n_sp, g_means[j], g_sd)
      }
      # for file name
      g_char = gsub(".", "", as.character(g_means[j]), fixed = TRUE)
      
      # running simulation
      z_list = CoevoMutNet2Sites(n_sp = n_sp, f = f, g = g, alpha = alpha, phi_A = phi_A, phi_B = phi_B,
                                 theta_A = theta_A, theta_B = theta_B, init_A = init_A, init_B = init_B, 
                                 m_A = m_A, m_B = m_B, epsilon = epsilon, t_max = t_max)
      
      # traits for each site
      z_A = z_list[[1]]
      z_B = z_list[[2]]
      
      # building data frames with results
      site_A = rbind(theta_A, z_A)
      colnames(site_A) = c(paste("R", 1:n_row, sep = ""), paste("C", 1:n_col, sep = ""))
      rownames(site_A) = c("theta_A", 1:nrow(z_A))
      site_B = rbind(theta_B, z_B)
      colnames(site_B) = c(paste("R", 1:n_row, sep = ""), paste("C", 1:n_col, sep = ""))
      rownames(site_B) = c("theta_B", 1:nrow(z_B))
      
      # saving species initial and final trait values 
      write.csv(site_A, file = paste(folder, "all_networks", "_mA", m_A_char, "_mB", m_B_char, "_g", 
                                     g_min_char,"-", g_max_char, "/", net_names[i], "/", net_names[i],
                                     "_mA", m_A_char, "_mB", m_B_char, "_g", g_char,
                                     "_alpha", alpha_char, "_phi", phi_char, 
                                     "_thetaA", theta_A_min_char, "-", theta_A_max_char,
                                     "_thetaB", theta_B_min_char, "-", theta_B_max_char,
                                     "_siteA", "_sim", k, ".csv", sep = ""))
      write.csv(site_B, file = paste(folder, "all_networks", "_mA", m_A_char, "_mB", m_B_char, "_g", 
                                     g_min_char,"-", g_max_char, "/", net_names[i], "/", net_names[i],
                                     "_mA", m_A_char, "_mB", m_B_char, "_g", g_char,
                                     "_alpha", alpha_char, "_phi", phi_char, 
                                     "_thetaA", theta_A_min_char, "-", theta_A_max_char,
                                     "_thetaB", theta_B_min_char, "-", theta_B_max_char,
                                     "_siteB", "_sim", k, ".csv", sep = ""))
    }
  }
}

#-----------------------------------------------------------------------------------------------------#