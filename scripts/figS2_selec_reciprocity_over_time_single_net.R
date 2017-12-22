#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads trait evolution through time data for one network, computes mean reciprocal selection and
#   plots reciprocal selection through time for different values of mutualistic selection
#
# Returns:
#   A plot.

# loading functions and packages
library(ggplot2)
library(cowplot)
library(viridis)
source("functions/MatchingMutNet.R")

# defining network
network = "FP_Galetti&Pizo_1996_bin"

# reading network
mat = as.matrix(read.table(paste("data/empirical_networks/binary/", network, ".txt", sep = ""), 
                           sep = " ", header = FALSE))
# defining number of rows, columns and species
n_row = nrow(mat)
n_col = ncol(mat)
n_sp = n_row + n_col
# building the square adjacency matrix f
f = rbind(cbind(matrix(0, n_row, n_row), mat), 
          cbind(t(mat), matrix(0, n_col, n_col)))

# defining m values
mA = c(0.1, 0.3, 0.5, 0.7, 0.9)
mB = 0.1
# defining mutualistic selection scenario
mA_char = paste("mA", gsub(".", "", mA, fixed = TRUE), sep = "")
mB_char = paste("mB", gsub(".", "", mB, fixed = TRUE), sep = "")
m_char = paste(mA_char, mB_char, sep = "_")
# defining gene flow value
g = 0
g_char = paste("g", g, "_", sep = "")
# defining site
site_char = "siteA"

# alpha value
alpha = 0.2

# vector to store mutualistic selection
m = c()
# vector to store times to equilibrium
time = c()
# vector to store reciprocal selection through time
mean_rec_sel = c()
# vector to store simulation number
simulation = c()

for (i in 1:length(m_char)) {
  # obtaining simulation files
  folders = dir("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_all_timesteps/")
  folder = folders[grep(m_char[i], folders)]
  files = dir(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_all_timesteps/",
                    folder, "/", network, sep = ""))
  files_g = files[grep(g_char, files)]
  files_g_site = files_g[grep(site_char, files_g)]
  
  for (j in 1:length(files_g_site)) {
    # read trait values through time
    traits_df = read.csv(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_all_timesteps/",
                               folder, "/", network, "/", files_g_site[j], sep = ""))
    # removing 1st column with labels
    traits_df = traits_df[ , 2:ncol(traits_df)]
    # theta values
    thetas = traits_df[1, ]
    # removing theta values
    traits_df = traits_df[2:nrow(traits_df), ]
    # time to equilibrium
    time = c(time, 1:nrow(traits_df))
    for (k in 1:nrow(traits_df)) {
      # traits at time j
      z = as.numeric(traits_df[k, ])
      # reciprocal selection
      z_dif = t(f*z) - f*z # matrix with all trait differences
      q = f*(exp(-alpha * (z_dif^2))) # calculating matrix q
      q_n = q / apply(q, 1, sum) # normalizing the matrix
      q_n_vec = q_n[f != 0]
      t_q_n_vec = t(q_n)[f != 0]
      rec_vec = log(q_n_vec) + log(t_q_n_vec) # calculating log of reciprocity
      mean_rec_sel = c(mean_rec_sel, mean(rec_vec)) # mean reciprocity for the network
    }
    
    simulation = c(simulation, rep(j, times = nrow(traits_df)))
    m = c(m, rep(mA[i], times = nrow(traits_df)))
    
  }
  
  results_df = data.frame(m = m,
                          simulation = simulation,
                          time = time,
                          mean_rec_sel = mean_rec_sel)
  
}

ggplot(data = subset(results_df, time < 30), 
       mapping = aes(x = time, y = mean_rec_sel, group = simulation)) +
  geom_line(size = 1, color = "black", alpha = 0.5) +
  xlab("Time") +
  ylab("Mean selection reciprocity (log)") +
  facet_grid( ~ m, labeller = label_both, scales = "free") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))

#-----------------------------------------------------------------------------------------------------#