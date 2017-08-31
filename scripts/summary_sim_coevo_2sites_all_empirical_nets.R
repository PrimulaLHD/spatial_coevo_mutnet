#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads coevolution simulation data from all 72 empirical networks with a specific combination of
#   m_A and m_B and many g values and calculates mean and standard deviation of (1) initial and final
#   trait matching, (2) initial and final environmental matching, and (3) initial and final geographic 
#   divergence.
#
# Returns:
#   Save the results as csv spreadsheets. 

# loading functions and packages
source("functions/MatchingMutNet.R")

# defining mutualistic selection values used in simulations
m_A = 0.9
m_B = 0.9
# for the file name
m_A_char = gsub(".", "", as.character(m_A), fixed = TRUE)
m_B_char = gsub(".", "", as.character(m_B), fixed = TRUE)

# defining gene flow range used in simulations
g_min = 0
g_max = 0.1
# for the file name
g_min_char = gsub(".", "", as.character(g_min), fixed = TRUE)
g_max_char = gsub(".", "", as.character(g_max), fixed = TRUE)

# defining folder with simulation results 
folder = paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/",
               "all_networks", "_mA", m_A_char, "_mB", m_B_char, "_g",
               g_min_char, "-", g_max_char, sep = "")

# network names
net_names = dir(folder)

# number of simulations per network
n_sim = length(dir(paste(folder, "/", net_names[1], sep = "")))

# creating data.frame to store results
summary_df = data.frame(network = rep(net_names, each = n_sim),
                        mutualism = rep(substr(net_names, start = 1, stop = 2), each = n_sim),
                        site = rep(NA, length(net_names)*n_sim),
                        g = rep(NA, length(net_names)*n_sim),
                        mean_init_mut_matching = rep(NA, length(net_names)*n_sim),
                        sd_init_mut_matching = rep(NA, length(net_names)*n_sim),
                        mean_final_mut_matching = rep(NA, length(net_names)*n_sim),
                        sd_final_mut_matching = rep(NA, length(net_names)*n_sim),
                        mean_init_env_matching = rep(NA, length(net_names)*n_sim),
                        sd_init_env_matching = rep(NA, length(net_names)*n_sim),
                        mean_final_env_matching = rep(NA, length(net_names)*n_sim),
                        sd_final_env_matching = rep(NA, length(net_names)*n_sim),
                        mean_init_divergence = rep(NA, length(net_names)*n_sim),
                        sd_init_divergence = rep(NA, length(net_names)*n_sim),
                        mean_final_divergence = rep(NA, length(net_names)*n_sim),
                        sd_final_divergence = rep(NA, length(net_names)*n_sim))

for (i in 1:length(net_names)) {
  
  print(net_names[i])
  
  # creating vectors to store results
  mean_init_mut_matching = c()
  sd_init_mut_matching = c()
  mean_final_mut_matching = c()
  sd_final_mut_matching = c()
  mean_init_env_matching = c()
  sd_init_env_matching = c()
  mean_final_env_matching = c()
  sd_final_env_matching = c()
  mean_init_divergence = c()
  sd_init_divergence = c()
  mean_final_divergence = c()
  sd_final_divergence = c()
  # reading network
  mat = as.matrix(read.table(paste("data/empirical_networks/", net_names[i], ".txt", sep = ""), 
                             sep = " ", header = FALSE))
  # defining number of rows, columns and species
  n_row = nrow(mat)
  n_col = ncol(mat)
  n_sp = n_row + n_col
  # building the square adjacency matrix f
  f = rbind(cbind(matrix(0, n_row, n_row), mat), 
            cbind(t(mat), matrix(0, n_col, n_col)))
  # simulation file names for the current network
  files = dir(paste(folder, "/", net_names[i], sep = ""))
  # extracting g values from the file names
  split1 = sapply(strsplit(files, split = "_site"), "[", 1)
  g_char = as.character(sub(".*_g", "", split1))
  g_char_point = paste("0", ".", substr(g_char, start = 2, stop = nchar(g_char)),
                       sep = "") 
  g = round(as.numeric(g_char_point), digits = 3)
  # extracting site information from the file names
  split2 = sapply(strsplit(files, split = "_sim"), "[", 1)
  site = sub(".*_site", "", split2)
  
  # calculating mutualistic and environmental matching, and divergence
  
  alpha = 0.2
  
  for (j in 1:length(files)) {
    traits = read.csv(paste(folder, "/", net_names[i], "/", files[j], sep = ""), row.names = 1)
    thetas = as.numeric(traits[1, ])
    init_traits = as.numeric(traits[2, ])
    final_traits = as.numeric(traits[3, ])
    # initial and final trait matching
    init_mut_mat = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                                  z = init_traits, method = "exponential", alpha = alpha)
    mean_init_mut_matching[j] = init_mut_mat[[1]]
    sd_init_mut_matching[j] = sd(init_mut_mat[[2]]$matching)
    final_mut_mat = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                                   z = final_traits, method = "exponential", alpha = alpha)
    mean_final_mut_matching[j] = final_mut_mat[[1]]
    sd_final_mut_matching[j] = sd(final_mut_mat[[2]]$matching)
    # initial and final environmental matching
    init_env_matching = exp(-alpha * (init_traits - thetas)^2)
    mean_init_env_matching[j] = mean(init_env_matching)
    sd_init_env_matching[j] = sd(init_env_matching)  
    final_env_matching = exp(-alpha * (final_traits - thetas)^2)
    mean_final_env_matching[j] = mean(final_env_matching)
    sd_final_env_matching[j] = sd(final_env_matching)
    
    if (grepl("siteA", files[j])) { # if the simulation results are from site A
      file_B = gsub(pattern = "siteA", replacement = "siteB", x = files[j])
      traits_B = read.csv(paste(folder, "/", net_names[i], "/", file_B, sep = ""), row.names = 1)
      thetas_B = as.numeric(traits_B[1, ])
      init_traits_B = as.numeric(traits_B[2, ])
      final_traits_B = as.numeric(traits_B[3, ])
      # initial and final geographical divergence 
      init_divergence = 1 - exp(-alpha * (init_traits - init_traits_B)^2)
      mean_init_divergence[j] = mean(init_divergence)
      sd_init_divergence[j] = sd(init_divergence)
      mean_init_divergence[which(files == file_B)] = mean_init_divergence[j] # adding results for site B as well
      sd_init_divergence[which(files == file_B)] = sd_init_divergence[j]
      final_divergence = 1 - exp(-alpha * (final_traits - final_traits_B)^2)
      mean_final_divergence[j] = mean(final_divergence)
      sd_final_divergence[j] = sd(final_divergence)
      mean_final_divergence[which(files == file_B)] = mean_final_divergence[j] # adding results for site B as well
      sd_final_divergence[which(files == file_B)] = sd_final_divergence[j]
    }
    
  }
  
  # adding g, site, trait matching, environmental matching and divergence for current network
  summary_df[summary_df$network == net_names[i], "g"] = g
  summary_df[summary_df$network == net_names[i], "site"] = site
  summary_df[summary_df$network == net_names[i], "mean_init_mut_matching"] = mean_init_mut_matching
  summary_df[summary_df$network == net_names[i], "sd_init_mut_matching"] = sd_init_mut_matching
  summary_df[summary_df$network == net_names[i], "mean_final_mut_matching"] = mean_final_mut_matching
  summary_df[summary_df$network == net_names[i], "sd_final_mut_matching"] = sd_final_mut_matching
  summary_df[summary_df$network == net_names[i], "mean_init_env_matching"] = mean_init_env_matching
  summary_df[summary_df$network == net_names[i], "sd_init_env_matching"] = sd_init_env_matching
  summary_df[summary_df$network == net_names[i], "mean_final_env_matching"] = mean_final_env_matching
  summary_df[summary_df$network == net_names[i], "sd_final_env_matching"] = sd_final_env_matching
  summary_df[summary_df$network == net_names[i], "mean_init_divergence"] = mean_init_divergence
  summary_df[summary_df$network == net_names[i], "sd_init_divergence"] = sd_init_divergence
  summary_df[summary_df$network == net_names[i], "mean_final_divergence"] = mean_final_divergence
  summary_df[summary_df$network == net_names[i], "sd_final_divergence"] = sd_final_divergence
  
}

# defining folder to store summary spreadsheet
folder_results = "~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/"

# saving results
write.csv(summary_df, row.names = FALSE, 
          file = paste(folder_results, "summary_coevo_results",
                       "_mA", m_A_char, "_mB", m_B_char,
                       "_g", g_min_char, "-", g_max_char, ".csv", sep = ""))

#-----------------------------------------------------------------------------------------------------#