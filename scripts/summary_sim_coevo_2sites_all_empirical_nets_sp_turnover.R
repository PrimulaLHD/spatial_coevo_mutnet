#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads coevolution simulation data from all 72 empirical networks with a specific combination of
#   m_A and m_B and many g values and calculates mean and standard deviation of (1) initial and final
#   trait matching, (2) initial and final trait convergence, and (3) initial and final environmental 
#   matching.
#
# Returns:
#   Save the results as csv spreadsheets. 

# loading functions and packages
source("functions/MatchingMutNet.R")
source("functions/ConvMutNet.R")

# defining mutualistic selection values used in simulations
m_A = 0.7
m_B = 0.7
# for the file name
m_A_char = gsub(".", "", as.character(m_A), fixed = TRUE)
m_B_char = gsub(".", "", as.character(m_B), fixed = TRUE)

# defining gene flow range used in simulations
g_min = 0
g_max = 0.3
# for the file name
g_min_char = gsub(".", "", as.character(g_min), fixed = TRUE)
g_max_char = gsub(".", "", as.character(g_max), fixed = TRUE)

# defining folder with simulation results 
folder = paste("~/Lucas/Projects/spatial_coevo_mutnet/output/data/simulations_empirical_networks/",
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
                        init_trait_var = rep(NA, length(net_names)*n_sim),
                        final_trait_var = rep(NA, length(net_names)*n_sim),
                        init_conv = rep(NA, length(net_names)*n_sim),
                        final_conv = rep(NA, length(net_names)*n_sim),
                        mean_init_env_matching = rep(NA, length(net_names)*n_sim),
                        sd_init_env_matching = rep(NA, length(net_names)*n_sim),
                        mean_final_env_matching = rep(NA, length(net_names)*n_sim),
                        sd_final_env_matching = rep(NA, length(net_names)*n_sim))

for (i in 1:length(net_names)) {
  
  print(net_names[i])
  
  # creating vectors to store results
  mean_init_mut_matching = c()
  sd_init_mut_matching = c()
  mean_final_mut_matching = c()
  sd_final_mut_matching = c()
  init_trait_var = c()
  final_trait_var = c()
  init_conv = c()
  final_conv = c()
  mean_init_env_matching = c()
  sd_init_env_matching = c()
  mean_final_env_matching = c()
  sd_final_env_matching = c()
  
  # reading network
  mat = as.matrix(read.table(paste("data/empirical_networks/sensitivity_analysis/", net_names[i], ".txt", sep = ""), 
                             sep = " ", header = FALSE))
  # simulation file names for the current network
  files = dir(paste(folder, "/", net_names[i], sep = ""))
  # extracting g values from the file names
  split1 = sapply(strsplit(files, split = "_alpha"), "[", 1)
  g_char = as.character(sub(".*_g", "", split1))
  g_char_point = paste("0", ".", substr(g_char, start = 2, stop = nchar(g_char)),
                       sep = "") 
  g = round(as.numeric(g_char_point), digits = 3)
  # extracting site information from the file names
  split2 = sapply(strsplit(files, split = "_sim"), "[", 1)
  site = sub(".*_site", "", split2)
  # extracting information on additional parameters
  split3 = sapply(strsplit(files, split = "_site"), "[", 1)
  add_par = paste("_alpha", as.character(sub(".*_alpha", "", split3)), sep = "")
  add_par_uni = unique(add_par)
  
  # calculating mutualistic and environmental matching
  
  alpha = 0.2
  
  for (j in 1:length(files)) {
    # building the square adjacency matrix f
    f = rbind(cbind(matrix(0, nrow(mat), nrow(mat)), mat), 
              cbind(t(mat), matrix(0, ncol(mat), ncol(mat))))
    # adding species labels to the matrix
    sp_labels = c(paste("R", 1:nrow(mat), sep = ""), paste("C", 1:ncol(mat), sep = ""))
    rownames(f) = sp_labels
    colnames(f) = sp_labels
    # reading simulation results
    traits = read.csv(paste(folder, "/", net_names[i], "/", files[j], sep = ""), row.names = 1)
    thetas = as.numeric(traits[1, ])
    init_traits = as.numeric(traits[2, ])
    final_traits = as.numeric(traits[3, ])
    # obtaining interaction matrix
    f = f[names(traits), names(traits)]
    # obtaining number of rows and columns of bipartite matrix
    n_row = length(grep("R", names(traits)))
    n_col = length(grep("C", names(traits)))
    n_sp = n_row + n_col
    # initial and final trait matching
    init_mut_mat = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                                  z = init_traits, method = "exponential", alpha = alpha)
    mean_init_mut_matching[j] = init_mut_mat[[1]]
    sd_init_mut_matching[j] = sd(init_mut_mat[[2]]$matching)
    final_mut_mat = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                                   z = final_traits, method = "exponential", alpha = alpha)
    mean_final_mut_matching[j] = final_mut_mat[[1]]
    sd_final_mut_matching[j] = sd(final_mut_mat[[2]]$matching)
    # initial and final trait convergence
    init_trait_var_both = ConvMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, 
                                     z = init_traits, method = "nuismer", alpha = alpha)
    init_trait_var[j] = mean(init_trait_var_both) 
    final_trait_var_both = ConvMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, 
                                      z = final_traits, method = "nuismer", alpha = alpha)
    final_trait_var[j] = mean(final_trait_var_both)
    init_conv_both = ConvMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, 
                                z = init_traits, method = "exponential", alpha = alpha)
    init_conv[j] = mean(init_conv_both)
    final_conv_both = ConvMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, 
                                 z = final_traits, method = "exponential", alpha = alpha)
    final_conv[j] = mean(final_conv_both)
    # initial and final environmental matching
    init_env_matching = exp(-alpha * (init_traits - thetas)^2)
    mean_init_env_matching[j] = mean(init_env_matching)
    sd_init_env_matching[j] = sd(init_env_matching)  
    final_env_matching = exp(-alpha * (final_traits - thetas)^2)
    mean_final_env_matching[j] = mean(final_env_matching)
    sd_final_env_matching[j] = sd(final_env_matching)
  }
  
  # adding g, site, trait matching, environmental matching and divergence for current network
  summary_df[summary_df$network == net_names[i], "g"] = g
  summary_df[summary_df$network == net_names[i], "site"] = site
  summary_df[summary_df$network == net_names[i], "mean_init_mut_matching"] = mean_init_mut_matching
  summary_df[summary_df$network == net_names[i], "sd_init_mut_matching"] = sd_init_mut_matching
  summary_df[summary_df$network == net_names[i], "mean_final_mut_matching"] = mean_final_mut_matching
  summary_df[summary_df$network == net_names[i], "sd_final_mut_matching"] = sd_final_mut_matching
  summary_df[summary_df$network == net_names[i], "init_trait_var"] = init_trait_var
  summary_df[summary_df$network == net_names[i], "final_trait_var"] = final_trait_var
  summary_df[summary_df$network == net_names[i], "init_conv"] = init_conv
  summary_df[summary_df$network == net_names[i], "final_conv"] = final_conv
  summary_df[summary_df$network == net_names[i], "mean_init_env_matching"] = mean_init_env_matching
  summary_df[summary_df$network == net_names[i], "sd_init_env_matching"] = sd_init_env_matching
  summary_df[summary_df$network == net_names[i], "mean_final_env_matching"] = mean_final_env_matching
  summary_df[summary_df$network == net_names[i], "sd_final_env_matching"] = sd_final_env_matching
  
}

# defining folder to store summary spreadsheet
folder_results = "~/Lucas/Projects/spatial_coevo_mutnet/output/data/simulations_empirical_networks/summary_coevo_results/"

# saving results
write.csv(summary_df, row.names = FALSE, 
          file = paste(folder_results, "summary_coevo_results",
                       "_mA", m_A_char, "_mB", m_B_char,
                       "_g", g_min_char, "-", g_max_char,
                       add_par_uni, ".csv", sep = ""))

#-----------------------------------------------------------------------------------------------------#