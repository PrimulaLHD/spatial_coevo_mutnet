#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads trait evolution data of all simulations with a given value of mutualistic selection and gene
#   flow and calculates the standardized degree and pairwise trait matching for each species.
#
# Returns:
#   A csv spreadsheet containing the results.

# loading functions and packages
source("functions/MatchingMutNet.R")

# defining mutualistic selection, gene flow, and site 
m_A = 0.7
m_B = 0.7
g = 0.3
site = "A"
# for the file names
m_A_char = paste("mA", gsub(".", "", as.character(m_A), fixed = TRUE), sep = "")
m_B_char = paste("mB", gsub(".", "", as.character(m_B), fixed = TRUE), sep = "")
g_char = paste("g", gsub(".", "", as.character(g), fixed = TRUE), sep = "")
  
# network names
net_names = dir(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/",
                      "all_networks_", m_A_char, "_", m_B_char, "_g0-03", sep = ""))

final_mut_mat_df_all_nets = data.frame()

for (i in 1:length(net_names)) {
  # simulation file names for the current network
  files = dir(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/",
                    "all_networks_", m_A_char, "_", m_B_char, "_g0-03/", net_names[i], sep = ""))
    
  # extracting g values from the file names
  split1 = sapply(strsplit(files, split = "_site"), "[", 1)
  g_char_all = as.character(sub(".*_g", "", split1))
  g_char_point_all = paste("0", ".", substr(g_char_all, start = 2, stop = nchar(g_char_all)),
                           sep = "") 
  g_all = as.numeric(g_char_point_all)
  # extracting site information from the file names
  split2 = sapply(strsplit(files, split = "_sim"), "[", 1)
  site_all = sub(".*_site", "", split2)
  
  # extracting files for simulations with no gene flow
  files_sub = files[which(g_all == g & site_all == site)]
  
  # reading network
  mat = as.matrix(read.table(paste("data/empirical_networks/binary/", net_names[i], ".txt", sep = ""), 
                             sep = " ", header = FALSE))
  # defining number of rows, columns and species
  n_row = nrow(mat)
  n_col = ncol(mat)
  n_sp = n_row + n_col
  # building the square adjacency matrix f
  f = rbind(cbind(matrix(0, n_row, n_row), mat), 
            cbind(t(mat), matrix(0, n_col, n_col)))
  
  # extracting species identities for interaction pairs
  row_names = paste("R", 1:n_row, sep = "") # codes for rows
  col_names = paste("C", 1:n_col, sep = "") # codes for columns
  rownames(f) = c(row_names, col_names)
  colnames(f) = c(row_names, col_names)
  names_rep = rep(c(row_names, col_names), n_sp) # getting codes for links
  names_node_1 = names_rep[which(f != 0)]
  k = apply(f, 1, sum) # node degrees
  names_node_2 = rep(c(row_names, col_names), k)
  n_links = sum(f)/2 # total number of links
  
  # calculating degree information for column species
  col = names_node_1[1:n_links]
  k_col = k[(n_row+1):(n_row+n_col)]
  k_col_rep = rep(k_col, k_col)
  k_col_rep_ord = k_col_rep[as.character(col)]
  std_k_col = (k_col_rep_ord - mean(k_col_rep_ord))/sd(k_col_rep_ord)
  # defining species topological role
  role_col = rep("specialist", n_links)
  role_col[std_k_col >= 1] = "generalist"
  
  # doing the same for row species
  row = names_node_2[1:n_links]
  k_row = k[1:n_row]
  k_row_rep = rep(k_row, k_row)
  k_row_rep_ord = k_row_rep[as.character(row)]
  std_k_row = (k_row_rep_ord - mean(k_row_rep_ord))/sd(k_row_rep_ord)
  role_row = rep("specialist", n_links)
  role_row[std_k_row > 1] = "generalist"
  
  # putting results in a single data frame
  final_mut_mat_df = data.frame(network = rep(net_names[i], n_links*length(files_sub)),
                                simulation = rep(1:length(files_sub), each = n_links),
                                interaction = rep(1:n_links, times = length(files_sub)),
                                col = rep(names_node_1[1:n_links], times = length(files_sub)),
                                row = rep(names_node_2[1:n_links], times = length(files_sub)),
                                matching = rep(NA, times = length(files_sub)*n_links),
                                k_col = rep(k_col_rep_ord, times = length(files_sub)),
                                std_k_col = rep(std_k_col, times = length(files_sub)),
                                role_col = rep(role_col, times = length(files_sub)),
                                k_row = rep(k_row_rep_ord, times = length(files_sub)),
                                std_k_row = rep(std_k_row, times = length(files_sub)),
                                role_row = rep(role_row, times = length(files_sub)))
  
  # calculating species trait matching
  for (j in 1:length(files_sub)) {
    # reading trait values
    traits = read.csv(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/",
                            "all_networks_", m_A_char, "_", m_B_char, "_g0-03/",
                            net_names[i], "/", files_sub[j], sep = ""), row.names = 1)
    # trait values at equilibrium
    final_traits = as.numeric(traits[3, ])
    # calculating trait matching
    final_mut_mat = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                                   z = final_traits, method = "exponential", alpha = 0.2)
    # paiwise trait matching
    final_mut_mat_pair = final_mut_mat[[2]]
    # adding values to data frame
    final_mut_mat_df$matching[final_mut_mat_df$simulation == j] = final_mut_mat_pair$matching
  }
  # putting results in a larger data frame
  final_mut_mat_df_all_nets = rbind(final_mut_mat_df_all_nets, final_mut_mat_df)
}

# saving results
write.csv(final_mut_mat_df_all_nets, row.names = FALSE, 
          file = paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/",
                       "coevo_results_pairwise_matching_sp_degree_",
                       m_A_char, "_", m_B_char, "_", g_char, "_site", site, ".csv", sep = ""))

#-----------------------------------------------------------------------------------------------------#