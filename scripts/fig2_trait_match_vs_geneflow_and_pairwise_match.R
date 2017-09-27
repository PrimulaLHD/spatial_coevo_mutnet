#-----------------------------------------------------------------------------------------------------#

# Description: 
#   
#
# Returns:
#   

# loading functions and packages
library(ggplot2)
library(cowplot)
library(viridis)
source("functions/MatchingMutNet.R")

# reading data
summ_final_mut_mat_df = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/summary_coevo_results_final_mut_matching.csv")

# changing mutualism names
summ_final_mut_mat_df$mutualism = as.character(summ_final_mut_mat_df$mutualism)
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AE"] = "ant-EFN"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AM"] = "ant-myrmecophyte"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "CC"] = "cleaner-client"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FA"] = "fish-anemone"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FP"] = "frugivore-plant"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "PP"] = "pollinator-plant"

#-------------------------------------------------------------------------#
# Fig. 2A: trait matching vs gene flow for one network with mA = mB = 0.7 #
#-------------------------------------------------------------------------#

# mean final trait matching at site A vs gene flow for m_A = m_B = 0.7
p_A_1 = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & m_A == 0.7 & m_B == 0.7 & network == "FP_Galetti&Pizo_1996_bin"),
               aes(x = g, y = mean_final_mut_mat)) +
  geom_errorbar(aes(ymin = lower_lim_final_mut_mat, ymax = upper_lim_final_mut_mat), 
                width = 0.004, color = "gray50", size = 0.7) +
  geom_point(size = 4, shape = 19, color = "#1AB245") + 
  xlab("Mean gene flow") +
  ylab("Network-level trait matching") +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.position = "none")

# reading matrix
mat = as.matrix(read.table("data/empirical_networks/binary/FP_Galetti&Pizo_1996_bin.txt", 
                           sep = " ", header = FALSE))
n_row = nrow(mat)
n_col = ncol(mat)
n_sp = n_row + n_col
f = rbind(cbind(matrix(0, n_row, n_row), mat),
          cbind(t(mat), matrix(0, n_col, n_col)))

# pairwise trait matching for one simulation with low gene flow
traits_sim = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/all_networks_mA07_mB07_g0-03/FP_Galetti&Pizo_1996_bin/FP_Galetti&Pizo_1996_bin_mA07_mB07_g0_siteA_sim10.csv",
                      row.names = 1)
# trait matching
matching = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                          z = as.numeric(traits_sim["z_A_final", ]),
                          method = "exponential", alpha = 0.2)
pair_mat = matching[[2]]
# changing order of species
pair_mat$row = n_row - as.numeric(substr(as.character(pair_mat$row), start = 2, 
                                 stop = nchar(as.character(pair_mat$row)))) + 1
pair_mat$column = as.numeric(substr(as.character(pair_mat$column), start = 2,
                                    stop = nchar(as.character(pair_mat$column))))
# plotting
p_A_2 = ggplot(data = pair_mat, aes(x = column, y = row, fill = matching,
                                    height = 0.9, width = 0.9)) + 
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(option = "plasma", direction = -1, limits = c(0, 1)) +
  ylab("Plant species") +
  xlab("Animal species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# pairwise trait matching for one simulation with high gene flow
traits_sim = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/all_networks_mA07_mB07_g0-03/FP_Galetti&Pizo_1996_bin/FP_Galetti&Pizo_1996_bin_mA07_mB07_g03_siteA_sim5.csv",
                      row.names = 1)
# trait matching
matching = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                          z = as.numeric(traits_sim["z_A_final", ]),
                          method = "exponential", alpha = 0.2)
pair_mat = matching[[2]]
# changing order of species
pair_mat$row = n_row - as.numeric(substr(as.character(pair_mat$row), start = 2, 
                                         stop = nchar(as.character(pair_mat$row)))) + 1
pair_mat$column = as.numeric(substr(as.character(pair_mat$column), start = 2,
                                    stop = nchar(as.character(pair_mat$column))))
# plotting
p_A_3 = ggplot(data = pair_mat, aes(x = column, y = row, fill = matching,
                                    height = 0.9, width = 0.9)) + 
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(option = "plasma", direction = -1, limits = c(0, 1)) +
  ylab("Plant species") +
  xlab("Animal species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#------------------------------------------------------------------------------#
# Fig. 2B: trait matching vs gene flow for one network with mA = 0.9, mB = 0.1 #
#------------------------------------------------------------------------------#

# mean final trait matching at site A vs gene flow for m_A = m_B = 0.7
p_B_1 = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & m_A == 0.9 & m_B == 0.1 & network == "FP_Galetti&Pizo_1996_bin"),
               aes(x = g, y = mean_final_mut_mat)) +
  geom_errorbar(aes(ymin = lower_lim_final_mut_mat, ymax = upper_lim_final_mut_mat), 
                width = 0.004, color = "gray50", size = 0.7) +
  geom_point(size = 4, shape = 19, color = "#1AB245") + 
  xlab("Mean gene flow") +
  ylab("Network-level trait matching") +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.position = "none")

# pairwise trait matching for one simulation with low gene flow
traits_sim = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/all_networks_mA09_mB01_g0-03/FP_Galetti&Pizo_1996_bin/FP_Galetti&Pizo_1996_bin_mA09_mB01_g0_siteA_sim7.csv",
                      row.names = 1)
# trait matching
matching = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                          z = as.numeric(traits_sim["z_A_final", ]),
                          method = "exponential", alpha = 0.2)
pair_mat = matching[[2]]
# changing order of species
pair_mat$row = n_row - as.numeric(substr(as.character(pair_mat$row), start = 2, 
                                         stop = nchar(as.character(pair_mat$row)))) + 1
pair_mat$column = as.numeric(substr(as.character(pair_mat$column), start = 2,
                                    stop = nchar(as.character(pair_mat$column))))
# plotting
p_B_2 = ggplot(data = pair_mat, aes(x = column, y = row, fill = matching,
                                    height = 0.9, width = 0.9)) + 
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(option = "plasma", direction = -1, limits = c(0, 1)) +
  ylab("Plant species") +
  xlab("Animal species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# pairwise trait matching for one simulation with high gene flow
traits_sim = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/all_networks_mA09_mB01_g0-03/FP_Galetti&Pizo_1996_bin/FP_Galetti&Pizo_1996_bin_mA09_mB01_g03_siteA_sim9.csv",
                      row.names = 1)
# trait matching
matching = MatchingMutNet(n_sp = n_sp, n_row = n_row, n_col = n_col, f = f,
                          z = as.numeric(traits_sim["z_A_final", ]),
                          method = "exponential", alpha = 0.2)
pair_mat = matching[[2]]
# changing order of species
pair_mat$row = n_row - as.numeric(substr(as.character(pair_mat$row), start = 2, 
                                         stop = nchar(as.character(pair_mat$row)))) + 1
pair_mat$column = as.numeric(substr(as.character(pair_mat$column), start = 2,
                                    stop = nchar(as.character(pair_mat$column))))
# plotting
p_B_3 = ggplot(data = pair_mat, aes(x = column, y = row, fill = matching,
                                    height = 0.9, width = 0.9)) + 
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(option = "plasma", direction = -1, limits = c(0, 1)) +
  ylab("Plant species") +
  xlab("Animal species") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = pair_mat, aes(x = row, y = column, fill = matching,
                                       height = 0.9, width = 0.9)) + 
  geom_tile(color = "black", size = 0.2) +
  scale_fill_viridis(option = "plasma", direction = -1, 
                     limits = c(0, 1), name = "Trait\nmatching") +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.55, "cm"))
  
legend = get_legend(p_legend)

# generating the whole figure
fig2 = ggdraw() +
  draw_plot(p_A_1, 0, 0.43, 0.5, 0.5) +
  draw_plot(p_A_2, 0, 0.02, 0.21, 0.38) +
  draw_plot(p_A_3, 0.22, 0.02, 0.21, 0.38) +
  draw_plot(p_B_1, 0.5, 0.43, 0.5, 0.5) +
  draw_plot(p_B_2, 0.47, 0.02, 0.21, 0.38) +
  draw_plot(p_B_3, 0.69, 0.02, 0.21, 0.38) +
  draw_plot(legend, 0.88, 0.15, 0.15, 0.15) +
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 1), size = 24)

# saving
save_plot("fig2.pdf", fig2, ncol = 2, nrow = 2, base_aspect_ratio = 1.5)

#-----------------------------------------------------------------------------------------------------#