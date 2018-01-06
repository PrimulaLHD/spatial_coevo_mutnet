#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads data on pairwise trait matching, summarizes data and plot pairwise trait matching as a function
#   of species degrees.
#
# Returns:
#   A plot.

# loading packages
library(plyr)
library(ggplot2)
library(cowplot)
library(viridis)

#----------------------------------#
# (A) Two hotspots (mA = mB = 0.7) #
#----------------------------------#

# reading data
final_mut_mat_df_all_nets_g0 = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA07_mB07_g0_siteA.csv")
final_mut_mat_df_all_nets_g03 = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA07_mB07_g03_siteA.csv")

# computing difference in trait matching between high and low gene flow
final_mut_mat_df_all_nets = final_mut_mat_df_all_nets_g0
final_mut_mat_df_all_nets$dif_matching = final_mut_mat_df_all_nets_g03$matching - final_mut_mat_df_all_nets_g0$matching

# reading network structure data
net_struct = read.csv("output/data/network_structure/network_structure.csv")
# adding mutualism type
n_sim = 100
final_mut_mat_df_all_nets$mutualism = rep(net_struct$mutualism, n_sim*net_struct$interactions)

# changing mutualism names
final_mut_mat_df_all_nets$mutualism = as.character(final_mut_mat_df_all_nets$mutualism)
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "AE"] = "ants-nectary-bearing plants"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "AM"] = "ants-myrmecophytes"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "CC"] = "marine cleaning"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "FA"] = "anemones-anemonefishes"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "FP"] = "seed dispersal"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "PP"] = "pollination"

# computing mean trait matching difference for each species pair across simulations
summ_final_dif_matching = ddply(final_mut_mat_df_all_nets, c("network", "interaction"),
                                summarize, mean_dif_matching = mean(dif_matching),
                                mutualism = unique(mutualism),
                                col = unique(col),
                                row = unique(row),
                                k_col = unique(k_col),
                                k_row = unique(k_row))

# plotting 
p_A = ggplot(data = subset(summ_final_dif_matching, k_col < 40 & k_row < 40), 
             aes(x = k_col, y = k_row, color = mean_dif_matching)) +
  geom_point(size = 1.8, alpha = 0.7) +
  facet_wrap( ~ mutualism) +
  xlab("First species standardized degree") +
  ylab("Second species\nstandardized degree") +
  scale_color_viridis(name = "Mean difference in\npairwise trait matching") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.title = element_text(size = 13))

#-------------------------------------------------------#
# (B) One hotspot and one coldspot (mA = 0.9, mB = 0.1) #
#-------------------------------------------------------#

# reading data
final_mut_mat_df_all_nets_g0 = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA09_mB01_g0_siteA.csv")
final_mut_mat_df_all_nets_g03 = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA09_mB01_g03_siteA.csv")

# computing difference in trait matching between high and low gene flow
final_mut_mat_df_all_nets = final_mut_mat_df_all_nets_g0
final_mut_mat_df_all_nets$dif_matching = final_mut_mat_df_all_nets_g03$matching - final_mut_mat_df_all_nets_g0$matching

# reading network structure data
net_struct = read.csv("output/data/network_structure/network_structure.csv")
# adding mutualism type
n_sim = 100
final_mut_mat_df_all_nets$mutualism = rep(net_struct$mutualism, n_sim*net_struct$interactions)

# changing mutualism names
final_mut_mat_df_all_nets$mutualism = as.character(final_mut_mat_df_all_nets$mutualism)
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "AE"] = "ants-nectary-bearing plants"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "AM"] = "ants-myrmecophytes"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "CC"] = "marine cleaning"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "FA"] = "anemones-anemonefishes"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "FP"] = "seed dispersal"
final_mut_mat_df_all_nets$mutualism[final_mut_mat_df_all_nets$mutualism == "PP"] = "pollination"

# computing mean trait matching difference for each species pair across simulations
summ_final_dif_matching = ddply(final_mut_mat_df_all_nets, c("network", "interaction"),
                                summarize, mean_dif_matching = mean(dif_matching),
                                mutualism = unique(mutualism),
                                col = unique(col),
                                row = unique(row),
                                k_col = unique(k_col),
                                k_row = unique(k_row))

# plotting 
p_B = ggplot(data = subset(summ_final_dif_matching, k_col < 40 & k_row < 40), 
             aes(x = k_col, y = k_row, color = mean_dif_matching)) +
  geom_point(size = 1.8, alpha = 0.7) +
  facet_wrap( ~ mutualism) +
  xlab("First species standardized degree") +
  ylab("Second species\nstandardized degree") +
  scale_color_viridis(name = "Mean difference in\npairwise trait matching",
                      direction = -1) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        legend.title = element_text(size = 13))

#---------------#
# Complete plot #
#---------------#

figSI = plot_grid(p_A, p_B, labels = c("A", "B"), label_size = 23, 
                  nrow = 2, align = "v")

save_plot("figSI.pdf", figSI,
          ncol = 1, nrow = 2,
          base_aspect_ratio = 2.35)

#-----------------------------------------------------------------------------------------------------#