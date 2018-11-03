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
final_mut_mat_df_all_nets_g0 = read.csv("output/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA07_mB07_g0_siteA.csv")
final_mut_mat_df_all_nets_g03 = read.csv("output/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA07_mB07_g03_siteA.csv")

# reading network structure data
net_struct = read.csv("output/data/network_structure/network_structure.csv")
# changing mutualism names
net_struct$mutualism = as.character(net_struct$mutualism)
net_struct$mutualism[net_struct$mutualism == "AE"] = "ants-nectary-bearing plants"
net_struct$mutualism[net_struct$mutualism == "AM"] = "ants-myrmecophytes"
net_struct$mutualism[net_struct$mutualism == "CC"] = "marine cleaning"
net_struct$mutualism[net_struct$mutualism == "FA"] = "anemones-anemonefishes"
net_struct$mutualism[net_struct$mutualism == "FP"] = "seed dispersal"
net_struct$mutualism[net_struct$mutualism == "PP"] = "pollination"

# computing mean trait matching difference for each species pair across simulations
summ_final_mut_matching_g0 = ddply(final_mut_mat_df_all_nets_g0, c("network", "interaction"),
                                   summarize, mean_matching = mean(matching),
                                   col = unique(col),
                                   row = unique(row),
                                   k_col = unique(k_col),
                                   k_row = unique(k_row))
summ_final_mut_matching_g03 = ddply(final_mut_mat_df_all_nets_g03, c("network", "interaction"),
                                    summarize, mean_matching = mean(matching),
                                    col = unique(col),
                                    row = unique(row),
                                    k_col = unique(k_col),
                                    k_row = unique(k_row))

# creating new data frame
summ_final_dif_matching = summ_final_mut_matching_g0[ , c("network", "interaction", "col", "row", "k_col", "k_row")]
# computing difference in trait matching with and without gene flow
summ_final_dif_matching$mean_dif_matching = summ_final_mut_matching_g03$mean_matching - summ_final_mut_matching_g0$mean_matching
# adding mutualism type
summ_final_dif_matching$mutualism = rep(net_struct$mutualism, net_struct$interactions)

# plotting 
p_A = ggplot(data = subset(summ_final_dif_matching, k_col < 40 & k_row < 40), 
             aes(x = k_col, y = k_row, color = mean_dif_matching)) +
  geom_point(size = 1.8, alpha = 0.7) +
  facet_wrap( ~ mutualism) +
  xlab("First species degree") +
  ylab("Second species degree") +
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
final_mut_mat_df_all_nets_g0 = read.csv("output/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA09_mB01_g0_siteA.csv")
final_mut_mat_df_all_nets_g03 = read.csv("output/data/simulations_empirical_networks/summary_coevo_results/coevo_results_pairwise_matching_sp_degree_mA09_mB01_g03_siteA.csv")

# computing mean trait matching difference for each species pair across simulations
summ_final_mut_matching_g0 = ddply(final_mut_mat_df_all_nets_g0, c("network", "interaction"),
                                   summarize, mean_matching = mean(matching),
                                   col = unique(col),
                                   row = unique(row),
                                   k_col = unique(k_col),
                                   k_row = unique(k_row))
summ_final_mut_matching_g03 = ddply(final_mut_mat_df_all_nets_g03, c("network", "interaction"),
                                    summarize, mean_matching = mean(matching),
                                    col = unique(col),
                                    row = unique(row),
                                    k_col = unique(k_col),
                                    k_row = unique(k_row))

# creating new data frame
summ_final_dif_matching = summ_final_mut_matching_g0[ , c("network", "interaction", "col", "row", "k_col", "k_row")]
# computing difference in trait matching with and without gene flow
summ_final_dif_matching$mean_dif_matching = summ_final_mut_matching_g03$mean_matching - summ_final_mut_matching_g0$mean_matching
# adding mutualism type
summ_final_dif_matching$mutualism = rep(net_struct$mutualism, net_struct$interactions)

# plotting 
p_B = ggplot(data = subset(summ_final_dif_matching, k_col < 40 & k_row < 40), 
             aes(x = k_col, y = k_row, color = mean_dif_matching)) +
  geom_point(size = 1.8, alpha = 0.7) +
  facet_wrap( ~ mutualism) +
  xlab("First species degree") +
  ylab("Second species degree") +
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

save_plot("output/figs/figS3.pdf", figSI,
          ncol = 1, nrow = 2,
          base_aspect_ratio = 2.35)

#-----------------------------------------------------------------------------------------------------#