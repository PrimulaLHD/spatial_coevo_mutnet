#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Choosing networks for the sensitivity analysis
#
# Returns:
#   

# loading functions and packages
library(ggplot2)
library(cowplot)
library(viridis)

# network structure data
network_structure = read.csv("output/data/network_structure/network_structure.csv")
network_structure$PC1 = -network_structure$PC1
network_structure$PC2 = -network_structure$PC2

# simulation data
summ_final_mut_mat_df = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/summary_coevo_results_final_mut_matching.csv")

# changing mutualism names
summ_final_mut_mat_df$mutualism = as.character(summ_final_mut_mat_df$mutualism)
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AE"] = "ant-EFN"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AM"] = "ant-myrmecophyte"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "CC"] = "cleaner-client"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FA"] = "fish-anemone"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FP"] = "frugivore-plant"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "PP"] = "pollinator-plant"

# adding PC1 and PC2 values
summ_final_mut_mat_df$PC1 = rep(rep(network_structure$PC1, each = 2), times = 15*31) 
summ_final_mut_mat_df$PC2 = rep(rep(network_structure$PC2, each = 2), times = 15*31)

set.seed(482)
select_sp = sample(1:nrow(network_structure),
                   nrow(network_structure)/2,
                   replace = FALSE)

sort(network_structure_select$network) = network_structure[select_sp, ]
table(network_structure_select$mutualism)/table(network_structure$mutualism)

# defining a color palette and ordering according to factor mutualism
my_palette = c("#00EEEE", "#CD2626", "#4876FF", "#FF7F00", "#A020F0", "#1AB245")

ggplot(network_structure, aes(x = PC1, y = PC2, fill = mutualism)) +
  scale_fill_manual(values = my_palette) +
  geom_point(size = 3, shape = 21) + 
  xlab("") +
  ylab("") +
  scale_x_continuous(limits = c(-3, 6.2)) +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

#-----------------------------------------------------------------------------------------------------#