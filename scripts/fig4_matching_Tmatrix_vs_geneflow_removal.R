#-----------------------------------------------------------------------------------------------------#

# Description: 
#   (Organize this code). Code for figure 4: gene flow perturbation, T matrix and trait matching
#
# Returns:
#   

library(plyr)
library(ggplot2)
library(cowplot)

net_struct = read.csv("output/data/network_structure/network_structure.csv")
result_files = dir("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_geneflow_removal/")

result_files = result_files[grep("g03-0.csv", result_files)]
# result_files = result_files[grep("g03-001.csv", result_files)]
# result_files = result_files[grep("g03-005.csv", result_files)]
# result_files = result_files[grep("g03-01.csv", result_files)]
# result_files = result_files[grep("g01-0.csv", result_files)]
# result_files = result_files[grep("g01-001.csv", result_files)]
# result_files = result_files[grep("g01-005.csv", result_files)]

net_names = gsub( "_geneflow.*$", "", result_files)
net_names_bin = gsub("wgt", "bin", net_names)

summ_df = c()

for (i in 1:length(result_files)) {
  current_df = read.csv(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_geneflow_removal/", 
                              result_files[i], sep = ""))
  
  curr_summ_df = ddply(current_df, c("pert", "m_A", "m_B"),
                       summarize,
                       network = unique(network),
                       matching_A = mean(matching_A),
                       matching_B = mean(matching_B),
                       dim_T = unique(dim_T),
                       eigen1_L = mean(eigen1_L),
                       eigen2_L = mean(eigen2_L),
                       sum_eigen_T = mean(sum_eigen_T),
                       mean_eigen_T = mean(mean_eigen_T),
                       sd_eigen_T = mean(sd_eigen_T),
                       eigen1_T = mean(eigen1_T),
                       eigen2_T = mean(eigen2_T),
                       det_T = mean(det_T),
                       mean_col_sd_T = mean(mean_col_sd_T))
  
  summ_df = rbind(summ_df, curr_summ_df)
}

summ_df$prop_eigen1_T = summ_df$eigen1_T/summ_df$sum_eigen_T

summ_df$mutualism = substr(summ_df$network, start = 1, stop = 2) 

# changing mutualism names
summ_df$mutualism[summ_df$mutualism == "AE"] = "ant-nectary bearing plant"
summ_df$mutualism[summ_df$mutualism == "AM"] = "ant-myrmecophyte"
summ_df$mutualism[summ_df$mutualism == "FP"] = "seed dispersal"
summ_df$mutualism[summ_df$mutualism == "PP"] = "pollination"

# defining a color palette and ordering according to factor mutualism
my_palette = c("#CD2626", "#00EEEE", "#A020F0", "#1AB245")

ggplot(data = subset(summ_df, m_A == 0.5 & m_B == 0.5), 
       aes(x = pert, y = matching_A, group = network, color = mutualism)) +
  geom_point(size = 2.4) +
  geom_line(size = 0.5) +
  scale_color_manual(values = my_palette) +
  xlab("Fraction of species without gene flow") +
  ylab("Mean trait matching") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

ggplot(data = summ_df,
       aes(x = mean_col_sd_T, y = matching_A, color = network)) +
  geom_line(size = 0.8, alpha = 0.8) +
  #scale_color_manual(values = my_palette) +
  xlab("Sd rows of T-matrix") +
  ylab("Mean trait matching") +
  facet_grid(m_A ~ m_B, labeller = label_both, scales = "free") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none")

#-----------------------------------------------------------------------------------------------------#