#-----------------------------------------------------------------------------------------------------#

# Description: 
#  Generates Fig. 4 of the paper, which shows how the gradual disruption of gene flow in mutualistic
#  networks would affect the evolution of trait matching using an analytical approximation of the 
#  coevolutionary model.
#
# Returns:
#  Saves the figure as a pdf file.    

# loading functions and packages
library(plyr)
library(ggplot2)
library(cowplot)

# defining folder with results
result_files = dir("output/data/simulations_empirical_networks/simulations_geneflow_removal")

# selecting the files with the desired gene flow values
result_files = result_files[grep("g03-0.csv", result_files)]
# result_files = result_files[grep("g03-001.csv", result_files)]
# result_files = result_files[grep("g03-005.csv", result_files)]
# result_files = result_files[grep("g03-01.csv", result_files)]
# result_files = result_files[grep("g01-0.csv", result_files)]
# result_files = result_files[grep("g01-001.csv", result_files)]
# result_files = result_files[grep("g01-005.csv", result_files)]

# getting network names
net_names = gsub( "_geneflow.*$", "", result_files)
net_names_bin = gsub("wgt", "bin", net_names)

# summarizing simulation results
summ_df = c()

for (i in 1:length(result_files)) {
  current_df = read.csv(paste("output/data/simulations_empirical_networks/simulations_geneflow_removal/", 
                              result_files[i], sep = ""))
  
  curr_summ_df = ddply(current_df, c("pert", "m_A", "m_B"),
                       summarize,
                       network = unique(network),
                       matching_A = mean(matching_A),
                       matching_B = mean(matching_B),
                       dim_T = unique(dim_T),
                       sum_eigen_T = mean(sum_eigen_T),
                       mean_eigen_T = mean(mean_eigen_T),
                       sd_eigen_T = mean(sd_eigen_T),
                       eigen1_T = mean(eigen1_T),
                       eigen2_T = mean(eigen2_T),
                       det_T = mean(det_T),
                       mean_row_corr_T = mean(mean_row_corr_T),
                       sd_g = mean(sd_g))
  
  summ_df = rbind(summ_df, curr_summ_df)
}

# magnitude of 1st eigenvalue in relation to all others
summ_df$prop_eigen1_T = summ_df$eigen1_T/summ_df$sum_eigen_T

# obtaining mutualism names
summ_df$mutualism = substr(summ_df$network, start = 1, stop = 2) 

# changing mutualism names
summ_df$mutualism[summ_df$mutualism == "AE"] = "ants - nectary-bearing plants"
summ_df$mutualism[summ_df$mutualism == "AM"] = "ants - myrmecophytes"
summ_df$mutualism[summ_df$mutualism == "FP"] = "seed dispersal"
summ_df$mutualism[summ_df$mutualism == "PP"] = "pollination"

# defining a color palette and ordering according to factor mutualism
my_palette = c("#CD2626", "#00EEEE", "#1AB245", "#A020F0")

# plotting
fig4 = ggplot(data = subset(summ_df, m_A == 0.5 & m_B == 0.5), 
           aes(x = pert, y = matching_A, group = network, color = mutualism)) +
  geom_point(size = 2) +
  geom_line(size = 0.4) +
  scale_color_manual(values = my_palette) +
  xlab("Fraction of species without gene flow") +
  ylab("Mean trait matching") +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_blank())

# saving
save_plot("output/figs/fig4.pdf", fig4, ncol = 1, nrow = 1, base_aspect_ratio = 2.2)

#-----------------------------------------------------------------------------------------------------#