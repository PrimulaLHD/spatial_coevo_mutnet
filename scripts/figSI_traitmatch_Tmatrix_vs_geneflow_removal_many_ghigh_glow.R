#-----------------------------------------------------------------------------------------------------#

# Description: 
#
#
# Returns:
#   

# loading functions and packages
library(ggplot2)
library(cowplot)
library(plyr)
library(viridis)

# defining folder with results
result_files = dir("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_geneflow_removal/")

# defining gene flow scenarios
g_scen = data.frame(g_high = c(0.3, 0.3, 0.3, 0.3, 0.1, 0.1, 0.1), 
                    g_low = c(0, 0.01, 0.05, 0.1, 0, 0.01, 0.05))
g_scen$g_high_char = sub(".", "", as.character(g_scen$g_high), fixed = TRUE)
g_scen$g_low_char = sub(".", "", as.character(g_scen$g_low), fixed = TRUE)

summ_df_all_g_scen = c()
for (j in 1:nrow(g_scen)) {
  # selecting the files with the desired gene flow values
  file_name = paste("g", g_scen$g_high_char[j], "-", g_scen$g_low_char[j], ".csv", sep = "")
  result_files_curr = result_files[grep(file_name, result_files)]
  # getting network names
  net_names = gsub( "_geneflow.*$", "", result_files_curr)
  net_names_bin = gsub("wgt", "bin", net_names)
  
  # summarizing simulation results
  summ_df = c()
  for (i in 1:length(result_files_curr)) {
    current_df = read.csv(paste("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/simulations_geneflow_removal/", 
                                result_files_curr[i], sep = ""))
    curr_summ_df = ddply(current_df, c("pert", "m_A", "m_B"),
                         summarize,
                         network = unique(network),
                         matching_A = mean(matching_A),
                         matching_B = mean(matching_B),
                         sd_g = mean(sd_g))
    summ_df = rbind(summ_df, curr_summ_df)
  }
  
  # obtaining mutualism names
  summ_df$mutualism = substr(summ_df$network, start = 1, stop = 2) 
  
  # changing mutualism names
  summ_df$mutualism[summ_df$mutualism == "AE"] = "ants - nectary-bearing plants"
  summ_df$mutualism[summ_df$mutualism == "AM"] = "ants - myrmecophytes"
  summ_df$mutualism[summ_df$mutualism == "FP"] = "seed dispersal"
  summ_df$mutualism[summ_df$mutualism == "PP"] = "pollination"
  
  # adding gene flow scenario
  summ_df$g_high = rep(g_scen$g_high[j], nrow(summ_df))
  summ_df$g_low = rep(g_scen$g_low[j], nrow(summ_df))
  
  summ_df_all_g_scen = rbind(summ_df_all_g_scen, summ_df)
}

# plotting
p = ggplot(data = subset(summ_df_all_g_scen, m_A == 0.5 & m_B == 0.5),
           aes(x = pert, y = matching_A, group = network, color = sd_g)) +
  geom_point(size = 2.7) +
  geom_line(size = 0.7) +
  scale_color_viridis(option = "plasma", 
                      name = "Standard deviation\nof gene flow\nacross species") +
  xlab("Fraction of species without gene flow") +
  ylab("Mean trait matching at site A") +
  facet_wrap(g_high ~ g_low, labeller = label_both) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 22),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.key.size = unit(0.8, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 17))

# saving
save_plot("figSI.pdf", p, ncol = 3, nrow = 3, base_aspect_ratio = 1.2)

#-----------------------------------------------------------------------------------------------------#