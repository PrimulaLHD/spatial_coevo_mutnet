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

# defining a color palette and ordering according to factor mutualism
my_palette = c("#CD2626", "#00EEEE", "#1AB245", "#A020F0")

# plotting
p_A = ggplot(data = subset(summ_df_all_g_scen, g_high == 0.3 & g_low == 0),
             aes(x = pert, y = matching_A, group = network, color = mutualism)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_line(size = 0.2, alpha = 0.5) +
  scale_color_manual(values = my_palette) +
  xlab("Fraction of species without gene flow") +
  ylab("Mean trait matching at site A") +
  facet_grid(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "none")

p_B = ggplot(data = subset(summ_df_all_g_scen, g_high == 0.3 & g_low == 0),
             aes(x = pert, y = matching_B, group = network, color = mutualism)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_line(size = 0.2, alpha = 0.5) +
  scale_color_manual(values = my_palette) +
  xlab("Fraction of species without gene flow") +
  ylab("Mean trait matching at site B") +
  facet_grid(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 13),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = subset(summ_df_all_g_scen, g_high == 0.3 & g_low == 0),
                  aes(x = pert, y = matching_B, group = network, color = mutualism)) +
  geom_point(size = 0.8) + 
  scale_color_manual(values = my_palette) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_blank())

legend = get_legend(p_legend)

# generating the whole figure
figSI = ggdraw() +
  draw_plot(p_A, 0, 0.5, 0.72, 0.5) +
  draw_plot(p_B, 0, 0, 0.72, 0.5) +
  draw_plot(legend, 0.85, 0.52, 0.01, 0.01) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.51), size = 22)

# saving
save_plot("figSI.pdf", figSI, ncol = 1, nrow = 2, base_aspect_ratio = 1.8)

#-----------------------------------------------------------------------------------------------------#