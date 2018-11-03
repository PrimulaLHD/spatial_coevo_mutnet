#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads networks structure data and summary from coevolution simulation and generates different plots.
#
# Returns:
#   Plots.

# loading functions and packages
library(ggplot2)
library(plyr)
library(cowplot)
library(viridis)

# network structure data
network_structure = read.csv("output/data/network_structure/network_structure.csv")

# summary of coevolution results
summaries = dir("output/data/simulations_empirical_networks/sensitivity_analysis/g_correlated_with_degree/summary_coevo_results/")

# extracting mA, mB and g values from files names
split_mA = sub(".*_mA", "", summaries)
mA_char = substr(split_mA, start = 1, stop = 2)
mA = as.numeric(sub("(.*0)(.*)", "\\1.\\2", mA_char))
split_mB = sub(".*_mB", "", summaries)
mB_char = substr(split_mB, start = 1, stop = 2)
mB = as.numeric(sub("(.*0)(.*)", "\\1.\\2", mB_char))

# calculating mean trait matching for many simulations
summ_final_mut_mat_df = c()
for (i in 1:length(summaries)) {
  current_df = read.csv(paste("output/data/simulations_empirical_networks/sensitivity_analysis/g_correlated_with_degree/summary_coevo_results/", 
                              summaries[i], sep = ""))
  current_df$m_A = rep(mA[i], nrow(current_df))
  current_df$m_B = rep(mB[i], nrow(current_df))
  curr_summ_final_mut_mat_df = ddply(current_df, c("g", "network", "site"),
                                     summarize, mean_final_mut_mat = mean(mean_final_mut_matching),
                                     lower_lim_final_mut_mat = quantile(mean_final_mut_matching, probs = 0.025),
                                     upper_lim_final_mut_mat = quantile(mean_final_mut_matching, probs = 0.975),
                                     mutualism = unique(mutualism),
                                     m_A = unique(m_A),
                                     m_B = unique(m_B))
  summ_final_mut_mat_df = rbind(summ_final_mut_mat_df, curr_summ_final_mut_mat_df)
}

# saving/reading results
write.csv(summ_final_mut_mat_df, row.names = FALSE, 
          file = "output/data/simulations_empirical_networks/sensitivity_analysis/g_correlated_with_degree/summary_coevo_results/summary_coevo_results_final_mut_matching.csv")
summ_final_mut_mat_df = read.csv("output/data/simulations_empirical_networks/sensitivity_analysis/g_correlated_with_degree/summary_coevo_results/summary_coevo_results_final_mut_matching.csv")

# changing mutualism names
summ_final_mut_mat_df$mutualism = as.character(summ_final_mut_mat_df$mutualism)
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AE"] = "ants-nectary-bearing plants"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AM"] = "ants-myrmecophytes"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "CC"] = "marine cleaning"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FA"] = "anemones-anemonefishes"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FP"] = "seed dispersal"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "PP"] = "pollination"

# defining a color palette and ordering according to factor mutualism
orange = "#FF7F00"
red = "#CD2626"
cyan = "#00EEEE"
blue = "#4876FF"
green = "#1AB245"
purple = "#A020F0"
my_palette = c(orange, red, cyan, blue, green, purple)

# Fig A: mean final trait matching as a function of gene flow for site A
p_A = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & 
                             ((m_A == 0.7 & m_B == 0.7) | (m_A == 0.9 & m_B == 0.1))),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.3, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  xlab("Mean gene flow") +
  ylab("Mean trait matching at site A") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig B: mean final trait matching as a function of gene flow for site B
p_B = ggplot(data = subset(summ_final_mut_mat_df, site == "B" & 
                             ((m_A == 0.7 & m_B == 0.7) | (m_A == 0.9 & m_B == 0.1))),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.3, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  xlab("Mean gene flow") +
  ylab("Mean trait matching at site B") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = subset(summ_final_mut_mat_df, site == "A"),
                  aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.3, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  guides(color = guide_legend(override.aes = list(size = 2.2))) +
  theme(legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

legend = get_legend(p_legend)

# generating the whole figure
figSI = ggdraw() +
  draw_plot(p_A, 0, 0.5, 0.72, 0.5) +
  draw_plot(p_B, 0, 0, 0.72, 0.5) +
  draw_plot(legend, 0.85, 0.52, 0.01, 0.01) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.51), size = 23)

# saving
save_plot("output/figs/figS6.pdf", figSI, ncol = 1, nrow = 2, base_aspect_ratio = 2.2)

#-----------------------------------------------------------------------------------------------------#