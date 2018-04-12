#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads summary from coevolution simulations and generates plots of environmental matching and 
#   geographical divergence as a function of gene flow.
#
# Returns:
#   Saves plots.

# loading functions and packages
library(ggplot2)
library(cowplot)

# reading data
summ_final_env_mat_df = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/summary_coevo_results_final_env_matching.csv")
summ_final_geo_div_df = read.csv("~/LUCAS/spatial_coevo_mutnet_results/data/simulations_empirical_networks/summary_coevo_results/summary_coevo_results_final_geo_div.csv")

# changing mutualism names
summ_final_env_mat_df$mutualism = as.character(summ_final_env_mat_df$mutualism)
summ_final_env_mat_df$mutualism[summ_final_env_mat_df$mutualism == "AE"] = "ants-nectary-bearing plants"
summ_final_env_mat_df$mutualism[summ_final_env_mat_df$mutualism == "AM"] = "ants-myrmecophytes"
summ_final_env_mat_df$mutualism[summ_final_env_mat_df$mutualism == "CC"] = "marine cleaning"
summ_final_env_mat_df$mutualism[summ_final_env_mat_df$mutualism == "FA"] = "anemones-anemonefishes"
summ_final_env_mat_df$mutualism[summ_final_env_mat_df$mutualism == "FP"] = "seed dispersal"
summ_final_env_mat_df$mutualism[summ_final_env_mat_df$mutualism == "PP"] = "pollination"
summ_final_geo_div_df$mutualism = as.character(summ_final_geo_div_df$mutualism)
summ_final_geo_div_df$mutualism[summ_final_geo_div_df$mutualism == "AE"] = "ants-nectary-bearing plants"
summ_final_geo_div_df$mutualism[summ_final_geo_div_df$mutualism == "AM"] = "ants-myrmecophytes"
summ_final_geo_div_df$mutualism[summ_final_geo_div_df$mutualism == "CC"] = "marine cleaning"
summ_final_geo_div_df$mutualism[summ_final_geo_div_df$mutualism == "FA"] = "anemones-anemonefishes"
summ_final_geo_div_df$mutualism[summ_final_geo_div_df$mutualism == "FP"] = "seed dispersal"
summ_final_geo_div_df$mutualism[summ_final_geo_div_df$mutualism == "PP"] = "pollination"

# defining a color palette and ordering according to factor mutualism
orange = "#FF7F00"
red = "#CD2626"
cyan = "#00EEEE"
blue = "#4876FF"
green = "#1AB245"
purple = "#A020F0"
my_palette = c(orange, red, cyan, blue, green, purple)

# Fig A: mean final environmental matching as a function of gene flow for site A
p_A = ggplot(data = subset(summ_final_env_mat_df, site == "A" & 
                             ((m_A == 0.7 & m_B == 0.7) | (m_A == 0.9 & m_B == 0.1))),
             aes(x = g, y = mean_final_env_mat, color = mutualism)) +
  geom_point(size = 1.3, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  xlab("Mean gene flow") +
  ylab("Mean environmental\nmatching at site A") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig B: mean final geographical divergence as a function of gene flow
p_B = ggplot(data = subset(summ_final_geo_div_df, 
                           (m_A == 0.7 & m_B == 0.7) | (m_A == 0.9 & m_B == 0.1)),
             aes(x = g, y = mean_final_geo_div, color = mutualism)) +
  geom_point(size = 1.3, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  xlab("Mean gene flow") +
  ylab("Mean geographical divergence") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = subset(summ_final_geo_div_df, 
                                (m_A == 0.7 & m_B == 0.7) | (m_A == 0.9 & m_B == 0.1)),
                  aes(x = g, y = mean_final_geo_div, color = mutualism)) +
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
save_plot("figSI.pdf", figSI, ncol = 1, nrow = 2, base_aspect_ratio = 2.2)

#-----------------------------------------------------------------------------------------------------#