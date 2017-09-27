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
library(GGally)
library(network)
library(sna)
library(scales)

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

#------------------------------------------------#
# Fig. 3A: PC1 and PC2 of all empirical networks #
#------------------------------------------------#

# defining a color palette and ordering according to factor mutualism
my_palette = c("#00EEEE", "#CD2626", "#4876FF", "#FF7F00", "#A020F0", "#1AB245")

p_A = ggplot(network_structure, aes(x = PC1, y = PC2, color = mutualism)) +
  scale_colour_manual(values = my_palette) +
  geom_point(size = 3, shape = 19) + 
  xlab("") +
  ylab("") +
  scale_x_continuous(limits = c(-3, 6.2)) +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

#----------------------------------------------------------------------------------------------#
# Fig. 3B: predicted model values of trait matching vs PC1 and PC2 for mA = mB = 0.7 and g = 0 #
#----------------------------------------------------------------------------------------------#

mod = lm(mean_final_mut_mat ~ PC1 + PC2, 
         data = subset(summ_final_mut_mat_df, m_A == 0.7 & m_B == 0.7 & g == 0))
summary(mod)
pred_values = seq(-3, 6.2, by = 0.1)
pred = predict(mod, data.frame(PC1 = rep(pred_values, times = length(pred_values)),
                               PC2 = rep(pred_values, each = length(pred_values))))

pred_df = data.frame(PC1 = rep(pred_values, times = length(pred_values)),
                     PC2 = rep(pred_values, each = length(pred_values)),
                     matching = pred)

p_B = ggplot(data = pred_df, aes(x = PC1, y = PC2)) +
  geom_tile(aes(fill = matching)) +
  scale_fill_viridis(limits = c(0.48, 0.91)) +
  geom_point(data = network_structure, aes(x = PC1, y = PC2), size = 3) +
  scale_x_continuous(limits = c(-3, 6.2)) +
  scale_y_continuous(limits = c(-3, 3)) +
  xlab("") +
  ylab("Principal component 2 (PC2)") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

#------------------------------------------------------------------------------------------------#
# Fig. 3C: predicted model values of trait matching vs PC1 and PC2 for mA = mB = 0.7 and g = 0.3 #
#------------------------------------------------------------------------------------------------#

mod = lm(mean_final_mut_mat ~ PC1 + PC2, 
         data = subset(summ_final_mut_mat_df, m_A == 0.7 & m_B == 0.7 & g == 0.3))
summary(mod)
pred_values = seq(-3, 6.2, by = 0.1)
pred = predict(mod, data.frame(PC1 = rep(pred_values, times = length(pred_values)),
                               PC2 = rep(pred_values, each = length(pred_values))))

pred_df = data.frame(PC1 = rep(pred_values, times = length(pred_values)),
                     PC2 = rep(pred_values, each = length(pred_values)),
                     matching = pred)

p_C = ggplot(data = pred_df, aes(x = PC1, y = PC2)) +
  geom_tile(aes(fill = matching)) +
  scale_fill_viridis(limits = c(0.48, 0.91)) +
  geom_point(data = network_structure, aes(x = PC1, y = PC2), size = 3) +
  scale_x_continuous(limits = c(-3, 6.2)) +
  scale_y_continuous(limits = c(-3, 3)) +
  xlab("Principal component 1 (PC1)") +
  ylab("") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = pred_df, aes(x = PC1, y = PC2)) +
  geom_tile(aes(fill = matching)) +
  scale_fill_viridis(limits = c(0.48, 0.91), direction = -1, 
                     name = "Trait\nmatching") +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.55, "cm"))

legend = get_legend(p_legend)

#---------------------------------#
# Plotting 2 networks as examples #
#---------------------------------#

# plotting Izzo's ant-myrmecophyte network
mat = as.matrix(read.table("data/empirical_networks/binary/AM_Izzo_mtcCOLpen_bin.txt"))
n_row = nrow(mat)
n_col = ncol(mat)
# creating network object for plotting
net = network(mat, matrix.type = "bipartite", directed = FALSE, ignore.eval = FALSE)
# plotting
set.seed(42)
p_net_izzo = ggnet2(net, node.size = 2.8, color = "black", edge.color = "gray70",
                    edge.size = 1, mode = "kamadakawai", 
                    shape = c(rep(15, n_row), rep(16, n_col)))

# plotting Galetti & Pizo's frugivore-plant network
mat = as.matrix(read.table("data/empirical_networks/binary/FP_Galetti&Pizo_1996_bin.txt"))
n_row = nrow(mat)
n_col = ncol(mat)
# creating network object for plotting
net = network(mat, matrix.type = "bipartite", directed = FALSE, ignore.eval = FALSE)
# plotting
set.seed(4)
p_net_galetti = ggnet2(net, node.size = 2.4, color = "black", edge.color = "gray70",
                       edge.size = 1, mode = "kamadakawai", 
                       shape = c(rep(15, n_row), rep(16, n_col)))

# generating the whole figure
fig3 = ggdraw() +
  draw_plot(p_A, 0.17, 0.64, 0.6, 0.35) +
  draw_plot(p_B, 0.17, 0.32, 0.6, 0.35) +
  draw_plot(p_C, 0.17, 0, 0.6, 0.35) +
  draw_plot(p_net_izzo, 0, 0.73, 0.2, 0.2) +
  draw_plot(p_net_galetti, 0.76, 0.71, 0.24, 0.24) +
  draw_plot(legend, 0.82, 0.3, 0.1, 0.1) +
  draw_plot_label(c("A", "B", "C"), c(0.1, 0.1, 0.1), c(1, 0.7, 0.38), size = 24)

# saving
save_plot("fig3.pdf", fig3, ncol = 1, nrow = 3, base_aspect_ratio = 2.2)

#-----------------------------------------------------------------------------------------------------#