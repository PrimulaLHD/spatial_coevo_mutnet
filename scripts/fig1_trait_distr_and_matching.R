#-----------------------------------------------------------------------------------------------------#

# Description: 
#   
#
# Returns:
#   
install.packages(c("ggplot2", "cowplot", "viridis", "GGally", "network", "sna", "scales", "ggnetwork"))

library(ggplot2)
library(cowplot)
library(viridis)
library(GGally)
library(network)
library(sna)
library(scales)

source("functions/MatchingMutNet.R")

# creating matrix with 3 species
n_row = 2
n_col = 1
n_sp = n_row + n_col
mat = matrix(c(1, 1), nrow = n_row, ncol = n_col)
f = rbind(cbind(matrix(0, n_row, n_row), mat),
          cbind(t(mat), matrix(0, n_col, n_col)))

theta = c(1.5, 9, 5)

# plotting graph
net = network(f, directed = FALSE, ignore.eval = FALSE)

set.seed(18)
p_net = ggnet2(net, node.size = 4, color = c("gray80", "black", "gray50"), shape = c(15, 15, 16),
               edge.color = c("black", "black"), edge.size = c(0.8, 0.8), mode = "spring")

#-----------------------------------------------------#
# Fig. 1A: intermediate trait matching - no gene flow #
#-----------------------------------------------------#

z = c(3.6, 6.1, 4.7)

# trait distributions
p_A_1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = z[1], sd = 0.7),
                color = "gray80", size = 1.2) +
  stat_function(fun = dnorm, args = list(mean = z[2], sd = 0.7),
                color = "black", size = 1.2) +
  stat_function(fun = dnorm, args = list(mean = z[3], sd = 0.7),
                color = "gray50", size = 1.2) +
  geom_vline(xintercept = theta[1], color = "gray80", size = 0.9,
             linetype = "dashed") +
  geom_vline(xintercept = theta[2], color = "black", size = 0.9,
             linetype = "dashed") +
  geom_vline(xintercept = theta[3], color = "gray50", size = 0.9,
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 10)) +
  xlab("Trait value") +
  ylab("Density") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.position = "none")

# trait matching heatmap
matching = MatchingMutNet(n_sp, n_row, n_col, f, z = z,
                          method = "exponential", alpha = 0.2)

p_A_2 = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", limits = c(0, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
#----------------------------------------------#
# Fig. 1B: low trait matching - with gene flow #
#----------------------------------------------#

z = c(4.2, 5.9, 8.2)

# trait distributions
p_B_1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = z[1], sd = 0.7),
                color = "gray80", size = 1.2) +
  stat_function(fun = dnorm, args = list(mean = z[2], sd = 0.7),
                color = "black", size = 1.2) +
  stat_function(fun = dnorm, args = list(mean = z[3], sd = 0.7),
                color = "gray50", size = 1.2) +
  geom_vline(xintercept = theta[1], color = "gray80", size = 0.9,
             linetype = "dashed") +
  geom_vline(xintercept = theta[2], color = "black", size = 0.9,
             linetype = "dashed") +
  geom_vline(xintercept = theta[3], color = "gray50", size = 0.9,
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 10)) +
  xlab("Trait value") +
  ylab("Density") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.position = "none")

# trait matching heatmap
matching = MatchingMutNet(n_sp, n_row, n_col, f, z = z,
                          method = "exponential", alpha = 0.2)

p_B_2 = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                         height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", limits = c(0, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#-----------------------------------------------#
# Fig. 1C: high trait matching - with gene flow #
#-----------------------------------------------#

z = c(5.3, 5.7, 6)

# trait distributions
p_C_1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = z[1], sd = 0.7),
                color = "gray80", size = 1.2) +
  stat_function(fun = dnorm, args = list(mean = z[2], sd = 0.7),
                color = "black", size = 1.2) +
  stat_function(fun = dnorm, args = list(mean = z[3], sd = 0.7),
                color = "gray50", size = 1.2) +
  geom_vline(xintercept = theta[1], color = "gray80", size = 0.9,
             linetype = "dashed") +
  geom_vline(xintercept = theta[2], color = "black", size = 0.9,
             linetype = "dashed") +
  geom_vline(xintercept = theta[3], color = "gray50", size = 0.9,
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 10)) +
  xlab("Trait value") +
  ylab("Density") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.position = "none")

# trait matching heatmap
matching = MatchingMutNet(n_sp, n_row, n_col, f, z = z,
                          method = "exponential", alpha = 0.2)

p_C_2 = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                         height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", limits = c(0, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                         height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", limits = c(0, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.6) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.55, "cm"))

legend = get_legend(p_legend)

# generating the whole figure
fig1 = ggdraw() +
  draw_plot(p_A_1, 0, 0.66, 0.8, 0.22) +
  draw_plot(p_A_2, 0.58, 0.87, 0.1, 0.1) +
  draw_plot(p_net, 0.22, 0.88, 0.17, 0.08) +
  draw_plot(p_B_1, 0, 0.33, 0.8, 0.22) +
  draw_plot(p_B_2, 0.58, 0.54, 0.1, 0.1) +
  draw_plot(p_net, 0.22, 0.55, 0.17, 0.08) +
  draw_plot(p_C_1, 0, 0, 0.8, 0.22) +
  draw_plot(p_C_2, 0.58, 0.21, 0.1, 0.1) +
  draw_plot(p_net, 0.22, 0.22, 0.17, 0.08) +
  draw_plot(legend, 0.8, 0.4, 0.25, 0.15) +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0), c(0.98, 0.65, 0.32), size = 24)

# saving
save_plot("fig1.pdf", fig1, ncol = 1, nrow = 3, base_aspect_ratio = 1.8)

#-----------------------------------------------------------------------------------------------------#