#-----------------------------------------------------------------------------------------------------#

# Description: 
#   
#
# Returns:
#   

library(ggplot2)
library(cowplot)
library(viridis)

source("functions/MatchingMutNet.R")

# creating matrix with 3 species
n_row = 2
n_col = 1
n_sp = n_row + n_col
mat = matrix(c(1, 1), nrow = n_row, ncol = n_col)
f = rbind(cbind(matrix(0, n_row, n_row), mat),
          cbind(t(mat), matrix(0, n_col, n_col)))

theta = c(1.5, 9, 5)

#-----------------------------------------------------#
# Fig. 1A: intermediate trait matching - no gene flow #
#-----------------------------------------------------#

z = c(3.7, 6, 4.7)

# trait distributions
p_A_1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = z[1], sd = 1.2),
                color = "gray80", size = 1) +
  stat_function(fun = dnorm, args = list(mean = z[2], sd = 1.2),
                color = "black", size = 1) +
  stat_function(fun = dnorm, args = list(mean = z[3], sd = 1.2),
                color = "gray50", size = 1) +
  geom_vline(xintercept = theta[1], color = "gray80", size = 0.8,
             linetype = "dashed") +
  geom_vline(xintercept = theta[2], color = "black", size = 0.8,
             linetype = "dashed") +
  geom_vline(xintercept = theta[3], color = "gray50", size = 0.8,
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 17)) +
  xlab("Trait value") +
  ylab("Density") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        legend.position = "none")

# trait matching heatmap
matching = MatchingMutNet(n_sp, n_row, n_col, f, z = z,
                          method = "exponential", alpha = 0.2)

p_A_2 = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", direction = -1,
                     limits = c(0.15, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.6) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        axis.ticks = element_blank(),
        legend.position = "none")
  
#----------------------------------------------#
# Fig. 1B: low trait matching - with gene flow #
#----------------------------------------------#

z = c(8, 13, 10)

# trait distributions
p_B_1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = z[1], sd = 1.2),
                color = "gray80", size = 1) +
  stat_function(fun = dnorm, args = list(mean = z[2], sd = 1.2),
                color = "black", size = 1) +
  stat_function(fun = dnorm, args = list(mean = z[3], sd = 1.2),
                color = "gray50", size = 1) +
  geom_vline(xintercept = theta[1], color = "gray80", size = 0.8,
             linetype = "dashed") +
  geom_vline(xintercept = theta[2], color = "black", size = 0.8,
             linetype = "dashed") +
  geom_vline(xintercept = theta[3], color = "gray50", size = 0.8,
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 17)) +
  xlab("Trait value") +
  ylab("Density") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        legend.position = "none")

# trait matching heatmap
matching = MatchingMutNet(n_sp, n_row, n_col, f, z = z,
                          method = "exponential", alpha = 0.2)

p_B_2 = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                         height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", direction = -1,
                     limits = c(0.15, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.6) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        axis.ticks = element_blank(),
        legend.position = "none")

#-----------------------------------------------#
# Fig. 1C: high trait matching - with gene flow #
#-----------------------------------------------#

z = c(9.8, 10.5, 10)

# trait distributions
p_C_1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, args = list(mean = z[1], sd = 1.2),
                color = "gray80", size = 1) +
  stat_function(fun = dnorm, args = list(mean = z[2], sd = 1.2),
                color = "black", size = 1) +
  stat_function(fun = dnorm, args = list(mean = z[3], sd = 1.2),
                color = "gray50", size = 1) +
  geom_vline(xintercept = theta[1], color = "gray80", size = 0.8,
             linetype = "dashed") +
  geom_vline(xintercept = theta[2], color = "black", size = 0.8,
             linetype = "dashed") +
  geom_vline(xintercept = theta[3], color = "gray50", size = 0.8,
             linetype = "dashed") +
  scale_x_continuous(limits = c(0, 17)) +
  xlab("Trait value") +
  ylab("Density") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        legend.position = "none")

# trait matching heatmap
matching = MatchingMutNet(n_sp, n_row, n_col, f, z = z,
                          method = "exponential", alpha = 0.2)

p_C_2 = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                         height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", direction = -1,
                     limits = c(0.15, 1), name = "Trait\nmatching") +
  geom_tile(color = "black", size = 0.7) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 17),
        axis.ticks = element_blank(),
        legend.position = "none")

# a dummie plot just to get the legend
p_legend = ggplot(data = matching[[2]], aes(x = column, y = row, fill = matching,
                                         height = 0.95, width = 0.9)) + 
  scale_fill_viridis(option = "plasma", direction = -1,
                     limits = c(0.15, 1), name = "trait\nmatching") +
  geom_tile(color = "black", size = 0.6) +
  xlab("") +
  ylab("") +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.key.size = unit(0.55, "cm"))

legend = get_legend(p_legend)

# generating the whole figure
fig1 = ggdraw() +
  draw_plot(p_A_1, 0, 0.66, 0.8, 0.22) +
  draw_plot(p_A_2, 0.57, 0.82, 0.2, 0.17) +
  draw_plot(p_B_1, 0, 0.33, 0.8, 0.22) +
  draw_plot(p_B_2, 0.57, 0.49, 0.2, 0.17) +
  draw_plot(p_C_1, 0, 0, 0.8, 0.22) +
  draw_plot(p_C_2, 0.57, 0.16, 0.2, 0.17) +
  draw_plot(legend, 0.8, 0.4, 0.25, 0.15) +
  draw_plot_label(c("A", "B", "C"), c(0, 0, 0), c(0.98, 0.65, 0.32), size = 24)

# saving
save_plot("fig1.pdf", fig1, ncol = 1, nrow = 3, base_aspect_ratio = 1.8)

#-----------------------------------------------------------------------------------------------------#