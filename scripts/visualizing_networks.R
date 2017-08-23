#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Different plots to visualize networks.
#
# Returns:
#   Several different plots.

library(RColorBrewer)
library(ggplot2)
library(GGally)
library(network)
library(sna)
library(scales)
library(reshape2)
library(cowplot)
library(fields)

source("functions/RemoveNodesNoLinks.R")

# creating adjacency matrix
set.seed(33)
n_row = 20
n_col = 20
n_sp = n_row + n_col
mat = matrix(data = sample(c(0, 1), size = n_row*n_col, replace = TRUE, prob = c(0.8, 0.2)),
             nrow = n_row, ncol = n_col)
# adding row and column names
rownames(mat) = paste("R", 1:n_row, sep = "")
colnames(mat) = paste("C", 1:n_col, sep = "")
mat = RemoveNodesNoLinks(mat)
n_row = nrow(mat)
n_col = ncol(mat)
n_sp = n_row + n_col

# creating network object for plotting
net = network(mat, matrix.type = "bipartite", directed = FALSE, ignore.eval = FALSE)

# network with black and grey nodes
set.seed(42)
ggnet2(net, node.size = 10, color = c(rep("black", n_row), rep("gray70", n_col)), 
       edge.color = "black", edge.size = 1, mode = "kamadakawai")

# network with yellow and blue nodes
set.seed(42)
ggnet2(net, node.size = 10, color = c(rep("#EBCC2A", n_row), rep("#78B7C5", n_col)), 
       edge.color = "black", edge.size = 1, mode = "kamadakawai")
# adding node labels
set.seed(42)
ggnet2(net, node.size = 10, color = c(rep("#EBCC2A", n_row), rep("#78B7C5", n_col)), 
       edge.color = "black", edge.size = 1, mode = "kamadakawai", label = TRUE)
# adding different shapes for different groups
set.seed(42)
ggnet2(net, node.size = 10, color = c(rep("#EBCC2A", n_row), rep("#78B7C5", n_col)), 
       edge.color = "black", edge.size = 1, mode = "kamadakawai", label = TRUE,
       shape = c(rep(15, n_row), rep(16, n_col)))
# changing line type for links
set.seed(42)
ggnet2(net, node.size = 10, color = c(rep("#EBCC2A", n_row), rep("#78B7C5", n_col)), 
       edge.color = "black", edge.size = 1, mode = "kamadakawai", label = TRUE, 
       shape = c(rep(15, n_row), rep(16, n_col)), edge.lty = "dashed")

# network with a different color for each node
pal_row = seq_gradient_pal(low = "gold1", 
                           high = "firebrick3")(seq(0, 1, length.out = n_row)) 
pal_col = seq_gradient_pal(low = "paleturquoise1", 
                           high = "forestgreen")(seq(0, 1, length.out = n_col))
set.seed(42)
ggnet2(net, node.size = 10, color = c(pal_row, pal_col), 
       edge.color = "black", edge.size = 1, mode = "kamadakawai", label = TRUE,
       shape = c(rep(15, n_row), rep(16, n_col)))

# network with edges colored according to quantitative attribute
matching = runif(sum(mat), 0, 1) # hypothetical trait matching values
heat = rev(heat.colors(100))[20:100]
bins = seq(0, 1, length = length(heat))
heat_pal = heat[as.numeric(cut(matching, breaks = bins))]
set.seed(42)
ggnet2(net, node.size = 10, color = "black", shape = c(rep(15, n_row), rep(16, n_col)),
       edge.color = heat_pal, edge.size = 2.5, mode = "kamadakawai")

# black and white adjacency matrix 
order_col = order(apply(mat, 2, sum), decreasing = TRUE)
order_row = order(apply(mat, 1, sum))
mat_ordered = mat[order_row, order_col]
melted_mat = melt(mat_ordered)
melted_mat$color[melted_mat$value == 1] = "black"
melted_mat$color[melted_mat$value == 0] = "white"
ggplot(data = melted_mat, aes(x = Var2, y = Var1, fill = value,
                              height = 0.75, width = 0.75)) + 
  geom_tile(fill = melted_mat$color) +
  xlab("Columns") +
  ylab("Rows")

# adjacency matrix colored according to quantitative attribute
melted_mat$matching = rep(0, nrow(melted_mat))
melted_mat$matching[melted_mat$value == 1] = matching
ggplot(data = subset(melted_mat, matching != 0), 
       aes(x = Var2, y = Var1, fill = matching,
           height = 0.75, width = 0.75)) + 
  geom_tile() +
  scale_fill_distiller(palette = "Spectral") +
  xlab("Columns") +
  ylab("Rows")

#-----------------------------------------------------------------------------------------------------#