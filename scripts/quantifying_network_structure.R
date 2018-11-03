#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Characterizes network structure for all networks in the dataset.
#
# Returns:
#   Calculates: species richness, connectance, and nestedness, and 
#   performs a principal components analysis with all network metrics.

# loading packages and functions 
library(bipartite)
library(plyr)

# getting network file names
net_files = list.files("data/empirical_networks/binary/")
net_names = substr(net_files, start = 1, stop = nchar(net_files) - 4) # removing the extension .txt

# mutualism types
mutualism = substr(net_names, start = 1, stop = 2)

# getting results from modular
results_modular = read.csv("output/data/network_structure/Modularity_Barber_SA_Null2.csv")
# checking if network names are in the same order
results_modular$network == net_names
modularity = results_modular$obs_mod
std_modularity = results_modular$std_mod
p_modularity = results_modular$p

# empty vectors
network = rep(NA, length(net_names))
rows = rep(NA, length(net_names))
columns = rep(NA, length(net_names))
richness = rep(NA, length(net_names))
connectance = rep(NA, length(net_names))
nestedness = rep(NA, length(net_names))

for (i in 1:length(net_names)) {
  net = as.matrix(read.table(paste("data/empirical_networks/binary/", net_files[i], sep = ""), 
                             sep = " ", header = FALSE))
  network[i] = net_names[i]
  rows[i] = nrow(net)
  columns[i] = ncol(net)
  richness[i] = nrow(net) + ncol(net)
  connectance[i] = sum(net)/(nrow(net)*ncol(net))
  nestedness[i] = nested(net, method = "NODF2")
}

# Principal components analysis
net_struct = data.frame(network, mutualism, rows, columns, richness, connectance,
                        nestedness, modularity)
mat = as.matrix(net_struct[ , c("richness", "connectance", "nestedness", "modularity")])
pca = prcomp(x = mat, scale = TRUE) # performing pca
pc_sd = pca[[1]]
pc_eigens = pc_sd^2
rel_PC_eigens = round(pc_eigens/sum(pc_eigens), 2) # percentage of variance explained by each PC
summary(pca) # same information as last lines

# adding PC1 and PC2 to data frame
net_struct$PC1 = as.data.frame(pca$x)$PC1
net_struct$PC2 = as.data.frame(pca$x)$PC2

# saving all results
write.csv(net_struct, file = "output/data/network_structure/network_structure.csv", row.names = FALSE)

#-----------------------------------------------------------------------------------------------------#