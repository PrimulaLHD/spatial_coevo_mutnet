#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Characterizes network structure for all networks in the dataset.
#
# Returns:
#   Calculates: species richness, connectance, degree variance,
#   nestedness and standardized nestedness (using null model 2), and performs a principal components
#   analysis with all network metrics.

# loading packages and functions 
library(bipartite)
library(plyr)
source("functions/NullModel2.R")

# getting network file names
net_files = list.files("~/Dropbox/spatial_coevo_mutnet/matrices/binary/")
net_names = substr(net_files, start = 1, stop = nchar(net_files) - 4) # removing the extension .txt

# mutualism types
mutualism = substr(net_names, start = 1, stop = 2)

# getting results from modular
results_modular = read.csv("~/Dropbox/spatial_coevo_mutnet/results/network_structure/Modularity_Barber_SA_Null2.csv")
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
interactions = rep(NA, length(net_names))
connectance = rep(NA, length(net_names))
degree_variance = rep(NA, length(net_names))
nestedness = rep(NA, length(net_names))
std_nestedness = rep(NA, length(net_names))
p_nestedness = rep(NA, length(net_names))

for (i in 1:length(net_names)) {
  net = as.matrix(read.table(paste("~/Dropbox/spatial_coevo_mutnet/matrices/binary/", net_files[i], sep = ""), 
                             sep = " ", header = FALSE))
  network[i] = net_names[i]
  rows[i] = nrow(net)
  columns[i] = ncol(net)
  richness[i] = nrow(net) + ncol(net)
  interactions[i] = sum(net)
  connectance[i] = sum(net)/(nrow(net)*ncol(net))
  
  # degree variance
  row_degree = apply(net, 1, sum)
  col_degree = apply(net, 2, sum)
  degrees = c(row_degree, col_degree)
  degree_variance[i] = var(degrees)
  
  # nestedness
  nestedness[i] = nested(net, method = "NODF2")
  null_nets = NullModel2(mat = net, mat_name = net_names[i],
                         write = FALSE, R = 100) # generating networks according to null model 2
  null_nestedness = ldply(null_nets, nested, method = "NODF2") # calculating nestedness for all networks
  mean_null = mean(null_nestedness$NODF2) # mean nestedness for theoretical matrices
  sd_null = sd(null_nestedness$NODF2) # nestedness sd for theoretical matrices
  std_nestedness[i] = (nestedness[i] - mean_null)/sd_null # standardized nestedness
  if (is.nan(std_nestedness[i])) # if sd is zero, std nestedness is zero
    std_nestedness[i] = 0
  p_nestedness[i] = sum(null_nestedness$NODF2 >= nestedness[i])/nrow(null_nestedness) # calculating p value
}

# Principal components analysis
net_struct = data.frame(network, mutualism, rows, columns, richness, interactions, connectance, degree_variance,
                        nestedness, std_nestedness, p_nestedness, modularity, std_modularity, p_modularity)
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
write.csv(net_struct, file = "~/Dropbox/spatial_coevo_mutnet/results/network_structure/network_structure.csv", row.names = FALSE)

#-----------------------------------------------------------------------------------------------------#