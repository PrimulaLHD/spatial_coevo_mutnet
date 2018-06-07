#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Creates a LaTeX table with the correlations between PC1, PC2 and network metrics
#
# Returns:
#   A tex file containing the table

# loading package
library(xtable)

# Principal components analysis
net_struct = read.csv("output/data/network_structure/network_structure.csv")
mat = as.matrix(net_struct[ , c("richness", "connectance", "nestedness", "modularity")])
pca = prcomp(x = mat, scale = TRUE) # performing pca
summary(pca)
table = pca[[2]]

# creating latex table
tex_table = xtable(table)

# sving file
print.xtable(tex_table, type = "latex", file = "Table_S3.tex")

#-----------------------------------------------------------------------------------------------------#