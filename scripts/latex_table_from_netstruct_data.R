#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Creates a LaTeX table using the data frame of network structure metrics
#
# Returns:
#   A tex file containing the table

# loading package
library(xtable)

# network structure data frame
network_structure = read.csv("output/data/network_structure/network_structure.csv")

# changing mutualism names
network_structure$mutualism = as.character(network_structure$mutualism)
network_structure$mutualism[network_structure$mutualism == "AE"] = "AN"
network_structure$mutualism[network_structure$mutualism == "CC"] = "MC"
network_structure$mutualism[network_structure$mutualism == "FA"] = "AA"
network_structure$mutualism[network_structure$mutualism == "FP"] = "SD"
network_structure$mutualism[network_structure$mutualism == "PP"] = "P"

# removing some variables
network_structure = network_structure[ , c(1, 2, 3, 4, 5, 7, 9, 12, 15, 16)]
# ordering according to mutualism type
network_structure = network_structure[order(network_structure$mutualism), ]

# creating latex table
tex_table = xtable(network_structure)

# sving file
print.xtable(tex_table, type = "latex", file = "Table_S1.tex")

#-----------------------------------------------------------------------------------------------------#