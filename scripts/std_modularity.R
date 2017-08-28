#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Calculates standardized modularity using results from program MODULAR.
#
# Returns:
#   Saves a csv file with modularity values and p values.

setwd("~/Dropbox/spatial_coevo_mutnet/results/network_structure/modularity_resultsSA_Barber/")
files = dir()

mod = read.table("OUT_MOD.txt", header = TRUE)
net_names = substr(as.character(mod$File), start = 1, stop = nchar(as.character(mod$File)) - 4)
mutualism = substr(net_names, start = 1, stop = 2)

out_2 = files[grep("OUT_2", files)]
net_names2 = substr(out_2, start = 6, stop = nchar(out_2) - 4)

net_names == net_names2

std_mod = rep(NA, length(net_names))

for (i in 1:length(net_names)) {
  nullmod_results = read.table(out_2[i], header = TRUE)
  obs_mod = mod$Modularity[i]
  mean_null_mod = mean(nullmod_results$Modularity)
  sd_null_mod = sd(nullmod_results$Modularity)
  std_mod[i] = (obs_mod - mean_null_mod)/sd_null_mod
  if (sd_null_mod == 0) 
    std_mod[i] = 0
}

results_df = data.frame(network = net_names, mutualism = mutualism, obs_mod = mod$Modularity,
                        std_mod = std_mod, p = mod$P.Null2)

write.csv(results_df, file = "~/Dropbox/spatial_coevo_mutnet/results/network_structure/Modularity_Barber_SA_Null2.csv",
          row.names = FALSE)

#-----------------------------------------------------------------------------------------------------#