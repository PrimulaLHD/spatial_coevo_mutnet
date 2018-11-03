#-----------------------------------------------------------------------------------------------------#

# Description: 
#   Reads summaries from coevolution simulations and generates a plot of trait matching as a function
#   of gene flow for several parameterizations of the model.
#
# Returns:
#   Plot.

# loading functions and packages
library(ggplot2)
library(plyr)
library(cowplot)
library(viridis)

# getting parameter names
pars = dir("output/data/simulations_empirical_networks/sensitivity_analysis/")

# to store results
summ_final_mut_mat_df = c()

for (i in 1:length(pars)) {
  # files for parameter i
  par_char = dir(paste("output/data/simulations_empirical_networks/sensitivity_analysis/",
                       pars[i], sep = ""))
  for (j in 1:length(par_char)) {
    # simulation result files for parameter value j
    summaries = dir(paste("output/data/simulations_empirical_networks/sensitivity_analysis/",
                          pars[i], "/", par_char[j], "/summary_coevo_results", sep = ""))
    # getting file names for mA = mB = 0.7 and mA = 0.9, mB = 0.1
    summaries_m = summaries[c(grep("mA07_mB07", summaries), grep("mA09_mB01", summaries))]
    # reading files and calculating mean trait matching for many simulations
    for (k in 1:length(summaries_m)) {
      # extracting parameter values from file name
      split_mA = sub(paste(".*", "mA", sep = ""), "", summaries_m[k])
      mA_char = sub("_.*", "", split_mA)
      mA = as.numeric(paste("0", ".", substr(mA_char, start = 2, stop = nchar(mA_char)), sep = ""))
      split_mB = sub(paste(".*", "mB", sep = ""), "", summaries_m[k])
      mB_char = sub("_.*", "", split_mB)
      mB = as.numeric(paste("0", ".", substr(mB_char, start = 2, stop = nchar(mB_char)), sep = ""))
      split_alpha = sub(paste(".*", "alpha", sep = ""), "", summaries_m[k])
      alpha_char = sub("_.*", "", split_alpha)
      alpha = as.numeric(paste("0", ".", substr(alpha_char, start = 2, stop = nchar(alpha_char)), sep = ""))
      split_phi = sub(paste(".*", "phi", sep = ""), "", summaries_m[k])
      phi_char = sub("_.*", "", split_phi)
      phi = as.numeric(phi_char)
      if (substring(phi_char, 1, 1) == "0")
        phi = as.numeric(paste("0", ".", substr(phi_char, start = 2, stop = nchar(phi_char)), sep = ""))
      split_thetaA = sub(paste(".*", "thetaA", sep = ""), "", summaries_m[k])
      thetaA_char = sub("_.*", "", split_thetaA)
      split_thetaB = sub(paste(".*", "thetaB", sep = ""), "", summaries_m[k])
      thetaB_char = substr(split_thetaB, start = 1, stop = nchar(split_thetaB) - 4)
      
      # reading file with simulation results
      current_df = read.csv(paste("output/data/simulations_empirical_networks/sensitivity_analysis/",
                                  pars[i], "/", par_char[j], "/summary_coevo_results/", summaries_m[k], sep = ""))
      
      # adding parameter values
      current_df$m_A = rep(mA, nrow(current_df))
      current_df$m_B = rep(mB, nrow(current_df))
      current_df$alpha = rep(alpha, nrow(current_df))
      current_df$phi = rep(phi, nrow(current_df))
      current_df$theta_A = rep(thetaA_char, nrow(current_df))
      current_df$theta_B = rep(thetaB_char, nrow(current_df))
      
      # calculating mean trait matching
      curr_summ_final_mut_mat_df = ddply(current_df, c("g", "network", "site"),
                                         summarize, mean_final_mut_mat = mean(mean_final_mut_matching),
                                         mutualism = unique(mutualism),
                                         m_A = unique(m_A),
                                         m_B = unique(m_B),
                                         alpha = unique(alpha),
                                         phi = unique(phi),
                                         theta_A = unique(theta_A),
                                         theta_B = unique(theta_B))
      
      # joining with previous results
      summ_final_mut_mat_df = rbind(summ_final_mut_mat_df, curr_summ_final_mut_mat_df)
    }
  }
}

# saving/reading results
write.csv(summ_final_mut_mat_df, row.names = FALSE, 
          file = "output/data/simulations_empirical_networks/summary_coevo_results/summary_coevo_results_final_mut_matching_sensitivity_analysis.csv")
summ_final_mut_mat_df = read.csv("output/data/simulations_empirical_networks/summary_coevo_results/summary_coevo_results_final_mut_matching_sensitivity_analysis.csv")

# changing mutualism names
summ_final_mut_mat_df$mutualism = as.character(summ_final_mut_mat_df$mutualism)
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AE"] = "ants - nectary-bearing plants"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "AM"] = "ants - myrmecophytes"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "CC"] = "marine cleaning"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FA"] = "anemones - anemonefishes"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "FP"] = "seed dispersal"
summ_final_mut_mat_df$mutualism[summ_final_mut_mat_df$mutualism == "PP"] = "pollination"

# defining a color palette and ordering according to factor mutualism
orange = "#FF7F00"
red = "#CD2626"
cyan = "#00EEEE"
blue = "#4876FF"
green = "#1AB245"
purple = "#A020F0"
my_palette = c(orange, red, cyan, blue, green, purple)

# Fig A: mean final trait matching as a function of gene flow for site A with theta_A ~ [0, 10], theta_B ~ [0, 10]  
p_A = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & theta_A == "0-10" & theta_B == "0-10"),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.9) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(theta[iA], " ~ U[0, 10], ", theta[iB], " ~ U[0, 10]"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig B: mean final trait matching as a function of gene flow for site A with theta_A ~ [0, 10], theta_B ~ [5, 15]  
p_B = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & theta_A == "0-10" & theta_B == "5-15"),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(theta[iA], " ~ U[0, 10], ", theta[iB], " ~ U[5, 15]"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig C: mean final trait matching as a function of gene flow for site A with theta_A ~ [0, 10], theta_B ~ [20, 30]  
p_C = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & theta_A == "0-10" & theta_B == "20-30"),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(theta[iA], " ~ U[0, 10], ", theta[iB], " ~ U[20, 30]"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig D: mean final trait matching as a function of gene flow for site A with theta_A ~ [0, 10], theta_B ~ theta_A + 10  
p_D = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & theta_A == "0-10" & theta_B == "10-20propA"),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(theta[iA], " ~ U[0, 10], ", theta[iB], " = ", theta[iA], " + N[10, 1]"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig E: mean final trait matching as a function of gene flow for site A with phi = 0.1
p_E = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & phi == 0.1),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(varphi[iA], ", ", varphi[iB], " ~ N[0.1, ", 10^-4, "]"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig F: mean final trait matching as a function of gene flow for site A with phi = 1
p_F = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & phi == 1),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(varphi[iA], ", ", varphi[iB], " ~ N[1, ", 10^-4, "]"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig G: mean final trait matching as a function of gene flow for site A with alpha = 0.05
p_G = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & alpha == 0.05),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(alpha, " = 0.05"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

# Fig H: mean final trait matching as a function of gene flow for site A with alpha = 0.8
p_H = ggplot(data = subset(summ_final_mut_mat_df, site == "A" & alpha == 0.8),
             aes(x = g, y = mean_final_mut_mat, color = mutualism)) +
  geom_point(size = 1.5, shape = 19, alpha = 0.6) + 
  scale_color_manual(values = my_palette) +
  ggtitle(expression(paste(alpha, " = 0.8"))) +
  xlab("") +
  ylab("") +
  facet_wrap(m_A ~ m_B, scales = "free", labeller = label_both) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.position = "none")

figSI = plot_grid(p_A, p_B, p_C, p_D, p_E, p_F, p_G, p_H, labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                   nrow = 4, ncol = 2, label_size = 20)

save_plot("output/figs/figS5.pdf", figSI, ncol = 2, nrow = 4, base_aspect_ratio = 1.4)

#-----------------------------------------------------------------------------------------------------#