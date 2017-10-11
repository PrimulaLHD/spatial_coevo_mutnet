require(plyr)
require(dplyr)
require(doMC)
registerDoMC(cores = 3)

source('functions/MatchingMutNet.R')

### load weighted networks
weighted.files <- dir('data/empirical_networks/weighted/')
weighted.networks <- list()

for(i in 1:length(weighted.files))
    weighted.networks [[i]] <-
        as.matrix(read.table(paste0('data/empirical_networks/weighted/', weighted.files [i])))

names(weighted.networks) <- gsub('_wgt.txt', '', weighted.files)

## fixed parameters
g <- 0.3
thetaAmin <- 0
thetaAmax <- 10
thetaBmin <- 10
thetaBmax <- 20

## parameters
mA_seq <- seq(0.1, 0.9, 0.2)
mB_seq <- seq(0.1, 0.9, 0.2)
pert_g_seq <- seq(0, 1, by = 0.05)

weight.factor <- factor(names(weighted.networks))

m.table <- as.matrix(expand.grid(mA_seq, mB_seq, pert_g_seq))
m.table <- m.table[m.table [, 1] >= m.table [, 2], ]

## run

perturb.results <-
    ldply(weighted.networks,
          function(L)
          {
              adply(m.table, 1, function(line)
              {
                  results.mat <- perturbEmpT(L, line [1], line [2], line [3],
                                         n_rep_pert = 30, n_theta = 30)
                  
                  par.tab <- matrix(rep(line, times = nrow(results.mat)),
                                    byrow = TRUE, ncol = 3)

                  colnames(par.tab) <- c('mA', 'mB', 'pert_g_seq')
                  
                  out <- cbind(par.tab, results.mat)

                  out
                  
              }, .parallel = TRUE)
          })

