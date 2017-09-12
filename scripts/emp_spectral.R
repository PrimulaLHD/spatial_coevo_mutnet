require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)

### load weighted networks
weighted.files <- dir('data/empirical_networks/weighted/')
weighted.networks <- list()

for(i in 1:length(weighted.files))
    weighted.networks [[i]] <-
        as.matrix(read.table(paste0('data/empirical_networks/weighted/', weighted.files [i])))

names(weighted.networks) <- gsub('_wgt.txt', '', weighted.files)

### load function
source('functions/empiricalSpectralT.R')
