require(doMC)
registerDoMC(cores = 32)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)

setwd('~/spatial_coevo_mutnet/')

Norm <- function(x) sqrt(sum (x * x))
Normalize <- function(x) x / Norm (x)

### load binary networks
binary.files <- dir('data/empirical_networks/binary/')
binary.networks <- list()

for(i in 1:length(binary.files))
    binary.networks [[i]] <-
        as.matrix(read.table(paste0('data/empirical_networks/binary/', binary.files [i])))

names(binary.networks) <- gsub('_bin.txt', '', binary.files)

### load functions
source('functions/twoSiteSpectral.R')

source('functions/assembleQ.R')
source('functions/convMutNet_2.R')

## choose network
net <- binary.networks [[45]] # Olesen 2002 (10 x 12)

## scenarios to be tested
flow <- seq(0, 0.4, by = 0.05)
hotA <- seq(0.2, 1, by = 0.2)
hotB <- hotA

## other pars
temp <- 0.2
h <- 1

### all combinations of simulation values
par.table <- as.matrix(expand.grid(flow, hotA, hotB))

par.table <- par.table [par.table [, 2] >= par.table [, 3], ]

## one hundred simulations per combination of g, mA, mB
it <- 100

set.seed(987)

## massive set of lists within lists
sim.spectral <-
    alply(1:it, 1, function(i)
    {
        iA <- runif(sum(dim(net)), 10, 30)
        iB <- runif(sum(dim(net)), 10, 30)

        thA <- runif(sum(dim(net)), 10, 20)
        thB <- runif(sum(dim(net)), 20, 30)
        
        out <- alply(1:nrow(par.table), 1, function(j)
        {
            tryCatch(expr = twoSiteSpectral(net, par.table[j, 1],
                                             h, temp,
                                             thA, thB,
                                             par.table[j, 2],
                                             par.table[j, 3], iA, iB),
                     error = function(cond) return(NA))
            
        }, .parallel = TRUE)

        print(i)

        return(list(out, thA, thB))
    })

save(sim.spectral, par.table, file = 'Arroyo1982.RData')

