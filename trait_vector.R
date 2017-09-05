require(doMC)
registerDoMC(cores = 32)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)


### load networks
binary.files <- dir('data/empirical_networks/binary/')
binary.networks <- list()

for(i in 1:length(binary.files))
    binary.networks [[i]] <-
        as.matrix(read.table(paste0('data/empirical_networks/binary/', binary.files [i])))

names(binary.networks) <- gsub('_bin.txt', '', binary.files)

### load functions
source('functions/twoSitesVectorMatch.R')
source('functions/assembleQ.R')

## choose network
net <- binary.networks [[61]] # Olesen 2002 (10 x 12)

## scenarios to be tested
flow <- seq(0.01, 0.1, by = 0.03)
hotA <- seq(0.2, 1, by = 0.2)
hotB <- hotA

## other pars
temp <- 0.1
h <- 0.1

thA <- rnorm(sum(dim(net)), 30, sqrt(40))
thB <- rnorm(sum(dim(net)), 70, sqrt(40))

### all combinations of simulation values
par.table <- as.matrix(expand.grid(flow, hotA, hotB))

## one hundred simulations per combination of g, mA, mB
it <- 100

## massive set of lists within lists
sim.spectral <-
    alply(1:it, 1, function(i)
    {
        iA <- runif(sum(dim(net)), 10, 100)
        iB <- runif(sum(dim(net)), 10, 100)

        out <- alply(1:nrow(par.table), 1, function(j)
        {
            tryCatch(expr = twoSitesVectorMatch(net, par.table[j, 1],
                                                h, temp,
                                                thA, thB,
                                                par.table[j, 2],
                                                par.table[j, 3], iA, iB),
                     error = function(cond) return(NA))
            
        }, .parallel = TRUE)

        print(i)

        out
    })

save(sim.spectral, par.table, file = 'Olesen2002_test.RData')

str(sim.spectral)
