require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(tidyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)


sim.nona <- sim.spectral [-c(5, 91)]

eval.df <-
    ldply(sim.nona, function(L1)
        ldply(L1, function(L2)
        {
            eval <- Re(L2 $ Q.eval)
            cbind(1:nrow(eval), eval)
        }))

colnames(eval.df) <- c('init', 'scen', paste('eval', 1:45))

eval.df <- data.frame(par.table[eval.df $ scen, ], eval.df)

colnames(eval.df) [1:3] <- c('gflow', 'hotA', 'hotB')

eval.gdf <- gather(eval.df)
