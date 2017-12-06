require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)

## load('Arroyo1982.RData')

colnames(par.table) <- c('g', 'mA', 'mB')

Arroyo.spec <- sim.spectral

Arroyo.Tmat <- 
    llply(Arroyo.spec,
          function(L1) llply(L1 [[1]],
                             function(L2) if(length(L2) > 1) return(L2 $ T.eq)))

llply(Arroyo.Tmat, function(L1)


Arroyo.Tmat <-
    lply(Arroyo.spec, function(L1)
        laply(L1 [[1]], function(L2) L2 $ T.eq))

dimnames(Arroyo.Tmat) [[3]] <- NULL

Ols.eval <- aaply(Arroyo.Tmat, c(1, 2), function(T) eigen(T) $ values)

Ols.eval [, , 2] #init, rep, eval

names(dimnames(Ols.eval)) <- c('init', 'rep', 'eval')

eval.df <- adply(Ols.eval, 1:3)

eval.df <- mutate(eval.df, 'Re' = Re(V1), 'Im' = Im(V1)),

eval.df <- mutate (eval.df, 'par' = par.table [as.numeric(eval.df $ rep), ])

plot(Im ~ Re, data = eval.df)

()
    
