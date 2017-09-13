require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)

load('Olesen2002.RData')

colnames(par.table) <- c('g', 'mA', 'mB')

Olesen.spec <- sim.spectral

llply(Olesen.spec [[1]], is.na)

Olesen.Tmat <-
    laply(Olesen.spec, function(L1)
        laply(L1 [[1]] [1:96], function(L2) L2 $ T.eq))

dimnames(Olesen.Tmat) [[3]] <- NULL

Ols.eval <- aaply(Olesen.Tmat, c(1, 2), function(T) eigen(T) $ values)

Ols.eval [, , 2] #init, rep, eval

names(dimnames(Ols.eval)) <- c('init', 'rep', 'eval')

eval.df <- adply(Ols.eval, 1:3)

eval.df <- mutate(eval.df, 'Re' = Re(V1), 'Im' = Im(V1)),

eval.df <- mutate (eval.df, 'par' = par.table [as.numeric(eval.df $ rep), ])

plot(Im ~ Re, data = eval.df)

()
    
