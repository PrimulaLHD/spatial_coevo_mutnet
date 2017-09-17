require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)

require(viridis)

## load('Arroyo1982.RData')

colnames(par.table) <- c('g', 'mA', 'mB')

Arroyo.spec <- sim.spectral

Arroyo.Tmat <- 
    llply(Arroyo.spec,
          function(L1) llply(L1 [[1]],
                             function(L2) if(length(L2) > 1) return(L2 $ T.eq)))

Arroyo.Tmat <- llply(Arroyo.Tmat, function(L1) L1 [!laply(L1, is.null)])
Arroyo.Tmat <- llply(Arroyo.Tmat, function(L1) laply(L1, identity))

Arroyo.Tmat <- laply(Arroyo.Tmat, identity)

Ayo.match <- 
    llply(Arroyo.spec,
          function(L1) llply(L1 [[1]],
                             function(L2) if(length(L2) > 1) return(L2 $ av.match)))

Ayo.match <- llply(Ayo.match, function(L1) L1 [!laply(L1, is.null)])
Ayo.match <- llply(Ayo.match, function(L1) laply(L1, identity))

Ayo.match <- laply(Ayo.match, identity)

dimnames(Arroyo.Tmat) [[3]] <- NULL

Ayo.eval <- aaply(Arroyo.Tmat, c(1, 2), function(T) eigen(T) $ values, .parallel = TRUE)

Ayo.evar <- aaply(Ayo.eval, c(1, 2), function(ev) var(Mod(ev)))

names(dimnames(Ayo.eval)) <- c('init', 'rep', 'eval')

names(dimnames(Ayo.evar)) <- c('init', 'rep')

eval2.df <- adply(Ayo.eval [, , 2], 1:2)

evar.df <- adply(Ayo.evar, 1:2)

match.df <- adply(Ayo.match, 1:2)

head(match.df)

metrics.df <- cbind(eval2.df, evar.df [, 3], match.df [, 3:4])

colnames(metrics.df) <- c('init', 'rep', 'eval2', 'evar', 'avmA', 'avmB')

metrics.df <- cbind(par.table [metrics.df $ rep, ], metrics.df)

metrics.df $ g <- factor(metrics.df $ g)
metrics.df $ mA <- factor(metrics.df $ mA)
metrics.df $ mB <- factor(metrics.df $ mB)

metrics.df %>%
    mutate('Mod2' = Mod(eval2),
           'mean.match' = rowMeans(cbind(avmA, avmB)),
           'scenario' = paste0('mA = ', mA, ', mB = ', mB)) %>%
    ggplot(.) +
    geom_point(aes(x = avmB, y = evar, color = mA, shape = mB)) +
    facet_wrap(~ g) +
    scale_color_viridis(option = 'C', direction = -1, discrete = TRUE) +
    theme_bw()

dev.off(dev.cur())
