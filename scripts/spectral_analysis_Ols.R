require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(RColorBrewer)
require(viridis)
load('Olesen2002.RData')

colnames(par.table) <- c('g', 'mA', 'mB')

Olesen.spec <- sim.spectral

### scenarios 97-100 with errors
Olesen.Tmat <-
    laply(Olesen.spec, function(L1)
        laply(L1 [[1]] [1:96], function(L2) L2 $ T.eq))

dimnames(Olesen.Tmat) [[3]] <- NULL

Ols.eval <- aaply(Olesen.Tmat, c(1, 2), function(T) eigen(T) $ values)

Ols.evar

names(dimnames(Ols.eval)) <- c('init', 'rep', 'eval')

## have some non zero im part

Ols.match <- 
    laply(Olesen.spec, function(L1)
        laply(L1 [[1]] [1:96], function(L2) L2 $ av.match))

match.df <- adply(Ols.match, 1:2)

eval2.df <- adply(Ols.eval [, , 2], 1:2)

eval.df <-  adply(Ols.eval, 1:3)

metrics.df <- cbind(eval2.df, match.df [, 2:3])

colnames(metrics.df) <- c('init', 'rep', 'eval2', 'avmA', 'avmB')

metrics.df <- cbind(par.table [metrics.df $ rep, ], metrics.df)

head(metrics.df)

metrics.df %>%
    mutate('Mod2' = Mod(eval2),
           'mean.match' = rowMeans(cbind(avmA, avmB)),
           'scenario' = paste0('mA = ', mA, ', mB = ', mB)) %>%
    filter(g > 0.09 & mA >= mB) %>%
    ggplot(.) +
    geom_point(aes(x = mean.match, y = Mod2)) +
    facet_wrap(~ scenario) +
    #scale_color_viridis() +
    theme_bw()

eval.df <- cbind(as.character(par.table [eval.df $ rep, ]), eval.df)

colnames(eval.df)[4:7] <- c('init', 'rep', 'eval', 'v')

eval.df %>%
    mutate('Reval' = Re(v), 'Imeval' = Im(v)) %>%
    filter(mA == '0.6' & mB == '0.2') %>%
    ggplot(., aes(x = Reval, y = Imeval, color = g)) +
    stat_density_2d() +
        ## facet_grid(mA ~ mB) +
    theme_bw()
