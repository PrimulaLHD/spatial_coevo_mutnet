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

Ols.evar <- aaply(Ols.eval, c(1, 2), function(ev) var(Mod(ev)))

Ols.det <- aaply(Olesen.Tmat, c(1, 2), function(T) det(T))

names(dimnames(Ols.eval)) <- c('init', 'rep', 'eval')

names(dimnames(Ols.evar)) <- c('init', 'rep')

names(dimnames(Ols.det)) <- c('init', 'rep')

## have some non zero im part

Ols.match <- 
    laply(Olesen.spec, function(L1)
        laply(L1 [[1]] [1:96], function(L2) L2 $ av.match))

match.df <- adply(Ols.match, 1:2)
eval2.df <- adply(Ols.eval [, , 2], 1:2)
evar.df <-  adply(Ols.evar, 1:2)
det.df <- adply(Ols.det, 1:2)

metrics.df <- cbind(eval2.df, evar.df [, 3], det.df [, 3], match.df [, 2:3])

colnames(metrics.df) <- c('init', 'rep', 'eval2', 'evar', 'det', 'avmA', 'avmB')

metrics.df <- cbind(par.table [metrics.df $ rep, ], metrics.df)

head(metrics.df)

metrics.df $ g <- factor(metrics.df $ g)
metrics.df $ mA <- factor(metrics.df $ mA)
metrics.df $ mB <- factor(metrics.df $ mB) 

var.avm.Ols <- 
    metrics.df %>%
    mutate('Mod2' = Mod(eval2),
           'mean.match' = rowMeans(cbind(avmA, avmB)),
           'scenario' = paste0('mA = ', mA, ', mB = ', mB)) %>%
    ggplot(.) +
    geom_point(aes(x = avmA, y = evar, color = mA, shape = mB)) +
    facet_wrap(~ g) +
    scale_color_viridis(option = 'C', direction = -1, discrete = TRUE) +
    theme_bw()

ggsave('var_avm_Ols.pdf', var.avm.Ols, width = 10, height = 10)

eval.df <- cbind(as.character(par.table [eval.df $ rep, ]), eval.df)

colnames(eval.df)[4:7] <- c('init', 'rep', 'eval', 'v')

eval.df %>%
    mutate('Reval' = Re(v), 'Imeval' = Im(v)) %>%
    filter(mA == '0.6' & mB == '0.2') %>%
    ggplot(., aes(x = Reval, y = Imeval, color = g)) +
    stat_density_2d() +
        ## facet_grid(mA ~ mB) +
    theme_bw()
