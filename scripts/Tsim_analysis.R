require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(viridis)
require(cowplot)

##

## load data
load('data/Tsims.RData')

source('functions/buildTmat.R')
source('functions/assembleQ.R')

## real columns as real
Tsims.re.trans <-
    aaply(
        as.matrix(Tsims.out[, laply(Tsims.out, is.complex)]), 2,
        function(C) !any(Im(C) > .Machine $ double.eps))


Tsims.re.df <- llply(Tsims.out[, laply(Tsims.out, is.complex)] [, Tsims.re.trans], Re)

Tsims.out[, laply(Tsims.out, is.complex)] [, Tsims.re.trans] <- as.data.frame(Tsims.re.df)

##
T.sims.sumdf <-
    Tsims.out %>%
    mutate('final.tm' = 0.5*(final.traitmatch.A + final.traitmatch.B)) %>%
    group_by(network, g, mA, mB, type) %>%
    summarise_each(funs(mean))

load('data/Tsims_gzero.RData')

## real columns as real
gzero.re.trans <-
    aaply(
        as.matrix(Tsims.gzero[, laply(Tsims.gzero, is.complex)]), 2,
        function(C) !any(Im(C) > .Machine $ double.eps))


Tsims.re.df <- llply(Tsims.gzero[, laply(Tsims.gzero, is.complex)] [, gzero.re.trans], Re)

Tsims.gzero[, laply(Tsims.gzero, is.complex)] [, gzero.re.trans] <- as.data.frame(Tsims.re.df)

gzero.sumdf <-
    Tsims.gzero %>%
    mutate('final.tm' = 0.5*(final.traitmatch.A + final.traitmatch.B)) %>%
    group_by(network, g, mA, type) %>%
    summarise_each(funs(mean))


## FIG gzero_all_sites_trait_match
pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_gzero_all_sites_trait_match_vs_corr.pdf', width = 8, height = 5)

ggplot(gzero.sumdf) +
    geom_line(aes(x = final.mcorr.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.mcorr.lines, y = final.tm, color = mA)) +
    #facet_wrap(~ type) +
    #facet_grid(mA ~ mB) +
    scale_color_viridis('Mutualism\nStrength', option = 'C', direction = -1) +
    theme_bw() +
    ylab('Trait Match') +
    xlab('Landscape Average Correlation')

dev.off(dev.cur())

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_gzero_all_sites_trait_match_vs_sd.pdf',
    width = 8, height = 5)

ggplot(gzero.sumdf) +
    geom_line(aes(x = final.sd.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.sd.lines, y = final.tm, color = mA)) +
    scale_color_viridis('Mutualism\nStrength', option = 'C', direction = -1) +
    theme_bw() +
    ylab('Trait Match') +
    xlab('Landscape Average Standard Deviation')

dev.off(dev.cur())


pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_gzero_all_sites_trait_match_vs_perc.pdf', width = 8, height = 5)

ggplot(gzero.sumdf) +
    geom_line(aes(x = final.perc, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.perc, y = final.tm, color = mA)) +
    scale_color_viridis('Mutualism\nStrength', option = 'C', direction = -1) +
    theme_bw() +
    ylab('Trait Match') +
    xlab('% First Eigenvalue')

dev.off(dev.cur())

## FIG S8

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_sitetype_m07_trait_match_vs_corr.pdf',
    width = 8, height = 8)

ggplot(subset(T.sims.sumdf, mA == '0.7' & mB == '0.7')) +
    geom_line(aes(x = final.mcorr.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.mcorr.lines, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ggtitle('mA = 0.7, mB = 0.7') +
    xlab('Average Landscape Correlation') +
    ylab('Trait Match') +
    theme_bw() 

dev.off(dev.cur())

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_sitetype_m07_trait_match_vs_sd.pdf',
    width = 8, height = 8)

ggplot(subset(T.sims.sumdf, mA == '0.7' & mB == '0.7')) +
    geom_line(aes(x = final.sd.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.sd.lines, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ggtitle('mA = 0.7, mB = 0.7') +
    xlab('Average Landscape Standard Deviation') +
    ylab('Trait Match') +
    theme_bw() 

dev.off(dev.cur())

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_sitetype_m07_trait_match_vs_perc.pdf',
    width = 8, height = 8)

ggplot(subset(T.sims.sumdf, mA == '0.7' & mB == '0.7')) +
    geom_line(aes(x = final.perc, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.perc, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ggtitle('mA = 0.7, mB = 0.7') +
    xlab('% First Eigenvalue') +
    ylab('Trait Match') +
    theme_bw() 

dev.off(dev.cur())

## FIG S8_2

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_sitetype_mA09_mB01_trait_match_vs_corr.pdf',
    width = 8, height = 8)

ggplot(subset(T.sims.sumdf, mA == '0.9' & mB == '0.1')) +
    geom_line(aes(x = final.mcorr.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.mcorr.lines, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ggtitle('mA = 0.7, mB = 0.7') +
    xlab('Average Landscape Correlation') +
    ylab('Trait Match') +
    theme_bw() 

dev.off(dev.cur())

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_sitetype_mA09_mB01_trait_match_vs_sd.pdf',
    width = 8, height = 8)

ggplot(subset(T.sims.sumdf, mA == '0.9' & mB == '0.1')) +
    geom_line(aes(x = final.sd.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.sd.lines, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ggtitle('mA = 0.7, mB = 0.7') +
    xlab('Average Landscape Standard Deviation') +
    ylab('Trait Match') +
    theme_bw() 

dev.off(dev.cur())

pdf('~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_sitetype_mA09_mB01_trait_match_vs_perc.pdf',
    width = 8, height = 8)

ggplot(subset(T.sims.sumdf, mA == '0.9' & mB == '0.1')) +
    geom_line(aes(x = final.perc, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.perc, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ggtitle('mA = 0.7, mB = 0.7') +
    xlab('% First Eigenvalue') +
    ylab('Trait Match') +
    theme_bw() 

dev.off(dev.cur())

### make some matrices

example.networks <- dir('example_sims', pattern = '.txt')

example.graph <-
    alply(1:2, 1,
          function(i)
              as.matrix(read.csv(paste0('example_sims/', example.networks [i]),
                                 header = FALSE, sep = ' ')))

names(example.graph) <- c('Carlo', 'Galetti')

example.sim <- dir('example_sims', pattern = '.csv')

example.pars <-
    aaply(example.sim, 1, function(fil)
    {
        params <- strsplit(fil, '[(?=[:alpha:]) (?=[:punct:])]') [[1]]
        ms <- params[params != ''] [2:4]
        ifelse(ms == '005', as.numeric(ms) * .01, as.numeric(ms) * .1)
    })

colnames(example.pars) <- c('mA', 'mB', 'g')

example.pars <- example.pars [c(1:8) * 2, ]

example.Z <-
    alply(example.sim, 1, function(fil)
        as.matrix(read.csv(paste0('example_sims/', fil), row.names = 1)))

ex.Z.A <- example.Z [(c(1:8) * 2) - 1]
ex.Z.B <- example.Z [(c(1:8) * 2)]

alpha <- 0.2
phi <- 0.5 # with 0.1 sd

## carlo

carloT <-
    alply(1:4, 1, function(i)
        buildTmat(example.graph [[1]],
                  g = example.pars [i, 'g'],
                  m.A = example.pars [i, 'mA'],
                  m.B = example.pars [i, 'mB'],
                  theta.A = ex.Z.A [[i]] ['theta_A', ],
                  theta.B = ex.Z.B [[i]] ['theta_B', ],
                  value.A = ex.Z.A [[i]] ['z_A_final', ],
                  value.B = ex.Z.B [[i]] ['z_B_final', ],
                  output.matrix = TRUE))

color2D.matplot(carloT [[1]])

carloT <-
    laply(carloT, function(mat)
    {
        rownames(mat) <- colnames(mat) <- NULL
        mat
    })

carloT <- aperm(carloT, c(2, 3, 1))

dimnames(carloT) <- list(1:72, 1:72, 1:4)

names(dimnames(carloT)) <- c('row', 'col', 'slice')

carlo.df <- adply(carloT, 1:3)

colnames(carlo.df) <- c('row', 'col', 'slice', 'value')

carlo.df <- cbind(carlo.df, example.pars [as.numeric(carlo.df $ slice), ])

carlo.df $ row <- factor(as.character(carlo.df $ row), levels = as.character(72:1))

carlo.df $ mA <- factor(paste('mA =', carlo.df $ mA),
                        levels = paste('mA =', unique(carlo.df $ mA)))

carlo.df $ mB <- factor(paste('mB =', carlo.df $ mA),
                        levels = paste('mB =', unique(carlo.df $ mA)))


## all matrices
ggplot(carlo.df) +
    geom_tile(aes(y = row, x = col, fill = value)) +
    facet_wrap(~ slice) +
    scale_fill_viridis(option = 'B', direction = -1, na.value = 'white') +
    theme_bw()

## submatrices no geneflow

pdf(file = '~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_Carlo_no_geneflow_m01_09.pdf',
    width = 15, height = 7)

ggplot(subset(carlo.df, row %in% as.character(1:36) & col %in% as.character(1:36) & g == 0)) +
    geom_tile(aes(y = row, x = col, fill = value)) +
    facet_wrap(~ mA, nrow = 1) +
    scale_fill_viridis('Interaction\nStrength',
                       option = 'B', direction = -1, na.value = 'white') +
    theme_bw() +
    scale_y_discrete(breaks = NULL) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL) +
    ylab(NULL)
    
dev.off(dev.cur())

## fullmatrices, geneflow effect

pdf(file = '~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_Carlo_g005_03_mA07_mB07.pdf',
    width = 15, height = 7)

ggplot(subset(carlo.df, g != 0)) +
    geom_tile(aes(y = row, x = col, fill = value)) +
    facet_wrap(~ g, nrow = 1) +
    scale_fill_viridis('Interaction\nStrength',
                       option = 'B', direction = -1, na.value = 'white') +
    theme_bw() +
    scale_y_discrete(breaks = NULL) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL) +
    ylab(NULL)

dev.off(dev.cur())

## galetti

galettiT <-
    alply(5:8, 1, function(i)
        buildTmat(example.graph [[2]],
                  g = example.pars [i, 'g'],
                  m.A = example.pars [i, 'mA'],
                  m.B = example.pars [i, 'mB'],
                  theta.A = ex.Z.A [[i]] ['theta_A', ],
                  theta.B = ex.Z.B [[i]] ['theta_B', ],
                  value.A = ex.Z.A [[i]] ['z_A_final', ],
                  value.B = ex.Z.B [[i]] ['z_B_final', ],
                  output.matrix = TRUE))

color2D.matplot(galettiT [[1]])

galettiT <-
    laply(galettiT, function(mat)
    {
        rownames(mat) <- colnames(mat) <- NULL
        mat
    })

galettiT <- aperm(galettiT, c(2, 3, 1))

dimnames(galettiT) <- list(1:128, 1:128, 1:4)

names(dimnames(galettiT)) <- c('row', 'col', 'slice')

galetti.df <- adply(galettiT, 1:3)

colnames(galetti.df) <- c('row', 'col', 'slice', 'value')

galetti.df <- cbind(galetti.df, example.pars [as.numeric(galetti.df $ slice), ])

galetti.df $ row <- factor(as.character(galetti.df $ row), levels = as.character(128:1))

galetti.df $ mA <- factor(paste('mA =', galetti.df $ mA),
                        levels = paste('mA =', unique(galetti.df $ mA)))

galetti.df $ mB <- factor(paste('mB =', galetti.df $ mA),
                        levels = paste('mB =', unique(galetti.df $ mA)))


## all matrices
ggplot(galetti.df) +
    geom_tile(aes(y = row, x = col, fill = value)) +
    facet_wrap(~ slice) +
    scale_fill_viridis(option = 'B', direction = -1, na.value = 'white') +
    theme_bw()

## submatrices no geneflow

pdf(file = '~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_Galetti_no_geneflow_m01_09.pdf',
    width = 15, height = 7)

ggplot(subset(galetti.df, row %in% as.character(1:64) &
                          col %in% as.character(1:64) & g == 0)) +
    geom_tile(aes(y = row, x = col, fill = value)) +
    facet_wrap(~ mA, nrow = 1) +
    scale_fill_viridis('Interaction\nStrength',
                       option = 'B', direction = -1, na.value = 'white') +
    theme_bw() +
    scale_y_discrete(breaks = NULL) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL) +
    ylab(NULL)
    
dev.off(dev.cur())

## fullmatrices, geneflow effect

pdf(file = '~/Dropbox/spatial_coevo_mutnet/results/figs/figSI_Galetti_g005_03_mA07_mB07.pdf',
    width = 15, height = 7)

ggplot(subset(galetti.df, g != 0)) +
    geom_tile(aes(y = row, x = col, fill = value)) +
    facet_wrap(~ g, nrow = 1) +
    scale_fill_viridis('Interaction\nStrength',
                       option = 'B', direction = -1, na.value = 'white') +
    theme_bw() +
    scale_y_discrete(breaks = NULL) +
    scale_x_discrete(breaks = NULL) +
    xlab(NULL) +
    ylab(NULL)

dev.off(dev.cur())
