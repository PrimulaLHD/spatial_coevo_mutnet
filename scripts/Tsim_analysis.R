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

