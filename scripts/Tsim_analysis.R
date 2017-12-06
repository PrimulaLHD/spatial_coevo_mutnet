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


## FIG S5
ggplot(gzero.sumdf) +
    geom_line(aes(x = final.mcorr.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.mcorr.lines, y = final.tm, color = mA)) +
    facet_wrap(~ type) +
    #facet_grid(mA ~ mB) +
    scale_color_viridis(option = 'C', direction = -1) +
    theme_bw() +
    #theme(legend.position = 'none') +
    scale_y_log10()

plot_grid(
    ggplot(subset(T.sims.sumdf, mA == '0.7' & mB == '0.7')) +
    geom_line(aes(x = final.mcorr.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.mcorr.lines, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ##scale_y_log10() +
    scale_x_log10() +
    theme_bw(), 

    ggplot(subset(T.sims.sumdf, mA == '0.9' & mB == '0.1')) +
    geom_line(aes(x = final.mcorr.lines, y = final.tm, group = network),
              alpha = 0.5) +
    geom_point(aes(x = final.mcorr.lines, y = final.tm, color = g), alpha = 0.5) +    
    facet_wrap(~ type, nrow = 3) +
    ##facet_grid(mA ~ mB) +
    scale_color_viridis(name = 'Gene\nFlow', option = 'A', direction = -1) +
    ##scale_y_log10() +
    scale_x_log10() +
    theme_bw(), 

    ncol = 2)
