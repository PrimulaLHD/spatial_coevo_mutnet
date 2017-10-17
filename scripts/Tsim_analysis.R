require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)
require(ggplot2)
require(viridis)

## load data
load('data/Tsims.RData')


## real columns as real
Tsims.re.trans <-
    aaply(
        as.matrix(Tsims.out[, laply(Tsims.out, is.complex)]), 2,
        function(C) !any(Im(C) > .Machine $ double.eps))


Tsims.re.df <- llply(Tsims.out[, laply(Tsims.out, is.complex)] [, Tsims.re.trans], Re)

Tsims.out[, laply(Tsims.out, is.complex)] [, Tsims.re.trans] <- as.data.frame(Tsims.re.df)


X11()
##
T.sims.sumdf <-
    Tsims.out %>%
    mutate('final.tm' = 0.5*(final.traitmatch.A + final.traitmatch.B)) %>%
    group_by(network, g, mA, mB, type) %>%
    summarise_each(funs(mean))

ggplot(T.sims.sumdf) +
    geom_point(aes(x = final.perc, y = final.tm,
                   color = g, shape = type),
               alpha = 0.5) +
    facet_wrap(~ type) +
    #facet_grid(mA ~ mB) +
    scale_color_viridis(option = 'C', direction = -1) +
    theme_bw() +
    scale_x_log10()
               
    
