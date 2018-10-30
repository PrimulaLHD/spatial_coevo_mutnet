require(doMC)
registerDoMC(cores = 3)
require(plotrix)
require(plyr)
require(dplyr)
require(tidyr)
require(magrittr)
require(ggplot2)
require(viridis)
require(cowplot)
require(RColorBrewer)

##

## load data
load('sim_analysis/Tsims.RData')

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

load('sim_analysis/Tsims_gzero.RData')

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


## Spectral
grad <- colorRampPalette(brewer.pal(11, 'Spectral'))

gzero.sumdf $ type <- factor(gzero.sumdf $ type)

gzero.sumdf $ type <- revalue(gzero.sumdf $ type,
                              c('AE' = 'Ants - Nectary-bearing Plants',
                                'AM' = 'Ants - Myrmecophytes',
                                'CC' = 'Marine Cleaning',
                                'FA' = 'Anemones - Anemonefishes',
                                'FP' = 'Seed Dispersal',
                                'PP' = 'Pollination'))
  

## FIG gzero_all_sites_trait_match
ggsave(
    '~/spatial_coevo_mutnet/results/figs/figSI_gzero_all_sites_trait_match_vs_corr.pdf',
    width = 12, height = 8,
    plot = 
        ggplot(gzero.sumdf) +
        geom_line(aes(x = final.mcorr.lines, y = final.traitmatch.A, group = network),
                  color = 'grey') +
        geom_point(aes(x = final.mcorr.lines, y = final.traitmatch.A, color = mA)) +
        facet_wrap(~ type) +
        ## facet_grid(mA ~ mB) +
        scale_color_viridis('Mutualistic\nSelection', option = 'D', direction = 1) +
        ## theme_bw() +
        ylab('Mean Trait Matching') +
        xlab('Mean Adaptive Landscape Correlation')
)

## FIG S8

T.sims.sumdf $ type <- factor(T.sims.sumdf $ type)

T.sims.sumdf $ type <- revalue(T.sims.sumdf $ type,
                              c('AE' = 'Ants - Nectary-bearing Plants',
                                'AM' = 'Ants - Myrmecophytes',
                                'CC' = 'Marine Cleaning',
                                'FA' = 'Anemones - Anemonefishes',
                                'FP' = 'Seed Dispersal',
                                'PP' = 'Pollination'))


ggsave(
    '~/spatial_coevo_mutnet/results/figs/figSI_sitetype_m07_trait_match_vs_corr.pdf',
    width = 12, height = 8,
    plot = 
        ggplot(subset(T.sims.sumdf, mA == '0.7' & mB == '0.7')) +
        geom_line(aes(x = final.mcorr.lines, y = final.traitmatch.A, group = network),
                  color = 'grey') +
        geom_point(aes(x = final.mcorr.lines, y = final.traitmatch.A, color = g)) +    
        facet_wrap(~ type, nrow = 2) +
        ##facet_grid(mA ~ mB) +
        scale_color_viridis(name = 'Gene\nFlow', option = 'D', direction = 1) +
        ggtitle('mA = 0.7, mB = 0.7') +
        ylab('Mean Trait Matching') +
        xlab('Mean Adaptive Landscape Correlation')
)

## FIG S8_2

ggsave(
   '~/spatial_coevo_mutnet/results/figs/figSI_sitetype_mA09_mB01_trait_match_vs_corr.pdf',
   width = 12, height = 8, 
   plot = 
        ggplot(subset(T.sims.sumdf, mA == '0.9' & mB == '0.1')) +
        geom_line(aes(x = final.mcorr.lines, y = final.traitmatch.A, group = network),
                  color = 'grey') +
        geom_point(aes(x = final.mcorr.lines, y = final.traitmatch.A, color = g)) +    
        facet_wrap(~ type, nrow = 2) +
        ##facet_grid(mA ~ mB) +
        scale_color_viridis(name = 'Gene\nFlow', option = 'D', direction = 1) +
        ggtitle('mA = 0.9, mB = 0.1') +
        ylab('Mean Trait Matching') +
        xlab('Mean Adaptive Landscape Correlation')
)

maxvalue = which.max(subset(T.sims.sumdf, mA == '0.7' &
                                          mB == '0.7' &
                                          type == 'Seed Dispersal') $ final.mcorr.lines)

subset(T.sims.sumdf, mA == '0.7' & mB == '0.7' & type == 'Seed Dispersal') [maxvalue, ]

### make some matrices

## FP_Sorensen_1981

sorensenFiles =
    dir('../sim_results', include.dirs = TRUE, full.names = TRUE,
        recursive = TRUE, pattern = 'FP_Sorensen_1981')

sorensenFiles = sorensenFiles[grepl('.csv', sorensenFiles)]

head(sorensenFiles)

sorensenSimOne = sorensenFiles[grepl('sim12\\.', sorensenFiles)]

m07Index =
    sorensenSimOne %>%
    strsplit('/') %>%
    llply(function(e) e[5]) %>%
    grepl(pattern = 'mA07_mB07', .)

g0Index = 
    sorensenSimOne %>%
    strsplit('/') %>%
    llply(function(e) e[5]) %>%
    grepl(pattern = 'g0_', .)

g005Index = 
    sorensenSimOne %>%
    strsplit('/') %>%
    llply(function(e) e[5]) %>%
    grepl(pattern = 'g005_', .)

g03Index = 
    sorensenSimOne %>%
    strsplit('/') %>%
    llply(function(e) e[5]) %>%
    grepl(pattern = 'g03_', .)

m01Index =
    sorensenSimOne %>%
    strsplit('/') %>%
    llply(function(e) e[5]) %>%
    grepl(pattern = 'mA01_mB01', .)

m09Index =
    sorensenSimOne %>%
    strsplit('/') %>%
    llply(function(e) e[5]) %>%
    grepl(pattern = 'mA09_mB09', .)



sorensenNetwork = read.table('data/empirical_networks/binary/FP_Sorensen_1981_bin.txt')
sorensenNetwork = as.matrix(sorensenNetwork)
dimnames(sorensenNetwork) = NULL

buildDf =
    function(index) {
        fs = sorensenSimOne[index] 

        contents = alply(fs, 1,function(f) {
            df = f %>% read.csv
            as.matrix(df [, -1])
        })

        params = strsplit(fs[1], '[(?=[:alpha:]) (?=[:punct:])]') [[1]]
        ms = params[params != '']
        ms = ms[length(ms) - 3:1]
        ms = ifelse(ms == '005', as.numeric(ms) * .01, as.numeric(ms) * .1)
        names(ms) = c('mA', 'mB', 'g')

        intMat = 
            buildTmat(
                sorensenNetwork,
                m.A = ms['mA'], m.B = ms['mB'], g = ms['g'],
                theta.A = contents[[1]] [1, ],
                theta.B = contents[[2]] [1, ],
                value.A = contents[[1]] [3, ],
                value.B = contents[[2]] [3, ],
                output.matrix = TRUE
            )
        dimnames(intMat) = NULL
        dimnames(intMat) = list('row' = 1:nrow(intMat), 'col' = 1:ncol(intMat))
        # names(dimnames(intMat)) = c('row', 'col')
        matDf = adply(intMat, 1:2)
        colnames(matDf) [3] = 'value'
        msMat = rep(ms, each = nrow(matDf)) %>% matrix(ncol = 3)
        colnames(msMat) = names(ms)
        matDf = cbind(matDf, msMat)

        matDf
    }

sorensenDf =
    ldply(list(m01Index & g0Index, 
               m07Index & g0Index,
               m09Index & g0Index,
               m07Index & g005Index,
               m07Index & g03Index), buildDf)

sorensenDf $ row = factor(as.character(sorensenDf $ row), levels = as.character(26:1))

sorensenDf $ mA %>% head

## submatrices no geneflow

ggsave(
    '~/spatial_coevo_mutnet/results/figs/figSI_Sorensen_no_geneflow_m01_09.pdf',
    width = 10, height = 5,
    plot = 
        ggplot(subset(sorensenDf,
                      row %in% as.character(1:13) &
                      col %in% as.character(1:13) &
                      g == 0 &
                      mA %in% c(0.1, 0.9) &
                      abs(value) > .Machine$double.eps)) +
        geom_tile(aes(y = row, x = col, fill = value)) +
        facet_wrap(~ mA, nrow = 1, labeller = labeller(mA = label_both)) +
        scale_fill_gradientn('Evolutionary\nEffect', colours = rev(grad(11))) +
        ##    scale_fill_viridis('Evolutionary\nEffect', option = 'B', direction = -1) +
        ##    theme_bw() +
        scale_y_discrete(breaks = NULL) +
        scale_x_discrete(breaks = NULL) +
        xlab(NULL) +
        ylab(NULL)
)
    
## full matrices, geneflow effect

ggsave(
    '~/spatial_coevo_mutnet/results/figs/figSI_Sorensen_g005_03_mA07_mB07.pdf',
    width = 15, height = 5,
    plot =
        ggplot(subset(sorensenDf,
                      mA == '0.7' & 
                      abs(value) > .Machine$double.eps)) +
        geom_tile(aes(y = row, x = col, fill = value)) +
        facet_wrap(~ g, nrow = 1, labeller = labeller(g = label_both)) +
        scale_fill_gradientn('Evolutionary\nEffect', colours = rev(grad(11))) +
        ##    scale_fill_viridis('Evolutionary\nEffect', option = 'A', direction = -1) +
        ##    theme_bw() +
        scale_y_discrete(breaks = NULL) +
        scale_x_discrete(breaks = NULL) +
        xlab(NULL) +
        ylab(NULL)
)
