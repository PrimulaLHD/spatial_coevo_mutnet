require(doMC)
registerDoMC(cores = 16)
require(plotrix)
require(plyr)
require(dplyr)
require(magrittr)

source('../functions/assembleQ.R')
source('../functions/buildTmat.R')
source('../functions/convMutNet_2.R')

## networks
binary.files <- dir('../data/empirical_networks/binary/')
binary.networks <- list()

for(i in 1:length(binary.files))
    binary.networks [[i]] <-
        as.matrix(read.table(paste0('../data/empirical_networks/binary/',
                                    binary.files [i])))

names(binary.networks) <- gsub('.txt', '', binary.files)

## unzip
filenames <- dir('../../sim_results/')
filenames <- filenames [grepl('all_networks', filenames)]

## fixed parameters
alpha <- 0.2
phi <- 0.5 # with 0.1 sd


Tsims.out <-
    adply(1:length(filenames), 1, function(i) ## run over zip files
    {

        system(paste0('unzip /home/guilherme/sim_results/', filenames [i]),
               ignore.stdout = TRUE)

        print('unzipped files')
        
        dirname <- gsub('.zip', '', filenames[i])

        ## get other values from filename
        params <- strsplit(filenames[i], '[(?=[:alpha:]) (?=[:punct:])]') [[1]]

        ms <- params[params != ''] [1:2]

        ms <- as.numeric(ms) * .1

        mA <- ms[1]
        mB <- ms[2]

        nets.output <-
            adply(1:length(binary.networks), 1, function(j)
            { ## run over networks
                
                print(names(binary.networks) [j])
        
                n.sps <- dim(binary.networks[[j]])

                path2csv <- paste(dirname, names(binary.networks)[j], sep = '/')
                csv.files <- dir(path2csv)
                ## need g and number and site

                csv.numbers <- strsplit(csv.files, '[(?=[:alpha:]) (?=[:punct:])]')

                csv.info <- laply(csv.numbers, function(S) S[S != ''])

                csv.info <- csv.info [, (-1:0) + ncol(csv.info)]
        
                csv.site <- c()
                csv.site[grepl('siteA', csv.files)] <- 'A'
                csv.site[grepl('siteB', csv.files)] <- 'B'

                csv.info <- data.frame('g' = as.numeric(csv.info [, 1]) *
                                           10^(-(nchar(csv.info [, 1]) - 1)),
                                       'site' = csv.site,
                                       'sim' = as.numeric(csv.info [, 2]))

                csv.info $ 'mA' <- rep(mA, times = nrow(csv.info))
                csv.info $ 'mB' <- rep(mB, times = nrow(csv.info))

                file.assoc <-
                    matrix(order(csv.info $ g, csv.info $ sim), ncol = 2, byrow = T)
        
                csv.info <- csv.info [file.assoc [, 1], -2]

                csv.info $ network <-
                    rep(gsub('_bin', '', names(binary.networks)[j]), nrow(csv.info))

                csv.info $ type <-
                    rep(strsplit(names(binary.networks)[j], '_') [[1]] [1], nrow(csv.info))

                sim.output <-
                    adply(1:nrow(file.assoc), 1, function(k)
                    {
                        siteAsim <- read.csv(paste(path2csv,
                                                   csv.files [file.assoc [k, 1]],
                                                   sep = '/'),
                                             row.names = 1)
                        siteBsim <- read.csv(paste(path2csv,
                                                   csv.files [file.assoc [k, 2]],
                                                   sep = '/'),
                                             row.names = 1)
        
                        siteAsim <- data.frame(t(siteAsim))
                        siteBsim <- data.frame(t(siteBsim))
        
                        T.initial <- buildTmat(binary.networks[[j]],
                                               csv.info $ g [k], phi, alpha,
                                               csv.info $ mA [k], csv.info $ mB [k],
                                               siteAsim $ theta_A, siteBsim $ theta_B,
                                               siteAsim $ z_A_init, siteBsim $ z_B_init)
        
                        T.final <- buildTmat(binary.networks[[j]],
                                             csv.info $ g [k], phi, alpha,
                                             csv.info $ mA [k], csv.info $ mB [k],
                                             siteAsim $ theta_A, siteBsim $ theta_B,
                                             siteAsim $ z_A_final, siteBsim $ z_B_final)
        
                        names(T.initial) <- paste('initial', names(T.initial), sep = '.')
                        names(T.final) <- paste('final', names(T.final), sep = '.')
        
                        c(T.initial, T.final)
                    }, .parallel = TRUE)

                cbind(csv.info, sim.output)
            })
        

        ## in the end (of i iteration)
        system(paste('rm -rfv', dirname), ignore.stdout = TRUE)
        system('rm -rfv __MACOSX', ignore.stdout = TRUE)

        print('removed files')
        
        nets.output
    })

save(Tsims.out, file = 'Tsims.RData')
