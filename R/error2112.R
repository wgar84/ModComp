require (plotrix)
require (ape)
require (Morphometrics)
require (geiger)
require (bbmle)
require (nlme)
require (plyr)
require (shapes)
require (geomorph)
require (numDeriv)
require (expm)
require (xtable)
require (ggplot2)
require (reshape2)
require (RColorBrewer)
require (grid)
require (gridExtra)
require (mvtnorm)

require (doMC)
registerDoMC (cores = 10)

attach ('../../Databases/ED.RData')
attach ('../../Databases/Sym.RData')
attach ('../../Databases/Aux.RData')
attach ('../../Databases/LOG.ED.RData')
attach ('../../Databases/NCS.ED.RData')
attach ('../../Databases/OneDef.RData')
attach ('../../Databases/Tree.RData')

source ('../../Func/ModTest.R')

modsim.Data <- list ()

modsim.Plots <- list()

modsim.Data $ cor.df <- melt (list ('Sym' = aaply (modcomp.Data $ Sym, 1, cov2cor),
                                    'ED' = aaply (modcomp.Data $ ED, 1, cov2cor),
                                    'Def' = aaply (modcomp.Data $ Def, 1, cov2cor)))

colnames (modsim.Data $ cor.df) <- c('otu', 't1', 't2', 'value', 'type')

modsim.Data $ cor.df <- subset (modsim.Data $ cor.df, t1 != t2)

modsim.Data $ cor.df <- subset (modsim.Data $ cor.df, t1 != 'logCS' & t2 != 'logCS')

modsim.Data $ var.df <- melt (list ('Sym' = aaply (modcomp.Data $ Sym, 1, diag),
                                    'ED' = aaply (modcomp.Data $ ED, 1, diag),
                                    'Def' = aaply (modcomp.Data $ Def, 1, diag)))

modsim.Data $ cor.df $ t1.num <- as.numeric (modsim.Data $ cor.df $ t1)
modsim.Data $ cor.df $ t2.num <- as.numeric (modsim.Data $ cor.df $ t2)

colnames (modsim.Data $ var.df) <- c('otu', 'trait', 'value', 'type')

modsim.Data $ var.df <- subset (modsim.Data $ var.df, trait != 'logCS')

### SIM PREP

source ('sim.func.R')
source ('error.R')

modsim.Data $ Sym.var <- subset (modsim.Data $ var.df, type == 'Sym') $ value
modsim.Data $ ED.var <- subset (modsim.Data $ var.df, type == 'ED') $ value
modsim.Data $ Def.var <- subset (modsim.Data $ var.df, type == 'Def') $ value

modsim.Data $ hyp.tensor.Sym <- outer (modcomp.Data $ Hyp $ Sym [-1, 1:6],
                                       modcomp.Data $ Hyp $ Sym [-1, 1:6])

modsim.Data $ hyp.tensor.EDef <- outer (modcomp.Data $ Hyp $ EDef [, 1:6],
                                        modcomp.Data $ Hyp $ EDef [, 1:6])

modsim.Data $ Sym.within <- c()
modsim.Data $ Sym.between <- c()

modsim.Data $ ED.within <- c()
modsim.Data $ ED.between <- c()

modsim.Data $ Def.within <- c()
modsim.Data $ Def.between <- c()

modsim.Data $ Sym.SR.within <- c()
modsim.Data $ Sym.SR.between <- c()

modsim.Data $ ED.SR.within <- c()
modsim.Data $ ED.SR.between <- c()

modsim.Data $ Def.SR.within <- c()
modsim.Data $ Def.SR.between <- c()


for (i in 1:21)
  for (j in 1:5)
    for (k in j:6)
      {
        tmp.cor <- cov2cor (modcomp.Data $ Sym [i, -1, -1]) [which (
          modsim.Data $ hyp.tensor.Sym [, j, , k] > 0, arr.ind = TRUE)]
        if (j == k)
          {
            n.diag <- sum (diag (modsim.Data $ hyp.tensor.Sym [, j, , k]))
            modsim.Data $ Sym.within [length (modsim.Data $ Sym.within) + 1] <-
              (sum (tmp.cor) - n.diag) / (length (tmp.cor) - n.diag)
          }
        else
          modsim.Data $ Sym.between [length (modsim.Data $ Sym.between) + 1] <-
            mean (tmp.cor)

        tmp.cor <- cov2cor (modcomp.Data $ ED [i, , ]) [which (
          modsim.Data $ hyp.tensor.EDef [, j, , k] > 0, arr.ind = TRUE)]
        if (j == k)
          {
            n.diag <- sum (diag (modsim.Data $ hyp.tensor.EDef [, j, , k]))
            modsim.Data $ ED.within [length (modsim.Data $ ED.within) + 1] <-
              (sum (tmp.cor) - n.diag) / (length (tmp.cor) - n.diag)
          }
        else
          modsim.Data $ ED.between [length (modsim.Data $ ED.between) + 1] <-
            mean (tmp.cor)
        
        tmp.cor <- cov2cor (modcomp.Data $ Def [i, -1, -1]) [which (
          modsim.Data $ hyp.tensor.EDef [, j, , k] > 0, arr.ind = TRUE)]
        if (j == k)
          {
           n.diag <- sum (diag (modsim.Data $ hyp.tensor.EDef [, j, , k]))
           modsim.Data $ Def.within [length (modsim.Data $ Def.within) + 1] <-
             (sum (tmp.cor) - n.diag) / (length (tmp.cor) - n.diag)
         }
        else
          modsim.Data $ Def.between [length (modsim.Data $ Def.between) + 1] <-
            mean (tmp.cor)
        
          
        tmp.cor <- cov2cor (modcomp.Data $ Sym.SR [i, , ]) [which (
          modsim.Data $ hyp.tensor.Sym [, j, , k] > 0, arr.ind = TRUE)]
        if (j == k)
          {
            n.diag <- sum (diag (modsim.Data $ hyp.tensor.Sym [, j, , k]))
            modsim.Data $ Sym.SR.within [length (modsim.Data $ Sym.SR.within) + 1] <-
              (sum (tmp.cor) - n.diag) / (length (tmp.cor) - n.diag)
          }
        else
          modsim.Data $ Sym.SR.between [length (modsim.Data $ Sym.SR.between) + 1] <-
            mean (tmp.cor)

        tmp.cor <- cov2cor (modcomp.Data $ ED.SR [i, , ]) [which (
          modsim.Data $ hyp.tensor.EDef [, j, , k] > 0, arr.ind = TRUE)]
        if (j == k)
         {
           n.diag <- sum (diag (modsim.Data $ hyp.tensor.EDef [, j, , k]))
           modsim.Data $ ED.SR.within [length (modsim.Data $ ED.SR.within) + 1] <-
             (sum (tmp.cor) - n.diag) / (length (tmp.cor) - n.diag)
         }
        else
          modsim.Data $ ED.SR.between [length (modsim.Data $ ED.SR.between) + 1] <-
            mean (tmp.cor)
        
        tmp.cor <- cov2cor (modcomp.Data $ Def.SR [i, , ]) [which (
         modsim.Data $ hyp.tensor.EDef [, j, , k] > 0, arr.ind = TRUE)]
        if (j == k)
          {
            n.diag <- sum (diag (modsim.Data $ hyp.tensor.EDef [, j, , k]))
            modsim.Data $ Def.SR.within [length (modsim.Data $ Def.SR.within) + 1] <-
              (sum (tmp.cor) - n.diag) / (length (tmp.cor) - n.diag)
          }
        else
          modsim.Data $ Def.SR.between [length (modsim.Data $ Def.SR.between) + 1] <-
           mean (tmp.cor)
      }

modsim.Data $ cor.wb.df <-
  with (modsim.Data,
        melt (list ('Def.within' = Def.within, 'Def.between' = Def.between,
                    'Def.SR.within' = Def.SR.within, 'Def.SR.between' = Def.SR.between,
                    'ED.within' = ED.within, 'ED.between' = ED.between,
                    'ED.SR.within' = ED.SR.within, 'ED.SR.between' = ED.SR.between,
                    'Sym.within' = Sym.within, 'Sym.between' = Sym.between,
                    'Sym.SR.within' = Sym.SR.within, 'Sym.SR.between' = Sym.SR.between)))

modsim.Data $ cor.wb.df $ size <-
  c('Retained', 'Removed') [grepl('.SR', modsim.Data $ cor.wb.df $ L1) + 1]

modsim.Data $ cor.wb.df $ type <-
  laply (strsplit(modsim.Data $ cor.wb.df $ L1, '\\.'), function (L) L [1])

modsim.Data $ cor.wb.df $ wb <-
  laply (strsplit(modsim.Data $ cor.wb.df $ L1, '\\.'), function (L) L [length (L)])

modsim.Data $ cor.wb.df <- modsim.Data $ cor.wb.df [, -2]

head (modsim.Data $ cor.wb.df)

modsim.Data $ cor.wb.df $ wb <- factor (modsim.Data $ cor.wb.df $ wb)
levels (modsim.Data $ cor.wb.df $ wb) <- c('Between', 'Within')

modsim.Data $ cor.wb.df $ wb <- relevel (modsim.Data $ cor.wb.df $ wb, 'Within')

modsim.Data $ cor.wb.df $ type <- factor (modsim.Data $ cor.wb.df $ type)
levels (modsim.Data $ cor.wb.df $ type) <-
  c('Local Shape Variables', 'Interlandmark Distances', 'Procrustes Residuals')

#modsim.Data $ cor.wb.df $ type <- relevel (modsim.Data $ cor.wb.df $ type,
#                                           'Procrustes Residuals')

modsim.Data $ cor.wb.df $ size <- factor (modsim.Data $ cor.wb.df $ size)
levels (modsim.Data $ cor.wb.df $ size) <- c('Size Removed', 'Size Retained')

modsim.Data $ cor.wb.df $ size <- relevel (modsim.Data $ cor.wb.df $ size,
                                           'Size Retained')

modsim.Plots $ cor.dist <- 
  ggplot (modsim.Data $ cor.wb.df) +
  geom_boxplot (aes (y = value, x = wb), outlier.shape = '+') +
  facet_grid(size ~ type) + theme_bw() +
  xlab('Correlation Type') + ylab('Correlation Value')
  
### SIMMMM

names (modsim.Data)

## modsim.Data <-
##   within (modsim.Data,
##           {
##             Def.12 <-
##               adply (1:10000, 1, function (i)
##                      {
##                        print(i)
##                        error2112(
##                          within = Def.within,
##                          between = Def.between,
##                          var = Def.var,
##                          n.char = 40,
##                          sample.sizes = seq (20, 100, 20),
##                          it = 1000)
##                      }, .parallel = TRUE)

##             Def.SR.12 <-
##               adply (1:10000, 1, function (i)
##                      {
##                        print (i)
##                        error2112(
##                          within = Def.SR.within,
##                          between = Def.SR.between,
##                          var = Def.var, 
##                          n.char = 40, sample.sizes = seq (20, 100, 20),
##                          samples.per.ss = 10, 
##                          it = 1000)
##                      }, .parallel = TRUE)
            
##             ED.12 <-
##               adply (1:10000, 1, function (i)
##                      {
##                        print (i)
##                        error2112(
##                          within = ED.within,
##                          between = ED.between,
##                          var = ED.var,
##                          n.char = 40,
##                          sample.sizes = seq (20, 100, 20),
##                          samples.per.ss = 10, 
##                          it = 1000)
##                      }, .parallel = TRUE)
            
##             ED.SR.12 <-
##               adply (1:10000, 1, function (i)
##                      {
##                        print (i)
##                        error2112(
##                          within = ED.SR.within,
##                          between = ED.SR.between,
##                          var = ED.var, 
##                          n.char = 40,
##                          sample.sizes = seq (20, 100, 20),
##                          samples.per.ss = 10, 
##                          it = 1000)
##                      }, .parallel = TRUE)
            
##             Sym.12 <-
##               adply (1:10000, 1, function (i)
##                      {
##                        print (i)
##                        error2112(
##                          within = Sym.within,
##                          between = Sym.between,
##                          var = Sym.var,
##                          n.char = 40,
##                          sample.sizes = seq (20, 100, 20),
##                          samples.per.ss = 10, 
##                          it = 1000)
##                      }, .parallel = TRUE)
            
##             Sym.SR.12 <-
##               adply (1:10000, 1, function (i)
##                      {
##                        print (i)
##                        error2112(
##                          within = Sym.SR.within,
##                          between = Sym.SR.between,
##                          var = Sym.var, 
##                          n.char = 40,
##                          sample.sizes = seq (20, 100, 20),
##                          samples.per.ss = 10, 
##                          it = 1000
##                          )
##                      }, .parallel = TRUE)

##           })

## save (modsim.Data, file = 'modsim12.RData')

load ('modsim12.RData')

head (modsim.Data $ Def.SR.12)

modsim.Data $ df.12 <- 
  ldply(list (modsim.Data $ Def.12, modsim.Data $ Def.SR.12,
              modsim.Data $ ED.12, modsim.Data $ ED.SR.12,
              modsim.Data $ Sym.12, modsim.Data $ Sym.SR.12))

head (modsim.Data $ df.12)

colnames (modsim.Data $ df.12) [1] <- 'replicates'

nrow(modsim.Data $ df.12)

modsim.Data $ df.12 $ type <-
  rep (c('Local Shape Variables',
         'Interlandmark Distances',
         'Procrustes Residuals'), each = nrow (modsim.Data $ df.12) / 3)

modsim.Data $ df.12 $ type <-
  factor(modsim.Data $ df.12 $ type,
         levels = 
         c('Local Shape Variables', 'Interlandmark Distances', 'Procrustes Residuals'))


modsim.Data $ df.12 $ size <-
  rep (c('Size Retained', 'Size Removed'), each = nrow (modsim.Data $ df.12) / 6,
       times = 3)

modsim.Data $ type1.df <-
  ddply (subset (modsim.Data $ df.12, sample.sizes != 0),
         .(type), summarize, 
         'effect.size' = cut (B12^2, breaks = c(0, summary (B12^2) [-1])),
         'pMIr' = pMIr,
         'pRVr' = pRVr,
         'size' = size,
         'sample.sizes' = sample.sizes)

modsim.Data $ ROC1.df <- 
  ddply (modsim.Data $ type1.df, .(sample.sizes, type, size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.1, by = 0.001),
         'Error.MI' = colSums (outer (pMIr, seq (0.01, 0.1, by = 0.001), '<')) /
         length (pMIr),
         'Error.RV' = colSums (outer (pRVr, seq (0.01, 0.1, by = 0.001), '<')) /
         length(pRVr))


modsim.Data $ ROC1.df $ sig.level <- factor (modsim.Data $ ROC1.df $ sig.level)
modsim.Data $ ROC1.df $ sample.sizes <-
  factor (paste ('n =', modsim.Data $ ROC1.df $ sample.sizes, sep = ' '),
          levels = paste ('n =', c(20, 40, 60, 80, 100), sep = ' '))

modsim.Data $ ROC1.df <- melt (modsim.Data $ ROC1.df)

modsim.Data $ ROC1.df $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC1.df $ sig.level)))

modsim.Data $ ROC1.df $ size <- factor (modsim.Data $ ROC1.df $ size,
                                       levels = c('Size Retained', 'Size Removed'))

head (modsim.Data $ ROC1.df)

modsim.Data $ type1.df <- 
  ddply (subset (modsim.Data $ df.12, sample.sizes != 0), 
         .(sample.sizes, type, size), 
         summarize,
         'Sig.Level' = seq (0.01, 0.1, by = 0.001),
         'Type1.MI' =
         colSums (outer (pMIr, seq (0.01, 0.1, by = 0.001), '<')) / length (pMIr),
         'Type1.RV' =
         colSums (outer (pRVr, seq (0.01, 0.1, by = 0.001), '<')) / length(pRVr))

modsim.Data $ ROC1.df $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC1.df $ sig.level)))

modsim.Data $ ROC1.df $ size <- factor (modsim.Data $ ROC1.df $ size,
                                       levels = c('Size Retained', 'Size Removed'))

modsim.Plots $ Type1.Sym <- 
  ggplot (subset (modsim.Data $ ROC1.df, type == 'Procrustes Residuals')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Type I Error Rate') +
  labs(title = 'Procrustes Residuals')

modsim.Plots $ Type1.ED <- 
  ggplot (subset (modsim.Data $ ROC1.df, type == 'Interlandmark Distances')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  scale_y_continuous(limits = c(0, 0.2)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Type I Error Rate') +
  labs(title = 'Interlandmark Distances')

modsim.Plots $ Type1.Def <- 
  ggplot (subset (modsim.Data $ ROC1.df, type == 'Local Shape Variables')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size, scales = 'free') +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Type I Error Rate') +
  labs(title = 'Local Shape Variables')



modsim.Data $ type1.df $ sample.sizes <- factor (modsim.Data $ type1.df $ sample.sizes)
modsim.Data $ type1.df $ Sig.Level <-
  factor (modsim.Data $ type1.df $ Sig.Level)

modsim.Data $ type1.df <- melt(modsim.Data $ type1.df)

modsim.Data $ type1.df $ Sig.Level <-
  as.numeric (as.character (modsim.Data $ type1.df $ Sig.Level))

modsim.Data $ type1.df $ size <- factor (modsim.Data $ type1.df $ size,
                                         levels = c('Size Retained', 'Size Removed'))

modsim.Data $ type1.df $ type <-
  factor(modsim.Data $ type1.df $ type,
         levels = 
         c('Local Shape Variables', 'Interlandmark Distances', 'Procrustes Residuals'))

head (modsim.Data $ type1.df)


modsim.Plots $ Type1 <- 
  ggplot (modsim.Data $ type1.df) +
  geom_line(aes(x = Sig.Level, y = value, linetype = variable, color = sample.sizes),
            stat = 'identity', position = 'dodge') +
  facet_grid(type ~ size) +
  theme_bw() +
  xlab ('Significance Level') + ylab('Type I Error Rate') +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Sample Size', palette = 'Set1') +
  geom_abline(intercept = 0, slope = 1)


### TYPE 2

modsim.Data $ type2.df <-
  ddply (subset (modsim.Data $ df.12, sample.sizes != 0),
         .(type), summarize, 
         'effect.size' = cut (B12^2, breaks = c(0, summary (B12^2) [-1])),
         'pMIs' = pMIs,
         'pRVs' = pRVs,
         'size' = size,
         'sample.sizes' = sample.sizes)

modsim.Data $ ROC.df <- 
  ddply (modsim.Data $ type2.df, .(sample.sizes, type, size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.1, by = 0.001),
         'Power.MI' = colSums (outer (pMIs, seq (0.01, 0.1, by = 0.001), '<')) /
         length (pMIs),
         'Power.RV' = colSums (outer (pRVs, seq (0.01, 0.1, by = 0.001), '<')) /
         length(pRVs))


modsim.Data $ ROC.df $ sig.level <- factor (modsim.Data $ ROC.df $ sig.level)
modsim.Data $ ROC.df $ sample.sizes <-
  factor (paste ('n =', modsim.Data $ ROC.df $ sample.sizes, sep = ' '),
          levels = paste ('n =', c(20, 40, 60, 80, 100), sep = ' '))

modsim.Data $ ROC.df <- melt (modsim.Data $ ROC.df)

modsim.Data $ ROC.df $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC.df $ sig.level)))

modsim.Data $ ROC.df $ size <- factor (modsim.Data $ ROC.df $ size,
                                       levels = c('Size Retained', 'Size Removed'))

modsim.Plots $ Type2.Sym <- 
  ggplot (subset (modsim.Data $ ROC.df, type == 'Procrustes Residuals')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Procrustes Residuals')

modsim.Plots $ Type2.ED <- 
  ggplot (subset (modsim.Data $ ROC.df, type == 'Interlandmark Distances')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Interlandmark Distances')

modsim.Plots $ Type2.Def <- 
  ggplot (subset (modsim.Data $ ROC.df, type == 'Local Shape Variables')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Local Shape Variables')


modsim.Data $ RV.df <-
  melt(modsim.Data $ df.12 [, c('type', 'sample.sizes', 'RVs', 'RVr')],
       id.vars = c('type', 'sample.sizes'))

modsim.Data $ MI.df <-
  melt (modsim.Data $ df.12 [, c('type', 'sample.sizes', 'MIs', 'MIr')],
        id.vars = c('type', 'sample.sizes'))

modsim.Data $ RV.df $ stat <- rep ('RV Coefficient', nrow (modsim.Data $ RV.df))
modsim.Data $ MI.df $ stat <- rep ('AVG Index', nrow (modsim.Data $ MI.df))

modsim.Data $ values.df <- rbind (subset (modsim.Data $ RV.df, sample.sizes != 0),
                                  subset (modsim.Data $ MI.df, sample.sizes != 0))

modsim.Data $ values.df $ variable <-
  c('Random', 'Structured')[grepl('s', modsim.Data $ values.df $ variable) + 1]

modsim.Data $ values.df $ variable <- factor(modsim.Data $ values.df $ variable,
                                             levels = c('Random', 'Structured'))

head (modsim.Data $ values.df)


modsim.Plots $ stat.dist.sim <- 
  ggplot (modsim.Data $ values.df) +
  geom_density(aes(x = value, fill = interaction(variable, type, sep = ' - '),
                   color = interaction (variable, type, sep = ' - ')),
                 alpha = 0.7, position = 'stack') +
  facet_wrap(~ stat, scales = 'free') +
  scale_fill_brewer(name = 'Simulated Matrix Type', palette = 'Paired') +
  scale_color_brewer(name = 'Simulated Matrix Type', palette = 'Paired') +
  xlab('Value') + ylab ('Density') +
  theme_bw()

save(modcomp.Data, modsim.Data,
     modcomp.Plots, modsim.Plots, file = '../../Tese/ModComp/modcomp.Results.RData')



### mean.df.12

mean.df.12 <-
  ddply (modsim.Data $ df.12, .(type, size, replicates, sample.sizes), summarize,
       'W1' = unique (W1),
       'W2' = unique (W2),
       'B12' = unique (B12),
       'pMIr' = sample(pMIr, 1),
       'pMIs' = sample(pMIs, 1),
       'pRVr' = sample(pRVr, 1),
       'pRVs' = sample(pRVs, 1),
       'nchar' = unique (n.char1))


alt.modsim.Plots <- list()

mean.type1.df <-
  ddply (subset (mean.df.12, sample.sizes != 0),
         .(type), summarize, 
         'effect.size' = cut (B12^2, breaks = c(0, summary (B12^2) [-1])),
         'pMIr' = pMIr,
         'pRVr' = pRVr,
         'size' = size,
         'sample.sizes' = sample.sizes)

mean.ROC1.df <- 
  ddply (mean.type1.df, .(sample.sizes, type, size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.1, by = 0.001),
         'Error.MI' = colSums (outer (pMIr, seq (0.01, 0.1, by = 0.001), '<')) /
         length (pMIr),
         'Error.RV' = colSums (outer (pRVr, seq (0.01, 0.1, by = 0.001), '<')) /
         length(pRVr))


mean.ROC1.df $ sig.level <- factor (mean.ROC1.df $ sig.level)
mean.ROC1.df $ sample.sizes <-
  factor (paste ('n =', mean.ROC1.df $ sample.sizes, sep = ' '),
          levels = paste ('n =', c(20, 40, 60, 80, 100), sep = ' '))

mean.ROC1.df <- melt (mean.ROC1.df)

mean.ROC1.df $ sig.level <-
  as.numeric (as.character (factor (mean.ROC1.df $ sig.level)))

mean.ROC1.df $ size <- factor (mean.ROC1.df $ size,
                                       levels = c('Size Retained', 'Size Removed'))

head (mean.ROC1.df)

mean.ROC1.df $ sig.level <-
  as.numeric (as.character (factor (mean.ROC1.df $ sig.level)))

mean.ROC1.df $ size <- factor (mean.ROC1.df $ size,
                                       levels = c('Size Retained', 'Size Removed'))

alt.modsim.Plots $ Type1.Sym <- 
  ggplot (subset (mean.ROC1.df, type == 'Procrustes Residuals')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.2)) +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Type I Error Rate') +
  labs(title = 'Procrustes Residuals')

alt.modsim.Plots $ Type1.ED <- 
  ggplot (subset (mean.ROC1.df, type == 'Interlandmark Distances')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  scale_y_continuous(limits = c(0, 0.2)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Type I Error Rate') +
  labs(title = 'Interlandmark Distances')

alt.modsim.Plots $ Type1.Def <- 
  ggplot (subset (mean.ROC1.df, type == 'Local Shape Variables')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size, scales = 'free') +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Type I Error Rate') +
  labs(title = 'Local Shape Variables')

mean.type1.df <- 
  ddply (subset (mean.df.12, sample.sizes != 0), 
         .(sample.sizes, type, size), 
         summarize,
         'Sig.Level' = seq (0.01, 0.1, by = 0.001),
         'Type1.MI' =
         colSums (outer (pMIr, seq (0.01, 0.1, by = 0.001), '<')) / length (pMIr),
         'Type1.RV' =
         colSums (outer (pRVr, seq (0.01, 0.1, by = 0.001), '<')) / length(pRVr))

mean.type1.df $ sample.sizes <- factor (mean.type1.df $ sample.sizes)
mean.type1.df $ Sig.Level <-
  factor (mean.type1.df $ Sig.Level)

mean.type1.df <- melt(mean.type1.df)

mean.type1.df $ Sig.Level <-
  as.numeric (as.character (mean.type1.df $ Sig.Level))

mean.type1.df $ size <- factor (mean.type1.df $ size,
                                         levels = c('Size Retained', 'Size Removed'))

mean.type1.df $ type <-
  factor(mean.type1.df $ type,
         levels = 
         c('Local Shape Variables', 'Interlandmark Distances', 'Procrustes Residuals'))


alt.modsim.Plots $ Type1 <- 
  ggplot (mean.type1.df) +
  geom_line(aes(x = Sig.Level, y = value, linetype = variable, color = sample.sizes),
            stat = 'identity', position = 'dodge') +
  facet_grid(type ~ size) +
  theme_bw() +
  xlab ('Significance Level') + ylab('Type I Error Rate') +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Sample Size', palette = 'Set1') +
  geom_abline(intercept = 0, slope = 1)


## alt type 2

mean.type2.df <-
  ddply (subset (mean.df.12, sample.sizes != 0),
         .(type), summarize, 
         'effect.size' = cut (B12^2, breaks = c(0, summary (B12^2) [-1])),
         'pMIs' = pMIs,
         'pRVs' = pRVs,
         'size' = size,
         'sample.sizes' = sample.sizes)

mean.ROC.df <- 
  ddply (mean.type2.df, .(sample.sizes, type, size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.1, by = 0.001),
         'Power.MI' = colSums (outer (pMIs, seq (0.01, 0.1, by = 0.001), '<')) /
         length (pMIs),
         'Power.RV' = colSums (outer (pRVs, seq (0.01, 0.1, by = 0.001), '<')) /
         length(pRVs))


mean.ROC.df $ sig.level <- factor (mean.ROC.df $ sig.level)
mean.ROC.df $ sample.sizes <-
  factor (paste ('n =', mean.ROC.df $ sample.sizes, sep = ' '),
          levels = paste ('n =', c(20, 40, 60, 80, 100), sep = ' '))

mean.ROC.df <- melt (mean.ROC.df)

mean.ROC.df $ sig.level <-
  as.numeric (as.character (factor (mean.ROC.df $ sig.level)))

mean.ROC.df $ size <- factor (mean.ROC.df $ size,
                                       levels = c('Size Retained', 'Size Removed'))

alt.modsim.Plots $ Type2.Sym <- 
  ggplot (subset (mean.ROC.df, type == 'Procrustes Residuals')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Procrustes Residuals')

alt.modsim.Plots $ Type2.ED <- 
  ggplot (subset (mean.ROC.df, type == 'Interlandmark Distances')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Interlandmark Distances')

alt.modsim.Plots $ Type2.Def <- 
  ggplot (subset (mean.ROC.df, type == 'Local Shape Variables')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  scale_linetype_discrete(name = 'Statistic', labels = c('AVGi', 'RV')) +
  scale_color_brewer(name = 'Squared Correlation\nBetween Modules', palette = 'Set1') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Local Shape Variables')


mean.RV.df <-
  melt(mean.df.12 [, c('type', 'sample.sizes', 'replicates', 'RVs', 'RVr')],
       id.vars = c('type', 'sample.sizes', 'replicates'))

mean.MI.df <-
  melt (mean.df.12 [, c('type', 'sample.sizes', 'replicates', 'MIs', 'MIr')],
        id.vars = c('type', 'sample.sizes', 'replicates'))

mean.RV.df $ stat <- rep ('RV Coefficient', nrow (mean.RV.df))
mean.MI.df $ stat <- rep ('AVG Index', nrow (mean.MI.df))

mean.values.df <- rbind (subset (mean.RV.df, sample.sizes != 0),
                                  subset (mean.MI.df, sample.sizes != 0))

mean.values.df $ variable <-
  c('Random', 'Structured')[grepl('s', mean.values.df $ variable) + 1]

mean.values.df $ variable <- factor(mean.values.df $ variable,
                                             levels = c('Random', 'Structured'))

head (mean.values.df)


alt.modsim.Plots $ stat.dist.sim <- 
  ggplot (mean.values.df) +
  geom_density(aes(x = value, fill = interaction(variable, type, sep = ' - '),
                   color = interaction (variable, type, sep = ' - ')),
                 alpha = 0.7, position = 'stack') +
  facet_wrap(~ stat, scales = 'free') +
  scale_fill_brewer(name = 'Simulated Matrix Type', palette = 'Paired') +
  scale_color_brewer(name = 'Simulated Matrix Type', palette = 'Paired') +
  xlab('Value') + ylab ('Density') +
  theme_bw()

save(modcomp.Data, modsim.Data,
     modcomp.Plots, modsim.Plots, file = '../../Tese/ModComp/modcomp.Results.RData')
