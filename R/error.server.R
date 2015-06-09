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
registerDoMC (cores = 100)

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

modsim.Data <-
  within (modsim.Data,
          {
            Def.12 <-
              adply (1:10000, 1, function (i)
                     {
                       print(i)
                       error2112(
                         within = Def.within,
                         between = Def.between,
                         var = Def.var,
                         n.char = 40,
                         sample.sizes = seq (20, 100, 20),
                         it = 1000)
                     }, .parallel = TRUE)

            Def.SR.12 <-
              adply (1:10000, 1, function (i)
                     {
                       print (i)
                       error2112(
                         within = Def.SR.within,
                         between = Def.SR.between,
                         var = Def.var, 
                         n.char = 40,
                         sample.sizes = seq (20, 100, 20),
                         it = 1000)
                     }, .parallel = TRUE)
            
            ED.12 <-
              adply (1:10000, 1, function (i)
                     {
                       print (i)
                       error2112(
                         within = ED.within,
                         between = ED.between,
                         var = ED.var,
                         n.char = 40,
                         sample.sizes = seq (20, 100, 20),
                         it = 1000)
                     }, .parallel = TRUE)
            
            ED.SR.12 <-
              adply (1:10000, 1, function (i)
                     {
                       print (i)
                       error2112(
                         within = ED.SR.within,
                         between = ED.SR.between,
                         var = ED.var, 
                         n.char = 40,
                         sample.sizes = seq (20, 100, 20),
                         it = 1000)
                     }, .parallel = TRUE)
            
            Sym.12 <-
              adply (1:10000, 1, function (i)
                     {
                       print (i)
                       error2112(
                         within = Sym.within,
                         between = Sym.between,
                         var = Sym.var,
                         n.char = 40,
                         sample.sizes = seq (20, 100, 20),
                         it = 1000)
                     }, .parallel = TRUE)
            
            Sym.SR.12 <-
              adply (1:10000, 1, function (i)
                     {
                       print (i)
                       error2112(
                         within = Sym.SR.within,
                         between = Sym.SR.between,
                         var = Sym.var, 
                         n.char = 40,
                         sample.sizes = seq (20, 100, 20),
                         it = 1000
                         )
                     }, .parallel = TRUE)

          })

save (modsim.Data, file = 'modsim12.RData')

