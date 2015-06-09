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

source ('../../Func/RV.new.R')

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

modsim.Plots $ cor.hist <- 
  ggplot (subset (modsim.Data $ cor.df, t1.num < t2.num)) +
  geom_histogram(aes(x = value, fill = type),
               alpha = 0.4, position = 'identity') +
  scale_fill_brewer(name = 'Data Type', palette = 'Paired') +
  theme_bw() + xlab ('Correlation') + ylab('Count')

modsim.Plots $ cormat.ED <- 
ggplot (subset (modsim.Data $ cor.df, type == 'ED')) +
  geom_tile(aes (x = t1, y = t2, fill = value)) +
  facet_wrap (~ otu) +
  scale_y_discrete(limits = rev(levels (modsim.Data $ cor.df $ t2)) [1:38]) +
  theme_bw() +
  theme(axis.text.x = element_text (size = 3, angle = 90),
        axis.text.y = element_text (size = 3)) +
  scale_fill_continuous(low = 'yellow', high = 'blue', space = 'Lab')

modsim.Plots $ cormat.Def <- 
ggplot (subset (modsim.Data $ cor.df, type == 'Def')) +
  geom_tile(aes (x = t1, y = t2, fill = value)) +
  facet_wrap (~ otu) +
  scale_y_discrete(limits = rev(levels (modsim.Data $ cor.df $ t2)) [1:38]) +
  theme_bw() +
  theme(axis.text.x = element_text (size = 3, angle = 90),
        axis.text.y = element_text (size = 3)) +
  scale_fill_continuous(low = 'yellow', high = 'blue', space = 'Lab')

modsim.Plots $ cormat.Sym <- 
ggplot (subset (modsim.Data $ cor.df, type == 'Sym')) +
  geom_tile(aes (x = t1, y = t2, fill = value)) +
  facet_wrap (~ otu) +
  scale_y_discrete(limits = rev(levels (modsim.Data $ cor.df $ t2)) [39:96]) +
  theme_bw() +
  theme(axis.text.x = element_text (size = 3, angle = 90),
        axis.text.y = element_text (size = 3)) +
  scale_fill_continuous(low = 'yellow', high = 'blue', space = 'Lab')


#print (arrangeGrob(modsim.Plots $ cormat.ED,
#            modsim.Plots $ cormat.Def,
#            modsim.Plots $ cormat.Sym, ncol = 3))


### SIM
source ('sim.func.R')
source ('error.R')

modsim.Data $ Sym.var <- subset (modsim.Data $ var.df, type == 'Sym') $ value
modsim.Data $ ED.var <- subset (modsim.Data $ var.df, type == 'ED') $ value
modsim.Data $ Def.var <- subset (modsim.Data $ var.df, type == 'Def') $ value

modsim.Data <-
  within(modsim.Data, 
         {

           sim.parm <- cbind (rep (c (2, 3, 6), times = 300),
                              rep (c (24, 48, 60), each = 300))

           colnames (sim.parm) <- c('n.hyp', 'n.char')
           
           type1.Def <-
             adply (sim.parm, 1, 
                    function (L) type1error (Def.var, n.char = L [2], n.hyp = L [1], 
                                             sample.sizes = seq(20, 200, 20)),
                    .parallel = TRUE)

           print ('Def done')

           type1.ED <-
             adply (sim.parm, 1, 
                    function (L) type1error (ED.var, n.char = L [2], n.hyp = L [1], 
                                             sample.sizes = seq(20, 200, 20)),
                    .parallel = TRUE)

           print ('ED done')

           type1.Sym <-
             adply (sim.parm, 1, 
                    function (L) type1error (Sym.var, n.char = L [2], n.hyp = L [1], 
                                             sample.sizes = seq(20, 200, 20)),
                    .parallel = TRUE)
           
         })

modsim.Data $ hyp.tensor.Sym <- outer (modcomp.Data $ Hyp $ Sym [, 1:6],
                                       modcomp.Data $ Hyp $ Sym [, 1:6])

modsim.Data $ hyp.tensor.EDef <- outer (modcomp.Data $ Hyp $ EDef [, 1:6],
                                        modcomp.Data $ Hyp $ EDef [, 1:6])

modsim.Data $ Sym.within <- c()
modsim.Data $ Sym.between <- c()

modsim.Data $ ED.within <- c()
modsim.Data $ ED.between <- c()

modsim.Data $ Def.within <- c()
modsim.Data $ Def.between <- c()

for (i in 1:21)
  for (j in 1:5)
    for (k in j:6)
      {
       tmp.cor <- cov2cor (modcomp.Data $ Sym [i, , ]) [which (
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
       
       tmp.cor <- cov2cor (modcomp.Data $ Def [i, , ]) [which (
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
      }

### matplot ggplot

modsim.Data <-
  within(modsim.Data, 
         {

           sim.parm <- rep (c (24, 48, 60), each = 300)
                              
           type2.Def <-
             adply (sim.parm, 1, 
                    function (i) type2error (Def.within, Def.between,
                                             Def.var, n.char = i,  
                                             sample.sizes = seq(20, 200, 20)),
                    .parallel = TRUE)

           print ('Def done')

           type2.ED <-
             adply (sim.parm, 1, 
                    function (i) type2error (ED.within, ED.between,
                                             ED.var, n.char = i,  
                                             sample.sizes = seq(20, 200, 20)),
                    .parallel = TRUE)

           print ('ED done')

           type2.Sym <-
             adply (sim.parm, 1, 
                    function (L) type2error (Sym.within, Sym.between,
                                             Sym.var, n.char = i, 
                                             sample.sizes = seq(20, 200, 20)),
                    .parallel = TRUE)
           
         })

# save (modsim.Data, file = 'type2.RData')

load ('type1.RData')

type1.df <- rbind (modsim.Data $ type1.Def, modsim.Data $ type1.ED,
                   modsim.Data $ type1.Sym)

nrow (type1.df)

type1.df $ type <- rep (c('Def', 'ED', 'Sym'), each = 18000)

head (type1.df)

load ('type2.RData')

type2.df <- rbind (modsim.Data $ type2.Def, modsim.Data $ type2.ED,
                   modsim.Data $ type2.Sym)

type2.df $ type <- rep (c('Def', 'ED', 'Sym'), each = 9000)

nrow (type2.df)

modsim.Data $ type1.df <- type1.df
modsim.Data $ type2.df <- type2.df

rm (type1.df, type2.df)

modsim.Data <-
  within(modsim.Data, 
         {

           sim.parm <- cbind (rep (c (2, 3, 6), times = 300),
                              rep (c (24, 48, 60), each = 300))
         })

modsim.Data $ type1.df $ n.hyp <-
  factor(modsim.Data $ sim.parm [as.numeric (modsim.Data $ type1.df $ X1), 1])

modsim.Data $ type1.df $ sample.size <-
  factor (modsim.Data $ type1.df $ sample.size)

modsim.Data $ type1.df $ n.traits <- 
  modsim.Data $ type1.df $ n.traits [(1:54000) + c(1, -1)]
  
modsim.Data $ type1.df.full <- subset (modsim.Data $ type1.df, hyp == 'Full')

nrow (modsim.Data $ type1.df.full)

modsim.Data $ type1.df.full $ n.traits <-
  factor (modsim.Data $ type1.df.full $ n.traits)

table(modsim.Data $ type1.df.full [, c('n.hyp', 'n.traits', 'sample.size', 'type')])

modsim.Data $ type1.psumm <- ddply (modsim.Data $ type1.df.full,
                                    .(sample.size, n.traits), summarize,
                                    'prop.MI' = (sum (pMI < 0.05) / length (pMI)) - 0.05,
                                    'prop.RV' = (sum (pRV < 0.05) / length (pRV)) - 0.05)

modsim.Data $ type1.psumm <- melt (modsim.Data $ type1.psumm)

modsim.Plots $ type1error <- 
  ggplot (modsim.Data $ type1.psumm) +
  geom_bar(aes (x = sample.size, y = value, fill = variable),
           stat = 'identity', position = 'dodge') +
  facet_wrap(~ n.traits) +
  theme_bw() + xlab('Sample Size') + ylab(expression(Delta ~ P(alpha))) +
  labs(title = 'Type I Error') +
  scale_fill_discrete(name = 'Test', labels = c('Mantel', 'RV'))

modsim.Data $ type2.df $ RV <- as.numeric(modsim.Data $ type2.df $ RV)

modsim.Data $ type2.df $ pRV <- as.numeric(modsim.Data $ type2.df $ pRV)

modsim.Data $ type2.df $ sample.size <-
  factor (modsim.Data $ type2.df $ sample.size)

modsim.Data $ type2.df $ n.char1 <-
  factor (modsim.Data $ type2.df $ n.char1)

modsim.Plots $ type2.Mantel.hyp1 <- 
  ggplot (modsim.Data $ type2.df) +
  geom_boxplot (aes(y = W1 - B, x = sample.size, fill = Hyp1_pMI < 0.05),
                width = 0.5, outlier.shape = '+') +
  facet_wrap(~ type) + theme_bw() +
  scale_fill_brewer(palette = 'Dark2')

modsim.Plots $ type2.RV.hyp1 <- 
ggplot (modsim.Data $ type2.df) +
  geom_boxplot (aes(y = B - W1, x = sample.size, fill = pRV < 0.05),
                width = 0.5, outlier.shape = '+') +
  facet_wrap(~ type) +
  scale_fill_brewer(palette = 'Dark2') + theme_bw()

modsim.Data $ type2.df.Sym <- within (subset (modsim.Data $ type2.df, type == 'Sym'),
                                              effect.size <- cut(B, breaks = 4))

modsim.Data $ ROC.df.Sym <- 
  ddply (modsim.Data $ type2.df.Sym, .(sample.size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.05, by = 0.001),
         'TPR.MI' = colSums (outer (Full_pMI, seq (0.01, 0.05, by = 0.001), '<')) /
         length (Full_pMI),
         'TPR.RV' = colSums (outer (pRV, seq (0.01, 0.05, by = 0.001), '<')) /
         length(Full_pMI))

modsim.Data $ ROC.df.Sym $ sig.level <- factor (modsim.Data $ ROC.df.Sym $ sig.level)

modsim.Data $ ROC.df.Sym <- melt (modsim.Data $ ROC.df.Sym)

modsim.Data $ ROC.df.Sym $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC.df.Sym $ sig.level)))

modsim.Plots $ ROC.Sym <- 
  ggplot (modsim.Data $ ROC.df.Sym) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_wrap ( ~ sample.size, ncol = 2) +
  geom_abline(intercept = 0, slope = 1) + theme_bw()

modsim.Data $ type2.df.Def <- within (subset (modsim.Data $ type2.df, type == 'Def'),
                                      effect.size <- cut(B, breaks = 4))

modsim.Data $ ROC.df.Def <- 
  ddply (modsim.Data $ type2.df.Def, .(sample.size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.05, by = 0.001),
         'TPR.MI' = colSums (outer (Full_pMI, seq (0.01, 0.05, by = 0.001), '<')) /
         length (Full_pMI),
         'TPR.RV' = colSums (outer (pRV, seq (0.01, 0.05, by = 0.001), '<')) /
         length(Full_pMI))

modsim.Data $ ROC.df.Def $ sig.level <- factor (modsim.Data $ ROC.df.Def $ sig.level)

modsim.Data $ ROC.df.Def <- melt (modsim.Data $ ROC.df.Def)

modsim.Data $ ROC.df.Def $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC.df.Def $ sig.level)))

modsim.Plots $ ROC.Def <- 
ggplot (modsim.Data $ ROC.df.Def) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
#  geom_point(aes (x = sig.level, y = value,
#                 color = effect.size, shape = variable)) +
  facet_wrap ( ~ sample.size, ncol = 2) +
  geom_abline(intercept = 0, slope = 1) + theme_bw()


modsim.Data $ type2.df.ED <- within (subset (modsim.Data $ type2.df, type == 'ED'),
                                              effect.size <- cut(B, breaks = 4))

modsim.Data $ ROC.df.ED <- 
  ddply (modsim.Data $ type2.df.ED, .(sample.size, effect.size), 
         summarize,
         'sig.level' = seq (0.01, 0.05, by = 0.001),
         'TPR.MI' = colSums (outer (Full_pMI, seq (0.01, 0.05, by = 0.001), '<')) /
         length (Full_pMI),
         'TPR.RV' = colSums (outer (pRV, seq (0.01, 0.05, by = 0.001), '<')) /
         length(Full_pMI))

modsim.Data $ ROC.df.ED $ sig.level <- factor (modsim.Data $ ROC.df.ED $ sig.level)

modsim.Data $ ROC.df.ED <- melt (modsim.Data $ ROC.df.ED)

modsim.Data $ ROC.df.ED $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC.df.ED $ sig.level)))

modsim.Plots $ ROC.ED <- 
ggplot (modsim.Data $ ROC.df.ED) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size, linetype = variable)) +
  facet_wrap ( ~ sample.size, ncol = 2) +
  geom_abline(intercept = 0, slope = 1) + theme_bw()

### Pela conversa com a Bárbara, acho que eu precisava casar os dois testes em uma coisa só, bem simétrica... ;)


