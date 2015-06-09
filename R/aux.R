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
require (gridExtra)

require (doMC)
registerDoMC (cores = 8)

attach ('../../Databases/ED.RData')
attach ('../../Databases/Sym.RData')
attach ('../../Databases/Aux.RData')
attach ('../../Databases/LOG.ED.RData')
attach ('../../Databases/NCS.ED.RData')
attach ('../../Databases/OneDef.RData')
attach ('../../Databases/Tree.RData')

source ('../../Func/ModTest.R')

### 2. selecionar dados e métodos
##### a. RV vs. Mantel (vs. Flex?)
##### b. ED, ED noSize, Sym, Def

data.frame (Aux $ sample.size)

### 1. selecionar OTUs

modcomp.Data <- list ()

modcomp.Data $ OTU <- c(names (OneDef) [Aux $ sample.size > 80], ### good samples
                        'Papio_anubis',
                        'Pan_troglodytes',
                        'Cercopithecus_ascanius',
                        'Macaca_fascicularis',
                        'Hylobates_lar',
                        'Saguinus_midas',
                        'Pithecia_pithecia',
                        'Cacajao_calvus',
                        'Chiropotes_chiropotes',
                        'Ateles_geoffroyi')

modcomp.Data $ OTU <- sort (modcomp.Data $ OTU)

modcomp.Data $ otu.table <- data.frame ('Espécie' = modcomp.Data $ OTU, 'n' = Aux $ sample.size [modcomp.Data $ OTU], 'Modelo' = Aux $ data.man.sp [modcomp.Data $ OTU, 2])

rownames (modcomp.Data $ otu.table) <- NULL

xtable(modcomp.Data $ otu.table, digits = 0)

### Matrizes


### Sym: só um alinhamento
modcomp.Data $ config <-
  ldply (Sym [modcomp.Data $ OTU], function (L) two.d.array (L $ sym))

modcomp.Data $ config.3d <- arrayspecs(modcomp.Data $ config [, -1], p = 36, k = 3)

dimnames (modcomp.Data $ config.3d) [1:2] <- dimnames (Sym [[1]] $ sym) [1:2]
                          
modcomp.Data $ sym.gpa <- gpagen(modcomp.Data $ config.3d, ShowPlot = FALSE)

modcomp.Data $ tan <-
  data.frame (
    'OTU' = modcomp.Data $ config [1],
    'logCS' = unlist (llply (Sym[modcomp.Data $ OTU], function (L) log (L $ cs))), 
    two.d.array (modcomp.Data $ sym.gpa $ coords)
  )

colnames (modcomp.Data $ tan) [-(1:2)] <- rownames (Aux $ sym.hyp [[1]])

modcomp.Data $ Sym <- daply (modcomp.Data $ tan, .(.id), 
                             function (D) var (D[, -1]))

mshape (modcomp.Data $ sym.gpa $ coords)


modcomp.Data $ ED <- laply (ED [modcomp.Data $ OTU],
                            function (L) L $ ed.vcv [-20, -20])
                             
modcomp.Data $ Def <- laply (OneDef [modcomp.Data $ OTU],
                             function (L) L $ ml.vcv)

modcomp.Data $ logED <- laply (LOG.ED [modcomp.Data $ OTU],
                               function (L) L $ ed.vcv [-20, -20])

modcomp.Data $ Hyp <- list ('Sym' = Aux $ sym.hyp [[1]], 
                            'EDef' = Aux $ ed.hyp [[1]] [-20, ])

### Vamos arrancar as redundancias das matriz de Sym

modcomp.Data $ Sym <- modcomp.Data $ Sym [, 1:67, 1:67]
modcomp.Data $ Sym <- modcomp.Data $ Sym [, -seq(3, 24, 3), -seq(3, 24, 3)]

dimnames (modcomp.Data $ Sym)

modcomp.Data $ Hyp $ Sym <- modcomp.Data $ Hyp $ Sym [1:66, ]
modcomp.Data $ Hyp $ Sym <- modcomp.Data $ Hyp $ Sym [-seq(2, 23, 3), ]

names (modcomp.Data)

### Size Removal

modcomp.Data <-
  within (modcomp.Data,
          {
### Sym
            dimnames (Sym) <- list (OTU,
                                    c('logCS', rownames (Hyp $ Sym)),
                                    c('logCS', rownames (Hyp $ Sym)))
            
            Sym.SR <- aaply (Sym, 1, RemoveCAC)

### ED
            dimnames (ED) [[1]] <- dimnames (logED) [[1]] <- OTU
            ED.SR <- aaply (ED, 1,
                            function (C) RemoveSize (cov2cor (C)))
            logED.SR <- aaply (logED, 1, function (C) RemoveSize (cov2cor (C)))
            
### Def
            dimnames (Def) [[1]] <- OTU
            Def.SR <- aaply (Def, 1, RemoveCAC)
          })

### Func

modcomp.Data <-
  within (modcomp.Data,
          {
### Sym
            rv.Sym <-
              adply (Sym [, -1, -1], 1,
                     RV, mod.hyp = Hyp $ Sym [, 1:6], proc.res = TRUE,
                     .parallel = TRUE)
            
            rv.Sym.SR <-
              adply (Sym.SR, 1,
                     RV, mod.hyp = Hyp $ Sym [, 1:6], proc.res = TRUE,
                     .parallel = TRUE)

            print ('Sym Done')
### ED
            rv.ED <-
              adply (ED, 1,
                     RV, mod.hyp = Hyp $ EDef [, 1:6], 
                     .parallel = TRUE)
            
            rv.ED.SR <-
              adply (ED.SR, 1,
                     RV, mod.hyp = Hyp $ EDef [, 1:6], 
                     .parallel = TRUE)
            
            print ('ED Done')

### Def
            rv.Def <-
              adply (Def[, -1, -1], 1,
                     RV, mod.hyp = Hyp $ EDef [, 1:6], 
                     .parallel = TRUE)
            
            rv.Def.SR <-
              adply (Def.SR, 1,
                     RV, mod.hyp = Hyp $ EDef [, 1:6], 
                     .parallel = TRUE)

          })

### MI

modcomp.Data <-
  within (modcomp.Data,
          {
### Sym
            mi.Sym <-
              adply (Sym [, -1, -1], 1,
                     AVGIndex, mod.hyp = Hyp $ Sym [, 1:6], proc.res = TRUE,
                     .parallel = TRUE)
            
            mi.Sym.SR <-
              adply (Sym.SR, 1,
                     AVGIndex, mod.hyp = Hyp $ Sym [, 1:6], proc.res = TRUE,
                     .parallel = TRUE)
            
### ED
            mi.ED <-
              adply (ED, 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 1:6], 
                     .parallel = TRUE)

            mi.ED.SR <-
              adply (ED.SR, 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 1:6],
                     .parallel = TRUE)
            
### Def
            mi.Def <-
              adply (Def[, -1, -1], 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 1:6],
                     .parallel = TRUE)
            
            mi.Def.SR <-
              adply (Def.SR, 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 1:6],
                     .parallel = TRUE)
          })

### Data framing

modcomp.Data $ RV.Summ <-
  melt (modcomp.Data [names(modcomp.Data) [grepl('rv', names (modcomp.Data))]])

colnames (modcomp.Data $ RV.Summ) <- c('otu', 'hyp', 'variable', 'value', 'data')

modcomp.Data $ RV.Summ $ data <- gsub('rv.', '', modcomp.Data $ RV.Summ $ data)

modcomp.Data $ RV.Summ $ size <-
  c('Size Retained', 'Size Removed') [grepl ('SR', modcomp.Data $ RV.Summ $ data) + 1]

modcomp.Data $ RV.Summ $ data <- gsub ('\\.SR', '', modcomp.Data $ RV.Summ $ data)

head (modcomp.Data $ RV.Summ)

modcomp.Data $ MI.Summ <-
  melt (modcomp.Data [names(modcomp.Data) [grepl('mi', names (modcomp.Data))]])

colnames (modcomp.Data $ MI.Summ) <- c('otu', 'hyp', 'variable', 'value', 'data')

modcomp.Data $ MI.Summ $ data <- gsub('mi.', '', modcomp.Data $ MI.Summ $ data)

modcomp.Data $ MI.Summ $ size <-
  c('Size Retained', 'Size Removed') [grepl ('SR', modcomp.Data $ MI.Summ $ data) + 1]

modcomp.Data $ MI.Summ $ data <- gsub ('\\.SR', '', modcomp.Data $ MI.Summ $ data)

modcomp.Data $ Summ <-
  dcast (modcomp.Data $ MI.Summ, otu + hyp + size + data ~ variable)
colnames (modcomp.Data $ Summ) [5] <- 'value'
colnames (modcomp.Data $ Summ) [6] <- 'p'

tmp <- dcast (modcomp.Data $ RV.Summ, otu + hyp + size + data ~ variable)
colnames (tmp) [5] <- 'value'
colnames (tmp) [6] <- 'p'

tmp $ type <- rep('RV', nrow(tmp))
modcomp.Data $ Summ $ type <- rep('MI', nrow(modcomp.Data $ Summ))

modcomp.Data $ Summ <-
  rbind (modcomp.Data $ Summ, tmp)

modcomp.Data $ Summ $ hyp <- factor (as.character (modcomp.Data $ Summ $ hyp),
                                     levels = c('Oral', 'Nasal', 'Zygomatic',
                                       'Orbit', 'Base', 'Vault'))

modcomp.Data $ Summ $ size <- factor (as.character (modcomp.Data $ Summ $ size),
                                     levels = c('Size Retained', 'Size Removed'))

modcomp.Data $ Summ $ data <- factor (as.character (modcomp.Data $ Summ $ data),
                                      levels = c('ED', 'Sym', 'Def'))

levels (modcomp.Data $ Summ $ data) <-  c('Interlandmark Distances',
                                          'Procrustes Residuals',
                                          'Local Shape Variables')

### Plots

modcomp.Plots <- list ()

modcomp.Plots $ RV.Func <- 
  ggplot (subset (modcomp.Data $ Summ, type == 'RV')) +
  geom_tile(aes(x = hyp, y = otu, fill = value)) +
  facet_grid(size ~ data) +
  theme_minimal() +
  scale_fill_continuous(name = 'RV', high = 'yellow', low = 'blue', space = 'Lab',
                        limits = c(0, 1), 
                        breaks = c(0.1, 0.5, 0.9)) +
  geom_point (aes (x = hyp, y = otu,
                   size = (p < 0.05) + (p < 0.01) + (p < 0.001),
                   alpha = c(0, 1) [1 + (p < 0.05)]), shape = 21) +
  scale_size_area(name = expression(P(alpha)),
                  labels = c('< 0.05', '< 0.01', '< 0.001'),
                  breaks = c(1, 2, 3)) +
  xlab ('Hypothesis') + ylab ('OTU') + labs(title = 'RV Coefficent') +
  scale_y_discrete(limits = rev(levels(modcomp.Data $ Summ $ otu))) +
  scale_alpha_continuous(limits = c(0, 1)) + guides(alpha = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0))

modcomp.Plots $ MI.Func <- 
  ggplot (subset (modcomp.Data $ Summ, type == 'MI')) +
  geom_tile(aes(x = hyp, y = otu, fill = value)) +
  facet_grid(size ~ data) +
  theme_minimal() +
  scale_fill_continuous(name = 'AVG Index',
                        high = 'blue', low = 'yellow', space = 'Lab',
                        breaks = c(-.3, 0, .3), limits = c(-.4, .4)) +
  geom_point (aes (x = hyp, y = otu,
                   size = (p < 0.05) + (p < 0.01) + (p < 0.001),
                   alpha = c(0, 1) [1 + (p < 0.05)]), shape = 21) +
  scale_size_area(name = expression(P(alpha)),
                  labels = c('< 0.05', '< 0.01', '< 0.001'),
                  breaks = c(1, 2, 3)) +
  xlab ('Hypothesis') + ylab ('OTU') + labs(title = 'AVG Index') +
  scale_alpha_continuous(limits = c(0, 1)) + guides(alpha = FALSE) +
  scale_y_discrete(limits = rev(levels(modcomp.Data $ Summ $ otu))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0))

modcomp.Plots $ Func <- arrangeGrob(modcomp.Plots $ RV.Func,
                                    modcomp.Plots $ MI.Func, ncol = 2)


### DEV

modcomp.Data <-
  within (modcomp.Data,
          {
### Sym
            mi.Sym.dev <-
              adply (Sym [, -1, -1], 1,
                     AVGIndex, mod.hyp = Hyp $ Sym [, 7:8], proc.res = TRUE,
                     hyp.sum = TRUE,
                     .parallel = TRUE)
            
### ED
            mi.ED.dev <-
              adply (ED, 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 7:8], hyp.sum = TRUE,
                     .parallel = TRUE)

### Def
            mi.Def.dev <-
              adply (Def[, -1, -1], 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 7:8], hyp.sum = TRUE,
                     .parallel = TRUE)

            mi.Sym.SR.dev <-
              adply (Sym.SR, 1,
                     AVGIndex, mod.hyp = Hyp $ Sym [, 7:8], proc.res = TRUE,
                     hyp.sum = TRUE,
                     .parallel = TRUE)
            
### ED
            mi.ED.SR.dev <-
              adply (ED.SR, 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 7:8], hyp.sum = TRUE,
                     .parallel = TRUE)

### Def
            mi.Def.SR.dev <-
              adply (Def.SR, 1,
                     AVGIndex, mod.hyp = Hyp $ EDef [, 7:8], hyp.sum = TRUE,
                     .parallel = TRUE)
          })

modcomp.Data <-
  within (modcomp.Data,
          {
### Sym
            rv.Sym.dev <-
              adply (Sym [, -1, -1], 1,
                     RV, mod.hyp = Hyp $ Sym [, 7:8], proc.res = TRUE,
                     .parallel = TRUE)

            rv.Sym.SR.dev <-
              adply (Sym.SR, 1,
                     RV, mod.hyp = Hyp $ Sym [, 7:8], proc.res = TRUE,
                     .parallel = TRUE)

            
            print ('Sym Done')
### ED
            rv.ED.dev <-
              adply (ED, 1,
                     RV, mod.hyp = Hyp $ EDef [, 7:8], 
                     .parallel = TRUE)

            rv.ED.SR.dev <-
              adply (ED.SR, 1,
                     RV, mod.hyp = Hyp $ EDef [, 7:8], 
                     .parallel = TRUE)

            
            print ('ED Done')

### Def
            rv.Def.dev <-
              adply (Def[, -1, -1], 1,
                     RV, mod.hyp = Hyp $ EDef [, 7:8], 
                     .parallel = TRUE)

            rv.Def.SR.dev <-
              adply (Def.SR, 1,
                     RV, mod.hyp = Hyp $ EDef [, 7:8], 
                     .parallel = TRUE)

            
          })

modcomp.Data $ RV.Dev <-
  melt (modcomp.Data [c('rv.Sym.dev', 'rv.ED.dev', 'rv.Def.dev',
                        'rv.Sym.SR.dev', 'rv.ED.SR.dev', 'rv.Def.SR.dev')])

head (modcomp.Data $ RV.Dev)

modcomp.Data $ RV.Dev <-
  modcomp.Data $ RV.Dev [seq(1, nrow (modcomp.Data $ RV.Dev), 2), -2]

colnames (modcomp.Data $ RV.Dev) <- c('otu', 'variable', 'value', 'data')

modcomp.Data $ RV.Dev $ data <- gsub('rv.', '', modcomp.Data $ RV.Dev $ data)
modcomp.Data $ RV.Dev $ data <- gsub('.dev', '', modcomp.Data $ RV.Dev $ data)        

modcomp.Data $ RV.Dev $ hyp <- rep ('Full', nrow (modcomp.Data $ RV.Dev))

modcomp.Data $ RV.Dev $ size <-
  c('Size Retained', 'Size Removed') [grepl ('SR', modcomp.Data $ RV.Dev $ data) + 1]

modcomp.Data $ RV.Dev $ data <- gsub ('\\.SR', '', modcomp.Data $ RV.Dev $ data)

modcomp.Data $ RV.Dev <- modcomp.Data $ RV.Dev [, c(1, 2, 4, 5, 6, 3)]

modcomp.Data $ MI.Dev <-
  melt (modcomp.Data [c('mi.Sym.dev', 'mi.ED.dev', 'mi.Def.dev',
                        'mi.Sym.SR.dev', 'mi.ED.SR.dev', 'mi.Def.SR.dev')])

head (modcomp.Data $ MI.Dev)

colnames (modcomp.Data $ MI.Dev) <- c('otu', 'hyp', 'variable', 'value', 'data')

modcomp.Data $ MI.Dev $ data <- gsub('mi.', '', modcomp.Data $ MI.Dev $ data)
modcomp.Data $ MI.Dev $ data <- gsub('.dev', '', modcomp.Data $ MI.Dev $ data)        

modcomp.Data $ MI.Dev $ size <-
  c('Size Retained', 'Size Removed') [grepl ('SR', modcomp.Data $ MI.Dev $ data) + 1]

modcomp.Data $ MI.Dev $ data <- gsub ('\\.SR', '', modcomp.Data $ MI.Dev $ data)

modcomp.Data $ MI.Dev <- modcomp.Data $ MI.Dev [, c(1, 3, 5, 2, 6, 4)]

modcomp.Data $ Summ.Dev <-
  dcast (modcomp.Data $ MI.Dev, otu + hyp + size + data ~ variable)
colnames (modcomp.Data $ Summ.Dev) [5] <- 'value'
colnames (modcomp.Data $ Summ.Dev) [6] <- 'p'

tmp <- dcast (modcomp.Data $ RV.Dev, otu + hyp + size + data ~ variable)
colnames (tmp) [5] <- 'value'
colnames (tmp) [6] <- 'p'

tmp $ type <- rep('RV', nrow(tmp))
modcomp.Data $ Summ.Dev $ type <- rep('MI', nrow(modcomp.Data $ Summ.Dev))

modcomp.Data $ Summ.Dev <-
  rbind (modcomp.Data $ Summ.Dev, tmp)

modcomp.Data $ Summ.Dev $ hyp <- factor (as.character (modcomp.Data $ Summ.Dev $ hyp),
                                         levels = c('Face', 'Neuro', 'Full'))

modcomp.Data $ Summ.Dev $ data <- factor (as.character (modcomp.Data $ Summ.Dev $ data),
                                      levels = c('ED', 'Sym', 'Def'))

levels (modcomp.Data $ Summ.Dev $ data) <-  c('Interlandmark Distances',
                                              'Procrustes Residuals',
                                              'Local Shape Variables')

levels (modcomp.Data $ Summ.Dev $ hyp) <- c('Face', 'Neuro', 'NeuroFace')

modcomp.Data $ Summ.Dev $ size <- factor (as.character (modcomp.Data $ Summ.Dev $ size),
                                     levels = c('Size Retained', 'Size Removed'))

modcomp.Plots $ RV.NeuroFace <- 
  ggplot (subset (modcomp.Data $ Summ.Dev, type == 'RV')) +
  geom_tile(aes(x = data, y = otu, fill = value)) +
  facet_wrap(~ size, nrow = 2) +
  theme_minimal() +
  scale_fill_continuous (name = 'RV', high = 'yellow', low = 'blue', space = 'Lab',
                         limits = c(0, 1), breaks = c(0.1, 0.5, 0.9)) +
  geom_point (aes (x = data, y = otu,
                   size = (p < 0.05) + (p < 0.01) + (p < 0.001),
                   alpha = c(0, 1) [1 + (p < 0.05)]), shape = 21) +
  scale_size_area(name = expression(P(alpha)),
                  labels = c('< 0.05', '< 0.01', '< 0.001'),
                  breaks = c(1, 2, 3)) +
  xlab ('Data Type') + ylab ('') + labs(title = 'Face/Neuro RV') +
  scale_y_discrete(limits = rev(levels(modcomp.Data $ Summ $ otu))) +
  scale_alpha_continuous(limits = c(0, 1)) + guides(alpha = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0))

modcomp.Plots $ MI.Dev <- 
  ggplot (subset (modcomp.Data $ Summ.Dev, type == 'MI')) +
  geom_tile(aes(x = hyp, y = otu, fill = value)) +
  facet_grid(size ~ data) +
  theme_minimal() +
  scale_fill_continuous(name = 'AVG Index',
                        high = 'blue', low = 'yellow', space = 'Lab',
                        breaks = c(-.3, 0, .3), limits = c(-.4, .4)) +
  geom_point (aes (x = hyp, y = otu,
                   size = (p < 0.05) + (p < 0.01) + (p < 0.001),
                   alpha = c(0, 1) [1 + (p < 0.05)]), shape = 21) +
  scale_size_area(name = expression(P(alpha)),
                  labels = c('< 0.05', '< 0.01', '< 0.001'),
                  breaks = c(1, 2, 3)) +
  xlab ('Hypothesis') + ylab ('') + labs(title = 'Face/Neuro AVG Index') +
  scale_alpha_continuous(limits = c(0, 1)) + guides(alpha = FALSE) +
  scale_y_discrete(limits = rev(levels(modcomp.Data $ Summ $ otu))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(size = 0))


modcomp.Plots $ Dev <- arrangeGrob(modcomp.Plots $ RV.NeuroFace,
                                   modcomp.Plots $ MI.Dev, ncol = 2, widths = c(3, 5))

