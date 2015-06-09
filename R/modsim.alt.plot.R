
modsim.Data $ ROC.df.RVes <- 
  ddply (modsim.Data $ type2.df, .(sample.sizes, type, size, effect.size.RV), 
         summarize,
         'sig.level' = seq (0.01, 0.1, by = 0.001),
         'Power.MI' = colSums (outer (pMIs, seq (0.01, 0.1, by = 0.001), '<')) /
         length (pMIs),
         'Power.RV' = colSums (outer (pRVs, seq (0.01, 0.1, by = 0.001), '<')) /
         length(pRVs))


modsim.Data $ ROC.df.RVes $ sig.level <- factor (modsim.Data $ ROC.df.RVes $ sig.level)
modsim.Data $ ROC.df.RVes $ sample.sizes <-
  factor (paste ('n =', modsim.Data $ ROC.df.RVes $ sample.sizes, sep = ' '),
          levels = paste ('n =', c(20, 40, 60, 80, 100), sep = ' '))

modsim.Data $ ROC.df.RVes <- melt (modsim.Data $ ROC.df.RVes)

modsim.Data $ ROC.df.RVes $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC.df.RVes $ sig.level)))

modsim.Data $ ROC.df.RVes $ size <- factor (modsim.Data $ ROC.df.RVes $ size,
                                       levels = c('Size Retained', 'Size Removed'))

modsim.Plots $ Type2.Sym.RV <- 
  ggplot (subset (modsim.Data $ ROC.df.RVes, type == 'Procrustes Residuals')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size.RV, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1, size = 2, alpha = 0.5) +
  theme_bw() +
  scale_linetype_discrete(name = 'Test', labels = c('Mantel', 'RV')) +
  scale_color_brewer(name = 'RV Value', palette = 'Dark2') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Procrustes Residuals')

modsim.Plots $ Type2.ED.RV <- 
  ggplot (subset (modsim.Data $ ROC.df.RVes, type == 'Interlandmark Distances')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size.RV, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1, size = 2, alpha = 0.5) +
  theme_bw() +
  scale_linetype_discrete(name = 'Test', labels = c('Mantel', 'RV')) +
  scale_color_brewer(name = 'RV Value', palette = 'Dark2') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Interlandmark Distances')

modsim.Plots $ Type2.Def.RV <- 
  ggplot (subset (modsim.Data $ ROC.df.RVes, type == 'Local Shape Variables')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size.RV, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1, size = 2, alpha = 0.5) +
  theme_bw() +
  scale_linetype_discrete(name = 'Test', labels = c('Mantel', 'RV')) +
  scale_color_brewer(name = 'RV Value', palette = 'Dark2') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Local Shape Variables')


modsim.Data $ ROC.df.MIes <- 
  ddply (modsim.Data $ type2.df, .(sample.sizes, type, size, effect.size.MI), 
         summarize,
         'sig.level' = seq (0.01, 0.1, by = 0.001),
         'Power.MI' = colSums (outer (pMIs, seq (0.01, 0.1, by = 0.001), '<')) /
         length (pMIs),
         'Power.RV' = colSums (outer (pRVs, seq (0.01, 0.1, by = 0.001), '<')) /
         length(pRVs))


modsim.Data $ ROC.df.MIes $ sig.level <- factor (modsim.Data $ ROC.df.MIes $ sig.level)
modsim.Data $ ROC.df.MIes $ sample.sizes <-
  factor (paste ('n =', modsim.Data $ ROC.df.MIes $ sample.sizes, sep = ' '),
          levels = paste ('n =', c(20, 40, 60, 80, 100), sep = ' '))

modsim.Data $ ROC.df.MIes <- melt (modsim.Data $ ROC.df.MIes)

modsim.Data $ ROC.df.MIes $ sig.level <-
  as.numeric (as.character (factor (modsim.Data $ ROC.df.MIes $ sig.level)))

modsim.Data $ ROC.df.MIes $ size <- factor (modsim.Data $ ROC.df.MIes $ size,
                                       levels = c('Size Retained', 'Size Removed'))

modsim.Plots $ Type2.Sym.MI <- 
  ggplot (subset (modsim.Data $ ROC.df.MIes, type == 'Procrustes Residuals')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size.MI, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1, size = 2, alpha = 0.5) +
  theme_bw() +
  scale_linetype_discrete(name = 'Test', labels = c('Mantel', 'RV')) +
  scale_color_brewer(name = 'MI Value', palette = 'Dark2') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Procrustes Residuals')

modsim.Plots $ Type2.ED.MI <- 
  ggplot (subset (modsim.Data $ ROC.df.MIes, type == 'Interlandmark Distances')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size.MI, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1, size = 2, alpha = 0.5) +
  theme_bw() +
  scale_linetype_discrete(name = 'Test', labels = c('Mantel', 'RV')) +
  scale_color_brewer(name = 'MI Value', palette = 'Dark2') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Interlandmark Distances')

modsim.Plots $ Type2.Def.MI <- 
  ggplot (subset (modsim.Data $ ROC.df.MIes, type == 'Local Shape Variables')) +
  geom_line(aes (x = sig.level, y = value,
                 color = effect.size.MI, linetype = variable)) +
  facet_grid (sample.sizes ~ size) +
  geom_abline(intercept = 0, slope = 1, size = 2, alpha = 0.5) +
  theme_bw() +
  scale_linetype_discrete(name = 'Test', labels = c('Mantel', 'RV')) +
  scale_color_brewer(name = 'MI Value', palette = 'Dark2') +
  xlab('Significance Level') + ylab ('Power') +
  labs(title = 'Local Shape Variables')

