
ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = RV1, y = MI1, color = type, shape = type),
              size = 3, alpha = 0.2) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

## First stuff: coherence

ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = W1, y = MI1, color = type, shape = type),
              size = 3, alpha = 0.2) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

## as expected
ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = W1, y = RV1, color = type, shape = type),
              size = 3, alpha = 0.2) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = B12, y = MI1,
                   color = type, shape = type),
              size = 3, alpha = 0.2) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = B12, y = RV1,
                   color = type, shape = type),
              size = 3, alpha = 0.2) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

print (arrangeGrob(
  ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = B12, y = RV1,
                   color = type, shape = type),
              size = 3, alpha = 0.4) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw(),
  
  ggplot (modsim.Data $ sim.df) +
  geom_boxplot (aes (y = B12, x = type, fill = type), alpha = 0.4, size = 1) +
  scale_fill_brewer(palette = 'Dark2') + coord_flip() +
  theme_bw(), nrow = 2), heights = c(3, 1))

ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = MI1, y = pMI1,
                   color = type, shape = type),
              size = 3, alpha = 0.2) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = RV1, y = pRV1,
                   color = type, shape = type),
              size = 3, alpha = 0.4) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()


ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = B12, y = MI.F,
                   color = type, shape = type),
              size = 3, alpha = 0.4) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = (B12 + B13 + B23)/3, y = pRV.F,
                   color = type, shape = type),
              size = 3, alpha = 0.4) +
  facet_wrap(~ type, scales = 'free_x') +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_boxplot (aes (y = W3, x = pMI3 < 0.05,
                   color = type),
              size = 1, alpha = 0.4) +
  #facet_wrap(~ type, scales = 'free_x') +
  #scale_color_continuous (high = 'green', low = 'red', trans = 'log') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_boxplot (aes (y = B12, x = pRV1 < 0.05, fill = type), size = 1, alpha = 0.4) +
  geom_boxplot (aes (y = ICV, x = pRV1 < 0.05, color = type), size = 1, alpha = 0.4) +
  #facet_wrap(~ type, scales = 'free_x') +
  #scale_color_continuous (high = 'green', low = 'red', trans = 'log') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_boxplot (aes (y = W1, x = type, fill = pMI1 < 0.05), size = 1, alpha = 1) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

ggplot (modsim.Data $ sim.df) +
  geom_boxplot (aes (y = B12, x = type, fill = pRV1 < 0.05), size = 1, alpha = 1) +
  theme_bw()


ggplot (modsim.Data $ sim.df) +
  geom_point (aes (x = (B12+B13+B23)/3, y = pRV.F,
                   color = ICV, shape = type),
              size = 3, alpha = 0.4) +
  #facet_wrap(~ type, scales = 'free_x') +
  scale_color_continuous(high = 'green', low = 'red', trans = 'log') +
  theme_bw()
