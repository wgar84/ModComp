require (plotrix)
require (ape)
require (Morphometrics)
require (geiger)
require (bbmle)
require (nlme)
require (plyr)
require (geomorph)
require (numDeriv)
require (expm)
require (xtable)
require (ggplot2)

require (doMC)
registerDoMC (cores = 3)

attach ('../../Databases/ED.RData')
attach ('../../Databases/Sym.RData')
attach ('../../Databases/Aux.RData')
attach ('../../Databases/LOG.ED.RData')
attach ('../../Databases/NCS.ED.RData')
attach ('../../Databases/Def.RData')
attach ('../../Databases/Tree.RData')


source ('../../Func/RV.new.R')

attach ('Mantel.RData')
attach ('RV.RData')
attach ('Flex.RData')

ICV <- function (cov.matrix, remove.size = FALSE)
  {
    if (remove.size)
      cov.matrix <- RemoveSize (cov.matrix)
    eval <- eigen (cov.matrix) $ values
    return (sd (eval) / mean (eval))
  }

CV <- function (x) sd (x) / mean (x)

fisherTrans <- function (x) 0.5 * (log (1 + x) - log (1 - x))

integration <- 
    cbind (laply (ED, function (x) ICV (x $ ed.vcv)),
           laply (LOG.ED, function (x) ICV (x $ ed.vcv)),
           laply (NCS.ED, function (x) ICV (x $ ed.vcv)),
           laply (Sym, function (x) CV (x $ cs)))

colnames (integration) <- c ('ICV', 'ICV.log', 'ICV.ncs', 'CSCV')

pdf ('../../10 Parcial 3/Fig/RV_sym_vs_ncs.pdf', width = 10, height = 15)
par (mfrow = c(4, 2))
for (i in 1:8)
    {
        plot (RV $ Sym [, i , 'RV'] ~ RV $ NCS.ED [, i, 'RV'],
              xlab = 'RV NCS ED', ylab = 'RV Sym',
              main = dimnames (RV $ Sym) [[2]] [i],
              col = ifelse (RV $ Sym [, i , 'pvalue'] < 0.05, 'blue', 'red'),
              pch = ifelse (RV $ NCS.ED [, i , 'pvalue'] < 0.05, 17, 21))
    }
legend ('bottomright', legend = c ('nsig', 'sig Sym', 'sig NCS ED', 'sig'),
        pch = c (21, 21, 17, 17), col = c('red', 'blue'), bty = 'n')
dev.off (dev.cur ())

pdf ('../../10 Parcial 3/Fig/RV_sym_vs_ed.pdf', width = 10, height = 15)
par (mfrow = c(4, 2))
for (i in 1:8)
    {
        plot (RV $ Sym [, i , 'RV'] ~ RV $ ED [, i, 'RV'],
              xlab = 'RV ED', ylab = 'RV Sym',
              main = dimnames (RV $ Sym) [[2]] [i],
              col = ifelse (RV $ Sym [, i , 'pvalue'] < 0.05, 'blue', 'red'),
              pch = ifelse (RV $ ED [, i , 'pvalue'] < 0.05, 17, 21))
    }
legend ('bottomright', legend = c ('nsig', 'sig Sym', 'sig ED', 'sig'),
        pch = c (21, 21, 17, 17), col = c('red', 'blue'), bty = 'n')
dev.off (dev.cur ())

pdf (file = '../../10 Parcial 3/Fig/RVSym.pdf', width = 16, height = 16)
pairs (RV $ Sym [, 1:6, 'RV'], upper.panel = NULL, pch = 20)
dev.off (dev.cur ())

pdf (file = '../../10 Parcial 3/Fig/RVED.pdf', width = 16, height = 16)
pairs (RV $ ED [, 1:6, 'RV'], upper.panel = NULL, pch = 20)
dev.off (dev.cur ())

pdf (file = '../../10 Parcial 3/Fig/RVNCS.pdf', width = 16, height = 16)
pairs (RV $ NCS.ED [, 1:6, 'RV'], upper.panel = NULL, pch = 20)
dev.off (dev.cur ())

pairs (Mantel $ ED $ Func [, 1:6, 'AVG+'])

wrapRV.vcv(list (ed.vcv = cov2cor (Sym [[1]] $ tan.vcv)), Aux $ sym.hyp [[1]])

sample.size <- laply (Sym, function (x) length (x $ cs))

par (mfrow = c(2, 4))
for (i in 1:8)
    {
        plot (abs (RV $ Sym [, i, 'RV'] - RV $ NCS.ED [, i, 'RV']) ~ sample.size)
    }



pdf (file = '../../10 Parcial 3/Fig/MantelSym.pdf', width = 16, height = 16)
pairs (Mantel $ Sym $ Func [, 1:6, 'ModIndex'], upper.panel = NULL, pch = 20)
dev.off (dev.cur ())

pdf (file = '../../10 Parcial 3/Fig/MantelED.pdf', width = 16, height = 16)
pairs (Mantel $ ED $ Func [, 1:6, 'ModIndex'], upper.panel = NULL, pch = 20)
dev.off (dev.cur ())

pdf (file = '../../10 Parcial 3/Fig/MantelNCS.pdf', width = 16, height = 16)
pairs (Mantel $ NCS.ED $ Func[, 1:6, 'ModIndex'], upper.panel = NULL, pch = 20)
dev.off (dev.cur ())

for (i in 1:8)
  print (cor (RV $ Sym [, i, 'RV'], RV $ ED [, i, 'RV']) ^ 2)

for (i in 1:8)
  print (cor (RV $ Sym [, i, 'RV'], RV $ NCS.ED [, i, 'RV']) ^ 2)

for (i in 1:8)
  print (summary (lm (RV $ Sym [, i, 'RV'] ~ RV $ ED [, i, 'RV'])))


pdf ('../../10 Parcial 3/Fig/Mantel_sym_vs_ed.pdf', width = 10, height = 15)
par (mfrow = c(4, 2))
for (i in 1:6)
    {
        plot (Mantel $ Sym $ Func [, i , 'ModIndex'] ~
              Mantel $ ED $ Func [, i, 'ModIndex'],
              xlab = 'MI ED', ylab = 'MI Sym',
              main = dimnames (Mantel $ Sym $ Func) [[2]] [i],
              col = ifelse (Mantel $ Sym $ Func[, i , 'Probability'] < 0.05, 'blue', 'red'),
              pch = ifelse (Mantel $ ED $ Func [, i , 'Probability'] < 0.05, 17, 21))
    }
for (i in 1:2)
    {
        plot (Mantel $ Sym $ Dev [, i , 'ModIndex'] ~
              Mantel $ ED $ Dev [, i, 'ModIndex'],
              xlab = 'MI ED', ylab = 'MI Sym',
              main = dimnames (Mantel $ Sym $ Dev) [[2]] [i],
              col = ifelse (Mantel $ Sym $ Dev[, i , 'Probability'] < 0.05, 'blue', 'red'),
              pch = ifelse (Mantel $ ED $ Dev [, i , 'Probability'] < 0.05, 17, 21))
    }
legend ('topright', legend = c ('nsig', 'sig Sym', 'sig NCS ED', 'sig'),
        pch = c (21, 21, 17, 17), col = c('red', 'blue'), bty = 'n')
dev.off (dev.cur ())


pdf ('../../10 Parcial 3/Fig/Mantel_sym_vs_ncs.pdf', width = 10, height = 15)
par (mfrow = c(4, 2))
for (i in 1:6)
    {
        plot (Mantel $ Sym $ Func [, i , 'ModIndex'] ~
              Mantel $ NCS.ED $ Func [, i, 'ModIndex'],
              xlab = 'MI ED', ylab = 'MI Sym',
              main = dimnames (Mantel $ Sym $ Func) [[2]] [i],
              col = ifelse (Mantel $ Sym $ Func[, i , 'Probability'] < 0.05, 'blue', 'red'),
              pch = ifelse (Mantel $ NCS.ED $ Func [, i , 'Probability'] < 0.05, 17, 21))
    }
for (i in 1:2)
    {
        plot (Mantel $ Sym $ Dev [, i , 'ModIndex'] ~
              Mantel $ NCS.ED $ Dev [, i, 'ModIndex'],
              xlab = 'MI ED', ylab = 'MI Sym',
              main = dimnames (Mantel $ Sym $ Dev) [[2]] [i],
              col = ifelse (Mantel $ Sym $ Dev[, i , 'Probability'] < 0.05, 'blue', 'red'),
              pch = ifelse (Mantel $ NCS.ED $ Dev [, i , 'Probability'] < 0.05, 17, 21))
    }
legend ('topright', legend = c ('nsig', 'sig Sym', 'sig NCS ED', 'sig'),
        pch = c (21, 21, 17, 17), col = c('red', 'blue'), bty = 'n')
dev.off (dev.cur ())

for (i in 1:8)
  print (sum (RV $ Sym [, i , 'pvalue'] < 0.05 & RV $ ED [, i , 'pvalue'] < 0.05) / 109 +
         sum (RV $ Sym [, i , 'pvalue'] > 0.05 & RV $ ED [, i , 'pvalue'] > 0.05) / 109)

for (i in 1:8)
  print (sum (RV $ Sym [, i , 'pvalue'] < 0.05 & RV $ NCS.ED [, i , 'pvalue'] < 0.05) / 109 +
         sum (RV $ Sym [, i , 'pvalue'] > 0.05 & RV $ NCS.ED [, i , 'pvalue'] > 0.05) / 109)

RV $ prop.symed <- NULL
for (i in 1:8)
    RV $ prop.symed <-
    cbind (RV $ prop.symed,
           c (sum (RV $ Sym [, i , 'pvalue'] < 0.05 & RV $ ED [, i , 'pvalue'] < 0.05),
              sum (RV $ Sym [, i , 'pvalue'] < 0.05 & RV $ ED [, i , 'pvalue'] > 0.05),
              sum (RV $ Sym [, i , 'pvalue'] > 0.05 & RV $ ED [, i , 'pvalue'] < 0.05),
              sum (RV $ Sym [, i , 'pvalue'] > 0.05 & RV $ ED [, i , 'pvalue'] > 0.05)))

RV $ prop.symncs <- NULL
for (i in 1:8)
    RV $ prop.symncs <-
    cbind (RV $ prop.symncs,
           c (sum (RV $ Sym [, i , 'pvalue'] < 0.05 & RV $ NCS.ED [, i , 'pvalue'] < 0.05),
              sum (RV $ Sym [, i , 'pvalue'] < 0.05 & RV $ NCS.ED [, i , 'pvalue'] > 0.05),
              sum (RV $ Sym [, i , 'pvalue'] > 0.05 & RV $ NCS.ED [, i , 'pvalue'] < 0.05),
              sum (RV $ Sym [, i , 'pvalue'] > 0.05 & RV $ NCS.ED [, i , 'pvalue'] > 0.05)))

Mantel $ prop.symed <- NULL
for (i in 1:6)
    Mantel $ prop.symed <-
    cbind (Mantel $ prop.symed,
           c (sum (Mantel $ Sym $ Func [, i , 'Probability'] < 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Func [, i , 'Probability'] < 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] > 0.05),
              sum (Mantel $ Sym $ Func [, i , 'Probability'] > 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Func [, i , 'Probability'] > 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] > 0.05)))
for (i in 1:2)
    Mantel $ prop.symed <-
    cbind (Mantel $ prop.symed,
           c (sum (Mantel $ Sym $ Dev [, i , 'Probability'] < 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Dev [, i , 'Probability'] < 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] > 0.05),
              sum (Mantel $ Sym $ Dev [, i , 'Probability'] > 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Dev [, i , 'Probability'] > 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] > 0.05)))

Mantel $ prop.symncs <- NULL
for (i in 1:6)
    Mantel $ prop.symncs <-
    cbind (Mantel $ prop.symncs,
           c (sum (Mantel $ Sym $ Func [, i , 'Probability'] < 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Func [, i , 'Probability'] < 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] > 0.05),
              sum (Mantel $ Sym $ Func [, i , 'Probability'] > 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Func [, i , 'Probability'] > 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] > 0.05)))
for (i in 1:2)
    Mantel $ prop.symncs <-
    cbind (Mantel $ prop.symncs,
           c (sum (Mantel $ Sym $ Dev [, i , 'Probability'] < 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Dev [, i , 'Probability'] < 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] > 0.05),
              sum (Mantel $ Sym $ Dev [, i , 'Probability'] > 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] < 0.05),
              sum (Mantel $ Sym $ Dev [, i , 'Probability'] > 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] > 0.05)))

dimnames (RV $ prop.symed) <- list (c ('sig', 'sig Sym', 'sig ED', 'nsig'),
                                    dimnames (RV $ Sym) [[2]])

dimnames (RV $ prop.symncs) <- list (c ('sig', 'sig Sym', 'sig NCS', 'nsig'),
                                    dimnames (RV $ Sym) [[2]])

dimnames (Mantel $ prop.symed) <- list (c ('sig', 'sig Sym', 'sig ED', 'nsig'),
                                    dimnames (RV $ Sym) [[2]])

dimnames (Mantel $ prop.symncs) <- list (c ('sig', 'sig Sym', 'sig NCS', 'nsig'),
                                    dimnames (RV $ Sym) [[2]])

dimnames (Mantel $ prop.symncs) [[2]] [3] <- 'Zygo'
dimnames (Mantel $ prop.symed) [[2]] [3] <- 'Zygo'
dimnames (RV $ prop.symncs) [[2]] [3] <- 'Zygo'
dimnames (RV $ prop.symed) [[2]] [3] <- 'Zygo'

Mantel $ prop.symncs <- Mantel $ prop.symncs [c(1, 4, 2, 3), ]
Mantel $ prop.symed <- Mantel $ prop.symed [c(1, 4, 2, 3), ]
RV $ prop.symncs <- RV $ prop.symncs [c(1, 4, 2, 3), ]
RV $ prop.symed <- RV $ prop.symed [c(1, 4, 2, 3), ]


pdf (file = '../../10 Parcial 3/Fig/sigcompmeth.pdf', width = 14, height = 12)
par (mfrow = c(2, 2))
barplot (RV $ prop.symed/109, legend.text = rownames (RV $ prop.symed),
         xlim = c (0, 11.5), main = '', xlab = 'Hipótese', ylab = 'Proporção',
         args.legend = list ('bty' = 'n'))
legend ('bottomright', legend = '(a)', cex = 2.5, bty = 'n')
barplot (RV $ prop.symncs/109, legend.text = rownames (RV $ prop.symncs),
         args.legend = list ('bty' = 'n'), xlab = 'Hipótese', ylab = 'Proporção',
         xlim = c (0, 11.5), main = '')
legend ('bottomright', legend = '(b)', cex = 2.5, bty = 'n')
barplot (Mantel $ prop.symed/109, legend.text = rownames (Mantel $ prop.symed),
         xlim = c (0, 11.5), main = '', xlab = 'Hipótese', ylab = 'Proporção',
         args.legend = list ('bty' = 'n'))
legend ('bottomright', legend = '(c)', cex = 2.5, bty = 'n')
barplot (Mantel $ prop.symncs/109, legend.text = rownames (Mantel $ prop.symncs),
         args.legend = list ('bty' = 'n'), xlab = 'Hipótese', ylab = 'Proporção',
         xlim = c (0, 11.5), main = '')
legend ('bottomright', legend = '(d)', cex = 2.5, bty = 'n')
dev.off (dev.cur ())         

RVvsMantel <- list ()
RVvsMantel $ prop.sym <- NULL
for (i in 1:6)
    RVvsMantel $ prop.sym <-
    cbind (RVvsMantel $ prop.sym,
           c (sum (RV $ Sym [, i , 2] < 0.05 &
                   Mantel $ Sym $ Func [, i , 'Probability'] < 0.05),
              sum (RV $ Sym [, i , 2] < 0.05 &
                   Mantel $ Sym $ Func [, i , 'Probability'] > 0.05),
              sum (RV $ Sym [, i , 2] > 0.05 &
                   Mantel $ Sym $ Func [, i , 'Probability'] < 0.05),
              sum (RV $ Sym [, i , 2] > 0.05 &
                   Mantel $ Sym $ Func [, i , 'Probability'] > 0.05)))
for (i in 1:2)
    RVvsMantel $ prop.sym <-
    cbind (RVvsMantel $ prop.sym,
           c (sum (RV $ Sym [, i+6 , 2] < 0.05 &
                   Mantel $ Sym $ Dev [, i , 'Probability'] < 0.05),
              sum (RV $ Sym [, i+6 , 2] < 0.05 &
                   Mantel $ Sym $ Dev [, i , 'Probability'] > 0.05),
              sum (RV $ Sym [, i+6 , 2] > 0.05 &
                   Mantel $ Sym $ Dev [, i , 'Probability'] < 0.05),
              sum (RV $ Sym [, i+6 , 2] > 0.05 &
                   Mantel $ Sym $ Dev [, i , 'Probability'] > 0.05)))

RVvsMantel $ prop.ed <- NULL
for (i in 1:6)
    RVvsMantel $ prop.ed <-
    cbind (RVvsMantel $ prop.ed,
           c (sum (RV $ ED [, i , 2] < 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] < 0.05),
              sum (RV $ ED [, i , 2] < 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] > 0.05),
              sum (RV $ ED [, i , 2] > 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] < 0.05),
              sum (RV $ ED [, i , 2] > 0.05 &
                   Mantel $ ED $ Func [, i , 'Probability'] > 0.05)))
for (i in 1:2)
    RVvsMantel $ prop.ed <-
    cbind (RVvsMantel $ prop.ed,
           c (sum (RV $ ED [, i+6 , 2] < 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] < 0.05),
              sum (RV $ ED [, i+6 , 2] < 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] > 0.05),
              sum (RV $ ED [, i+6 , 2] > 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] < 0.05),
              sum (RV $ ED [, i+6 , 2] > 0.05 &
                   Mantel $ ED $ Dev [, i , 'Probability'] > 0.05)))

RVvsMantel $ prop.ncs.ed <- NULL
for (i in 1:6)
    RVvsMantel $ prop.ncs.ed <-
    cbind (RVvsMantel $ prop.ncs.ed,
           c (sum (RV $ NCS.ED [, i , 2] < 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] < 0.05),
              sum (RV $ NCS.ED [, i , 2] < 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] > 0.05),
              sum (RV $ NCS.ED [, i , 2] > 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] < 0.05),
              sum (RV $ NCS.ED [, i , 2] > 0.05 &
                   Mantel $ NCS.ED $ Func [, i , 'Probability'] > 0.05)))
for (i in 1:2)
    RVvsMantel $ prop.ncs.ed <-
    cbind (RVvsMantel $ prop.ncs.ed,
           c (sum (RV $ NCS.ED [, i+6 , 2] < 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] < 0.05),
              sum (RV $ NCS.ED [, i+6 , 2] < 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] > 0.05),
              sum (RV $ NCS.ED [, i+6 , 2] > 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] < 0.05),
              sum (RV $ NCS.ED [, i+6 , 2] > 0.05 &
                   Mantel $ NCS.ED $ Dev [, i , 'Probability'] > 0.05)))


RVvsMantel <- llply (RVvsMantel, function (x)
    {
        dimnames (x) <- list (c ('sig', 'sig RV', 'sig Mantel', 'nsig'),
                              dimnames (RV $ Sym) [[2]])
        dimnames (x) [[2]] [3] <- 'Zygo'
        x <- x [c (1, 4, 2, 3),]
        x
    })

pdf (file = '../../10 Parcial 3/Fig/sigcomprep.pdf', width = 14, height = 12)
par (mfrow = c(2, 2))
barplot (RVvsMantel $ prop.sym/109, legend.text = rownames (RVvsMantel $ prop.sym),
         xlim = c (0, 11.5), main = '', xlab = 'Hipótese', ylab = 'Proporção',
         args.legend = list ('bty' = 'n'))
legend ('bottomright', legend = '(a)', cex = 2.5, bty = 'n')
barplot (RVvsMantel $ prop.ed/109, legend.text = rownames (RVvsMantel $ prop.ed),
         args.legend = list ('bty' = 'n'), xlab = 'Hipótese', ylab = 'Proporção',
         xlim = c (0, 11.5), main = '')
legend ('bottomright', legend = '(b)', cex = 2.5, bty = 'n')
barplot (RVvsMantel $ prop.ncs.ed/109, legend.text = rownames (RVvsMantel $ prop.ncs.ed),
         xlim = c (0, 11.5), main = '', xlab = 'Hipótese', ylab = 'Proporção',
         args.legend = list ('bty' = 'n'))
legend ('bottomright', legend = '(c)', cex = 2.5, bty = 'n')
dev.off (dev.cur ())

PTable <- function (values, probs, italic = 0.1, bold = 0.05)
  {
    values <- round (values, 4)
    values <- ifelse (probs < bold, paste ('\\textbf{', values, '}', sep = ''),
                      ifelse (probs < italic, paste ('\\textit{', values, '}', sep = ''),
                              values))
    values <- data.frame (values)
    rownames (values) <- paste ('\\textit{', rownames (values), '}', sep = '')
    print (xtable (values), sanitize.text.function = function (x) {x})
  }

PTable (RV $ Sym [, 1:7, 'RV'], RV $ Sym [, 1:7, 'pvalue'])
PTable (RV $ ED [, , 'RV'], RV $ ED [, , 'pvalue'])
PTable (RV $ NCS.ED [, , 'RV'], RV $ NCS.ED [, , 'pvalue'])


dimnames (RV $ Sym) [[1]] <- names (Sym)
dimnames (RV $ ED) [[1]] <- names (Sym)
dimnames (RV $ NCS.ED) [[1]] <- names (Sym)


tmp <- cbind (Mantel $ Sym $ Func [, 1:6, 'ModIndex'],
              Mantel $ Sym $ Dev [, 1:3, 'ModIndex'])
tmp.p <- cbind (Mantel $ Sym $ Func [, 1:6, 'Probability'],
              Mantel $ Sym $ Dev [, 1:3, 'Probability'])
colnames (tmp) [9] <- 'NeuroFace'

rownames (tmp) <- names (Sym)
rownames (tmp.p) <- names (Sym)

PTable (tmp, tmp.p)

tmp <- cbind (Mantel $ ED $ Func [, 1:6, 'ModIndex'],
              Mantel $ ED $ Dev [, 1:3, 'ModIndex'])
tmp.p <- cbind (Mantel $ ED $ Func [, 1:6, 'Probability'],
              Mantel $ ED $ Dev [, 1:3, 'Probability'])
colnames (tmp) [9] <- 'NeuroFace'

rownames (tmp) <- names (ED)
rownames (tmp.p) <- names (ED)

PTable (tmp, tmp.p)

tmp <- cbind (Mantel $ NCS.ED $ Func [, 1:6, 'ModIndex'],
              Mantel $ NCS.ED $ Dev [, 1:3, 'ModIndex'])
tmp.p <- cbind (Mantel $ NCS.ED $ Func [, 1:6, 'Probability'],
              Mantel $ NCS.ED $ Dev [, 1:3, 'Probability'])
colnames (tmp) [9] <- 'NeuroFace'

rownames (tmp) <- names (NCS.ED)
rownames (tmp.p) <- names (NCS.ED)

PTable (tmp, tmp.p)
