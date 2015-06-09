# require (plotrix)
require (ape)
require (Morphometrics)
require (geiger)
require (bbmle)
require (nlme)
require (plyr)
require (geomorph)
require (numDeriv)
require (expm)
# require (xtable)

require (doMC)
registerDoMC (cores = 10)

attach ('../../Databases/ED.RData')
attach ('../../Databases/Sym.RData')
attach ('../../Databases/Aux.RData')
attach ('../../Databases/LOG.ED.RData')
attach ('../../Databases/NCS.ED.RData')
attach ('../../Databases/Def.RData')

source ('../../Func/RV.new.R')

ICV <- function (cov.matrix, remove.size = FALSE)
  {
    if (remove.size)
      cov.matrix <- RemoveSize (cov.matrix)
    eval <- eigen (cov.matrix) $ values
    return (sd (eval) / mean (eval))
  }

CV <- function (x) sd (x) / mean (x)

fisherTrans <- function (x) 0.5 * (log (1 + x) - log (1 - x))

## rv.comp <- list ()

## rv.comp $ gm.func <- 
##   laply (Sym, wrapRV, hyps = Aux $ simple.hyp $ Func, .parallel = TRUE)

## rv.comp $ gm.dev <- 
##   laply (Sym, wrapRV, hyps = Aux $ simple.hyp $ Dev, .parallel = TRUE)

## rv.comp $ ed <- 
##   laply (ED, wrapRV.vcv,
##          hyps = Aux $ ed.hyp [[1]], n.it = 10000, .parallel = TRUE)

## rv.comp $ log.ed <- 
##   laply (LOG.ED, wrapRV.vcv,
##          hyps = Aux $ ed.hyp [[1]], n.it = 10000, .parallel = TRUE)

## rv.comp $ ncs.ed <- 
##   laply (NCS.ED, wrapRV.vcv,
##          hyps = Aux $ ed.hyp [[1]], n.it = 10000, .parallel = TRUE)

## save (rv.comp, file = 'rv.comp.RData')

attach ('rv.comp.RData')
detach (file:rv.comp.RData)

integration <- 
    cbind (laply (ED, function (x) ICV (x $ ed.vcv)),
           laply (LOG.ED, function (x) ICV (x $ ed.vcv)),
           laply (NCS.ED, function (x) ICV (x $ ed.vcv)),
           laply (Sym, function (x) CV (x $ cs)))

colnames (integration) <- c ('ICV', 'ICV.log', 'ICV.ncs', 'CSCV')

dim (rv.comp $ log.ed)

par (mfrow = c(2, 4))
for (i in 1:8)
plot (fisherTrans (rv.comp $ gm [, i , 'RV']) ~ integration [, 'CSCV'], xlab = 'CSCV', ylab = 'RV',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (rv.comp $ ed [, i , 'RV'] ~ log (integration [, 'ICV']), xlab = 'ICV', ylab = 'RV',
      main = dimnames (rv.comp $ ed) [[2]] [i], pch = ifelse (rv.comp $ ed [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (log (rv.comp $ gm [, i , 'Between']) ~ log (rv.comp $ gm [, i , 'Within1']), xlab = 'Within.Out', ylab = 'Between',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (log (rv.comp $ gm [, i , 'Between']) ~ log (rv.comp $ gm [, i , 'Within2']), xlab = 'Within.In', ylab = 'Between',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (log (rv.comp $ gm [, i , 'Between']) ~ log (sqrt (rv.comp $ gm [, i , 'Within1'] * rv.comp $ gm [, i , 'Within2'])),
      xlab = 'Within', ylab = 'Between',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (log (rv.comp $ gm [, i , 'Between']) ~ log (integration [, 'CSCV']), xlab = 'CSCV', ylab = 'RV',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (log (rv.comp $ gm [, i , 'Within2']) ~ log (integration [, 'CSCV']), xlab = 'CSCV', ylab = 'RV',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 4))
for (i in 1:8)
plot (log (rv.comp $ gm [, i , 'RV']) ~ log (rv.comp $ ncs.ed [, i, 'RV']), xlab = 'RV ED', ylab = 'RV GM',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

Sym <- llply (Sym, function (x)
              {
                  x <- within (x,
                               {
                                   gpa <- gpagen (sym, ShowPlot = FALSE)
                                   tan <- t (t (two.d.array (gpa $ coords)) - c (mshape (gpa $ coords)))
                                   tan.vcv <- cov (tan)
                                   dimnames (tan.vcv) [[1]] <- dimnames (tan.vcv) [[1]] <-
                                       paste (rep (rownames (sym), each = 3), rep (c ('X', 'Y', 'Z'),
                                                                       times = nrow (sym)),
                                              sep = '.')
                               })
                  x
              })

save (Sym, file = '../../Databases/Sym.RData')

Sym.Mantel <- list ();
Sym.Mantel $ Func <-
    laply (Sym, function (x) ModNew (x $ tan.vcv [1:66, 1:66], Aux $ sym.hyp [[1]] [1:66, 1:6], n.it = 10000));
Sym.Mantel $ Dev <-
    laply (Sym, function (x) ModNew (x $ tan.vcv, Aux $ sym.hyp [[1]] [1:66, 7:8], n.it = 10000));

ED.Mantel <- list ();
ED.Mantel $ Func <-
    laply (ED, function (x) ModNew (x $ ed.vcv, Aux $ ed.hyp [[1]] [, 1:6], n.it = 10000));
ED.Mantel $ Dev <-
    laply (ED, function (x) ModNew (x $ ed.vcv, Aux $ ed.hyp [[1]] [, 7:8], n.it = 10000));



attach ('Mantel.RData')

par (mfrow = c(2, 3))
for (i in 1:6)
    plot (Sym.Mantel $ Func [, i, 'ModIndex'] ~
          I (rv.comp $ gm [, i, 'Between'] / rv.comp $ gm [, i, 'Within2']))

par (mfrow = c(2, 3))
for (i in 1:6)
    plot (I (ED.Mantel $ Func [, i, 'AVG+'] - ED.Mantel $ Func [, i, 'AVG-']) ~
          rv.comp $ ed [, i, 'RV'])

par (mfrow = c(2, 3))
for (i in 1:6)
    plot (ED.Mantel $ Func [, i, 'Gamma'] ~ rv.comp $ ed [, i, 'RV'])


attach ('def.comp.RData')

par (mfrow = c(2, 4))
for (i in 1:8)
    {
        tmp <- cbind (def.comp $ RV [, i, 'pvalue'] < 0.05, rv.comp $ ed [, i, 'pvalue'] < 0.05)
        plot (def.comp $ RV [, i, 'RV'] ~ rv.comp $ ed [, i, 'RV'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'))
    }

par (mfrow = c(2, 4))
for (i in 1:8)
    {
        tmp <- cbind (def.comp $ RV [, i, 'pvalue'] < 0.05, rv.comp $ gm [, i, 'pvalue'] < 0.05)
        plot (def.comp $ RV [, i, 'RV'] ~ rv.comp $ gm [, i, 'RV'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'))
    }

par (mfrow = c(2, 4))
for (i in 1:8)
    {
        tmp <- cbind (rv.comp $ ed [, i, 'pvalue'] < 0.05, rv.comp $ gm [, i, 'pvalue'] < 0.05)
        plot (rv.comp $ ed [, i, 'RV'] ~ rv.comp $ gm [, i, 'RV'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'))
    }


par (mfrow = c(2, 5))
for (i in 1:7)
    {
        tmp <- cbind (def.comp $ Mantel.Func [, i, 'Probability'] < 0.05,
                      Re (Sym.Mantel $ Func [, i, 'Probability']) < 0.05)
        plot (def.comp $ Mantel.Func [, i, 'ModIndex'] ~ Sym.Mantel $ Func [, i, 'ModIndex'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'),
              main = dimnames (def.comp $ Mantel.Func) [[2]] [i])
    }
for (i in 1:3)
    {
        tmp <- cbind (def.comp $ Mantel.Dev [, i, 'Probability'] < 0.05,
                      Re (Sym.Mantel $ Dev [, i, 'Probability']) < 0.05)
        plot (def.comp $ Mantel.Dev [, i, 'ModIndex'] ~ Sym.Mantel $ Dev [, i, 'ModIndex'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'),
              main = dimnames (def.comp $ Mantel.Dev) [[2]] [i])
    }

par (mfrow = c(2, 5))
for (i in 1:7)
    {
        tmp <- cbind (Re (Sym.Mantel $ Func [, i, 'Probability']) < 0.05,
                      ED.Mantel $ Func [, i, 'Probability'] < 0.05)
        plot (Sym.Mantel $ Func [, i, 'ModIndex'] ~ ED.Mantel $ Func [, i, 'ModIndex'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'),
              main = dimnames (Sym.Mantel $ Func) [[2]] [i])
    }
for (i in 1:3)
    {
        tmp <- cbind (Re (Sym.Mantel $ Dev) [, i, 'Probability'] < 0.05,
                      ED.Mantel $ Dev [, i, 'Probability'] < 0.05)
        plot (Sym.Mantel $ Dev [, i, 'ModIndex'] ~ ED.Mantel $ Dev [, i, 'ModIndex'],
              col = ifelse (apply (tmp, 1, all), 'darkgreen', 'red'),
              main = dimnames (def.comp $ Mantel.Dev) [[2]] [i])
    }



par (mfrow = c(2, 4))
for (i in 1:8)
plot (rv.comp $ gm [, i , 'RV'] ~ log (integration [, 'CSCV']), xlab = 'CSCV', ylab = 'RV',
      main = dimnames (rv.comp $ gm) [[2]] [i], pch = ifelse (rv.comp $ gm [, i , 'pvalue'] < 0.05, 17, 21))

par (mfrow = c(2, 5))
for (i in 1:7)
plot (ED.Mantel $ Func [, i , 'ModIndex'] ~ log (integration [, 'ICV']), xlab = 'CSCV', ylab = 'ED.ModIndex',
      main = dimnames (ED.Mantel $ Func) [[2]] [i], pch = ifelse (ED.Mantel $ Func [, i , 2] < 0.05, 17, 21))
for (i in 1:3)
plot (ED.Mantel $ Dev [, i , 'ModIndex'] ~ log (integration [, 'ICV']), xlab = 'CSCV', ylab = 'ED.ModIndex',
      main = dimnames (ED.Mantel $ Dev) [[2]] [i], pch = ifelse (ED.Mantel $ Dev [, i , 2] < 0.05, 17, 21))

