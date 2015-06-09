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

load ('Alt.Mantel.RData')
load ('Mantel.RData')
load ('rv.comp.RData')
load ('def.comp.RData')

Mantel <- list ()
Mantel $ ED <- list ()
Mantel $ ED $ Func <- ED.Mantel $ Func
Mantel $ ED $ Dev <- ED.Mantel $ Dev

Mantel $ LOG.ED <- list ()
Mantel $ LOG.ED $ Func <- Log.Mantel $ Func
Mantel $ LOG.ED $ Dev <- Log.Mantel $ Dev

Mantel $ NCS.ED <- list ()
Mantel $ NCS.ED $ Func <- Ncs.Mantel $ Func
Mantel $ NCS.ED $ Dev <- Ncs.Mantel $ Dev

Mantel $ Sym <- list ()
Mantel $ Sym $ Func <- Sym.Mantel $ Func
Mantel $ Sym $ Dev <- Sym.Mantel $ Dev

Mantel $ Def <- list ()
Mantel $ Def $ Func <- def.comp $ Mantel.Func
Mantel $ Def $ Dev <- def.comp $ Mantel.Dev

RV <- list ()
RV $ ED <- rv.comp $ ed
RV $ LOG.ED <- rv.comp $ log.ed
RV $ NCS.ED <- rv.comp $ ncs.ed
RV $ Sym <- rv.comp $ gm
RV $ Def <- def.comp $ RV

save (Mantel, file = 'Mantel.RData')
save (RV, file = 'RV.RData')


