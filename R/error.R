error2112 <- function (within, between, var,
                       n.char = 20, n.hyp = 2, sample.sizes = c(10, 20, 40),
                       it = 1000)
  {
    ### build hypothesis
    repeat {
      base.hyp <- NULL
      for (i in 1:n.char)
        base.hyp <- rbind (base.hyp, sample (c(rep (0, times = n.hyp - 1), 1), n.hyp))

      colnames (base.hyp) <- 1:ncol (base.hyp)
      
      if(all (colSums (base.hyp) > 5)) break
      else (print ('hyp with few traits, sample it again'))
    }

    ### build structured matrix
    base.tensor <- outer (base.hyp, base.hyp)
      
    within.dist <- atanh (within)
    between.dist <- atanh (between)
    
    repeat {
      
      within.cor <- tanh (rnorm (n.hyp, mean (within.dist), sd (within.dist)))
      between.cor <- tanh (rnorm (choose (n.hyp, 2), mean (between.dist),
                                  sd(between.dist)))
      
      names (within.cor) <- paste ('W', 1:n.hyp, sep = '')
    
      cor.mat <- array (0, c(n.char, n.char))
      for (i in 1:n.hyp)
        cor.mat <- cor.mat + within.cor [i] * base.tensor [, i, , i]
      
      k <- 1
      for (i in 1:(n.hyp-1))
        for (j in (i+1):n.hyp)
          {
            cor.mat <- cor.mat + between.cor[k] * (base.tensor[, i, , j] +
                                                   base.tensor[, j, , i])
            names (between.cor) [k] <- paste ('B', i, j, sep = '')
            k <- k + 1
          }

      out.cor <- NULL
      
      for(i in 1:(length (sample.sizes) + 1))
        out.cor <- rbind (out.cor, c(within.cor, between.cor))

      colnames (out.cor) <- c(names (within.cor), names (between.cor))
      
      diag(cor.mat) <- 1
      
      tryCatch (expr = {chol(cor.mat); break},
                error = function (cond) print ('not PD. do it again'))
      
    }
    
    stdev <- sqrt (sample (var, n.char))

    vcv.mat.str <- diag(stdev) %*% cor.mat %*% diag(stdev)

    ### build (pseudo) random matrix
    repeat {
      shuffle <- sample (n.char)
     
      vcv.mat.ran <- vcv.mat.str[shuffle, shuffle]

      tryCatch (expr = {chol(vcv.mat.ran); break},
                error = function (cond) print ('not PD 2. do it again'))
    }

    MI.real.str <- AVGIndex (vcv.mat.str, mod.hyp = base.hyp,
                             n.it = it, hyp.sum = T) [n.hyp + 1, -1]

    RV.real.str <- RV (vcv.mat.str, base.hyp, n.it = it) [1, -1]

    MI.real.ran <- AVGIndex (vcv.mat.ran,
                             base.hyp, n.it = it, hyp.sum = T) [n.hyp + 1, -1]

    RV.real.ran <- RV (vcv.mat.ran, base.hyp, n.it = it) [1, -1]

    ICV.real.str <- ICV (vcv.mat.str)
    ICV.real.ran <- ICV (vcv.mat.ran)

    ### sample
    sampled.vcv.str <- 
      alply (sample.sizes, 1, function (i)
             var (rmvnorm(i, sigma = vcv.mat.str)))
    
    sampled.vcv.ran <-
      alply (sample.sizes, 1, function (i)
             var (rmvnorm(i, sigma = vcv.mat.ran)))


    ### mod tests (full integration only)
    MI.str <-
      ldply (sampled.vcv.str,
             function (C)
             AVGIndex (C, base.hyp,
                       n.it = it, hyp.sum = TRUE) [n.hyp + 1, -1]) [, -1]
    
    MI.ran <-
      ldply (sampled.vcv.ran, function (C)
             AVGIndex (C, base.hyp,
                       n.it = it, hyp.sum = TRUE) [n.hyp + 1, -1]) [, -1]

    RV.str <-
      ldply (sampled.vcv.str, function (C)
             RV (C, base.hyp, n.it = it) [1, -1]) [, -1]

    RV.ran <-
      ldply (sampled.vcv.ran, function (C)
             RV (C, base.hyp, n.it = it) [1, -1]) [, -1]
    
    icv.str <- laply (sampled.vcv.str, ICV)
    icv.ran <- laply (sampled.vcv.ran, ICV)

    out <- cbind (rep(sum(base.hyp [, 1]), length (sample.sizes)), 
                  sample.sizes, MI.str, MI.ran,
                  RV.str, RV.ran, icv.str, icv.ran)
               
    colnames (out) <-
      c('n.char1', 'sample.sizes', 'MIs', 'pMIs', 'MIr','pMIr',
        'RVs', 'pRVs', 'RVr', 'pRVr', 'ICVs', 'ICVr')

    real.out <- c(sum(base.hyp [, 1]), 0, MI.real.str, MI.real.ran,
                  RV.real.str, RV.real.ran, ICV.real.str, ICV.real.ran)

    names (real.out) <-
      c('n.char1', 'sample.sizes', 'MIs', 'pMIs', 'MIr','pMIr',
        'RVs', 'pRVs', 'RVr', 'pRVr', 'ICVs', 'ICVr')
    
    out <- rbind (out, real.out)
    
    out <- cbind (out.cor, out)
               
    out
    
  }


