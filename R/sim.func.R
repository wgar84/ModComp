GenerateMatrix <- function (within, between, var, base.hyp)
  {
    n.hyp <- ncol (base.hyp)
    n.char <- nrow (base.hyp)
    base.tensor <- outer (base.hyp, base.hyp)

    within.dist <- atanh (within)
    between.dist <- atanh (between)
    
    within.cor <- tanh (rnorm (n.hyp, mean (within.dist), sd (within.dist)))
    between.cor <- tanh (rnorm (choose (n.hyp, 2), mean (between.dist), sd (between.dist)))
        
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

    diag(cor.mat) <- 1

    stdev <- sqrt (sample (var, n.char))

    var.mat <- diag(stdev) %*% cor.mat %*% diag(stdev)

    icv <- ICV(var.mat) [1]
    
    RV <- wrapRV.vcv (var.mat, base.hyp, n.it = 10000)
    MI <- ModIndex (var.mat, base.hyp, iterations = 10000)[, c('ModIndex', 'Probability')]
    
    out <- c(within.cor, between.cor, icv, RV [, 1], RV [, 2], MI [, 1], MI [, 2])
    names (out) <- c(names (within.cor), names (between.cor), 'ICV',
                     paste ('RV', 1:n.hyp, sep = ''), 'RV.F',
                     paste ('pRV', 1:n.hyp, sep = ''), 'pRV.F',
                     paste ('MI', 1:n.hyp, sep = ''), 'MI.F',
                     paste ('pMI', 1:n.hyp, sep = ''), 'pMI.F')

    unlist (out)
  }

MatrixSim <- function (within, between, var, base.hyp,
                       n.mat = 1000, parallel = TRUE)
  adply (1:n.mat, 1, function (i)
         GenerateMatrix (within, between, var, base.hyp), .parallel = parallel)
   
isPD <- function (WB.cor, base.hyp)
  {

    n.char <- nrow (base.hyp)
    n.hyp <- ncol (base.hyp)

    within.cor <- WB.cor[1:n.hyp]
    between.cor <- WB.cor[-(1:n.hyp)]
    
    base.tensor <- outer (base.hyp, base.hyp)
    
    cor.mat <- array (0, c(n.char, n.char))
    for (i in 1:n.hyp)
      cor.mat <- cor.mat + within.cor [i] * base.tensor [, i, , i]

    k <- 1
    for (i in 1:(n.hyp-1))
      for (j in (i+1):n.hyp)
        {
          cor.mat <- cor.mat + between.cor[k] * (base.tensor[, i, , j] +
                                                 base.tensor[, j, , i])
          names (between.cor) [k] <- paste ('B', i, '-', j, sep = '')
          k <- k + 1
        }

    diag(cor.mat) <- 1
    
    all (eigen (cor.mat) $ values > 0) #?

  }

