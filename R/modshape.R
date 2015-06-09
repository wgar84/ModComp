ModIndex <- function (cov.ss, mod.hyp, hyp.sum = FALSE, proc.res = FALSE, n.it = 1000)
  {
    n.traits <- nrow (mod.hyp)
    if (any (rownames (cov.ss) == 'logCS'))
      cov.ss <- cov.ss [- which (rownames (cov.ss) == 'logCS'),
                        - which (rownames (cov.ss) == 'logCS')]
    cov.ss.eval <- eigen (cov.ss) $ values
    icv <- sd (cov.ss.eval) / mean (cov.ss.eval)
    if (proc.res)
      {
        landmarks <- do.call (rbind, strsplit (rownames (cov.ss), '\\.')) [, 1]
        landmarks <- as.numeric (as.factor (landmarks))
      }
    cor.ss <- cov2cor (cov.ss)
    cor.ss.vec <- cor.ss [lower.tri (cor.ss)]
    if (hyp.sum)
      {
        mod.hyp <- cbind (mod.hyp, rep (0, times = nrow (mod.hyp)))
        colnames (mod.hyp) [ncol (mod.hyp)] <- 'Full'
      }
    hyp.mat <- aaply (mod.hyp, 2, function (H) H %*% t(H))
    if (hyp.sum)
      {
        hyp.mat[ncol(mod.hyp), , ] <- aaply (hyp.mat[- ncol(mod.hyp), , ], c(2, 3), sum)
        hyp.mat[ncol(mod.hyp), , ] <- (hyp.mat[ncol(mod.hyp), , ] >= 1) * 1
      }
    real.mod <-
      aaply (hyp.mat, 1,
             function (H.mat)
             {
               H.mat.vec <- H.mat [lower.tri (H.mat)]
               avg.plus <- mean (cor.ss.vec [H.mat.vec == 1])
               avg.minus <- mean (cor.ss.vec [H.mat.vec != 1])
               c(avg.plus, avg.minus, (avg.plus - avg.minus) / icv)
             })
    perm.mod.single <- function (H.mat)
      {
        if (proc.res)
          {
            n.land <- max (landmarks)
            perm.lm <- sample(1:n.land)
            use.as.index <-
              do.call(c,
                      alply (perm.lm, 1,
                             function (i) (1:n.traits) [which (i == landmarks)]))
          }
        else
          use.as.index <- sample (1:n.traits)
        H.perm.mat <- H.mat [use.as.index, use.as.index]
        H.perm.vec <- H.perm.mat [lower.tri (H.perm.mat)]
        avg.plus <- mean (cor.ss.vec [H.perm.vec == 1])
        avg.minus <- mean (cor.ss.vec [H.perm.vec != 1])
        (avg.plus - avg.minus) / icv
      }
    perm.dist <- aaply (1:n.it, 1, function (i) aaply (hyp.mat, 1, perm.mod.single))
    probs <-
      aaply (rbind (real.mod [, 3], perm.dist), 2,
             function (V) sum (V [1] < V [-1]) / length (V [-1]))
    out <- data.frame (colnames (mod.hyp), real.mod, probs)
    colnames (out) <- c('Hyp', 'AVG+', 'AVG-', 'AVG.Index', 'P')
    out
  }


RemoveCAC <- function (cov.ss)
  {
    logCS.pos <- which(rownames(cov.ss) == 'logCS')
    cov.shape <- cov.ss [- logCS.pos, - logCS.pos]
    CAC <- Normalize (cov.ss [logCS.pos, -logCS.pos] / cov.ss [logCS.pos, logCS.pos])
    Iaa <- diag (nrow (cov.shape)) - (CAC %*% t (CAC))
    cov.res <- t(Iaa) %*% cov.shape %*% Iaa
    cov.res
  }

ModIndex <- function (cov.matrix, mod.hyp, remove.size = FALSE, ...)
  {
    cor.matrix <- cov2cor (cov.matrix)
    if (remove.size)
      {
        cov.matrix <- RemoveSize (cov.matrix)
        cor.matrix <- RemoveSize (cor.matrix)
      }
    icv <- ICV (cov.matrix)
    mod.test <- TestModularity (cor.matrix, mod.hyp, ...)
    mod.index <- (mod.test [, 'AVG+'] - mod.test [, 'AVG-']) / icv
    out <- cbind (mod.test, mod.index)
    colnames (out) <- c (colnames (mod.test), 'ModIndex')
    return (out)
  }


