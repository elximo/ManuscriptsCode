## PERMANOVA with whole-farm (restricted) permutation for UNBALANCED designs
## Works when farms have unequal numbers of cores (e.g. some have 4, most have 5).
##
## d     : dist object (e.g. vegdist(otu, "bray")), rows in the same order as meta
## meta  : data.frame with columns Location and Soil_Amendment
## Tests Soil_Amendment using farms as the experimental unit, by enumerating
## (or sampling) all reassignments of the amendment labels among farms.

#Name: Hazem Sharaf
#Location: Blacksburg, VA
#Date: November 14, 2026
#v0.0.1
#This initial version was generated using Claude, Opus model in High effort mode.
#nested.npmanova was very slow as it has to make a call to adonis2(). This is a much faster version on 2 orders of magnitude improvement.

farm_permanova <- function(d, meta, farm = "Location", trt = "Soil_Amendment",
                           exact = TRUE, nperm = 9999, seed = 42) {
  stopifnot(inherits(d, "dist"), attr(d, "Size") == nrow(meta))
  D2   <- as.matrix(d)^2
  f    <- factor(meta[[farm]])
  N    <- nrow(meta)

  ## farm-level treatment (must be constant within farm)
  tf <- tapply(as.character(meta[[trt]]), f, unique)
  if (is.list(tf)) stop("Treatment varies within a farm")
  lev   <- sort(unique(tf))
  isA   <- tf == lev[2]                 # TRUE for farms in the 2nd level (e.g. "Yes")
  nF    <- nlevels(f)
  m     <- as.vector(table(f))          # cores per farm

  ## 17 x 17 block sums of squared distances between farms
  Z <- model.matrix(~ f - 1)
  B <- t(Z) %*% D2 %*% Z

  SS_T      <- sum(D2) / (2 * N)
  SS_within <- sum(diag(B) / (2 * m))   # pooled within-farm SS (core-to-core)
  SS_farms  <- SS_T - SS_within         # between-farm SS (invariant under permutation)

  ## treatment SS for any farm-level labelling (logical vector over farms)
  ss_trt <- function(s) {
    nA <- sum(m[s]); nB <- N - nA
    SS_W <- sum(B[s, s]) / (2 * nA) + sum(B[!s, !s]) / (2 * nB)
    SS_T - SS_W
  }

  SS_A   <- ss_trt(isA)
  df_trt <- 1
  df_err <- nF - 2                      # farms within treatment
  F_obs  <- (SS_A / df_trt) / ((SS_farms - SS_A) / df_err)

  ## null distribution: reassign the k "Yes" labels among the nF farms
  k <- sum(isA)
  if (exact) {
    combos <- combn(nF, k)              # choose(17, 8) = 24,310 assignments
    null <- apply(combos, 2, function(idx) { s <- logical(nF); s[idx] <- TRUE; ss_trt(s) })
  } else {
    set.seed(seed)
    null <- replicate(nperm, ss_trt(sample(isA)))
    null <- c(SS_A, null)               # include observed
  }
  ## F is monotone in SS_A (SS_farms fixed), so compare SS directly
  p <- mean(null >= SS_A - 1e-12)

  data.frame(
    term   = c(trt, paste0(farm, " (error)"), "Residual (within farm)", "Total"),
    Df     = c(df_trt, df_err, N - nF, N - 1),
    SumOfSqs = c(SS_A, SS_farms - SS_A, SS_within, SS_T),
    R2     = c(SS_A, SS_farms - SS_A, SS_within, SS_T) / SS_T,
    F      = c(F_obs, NA, NA, NA),
    `Pr(>F)` = c(p, NA, NA, NA),
    check.names = FALSE
  )
}
