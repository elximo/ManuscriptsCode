#Name: Hazem Sharaf
#Location: Blacksburg, VA
#Date: November 15, 2026
#v0.0.2
#I wanted to get a p-value for the permuted location effect similar to nested.npmanova()
#I also changed function name. Hopefully this will go into use in other manuscripts or codes.
#I will continue to modify
#Output is still crude and will be further improved (i.e. format for significance
#Also seed updated to match the rest of the code

## Nested PERMANOVA with restricted permutations for UNBALANCED designs
## Works when farms have unequal numbers of cores (e.g. some have 4, most have 5).
##
## d     : dist object (e.g. vegdist(otu, "bray")), rows in the same order as meta.
## meta  : data.frame with columns Location and Soil_Amendment
#
#I used the one from QIIME2 as imported using qiime2R. More details on README
##
## Two tests, each with its own correct error term and permutation scheme:
##  1. Soil_Amendment: F = MS_amendment / MS_farms-within-amendment, df = (1, nFarms - 2).
##     Null distribution: amendment labels reassigned among whole farms. With
##     exact = TRUE, all choose(17, 8) = 24,310 assignments are enumerated (exact p).
##  2. Location (farms within amendment): F = MS_farms / MS_samples-within-farms,
##     df = (nFarms - 2, N - nFarms). Null distribution: cores shuffled among farms
##     that share the same amendment status (farm sizes kept), nperm_loc permutations.
##
## Both results match BiodiversityR::nested.npmanova(d ~ Soil_Amendment + Location),
## which uses the same tests with random permutations.


permanova_unbalanced <- function(d, meta, farm = "Location", trt = "Soil_Amendment",
                           exact = TRUE, nperm = 9999, nperm_loc = 9999, seed = 8112016) {
  stopifnot(inherits(d, "dist"), attr(d, "Size") == nrow(meta))
  ## reproducible permutations without disturbing the caller's random-number stream
  if (!is.null(seed)) {
    genv <- globalenv()
    old  <- if (exists(".Random.seed", envir = genv)) get(".Random.seed", envir = genv) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old, envir = genv),
            add = TRUE)
    set.seed(seed)
  }
  D2 <- as.matrix(d)^2
  f  <- factor(meta[[farm]])
  tr <- as.character(meta[[trt]])
  N  <- nrow(meta)

  ## farm-level treatment (must be constant within farm)
  tf <- tapply(tr, f, unique)
  if (is.list(tf)) stop("Treatment varies within a farm")
  lev <- sort(unique(tf))
  isA <- tf == lev[2]                   # TRUE for farms in the 2nd level (e.g. "Yes")
  nF  <- nlevels(f)
  m   <- as.vector(table(f))            # cores per farm

  ## nF x nF block sums of squared distances between farms
  Z <- model.matrix(~ f - 1)
  B <- t(Z) %*% D2 %*% Z

  SS_T      <- sum(D2) / (2 * N)
  SS_within <- sum(diag(B) / (2 * m))   # pooled within-farm SS (core-to-core)
  SS_farms  <- SS_T - SS_within         # between-farm SS

  df_trt <- 1
  df_loc <- nF - 2                      # farms within amendment
  df_res <- N - nF                      # cores within farms

  ## ---------- Test 1: Soil_Amendment, farms as the unit ----------
  ss_trt <- function(s) {               # s: logical vector over farms
    nA <- sum(m[s]); nB <- N - nA
    SS_T - (sum(B[s, s]) / (2 * nA) + sum(B[!s, !s]) / (2 * nB))
  }
  SS_A  <- ss_trt(isA)
  SS_L  <- SS_farms - SS_A
  F_trt <- (SS_A / df_trt) / (SS_L / df_loc)

  k <- sum(isA)
  if (exact) {
    combos   <- combn(nF, k)
    null_trt <- apply(combos, 2, function(idx) { s <- logical(nF); s[idx] <- TRUE; ss_trt(s) })
    n_trt    <- ncol(combos)
  } else {
    null_trt <- c(SS_A, replicate(nperm, ss_trt(sample(isA))))   # include observed
    n_trt    <- nperm
  }
  ## SS_farms is fixed under whole-farm shuffling, so F is monotone in SS_A
  p_trt <- mean(null_trt >= SS_A - 1e-12)

  ## ---------- Test 2: Location (farms within amendment), cores as the unit ----------
  ## Shuffling cores within amendment groups leaves SS_T and SS_A unchanged,
  ## so F_loc is monotone (decreasing) in the within-farm SS.
  ss_within_of <- function(lab) {
    sum(vapply(split(seq_len(N), lab),
               function(g) sum(D2[g, g]) / (2 * length(g)), numeric(1)))
  }
  F_loc  <- (SS_L / df_loc) / (SS_within / df_res)
  groups <- split(seq_len(N), tr)
  fc     <- as.character(f)
  null_w <- replicate(nperm_loc, {
    lab <- fc
    for (g in groups) lab[g] <- fc[g][sample(length(g))]
    ss_within_of(lab)
  })
  p_loc <- (1 + sum(null_w <= SS_within + 1e-12)) / (nperm_loc + 1)

  out <- data.frame(
    Df       = c(df_trt, df_loc, df_res, N - 1),
    SumOfSqs = c(SS_A, SS_L, SS_within, SS_T),
    R2       = c(SS_A, SS_L, SS_within, SS_T) / SS_T,
    F        = c(F_trt, F_loc, NA, NA),
    N.Perm   = c(n_trt, nperm_loc, NA, NA),
    `Pr(>F)` = c(p_trt, p_loc, NA, NA),
    row.names = c(trt, paste(farm, "(within", paste0(trt, ")")),
                  "Residual (cores within farm)", "Total"),
    check.names = FALSE
  )
  attr(out, "note") <- paste0(
    trt, " tested against ", farm, " (df ", df_trt, ", ", df_loc, "), ",
    if (exact) "exact p from all " else "p from ", n_trt, " farm-level assignments. ",
    farm, " tested against residual (df ", df_loc, ", ", df_res, "), ",
    nperm_loc, " permutations of cores within ", trt, " groups.")
  out
}
