## ==========================================================================
##  Bates, Candes, Lei, Romano, Sesia
##  "Testing for outliers with conformal p-values", Section 5.2
##  Outlier detection on simulated data -- replication in R.
##
##  FOUR p-value calibrations, as in Figure 7:
##      Marginal                 plain conformal p-values
##      Conditional (Simes)      betainv_simes
##      Conditional (Asymptotic) betainv_asymptotic
##      Conditional (Monte Carlo) betainv_mc
##  The three conditional ones are ports of utils_calib.py from the authors'
##  Python code. They inflate the marginal p-values so that validity holds
##  CONDITIONALLY on the calibration set, with confidence 1 - delta_ccv.
##
##  Multiple-testing step: mFDP (Hemerik) instead of BH.
##      BH   guarantees  E[FDP] <= alpha         (mean)
##      mFDP guarantees  P(FDP <= alpha) >= 0.5  (median)
##
##  Two distinct uses of "delta" -- do not conflate:
##      delta_ccv   confidence for the calibration-conditional correction
##      (mFDP has no delta; it is a median guarantee)
##
##  Output tree (section 8):
##      mfdp_results/
##        all_results.csv / .rds     every (practitioner, a, alpha, method)
##        summary.csv                collapsed over practitioners
##        by_method/Marginal.csv ...
##        by_a/a_1.00.csv ...
## ==========================================================================

suppressPackageStartupMessages({
  library(e1071)
  library(mFDP)
})

RNGkind("L'Ecuyer-CMRG")
set.seed(2024)

## --------------------------------------------------------------------------
## 1.  Parameters
## --------------------------------------------------------------------------

d        <- 50        # ambient dimension
n_atom   <- 50        # |W| : number of mixture components
n_train  <- 1000      # fit the one-class SVM
n_cal    <- 1000      # calibrate conformal p-values
n_test   <- 1000      # size of each test set
prop_out <- 0.10      # fraction of outliers in each test set

## Paper: J = 100, L = 100.
J <- 20
L <- 20

a_grid     <- c(1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3)
alpha_grid <- 0.1

## Calibration-conditional parameters (paper Section 5.2.2: delta = 0.1,
## k = n_cal/2).
delta_ccv  <- 0.1
simes_kden <- 2
k_ccv      <- as.integer(n_cal / simes_kden)
n_mc       <- 10000   # MC replicates for the finite-sample correction

## How the L test sets are collapsed for each practitioner.
AGG <- median         # `mean` reproduces the paper's cFDR

## mFDP.adjust tuning. s2 must cover the largest alpha thresholded at.
mfdp_c  <- "1/(2m)"   # "1/m" is looser; try it if power looks low
mfdp_s1 <- 0
mfdp_s2 <- max(alpha_grid)

PARALLEL <- FALSE
n_cores  <- max(1L, parallel::detectCores() - 1L)

OUTDIR <- "mfdp_results"

## --------------------------------------------------------------------------
## 2.  Data generating model  P_X^a
##
##     X_i = sqrt(a) * V_i + W_i,   V_i ~ N(0, I_d),   W_i ~ Unif(W-set)
##     p(x) = (1/50) sum_k N(x | w^(k), a * I_d)
## --------------------------------------------------------------------------

Wset <- matrix(runif(n_atom * d, min = -3, max = 3), nrow = n_atom, ncol = d)

rP <- function(n, a) {
  V <- matrix(rnorm(n * d), n, d)
  W <- Wset[sample.int(n_atom, n, replace = TRUE), , drop = FALSE]
  sqrt(a) * V + W
}

## --------------------------------------------------------------------------
## 3.  Calibration-conditional machinery -- port of utils_calib.py
##
##     Each method builds a simultaneous confidence band `aseq` for the
##     uniform order statistics; betainv_generic then maps a marginal
##     p-value through it. Crucially aseq depends ONLY on (n_cal, delta,
##     k), never on the data, so all three are built once here and reused
##     for every test set. Per test set the correction is a lookup.
## --------------------------------------------------------------------------

compute_cn <- function(delta, n) {
  cn <- -log(-log(1 - delta)) + 2 * log(log(n)) +
    0.5 * log(log(log(n))) - 0.5 * log(pi)
  cn / sqrt(2 * log(log(n)))
}

## rolling mean of width k, "valid" mode -> length n - k + 1
moving_average <- function(x, k) {
  cs <- cumsum(c(0, x))
  (cs[(k + 1):(length(x) + 1)] - cs[1:(length(x) - k + 1)]) / k
}

compute_aseq <- function(n, k, delta) {
  k    <- as.integer(k)
  fac1 <- log(delta) / k - mean(log((n - k + 1):n))
  fac2 <- moving_average(log(1:n), k)
  c(rep(0, k - 1), exp(fac2 + fac1))
}

## maps marginal p-values through a band; aseq must have length n_cal
betainv_generic <- function(pvals, aseq) {
  n   <- length(aseq)
  idx <- pmax(1L, pmin(n, as.integer(floor((n + 1) * (1 - pvals)))))
  1 - aseq[idx]
}

aseq_simes <- function(n, k, delta) compute_aseq(n, k, delta)

aseq_asymptotic <- function(n, delta) {
  iseq <- 1:n
  cn   <- compute_cn(delta, n)
  b    <- iseq / n + cn * sqrt(iseq * (n - iseq)) / (n * sqrt(n))
  1 - pmin(1, rev(b))
}

compute_hybrid_bound <- function(delta, n, gamma) {
  i     <- 1:n
  cna   <- compute_cn(delta - gamma, n)
  bound <- i / n + cna * sqrt(i * (n - i)) / (n * sqrt(n))
  
  ## linearise past n/2
  kl    <- as.integer(n / 2)
  slope <- bound[kl] - bound[kl - 1]
  idx   <- (kl + 1):n
  bound[idx] <- bound[kl] + slope * (i[idx] - kl)
  
  ## and cap by the Simes band
  bound_s <- 1 - rev(compute_aseq(n, as.integer(n / 2), delta))
  pmin(bound_s, bound)
}

## Binary search for the finite-sample correction gamma: the value making
## the simulated probability of the empirical process crossing the hybrid
## bound equal to delta. This is the expensive step -- run once.
estimate_fs_correction <- function(delta, n, n_mc = 10000) {
  ## columns are sorted uniform samples, so `Ut > bound` recycles bound
  ## down each column (column-major) without materialising a big matrix
  Ut <- matrix(runif(n_mc * n), nrow = n, ncol = n_mc)
  Ut <- apply(Ut, 2, sort)
  
  prob_crossing <- function(gamma) {
    bound <- compute_hybrid_bound(delta, n, gamma)
    mean(colSums(Ut > bound) > 0)
  }
  
  gamma0 <- -(1 - 1e-6 - delta)
  gamma1 <- delta - 1e-6
  gamma  <- gamma1
  while (abs(gamma1 - gamma0) > 1e-6) {
    gamma <- (gamma0 + gamma1) / 2
    if (prob_crossing(gamma) > delta) gamma0 <- gamma else gamma1 <- gamma
  }
  gamma
}

cat("building calibration-conditional bands ...\n")
t_band <- Sys.time()
fs_correction <- estimate_fs_correction(delta_ccv, n_cal, n_mc)
ASEQ <- list(
  "Conditional (Simes)"       = aseq_simes(n_cal, k_ccv, delta_ccv),
  "Conditional (Asymptotic)"  = aseq_asymptotic(n_cal, delta_ccv),
  "Conditional (Monte Carlo)" = 1 - pmin(1, rev(compute_hybrid_bound(delta_ccv, n_cal, fs_correction)))
)
cat(sprintf("  fs_correction = %.6f   (%.1f s)\n", fs_correction,
            as.numeric(difftime(Sys.time(), t_band, units = "secs"))))

METHODS <- c("Marginal", names(ASEQ))
n_meth  <- length(METHODS)

## applies the calibration named by `m` to marginal p-values
calibrate <- function(p, m) {
  if (m == "Marginal") p else pmin(1, betainv_generic(p, ASEQ[[m]]))
}

## --------------------------------------------------------------------------
## 4.  Scoring function and orientation
##
##     e1071 one-class SVM: f(x) positive for inliers, so LOWER = more
##     outlying. Checked rather than assumed; the flip fixes the software's
##     sign convention, not the method.
## --------------------------------------------------------------------------

score <- function(fit, X) {
  as.numeric(attr(predict(fit, X, decision.values = TRUE), "decision.values"))
}

check_orientation <- function() {
  Xtr <- rP(n_train, a = 1)
  f <- svm(Xtr, type = "one-classification", kernel = "radial",
           nu = 0.05, gamma = 1 / (d * var(as.vector(Xtr))), scale = FALSE)
  if (mean(score(f, rP(2000, a = 3))) > mean(score(f, rP(2000, a = 1)))) {
    message("Decision values inverted on this build; flipping sign.")
    return(-1)
  }
  1
}
SIGN <- check_orientation()

## --------------------------------------------------------------------------
## 5.  Conformal p-values and the multiple-testing step
##
##     u_hat(X) = (1 + #{i : s(X_cal_i) <= s(X)}) / (n_cal + 1)
##     Granularity floor: min p = 1/(n_cal + 1) ~ 9.99e-4.
##
##     mFDP.adjust can return Inf, to be read as 1; `padj <= alpha` handles
##     that correctly.
## --------------------------------------------------------------------------

conformal_p <- function(s_test, s_cal_sorted) {
  (1 + findInterval(s_test, s_cal_sorted)) / (length(s_cal_sorted) + 1)
}

## SWAP POINT. mFDP.adjust(P, c, s1, s2) takes no alpha, so adjusted
## p-values serve the whole alpha grid from one call.
##ADJUST <- function(p) mFDP.adjust(p, c = mfdp_c, s1 = mfdp_s1, s2 = mfdp_s2)
ADJUST <- function(p) p.adjust(p, "BH")   # to reproduce the paper

## --------------------------------------------------------------------------
## 6.  One practitioner
##
##     D_j is drawn with a = 1 regardless of test-set signal strength, so
##     the SVM and calibration scores are fitted once and reused across the
##     whole a_grid. Only test sets are regenerated.
## --------------------------------------------------------------------------

run_practitioner <- function(j) {
  
  D       <- rP(n_train + n_cal, a = 1)
  X_train <- D[seq_len(n_train), , drop = FALSE]
  X_cal   <- D[n_train + seq_len(n_cal), , drop = FALSE]
  
  fit <- svm(X_train,
             type   = "one-classification",
             kernel = "radial",
             nu     = 0.05,
             gamma  = 1 / (d * var(as.vector(X_train))),   # sklearn gamma='scale'
             scale  = FALSE)
  
  s_cal_sorted <- sort(SIGN * score(fit, X_cal))
  
  n_out   <- round(prop_out * n_test)
  n_in    <- n_test - n_out
  out_idx <- n_in + seq_len(n_out)
  n_al    <- length(alpha_grid)
  
  out <- vector("list", length(a_grid))
  
  for (ai in seq_along(a_grid)) {
    a <- a_grid[ai]
    
    ## all L test sets for this signal strength, scored in one call
    Xt <- rbind(rP(n_in * L, a = 1), rP(n_out * L, a = a))
    st <- SIGN * score(fit, Xt)
    S_in  <- matrix(st[seq_len(n_in * L)],  nrow = n_in,  ncol = L)
    S_out <- matrix(st[-seq_len(n_in * L)], nrow = n_out, ncol = L)
    
    ## at a = 1 the "outlier" block is drawn from a = 1 too: no signal,
    ## every rejection is false, power undefined
    truth <- if (a > 1) out_idx else integer(0)
    
    fdp <- array(NA_real_, c(L, n_al, n_meth))
    pwr <- array(NA_real_, c(L, n_al, n_meth))
    
    for (l in seq_len(L)) {
      p_marg <- conformal_p(c(S_in[, l], S_out[, l]), s_cal_sorted)
      
      for (mi in seq_len(n_meth)) {
        padj <- ADJUST(calibrate(p_marg, METHODS[mi]))
        
        for (k in seq_len(n_al)) {
          R <- which(padj <= alpha_grid[k])
          fdp[l, k, mi] <- if (length(R)) length(setdiff(R, truth)) / length(R) else 0
          pwr[l, k, mi] <- if (length(truth)) length(intersect(R, truth)) / length(truth) else NA_real_
        }
      }
    }
    
    ## collapse the L test sets -> one value per (alpha, method)
    grid <- expand.grid(alpha = alpha_grid, method = METHODS,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    out[[ai]] <- data.frame(practitioner = j,
                            a            = a,
                            alpha        = grid$alpha,
                            method       = grid$method,
                            cFDR         = as.vector(apply(fdp, c(2, 3), AGG)),
                            cPower       = as.vector(apply(pwr, c(2, 3), AGG)))
  }
  
  do.call(rbind, out)
}

## --------------------------------------------------------------------------
## 7.  Run -- every practitioner recomputed on every run
## --------------------------------------------------------------------------

t0 <- Sys.time()
if (PARALLEL) {
  res <- parallel::mclapply(seq_len(J), run_practitioner, mc.cores = n_cores)
} else {
  res <- lapply(seq_len(J), function(j) {
    cat(sprintf("practitioner %3d / %d\n", j, J)); flush.console()
    run_practitioner(j)
  })
}
results <- do.call(rbind, res)
results$method <- factor(results$method, levels = METHODS)
results <- results[order(results$method, results$a,
                         results$alpha, results$practitioner), ]
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nelapsed: %.1f min   (rows: %d, expected %d)\n", elapsed,
            nrow(results), J * length(a_grid) * length(alpha_grid) * n_meth))

## --------------------------------------------------------------------------
## 8.  Summary
##
##     mFDR is the MEDIAN of cFDR across practitioners (matching the median
##     guarantee). frac_exceed compares each row against ITS OWN alpha;
##     under a median guarantee, near 0.5 is expected, well below means the
##     procedure is conservative.
## --------------------------------------------------------------------------

key  <- interaction(results$method, results$a, results$alpha, drop = TRUE)
summ <- do.call(rbind, lapply(split(results, key), function(g) data.frame(
  method      = g$method[1],
  a           = g$a[1],
  alpha       = g$alpha[1],
  mFDR        = median(g$cFDR),
  cFDR_q90    = quantile(g$cFDR, 0.90, names = FALSE),
  frac_exceed = mean(g$cFDR > g$alpha[1]),
  mPower      = median(g$cPower)
)))
summ <- summ[order(summ$method, summ$a, summ$alpha), ]
print(summ, row.names = FALSE, digits = 3)

## --------------------------------------------------------------------------
## 9.  Save
## --------------------------------------------------------------------------

dir.create(file.path(OUTDIR, "by_method"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "by_a"),      recursive = TRUE, showWarnings = FALSE)

write.csv(results, file.path(OUTDIR, "all_results.csv"), row.names = FALSE)
write.csv(summ,    file.path(OUTDIR, "summary.csv"),     row.names = FALSE)

saveRDS(list(results = results, summary = summ, aseq = ASEQ,
             params = list(J = J, L = L, d = d, n_train = n_train,
                           n_cal = n_cal, n_test = n_test,
                           prop_out = prop_out, a_grid = a_grid,
                           alpha_grid = alpha_grid, methods = METHODS,
                           delta_ccv = delta_ccv, k_ccv = k_ccv,
                           n_mc = n_mc, fs_correction = fs_correction,
                           mfdp_c = mfdp_c, mfdp_s1 = mfdp_s1, mfdp_s2 = mfdp_s2,
                           agg = deparse(substitute(AGG)),
                           Wset = Wset, seed = 2024, sign = SIGN,
                           elapsed_min = elapsed,
                           R = R.version.string,
                           sessionInfo = sessionInfo())),
        file.path(OUTDIR, "all_results.rds"))

for (m in METHODS) {
  g <- results[results$method == m, c("practitioner", "a", "alpha", "cFDR", "cPower")]
  write.csv(g, file.path(OUTDIR, "by_method",
                         sprintf("%s.csv", gsub("[^A-Za-z]+", "_", m))),
            row.names = FALSE)
}

for (aa in a_grid) {
  g <- results[results$a == aa, c("practitioner", "method", "alpha", "cFDR", "cPower")]
  write.csv(g, file.path(OUTDIR, "by_a", sprintf("a_%.2f.csv", aa)), row.names = FALSE)
}

cat(sprintf("\nwritten to %s\n", normalizePath(OUTDIR)))
print(list.files(OUTDIR, recursive = TRUE))

## --------------------------------------------------------------------------
## 10.  Plots
##
##      Distribution across the J practitioners, one box per signal
##      strength. The box's centre line is the median across practitioners,
##      i.e. the mFDR reported in section 8.
##
##      Two PNGs per method, written to OUTDIR/by_method/:
##          <Method>_cfdr.png
##          <Method>_power.png
##      plus combined faceted views at the top level. All per-method plots
##      share a y-range so they can be compared by eye.
## --------------------------------------------------------------------------

al       <- alpha_grid[1]
ylim_fdr <- range(results$cFDR,   na.rm = TRUE)
ylim_pwr <- range(results$cPower, na.rm = TRUE)
slug     <- function(m) gsub("[^A-Za-z]+", "_", m)

dir.create(file.path(OUTDIR, "by_method"), recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  ## ---- one pair of PNGs per method ---------------------------------------
  for (m in METHODS) {
    g <- results[results$method == m, ]
    
    pm_fdr <- ggplot(g, aes(factor(a), cFDR)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      geom_hline(yintercept = al, linetype = "dashed", colour = "grey40") +
      coord_cartesian(ylim = ylim_fdr) +
      labs(x = "Signal strength", y = "cFDR",
           title = sprintf("%s -- cFDP across %d practitioners", m, J),
           subtitle = sprintf("mFDP, alpha = %.2f", al)) +
      theme_bw()
    
    pm_pwr <- ggplot(g, aes(factor(a), cPower)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      coord_cartesian(ylim = ylim_pwr) +
      labs(x = "Signal strength", y = "cPower",
           title = sprintf("%s -- power across %d practitioners", m, J),
           subtitle = sprintf("mFDP, alpha = %.2f", al)) +
      theme_bw()
    
    ggsave(file.path(OUTDIR, "by_method", sprintf("%s_cfdr.png",  slug(m))),
           pm_fdr, width = 6.5, height = 4.5, dpi = 150)
    ggsave(file.path(OUTDIR, "by_method", sprintf("%s_power.png", slug(m))),
           pm_pwr, width = 6.5, height = 4.5, dpi = 150)
  }
  
  ## ---- combined: one panel per method -------------------------------------
  p_fdr <- ggplot(results, aes(factor(a), cFDR)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_hline(yintercept = al, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ method, nrow = 1) +
    labs(x = "Signal strength", y = "cFDR",
         title = sprintf("cFDP by calibration (mFDP, alpha = %.2f)", al)) +
    theme_bw()
  
  p_pwr <- ggplot(results, aes(factor(a), cPower)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    facet_wrap(~ method, nrow = 1) +
    labs(x = "Signal strength", y = "cPower",
         title = sprintf("Power by calibration (mFDP, alpha = %.2f)", al)) +
    theme_bw()
  
  ## ---- combined: methods dodged at each signal strength -------------------
  p_fdr_d <- ggplot(results, aes(factor(a), cFDR, fill = method)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    geom_hline(yintercept = al, linetype = "dashed", colour = "grey40") +
    labs(x = "Signal strength", y = "cFDR", fill = "Calibration",
         title = sprintf("cFDP across %d practitioners (mFDP, alpha = %.2f)", J, al)) +
    theme_bw()
  
  p_pwr_d <- ggplot(results, aes(factor(a), cPower, fill = method)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8)) +
    labs(x = "Signal strength", y = "cPower", fill = "Calibration",
         title = sprintf("Power across %d practitioners (mFDP, alpha = %.2f)", J, al)) +
    theme_bw()
  
  ggsave(file.path(OUTDIR, "cfdr_by_method.png"),  p_fdr,   width = 13, height = 4,   dpi = 150)
  ggsave(file.path(OUTDIR, "power_by_method.png"), p_pwr,   width = 13, height = 4,   dpi = 150)
  ggsave(file.path(OUTDIR, "cfdr_boxplot.png"),    p_fdr_d, width = 10, height = 4.5, dpi = 150)
  ggsave(file.path(OUTDIR, "power_boxplot.png"),   p_pwr_d, width = 10, height = 4.5, dpi = 150)
  
  print(p_fdr); print(p_pwr)
  
} else {
  
  ## ---- base R fallback: same files, no ggplot2 ----------------------------
  for (m in METHODS) {
    g <- results[results$method == m, ]
    
    png(file.path(OUTDIR, "by_method", sprintf("%s_cfdr.png", slug(m))),
        width = 1000, height = 750, res = 150)
    boxplot(cFDR ~ a, g, ylim = ylim_fdr, main = sprintf("%s -- cFDP", m),
            xlab = "Signal strength", ylab = "cFDR")
    abline(h = al, lty = 2)
    dev.off()
    
    png(file.path(OUTDIR, "by_method", sprintf("%s_power.png", slug(m))),
        width = 1000, height = 750, res = 150)
    boxplot(cPower ~ a, g, ylim = ylim_pwr, main = sprintf("%s -- power", m),
            xlab = "Signal strength", ylab = "cPower")
    dev.off()
  }
  
  png(file.path(OUTDIR, "cfdr_by_method.png"), width = 1800, height = 500, res = 130)
  op <- par(mfrow = c(1, n_meth), mar = c(4, 4, 3, 1))
  for (m in METHODS) {
    boxplot(cFDR ~ a, results[results$method == m, ], ylim = ylim_fdr,
            main = m, xlab = "Signal strength", ylab = "cFDR")
    abline(h = al, lty = 2)
  }
  par(op); dev.off()
  
  png(file.path(OUTDIR, "power_by_method.png"), width = 1800, height = 500, res = 130)
  op <- par(mfrow = c(1, n_meth), mar = c(4, 4, 3, 1))
  for (m in METHODS) {
    boxplot(cPower ~ a, results[results$method == m, ], ylim = ylim_pwr,
            main = m, xlab = "Signal strength", ylab = "cPower")
  }
  par(op); dev.off()
}

cat("\nplots written:\n")
print(list.files(OUTDIR, pattern = "\\.png$", recursive = TRUE))

## --------------------------------------------------------------------------
## 11.  Sanity checks -- run these once by hand
## --------------------------------------------------------------------------
if (FALSE) {
  ## each band must have length n_cal and be non-decreasing
  stopifnot(all(sapply(ASEQ, length) == n_cal))
  sapply(ASEQ, function(v) all(diff(v) >= -1e-12))
  
  ## conditional p-values must dominate the marginal ones everywhere
  pt <- seq(1 / (n_cal + 1), 1, length.out = 500)
  sapply(names(ASEQ), function(m) all(calibrate(pt, m) >= pt - 1e-12))
  
  ## and the smallest attainable value under each calibration
  sapply(METHODS, function(m) min(calibrate(1 / (n_cal + 1), m)))
}

