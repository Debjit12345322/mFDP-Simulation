
# Conformal outlier detection — Section 5.2 replication (R)

An R replication of the simulated outlier-detection experiment in Bates, Candès, Lei,
Romano & Sesia, *Testing for outliers with conformal p-values* (Section 5.2, Figure 7),
with the multiple-testing step swappable between Benjamini–Hochberg and mFDP.

The script fits a one-class SVM, turns its scores into conformal p-values, applies four
different calibrations to those p-values, and reports the false-discovery proportion and
power that each calibration achieves across many simulated "practitioners".

---

## Requirements

| Package | Role |
| --- | --- |
| `e1071` | one-class SVM (the scoring function) |
| `mFDP` | Hemerik's median-FDP adjustment |
| `ggplot2` | optional — base R boxplots are used if it's absent |
| `parallel` | optional — only if `PARALLEL <- TRUE` |

```r
install.packages(c("e1071", "ggplot2"))
# mFDP: install from its source repository
```

## Running it

```r
source("outlier_mfdp.R")
```

Everything runs top to bottom; no arguments, no cached state. `RNGkind("L'Ecuyer-CMRG")`
and `set.seed(2024)` make the run reproducible, including under `mclapply`.

At the default `J = 20`, `L = 20` this is a few minutes of work. The paper uses
`J = 100`, `L = 100`, which is roughly 25× the compute.

---

## What it does

**Data model.** Each observation is drawn from a 50-component Gaussian mixture in
`d = 50` dimensions:

```
X = sqrt(a) * V + W,   V ~ N(0, I_d),   W ~ Unif{w_1..w_50}
```

The mixture centres `Wset` are drawn once and stored with the results. Inliers use
`a = 1`; outliers use a larger `a`, so `a` is the signal strength and `a = 1` means *no
signal at all*.

**Per practitioner** (`run_practitioner`):

1. Draw `n_train + n_cal = 2000` points at `a = 1`; fit a radial one-class SVM on the
   first 1000 and score the other 1000 for calibration.
2. For each `a` in the grid, generate `L` test sets of 1000 points, 10% of which are
   outliers drawn at that `a`.
3. Compute marginal conformal p-values
   `u(X) = (1 + #{i : s(X_cal_i) <= s(X)}) / (n_cal + 1)`.
4. Apply each of the four calibrations, adjust for multiplicity, threshold at `alpha`,
   and record FDP and power.
5. Collapse the `L` test sets with `AGG` (default `median`; `mean` reproduces the paper's
   cFDR).

The training/calibration split is drawn at `a = 1` regardless of the test signal
strength, so the SVM and the calibration scores are fitted **once** per practitioner and
reused across the whole `a_grid`. Only test sets are regenerated.

---

## The four calibrations

| Name | What it does |
| --- | --- |
| `Marginal` | plain conformal p-values, valid on average over calibration sets |
| `Conditional (Simes)` | Simes-type simultaneous band |
| `Conditional (Asymptotic)` | DKW/normal-approximation band |
| `Conditional (Monte Carlo)` | hybrid band with a simulated finite-sample correction |

The three conditional variants are ports of `utils_calib.py` from the authors' Python
code. They inflate the marginal p-values so validity holds *conditionally* on the
calibration set, with confidence `1 - delta_ccv`.

Each is a confidence band `aseq` over the uniform order statistics; `betainv_generic`
maps a marginal p-value through it. Because `aseq` depends only on `(n_cal, delta, k)`
and never on the data, all three bands are built once at startup and reused — per test
set, the correction is a table lookup.

The Monte Carlo band needs `estimate_fs_correction`, a binary search over `gamma` that
simulates `n_mc = 10000` sorted uniform samples at each step. This is the expensive
setup step; it also runs only once, and the resulting `fs_correction` is saved with the
results.

---

## Two different deltas

Don't conflate them:

- **`delta_ccv`** — confidence level for the calibration-conditional correction (0.1).
- mFDP has **no** delta. It is a median guarantee, not a tail guarantee.

## BH vs mFDP

```
BH    guarantees  E[FDP] <= alpha          (mean)
mFDP  guarantees  P(FDP <= alpha) >= 0.5   (median)
```

The swap point is a single line in Section 5:

```r
##ADJUST <- function(p) mFDP.adjust(p, c = mfdp_c, s1 = mfdp_s1, s2 = mfdp_s2)
ADJUST <- function(p) p.adjust(p, "BH")   # to reproduce the paper
```

**As shipped, the script runs BH**, not mFDP — uncomment the first line and comment the
second to switch. Note that the plot titles and subtitles are hard-coded to say "mFDP"
either way, so they will mislabel a BH run.

`mFDP.adjust` takes no `alpha`, so one call serves the entire `alpha` grid. It can return
`Inf`, which should be read as 1; the `padj <= alpha` comparison handles that correctly.
`mfdp_s2` must cover the largest `alpha` you threshold at.

---

## Parameters

| Variable | Default | Meaning |
| --- | --- | --- |
| `d` | 50 | ambient dimension |
| `n_atom` | 50 | mixture components |
| `n_train` / `n_cal` / `n_test` | 1000 each | SVM fit / calibration / test set size |
| `prop_out` | 0.10 | outlier fraction per test set |
| `J` | 20 | practitioners (paper: 100) |
| `L` | 20 | test sets per practitioner (paper: 100) |
| `a_grid` | 1 → 3 by 0.25 | signal strength |
| `alpha_grid` | 0.1 | nominal level |
| `delta_ccv` | 0.1 | calibration-conditional confidence |
| `k_ccv` | `n_cal / 2` | Simes band parameter |
| `n_mc` | 10000 | MC replicates for the finite-sample correction |
| `AGG` | `median` | how the `L` test sets collapse |
| `mfdp_c` | `"1/(2m)"` | mFDP tuning; `"1/m"` is looser |
| `PARALLEL` | `FALSE` | use `mclapply` over practitioners |
| `OUTDIR` | `"mfdp_results"` | output root |

---

## Output

```
mfdp_results/
  all_results.csv / .rds     every (practitioner, a, alpha, method)
  summary.csv                collapsed over practitioners
  cfdr_by_method.png         faceted, one panel per calibration
  power_by_method.png
  cfdr_boxplot.png           calibrations dodged at each signal strength
  power_boxplot.png
  by_method/
    Marginal.csv, Marginal_cfdr.png, Marginal_power.png, ...
  by_a/
    a_1.00.csv, a_1.25.csv, ...
```

`all_results.rds` also stores the bands, the mixture centres, `fs_correction`, the full
parameter list, the seed, elapsed time, and `sessionInfo()` — enough to reconstruct the
run.

**Reading `summary.csv`:** `mFDR` is the *median* cFDR across practitioners, matching the
median guarantee. `frac_exceed` compares each row against its own `alpha`; under a median
guarantee you expect something near 0.5, and well below that means the procedure is
conservative.

---

## Notes and gotchas

- **SVM sign convention.** `e1071` returns decision values that are positive for
  inliers, so lower means more outlying. `check_orientation()` tests this empirically at
  startup and flips the sign if the local build disagrees. The flip corrects the
  software's convention, not the method.
- **`gamma = 1 / (d * var(as.vector(X_train)))`** mirrors scikit-learn's `gamma='scale'`.
- **p-value granularity.** The smallest attainable marginal p-value is
  `1 / (n_cal + 1) ≈ 9.99e-4`. With `alpha = 0.1` and 1000 tests this is a real
  constraint on power, not a rounding detail.
- **`a = 1` is a null column.** The "outlier" block is drawn at `a = 1` too, so every
  rejection is false and power is undefined (`NA`). The FDP boxplot at `a = 1` is the
  type-I error check.
- **Boxplots show spread across practitioners**, not across test sets — the `L` test sets
  were already collapsed by `AGG`. The centre line of each box is the `mFDR` reported in
  the summary table. Per-method PNGs share a y-range so they can be compared by eye.
- **`agg` in the saved RDS** is recorded via `deparse(substitute(AGG))` at top level,
  which yields the literal string `"AGG"` rather than `"median"`. Harmless, but don't
  trust that field to tell you which aggregator ran.
- **Sanity checks** sit in an `if (FALSE)` block at the end of the script — band lengths,
  monotonicity, conditional-dominates-marginal, and the smallest attainable calibrated
  p-value. Run them by hand after changing anything in Section 3.
