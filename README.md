# cytokineR

<!-- badges: start -->
[![R-CMD-check](https://github.com/andjar/cytokineR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/andjar/cytokineR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Transparent, modular pre-processing of multiplex cytokine data** from
bead-based immunoassays (Bio-Rad Bio-Plex, Luminex xMAP): import, quality
control, plate and batch correction, fluorescence-to-concentration
calibration and imputation of censored values -- plus a benchmarking
framework to find out which combination works best for *your* data.

```r
# install.packages("remotes")
remotes::install_github("andjar/cytokineR")
```

## Why

Standard practice converts fluorescence intensity (FI) to concentration with
one five-parameter logistic curve per plate, then corrects batches somehow.
Those per-plate curves are fragile: few points, a flat low end where many
cytokines live, sensitivity to a single bad standard well, and plate-dependent
"out of range" censoring. Alternatives exist -- working on FI directly,
multiples of blank, standard or control, curves shared across plates -- but
it is rarely clear which is best for a given assay and layout.

cytokineR makes every step a small, swappable, logged method, and provides
the tools to compare them honestly.

## Quick start

```r
library(cytokineR)

x <- read_bioplex(c("plate1.xlsx", "plate2.xlsx"), batch = "lot_A")   # or read_xponent()
x <- x |>
  flag_beads(min_count = 35) |>
  cyto_calibrate("shared") |>          # shared-shape 5PL across plates
  cyto_impute("rank_tobit") |>         # censored values, ordered by their FI
  cyto_correct_batch("bridge")         # bridging samples between lots

cyto_steps(x)        # what was done, with which arguments
cyto_fit(x, "calibrate")$curves   # what was estimated
w <- cyto_wide(x)    # one row per sample, one column per analyte
```

## What is inside

| stage | function | methods |
|---|---|---|
| import | `read_bioplex()`, `read_xponent()`, `cyto_from_long()`, `cyto_from_wide()` | Bio-Plex Manager multi-analyte and per-analyte exports (Excel/CSV), xPONENT CSV, any table |
| QC | `flag_beads()`, `flag_replicates()`, `qc_pca()`, `qc_plates()`, `qc_controls()`, `qc_variance()` | bead counts, replicate CVs, PCA by plate/batch, plate z-scores, Levey-Jennings + Westgard rules, variance decomposition |
| signal | `cyto_transform()` | `log2`, `bgsub`, `mob` (multiple of blank, EB-shrunk), `mos` (multiple of standard or reference concentration), `moc` (multiple of control), `vendor` |
| plate / batch | `cyto_normalize()`, `cyto_correct_batch()` | `median` (within strata), `quantile`, `anchor` (offset or affine; controls, standards), `bridge` (offset or Deming), `regress` (covariate-protected), `combat` |
| calibration | `cyto_calibrate()` | `5pl`, `4pl` (optionally robust), `spline` (monotone, no extrapolation), `shared` (robust shared-shape 5PL with shrunken gain/potency/background shifts); per plate, per batch, or one curve per study |
| censored values | `cyto_impute()` | `fraction`, `tobit` (covariates, draws), `rank_tobit` |
| comparison | `cyto_pipeline()`, `expand_pipelines()`, `benchmark_pipelines()`, `benchmark_simulations()`, `simulate_cyto()` | hold-out benchmarking with scale-free metrics; simulator with known truth |
| attribution | `waterfall_lattice()`, `metric_effect()` | every on/off combination of the steps; each step's waterfall change, order-free Shapley share and range across contexts |

New methods are registered with `register_cyto_method()` and are immediately
usable in pipelines and benchmarks (`vignette("extending")`).

## Finding better methods

```r
grid <- expand_pipelines(
  signal = list(log2 = pipe_transform("log2"), MoB = pipe_transform("mob")),
  plate  = list(none = NULL, median = pipe_normalize("median", within = "time"),
                ctrl = pipe_normalize("anchor", anchors = "pool_A", scope = "batch")),
  batch  = list(bridge = pipe_correct_batch("bridge"))
)
bm <- benchmark_pipelines(x_raw, grid,
                          holdout = c(list(pool = "pool_B"),
                                      holdout_folds(x_raw, "bridge", k = 2)),
                          metrics = default_metrics(bio = ~ time, adjust = ~ subject))
bm; plot(bm)
```

* **Hold-out:** material used to judge a pipeline (a control pool, half of the
  bridging samples) is hidden from every fitting step, so no method is graded
  on what it was fitted to.
* **Scale-free metrics:** technical reliability (an intra-class correlation of
  repeated material), residual plate R², retained biological signal,
  detection rate, and -- in simulations -- agreement with the truth. Pipelines
  producing log2 FI, ratios or concentrations are compared on equal terms.
* **Simulation with known truth:** `simulate_cyto()` models plate gain,
  potency and background effects, kit-lot errors in the standards, bad
  standard wells, low-abundance analytes and confounded layouts;
  `benchmark_simulations()` runs pipelines across such scenarios.

## Which step did it?

A table of results "after each step", or a waterfall chart, credits each step
with the change it caused *in that order*. Steps that partly do the same job
split the credit by order: the first looks decisive, the next one useless.
`waterfall_lattice()` runs every on/off combination of the optional steps
(2^K pipelines) and computes a statistic for each:

```r
steps <- list(
  log2   = pipe_transform("log2"),
  median = pipe_normalize("median", within = "time"),
  ctrl   = pipe_normalize("anchor", anchors = "ctrl_low", scope = "batch"),
  bridge = pipe_correct_batch("bridge")
)
lat <- waterfall_lattice(x_raw, steps, fixed = "log2",
                         statistic = default_metrics()$reliability_batch,
                         holdout = "ctrl_high")
lat                        # waterfall change vs Shapley share per step
plot(lat)                  # every pipeline a node, the waterfall route in bold
plot(lat, "contributions") # each step's change in every context
```

In the example data, the waterfall credits median plate normalisation with
almost all of the gain in batch-to-batch reliability and bridging correction
with almost nothing. Their order-free Shapley shares are about equal, because
centring every plate also removes most of the batch shift. Any statistic
works, for example the effect you will report
(`metric_effect(x, ~ time, adjust = ~ subject)`, with p-values marking which
pipelines reach significance). See `vignette("waterfall-lattice")`.

## What the simulations suggest

Eleven representative pipelines, six scenarios (default, strong plate
effects, confounded layout, 15% bad standard wells, low abundance, noisy
blanks), five simulated studies each
(`inst/scripts/simulation_study.R`). Mean rank over the six scenarios
(1 = best), and two of the metrics averaged over scenarios:

| pipeline | mean rank | worst scenario rank | truth ρ | reliability (plate) |
|---|---|---|---|---|
| shared 5PL > rank-Tobit > bridge | 2.8 | 3.1 | 0.92 | 0.94 |
| log2 FI > sample median (within time) > bridge | 3.2 | 3.6 | 0.91 | 0.92 |
| harmonise on standards > one curve (reference lot) > rank-Tobit | 3.6 | 4.0 | 0.91 | 0.94 |
| multiple of control > bridge | 3.9 | 4.9 | 0.91 | 0.90 |
| log2 FI > standards anchor (affine) > bridge | 4.1 | 4.8 | 0.90 | 0.92 |
| multiple of standard (ref. conc.) > bridge | 6.2 | 7.6 | 0.89 | 0.86 |
| spline per plate > Tobit > bridge | 7.1 | 8.1 | 0.84 | 0.76 |
| multiple of blank > bridge | 7.3 | 9.1 | 0.84 | 0.76 |
| 5PL per plate > Tobit > bridge | 8.1 | 8.6 | 0.80 | 0.70 |
| log2 FI, uncorrected | 8.5 | 9.8 | 0.81 | 0.71 |
| 5PL per plate > Tobit | 8.7 | 9.4 | 0.80 | 0.70 |

Take-aways, *under the simulator's model*:

* Per-plate 5PL curves were near the bottom in every scenario; with 15% bad
  standard wells their reliability dropped to 0.65. Sharing the curve shape
  across the plates of a kit lot (robust fit, shrunken per-plate shifts) was
  the most accurate route to concentrations (bias-corrected RMSE 0.30--0.51
  log2 vs 0.45--0.76 for per-plate 5PL).
* Without any curve, simple FI-based routes (sample-median centring within
  time point, control anchoring) came close -- attractive for labs that only
  need relative values.
* Multiple-of-blank was weak whenever blanks do not track the plate gain or
  are few and noisy; it remains useful as a detectability scale.
* A plate correction applied before *per-plate* calibration is undone by it;
  harmonise first and then calibrate once (`curve_by = NULL`).

These are hypotheses to test, not conclusions: the simulator encodes one
model of plate effects. Run `benchmark_pipelines()` with hold-out material on
real data before choosing.

## Vignettes

* `vignette("cytokineR")` -- getting started
* `vignette("signal-scales")` -- plate effects, signal representations, why
  standard curves are fragile, harmonise-then-calibrate
* `vignette("benchmarking")` -- comparing pipelines on your own data
* `vignette("waterfall-lattice")` -- how much each step contributes, in
  every combination of the other steps
* `vignette("extending")` -- your own methods and metrics

## Status

Version 0.1.0. The importers have been tested on synthetic files that follow
the documented Bio-Plex Manager and xPONENT layouts; exports from other
software versions may need small adjustments -- please open an issue with an
(anonymised) example file.

## License

GPL (>= 3)
