# The `tsd` Discrepancy Between pmsimstats-orig and pmsimstats2025

## Summary

The `tsd` (time since discontinuation) variable in
`carryover_power_simulation.R` (pmsimstats2025) is computed differently
from the original `buildtrialdesign()` in pmsimstats-orig. The 2025
code produces `tsd = 0` at every first off-drug measurement, while
the original code produces `tsd >= 1` (the interval width leading into
that measurement). This discrepancy is the sole cause of the
$t_{\text{sd}} = 0$ discontinuity documented in
`predictor_comparison.tex`. The discontinuity does not exist in
Hendrickson's original implementation.

## The Two Algorithms

### Original: `buildtrialdesign()` (pmsimstats-orig)

**File:** `01c-pmsimstats-orig/pmsimstats/R/buildtrialdesign.R`,
lines 86--106.

The original computes `tsd` as the cumulative duration of off-drug
intervals, measured in interval widths (`t_wk`):

```r
t_wk <- c(4, 4, 1, 1, 1, 1, 4, 4)
tsd  <- t_wk * (od != 1)
for (iTP in 2:nP) {
  if (od_intervals[iTP] > 0)
    tod[iTP] <- tod[iTP] + tod[iTP-1]
  else
    tsd[iTP] <- tsd[iTP] + tsd[iTP-1]
}
tsd <- everondrug * tsd
```

At the first off-drug timepoint, `tsd` equals the interval width
leading into that timepoint (1 week for weeks 10--12 in the clustered
phase, 4 weeks for the spaced phase). It is never zero.

The implicit assumption: the drug was discontinued at the *start* of
the interval, so by the measurement at the end of the interval, the
participant has been off drug for the full interval duration.

### pmsimstats2025: `carryover_power_simulation.R`

**File:** `carryover_focus/carryover_power_simulation.R`,
lines 332--343.

The 2025 code computes `tsd` as the difference between the current
calendar week and the week at which discontinuation was detected:

```r
just_discontinued = treatment == 0 &
  lag(treatment, default = 0) == 1
discontinuation_week = if_else(
  just_discontinued, week, NA_real_
)
discontinuation_week = zoo::na.locf(
  discontinuation_week, na.rm = FALSE
)
tsd = if_else(
  treatment == 0 & !is.na(discontinuation_week),
  pmax(week - discontinuation_week, local_tsd_min),
  0
)
```

At the first off-drug timepoint, `week == discontinuation_week`, so
`tsd = 0`. The implicit assumption: discontinuation occurs *at* the
measurement week, not before it.

## The Hybrid Design

Both implementations use the same measurement weeks and path
definitions. The Hybrid design (called 'primary N-of-1 design' in
Hendrickson) has:

- **Measurement weeks:** 4, 8, 9, 10, 11, 12, 16, 20
- **Interval widths (`t_wk`):** 4, 4, 1, 1, 1, 1, 4, 4

The four `ondrug` vectors from the vignette
(`Produce_Publication_Results_1_generate_data.Rmd`, lines 54--60):

```
Path A: c(1, 1, 1, 1, 0, 0, 1, 0)
Path B: c(1, 1, 1, 1, 0, 0, 0, 1)
Path C: c(1, 1, 1, 0, 0, 0, 1, 0)
Path D: c(1, 1, 1, 0, 0, 0, 0, 1)
```

These correspond to pmsimstats2025 Paths 1--4 respectively, confirmed
by matching the `case_when` logic in `create_hybrid_design()` at the
three decision points (weeks 10, 16, 20).

## Path-by-Path Comparison

In the tables below, 'Orig' refers to the original
`buildtrialdesign()` output and '2025' refers to
`carryover_power_simulation.R`. Only off-drug rows are shown (on-drug
rows have `tsd = 0` in both implementations).

### Path A / Path 1: `ondrug = c(1, 1, 1, 1, 0, 0, 1, 0)`

| Week | ondrug | Orig tsd | 2025 tsd | Match |
|-----:|:------:|---------:|---------:|:-----:|
|    4 |      1 |        0 |        0 | Yes   |
|    8 |      1 |        0 |        0 | Yes   |
|    9 |      1 |        0 |        0 | Yes   |
|   10 |      1 |        0 |        0 | Yes   |
|   11 |      0 |    **1** |    **0** | **No**|
|   12 |      0 |    **2** |    **1** | **No**|
|   16 |      1 |        0 |        0 | Yes   |
|   20 |      0 |    **4** |    **0** | **No**|

### Path B / Path 2: `ondrug = c(1, 1, 1, 1, 0, 0, 0, 1)`

| Week | ondrug | Orig tsd | 2025 tsd | Match |
|-----:|:------:|---------:|---------:|:-----:|
|    4 |      1 |        0 |        0 | Yes   |
|    8 |      1 |        0 |        0 | Yes   |
|    9 |      1 |        0 |        0 | Yes   |
|   10 |      1 |        0 |        0 | Yes   |
|   11 |      0 |    **1** |    **0** | **No**|
|   12 |      0 |    **2** |    **1** | **No**|
|   16 |      0 |    **6** |    **5** | **No**|
|   20 |      1 |        0 |        0 | Yes   |

### Path C / Path 3: `ondrug = c(1, 1, 1, 0, 0, 0, 1, 0)`

| Week | ondrug | Orig tsd | 2025 tsd | Match |
|-----:|:------:|---------:|---------:|:-----:|
|    4 |      1 |        0 |        0 | Yes   |
|    8 |      1 |        0 |        0 | Yes   |
|    9 |      1 |        0 |        0 | Yes   |
|   10 |      0 |    **1** |    **0** | **No**|
|   11 |      0 |    **2** |    **1** | **No**|
|   12 |      0 |    **3** |    **2** | **No**|
|   16 |      1 |        0 |        0 | Yes   |
|   20 |      0 |    **4** |    **0** | **No**|

### Path D / Path 4: `ondrug = c(1, 1, 1, 0, 0, 0, 0, 1)`

| Week | ondrug | Orig tsd | 2025 tsd | Match |
|-----:|:------:|---------:|---------:|:-----:|
|    4 |      1 |        0 |        0 | Yes   |
|    8 |      1 |        0 |        0 | Yes   |
|    9 |      1 |        0 |        0 | Yes   |
|   10 |      0 |    **1** |    **0** | **No**|
|   11 |      0 |    **2** |    **1** | **No**|
|   12 |      0 |    **3** |    **2** | **No**|
|   16 |      0 |    **7** |    **6** | **No**|
|   20 |      1 |        0 |        0 | Yes   |

## The Pattern

The 2025 `tsd` is systematically shorter than the original by exactly
one interval width at the first off-drug measurement after each
discontinuation, and by the same offset at all subsequent off-drug
measurements in that sequence.

The rule:

- **Original:** First off-drug `tsd` = interval width into that
  timepoint (1 week in the clustered phase, 4 weeks in the spaced
  phase). Subsequent off-drug measurements accumulate from there.
- **2025:** First off-drug `tsd` = 0 (because
  `week - discontinuation_week = 0`). Subsequent off-drug measurements
  are offset by the same initial deficit.

The offset is constant within each off-drug run but varies across
discontinuation events depending on the interval width:

- Weeks 10--12 (clustered phase, `t_wk = 1`): offset = 1 week
- Week 16 (spaced phase, `t_wk = 4`): offset = 4 weeks
- Week 20 (spaced phase, `t_wk = 4`): offset = 4 weeks

## Consequence for `Dbc`

The `Dbc` predictor is computed identically in both implementations:

$$
\text{Dbc} = \left(\frac{1}{2}\right)^{\text{tsd} / t_{1/2}}
$$

The different `tsd` values produce dramatically different `Dbc` at
boundary measurements. For a Path C / Path 3 participant at
$t_{1/2} = 0.1$ weeks:

| Week | Orig tsd | Orig Dbc          | 2025 tsd | 2025 Dbc |
|-----:|---------:|------------------:|---------:|---------:|
|   10 |        1 | $(1/2)^{10} \approx 0.001$ |        0 |    1.000 |
|   11 |        2 | $(1/2)^{20} \approx 0$     |        1 |    0.001 |
|   12 |        3 | $(1/2)^{30} \approx 0$     |        2 | $\approx 0$ |
|   20 |        4 | $(1/2)^{40} \approx 0$     |        0 |    1.000 |

At the longer half-life $t_{1/2} = 0.2$ weeks:

| Week | Orig tsd | Orig Dbc          | 2025 tsd | 2025 Dbc |
|-----:|---------:|------------------:|---------:|---------:|
|   10 |        1 | $(1/2)^{5} = 0.031$        |        0 |    1.000 |
|   11 |        2 | $(1/2)^{10} \approx 0.001$ |        1 |    0.031 |
|   12 |        3 | $(1/2)^{15} \approx 0$     |        2 |    0.001 |
|   20 |        4 | $(1/2)^{20} \approx 0$     |        0 |    1.000 |

The original produces near-zero carryover at boundary measurements for
short half-lives. The 2025 code produces full carryover (`Dbc = 1.000`)
at every first off-drug timepoint regardless of half-life.

## Two Separate Issues

The `tsd` discrepancy and the published power behavior are related
but distinct problems. They must be addressed in sequence.

### Issue 1: The 2025 reimplementation does not match the original

The $t_{\text{sd}} = 0$ step-function power cliff documented in
`predictor_comparison.tex` is caused by the 2025 `tsd` computation,
not by Hendrickson's original logic. Under the original computation:

- Boundary `tsd` is always at least 1 week (clustered phase) or
  4 weeks (spaced phase).
- For $t_{1/2} = 0.1$ weeks, boundary `Dbc` $\approx 0.001$ (not
  1.000).
- The step-function artifact does not arise because the boundary
  measurements carry negligible carryover contamination.

**Required fix:** Adopt the original interval-based `tsd`
computation so that pmsimstats2025 faithfully reproduces the
Hendrickson simulation as a baseline. This is a prerequisite for
any further analysis of power behavior.

### Issue 2: The original itself shows excessive power dropoff

Independent of the 2025 bug, the published Hendrickson heatmaps
show power declining more rapidly with increasing carryover
half-life than pharmacological reasoning would predict. This is
the motivating observation for the entire carryover_focus analysis.

Once the 2025 `tsd` is corrected to match the original, the
excessive power dropoff (if it persists) must have a different
mechanism than the $t_{\text{sd}} = 0$ discontinuity. Candidate
explanations include:

1. **Scale factor mismatch in the original code.**
   `generateData()` applies `scalefactor = 2` to the carryover
   decay in the DGP (line 96: `(1/2)^(scalefactor * tsd /
   t_half)`), while `lme_analysis()` uses
   `carryover_scalefactor = 1` when computing `Dbc` (line 109:
   `(1/2)^(op$carryover_scalefactor * tsd / op$carryover_t1half)`).
   This means the DGP generates faster decay than the analysis
   model assumes, creating a systematic misspecification between
   the true carryover and the modeled carryover.

2. **Cumulative carryover contamination of biological response
   means.** In `generateData()` (line 96), carryover is added to
   `brmeans` cumulatively across timepoints:
   `brmeans[p] <- brmeans[p] + brmeans[p-1] * (1/2)^(sf * tsd /
   t_half)`. This recursive structure means that the biological
   response at off-drug timepoints depends on the entire preceding
   on-drug trajectory, not just the most recent on-drug value.

3. **Interaction between carryover and the correlation structure.**
   The original DGP generates correlated random components across
   timepoints, and the carryover adjustment modifies the mean
   structure but not the correlation matrix. This could produce
   model misspecification that grows with the half-life.

These candidates require empirical investigation by running the
corrected 2025 simulation (with original-matching `tsd`) and
comparing its power curves against the published results.

### The role of `tsd_min` going forward

The `tsd_min` parameter was introduced in pmsimstats2025 as a
correction for the $t_{\text{sd}} = 0$ discontinuity. After fixing
the `tsd` computation to match the original, the `tsd_min`
parameter is no longer needed for that purpose.

However, `tsd_min` (or something like it) may still have value as
a *design improvement* on the original. If the original's own
power behavior proves to be excessive, a minimum washout interval
represents a pharmacologically defensible modification: in clinical
practice, there is always a nonzero gap between the last dose and
the next assessment. Whether such a modification is appropriate
depends on the results of the corrected simulation, and should be
evaluated as a deliberate modeling choice rather than a bug fix.

## Recommended Next Steps

1. **Fix `tsd` in `carryover_power_simulation.R`** to use the
   original interval-based accumulation logic. Verify that the
   corrected code produces identical `tsd` values to
   `buildtrialdesign()` for all four paths.

2. **Re-run the two-scenario comparison** with the corrected `tsd`
   (and `tsd_min = 0`) to establish a faithful reproduction of
   Hendrickson's simulation.

3. **Compare power curves** from the corrected 2025 simulation
   against the published heatmaps. If power still drops too
   rapidly, investigate the scale factor mismatch and cumulative
   carryover mechanisms listed above.

4. **Evaluate `tsd_min` as a design modification**, not a bug fix,
   by running a separate set of simulations with `tsd_min > 0` and
   comparing against the baseline.

## Source Files

- **Original:** `~/prj/alz/01c-pmsimstats-orig/pmsimstats/R/buildtrialdesign.R`
- **Original vignette:** `~/prj/alz/01c-pmsimstats-orig/pmsimstats/vignettes/Produce_Publication_Results_1_generate_data.Rmd`
- **Original analysis:** `~/prj/alz/01c-pmsimstats-orig/pmsimstats/R/lme_analysis.R`
- **2025 simulation:** `~/prj/alz/01-pmsimstats/pmsimstats2025/analysis/scripts/carryover_focus/carryover_power_simulation.R`
- **2025 analysis document:** `~/prj/alz/01-pmsimstats/pmsimstats2025/analysis/scripts/carryover_focus/predictor_comparison.tex`

---
*Rendered on 2026-03-18 at 09:23 PDT.*
*Source: ~/prj/alz/01-pmsimstats/pmsimstats2025/analysis/scripts/carryover_focus/tsd_comparison.md*
