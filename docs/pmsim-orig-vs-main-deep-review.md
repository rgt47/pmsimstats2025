# Deep Code Review: pmsimstats-orig vs. pmsimstats2025

**Systematic comparison of every element of the
data-generating process, covariance structure, and
analysis model**

---

## 1. Biomarker-Treatment Interaction Mechanism

**This is the most consequential difference between the
two codebases.**

### 1.1 orig: Correlation-based interaction

The biomarker-treatment interaction is encoded entirely
within the covariance matrix. There is no biomarker term
in the mean structure of the response components.

The mechanism:

1. The biomarker (`bm`) and biological response
   components (`tp.br`) are drawn jointly from a single
   MVN.
2. The correlation between `bm` and each `tp.br` is set
   to `c.bm` when on-drug, and scaled down (or to zero)
   when off-drug, via the Ron Thomas logic.
3. `MASS::mvrnorm()` generates the joint sample.
   Participants with higher biomarker values tend to have
   higher BR values (on-drug) purely through the
   multivariate correlation.

The code (`generateData.R`, lines 146-162):

```r
for (p in 1:nP) {
  n1 <- paste(trialdesign$timeptname[p], "br", sep=".")
  if (p > 1) {
    n0 <- paste(trialdesign$timeptname[p-1], "br",
                sep=".")
    mm1 <- means[which(n1 == labels)]
    mm0 <- means[which(n0 == labels)]
    correlations["bm", n1] <-
      correlations[n1, "bm"] <-
        ifelse(brtest[p],
          ifelse(brmeans[p] == 0, 0,
                 (mm1/mm0) * modelparam$c.bm),
          modelparam$c.bm)
  }
}
```

Key features:

- `brtest[p]` is TRUE when the *raw* (pre-carryover) BR
  mean at timepoint p is zero (off-drug per the Gompertz
  of `tod`)
- When off-drug with carryover: correlation is scaled by
  `mm1/mm0` (ratio of carryover-adjusted BR means)
- When off-drug without carryover: `brmeans[p] == 0`, so
  correlation is set to 0
- When on-drug: correlation is `c.bm` (unscaled)
- **Timepoint p=1 is never assigned a bm correlation**
  (the `if (p > 1)` guard skips it)

The interaction strength is controlled by `c.bm`
(a correlation, range 0-1). No explicit biomarker term
appears in the mean structure. The BR mean at each
timepoint is `modgompertz(tod, ...)` regardless of
biomarker value.

### 1.2 main: Mean-structure interaction

The biomarker-treatment interaction is encoded directly
in the mean structure of the BR component:

```r
effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered)
br_on_drug = weeks_on_drug * effective_BR_rate
```

High-biomarker participants get a deterministically larger
BR mean. The interaction strength is controlled by
`biomarker_moderation` (range 0-1).

The covariance matrix uses `c.bm` to set the correlation
between biomarker and BR components, but this is
**constant across all timepoints** and **does not vary by
treatment status**:

```r
Sigma_12[br_idx, 1] <- effective_c.bm *
  within_subject_sd * biomarker_sd
```

There is no Ron Thomas scaling. No treatment-dependent
correlation adjustment.

### 1.3 Consequence

The orig produces the interaction via *differential
correlation*: the on/off-drug covariance contrast is the
signal. This is stochastic and requires larger samples to
detect.

The main produces the interaction via *differential
means*: each participant carries a deterministic
biomarker-modulated treatment effect. This is more
detectable per participant.

**These are fundamentally different statistical
mechanisms.** Power estimates from the two are not
comparable for the same nominal parameter values.

---

## 2. Per-Path vs. Shared Sigma Generation

### 2.1 orig: Path-specific sigma matrices

Data is generated **separately for each path** in
`generateSimulatedResults.R` (lines 133-149):

```r
nP <- length(td)
for (iP in 1:nP) {
  dat <- generateData(modelparam, respparam, blparam,
                      td[[iP]], FALSE, TRUE)
  dat[, path := iP]
}
```

Each call to `generateData()` builds a **path-specific**
sigma matrix. The sigma depends on the path's `ondrug`
vector because:

- BR means change (Gompertz of path-specific `tod`)
- BM-BR correlations change (Ron Thomas scaling depends
  on each timepoint's on/off-drug status for that path)
- ER SDs change (scaled by path-specific expectancy)

Path A (drug-first in CO) has `c.bm` correlations at
timepoints 1-4 and zero at 5-8. Path B has the reverse.
Participants on different paths are drawn from
**different multivariate distributions**.

### 2.2 main: Single sigma per design

All participants in a given design share one cached
sigma matrix:

```r
sigma_cache[[cache_key]] <-
  build_sigma_guaranteed_pd(weeks, c.bm, params)
```

Path-specific treatment effects are applied **after** the
MVN draw, in the mean-structure computation:

```r
BR_mean = if_else(treatment == 1,
                  br_on_drug, br_off_drug)
BR = BR_mean + br_random
```

The `br_random` values (from the shared MVN) represent
the stochastic component. `BR_mean` is path-specific
and deterministic.

### 2.3 Consequence

In the orig, different paths induce different covariance
structures. This is statistically necessary for the
correlation-based interaction to work: the biomarker must
correlate with BR differently depending on treatment
status. In the main, the shared sigma means the biomarker
correlates with BR identically regardless of treatment
status.

---

## 3. Response Growth Model

| Aspect                | orig                           | main                      |
|-----------------------|--------------------------------|---------------------------|
| Growth function       | Modified Gompertz              | Linear rate               |
| BR on-drug            | `modgompertz(tod, max, d, r)`  | `weeks_on_drug * BR_rate` |
| ER                    | `modgompertz(tpb, max, d, r)`  | `weeks_w_expect * ER_rate`|
|                       | `* expectancy`                 |                           |
| TR                    | `modgompertz(t_cum, max, d, r)`| `weeks_in_trial * TR_rate`|
| Parameters per factor | 4 (max, disp, rate, sd)        | 1 (rate)                  |
| Saturation            | Yes (asymptote)                | No (unbounded linear)     |

The orig's modified Gompertz (`utilities.R`, lines 34-40):

```r
modgompertz <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y <- y * (maxr / (maxr - vert_offset))
  return(y)
}
```

This produces a sigmoidal curve that starts at 0 (t=0)
and saturates at `maxr`. The main uses a simple linear
model with no saturation.

---

## 4. Cumulative Drug Exposure (`tod` vs. `weeks_on_drug`)

### 4.1 Units

The orig computes `tod` as cumulative time-weighted
exposure in **weeks**:

```r
t_wk = c(4, 4, 1, 1, 1, 1, 4, 4)  # interval durations
od_intervals = t_wk * od
tod[i] = tod[i] + tod[i-1]  # when on-drug
```

The main counts on-drug **measurement occasions**:

```r
weeks_on_drug = cumsum(treatment)
```

For the Hybrid design Path A:

```
Week:    4    8    9   10   11   12   16   20
orig:    4    8    9   10    0    0    4    0
main:    1    2    3    4    4    4    5    5
```

### 4.2 Off-drug behavior

The orig's `tod` drops to 0 during off-drug periods
and **resets** when drug resumes. At week 16 (drug
re-initiated after off-drug at 11-12), `tod = 4` (just
the 4-week interval), not 14 (total lifetime exposure).

The main's `cumsum` is monotonically non-decreasing and
**plateaus** at the last on-drug count.

### 4.3 Consequence

The orig treats drug re-initiation as a partial restart
(the Gompertz re-enters its steep early phase). The main
treats it as continuation of cumulative exposure.

---

## 5. Time Since Discontinuation (`tsd`)

### 5.1 Computation

The orig accumulates off-drug interval durations:

```r
tsd <- t_wk * (od != 1)
for (iTP in 2:length) {
  if (od_intervals[iTP] == 0)
    tsd[iTP] <- tsd[iTP] + tsd[iTP-1]
}
everondrug <- (cumulative(od) > 0)
tsd <- everondrug * tsd
```

The main measures elapsed time since the first off-drug
observation:

```r
just_discontinued = treatment == 0 &
  lag(treatment, default = 0) == 1
discontinuation_week = if_else(
  just_discontinued, week, NA_real_)
discontinuation_week = zoo::na.locf(...)
tsd = if_else(treatment == 0 & !is.na(disc_week),
              week - disc_week, 0)
```

### 5.2 Comparison (Hybrid Path A)

```
Week:    4    8    9   10   11   12   16   20
orig:    0    0    0    0    1    2    0    4
main:    0    0    0    0    0    1    0    0
```

Systematic offset: orig's tsd is one measurement interval
larger at every off-drug timepoint. The orig measures
"how long off-drug by measurement time." The main
measures "time since the transition event."

---

## 6. Carryover Formula

### 6.1 Scale factor

The orig uses `scalefactor = 2`:

```r
brmeans[p] <- brmeans[p] +
  brmeans[p-1] * (1/2)^(scalefactor * tsd / t1half)
```

The main omits the scale factor:

```r
carryover_effect = (1/2)^(tsd / carryover_t_half)
```

### 6.2 Where carryover is applied

The orig applies carryover to the **BR mean vector**
before MVN generation. Off-drug BR means are adjusted
upward by the decayed previous BR mean. This carryover
then feeds into the Ron Thomas correlation scaling.

The main applies carryover **after** MVN generation,
in the deterministic mean-structure computation.

### 6.3 Combined impact (Hybrid Path A, t1/2 = 0.2)

```
Week:         11       12       20
orig:         0.00098  ~0       ~0
main:         1.000    0.031    1.000
```

Three-order-of-magnitude difference at the first
off-drug timepoint.

---

## 7. Covariance Structure

### 7.1 Within-factor autocorrelation

| orig                         | main                            |
|------------------------------|---------------------------------|
| Compound symmetry:           | Time-based AR(1):               |
| `corr[i,j] = c.br`          | `Cov[i,j] = sigma^2 *`         |
| (constant for all lags)      | `rho^abs(t_i - t_j)`           |

### 7.2 Cross-factor correlation

| orig                         | main                            |
|------------------------------|---------------------------------|
| Same-time: `c.cf1t`         | Same-time: `c.cf1t * sigma^2`  |
| Diff-time: `c.cfct`         | Diff-time: `c.cfct * sigma^2`  |
|                              | `* 0.9^time_lag`                |

The orig uses constant cross-factor correlations; the
main decays them with time lag.

### 7.3 BM-BR correlation

| orig                         | main                            |
|------------------------------|---------------------------------|
| Treatment-dependent:         | Constant across all TPs:        |
| On-drug: `c.bm`             | `c.bm` for all BR timepoints    |
| Off-drug: scaled by          |                                 |
| `mm1/mm0 * c.bm`            |                                 |
| (Ron Thomas adjustment)      |                                 |
| Timepoint 1: always 0        | Timepoint 1: `c.bm`            |

### 7.4 BM-ER and BM-TR correlation

The orig sets **no** correlation between the biomarker
and ER or TR components (only BR gets `c.bm`).

The main sets BM-ER and BM-TR correlations at
`c.bm * 0.5`:

```r
Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 *
  within_subject_sd * biomarker_sd
Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 *
  within_subject_sd * biomarker_sd
```

This means the main's biomarker correlates with ER and
TR components (at half strength), while the orig's
biomarker correlates only with BR.

### 7.5 BM-Baseline and Baseline-Response correlation

The orig has no explicit baseline-response or
biomarker-baseline correlation beyond what emerges from
the shared MVN. The baseline (`BL`) and biomarker (`bm`)
are dimensions 1-2 of the MVN, but their correlations
with other components are not explicitly set (they remain
at zero in the correlation matrix except for `bm`-`br`).

The main explicitly sets:

```r
# BM-Baseline cross-covariance:
Sigma_22[1,2] = c.bm_baseline * biomarker_sd *
                between_subject_sd

# Baseline-Response cross-covariance:
Sigma_12[all_idx, 2] = c.baseline_resp *
  within_subject_sd * between_subject_sd
```

### 7.6 Sigma matrix dimensions

Both use `(2 + 3*nTP)` dimensions for 2 baseline
variables (bm, BL) and 3 factors (BR, ER, TR) at each
of nTP timepoints. For 8 timepoints: 26x26.

### 7.7 Positive definiteness enforcement

Both use `corpcor::make.positive.definite(sigma, tol=1e-3)`.

The main additionally pre-validates all sigma matrices
before simulation via eigenvalue checking and condition
number reporting.

---

## 8. Expectancy Response SD Scaling

### 8.1 orig

ER standard deviations are scaled by expectancy
(`generateData.R`, line 73):

```r
sds <- c(sds,
  rep(respparam[cat=="pb"]$sd, nP) * trialdesign$e)
```

At blinded timepoints (`e = 0.5`), ER SD is halved. At
open-label (`e = 1.0`), it is at full value.

### 8.2 main

No ER SD scaling. The sigma uses `within_subject_sd`
uniformly for the ER block. Expectancy enters only through
the ER mean:

```r
ER_mean = weeks_with_expectancy * ER_rate
```

### 8.3 Consequence

The orig reduces ER *variance* during blinded phases;
the main reduces only the ER *mean*. The orig produces
tighter ER distributions when blinded.

---

## 9. Outcome Computation

### 9.1 Sign convention

| orig                        | main                            |
|-----------------------------|---------------------------------|
| `Sx = BL - (tv + pb + br)` | `response = BL + BR + ER + TR`  |
| Symptom score (high = bad)  | Improvement (high = good)       |

The sign reversal does not affect the interaction test
(p-value is identical; coefficient flips sign).

### 9.2 Component decomposition

In the orig, the entire value of each component (signal
+ noise) comes from the MVN draw. The Gompertz-derived
values are the *means* of the MVN.

In the main, each component is explicitly decomposed:

```r
BR = BR_mean + br_random
```

where `BR_mean` is deterministic (from treatment schedule
and biomarker) and `br_random` is from the MVN.

---

## 10. Analysis Model

| Aspect              | orig                          | main                          |
|---------------------|-------------------------------|-------------------------------|
| Engine              | `lme4::lmer()`                | `nlme::lme()` or `lmer()`    |
| Residual corr.      | None (iid residuals)          | `corCAR1` (lme only)         |
| Treatment predictor | `Dbc` (0-1 continuous)        | `effective_drug_weeks`        |
|                     |                               | (0 to max_weeks)              |
| Interaction term    | `bm:Dbc`                      | `eff_drug_wks:bm_centered`   |
| Time covariate      | `t` (cumulative)              | `week`                        |
| BM predictor        | Raw `bm`                      | `bm_centered` (standardized) |
| Formula building    | `eval(parse(...))`            | Static `case_when`            |

### 10.1 Dbc vs. effective_drug_weeks

The orig's `Dbc` (`lme_analysis.R`, lines 107-109):

```r
data.m2[Db == TRUE, Dbc := 1]
data.m2[Db == FALSE,
  Dbc := (1/2)^(sf * tsd / t1half)]
```

Range [0, 1]. On-drug = 1; off-drug = exponential decay.

The main's `effective_drug_weeks` (line 697-701):

```r
effective_drug_weeks = case_when(
  treatment == 1 ~ weeks_on_drug,
  treatment == 0 ~ weeks_at_discont * carryover_effect,
  TRUE ~ 0
)
```

Range [0, max(weeks_on_drug)]. On-drug = cumulative count;
off-drug = count * decay.

These produce coefficients on different scales but test
the same null (interaction = 0).

### 10.2 Dynamic formula selection (orig only)

The orig dynamically selects the model formula based on
data characteristics:

```r
modelbase = "Sx~bm+t"
if (varInDb)
  modelbase = paste0(modelbase, "+Dbc+bm*Dbc")
else
  modelbase = paste0(modelbase, "+bm*t")
```

The main uses fixed formulas per design (OL vs. others).

---

## 11. Parameter Defaults

```
+----------------------------+----------+----------+
| Parameter                  | orig     | main     |
+----------------------------+----------+----------+
| c.bm (interaction param)  | 0, .3, .6| 0, .3    |
| c.br, c.pb, c.tv          | 0.8      | 0.75     |
| c.cf1t                     | 0.2      | 0.12     |
| c.cfct                     | 0.1      | 0.05     |
| carryover_t1half (weeks)   | 0, .1, .2| 0, .1, .2|
| scalefactor                | 2        | 1*       |
| N                          | 35, 70   | 35, 70   |
| BR growth                  | Gompertz | Linear   |
|                            | (4 param)| (1 param)|
| BM mean                    | ~125     | 5.0      |
| BM sd                      | ~15      | 2.0      |
| BL mean                    | ~35      | 10.0     |
| BL sd                      | ~16      | 2.0      |
| Iterations                 | 1000     | 500      |
+----------------------------+----------+----------+
* simulation_clustered.R only; pm_functions.R uses 2
```

The orig's baseline/biomarker values are extracted from
real clinical data (CAPS5 scores, standing blood
pressure). The main uses smaller synthetic values.

---

## 12. Architectural Differences

| Aspect                | orig                         | main                       |
|-----------------------|------------------------------|----------------------------|
| Data structures       | `data.table`                 | `tibble` / `tidyverse`     |
| Package structure     | Installable R package        | Research compendium        |
| Exported functions    | 10 (NAMESPACE)               | 0 (empty NAMESPACE)        |
| Parallelism           | Sequential                   | `furrr::future_map_dfr()`  |
| Sigma caching         | None (rebuilt each call)     | Pre-built and cached       |
| Two-stage MVN         | No (full joint draw)         | Yes (conditional draw)     |
| PD pre-validation     | No                           | Yes (eigenvalue check)     |
| Progressive save      | Yes (RDS checkpointing)      | No (in-memory)             |
| Reproducibility       | None (no Docker/renv)        | Docker + renv + Makefile   |
| Test coverage         | 0 tests                      | 1 trivial test             |

---

## 13. Priority for Faithful Replication

To make pmsimstats2025 reproduce the orig's results:

| # | Change                                     | Priority |
|---|--------------------------------------------|:--------:|
| 1 | Replace mean-based BM interaction with      | Critical |
|   | correlation-based interaction (differential |          |
|   | c.bm by treatment status, Ron Thomas       |          |
|   | scaling)                                    |          |
| 2 | Generate data per-path with path-specific   | Critical |
|   | sigma matrices                              |          |
| 3 | Fix tsd (anchor to last on-drug week)       | Critical |
| 4 | Fix tod (time-weighted with gap reset)      | Critical |
| 5 | Restore scalefactor = 2                     | Critical |
| 6 | Implement modified Gompertz                 | High     |
| 7 | Scale ER SD by expectancy                   | Moderate |
| 8 | Skip BM-BR correlation at timepoint 1       | Low      |
| 9 | Remove BM-ER and BM-TR correlations         | Moderate |
|10 | Use orig's clinical parameter values        | Low      |
|   | (or verify results are scale-invariant)     |          |

Items 1-2 are the structural changes that drive the
largest differences in power estimates. Items 3-5 are
the numerical corrections that compound to produce
order-of-magnitude differences in carryover behavior.
Items 6-10 are refinements that improve fidelity.

The autocorrelation structure (compound symmetry vs.
AR(1)) is a deliberate improvement that may be retained,
but its impact should be quantified by running both
structures and comparing power curves.

---

*Rendered on 2026-03-18 at 09:36 PDT.*
*Source: ~/prj/alz/pmsim-orig-vs-main-deep-review.md*
