# Technical Audit: pmsimstats Project Family

**Comprehensive Review of Three Sibling Repositories for
Precision Medicine Trial Simulation**

Author: Technical Audit (Claude)
Date: 2026-03-18

---

---

## 1. Structural Inventory

### 1.1 Repository Overview

| Attribute       | orig                  | main                  | simple             |
|-----------------|-----------------------|-----------------------|--------------------|
| Package name    | pmsimstats            | pmsimstats2025        | pmsimstatssimple   |
| Version         | 0.0.0.900             | 0.1.0                 | 0.0.0.9000         |
| License         | GPL 2.0               | GPL-3                 | GPL-3              |
| Total R lines   | 3,142                 | 15,584                | 562                |
| R files (src)   | 11                    | 44                    | 1                  |
| Test files      | 0                     | 5                     | 1                  |
| Doc files       | 4 (README + 3 Rmd)    | 40+ (md, tex, Rmd)    | 4 (md)             |
| Docker          | No                    | Yes                   | Yes                |
| renv            | No                    | Yes                   | Yes                |
| Makefile        | No                    | Yes                   | Yes                |

### 1.2 Directory Trees

#### 1.2.1 orig (`pmsimstats`)

```
pmsimstats/
+-- DESCRIPTION
+-- NAMESPACE
+-- README.md
+-- R/
|   +-- buildtrialdesign.R      (119 lines)
|   +-- censordata.R            ( 89 lines)
|   +-- comb.R                  (1571 lines)
|   +-- datadocumentation.R     ( 61 lines)
|   +-- generateData.R          (223 lines)
|   +-- generateSimulatedResults.R (255 lines)
|   +-- lme_analysis.R          (217 lines)
|   +-- packagedocumentation.R  ( 19 lines)
|   +-- plottingfunctions.R     (341 lines)
|   +-- scratch_WorkingScript.R (142 lines)
|   +-- utilities.R             (105 lines)
+-- vignettes/
    +-- Produce_Publication_Results_1_generate_data.Rmd
    +-- Produce_Publication_Results_2_Look_at_data.Rmd
    +-- Produce_Publication_Results_3_apply_to_actual_clinical_trial_data.Rmd
```

#### 1.2.2 main (`pmsimstats2025`)

```
pmsimstats2025/
+-- DESCRIPTION
+-- NAMESPACE
+-- CLAUDE.md
+-- Dockerfile
+-- Makefile
+-- renv.lock
+-- R/
|   +-- utils.R
+-- analysis/
|   +-- scripts/
|   |   +-- pm_functions.R               (1301 lines)
|   |   +-- simulation_clustered.R       (1233 lines)
|   |   +-- simulationplus.R             ( 995 lines)
|   |   +-- simulation_evenly_spaced.R   ( 649 lines)
|   |   +-- simulation_carryover_robustness.R
|   |   +-- simulation_carryover_spectrum.R
|   |   +-- simulation_censoring.R
|   |   +-- simulation_clustered_exponential.R
|   |   +-- simulation_clustered_flexible_correlation.R
|   |   +-- simulation_clustered_twoweek.R
|   |   +-- simulation_evenly_spaced_flexible_correlation.R
|   |   +-- simulation_evenly_spaced_powerlaw.R
|   |   +-- simulation_response_params.R
|   |   +-- visualize_heatmaps.R
|   |   +-- visualize_hendrickson_style.R
|   |   +-- visualize_power_heatmaps_4panel.R
|   |   +-- plot_results.R
|   |   +-- plot_censoring.R
|   |   +-- summarize_results.R
|   |   +-- check_model_specification.R (x3)
|   |   +-- debug_effect_recovery.R
|   |   +-- diagnose_nonmonotonic_power.R
|   |   +-- example_parameter_validation.R
|   |   +-- profile_simulation.R
|   |   +-- test_*.R (5 ad-hoc test scripts)
|   |   +-- 00_*.R (setup scripts)
|   |   +-- 02_data_validation.R
|   |   +-- 99_reproducibility_check.R
|   |   +-- carryover_focus/
|   |       +-- carryover_power_simulation.R
|   +-- data/
|   +-- output/
+-- docs/
|   +-- (40+ markdown, tex, and Rmd files)
+-- tests/
|   +-- testthat.R
|   +-- testthat/
|   |   +-- test-utils.R
|   +-- integration/
|       +-- test-analysis-scripts.R
|       +-- test-data-pipeline.R
|       +-- test-report-rendering.R
+-- vignettes/
    +-- 01_generate_simulation_data.Rmd
    +-- quickstart.Rmd
```

#### 1.2.3 simple (`pmsimstats-simple`)

```
pmsimstats-simple/
+-- DESCRIPTION
+-- NAMESPACE
+-- CLAUDE.md
+-- Dockerfile
+-- Makefile
+-- renv.lock
+-- zzcollab.yaml
+-- .github/workflows/
|   +-- r-package.yml
+-- analysis/
|   +-- scripts/
|   |   +-- simulation.R            (562 lines)
|   +-- data/
|   |   +-- README.md
|   +-- output/
|       +-- comparison_whitepaper.md
+-- docs/
|   +-- simplification-plan.md
+-- tests/
    +-- testthat.R
    +-- testthat/
        +-- test-basic.R
```

### 1.3 Source File Manifest

#### 1.3.1 orig -- Exported Functions (from NAMESPACE)

| Function                   | File                        | Lines | Purpose                                           |
|----------------------------|-----------------------------|------:|---------------------------------------------------|
| `buildtrialdesign`         | buildtrialdesign.R          |   119 | Construct trial design from user-friendly inputs   |
| `generateData`             | generateData.R              |   223 | Core MVN data-generating process                   |
| `generateSimulatedResults` | generateSimulatedResults.R  |   255 | Parameter-space iteration wrapper                  |
| `lme_analysis`             | lme_analysis.R              |   217 | LME model fitting and coefficient extraction       |
| `censordata`               | censordata.R                |    89 | Simulate dropout with outcome-dependent censoring  |
| `cumulative`               | utilities.R                 |   105 | Convert intervals to cumulative time (shared)      |
| `modgompertz`              | utilities.R                 |       | Modified Gompertz growth curve                     |
| `reknitsimresults`         | utilities.R                 |       | Recombine chunked simulation outputs               |
| `PlotModelingResults`      | plottingfunctions.R         |   341 | Heatmap visualization of power results             |
| `plotfactortrajectories`   | plottingfunctions.R         |       | Mean trajectory plots by factor                    |

Non-exported / legacy files:

| File                       | Lines | Purpose                                        |
|----------------------------|------:|------------------------------------------------|
| comb.R                     | 1,571 | Combined scratch: duplicated `lme_analysis`,   |
|                            |       | `generateData`, assorted utilities, dead code  |
| scratch_WorkingScript.R    |   142 | Development notes and parameter exploration    |
| datadocumentation.R        |    61 | roxygen2 stubs for embedded datasets           |
| packagedocumentation.R     |    19 | Package-level documentation                    |

#### 1.3.2 main -- Key Functions (in pm_functions.R)

| Function                              | Purpose                                           |
|---------------------------------------|---------------------------------------------------|
| `mod_gompertz`                        | Modified Gompertz (reformulated)                  |
| `calculate_carryover`                 | Multi-model decay (exp, linear, Weibull)          |
| `calculate_bio_response_with_interaction` | BR means with biomarker interaction           |
| `calculate_component_halflives`       | Component-specific half-life mapping              |
| `apply_carryover_to_component`        | Apply decay to BR component                       |
| `build_correlation_matrix`            | Full correlation structure with BM interaction    |
| `generate_data`                       | Primary MVN data generation                       |
| `build_sigma_matrix`                  | Standalone sigma builder for caching              |
| `validate_correlation_structure`      | Eigenvalue diagnostics for sigma                  |
| `validate_parameter_grid`             | Pre-simulation PD validation                      |
| `lme_analysis`                        | Mixed-effects model fitting                       |

#### 1.3.3 simple -- Functions (in simulation.R)

| Function                | Purpose                                              |
|-------------------------|------------------------------------------------------|
| `create_design`         | Trial design specification (switch-based)            |
| `randomize_assignment`  | Per-participant arm randomization                    |
| `calc_carryover`        | Simple exponential decay                             |
| `generate_response`     | Single-timepoint signal computation                  |
| `generate_participant`  | Full participant trajectory                          |
| `generate_trial`        | N-participant trial generation                       |
| `analyze_trial`         | Design-specific ANCOVA analysis                      |

### 1.4 Dependency Comparison

| Package        | orig | main | simple | Notes                        |
|----------------|:----:|:----:|:------:|------------------------------|
| data.table     |  D   |      |        | Core data structure (orig)   |
| lme4           |  D   |  I   |        | `lmer()` fitting             |
| lmerTest       |  D   |  I   |        | p-values for lmer            |
| corpcor        |  D   |  I   |        | PD checking/enforcement      |
| ggplot2        |  D   |  I   |        | Visualization                |
| MASS           |  D   |  I   |        | `mvrnorm()` for MVN draws    |
| svMisc         |  D   |      |        | Progress bar (legacy)        |
| tictoc         |  D   |      |        | Timing utilities             |
| ggpubr         |  D   |      |        | Publication-quality plots    |
| gridExtra      |  D   |      |        | Multi-panel layouts          |
| nlme           |      |  I   |        | `lme()` with `corCAR1`       |
| dplyr          |      |  I   |  *     | Data manipulation            |
| tidyr          |      |  *   |  *     | (via tidyverse)              |
| readr          |      |  I   |  *     | CSV I/O                      |
| tidyverse      |      |  I   |  *     | Meta-package (simple uses)   |
| furrr / future |      |  I   |        | Parallel iteration           |
| foreach        |      |  I   |        | Parallel loops (alternative) |
| doParallel     |      |  I   |        | Parallel backend             |
| zoo            |      |  I   |        | `na.locf`, time-series       |
| patchwork      |      |  I   |        | Plot composition             |
| viridis        |      |  I   |        | Color scales                 |
| scales         |      |  I   |        | Axis formatting              |
| here           |      |  I   |  *     | Working directory management |
| renv           |      |  I   |  I     | Package management           |
| testthat       |      |  I   |  I     | Testing framework            |
| knitr          |  S   |  I   |        | Document rendering           |
| rmarkdown      |  S   |  I   |        | Rmarkdown support            |
| bookdown       |      |  I   |        | Long-form documents          |
| DBI/odbc/etc.  |      |  I   |        | Database (unused in sims)    |
| palmerpenguins |      |  I   |        | Example data (template)      |
| skimr          |      |  I   |        | Data summaries               |
| naniar/visdat  |      |  I   |        | Missing data visualization   |
| emmeans        |      |  I   |        | Estimated marginal means     |
| jsonlite       |      |  I   |        | JSON processing              |
| tibble         |      |  I   |  *     | (via tidyverse)              |

Key: D = Depends, I = Imports, S = Suggests, * = implicit
via tidyverse

### 1.5 Test Inventory

| Repository | Test File                       | Coverage                    | Assertions |
|------------|---------------------------------|-----------------------------|:----------:|
| orig       | (none)                          | --                          |     0      |
| main       | test-utils.R                    | Trivial pass-through        |     1      |
| main       | test-analysis-scripts.R         | Script executability        |    ~3      |
| main       | test-data-pipeline.R            | Data workflow pipeline      |    ~3      |
| main       | test-report-rendering.R         | Rmd rendering               |    ~2      |
| simple     | test-basic.R                    | Trivial pass-through        |     1      |

All three repositories have effectively zero unit-test coverage
of the core simulation logic (DGP, carryover, correlation
matrix, analysis model). The main repository has integration
tests that verify scripts execute without error, but do not
validate statistical correctness.

---

## 2. Statistical Model Comparison

### 2.1 Response Decomposition

All three codebases decompose the outcome as:

```
Y_i(t) = BL_i - [BR_i(t) + ER_i(t) + TR_i(t)]
```

where BL is baseline symptom severity and the bracketed terms
represent improvement (higher = more improvement = lower
symptoms). The sign convention treats the outcome as a symptom
score that decreases with effective treatment.

#### 2.1.1 Biological Response (BR)

| Aspect            | orig                                  | main                                    | simple                                |
|-------------------|---------------------------------------|-----------------------------------------|---------------------------------------|
| Growth curve      | Modified Gompertz                     | Modified Gompertz (reformulated)        | Linear rate                           |
| Formula (on-drug) | `modgompertz(tod, max, disp, rate)`   | `mod_gompertz(tod, max, disp, rate)`    | `br_rate * (1 + bm_mod * bm)`        |
| Gompertz form     | `max*exp(-d*exp(-r*t)) - offset`      | `max * (1 - exp(-d*exp(-r*t)))`         | N/A (linear)                          |
| Biomarker mod.    | Via correlation (c.bm in sigma)       | Via correlation (c.bm in sigma)         | Via mean structure (explicit term)    |
| Time dependence   | Cumulative time-on-drug               | Cumulative time-on-drug                 | Per-week (not cumulative)             |

**Critical difference in Gompertz formulation.** The orig
uses a traditional Gompertz with a vertical offset subtraction
and rescaling:

```r
# orig: modgompertz
y <- maxr * exp(-disp * exp(-rate * t))
vert_offset <- maxr * exp(-disp * exp(-rate * 0))
y <- y - vert_offset
y <- y * (maxr / (maxr - vert_offset))
```

The main reformulates this as:

```r
# main: mod_gompertz
y <- max_value * (1 - exp(-displacement * exp(-rate * time)))
```

These two formulations are NOT algebraically identical. The
orig produces a curve that passes through zero at t=0 and
approaches `maxr` asymptotically, but with different curvature
than the main's version. The main's formulation is a
complementary Gompertz that also passes through zero at t=0
but with a different displacement parameterization. Whether
the `disp` and `rate` parameter values are recalibrated to
produce equivalent trajectories is not documented.

The simple version abandons the Gompertz entirely and uses a
constant rate model where BR is proportional to `br_rate` with
no saturation, yielding a fundamentally different dose-response
shape.

#### 2.1.2 Expectancy Response (ER)

| Aspect          | orig                                 | main                                    | simple                       |
|-----------------|--------------------------------------|-----------------------------------------|------------------------------|
| Growth curve    | Modified Gompertz                    | Modified Gompertz (reformulated)        | Constant rate                |
| Input variable  | `tpb` (time in belief period)        | `tpb` (time in belief period)           | Binary (unblinded + on-drug) |
| Expectancy wt.  | Scaled by `e` (0-1 probability)      | Scaled by expectancy parameter          | Binary: `er_rate` or 0       |
| SD scaling      | `sd * e` (reduces variance too)      | Similar expectancy scaling              | Fixed `sigma_within`         |

#### 2.1.3 Time-Variant Response (TR)

| Aspect          | orig                                 | main                                    | simple                       |
|-----------------|--------------------------------------|-----------------------------------------|------------------------------|
| Growth curve    | Modified Gompertz                    | Modified Gompertz (reformulated)        | **Absent**                   |
| Input variable  | `t_wk_cumulative` (absolute time)    | Cumulative elapsed time                 | --                           |
| Purpose         | Natural disease progression / drift  | Natural disease progression / drift     | Assumed stable baseline      |

The simple version omits TR entirely, assuming no natural
progression over the trial duration. This is a meaningful
simplification that removes one source of confounding from
the simulation.

### 2.2 Covariance Structure

| Aspect                 | orig                                  | main                                          | simple                       |
|------------------------|---------------------------------------|-----------------------------------------------|------------------------------|
| Matrix dimension       | (2 + 3*nTP) x (2 + 3*nTP)            | (2 + 3*nTP) x (2 + 3*nTP)                    | N/A                          |
|                        | e.g., 26x26 for 8 timepoints         | e.g., 26x26 for 8 timepoints                 |                              |
| Autocorrelation (BR)   | Constant: `c.br`                      | Time-based AR(1): `c.br^abs(t_i - t_j)`      | None (iid)                   |
| Autocorrelation (ER)   | Constant: `c.pb`                      | Time-based AR(1): `c.pb^abs(t_i - t_j)`      | None (iid)                   |
| Autocorrelation (TR)   | Constant: `c.tv`                      | Time-based AR(1): `c.tv^abs(t_i - t_j)`      | None                         |
| Cross-factor (same t)  | `c.cf1t`                              | `c.cf1t`                                      | None                         |
| Cross-factor (diff t)  | `c.cfct`                              | `c.cfct * 0.9^time_lag`                       | None                         |
| BM-BR correlation      | `c.bm` (with Ron Thomas adj.)        | `c.bm` (on-drug only, 0 off-drug)            | N/A (mean structure)         |
| PD enforcement         | `corpcor::make.positive.definite`     | `corpcor::make.positive.definite`             | Not needed                   |
| PD tolerance           | `tol = 1e-3`                          | `tol = 1e-3`                                 | --                           |
| Data generation        | `MASS::mvrnorm`                       | `MASS::mvrnorm`                               | `rnorm` (indep. components)  |
| Variance model         | Full MVN (all components jointly)     | Full MVN (all components jointly)             | Random intercept + iid       |

**Key structural differences:**

1. The orig uses **constant** autocorrelation across all
   time lags (compound symmetry within each factor). The main
   uses **time-based AR(1)** where correlation decays with
   actual elapsed time. This is a substantive methodological
   improvement in the main.

2. The orig's biomarker-BR correlation uses a Ron Thomas
   modification that scales `c.bm` by a ratio of means
   (`mm1/mm0`) when `brtest[p]` is TRUE (i.e., when the raw
   BR mean is zero at a timepoint). The main simplifies this
   to a binary rule: `c.bm` when on treatment, 0 when off
   treatment.

3. The simple bypasses the covariance structure entirely,
   placing the biomarker interaction directly in the mean
   structure (`br_rate * (1 + bm_mod * bm)`) rather than
   encoding it as a correlation. This is a fundamentally
   different statistical mechanism for the biomarker-treatment
   interaction.

### 2.3 Carryover Model

| Aspect               | orig                                          | main                                                   | simple                               |
|-----------------------|-----------------------------------------------|--------------------------------------------------------|--------------------------------------|
| Decay function        | Exponential only                              | Exponential, Linear, Weibull                           | Exponential only                     |
| Formula (exponential) | `BR[p-1] * (0.5)^(sf * tsd / t_half)`        | `prev_effect * (0.5)^(sf * tsd / t_half)`             | `0.5^(tsd / t_half)`                |
| Scale factor default  | 2                                             | 2                                                      | 1 (implicit)                         |
| Half-life values      | Not specified in code (set externally)        | 0, 0.1, 0.2 weeks                                     | 0, 0.5, 1, 2 weeks                  |
| Components affected   | BR only                                       | BR only                                                | BR only                              |
| Boundary handling     | Skipped if `tsd == 0`                         | Skipped if off-drug time is zero                       | Returns 0 if `halflife <= 0`         |
| Where applied         | DGP means (generateData)                      | DGP means (pm_functions) AND analysis model            | DGP (generate_response) AND analysis |

**Scale factor discrepancy.** The orig and main use
`scalefactor = 2` by default, which doubles the effective
decay rate. At `t_half = 0.1` weeks with `scalefactor = 2`,
the effective half-life is 0.05 weeks (~8.4 hours). The
simple omits the scale factor (effectively `scalefactor = 1`),
so its decay at `t_half = 0.5` weeks gives a true half-life
of 0.5 weeks. These are not comparable without adjustment.

**Half-life value divergence.** The main's CLAUDE.md explicitly
warns against using the values 0.5, 1.0, 2.0 weeks, calling
them "legacy" and "incorrect." Yet the simple uses exactly
these values. The simple's CLAUDE.md inherits this warning
but the code contradicts it.

### 2.4 Trial Designs Supported

| Design                         | orig | main | simple |
|--------------------------------|:----:|:----:|:------:|
| Open-Label (OL)                | Yes  | Yes  | Yes*   |
| Open-Label + Blinded Disc.     | Yes  | Yes  | Yes    |
| Traditional Crossover (CO)     | Yes  | Yes  | Yes    |
| Hybrid (N-of-1)                | Yes  | Yes  | Yes    |
| Parallel                       | Yes  | Yes  | No     |

*OL is defined in simple but excluded from the main
simulation because it has no drug comparison arm.

**Design specification mechanism:**

- orig: User-specified via `buildtrialdesign()` with
  arbitrary timepoints, expectancy vectors, and ondrug
  path lists. Highly flexible.
- main: Hard-coded design functions in simulation scripts
  (e.g., `create_hybrid_design()`). Less flexible but
  more reproducible for the specific comparison.
- simple: Switch-based `create_design()` with four
  pre-defined designs. Simplest but least extensible.

### 2.5 Analysis Model

| Aspect                | orig                                             | main                                              | simple                                 |
|-----------------------|--------------------------------------------------|---------------------------------------------------|----------------------------------------|
| Fitting function      | `lme4::lmer()`                                   | `lme4::lmer()` or `nlme::lme()`                  | `stats::lm()`                          |
| Correlation structure | None (random effects only)                       | Optional `corCAR1` (nlme only)                    | None                                   |
| Random effects        | `(1|ptID)` or `(1+t|ptID)`                       | `(1|participant_id)`                              | None (ANCOVA on means)                 |
| Fixed effects (CO/HY) | `Sx ~ bm + t + Dbc + bm*Dbc`                    | `response ~ eff_drug_wks * bm_centered + week`    | `drug_effect ~ biomarker`              |
| Fixed effects (OL)    | `Sx ~ bm + t + bm*t`                            | `response ~ week * bm_centered`                   | `mean_response ~ biomarker`            |
| Interaction term      | `bm:Dbc` or `bm:t`                              | `effective_drug_weeks:bm_centered`                | `biomarker` (main effect on diff)      |
| Carryover in analysis | `Dbc = (0.5)^(sf*tsd/t_half)` (continuous)      | `effective_drug_weeks` (decayed cumulative)        | Optional subtraction from response     |
| Two-scenario design   | Single model per condition                       | Yes: model_carryover TRUE vs FALSE                | Yes: carryover_adjust TRUE vs FALSE    |

**Key analytical differences:**

1. The orig constructs a continuous drug indicator `Dbc` that
   equals 1 on-drug and decays exponentially off-drug. The
   biomarker interaction is then `bm * Dbc`. This is an
   elegant approach that embeds carryover directly into the
   predictor variable.

2. The main constructs `effective_drug_weeks` (cumulative
   drug exposure with decay) and tests
   `effective_drug_weeks * bm_centered`. It also compares
   this against a misspecified model using `weeks_on_drug`
   (no decay) to quantify the cost of ignoring carryover.

3. The simple collapses repeated measures to phase means
   and runs simple linear regression. This sacrifices
   within-subject temporal information but eliminates all
   convergence issues associated with mixed-effects models.

### 2.6 Censoring / Dropout

| Aspect                | orig                                   | main                                 | simple                        |
|-----------------------|----------------------------------------|--------------------------------------|-------------------------------|
| Implemented           | Yes (`censordata()`)                   | Yes (simulation_censoring.R)         | Yes (per-week probability)    |
| Mechanism             | Time-dependent + outcome-dependent     | Multiple patterns                    | MCAR only                     |
| Parameters            | `beta0` (rate), `beta1` (outcome       | 5 censoring scenarios                | `dropout_prob` (0-0.10/week)  |
|                       | fraction), `eb1` (exponent)            |                                      |                               |
| Time dependence       | Proportional to interval duration      | Various                              | Constant per week             |
| Outcome dependence    | Yes (weighted by delta scores)         | Yes (some scenarios)                 | No (MCAR)                     |
| Implementation        | Post-hoc NA insertion                  | Post-hoc NA insertion                | Row filtering                 |

### 2.7 Side-by-Side Summary

```
+--------------------+-------------------+---------------------+-----------------+
| Component          | orig              | main                | simple          |
+--------------------+-------------------+---------------------+-----------------+
| Growth curve       | Modified Gompertz | Modified Gompertz   | Linear rate     |
|                    | (offset-rescale)  | (complementary)     |                 |
+--------------------+-------------------+---------------------+-----------------+
| Factors modeled    | BR, ER, TR        | BR, ER, TR          | BR, ER          |
+--------------------+-------------------+---------------------+-----------------+
| BM interaction     | Correlation-based | Correlation-based   | Mean-based      |
+--------------------+-------------------+---------------------+-----------------+
| Covariance         | Compound symmetry | AR(1), time-based   | Indep. + RE     |
+--------------------+-------------------+---------------------+-----------------+
| Carryover models   | Exponential       | Exp, Linear, Weibull| Exponential     |
+--------------------+-------------------+---------------------+-----------------+
| Scale factor       | 2 (default)       | 2 (default)         | 1 (implicit)    |
+--------------------+-------------------+---------------------+-----------------+
| Half-lives (wk)    | External          | 0, 0.1, 0.2        | 0, 0.5, 1, 2   |
+--------------------+-------------------+---------------------+-----------------+
| Analysis model     | lmer              | lmer or lme+corCAR1 | lm (ANCOVA)     |
+--------------------+-------------------+---------------------+-----------------+
| Censoring          | Time + outcome    | Multiple patterns   | MCAR only       |
+--------------------+-------------------+---------------------+-----------------+
| Parallelism        | Sequential        | furrr/future        | Sequential      |
+--------------------+-------------------+---------------------+-----------------+
| PD enforcement     | corpcor           | corpcor             | Not needed      |
+--------------------+-------------------+---------------------+-----------------+
```

---

## 3. Code Quality and Architecture Audit

### 3.1 Design Patterns

#### 3.1a orig

- **Modularity:** Moderate. Functions are separated into
  logical files (design, data gen, analysis, plotting), but
  `comb.R` (1,571 lines) duplicates much of the other code
  and contains scratch work.
- **Parameter specification:** Data.table-based parameter
  structures with named columns. No formal validation.
- **Parallelism:** Sequential only. The `svMisc` progress
  bar and `tictoc` timing are the only concessions to long
  runtimes.
- **Output management:** Progressive RDS saving with
  numbered suffixes (`_save1`, `_save2`, etc.). Fragile
  naming convention.
- **Reproducibility:** No Docker, renv, or Makefile. Depends
  on the user having all packages installed at compatible
  versions.

#### 3.1b main

- **Modularity:** Good separation between shared functions
  (`pm_functions.R`) and simulation scripts. However, the
  simulation scripts are large monoliths (1,233 lines for
  `simulation_clustered.R`) that inline much logic that
  could be factored out.
- **Parameter specification:** Explicit parameter grids
  with `expand.grid()` or `expand_grid()`. Pre-simulation
  validation via `validate_parameter_grid()`.
- **Parallelism:** `furrr::future_map_dfr()` for
  iteration-level parallelism. Configurable via
  `future::plan()`.
- **Output management:** Timestamped filenames, CSV export,
  structured logging to file.
- **Reproducibility:** Docker + renv + Makefile. The most
  robust of the three.

#### 3.1c simple

- **Modularity:** Single file by design. All functions,
  parameters, and execution in one script. Readable but
  not reusable.
- **Parameter specification:** Hard-coded constants at
  script top. Clean and transparent.
- **Parallelism:** None. Acceptable given ~4 minute
  runtime.
- **Output management:** Timestamped CSV and PNG output.
- **Reproducibility:** Docker + renv + Makefile + CI/CD
  via GitHub Actions.

### 3.2 Code Issues

#### 3.2a orig -- Issues

1. **Dead code: `comb.R` (1,571 lines).** This file
   contains duplicated versions of `lme_analysis()`,
   `generateData()`, and various utilities alongside
   scratch work and commented-out experiments. It is
   actively loaded by the package (listed in R/) and
   creates namespace pollution and confusion about which
   version of a function is authoritative.

2. **`eval(parse(...))` in `generateData.R` (lines
   204-214).** The outcome computation uses string
   construction and `eval(parse())` to dynamically build
   data.table expressions:
   ```r
   evalstring <- paste(
     "dat[,D_", n, ":=sum(",
     paste(paste(n, cl, sep='.', collapse=',')),
     "),by='ptID']", sep="")
   eval(parse(text=evalstring))
   ```
   This is fragile, hard to debug, and a potential
   injection vector if inputs are not sanitized. It should
   be replaced with programmatic data.table column
   operations.

3. **Hard-coded indices in `generateData.R` (line 89).**
   `names(brtest) <- labels[19:26]` assumes exactly 8
   timepoints and 2 baseline variables. This breaks for
   any other trial design length.

4. **Missing error handling.** No `tryCatch()` around
   `lmer()` calls (the code comments note "should insert
   tryCatch here, but not working!"). No validation of
   input parameters. No check for degenerate covariance
   matrices before `mvrnorm()`.

5. **`scratch_WorkingScript.R` in R/.** Development notes
   do not belong in the package source directory.

6. **Unqualified function calls.** Functions from
   `data.table`, `MASS`, `lme4`, and `corpcor` are called
   without namespace qualification, relying on `Depends:`
   to attach them. This is fragile and against modern R
   packaging conventions.

7. **No tests whatsoever.** Zero test files.

8. **Symmetry bug (fixed in-place).** Line 134 of
   `generateData.R` contains a note about fixing a typo
   where both `correlations[n1,n2]` and `correlations[n2,n1]`
   were set with the same indices, breaking matrix symmetry.

#### 3.2b main -- Issues

1. **Massive DESCRIPTION dependency list.** 38 packages in
   Imports, including database drivers (`DBI`, `odbc`,
   `RMySQL`, `RPostgres`, `RSQLite`), example datasets
   (`palmerpenguins`), and tools (`conflicted`, `sessioninfo`,
   `skimr`, `naniar`, `visdat`) that are never used in the
   simulation code. This is a zzcollab template artifact
   that inflates the dependency footprint.

2. **Empty NAMESPACE.** The package exports nothing.
   All functional code lives in `analysis/scripts/` rather
   than `R/`. This means the "package" structure is
   vestigial; it functions as a research compendium, not an
   installable R package.

3. **Proliferation of simulation script variants.** There
   are 13+ simulation scripts (`simulation_clustered.R`,
   `simulation_clustered_exponential.R`,
   `simulation_clustered_flexible_correlation.R`,
   `simulation_clustered_twoweek.R`,
   `simulation_evenly_spaced.R`,
   `simulation_evenly_spaced_flexible_correlation.R`,
   `simulation_evenly_spaced_powerlaw.R`,
   `simulationplus.R`, etc.) with substantial code
   duplication. Shared logic (design creation, parameter
   grid setup, output formatting) is copied across files
   rather than factored into `pm_functions.R`.

4. **Ad-hoc test scripts in analysis/scripts/.** Files like
   `test_autocorr_for_cbm_06.R`,
   `test_boundary_autocorr_06.R`, `test_cbm_06.R`,
   `test_new_params.R`, `test_pd_boundaries.R` are not
   integrated into the testthat framework. They are
   standalone exploration scripts that could silently break.

5. **Trivial test coverage.** The single testthat test is:
   ```r
   test_that("package loads correctly", {
     expect_true(TRUE)
   })
   ```
   This provides zero coverage of the 1,301-line
   `pm_functions.R` or any simulation logic.

6. **Documentation volume vs. coherence.** Over 40
   documentation files exist, but many overlap, some
   contradict each other, and there is no clear canonical
   reference. A new collaborator would struggle to identify
   which document is authoritative.

7. **Deprecated function still present.**
   `calculate_carryover_adjusted_correlations()` is marked
   deprecated in comments but may still be called from
   older simulation scripts.

#### 3.2c simple -- Issues

1. **Carryover half-life values contradict project
   convention.** The CLAUDE.md (inherited from main) states
   "Do NOT use the legacy values (0.5, 1.0, 2.0 weeks)"
   but `simulation.R` line 37 uses exactly
   `c(0, 0.5, 1, 2)`.

2. **Missing DESCRIPTION dependency.** The script requires
   `tidyverse` via `library(tidyverse)` but DESCRIPTION
   lists only `renv` and `testthat` as imports. This would
   fail R CMD check.

3. **Trivial test coverage.** Same `expect_true(TRUE)`
   placeholder as main.

4. **Placeholder author in DESCRIPTION.** The author field
   contains "Your Name" and "your.email@example.com" from
   the zzcollab template.

5. **`here::here()` used without `here` in Imports.**
   Works in the Docker container (via tidyverse transitive
   dependency) but is not formally declared.

6. **Magic number thresholds.** Minimum sample sizes for
   analysis (10 for OL+BDC, 5 for CO/N-of-1) are
   hard-coded without justification.

### 3.3 Documentation Quality

#### 3.3a orig

- **Completeness:** The vignettes describe the workflow
  (generate, visualize, apply to real data) but do not
  fully specify the DGP mathematics. The modified Gompertz,
  carryover formula, and correlation structure are
  documented only in code comments and roxygen headers.
- **Accuracy:** The code is the ground truth. The
  `scalefactor` parameter is documented as "TODO update
  when understand what this does?" (line 43 of
  generateData.R), indicating incomplete understanding
  even by the original author.
- **Accessibility:** A new collaborator would need to
  read the code to understand the simulation. The vignettes
  provide usage examples but not mathematical derivations.

#### 3.3b main

- **Completeness:** Extensive. Over 40 documents cover
  the mathematical foundations, carryover theory, sigma
  matrix derivation, positive definiteness constraints,
  Hendrickson alignment, and mixed model comparison. LaTeX
  documents provide formal mathematical treatment.
- **Accuracy:** Generally good, but the volume creates
  risk of contradictions. The critical finding about
  carryover half-life scale is well-documented in
  `CRITICAL_FINDING_CARRYOVER_HALFLIFE.md`.
- **Accessibility:** Mixed. The sheer volume (40+ files)
  is overwhelming. There is no single entry-point document
  that provides a complete mathematical specification.
  `docs/readme.md` exists but may not index all documents.
  A new collaborator could reproduce the simulation
  from documentation, but would need significant time to
  navigate the corpus.

#### 3.3c simple

- **Completeness:** The `simplification-plan.md` document
  thoroughly explains what was simplified and why. The
  `comparison_whitepaper.md` compares all three codebases.
  However, the DGP is documented only in code comments.
- **Accuracy:** The code is straightforward enough that
  the documentation matches. The carryover half-life
  mismatch with CLAUDE.md is the notable exception.
- **Accessibility:** Good. A single 562-line script with
  inline documentation is easy to follow. A new
  collaborator could understand the simulation quickly.

---

## 4. Evolution Analysis

### 4.1 Changes from orig to main

| Category              | Change                                            | Significance |
|-----------------------|---------------------------------------------------|:------------:|
| **Architecture**      | Package functions to research compendium scripts  | Major        |
| **Data structures**   | data.table to tibble/tidyverse                    | Moderate     |
| **Gompertz formula**  | Offset-rescale to complementary form              | Moderate     |
| **Autocorrelation**   | Compound symmetry to time-based AR(1)             | Major        |
| **Cross-correlation** | Constant to time-decaying (`0.9^lag`)             | Moderate     |
| **BM correlation**    | Ron Thomas mean-ratio adj. to binary on/off       | Major        |
| **Carryover models**  | Exponential only to Exp + Linear + Weibull        | Moderate     |
| **Carryover values**  | Unspecified to 0, 0.1, 0.2 weeks (corrected)     | Critical     |
| **Analysis engine**   | lmer only to lmer + lme with corCAR1              | Major        |
| **Two-scenario**      | Single model to WITH/WITHOUT carryover comparison | Major        |
| **Parallelism**       | Sequential to furrr/future                        | Moderate     |
| **Sigma validation**  | None to pre-simulation eigenvalue diagnostics     | Major        |
| **Sigma caching**     | Rebuilt per iteration to pre-computed cache        | Moderate     |
| **Dependencies**      | 10 packages to 38+ packages                       | Negative     |
| **Dead code**         | comb.R present to (mostly) removed                | Positive     |
| **eval(parse)**       | Present to eliminated                             | Positive     |
| **Docker/renv**       | Absent to present                                 | Major        |
| **Documentation**     | 4 files to 40+ files                              | Major        |

### 4.2 Simplifications from main to simple

| Simplification                        | Statistical validity preserved?        |
|---------------------------------------|----------------------------------------|
| MVN covariance to random intercept    | Partially. Loses temporal correlation  |
|                                       | structure; overstates effective N and   |
|                                       | thus overstates power.                 |
| Modified Gompertz to linear rate      | Partially. Valid for short trials where |
|                                       | saturation is negligible; misrepresents |
|                                       | dose-response shape for longer trials.  |
| Correlation-based BM interaction      | Different mechanism. The mean-based     |
| to mean-based interaction             | approach is more transparent but        |
|                                       | produces a different effect-size        |
|                                       | interpretation than the MVN approach.   |
| lmer/lme to lm on phase means        | Partially. Loss of temporal resolution  |
|                                       | reduces power but avoids convergence    |
|                                       | issues. The estimand changes from       |
|                                       | per-timepoint to phase-mean difference. |
| TR removed                            | Valid if natural progression is         |
|                                       | negligible over trial duration (~20 wk).|
| Outcome-dependent dropout to MCAR     | Partially. MCAR preserves unbiasedness  |
|                                       | but misses the bias introduced by       |
|                                       | informative dropout.                    |
| Carryover half-lives 0.1/0.2 to      | Invalid for comparison. The values      |
| 0.5/1/2 weeks                         | are 5-20x longer, producing much        |
|                                       | slower decay and larger carryover       |
|                                       | effects than Hendrickson.               |
| Scale factor 2 to 1 (implicit)        | Changes the effective decay rate.       |
|                                       | Combined with longer half-lives, this   |
|                                       | produces a categorically different       |
|                                       | carryover regime.                       |

### 4.3 Backport Candidates

| Improvement                            | From    | To       | Priority |
|----------------------------------------|---------|----------|:--------:|
| Pre-simulation sigma validation        | main    | orig     | High     |
| Time-based AR(1) autocorrelation       | main    | orig     | High     |
| `eval(parse())` elimination            | main    | orig     | High     |
| Docker + renv infrastructure           | main    | orig     | Medium   |
| Carryover half-life correction (0.1/2) | main    | simple   | Critical |
| CI/CD pipeline (GitHub Actions)        | simple  | main     | Medium   |
| Single-file pedagogical script         | simple  | main     | Low      |
| Weibull/linear carryover options       | main    | simple   | Low      |

### 4.4 Unresolved Divergences

#### 4.4.1 Gompertz formulation

The orig and main use different mathematical forms of the
modified Gompertz function with the same parameter names
(`max`, `disp`/`displacement`, `rate`). Unless the parameter
values are recalibrated, the two produce different response
trajectories for the same inputs. **The main's formulation
appears correct** as a standard complementary Gompertz, while
the orig's offset-rescale approach is non-standard and harder
to interpret.

#### 4.4.2 Biomarker-treatment interaction mechanism

The orig and main encode the interaction via differential
correlation in the MVN structure (c.bm > 0 on-drug,
reduced/zero off-drug). The simple encodes it as an additive
term in the mean structure (`br_rate * bm_mod * bm`). These
produce different estimands: the MVN approach creates a
correlation-mediated association that manifests as a
regression coefficient, while the mean-based approach creates
a direct moderation effect. **Both are defensible** but they
answer slightly different statistical questions and produce
non-comparable power estimates.

#### 4.4.3 Carryover half-life values

- main: 0, 0.1, 0.2 weeks (with scale factor 2)
- simple: 0, 0.5, 1, 2 weeks (with implicit scale factor 1)

Effective half-lives:
- main: 0.05, 0.10 weeks (~8h, ~17h)
- simple: 0.5, 1.0, 2.0 weeks (~84h, ~168h, ~336h)

The simple's carryover persists 10-20x longer than the
main's. **The main's values match Hendrickson et al. (2020)
and are correct for the intended application.** The simple's
values may be deliberately chosen to demonstrate carryover
effects more visibly, but they are not comparable to the
published methodology.

#### 4.4.4 Autocorrelation structure

- orig: Compound symmetry (constant correlation across lags)
- main: AR(1) with time-based decay
- simple: Independence (no temporal correlation)

**The main's AR(1) structure is most realistic** for
longitudinal biomarker data. The orig's compound symmetry
is an overly strong assumption. The simple's independence
assumption is appropriate for its pedagogical purpose but
would overstate power in realistic applications.

---

## 5. Recommendations

### 5.1 Statistical Correctness (Priority: Critical)

1. **Harmonize carryover half-life values in simple.**
   Replace `c(0, 0.5, 1, 2)` with `c(0, 0.1, 0.2)` to
   match Hendrickson et al. (2020) and the main codebase.
   If the longer values are retained for pedagogical
   purposes, they must be clearly labeled as
   non-Hendrickson exploratory values, and the scale
   factor difference must be documented.

2. **Verify Gompertz equivalence between orig and main.**
   Run both formulations with identical parameters and
   compare response trajectories. If they diverge, document
   which is authoritative and recalibrate parameters in
   the non-canonical version.

3. **Document the biomarker interaction mechanism
   difference.** The correlation-based (orig, main) and
   mean-based (simple) approaches are not equivalent. A
   brief technical note comparing the two mechanisms and
   their implications for power estimation would prevent
   inappropriate cross-codebase comparisons.

4. **Validate type I error rates.** None of the codebases
   includes a null-case validation (biomarker moderation =
   0, expected power = alpha). This is essential for
   confirming that the simulation and analysis pipeline
   does not produce inflated false-positive rates.

### 5.2 Software Engineering (Priority: High)

5. **Add unit tests for core simulation functions.** At
   minimum, test:
   - Gompertz function: known input/output pairs
   - Carryover decay: boundary conditions (t=0, t=inf,
     halflife=0)
   - Correlation matrix: symmetry, PD, correct dimensions
   - Analysis model: coefficient recovery under known
     effect sizes
   - Null case: type I error rate within binomial bounds

6. **Eliminate `comb.R` from orig.** This 1,571-line file
   is dead weight. Archive it outside the R/ directory or
   delete it entirely.

7. **Eliminate `eval(parse())` from orig.** Replace with
   programmatic data.table column operations using
   `.SDcols` or `set()`.

8. **Reduce main's DESCRIPTION dependencies.** Remove
   database drivers, palmerpenguins, skimr, naniar, visdat,
   and other packages not used in any simulation or
   analysis script.

9. **Consolidate simulation script variants in main.**
   Factor shared logic (design creation, parameter grid
   setup, model fitting, output formatting) into
   `pm_functions.R` and reduce the 13+ scripts to 2-3
   canonical entry points with configuration parameters.

10. **Fix simple's DESCRIPTION.** Add `tidyverse` (or
    specific tidyverse packages) to Imports. Replace
    placeholder author information.

### 5.3 Documentation (Priority: Medium)

11. **Create a canonical mathematical specification.**
    A single LaTeX document (or structured Rmd) that fully
    specifies the DGP, analysis model, and power
    calculation. This should be the authoritative reference
    that all three codebases implement, with a
    correspondence table mapping mathematical notation
    to variable names in each codebase.

12. **Create a documentation index for main.** The 40+
    documents need a clear hierarchy: which is
    authoritative, which is historical, which is
    exploratory. A single `docs/INDEX.md` with status
    annotations would suffice.

13. **Document the scale factor.** The orig's roxygen
    for `scalefactor` says "TODO update when understand
    what this does?" This parameter fundamentally changes
    the carryover regime and must be properly documented
    in all three codebases.

### 5.4 Consolidation Strategy (Priority: Strategic)

The three repositories serve distinct purposes and should
be maintained as separate artifacts, but with tighter
alignment:

- **orig**: Archive. Freeze as the published reference
  implementation for Hendrickson et al. (2020). Do not
  develop further. Apply only the minimal fixes needed
  for reproducibility (renv, Docker, elimination of
  `comb.R`).

- **main**: Primary research codebase. Continue active
  development. Consolidate simulation variants, expand
  test coverage, and produce a canonical mathematical
  specification document. This is the codebase that
  should be used for publication-quality results.

- **simple**: Pedagogical companion. Align carryover
  parameters with Hendrickson. Add a comparison mode
  that runs both mean-based and correlation-based
  biomarker interaction mechanisms side by side. This
  codebase should prioritize clarity and be usable as
  teaching material or a supplementary methods appendix.

A shared `pmsim-common` package containing only the
mathematical primitives (Gompertz, carryover decay,
correlation matrix builder) would reduce drift across
the three codebases, but the overhead of maintaining a
fourth repository may not be justified given the current
team size.

---

## Appendix A: Function Signature Cross-Reference

```
+---------------------------+--------------------------+-----------------------+
| orig                      | main                     | simple                |
+---------------------------+--------------------------+-----------------------+
| modgompertz(t, maxr,      | mod_gompertz(time,       | (linear: br_rate *    |
|   disp, rate)             |   max_value, disp, rate) |   (1 + bm_mod * bm)) |
+---------------------------+--------------------------+-----------------------+
| (inline in generateData)  | calculate_carryover(     | calc_carryover(       |
|                           |   tsd, prev_effect,      |   weeks_since_drug,   |
|                           |   model, params)         |   halflife)           |
+---------------------------+--------------------------+-----------------------+
| generateData(modelparam,  | generate_data(           | generate_participant( |
|   respparam, blparam,     |   model_param, resp_param|   id, design,         |
|   trialdesign, empirical, |   baseline_param,        |   biomarker, ...)     |
|   makePositiveDefinite,   |   trial_design,          |                       |
|   seed, scalefactor)      |   empirical, make_pd,    |                       |
|                           |   seed, scale_factor)    |                       |
+---------------------------+--------------------------+-----------------------+
| buildtrialdesign(name,    | create_*_design(         | create_design(        |
|   timepoints, expectancy, |   n_participants,        |   design_name)        |
|   ondrug)                 |   measurement_weeks)     |                       |
+---------------------------+--------------------------+-----------------------+
| lme_analysis(td_set,      | lme_analysis(td_set,     | analyze_trial(        |
|   dat, op)                |   data, options)         |   trial_data,         |
|                           |                          |   design_name,        |
|                           |                          |   adjust_carryover)   |
+---------------------------+--------------------------+-----------------------+
| censordata(dat, td,       | (in simulation scripts)  | (inline in            |
|   censorparam)            |                          |   generate_participant|
+---------------------------+--------------------------+-----------------------+
| generateSimulatedResults( | (in simulation scripts)  | (inline main loop)    |
|   trialdesigns, resp,     |                          |                       |
|   bl, censor, model,      |                          |                       |
|   sim, analysis)          |                          |                       |
+---------------------------+--------------------------+-----------------------+
```

## Appendix B: Parameter Defaults

```
+----------------------------+----------+----------+----------+
| Parameter                  | orig     | main     | simple   |
+----------------------------+----------+----------+----------+
| N (participants)           | external | 35, 70   | 35, 70   |
| Carryover half-life (wk)   | external | 0/0.1/0.2| 0/.5/1/2 |
| Scale factor               | 2        | 2        | 1        |
| c.bm (BM correlation)      | external | 0.3      | N/A      |
| c.br (BR autocorrelation)  | external | 0.75     | N/A      |
| c.pb (ER autocorrelation)  | external | 0.75     | N/A      |
| c.tv (TR autocorrelation)  | external | 0.75     | N/A      |
| c.cf1t (cross-factor)      | external | 0.12     | N/A      |
| c.cfct (cross-time)        | external | 0.05     | N/A      |
| sigma_between              | external | (in MVN) | 1.0      |
| sigma_within               | external | (in MVN) | 0.5      |
| Iterations                 | external | 500      | 500      |
| alpha                      | 0.05     | 0.05     | 0.05     |
| Designs                    | flexible | 4 fixed  | 3 fixed  |
+----------------------------+----------+----------+----------+
```

---

*Rendered on 2026-03-18 at 08:55 PDT.*
*Source: ~/prj/alz/pmsim-audit.md*
