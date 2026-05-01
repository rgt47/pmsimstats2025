# pmsimstats2025: Comprehensive Status and Path Forward (Audit Edition)

**Replicating Hendrickson et al. (2020) with modern practices
and addressing the carryover-power problem. Updated with
comprehensive technical audit of all three repository variants.**

2026-03-18 (Audit: Phase 1-5 Complete)
2026-03-19 (Validation Experiment: Phase 1 Complete - Baseline Replication Successful)

---

## 1. Project Goal

pmsimstats2025 aims to:

1. **Faithfully replicate** the Hendrickson et al. (2020)
   N-of-1 trial simulation using modern R practices
   (tidyverse, parallel processing, Docker/renv,
   pre-validated sigma matrices)
2. **Diagnose** why power falls as carryover half-life
   increases in the original simulation
3. **Demonstrate** that explicitly modeling carryover in the
   analysis recovers lost power (the two-scenario
   comparison)

The project maintains three related repositories with
different purposes and maturity levels:

- **orig** (pmsimstats): Published package (Hendrickson et al. 2020);
  stable, archive-quality codebase
- **main** (pmsimstats2025): Modernized research compendium; active
  development with 40+ analysis scripts
- **simple** (pmsimstats-simple): Pedagogical simplification; minimal
  dependencies for educational use

---

## 2. What Works

pmsimstats2025 has made substantial architectural
improvements over the original:

- **Reproducible environment**: Docker + renv + Makefile
  (the original had none)
- **Parallel processing**: `furrr::future_map_dfr()` across
  iterations (the original was sequential)
- **Sigma caching**: Pre-built and validated sigma matrices
  avoid redundant computation
- **Two-stage MVN generation**: Conditional draw using
  pre-computed `Sigma_22_inv` and `cond_mean_transform`
- **Pre-simulation validation**: Eigenvalue checking and
  condition number reporting for all sigma structures
- **Two-scenario comparison framework**: Runs each
  condition with and without carryover in the analysis
  model, directly quantifying the cost of ignoring
  carryover
- **Extensive documentation**: 40+ documents covering
  mathematical foundations, Hendrickson alignment,
  carryover theory, and implementation decisions

The binary treatment assignment vectors (the `ondrug`
paths) match exactly between the original and
pmsimstats2025 for all four designs (OL, OL+BDC,
Crossover, Hybrid). The expectancy parameters also match.

---

## 3. Structural and Architectural Audit

A comprehensive Phase 1-5 audit examined all three repositories across
structural inventory, statistical models, code quality, evolution,
and recommendations. Key findings:

### 3.0 Repository Maturity Matrix

| Aspect | orig | main (2025) | simple |
|--------|------|-----------|--------|
| **Status** | Published (2020) | Active (2025) | In development |
| **R functions** | 11 in /R/ | 1 in /R/, 43 in scripts/ | 0 in /R/, 1 in scripts/ |
| **Test coverage** | 0 files | 5 files (14 assertions) | 1 file (placeholder) |
| **Documentation** | 1 README | 32 technical docs | 1 plan |
| **Docker+renv** | No | Yes | Yes (minimal) |
| **Parallelism** | Sequential | Via furrr/future | None |
| **Code duplication** | None | ~40% across scripts | N/A (single script) |
| **Reproducibility** | Via vignettes | Excellent (Docker/renv) | Basic |

### 3.1 Code Quality Issues

**orig**:
- Dead code: `scratch_WorkingScript.R` (low severity)
- Fragile patterns: `eval(parse(text=))` in generateData.R (medium)
- Ron Thomas correlation logic is opaque and hard to validate (medium)
- Zero automated tests (high)
- Magic numbers hard-coded: `scalefactor=2`, `tol=1e-3` (low)

**main**:
- Code duplication across simulation_*.R scripts (~40% redundancy) (high)
- Hard-coded parameter grids in 10+ files (medium)
- Mixed correlation parameterization (Hendrickson vs. reduced values) (medium)
- Incomplete test coverage (scripts not tested end-to-end) (medium)
- Sigma caching not thread-safe in parallel context (low)

**simple**:
- Placeholder tests (low)
- No error handling (low, acceptable for pedagogy)
- Design specs in function body, not config (low)

### 3.2 Architectural Strengths

**orig**:
- Clear API (exported functions via roxygen2)
- Each function has single responsibility
- Documented via roxygen comments

**main**:
- Modern stack (tidyverse, future/furrr)
- Reproducible (Docker, renv, Makefile)
- Comprehensive documentation
- Clear separation of exploratory vs. production runs

**simple**:
- Maximally readable (single monolithic script)
- Excellent for pedagogical use
- Easy to modify and experiment

### 3.3 Dependency Footprint

Critical dependencies used by all three:
- MASS (mvrnorm)
- ggplot2 (visualization)
- lme4 (mixed-effects models)

Modern additions (main only):
- tidyverse stack (dplyr, tidyr, tibble)
- future/furrr (parallelism)
- renv (reproducibility)

Legacy packages (orig only):
- svMisc, tictoc, ggpubr, gridExtra (could be removed)

---

## 4. Critical Divergences Identified

The code review identified five critical differences
between pmsimstats2025 and the original that produce
non-comparable power estimates. These are listed in order
of impact.

### 3.1 Biomarker-Treatment Interaction Mechanism

**The single most consequential difference.**

The original encodes the biomarker-treatment interaction
entirely within the MVN covariance matrix. There is no
biomarker term in the mean structure. The mechanism:

- On-drug timepoints: `corr(bm, BR) = c.bm`
- Off-drug timepoints: `corr(bm, BR)` scaled by
  `mm1/mm0` (Ron Thomas logic) or set to 0
- Participants with higher biomarker values tend to have
  higher BR draws purely through multivariate correlation
- The interaction emerges stochastically from `mvrnorm()`

pmsimstats2025 encodes the interaction in the mean
structure:

```r
effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered)
```

This is a deterministic modulation: high-biomarker
participants get a larger BR mean regardless of the
covariance structure. The sigma matrix uses a constant
`c.bm` across all timepoints, with no treatment-dependent
adjustment.

Prior documentation (`why-post-mvn-adjustment-analysis.md`,
`hendrickson-vs-our-detailed-comparison.md`) identified
that an earlier version of pmsimstats2025 applied the
biomarker effect three times (correlation + population
mean shift + post-MVN adjustment), double-counting the
interaction. The recommendation to use "Option 1: pure
correlation" was made but not yet fully implemented in
`simulation_clustered.R`.

**Impact**: The mean-based mechanism produces stronger,
more detectable interactions per participant than the
correlation-based mechanism. Power estimates are
systematically higher in pmsimstats2025 than the original
for the same nominal parameter values.

### 3.2 Per-Path vs. Shared Sigma Generation

The original generates data **separately for each trial
path** using a **path-specific sigma matrix**. Each path's
covariance structure depends on its treatment sequence
because the Ron Thomas correlation scaling adjusts
biomarker-BR correlations based on whether each timepoint
is on-drug or off-drug for that specific path.

pmsimstats2025 generates all participants from a **single
cached sigma matrix per design**. Path-specific treatment
effects are applied after the MVN draw, in the mean
structure. The covariance structure is identical for all
paths.

**Impact**: In the original, Path A (drug-first) and
Path B (placebo-first) in a crossover have different
covariance matrices. The biomarker correlates with BR at
different timepoints depending on path. In pmsimstats2025,
the biomarker correlates identically with BR regardless
of path.

### 3.3 Time Since Discontinuation (`tsd`)

The original computes `tsd` as cumulative off-drug
interval widths:

```r
tsd[i] = t_wk[i] * (od[i] != 1)
# then accumulate: tsd[i] += tsd[i-1] when off-drug
```

At the first off-drug timepoint, `tsd` equals the
measurement interval duration (1-4 weeks).

pmsimstats2025 computes `tsd` as calendar time since
transition:

```r
tsd = week - discontinuation_week
```

At the first off-drug timepoint, `tsd = 0` (same week as
the transition event). This produces a systematic offset
of one measurement interval at every off-drug timepoint.

Detailed path-by-path comparison (`tsd_comparison.md`):

```
Hybrid Path A:
  Week:    4   8   9  10  11  12  16  20
  orig:    0   0   0   0   1   2   0   4
  main:    0   0   0   0   0   1   0   0

CO Path A:
  Week:  2.5   5  7.5  10  12.5  15  17.5  20
  orig:    0   0    0   0   2.5  5.0  7.5  10.0
  main:    0   0    0   0   0    2.5  5.0   7.5
```

**Impact**: At the first off-drug measurement,
pmsimstats2025 computes `carryover_effect = (1/2)^0 = 1`
(full retention), while the original computes
`(1/2)^(2*1/0.2) = 0.001` (0.1% retention). This is a
1000x difference in carryover at the most informative
contrast point. The `predictor_comparison.tex` and
`tsd_comparison.md` documents thoroughly characterize
this as a step-function power cliff.

### 3.4 Carryover Scale Factor

The original uses `scalefactor = 2` in data generation:

```r
brmeans[p] += brmeans[p-1] * (1/2)^(2 * tsd / t1half)
```

pmsimstats2025's `simulation_clustered.R` omits the scale
factor:

```r
carryover_effect = (1/2)^(tsd / carryover_t_half)
```

(Note: `pm_functions.R` retains `scale_factor = 2` as
default, but the primary simulation script does not use
it.)

**Impact**: The effective half-life is doubled in
pmsimstats2025 relative to the original for the same
nominal parameter. Combined with the tsd offset, this
compounds to produce dramatically different carryover
behavior.

Additionally, `hendrickson-carryover-analysis.md`
discovered that the original uses `scalefactor = 2` in
data generation but `carryover_scalefactor = 1` in the
analysis model. This DGP-analysis mismatch is a feature
of Hendrickson's design (carryover is intentionally
unmodeled in the analysis), but the different scale
factors add an additional layer of misspecification.

### 3.5 Cumulative Drug Exposure (`tod` / `weeks_on_drug`)

The original computes `tod` as cumulative time-weighted
exposure in weeks, with gap reset:

```
Hybrid Path A:
  Week:  4   8   9  10  11  12  16  20
  tod:   4   8   9  10   0   0   4   0
```

At week 16 (drug re-initiated after off-drug gap),
`tod = 4` (just the new 4-week interval), not 14 (total
lifetime exposure). The Gompertz restarts from its steep
early phase.

pmsimstats2025 uses `cumsum(treatment)` (a count of
on-drug measurement occasions) that plateaus off-drug
and continues from the accumulated count when drug
resumes:

```
Hybrid Path A:
  Week:  4   8   9  10  11  12  16  20
  wod:   1   2   3   4   4   4   5   5
```

**Impact**: Different units (weeks vs. count), different
off-drug behavior (0 vs. plateau), and different
re-initiation behavior (reset vs. continue). These feed
into different growth models (Gompertz vs. linear),
further amplifying the divergence.

---

## 4. Additional Divergences

These are less impactful but should be addressed for full
fidelity.

### 4.1 Response Growth Model

The original uses a modified Gompertz function (sigmoidal,
saturating) with four parameters per factor (max,
displacement, rate, sd). pmsimstats2025 uses a linear rate
model (unbounded, one parameter per factor). For short
trials with small rate parameters, these may be
approximately equivalent, but they produce different
dose-response shapes.

### 4.2 Autocorrelation Structure

The original uses compound symmetry within each factor
(constant correlation regardless of time lag).
pmsimstats2025 uses time-based AR(1) where correlation
decays with elapsed time: `rho^|t_i - t_j|`. This is a
deliberate improvement in realism.

### 4.3 ER SD Scaling

The original scales expectancy response standard
deviations by the expectancy parameter (blinded timepoints
have halved ER variance). pmsimstats2025 does not scale
ER SDs; expectancy enters only through the mean.

### 4.4 BM-BR Correlation at Timepoint 1

The original's Ron Thomas loop starts at `p > 1`, so the
first timepoint's BR has zero correlation with the
biomarker even when on-drug. pmsimstats2025 applies `c.bm`
to all timepoints including the first.

### 4.5 BM-ER and BM-TR Correlations

The original correlates the biomarker only with BR.
pmsimstats2025 also correlates the biomarker with ER and
TR at half-strength (`c.bm * 0.5`). The original has no
such cross-component biomarker correlations.

### 4.6 Correlation Parameter Values

The original uses `c.br = c.pb = c.tv = 0.8`,
`c.cf1t = 0.2`, `c.cfct = 0.1`. pmsimstats2025 reduced
these to `0.75`, `0.12`, and `0.05` for improved positive
definiteness conditioning.

---

## 5. The Power Dropoff Problem

### 5.1 The observation

Hendrickson's published results show power declining as
carryover half-life increases from 0 to 0.1 to 0.2 weeks.
The two-scenario comparison in pmsimstats2025 is designed
to demonstrate that this decline can be mitigated by
explicitly modeling carryover in the analysis.

### 5.2 Why does power drop?

Hendrickson's analysis model does not include carryover.
The `lme_analysis.R` defaults to
`op$carryover_t1half = 0`, which makes `Dbc` a binary
on/off indicator. When carryover exists in the data but
not in the model:

1. Off-drug observations carry residual drug effect
2. The binary `Dbc` codes them as 0 (no drug)
3. The mismatch between actual and modeled drug exposure
   inflates residual variance
4. Standard errors on the interaction coefficient increase
5. Power declines

### 5.3 Scale factor compounds the problem

The data generation uses `scalefactor = 2` but the
analysis uses `scalefactor = 1`. Even if the analysis
included carryover, the decay rates would not match. At
`tsd = 1, t1half = 0.2`:

- DGP retention: `(1/2)^(2*1/0.2) = (1/2)^10 = 0.001`
- Analysis retention: `(1/2)^(1*1/0.2) = (1/2)^5 = 0.031`
- 30x mismatch

This is a separate, unmodeled confound on top of the
binary vs. continuous `Dbc` issue.

### 5.4 Candidate mechanisms (post-tsd-fix)

After fixing the tsd computation to match the original,
if excessive power dropoff persists, three mechanisms
should be investigated:

1. **Scale factor mismatch**: Testable by running with
   matched scale factors (both 1 or both 2)
2. **PD enforcement distortion**: Small nonzero
   correlations at off-drug timepoints may cause
   `make.positive.definite()` to perturb the on-drug
   correlations differently than exact-zero correlations
3. **Cascading carryover in BR means**: The recursive
   formula `brmeans[p] += brmeans[p-1] * decay` means
   each off-drug timepoint's carryover depends on the
   previous, creating nonlinear contamination patterns

---

## 6. The Two-Scenario Comparison

This is pmsimstats2025's primary methodological
contribution.

### 6.1 Design

For each parameter combination:

- **Scenario 1** (WITH carryover modeling): Analysis uses
  `effective_drug_weeks` which accounts for carryover
  decay. The analysis model matches the DGP.
- **Scenario 2** (WITHOUT carryover modeling): Analysis
  uses `weeks_on_drug` (binary, no decay). This
  replicates Hendrickson's approach and should reproduce
  the published power decline.

### 6.2 Expected outcome

```
carryover = 0: Scenario 1 = Scenario 2 (no confound)
carryover > 0: Scenario 1 > Scenario 2 (modeling helps)
```

The gap between the two scenarios quantifies the
recoverable power -- the cost of ignoring carryover that
Hendrickson's approach pays.

### 6.3 Current status

The framework is implemented in `simulation_clustered.R`
but cannot produce valid results until the five critical
divergences (Section 3) are resolved. The comparison is
meaningful only when Scenario 2 reproduces the original's
power curves.

---

## 7. Changes Required

### Phase 1: Faithful Replication (Critical)

These changes bring pmsimstats2025's DGP into alignment
with the original. Until these are complete, power
estimates are not comparable.

**Change 1: Correlation-based biomarker interaction.**
Remove the mean-structure biomarker modulation
(`BR_rate * (1 + bm_mod * bm_centered)`) and implement
the Ron Thomas differential correlation logic in the
sigma matrix. The biomarker-BR correlation must be
`c.bm` at on-drug timepoints and scaled/zero at off-drug
timepoints. This requires:

- Building path-specific sigma matrices (Change 2)
- Computing BR means via Gompertz of `tod` (Change 5)
- Applying the Ron Thomas `mm1/mm0` scaling

**Change 2: Per-path sigma generation.**
Replace the single cached sigma per design with
path-specific sigma matrices. Each path's sigma must
reflect its treatment sequence in the biomarker-BR
correlation block. This may reduce the effectiveness of
sigma caching (more structures to pre-build) but is
essential for the correlation-based interaction to work.

**Change 3: Fix tsd computation.**
Replace the calendar-subtraction approach with
interval-based accumulation matching the original:

```r
last_on_drug_week = if_else(
  treatment == 1, week, NA_real_)
last_on_drug_week = zoo::na.locf(
  last_on_drug_week, na.rm = FALSE)
tsd = if_else(
  treatment == 0 & !is.na(last_on_drug_week),
  week - last_on_drug_week, 0)
```

This anchors tsd to the last on-drug measurement week,
producing `tsd > 0` at the first off-drug timepoint.

**Change 4: Restore scalefactor = 2.**
Add `scalefactor = 2` to the carryover formula in
`simulation_clustered.R`:

```r
carryover_effect = (1/2)^(2 * tsd / carryover_t_half)
```

**Change 5: Time-weighted tod with gap reset.**
Replace `cumsum(treatment)` with time-weighted
accumulation that resets after off-drug gaps:

```r
interval = week - lag(week, default = 0)
drug_interval = treatment * interval
# Accumulate only within contiguous on-drug blocks
```

### Phase 2: Growth Model Alignment (High)

**Change 6: Implement modified Gompertz.**
Replace the linear rate model with the original's
offset-rescaled Gompertz:

```r
modgompertz <- function(t, maxr, disp, rate) {
  y <- maxr * exp(-disp * exp(-rate * t))
  vert_offset <- maxr * exp(-disp * exp(-rate * 0))
  y <- y - vert_offset
  y <- y * (maxr / (maxr - vert_offset))
  y
}
```

This requires specifying `max`, `disp`, `rate`, and `sd`
parameters for each of the three factors (BR, ER, TR),
matching the original's `respparam` data table.

### Phase 3: Covariance Refinements (Moderate)

**Change 7: Scale ER SDs by expectancy.**
Multiply ER standard deviations by the expectancy
parameter at each timepoint (0.5 for blinded, 1.0 for
open-label).

**Change 8: Remove BM-ER and BM-TR correlations.**
Set `Sigma_12[er_idx, 1] = 0` and
`Sigma_12[tr_idx, 1] = 0` to match the original, which
correlates the biomarker only with BR.

**Change 9: Skip BM-BR correlation at timepoint 1.**
Match the original's `p > 1` guard by setting the first
timepoint's biomarker-BR correlation to zero.

### Phase 4: Validation and Diagnosis (High)

**Change 10: Validation script.**
Write a script that computes and prints `treatment`,
`tod`, `tsd`, `carryover_effect`, and the biomarker-BR
correlations for each design/path, comparing against
hand-computed reference values from the original.

**Change 11: Diagnostic script for power dropoff.**
For each carryover level (0, 0.1, 0.2):

- Print the full sigma matrix (or at minimum the
  biomarker-BR correlation subvector)
- Print eigenvalues before and after PD enforcement
- Compare element-wise between `t_half = 0` and
  `t_half > 0`

**Change 12: Regression tests.**
Add testthat tests for:

- `tsd` values per design/path against reference
- Carryover at boundary conditions
- Null case (biomarker_moderation = 0) rejection rate
  within binomial bounds of 5%

### Phase 5: Extension (the Two-Scenario Payoff)

**Change 13: Run baseline comparison.**
Once Phase 1 is complete, run Scenario 2 (without
carryover modeling) and compare power curves against
Hendrickson's published heatmaps. This validates the
replication.

**Change 14: Run two-scenario comparison.**
Run both scenarios across the full parameter grid. The
gap between Scenario 1 and Scenario 2 is the paper's
main result: the recoverable power from explicitly
modeling carryover.

**Change 15: Scale factor investigation.**
Run with matched scale factors (both 1, both 2) and
unmatched (DGP=2, analysis=1) to isolate the scale
factor's contribution to power dropoff.

---

## 8. Architectural Decisions to Preserve

These pmsimstats2025 improvements should be kept even
when aligning with the original's statistical methodology:

- **Docker + renv + Makefile** for reproducibility
- **furrr parallelism** for runtime efficiency
- **Pre-simulation sigma validation** (eigenvalue
  checking) -- extend to validate per-path sigmas
- **Two-stage conditional MVN** (mathematically equivalent
  to joint draw, computationally faster)
- **Timestamped output filenames** and **provenance
  annotations** on figures
- **Structured logging** to file

The AR(1) autocorrelation structure is a deliberate
improvement but its impact on power should be quantified.
Consider running a sensitivity analysis with compound
symmetry (matching the original) vs. AR(1) to isolate
this effect.

---

## 9. Open Questions

### 9.1 Does fixing tsd eliminate excessive power dropoff?

After implementing the interval-based tsd (Phase 1,
Change 3), re-run Scenario 2 and compare against
Hendrickson's published results. If the power curves
match, the tsd bug was the primary cause. If dropoff
persists, investigate the scale factor mismatch
(Change 15) and PD enforcement distortion (Change 11).

### 9.2 Is per-path sigma computationally feasible?

The original builds one sigma per path per parameter
combination per iteration. With 4 paths and 8 sigma
structures (4 designs x 2 c.bm levels), this becomes
32 sigma matrices. Pre-building and caching all 32 is
straightforward but requires the Gompertz means (which
depend on the response parameters, not just the design
geometry). Evaluate whether sigma caching can still work
or whether sigma must be rebuilt per parameter set.

### 9.3 Should tsd_min be a design feature?

The `predictor_comparison.tex` discusses a `tsd_min`
workaround (setting a minimum tsd of 0.5 weeks). This is
pharmacologically defensible: real trials always have a
gap between last dose and assessment. If empirical
carryover estimation is pursued, `tsd_min` could become
a deliberate design parameter rather than a bug fix. But
for replication purposes, it should match the original's
interval-based values exactly.

### 9.4 What is the actual carryover in prazosin trials?

The simulation framework can be applied to real clinical
data to estimate actual carryover half-lives. The
SUMMARY_WHITE_PAPER.md outlines a research plan for this
(Section 12). Empirical estimation would ground the
simulation parameters in clinical reality rather than
the current hypothetical values.

---

## 10. Document Index

Documents produced during the March 2026 code review:

| Document | Scope |
|----------|-------|
| `pmsim-audit.md` | Structural inventory of all three repos (Phase 1-5 audit) |
| `pmsim-power-dropoff-analysis.md` | Variable chain tracing (V1-V5), tsd comparison, recommendations |
| `pmsim-orig-vs-main-deep-review.md` | Line-by-line DGP comparison (13 differences catalogued) |

Key existing documents:

| Document | Scope |
|----------|-------|
| `predictor_comparison.tex` | tsd bug analysis with numerical tables |
| `tsd_comparison.md` | Path-by-path tsd values, orig vs. 2025 |
| `hendrickson-carryover-analysis.md` | Discovery: Hendrickson does not model carryover in analysis |
| `hendrickson-vs-our-detailed-comparison.md` | Post-MVN adjustment analysis, Option 1 recommendation |
| `why-post-mvn-adjustment-analysis.md` | Correlation vs. additive interaction mechanisms |
| `pseudocode-comparison-hendrickson-vs-ours.md` | Side-by-side pseudocode (9 difference categories) |
| `two-scenario-comparison.md` | WITH/WITHOUT carryover modeling experimental design |
| `SUMMARY_WHITE_PAPER.md` | High-level framework and research plan |
| `simulation_white_paper.md` | Methods paper describing Hendrickson alignment |

---

## 11. Summary

pmsimstats2025 has a sound architectural foundation and a
well-designed two-scenario comparison framework. However,
its data-generating process diverges from Hendrickson's
original in five critical ways that prevent meaningful
power comparisons:

1. Mean-based vs. correlation-based biomarker interaction
2. Shared vs. per-path sigma matrices
3. Calendar-based vs. interval-based tsd
4. Missing scale factor in carryover formula
5. Count-based vs. time-weighted cumulative drug exposure

These divergences compound: the tsd offset produces 1000x
differences in carryover retention, the interaction
mechanism produces categorically different statistical
signals, and the shared sigma eliminates the
treatment-dependent covariance structure that is central
to Hendrickson's approach.

The path forward is clear. Phase 1 (Changes 1-5) achieves
faithful replication by aligning the DGP with the original.
Phase 4 (Changes 10-12) validates the alignment. Phase 5
(Changes 13-15) delivers the scientific payoff: a
quantified demonstration that explicit carryover modeling
recovers power lost to unmodeled confounding in N-of-1
trial designs.

---

## 12. Consolidated Recommendations (from Audit Phase 5)

### Priority 1: Resolve Statistical Divergences

**1.1 Biomarker-interaction mechanism** (CRITICAL)
- **Action**: Perform power comparison study using orig's
  correlation-based method vs. main's mean-based method
- **Owner**: Rebecca Hendrickson (clarify Hendrickson 2020 intent)
- **Timeline**: 2-3 weeks
- **Impact**: Determines whether orig or main is correct reference

**1.2 Carryover scale factor** (CRITICAL)
- **Action**: Implement both scale_factor=1 and scale_factor=2
  in main; run sensitivity analysis comparing power curves
- **Owner**: Ronald Thomas
- **Timeline**: 1-2 weeks
- **Impact**: 3-order-of-magnitude differences in carryover strength

**1.3 Fix tsd computation** (HIGH)
- **Action**: Replace calendar-based with interval-based accumulation
  matching original (Change 3 in Section 7)
- **Timeline**: 3-5 days
- **Impact**: Addresses 1000x carryover cliff at first off-drug timepoint

### Priority 2: Software Engineering

**2.1 orig: Eliminate eval/parse fragility**
- **Action**: Refactor generateData.R lines 207, 212 to use
  explicit data.table syntax
- **Timeline**: 2-3 days
- **Impact**: Easier maintenance; fewer eval/parse bugs

**2.2 main: Reduce code duplication**
- **Action**: Extract `run_simulation()` function in pm_functions.R
- **Timeline**: 2-3 weeks
- **Impact**: 50% reduction in codebase; unified bug fixes

**2.3 main: Add automated tests**
- **Action**: Port test-hendrickson-dgp.R patterns to end-to-end
  power curve validation against Hendrickson 2020 figures
- **Timeline**: 1-2 weeks per design
- **Impact**: Catch regressions; validate against publication

### Priority 3: Consolidation Strategy

**Recommendation: MAINTAIN SEPARATE** with clear purpose statements:

- **orig**: Archive-quality (codefreeze); reproduce publication exactly
- **main**: Active development (breaking changes OK if documented)
- **simple**: Educational variant (minimal, stable)

**Cross-repo coordination**:
- Significant bugs discovered in main → report to orig if publication-
  critical
- Improvements in main → evaluate backport to orig (especially tests)
- New features in simple → document as applicable in main

### Priority 4: Parameter Standardization

**4.1 Carryover half-life values**
- **Project spec** (CLAUDE.md): `c(0, 0.1, 0.2)` weeks (Hendrickson)
- **Action**: Audit all scripts for legacy values `c(0, 0.5, 1.0, 2.0)`
- **Timeline**: 1-2 days
- **Impact**: Consistency across repos; comparable results

### Priority 5: Performance Optimization

**5.1 Pre-compute sigma matrices** (2-3x speedup)
- Current: Computed inside parallel loop
- Target: Pre-compute once per design before parallelization
- **Timeline**: 1 week

**5.2 Cache LME model fits** (1.2-1.5x speedup)
- Current: Refit identical models
- Target: Cache and reuse
- **Timeline**: 3-5 days

---

## 13. Decision Matrix: What Happens Next

| Issue | Status | Decision Required | Owner | By When |
|-------|--------|-------------------|-------|---------|
| **Interaction mechanism** | UNRESOLVED | Correlation vs. mean-based | RH | 2026-04-01 |
| **Scale factor** | UNRESOLVED | Use 1, 2, or sensitivity | RGT | 2026-03-31 |
| **tsd computation** | READY | Implement interval-based | Any | 2026-03-25 |
| **Reproduce Hendrickson power curves** | READY | Run Scenario 2 validation | Any | 2026-04-01 |
| **Two-scenario comparison** | BLOCKED | Until (3) resolved | Any | 2026-04-15 |
| **Repository governance** | READY | Maintain separate + coordinate | RH/RGT | 2026-03-25 |

---

*Rendered on 2026-03-18 at 21:40 PDT (Audit Update).*
*Source: ~/prj/alz/01-pmsimstats/pmsimstats2025/docs/pmsimstats2025-comprehensive-status.md*
