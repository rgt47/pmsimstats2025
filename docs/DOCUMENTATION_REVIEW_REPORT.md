# Documentation Review Report: pmsimstats2025

## Comprehensive Assessment of Documentation Accuracy, Consistency, and Publication Readiness

**Review Date:** February 2026
**Reviewer:** Claude Code Analysis
**Project:** N-of-1 Trial Simulation for Biomarker-Treatment Interaction Detection

---

## Executive Summary

This report evaluates the pmsimstats2025 repository documentation against current
best practices for N-of-1 trial simulation, carryover effect modeling, and
biomarker-treatment interaction detection. The analysis identifies areas of
strength, inconsistencies requiring correction, and enhancements needed for
publication readiness.

**Overall Assessment:** The project has a strong methodological foundation but
requires documentation updates to reflect recent code changes (particularly the
switch from `lmer` to `lme` for proper AR(1) correlation) and alignment with
contemporary best practices.

---

## 1. Documentation Inventory

### 1.1 Primary Documentation Files Reviewed

| File | Status | Priority |
|:-----|:-------|:---------|
| `README.md` | Needs Update | High |
| `docs/simulation_white_paper.md` | Needs Update | High |
| `vignettes/quickstart.Rmd` | Needs Update | Medium |
| `analysis/docs/carryover_variance_whitepaper.tex` | Current | Low |
| `analysis/scripts/README.md` | Needs Review | Medium |

### 1.2 Documentation Statistics

- Total markdown files: 22
- Total LaTeX documents: 10
- Total vignettes: 4
- Total PDF outputs: 16

---

## 2. Critical Inconsistencies Identified

### 2.1 Model Engine Discrepancy (HIGH PRIORITY)

**Issue:** Documentation references `lmer` (lme4) throughout, but the current
implementation uses `lme` (nlme) with `corCAR1` for proper AR(1) correlation.

**Affected Files:**

- `README.md` (lines 27, 154-162): References `lmer` models
- `docs/simulation_white_paper.md` (lines 237-251): Shows `lmer` code examples
- `vignettes/quickstart.Rmd` (lines 131-147): Uses `lmer` examples

**Current Implementation** (`simulation_clustered.R`, line 54):
```r
model_engine <- "lme"  # Uses nlme::lme with corCAR1
```

**Correction Required:** All documentation should reference `lme` with `corCAR1`
for continuous-time AR(1) correlation structure that matches the DGP.

**Why This Matters:** The web consensus from [Schmid & Yang (2014)](https://effectivehealthcare.ahrq.gov/products/n-1-trials/research-2014-1)
emphasizes that "models that adjust for autocorrelation take two main forms" and
proper specification is critical for valid inference. The `lme4::lmer` package
does not support continuous-time AR(1) correlation via `corCAR1`; only `nlme::lme`
does.

### 2.2 Correlation Parameter Values (MEDIUM PRIORITY)

**Issue:** Documentation cites Hendrickson's original values, but implementation
uses reduced values for numerical stability.

**Documentation States** (`docs/simulation_white_paper.md`, lines 73-78):
```
c.tv = c.pb = c.br = 0.8
c.cf1t = 0.2
c.cfct = 0.1
```

**Current Implementation** (`simulation_clustered.R`, lines 70-76):
```r
c.br <- 0.75  # Reduced from 0.8
c.er <- 0.75  # Reduced from 0.8
c.tr <- 0.75  # Reduced from 0.8
c.cf1t <- 0.12  # Reduced from 0.2
c.cfct <- 0.05  # Reduced from 0.1
```

**Correction Required:** Document the rationale for parameter reduction
(positive definiteness constraints) and update stated values.

### 2.3 Carryover Parameterization (MEDIUM PRIORITY)

**Issue:** Documentation describes carryover as `carryover_decay_rate` (proportion),
but implementation uses `carryover_halflife` (half-life in weeks).

**Documentation** (`docs/simulation_white_paper.md`, lines 284-285):
```
Carryover decay rate | 0, 0.5 | Proportion of effect persisting
```

**Current Implementation** (`simulation_clustered.R`, lines 95-96):
```r
carryover_mode <- "halflife"  # Not "proportion"
carryover_halflife = c(0, 0.5, 1.0, 2.0)  # Half-lives in weeks
```

**Correction Required:** Update documentation to use half-life parameterization
consistently, which is more interpretable pharmacologically.

### 2.4 Trial Design Coverage (LOW PRIORITY)

**Issue:** README mentions 5 designs, but simulation_clustered.R only implements
Hybrid and OL+BDC.

**Documentation** (`README.md`, lines 19-20):
```
Multiple Trial Designs: Hybrid, OL+BDC, Open-Label, Crossover, Parallel
```

**Current State:** `simulation_clustered.R` implements only Hybrid and OL+BDC.
Other designs are in `simulation_evenly_spaced.R`.

**Correction Required:** Clarify which script implements which designs.

---

## 3. Alignment with Contemporary Best Practices

### 3.1 Carryover Effect Modeling

**Web Consensus:** According to [Wienke et al. (2023)](https://link.springer.com/article/10.1186/s12874-023-02012-5),
"best-practice recommendations are not available" for carryover modeling in
aggregated N-of-1 trials, and researchers have introduced methods like the
"carry-over adjusted parametric model (COAPM)."

**Project Status:** The pmsimstats2025 approach using exponential decay with
configurable half-life is well-aligned with current methodology. The variance-power
tradeoff whitepaper provides rigorous justification that is ahead of published
literature.

**Recommendation:** Add citation to Wienke et al. (2023) and position the
project as addressing the identified gap in best-practice recommendations.

### 3.2 Biomarker-Treatment Interaction Detection

**Web Consensus:** [Tian et al. (2018)](https://trialsjournal.biomedcentral.com/articles/10.1186/s13063-018-2491-0)
demonstrate that "tests of interactions are often lacking statistical power" and
that "the probability of detecting a biomarker-treatment interaction can be
increased by including prognostic variables."

**Project Status:** The simulation includes biomarker-response correlation in the
covariance structure, which is consistent with best practice. The power analysis
across biomarker moderation levels (0 to 0.65) provides valuable guidance.

**Recommendation:** Add discussion comparing power results to literature
benchmarks for interaction detection.

### 3.3 Type I Error Evaluation

**Web Consensus:** [Nolan et al. (2019)](https://www.mdpi.com/2227-9032/7/4/137)
found that "N-of-1 designs resulted in a higher type-I error probability than
parallel RCT and crossover designs when moderate-to-strong carryover effects
were not considered."

**Project Status:** The simulation now includes null conditions (bm_mod = 0) for
Type I error evaluation. Recent results show error rates of 3-9%, consistent with
the nominal 5% level.

**Recommendation:** Highlight Type I error control as a key finding; this
addresses a documented concern in the literature.

### 3.4 Autocorrelation Structure

**Web Consensus:** [AHRQ Design Guide (2014)](https://effectivehealthcare.ahrq.gov/products/n-1-trials/research-2014-1)
describes autoregressive models where "εt = δεt-1 + ut" with "δ the correlation
between consecutive errors."

**Project Status:** The implementation uses continuous-time AR(1) via `corCAR1`
in the DGP and analysis model, which is the correct approach for irregularly
spaced measurements. This is superior to discrete-time AR(1).

**Recommendation:** Document the choice of continuous-time AR(1) and its
advantages for clustered measurement schedules.

---

## 4. Strengths of Current Documentation

### 4.1 Theoretical Foundations

The `carryover_variance_whitepaper.tex` provides rigorous mathematical
explanation of the variance-power tradeoff that is not found in published
literature. Key strengths:

- Clear derivation of SE dependence on predictor variance
- Numerical examples with realistic parameter values
- Practical implications for trial design
- Proper acknowledgment that power loss is not model misspecification

### 4.2 DGP-Analysis Model Matching

Section 2.3 of the whitepaper correctly explains why exponential decay in the DGP
with linear terms in the analysis model is not misspecification. This addresses
a common source of confusion.

### 4.3 Covariance Matrix Construction

The partitioned construction approach guarantees positive definiteness and is
well-documented in `docs/sigma_matrix_derivation.tex`.

---

## 5. Required Updates for Publication Readiness

### 5.1 High Priority (Must Fix)

1. **Update model references from `lmer` to `lme`**
   - All documentation files
   - Explain why: `corCAR1` support for proper AR(1)

2. **Add Type I error results to main findings**
   - Update README.md with Type I error rates
   - Add "Null" row to all summary tables

3. **Synchronize parameter values**
   - Document actual correlation values used
   - Explain reduction from Hendrickson's values

### 5.2 Medium Priority (Should Fix)

4. **Update carryover parameterization**
   - Switch all documentation to half-life notation
   - Provide conversion formulas

5. **Add contemporary references**
   - Wienke et al. (2023) for carryover methods
   - Chen & Ho (2023) for Bayesian approaches
   - JMIR (2019) for N-of-1 simulation methodology

6. **Clarify script-design mapping**
   - Which designs in which scripts
   - Why designs are grouped by measurement schedule

### 5.3 Low Priority (Nice to Have)

7. **Add power comparison to literature benchmarks**

8. **Expand discussion of washout period alternatives**

9. **Add sensitivity analysis guidance**

---

## 6. Specific File Updates Required

### 6.1 README.md Updates

```markdown
## Line 27: Change
- **Mixed-Effects Analysis**: `lmer` models with and without explicit carryover
+ **Mixed-Effects Analysis**: `lme` models (nlme) with corCAR1 for AR(1)
  correlation and with/without explicit carryover terms

## Lines 154-162: Replace lmer examples
```r
# Change from:
lmer(response ~ treatment * biomarker + week + carryover_effect +
       (1 | participant_id))

# To:
lme(response ~ effective_drug_weeks * bm_centered + week,
    random = ~ 1 | participant_id,
    correlation = corCAR1(form = ~ week | participant_id),
    data = trial_data)
```

## Add new section after line 180:
### Type I Error Control

Under null conditions (biomarker_moderation = 0), empirical Type I error rates
range from 3-9% across designs and carryover levels, consistent with the
nominal 5% level (100 iterations, 95% CI ≈ 0.01-0.11).
```

### 6.2 simulation_white_paper.md Updates

Major revision needed:

1. Change all `lmer` references to `lme`
2. Update parameter table (Section 3.5) with actual values
3. Change `carryover_decay_rate` to `carryover_halflife`
4. Add section on AR(1) correlation structure choice
5. Update pseudocode (Section 7) to reflect current implementation

### 6.3 quickstart.Rmd Updates

```markdown
## Lines 131-147: Replace model examples

### With Explicit Carryover (Current Implementation)

```r
library(nlme)

fit_with <- lme(
  response ~ effective_drug_weeks * bm_centered + week,
  random = ~ 1 | participant_id,
  correlation = corCAR1(form = ~ week | participant_id),
  data = trial_data
)
```

### Without Carryover Term

```r
fit_without <- lme(
  response ~ treatment * bm_centered + week,
  random = ~ 1 | participant_id,
  correlation = corCAR1(form = ~ week | participant_id),
  data = trial_data
)
```
```

---

## 7. New Documentation Recommendations

### 7.1 Create: METHODS_COMPARISON.md

Document how pmsimstats2025 compares to:

- Hendrickson et al. (2020) original implementation
- Wienke et al. (2023) COAPM approach
- Chen & Ho (2023) Bayesian distributed lag model
- AHRQ Design Guide (2014) recommendations

### 7.2 Create: PARAMETER_SELECTION_RATIONALE.md

Explain:

- Why correlation values were reduced from Hendrickson
- Positive definiteness constraints
- Half-life vs proportion parameterization choice
- Measurement schedule design rationale

### 7.3 Update: CITATION.cff

Add preprint/publication information when available.

---

## 8. Conclusion

The pmsimstats2025 project has strong methodological foundations and addresses
an identified gap in the literature regarding best practices for carryover
modeling in N-of-1 trials. However, documentation has drifted from the current
implementation, particularly regarding:

1. Model engine (`lme` vs `lmer`)
2. Correlation parameter values
3. Carryover parameterization

Addressing the high-priority updates identified in Section 5.1 is essential
before publication. The variance-power tradeoff whitepaper represents a novel
contribution to the field and should be highlighted as a key deliverable.

---

## References

- [Chen & Ho (2023). Analysis of N-of-1 trials using Bayesian distributed lag model. Statistics in Medicine](https://link.springer.com/article/10.1186/s12874-023-02012-5)
- [Nolan et al. (2019). Comparison of Aggregated N-of-1 Trials with RCTs. Healthcare](https://www.mdpi.com/2227-9032/7/4/137)
- [Percha & Altman (2019). Designing Robust N-of-1 Studies. JMIR](https://www.jmir.org/2019/4/e12641/)
- [Schmid & Yang (2014). Statistical Design for N-of-1 Trials. AHRQ](https://effectivehealthcare.ahrq.gov/products/n-1-trials/research-2014-1)
- [Tian et al. (2018). Biomarker-treatment interaction in RCTs. Trials](https://trialsjournal.biomedcentral.com/articles/10.1186/s13063-018-2491-0)
- [Wienke et al. (2023). Comparison of methods for N-of-1 carry-over. BMC Med Res Methodol](https://link.springer.com/article/10.1186/s12874-023-02012-5)

---

*Report generated: February 2026*
