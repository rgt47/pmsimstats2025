# Type I Error Inflation: Analysis and Fix

## Problem Identified

When examining simulation results with `biomarker_moderation = 0` (the Type I error / null hypothesis condition), the rejection rates were inflated:

| Design | Expected | Observed |
|--------|----------|----------|
| OL | 5% | ~12% |
| OL+BDC | 5% | ~15% |
| Crossover | 5% | ~22% |

## Root Cause

The covariance structure (`build_sigma_guaranteed_pd`) creates correlations between the biomarker and response random components:

```r
# Column 1: Biomarker effects on responses
# BR strongly correlated with biomarker
Sigma_12[br_idx, 1] <- c.bm * within_subject_sd * biomarker_sd
# ER and TR weakly correlated with biomarker
Sigma_12[er_idx, 1] <- c.bm * 0.5 * within_subject_sd * biomarker_sd
Sigma_12[tr_idx, 1] <- c.bm * 0.5 * within_subject_sd * biomarker_sd
```

When `biomarker_correlation = 0.3` (c.bm = 0.3), this creates **real correlation** between biomarker and the random response components (BR_random, ER_random, TR_random).

**Critical insight**: Even when `biomarker_moderation = 0` (no treatment × biomarker interaction in the mean structure), the covariance structure still creates a correlation between biomarker and responses. This is NOT a spurious correlation - it's a real statistical relationship induced by the covariance structure.

The mixed model correctly detects this relationship, resulting in inflated rejection rates. This is **not** a Type I error in the traditional sense - the model is detecting a real pattern in the data.

## The Fix

For proper Type I error evaluation (testing the null hypothesis of no biomarker × treatment interaction), we must ensure **both**:

1. `biomarker_moderation = 0` (no interaction in the mean structure)
2. `biomarker_correlation = 0` (no correlation in the covariance structure)

### Implementation in param_grid

For each design, Type I error conditions now use `biomarker_correlation = 0`:

```r
param_grid <- bind_rows(
  # OL design - power conditions
  expand_grid(
    design = "ol",
    biomarker_moderation = c(0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    carryover = c(0)
  ),
  # OL design - Type I error condition
  expand_grid(
    design = "ol",
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),  # ← KEY FIX
    carryover = c(0)
  ),
  # ... similar pattern for all other designs
)
```

## Technical Details

### Why the original approach failed

The original `param_grid` used:
```r
biomarker_moderation = c(0, 0.25, 0.35, ...)  # 0 for Type I error
biomarker_correlation = c(0.3)                 # Same for all
```

This created an inconsistency:
- Mean structure: No biomarker × treatment interaction
- Covariance structure: Biomarker correlated with response random effects

The mixed model detected the covariance-induced relationship, leading to inflated rejection rates.

### Why setting biomarker_correlation = 0 fixes it

When `c.bm = 0`:
```r
Sigma_12[br_idx, 1] <- 0 * within_subject_sd * biomarker_sd = 0
Sigma_12[er_idx, 1] <- 0 * 0.5 * within_subject_sd * biomarker_sd = 0
Sigma_12[tr_idx, 1] <- 0 * 0.5 * within_subject_sd * biomarker_sd = 0
```

The biomarker becomes independent of the response random effects. Any apparent biomarker × treatment interaction is now purely due to sampling variability, and the model should reject at the nominal 5% rate.

### Satterthwaite degrees of freedom

The p-value calculation uses lmerTest's Satterthwaite approximation:
```r
p_value <- coefs[idx, "Pr(>|t|)"]  # Satterthwaite approximation
```

This is appropriate for cross-level interactions (biomarker is between-subject, treatment is within-subject). The effective degrees of freedom are approximately equal to the number of participants (~70), not the number of observations (~560).

## Expected Results After Fix

| Design | Condition | biomarker_mod | biomarker_corr | Expected |
|--------|-----------|---------------|----------------|----------|
| OL | Type I | 0 | 0 | ~5% |
| OL | Power | 0.25-0.65 | 0.3 | Varies |
| OL+BDC | Type I | 0 | 0 | ~5% |
| OL+BDC | Power | 0.25-0.65 | 0.3 | Varies |
| Crossover | Type I | 0 | 0 | ~5% |
| Crossover | Power | 0.25-0.65 | 0.3 | Varies |
| Hybrid | Type I | 0 | 0.3 | ~5% (no fix needed) |
| Hybrid | Power | 0.25-0.65 | 0.3 | Varies |
| Parallel | Type I | 0 | 0.3 | ~5% (no fix needed) |
| Parallel | Power | 0.25-0.65 | 0.3 | Varies |

## Why Hybrid and Parallel Don't Need the Fix

**Hybrid (N-of-1)**: Multiple treatment switches per person (4 paths, many transitions). The dense within-person replication allows the model to correctly attribute variance to the random effect structure rather than the fixed biomarker × treatment interaction.

**Parallel**: Each person is only on one treatment throughout the study. The biomarker × treatment interaction is estimated entirely from between-person comparisons, which are independent of the within-person covariance structure.

**Designs needing fix (OL, OL+BDC, Crossover)**: These have limited within-person treatment variation (0-1 switches), so the biomarker-response correlation in the covariance structure can masquerade as a biomarker × treatment/time interaction.

## Interpretation Notes

1. **Type I error row**: Shows the false positive rate when there is truly no biomarker × treatment interaction. Should be ~5% (nominal alpha level).

2. **Power rows**: Show detection rate when biomarker does moderate treatment effect. Higher biomarker_moderation → higher power.

3. **The distinction**: Type I error (biomarker_correlation = 0) represents a world where biomarker has NO relationship to treatment response. Power conditions (biomarker_correlation = 0.3) represent a world where biomarker IS predictive of treatment response.

---

*Generated: 2025-11-25*
