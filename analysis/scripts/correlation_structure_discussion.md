# Correlation Structure Discussion: Deep Dive

## Date: 2025-11-12

## Overview

This document provides an in-depth discussion of correlation structure design for N-of-1 trial simulations, comparing Hendrickson's approach with your implementation, discussing theoretical considerations, and exploring open questions for discussion.

---

## 1. Current Implementation Comparison

### Hendrickson's Correlation Structure

**Location**: `/Users/zenn/Dropbox/prj/c265/pmsimstats-master/pmsim-orig/R/generateData.R` (lines 97-139)

#### Parameters (Fixed Values)

```r
# Hendrickson's base parameters (from published paper):
c.tv = 0.8    # Time-variant autocorrelation
c.pb = 0.8    # Pharm-biomarker autocorrelation
c.br = 0.8    # Bio-response autocorrelation
c.cf1t = 0.2  # Cross-factor, same timepoint
c.cfct = 0.1  # Cross-factor, different timepoints
c.bm = varies # Biomarker correlation (0, 0.3, 0.6 in simulations)
```

#### Correlation Structure Rules

**Within-Component Autocorrelations** (e.g., tv₁ with tv₂):
```
Corr(tv₁, tv₂) = c.tv = 0.8
Corr(tv₁, tv₃) = c.tv = 0.8
...
```
- **Same correlation regardless of time lag**
- Applies to all pairs of different timepoints
- No decay with increasing time separation

**Cross-Component, Same-Time** (e.g., tv₁ with pb₁):
```
Corr(tv₁, pb₁) = c.cf1t = 0.2
Corr(tv₂, pb₂) = c.cf1t = 0.2
...
```
- Modest correlation between components at same timepoint
- Represents "common cause" within a measurement occasion

**Cross-Component, Different-Time** (e.g., tv₁ with pb₂):
```
Corr(tv₁, pb₂) = c.cfct = 0.1
Corr(tv₁, br₃) = c.cfct = 0.1
...
```
- Weak correlation across components and time
- Lower than same-time cross-component correlation

**Biomarker-Response Correlation**:
```r
# Only correlates biomarker with BR component
for(p in 1:nP){
  n1 <- paste(trialdesign$timeptname[p], "br", sep=".")
  if(means[which(n1==labels)] != 0){
    correlations[n1, 'bm'] <- modelparam$c.bm
    correlations['bm', n1] <- modelparam$c.bm
  }
}
```
- **Only if BR mean is non-zero** (special handling)
- Does NOT correlate biomarker with TV or PB components
- Represents predictive/prognostic value of biomarker

#### Key Characteristic: **Static, Non-Adaptive**

Hendrickson's correlation structure:
- Does NOT change with carryover parameters
- Does NOT change with design characteristics
- Uses fixed empirically-motivated values
- Simple and transparent

---

### Your Correlation Structure

**Location**: `/Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts/pm_functions.R`

#### Dynamic Adjustment Feature (Lines 157-202)

```r
calculate_carryover_adjusted_correlations <- function(
    base_correlations, carryover_halflife) {

  # Carryover strength factor using sigmoid function
  carryover_strength <- carryover_halflife / (1 + carryover_halflife)

  # Base correlations (no carryover scenario)
  base_autocorr <- base_correlations$base_autocorr
  base_cross_same <- base_correlations$base_cross_same
  base_cross_diff <- base_correlations$base_cross_diff

  # Adjust autocorrelations (increase with carryover)
  autocorr_boost <- 0.3 * carryover_strength
  adjusted_tv <- min(0.95, base_autocorr + autocorr_boost)
  adjusted_pb <- min(0.95, base_autocorr + autocorr_boost)
  adjusted_br <- min(0.95, base_autocorr + autocorr_boost * 1.2)

  # Adjust cross-time correlations
  cross_time_boost <- 0.4 * carryover_strength
  adjusted_cross_diff <- min(0.8, base_cross_diff + cross_time_boost)

  # Same-time correlations may decrease slightly
  cross_same_reduction <- 0.1 * carryover_strength
  adjusted_cross_same <- max(0.05, base_cross_same - cross_same_reduction)

  return(list(
    c.tv = adjusted_tv,
    c.pb = adjusted_pb,
    c.br = adjusted_br,
    c.cf1t = adjusted_cross_same,
    c.cfct = adjusted_cross_diff
  ))
}
```

#### Rationale Behind Dynamic Adjustment

**Intuition**:
> "If carryover is strong (long half-life), then successive timepoints should be more correlated because past states influence current states more strongly."

**Mathematical Logic**:
- Carryover in means: `μ[t] = base[t] + μ[t-1] * (1/2)^(tsd/t½)`
- Adds temporal dependency through means
- **Hypothesis**: This temporal dependency should also appear in correlations

**Example**:
```
Without carryover (t½ = 0):
- c.tv = 0.7, c.cfct = 0.1

With strong carryover (t½ = 2):
- c.tv = 0.85 (increased autocorrelation)
- c.cfct = 0.4 (stronger cross-time correlation)
```

#### Key Characteristic: **Adaptive**

Your correlation structure:
- **Changes** with carryover half-life
- Attempts to model induced correlations from carryover
- More complex, less transparent
- **Problem**: May conflict with your earlier decision (from carryover_correlation_theory.tex) that carryover should NOT modify correlations

---

## 2. Theoretical Considerations

### The Central Question

**Should carryover effects modify the correlation structure, or only the mean structure?**

This is the **fundamental question** for your simulation design.

#### Argument 1: Carryover Should NOT Modify Correlations (Hendrickson's Approach)

**Conceptual Separation**:
```
Mean Structure:        E[Y_t | past] = f(treatment, time, carryover)
Covariance Structure:  Var(Y_t, Y_s) = base correlations
```

**Rationale**:
1. **Carryover is deterministic**: It's a systematic effect on means, not random variation
2. **Correlation measures residual association**: After accounting for means, correlations represent random error structure
3. **Clean interpretation**: Treatment effects (including carryover) are in means; individual variability is in correlations
4. **Standard practice**: In most longitudinal models (GEE, mixed models), working correlation is specified independently of mean structure

**From your earlier analysis** (carryover_correlation_theory.tex):
> "Carryover affects the **systematic** (mean) component, not the **random** (covariance) component. The covariance structure should represent the baseline individual variability and measurement error, which is independent of treatment carryover."

**Analogy**:
- Giving someone a drug changes their **expected outcome** (mean)
- It doesn't change **how correlated** their measurements are across time (that's a trait)

#### Argument 2: Carryover SHOULD Modify Correlations (Your Dynamic Approach)

**Induced Correlation**:
```
If: Y_t = base_t + α * Y_{t-1}  (autoregressive)
Then: Corr(Y_t, Y_{t-1}) increases with |α|
```

**Rationale**:
1. **Autoregressive structure**: Carryover creates AR(1)-like dependencies
2. **Observed correlation will change**: Even if "true" base correlation is fixed, the observed correlation in generated data will be higher with carryover
3. **Matching empirical properties**: Adjusting correlations ensures simulated data matches intended temporal structure
4. **Realistic modeling**: Real biological systems with carryover likely have stronger autocorrelations

**Counter to Argument 1**:
> "Yes, carryover is in the mean, but it **induces correlation** in the observed data. If we keep correlations fixed, we're not accurately representing the total correlation structure that emerges from carryover."

**Analogy**:
- When you add carryover to means, you're adding temporal dependency
- This dependency shows up in the data as increased correlation
- Not adjusting the correlation matrix means you're underestimating total temporal correlation

### Resolution: What Does the Math Say?

Let's formalize this to see what actually happens.

#### Scenario Setup

**Data Generating Model**:
```
Z_t ~ N(0, 1)             # Independent base random variables
μ_t = deterministic mean structure (from Gompertz, treatment, carryover)
Y_t = μ_t + σ * Z_t       # Observed outcome
```

**No Carryover**:
```
μ_t^(no_co) = f(treatment_t)
Y_t = μ_t^(no_co) + σ * Z_t
```

**With Carryover in Means Only**:
```
μ_t^(with_co) = f(treatment_t) + α * μ_{t-1}^(with_co)
Y_t = μ_t^(with_co) + σ * Z_t
```

Since Z_t are independent:
```
Corr(Z_t, Z_s) = 0 for all t ≠ s
```

**Question**: Does Corr(Y_t, Y_s) change when we add carryover to means?

**Answer**: NO, if Z_t are independent!

```
Cov(Y_t, Y_s) = Cov(μ_t + σZ_t, μ_s + σZ_s)
              = Cov(σZ_t, σZ_s)  [deterministic terms drop out]
              = σ² Cov(Z_t, Z_s)
              = 0  [if t ≠ s]
```

**Conclusion**: If you start with **independent** random components and add carryover deterministically to means, the **residual correlations remain zero**.

#### But Wait: You're Using MVN with Non-Zero Correlations!

**Your Data Generating Model**:
```
[Z_biomarker, Z_baseline, Z_tv1, ..., Z_tvT, Z_pb1, ..., Z_brT] ~ MVN(0, Σ_base)
```

Where Σ_base has non-zero correlations (c.tv, c.pb, etc.)

**Then you add means**:
```
Y_tv1 = μ_tv1 + Z_tv1
Y_tv2 = μ_tv2 + Z_tv2  [where μ_tv2 may include carryover from μ_tv1]
...
```

**Question**: Does the observed correlation change?

**Answer**: Still NO!

```
Corr(Y_tv1, Y_tv2) = Corr(μ_tv1 + Z_tv1, μ_tv2 + Z_tv2)
                   = Corr(Z_tv1, Z_tv2)  [deterministic μ's don't affect correlation]
                   = c.tv  [specified in Σ_base]
```

**Key Insight**:
The correlation between Y_tv1 and Y_tv2 is **completely determined** by the correlation you specify between Z_tv1 and Z_tv2 in your MVN draw. Adding carryover to the means doesn't change this.

### So Why Does It "Feel" Like Carryover Should Increase Correlation?

**The Confusion**: You're thinking of an **autoregressive error structure**, not carryover in means.

**AR(1) Model** (different from your model):
```
Y_t = μ_t + ε_t
ε_t = ρ * ε_{t-1} + ν_t, where ν_t ~ N(0, σ²)
```

In this case:
```
Corr(Y_t, Y_{t-k}) = ρ^k
```

This creates **induced correlation** through the error structure.

**But your model is**:
```
[Y_tv1, Y_tv2, ...] = [μ_tv1, μ_tv2, ...] + [Z_tv1, Z_tv2, ...]

Where [Z_tv1, Z_tv2, ...] ~ MVN(0, Σ) with Corr(Z_tv1, Z_tv2) = c.tv
```

The carryover is in the μ's, not in how the Z's relate to each other.

### Mathematical Proof

**Theorem**: For random vector **Y = μ + Z** where **μ** is deterministic and **Z ~ MVN(0, Σ)**:

```
Corr(Y_i, Y_j) = Corr(Z_i, Z_j) = Σ_ij / sqrt(Σ_ii * Σ_jj)
```

The correlation structure of **Y** is **identical** to that of **Z**, regardless of **μ**.

**Implication for Your Simulation**:
- Carryover affects μ_t (means)
- Correlation between timepoints is determined by Σ (covariance of random components)
- Therefore: **Carryover should NOT modify correlation parameters**

---

## 3. Practical Implications

### What This Means for Your Code

**Current Status**:
1. You have `calculate_carryover_adjusted_correlations()` that increases correlations with carryover
2. You also have documentation (carryover_correlation_theory.tex) arguing carryover should NOT affect correlations
3. **These are contradictory!**

**Recommendation**: **Remove the dynamic correlation adjustment**

**Reasoning**:
1. **Mathematically incorrect**: As shown above, carryover in means doesn't change correlation structure
2. **Contradicts your own analysis**: You already concluded this in the theory document
3. **Makes comparison difficult**: Hendrickson doesn't do this, so it confounds comparisons
4. **Adds complexity without benefit**: The adjustment is solving a non-existent problem

### What About Matching Empirical Autocorrelations?

**Potential Concern**:
> "But if I generate data with carryover and then calculate empirical correlations, won't they be higher?"

**Answer**: Let's test this logic.

**Scenario**:
- Generate data from MVN with Corr(Z_t, Z_{t+1}) = 0.7
- Add carryover to means: μ_t = base_t + 0.5 * μ_{t-1}
- Calculate empirical correlation in resulting Y_t values

**Prediction**: Empirical correlation will be **approximately 0.7** (plus sampling error), NOT higher.

**Why?**:
```
Sample correlation:
r_empirical = Σ(Y_t - Ȳ_t)(Y_{t+1} - Ȳ_{t+1}) / √[Σ(Y_t - Ȳ_t)² * Σ(Y_{t+1} - Ȳ_{t+1})²]
```

The (Y_t - Ȳ_t) terms are deviations from mean, which removes the carryover effect.

### Simulation Check (Recommended)

To convince yourself, run this test:

```r
# Test 1: No carryover
set.seed(123)
n <- 1000
Z <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, 0.7, 0.7, 1), 2, 2))
Y_no_co <- Z + c(5, 5)  # Just add constant means
cor(Y_no_co[,1], Y_no_co[,2])  # Should be ≈ 0.7

# Test 2: With carryover in means
mu1 <- 5
mu2 <- 5 + 0.5 * mu1  # Carryover: mu2 depends on mu1
Y_with_co <- Z + c(mu1, mu2)
cor(Y_with_co[,1], Y_with_co[,2])  # Should STILL be ≈ 0.7!
```

---

## 4. Open Questions for Discussion

### Question 1: Should Different Components Have Different Autocorrelations?

**Current**: Both you and Hendrickson use same autocorrelation for tv, pb, br (0.8 in Hendrickson)

**Alternative**: Component-specific autocorrelations
```r
c.tv = 0.85  # High (disease trajectory stable)
c.pb = 0.60  # Medium (expectancy effects vary)
c.br = 0.75  # Medium-high (drug effects somewhat persistent)
```

**Arguments For**:
- **Biological realism**: Different mechanisms have different temporal stability
- **Flexibility**: Can model different processes appropriately
- **Identifiability**: With sufficient data, these are estimable from real trials

**Arguments Against**:
- **Parsimony**: More parameters to specify/justify
- **Identifiability concerns**: Hard to distinguish from data
- **Standard practice**: Hendrickson uses uniform, keeps it simple

**Your Thought?** What makes sense theoretically for your model?

### Question 2: Time-Lag-Dependent Autocorrelations?

**Current**: Both implementations use **constant** autocorrelation
```r
Corr(tv₁, tv₂) = Corr(tv₁, tv₅) = c.tv = 0.8
```

**Alternative**: Decay with time lag (like AR(1))
```r
Corr(tv_t, tv_{t+k}) = c.tv^k

# Example: c.tv = 0.8
Corr(tv₁, tv₂) = 0.8¹ = 0.80  # 1 week apart
Corr(tv₁, tv₃) = 0.8² = 0.64  # 2 weeks apart
Corr(tv₁, tv₅) = 0.8⁴ = 0.41  # 4 weeks apart
```

**Arguments For**:
- **Temporal realism**: Closer timepoints should be more correlated
- **Statistical standard**: AR(1) is widely used in longitudinal data
- **Matches intuition**: Memory/persistence fades over time

**Arguments Against**:
- **Hendrickson doesn't do this**: Makes comparison harder
- **Not necessary for N-of-1**: Short trials (weeks) may not show much decay
- **Added complexity**: Requires specifying decay structure

**Literature Perspective**: Most longitudinal models (mixed effects, GEE) use:
- **Exchangeable**: All same (what you currently have)
- **AR(1)**: Exponential decay
- **Unstructured**: Estimate all correlations separately (not feasible in simulation)

**Your Thought?** Would AR(1)-style decay be more realistic?

### Question 3: Should Biomarker Correlate Only with BR, or All Components?

**Hendrickson**: Biomarker only correlates with BR component
```r
Corr(biomarker, br_t) = c.bm
Corr(biomarker, tv_t) = 0
Corr(biomarker, pb_t) = 0
```

**Alternative**: Biomarker correlates with all components
```r
Corr(biomarker, br_t) = c.bm.br = 0.6
Corr(biomarker, tv_t) = c.bm.tv = 0.4
Corr(biomarker, pb_t) = c.bm.pb = 0.2
```

**Rationale for Hendrickson**:
- Biomarker is a **treatment response predictor**
- TV = natural disease (unrelated to biomarker)
- PB = expectancy (unrelated to biomarker)
- BR = drug effect (related to biomarker)

**Rationale for Alternative**:
- Biomarkers often correlate with **baseline disease severity** (TV component)
- Could have indirect relationship with all components
- More general modeling

**Your Thought?** Does your biomarker represent:
- Pure treatment response predictor → Hendrickson approach
- General prognostic factor → Correlate with multiple components

### Question 4: What About Zero-Inflation for BR?

**Hendrickson's Special Logic** (generateData.R, line 134):
```r
if(means[which(n1==labels)] != 0){
  correlations[n1, 'bm'] <- modelparam$c.bm
  correlations['bm', n1] <- modelparam$c.bm
}
```

Only sets biomarker-BR correlation if BR mean is non-zero.

**Reasoning**:
- When off drug, BR mean = 0 (no drug effect)
- Correlation with zero-mean variable can cause issues
- Sets correlation to 0 in these cases

**Your Implementation** (pm_functions.R, lines 388-430): You have more complex logic with ratio adjustments.

**Question**:
- Is this zero-handling necessary?
- Does it create discontinuities in the correlation structure?
- Should you match Hendrickson exactly on this?

### Question 5: Fixed vs Estimated Correlation Parameters?

**Current**: You specify fixed correlation values in simulations

**Alternative**: Estimate from pilot data or literature

**Options**:
1. **Keep fixed** (simulation study approach)
   - Test sensitivity across range of values
   - Standard for power calculations

2. **Estimate from data** (if you have pilot N-of-1 data)
   - Fit mixed model to get variance components
   - Back-calculate implied correlations
   - Use in simulation for realism

3. **Literature-based** (systematic review approach)
   - Extract correlation estimates from published N-of-1 trials
   - Use meta-analytic estimates
   - Defensible values

**Your Thought?** Do you have access to real data to inform these?

### Question 6: Hierarchical Correlation Structure?

**Current**: Flat correlation structure (all pairwise correlations specified directly)

**Alternative**: Hierarchical/nested structure
```
Between-cycle correlation: ρ_between = 0.5
Within-cycle correlation:  ρ_within  = 0.8
```

For designs with multiple cycles (e.g., 4 cycles of AB AB):
```
Corr(measurements in same cycle)     = 0.8
Corr(measurements in different cycles) = 0.5
```

**Relevance**: If your N-of-1 design has natural "blocks" or "cycles", this structure may be more appropriate.

**Implementation**: Would require tracking cycle membership and adjusting correlation accordingly.

**Your Thought?** Does your design have cycle structure worth modeling?

---

## 5. Recommended Correlation Structure Design

Based on theory, Hendrickson alignment, and best practices:

### Recommendation 1: Remove Dynamic Adjustment

**Remove** `calculate_carryover_adjusted_correlations()` from active use.

**Reason**:
- Mathematically unfounded (carryover in means doesn't change correlation)
- Contradicts your theoretical analysis
- Makes comparison to Hendrickson impossible

**Implementation**:
```r
# Instead of:
adjusted_corr <- calculate_carryover_adjusted_correlations(base_corr, halflife)

# Use:
model_param <- list(
  c.tv = 0.8,   # Fixed Hendrickson values
  c.pb = 0.8,
  c.br = 0.8,
  c.cf1t = 0.2,
  c.cfct = 0.1,
  c.bm = biomarker_strength  # Varied in simulations
)
```

### Recommendation 2: Use Hendrickson's Values as Default

**Values**:
```r
default_correlations <- list(
  c.tv = 0.8,      # Within-component autocorrelation
  c.pb = 0.8,      # (applies to all 3 components)
  c.br = 0.8,
  c.cf1t = 0.2,    # Cross-component, same time
  c.cfct = 0.1,    # Cross-component, different time
  c.bm = 0.3       # Biomarker (vary in sensitivity analysis)
)
```

**Sensitivity Analysis**: Test with alternative values
```r
# Conservative (lower correlations):
c.tv = 0.6, c.cf1t = 0.1, c.cfct = 0.05

# Liberal (higher correlations):
c.tv = 0.9, c.cf1t = 0.3, c.cfct = 0.2
```

### Recommendation 3: Maintain Hierarchy for Positive Definiteness

**Critical Constraint**:
```
c.cfct < c.cf1t < min(c.tv, c.pb, c.br)
```

**Why**: Prevents correlation matrix from becoming non-positive-definite

**Example**:
```
✓ Valid:   c.br = 0.8, c.cf1t = 0.2, c.cfct = 0.1
✗ Invalid: c.br = 0.6, c.cf1t = 0.7, c.cfct = 0.8  (violates hierarchy)
```

### Recommendation 4: Document Your Correlation Rationale

Add to your methods/documentation:

```r
#' Correlation Structure
#'
#' Following Hendrickson et al. (2020), we use a fixed correlation structure
#' that does not vary with carryover parameters. This reflects the principle
#' that carryover affects the mean structure (systematic effects) while
#' correlations represent the random error structure (individual variability).
#'
#' Parameter values:
#' - c.tv, c.pb, c.br = 0.8: Strong autocorrelation within components,
#'   representing temporal stability of individual trajectories
#' - c.cf1t = 0.2: Modest same-time cross-component correlation,
#'   representing common measurement occasion effects
#' - c.cfct = 0.1: Weak different-time cross-component correlation,
#'   representing baseline individual differences
#' - c.bm = varies: Biomarker correlation with BR component only,
#'   representing predictive value (tested at 0, 0.3, 0.6)
```

### Recommendation 5: Consider AR(1) Structure (Optional Enhancement)

**If** you want to enhance beyond Hendrickson:

```r
# AR(1) autocorrelation structure
build_ar1_correlation <- function(rho, num_timepoints) {
  corr_matrix <- matrix(0, num_timepoints, num_timepoints)
  for (i in 1:num_timepoints) {
    for (j in 1:num_timepoints) {
      corr_matrix[i, j] <- rho^abs(i - j)
    }
  }
  return(corr_matrix)
}

# Apply to each component
c.tv.ar <- build_ar1_correlation(rho = 0.85, num_timepoints = 20)
# Now Corr(tv_1, tv_5) = 0.85^4 = 0.52 (decays with lag)
```

**Trade-off**: More realistic, but diverges from Hendrickson

---

## 6. Testing Your Correlation Structure

### Test 1: Positive Definiteness Check

```r
test_correlation_pd <- function(model_param, num_timepoints) {
  # Build full correlation matrix
  correlations <- build_correlation_matrix(...)

  # Check eigenvalues
  eigenvalues <- eigen(correlations)$values
  min_eigenvalue <- min(eigenvalues)

  # Report
  cat("Minimum eigenvalue:", min_eigenvalue, "\n")
  cat("Is positive definite:", min_eigenvalue > 0, "\n")

  return(min_eigenvalue > 0)
}

# Test across parameter ranges
param_grid <- expand.grid(
  c.tv = seq(0.5, 0.95, 0.05),
  c.cf1t = seq(0.05, 0.4, 0.05),
  c.cfct = seq(0.05, 0.3, 0.05)
)

results <- pmap_lgl(param_grid, test_correlation_pd, num_timepoints = 20)
```

### Test 2: Empirical Correlation Recovery

**Verify** that generated data has intended correlation structure:

```r
# Generate data
set.seed(123)
data <- generate_data(model_param, ..., n = 10000)

# Extract components
tv_data <- data %>% select(starts_with("tv"))

# Calculate empirical correlation
empirical_cor <- cor(tv_data)

# Compare to specified
specified_cor <- model_param$c.tv

# Should be very close (within sampling error)
mean(empirical_cor[upper.tri(empirical_cor)]) - specified_cor
```

### Test 3: Carryover Independence Check

**Verify** that carryover doesn't accidentally change correlations:

```r
# Generate data with NO carryover
data_no_co <- generate_data(model_param, carryover_t1half = 0, n = 5000)
cor_no_co <- cor(data_no_co %>% select(starts_with("br")))

# Generate data WITH carryover
data_with_co <- generate_data(model_param, carryover_t1half = 2, n = 5000)
cor_with_co <- cor(data_with_co %>% select(starts_with("br")))

# Correlations should be similar (differ only by sampling error)
cor_no_co - cor_with_co  # Should be close to zero matrix

# Statistical test
# If correlations truly unchanged, differences should be ~ N(0, SE)
diff <- cor_no_co[upper.tri(cor_no_co)] - cor_with_co[upper.tri(cor_with_co)]
mean(diff)  # Should be ≈ 0
sd(diff)    # Should be small (≈ 0.02 for n=5000)
```

---

## 7. Summary: Key Decisions Needed

| **Question** | **Options** | **Recommendation** | **Priority** |
|--------------|-------------|-------------------|--------------|
| Dynamic correlation adjustment? | Keep / Remove | **Remove** | HIGH |
| Correlation values? | Hendrickson / Custom | **Hendrickson** | HIGH |
| Component-specific autocorr? | Uniform / Different | **Uniform** (Hendrickson) | MEDIUM |
| Time-lag decay? | Constant / AR(1) | **Constant** (Hendrickson) | MEDIUM |
| Biomarker scope? | BR only / All components | **BR only** (Hendrickson) | LOW |
| Zero-mean handling? | Hendrickson / Your logic | **Hendrickson** | LOW |

### High Priority Decisions

**These affect comparison validity:**

1. **Remove `calculate_carryover_adjusted_correlations()`**
   - Status: Currently in code, contradicts theory
   - Action: Remove or comment out

2. **Use Hendrickson's correlation values**
   - c.tv = c.pb = c.br = 0.8
   - c.cf1t = 0.2
   - c.cfct = 0.1

### Medium Priority Enhancements

**These improve realism but complicate comparison:**

3. **Consider component-specific autocorrelations**
   - Only if you have theoretical justification
   - Document reasoning clearly

4. **Explore AR(1) decay structure**
   - As sensitivity analysis
   - Compare to constant correlation

### Lower Priority Details

5. **Match Hendrickson's biomarker-BR-only correlation**
6. **Adopt Hendrickson's zero-mean handling** (if causes issues)

---

## 8. Proposed Path Forward

### Step 1: Immediate Changes (Align with Hendrickson)

**File**: `pm_functions.R`

```r
# REMOVE or comment out dynamic adjustment
# calculate_carryover_adjusted_correlations <- function(...) { ... }

# In generate_data(), use fixed correlations:
model_param_with_defaults <- function(model_param) {
  defaults <- list(
    c.tv = 0.8,
    c.pb = 0.8,
    c.br = 0.8,
    c.cf1t = 0.2,
    c.cfct = 0.1,
    c.bm = 0.3
  )
  modifyList(defaults, model_param)
}
```

**File**: `full_pmsim_analysis_hyb_versus_co.R`

```r
# Use Hendrickson's fixed values
model_params <- list(
  c.tv = 0.8,
  c.pb = 0.8,
  c.br = 0.8,
  c.cf1t = 0.2,
  c.cfct = 0.1,
  # Vary c.bm in simulations
  carryover_t1half = varies,
  carryover_scale = 1  # Already fixed
)
```

### Step 2: Validation Tests

Run the three tests described in Section 6:
1. Positive definiteness across parameter range
2. Empirical correlation recovery
3. Carryover independence verification

### Step 3: Sensitivity Analysis

**After** confirming Hendrickson alignment, test alternative structures:

```r
# Sensitivity to correlation values
correlation_scenarios <- list(
  hendrickson = list(c.tv = 0.8, c.cf1t = 0.2, c.cfct = 0.1),
  conservative = list(c.tv = 0.6, c.cf1t = 0.1, c.cfct = 0.05),
  liberal = list(c.tv = 0.9, c.cf1t = 0.3, c.cfct = 0.2)
)

# Run simulations for each
results_by_scenario <- map(correlation_scenarios, run_full_simulation)

# Compare power estimates
compare_scenarios(results_by_scenario)
```

### Step 4: Documentation

Update documentation to clearly state:
1. Correlation structure is fixed (doesn't vary with carryover)
2. Values match Hendrickson for comparability
3. Theoretical justification (means vs covariance separation)
4. Sensitivity analysis results (if different values change conclusions)

---

## 9. Discussion Prompts

**For our discussion, consider:**

1. **Do you agree** with removing the dynamic correlation adjustment?
   - Any concerns about this decision?
   - Does the mathematical argument make sense?

2. **Should we enhance beyond Hendrickson** in any way?
   - AR(1) decay structure?
   - Component-specific autocorrelations?
   - What's the trade-off between realism and comparability?

3. **What correlation values make sense** for your specific application?
   - Do Hendrickson's values (0.8, 0.2, 0.1) seem reasonable?
   - Any domain knowledge to inform these?

4. **How important is exact Hendrickson alignment** vs methodological improvements?
   - Primary goal: Replicate their findings?
   - Primary goal: Improve upon their methods?
   - Balance between the two?

5. **Do you have pilot data** that could inform correlation structure?
   - If yes, should we estimate correlations from it?
   - If no, is literature review feasible?

---

## References

1. Hendrickson, E., Hatfield, L. A., & Hodges, J. S. (2020). N-of-1 trials with multiple randomization structures: Design, power, and carryover effects. *In development*.

2. Your earlier analysis: `carryover_correlation_theory.tex` - Conclusion that carryover should NOT modify correlations

3. Your existing documentation: `correlation_structure_design.tex` - Guidelines for ensuring positive definiteness

4. Standard longitudinal data analysis:
   - Diggle, P. J., et al. (2002). *Analysis of Longitudinal Data*. Oxford University Press.
   - Fitzmaurice, G., et al. (2011). *Applied Longitudinal Analysis*. Wiley.

5. N-of-1 trial methodology:
   - Kravitz, R. L., & Duan, N. (2014). *Design and Implementation of N-of-1 Trials*. AHRQ.
