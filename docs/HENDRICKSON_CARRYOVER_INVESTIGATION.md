# Investigation: Why Power Doesn't Decrease with Carryover

## Problem Statement

In Hendrickson et al. (2020) Figure 4, power **decreases** as carryover half-life increases (when carryover is not modeled). Our simulations show **no decrease** - power stays roughly constant across carryover levels.

## Hendrickson's Expected Pattern

**Figure 4 (WITHOUT carryover in model)**:
- No carryover (t½=0): Power ~80%
- t½=1 week: Power ~70%
- t½=2 weeks: Power ~60%
- t½=4 weeks: Power ~50%

**Clear downward trend**: More carryover → Less power

## Our Current Results

**WITHOUT carryover model (c.bm = 0.48)**:
- t½=0: Power = 15%
- t½=1: Power = 20%
- t½=2: Power = 15%

**No clear pattern**: Power varies randomly, no systematic decrease

## Theoretical Expectation

When carryover exists but is NOT modeled:

1. **Systematic bias**: Control periods (off treatment) have elevated responses due to carryover from previous treatment periods

2. **Increased residual variance**: Carryover varies by:
   - Individual response magnitude (related to biomarker)
   - Time since treatment stopped
   - Creates heteroscedasticity

3. **Confounding with interaction**:
   - High biomarker → Strong treatment response → Strong carryover
   - Carryover noise is correlated with the biomarker×treatment interaction
   - Makes interaction harder to detect

4. **Power reduction mechanism**:
   - Unmodeled systematic effect increases residual variance
   - Reduces signal-to-noise ratio
   - Decreases power to detect interaction

## Hypothesis: What We Might Be Missing

### Hypothesis 1: Carryover Not Actually in Generated Data

**Check**: Does `generate_data()` actually add carryover to the response values?

**Location**: `pm_functions.R`, `apply_carryover_to_component()` function

**Test**: Generate data with/without carryover using same seed, compare values

### Hypothesis 2: Time Effect Absorbing Carryover

**Issue**: We include `week` as a covariate in the model

**Possibility**: In alternating design (A-B-A-B), carryover creates a systematic pattern that's collinear with time

**Effect**: Time covariate "mops up" the carryover effect even when not explicitly modeled

**Test**: Run model without `week` covariate, see if power pattern changes

### Hypothesis 3: Carryover Scale Too Small

**Current**: Carryover added to bio_response means only

**Issue**: If treatment effect is small relative to noise, carryover might be negligible

**Check**: Calculate actual magnitude of carryover vs treatment effect

### Hypothesis 4: Wrong Carryover Model

**Hendrickson uses**: Simple exponential decay

**We use**: `apply_carryover_to_component()` with specific logic

**Difference**: Need to verify we're implementing exact same decay

### Hypothesis 5: Biomarker-Carryover Interaction Missing

**Critical**: Carryover should be **person-specific**
- High biomarker → Strong response → Strong carryover
- Low biomarker → Weak response → Weak carryover

**Check**: Is carryover currently person-specific or uniform?

**If uniform**: This is the problem! Uniform carryover can be absorbed by time effect

**If person-specific**: Should create the noise pattern that reduces power

## Investigation Steps

### Step 1: Verify Carryover in Data ✓ (In Progress)

```r
# Generate same data with/without carryover
# Same seed → Same MVN draws
# Difference = carryover effect
```

**Expected**: Non-zero differences at OFF timepoints following ON periods

### Step 2: Check Carryover Magnitude

```r
# Compare:
# - Treatment effect size
# - Carryover effect size
# - Residual SD
```

**Expected**: Carryover should be non-trivial relative to noise

### Step 3: Test Time Effect Hypothesis

```r
# Run model without week covariate
lmer(response ~ treatment * bm + (1 | participant_id))
```

**Prediction**: If time absorbs carryover, removing it should show power decrease pattern

### Step 4: Check Person-Specificity

```r
# Is carryover:
# A) Same for all participants? (WRONG - can be absorbed by time)
# B) Proportional to individual's response? (CORRECT - creates noise)
```

### Step 5: Review Hendrickson Implementation

Need to verify:
1. Exactly how Hendrickson adds carryover to data
2. Whether it's person-specific or uniform
3. Whether they include time effect in analysis model

## Key Code Sections to Review

### 1. Carryover Application
**File**: `pm_functions.R`
**Lines**: 221-272 (`apply_carryover_to_component`)

```r
component_means[idx] <- component_means[idx] +
  component_means[prev_idx] * decay_factor
```

**Question**: Is this adding previous INDIVIDUAL's response or population mean?

### 2. Data Generation
**File**: `pm_functions.R`
**Function**: `generate_data()`

**Question**: When MVN data is drawn, does each person get their own response that carries over, or is carryover applied uniformly?

### 3. Analysis Model
**File**: `full_pmsim_analysis_hyb_versus_co.R`
**Lines**: 217-224

```r
# WITH carryover model
response ~ treatment * bm + week + carryover_effect + (1 | participant_id)

# WITHOUT carryover model
response ~ treatment * bm + week + (1 | participant_id)
```

**Question**: Does including `week` prevent us from seeing carryover effect on power?

## Critical Insight

Looking at the carryover application code:

```r
component_means[idx] <- component_means[idx] +
  component_means[prev_idx] * decay_factor
```

This adds **previous period's MEAN** to current period. But these are POPULATION means, not individual-specific values!

**If this is the problem**: Carryover is being added to means BEFORE MVN sampling, making it a systematic shift that affects all participants equally.

**What we need**: Carryover should be added AFTER MVN sampling, so each person's carryover depends on THEIR response magnitude.

## Proposed Fix

**Current flow**:
1. Calculate means (including carryover)
2. Draw MVN data using those means
3. Result: Everyone gets same carryover (systematic)

**Correct flow**:
1. Calculate means (without carryover)
2. Draw MVN data
3. Add person-specific carryover based on individual's previous response
4. Result: Carryover varies by person (creates noise)

## Next Actions

1. ✓ Increase iterations to 100
2. Verify carryover is in generated data
3. Check if carryover is person-specific or population-level
4. If population-level, modify to be person-specific
5. Re-run simulation and check for power decrease pattern

---

*Investigation started: 2025-11-18*
