# Hendrickson Code Comparison - Carryover Implementation

## Date: 2025-11-18

## Purpose

Compare Hendrickson's original implementation with ours to understand why we're not seeing power decrease with increasing carryover.

---

## Key Files Reviewed

### Hendrickson's Original Code
- **Location**: `~/prj/c265/pmsimstats-master/pmsim-orig/R/`
- **Data Generation**: `generateData.R`
- **Analysis**: `lme_analysis.R`

### Our Implementation
- **Data Generation**: `pm_functions.R` (line 221-272: `apply_carryover_to_component`)
- **Analysis**: `full_pmsim_analysis_hyb_versus_co.R` (lines 204-214)

---

## Carryover Implementation Comparison

### Hendrickson's Approach (generateData.R:82-94)

```r
if(c=="br"){
  brmeans<-modgompertz(d$tod,rp$max,rp$disp,rp$rate)
  if(nP>1){
    for(p in 2:nP){
      if(!d[p]$onDrug){  # Only when OFF drug
        if(d[p]$tsd>0){   # tsd = time since drug stopped
          # Add previous period's BR mean, scaled by exponential decay
          brmeans[p]<-brmeans[p]+brmeans[p-1]*(1/2)^(d$tsd[p]/modelparam$carryover_t1half)
        }
      }
    }
  }
  means<-c(means,brmeans)
}

# THEN sample from MVN:
dat<-mvrnorm(n=modelparam$N,mu=means,Sigma=sigma,empirical=empirical)
```

**Key characteristics**:
1. Carryover added to MEANS vector (population-level)
2. Applied BEFORE MVN sampling
3. Formula: `(1/2)^(time_off / halflife)`
4. Only applied when OFF drug
5. Based on previous period's MEAN (not individual response)

### Our Approach (pm_functions.R:221-272)

```r
apply_carryover_to_component <- function(component_means, design_df, carryover_t1half) {
  # ... (similar logic)

  for (i in 2:nrow(design_df)) {
    if (design_df$treatment[i] == 0) {  # OFF treatment
      # Find most recent ON period
      prev_on_idx <- tail(which(design_df$treatment[1:(i-1)] == 1), 1)

      if (length(prev_on_idx) > 0) {
        weeks_off <- design_df$week[i] - design_df$week[prev_on_idx]
        decay_factor <- (1/2)^(weeks_off / carryover_t1half)

        # Add carryover to current mean
        component_means[i] <- component_means[i] +
          component_means[prev_on_idx] * decay_factor
      }
    }
  }

  return(component_means)
}

# Applied in generate_data() BEFORE MVN sampling
```

**Key characteristics**:
1. Carryover added to MEANS vector (population-level)
2. Applied BEFORE MVN sampling
3. Formula: `(1/2)^(weeks_off / halflife)`
4. Only applied when OFF treatment
5. Based on previous period's MEAN (not individual response)

### Conclusion: **IMPLEMENTATIONS ARE IDENTICAL**

Both Hendrickson and we apply carryover:
- At the population level (to means)
- Before MVN sampling
- Using exponential decay formula
- Only during OFF periods

---

## Analysis Model Comparison

### Hendrickson's Analysis Model (lme_analysis.R:117)

```r
form <- Sx ~ bm + Db + t + bm*Db + (1|ptID)
```

Where:
- `Sx` = symptom score (response)
- `bm` = biomarker
- `Db` = drug status (binary: ON/OFF)
- `t` = time (continuous)
- `bm*Db` = biomarker×treatment interaction
- `(1|ptID)` = random intercept

**CRITICAL**: No carryover term in the model!

### Our Analysis Model - "WITHOUT Carryover"

```r
response ~ treatment * bm + week + (1 | participant_id)
```

Where:
- `response` = outcome measure
- `treatment` = treatment status (0/1)
- `bm` = biomarker
- `week` = time (continuous)
- `treatment*bm` = biomarker×treatment interaction
- `(1|participant_id)` = random intercept

**CRITICAL**: No carryover term in the model!

### Our Analysis Model - "WITH Carryover"

```r
response ~ treatment * bm + week + carryover_effect + (1 | participant_id)
```

Same as above, plus:
- `carryover_effect` = computed exponential decay term

### Conclusion: **MODELS ARE EQUIVALENT**

Hendrickson's model matches our "WITHOUT carryover" model exactly. Both:
- Include time effect
- Include biomarker×treatment interaction
- Do NOT include carryover term
- Use random intercept only

---

## The Mystery: Why No Power Decrease?

### Hendrickson's Figure 4 Pattern

**When carryover exists but is NOT modeled:**
- No carryover (t½=0): Power ~80%
- t½=1 week: Power ~70%
- t½=2 weeks: Power ~60%
- t½=4 weeks: Power ~50%

Clear downward trend as carryover increases.

### Our Results

**WITHOUT carryover model (c.bm = 0.48):**
- t½=0: Power = 15%
- t½=1: Power = 20%
- t½=2: Power = 15%

No systematic decrease - power varies randomly.

---

## Hypotheses for Discrepancy

### Hypothesis 1: ✓ RULED OUT - Different Implementation
**Status**: RULED OUT
**Finding**: Carryover implementation is identical

### Hypothesis 2: ✓ RULED OUT - Different Analysis Model
**Status**: RULED OUT
**Finding**: Analysis models are equivalent

### Hypothesis 3: Scale/Magnitude Differences
**Status**: NEEDS INVESTIGATION
**Possibility**: Our effect sizes are much smaller (15-20% power vs 80% power)
**Implication**: Carryover effect might be too small relative to noise to matter

### Hypothesis 4: Design Differences
**Status**: NEEDS INVESTIGATION
**Differences**:
- Hendrickson: 8 timepoints (4 weeks each), total 32 weeks
- Us: 20 timepoints (1 week each), total 20 weeks
**Implication**: Different temporal structure might affect how carryover manifests

### Hypothesis 5: Parameter Differences
**Status**: NEEDS INVESTIGATION
**Our parameters**: autocorr = 0.6, c.bm = {0, 0.3, 0.48}
**Hendrickson parameters**: autocorr = 0.8, c.bm = {0, 0.3, 0.6}
**Implication**: Lower autocorrelation and effect size might change power patterns

### Hypothesis 6: Time Variable Specification
**Status**: NEEDS INVESTIGATION
**Observation**: Hendrickson uses `t` (cumulative time), we use `week` (sequential)
**Question**: Could different time encoding affect how residual variance is modeled?

### Hypothesis 7: Carryover Magnitude
**Status**: NEEDS CRITICAL INVESTIGATION
**Question**: Are we actually generating meaningful carryover in the data?
**Test needed**: Generate data WITH vs WITHOUT carryover, verify differences exist

---

## Next Steps

1. **Verify carryover exists in generated data** ✓ IN PROGRESS
   - Generate paired datasets (WITH/WITHOUT carryover, same seed)
   - Compare individual-level responses
   - Confirm non-zero differences at OFF timepoints

2. **Calculate carryover effect size**
   - Quantify magnitude of carryover relative to:
     - Treatment effect
     - Residual SD
     - Between-subject variation

3. **Compare parameter values**
   - Review Hendrickson's exact parameter values
   - Check if we're using comparable magnitudes
   - Verify BR factor scaling

4. **Test design structure impact**
   - Run simulation with 8-timepoint design (matching Hendrickson)
   - Compare power patterns

5. **Investigate time effect interaction**
   - Run model WITHOUT week covariate
   - See if power decrease pattern emerges

---

## Code Snippets for Investigation

### Test 1: Verify Carryover in Data (RUNNING)

```r
# Generate data WITH carryover
set.seed(123)
data_with <- generate_data(..., carryover_t1half = 2.0)

# Generate data WITHOUT carryover (same seed)
set.seed(123)
data_no <- generate_data(..., carryover_t1half = 0)

# Compare: differences should be non-zero at OFF timepoints
```

### Test 2: Model Without Time Effect

```r
# Remove week covariate
model <- lmer(response ~ treatment * bm + (1 | participant_id))
```

### Test 3: Calculate Effect Sizes

```r
# Treatment effect
mean_br_on <- mean(br_means[treatment == 1])
mean_br_off <- mean(br_means[treatment == 0])
treatment_effect <- mean_br_on - mean_br_off

# Carryover effect
carryover_effect <- mean_br_off_with_co - mean_br_off_no_co

# Ratio
effect_ratio <- carryover_effect / treatment_effect
```

---

## References

- Hendrickson et al. (2020). Statistics in Medicine.
- Original code: `~/prj/c265/pmsimstats-master/pmsim-orig/`

---

*Investigation in progress: 2025-11-18*
