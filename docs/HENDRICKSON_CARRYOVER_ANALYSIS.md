# Hendrickson Carryover Implementation Analysis

## Summary

After examining the original Hendrickson repository (`~/prj/c265/pmsimstats-master`), I found that:

1. **Hendrickson DID vary carryover in data generation** - but **DID NOT model it in the statistical analysis**
2. This is a **critical difference** from our current implementation
3. This likely explains why power patterns differ from my initial (incorrect) intuition

---

## Detailed Findings

### 1. Data Generation - Carryover IS Included

**File**: `vignettes/Produce_Publication_Results_1_generate_data.Rmd` (lines 315-321)

```r
coremodelparams<-expand.grid(
  N=c(35,70),
  c.bm=c(0,.3,.6),
  carryover_t1half=c(0,.1,.2),  # <-- Carryover VARIED in data generation
  c.tv=.8,c.pb=.8,c.br=.8,
  c.cf1t=.2,c.cfct=.1
)
```

Hendrickson tested three carryover conditions:
- `carryover_t1half = 0` (no carryover)
- `carryover_t1half = 0.1` (short half-life)
- `carryover_t1half = 0.2` (longer half-life)

These parameters are passed to `generateData()` which implements carryover effects in the biological response means.

---

### 2. Statistical Analysis - Carryover is NOT Modeled

**File**: `vignettes/Produce_Publication_Results_1_generate_data.Rmd` (line 350)

```r
# Parameters on the analysis side
analysisparams<-expand.grid(useDE=FALSE,t_random_slope=FALSE,full_model_out=FALSE)
```

**Notice**: `carryover_t1half` is **NOT included** in `analysisparams`!

---

### 3. What Happens in lme_analysis()

**File**: `R/lme_analysis.R` (lines 62-64, 107-109)

```r
# Default carryover parameters when not specified
if(!("carryover_t1half"%in%names(op))){
  op$carryover_t1half=0
  op$carryover_scalefactor=1
}

# Create Dbc variable (continuous drug indicator with carryover)
data.m2[,Dbc:=as.numeric(NA)]
data.m2[Db==TRUE,Dbc:=1]  # On drug: Dbc = 1
data.m2[Db==FALSE,Dbc:=((1/2)^(op$carryover_scalefactor*tsd/op$carryover_t1half))]
```

**When `op$carryover_t1half = 0` (the default)**:
- `tsd/0 = Infinity`
- `(1/2)^Infinity = 0`
- Therefore: `Dbc = 0` when off drug

**Result**: `Dbc` becomes a **binary 0/1 variable** (same as `Db`), with NO gradient for carryover decay.

---

### 4. Statistical Model Used

**File**: `R/lme_analysis.R` (line 135)

```r
if(varInDb){
  modelbase=paste0(modelbase,"+Dbc+bm*Dbc")
}
```

**Final model** (typical case):
```r
Sx ~ bm + t + Dbc + bm*Dbc + (1|ptID)
```

Where:
- `Sx` = symptom severity (outcome)
- `bm` = biomarker
- `t` = time
- `Dbc` = drug indicator (binary, despite the name suggesting continuous)
- `bm*Dbc` = **biomarker × treatment interaction** (primary effect of interest)
- `(1|ptID)` = random intercept

**Crucially**: `Dbc` does NOT capture carryover decay because `op$carryover_t1half = 0`.

---

## Comparison: Hendrickson vs. Our Implementation

| Aspect | Hendrickson (2020) | Our Implementation |
|--------|-------------------|-------------------|
| **Data Generation** | | |
| Carryover in data? | ✓ Yes (t1half = 0, 0.1, 0.2) | ✓ Yes (t1half = 0, 1.0, 2.0) |
| Carryover affects means? | ✓ Yes | ✓ Yes |
| Carryover affects correlations? | ✗ No (fixed) | ✗ No (fixed) |
| **Statistical Analysis** | | |
| Carryover modeled? | ✗ **NO** (defaults to 0) | ✓ **YES** (explicitly included) |
| Model includes Dbc? | ✓ Yes (but binary) | ✓ Yes (continuous with decay) |
| Model includes carryover_effect? | ✗ No | ✓ **YES** |
| Model includes bm*Dbc interaction? | ✓ Yes | ✓ Yes (as treatment*bm) |

---

## Implications for Power Analysis

### Hendrickson's Approach (Carryover NOT Modeled)

When carryover exists in data but is NOT modeled:

1. **Increased Noise**: Carryover adds unexplained variance
2. **Confounding**: Off-drug periods show residual effects not captured by model
3. **Expected Pattern**: Power should **DECLINE** as carryover increases

### Our Approach (Carryover IS Modeled)

When carryover exists in data AND is modeled:

1. **Controlled Confound**: `carryover_effect` removes confounding
2. **Reduced Noise**: Model accounts for decay pattern
3. **Expected Pattern**: Power should **REMAIN STABLE** across carryover levels

---

## Why Our Results Show Stable Power

Our simulation shows power does NOT decline with increasing carryover because:

1. **Perfect Control**: The `carryover_effect` variable in our model uses the **exact same** decay function as data generation:
   ```r
   # Data generation (pm_functions.R:262)
   decay_factor <- (1/2)^(scale_factor * time_lag / component_halflives$br)

   # Analysis model (full_pmsim_analysis_hyb_versus_co.R:183-185)
   carryover_effect = (1/2)^(params$carryover_scale *
                             time_since_discontinuation /
                             params$carryover_t1half)
   ```

2. **Model Specification** (full_pmsim_analysis_hyb_versus_co.R:205-206):
   ```r
   model <- lmer(
     response ~ treatment * bm + week + carryover_effect + (1 | participant_id)
   )
   ```

3. **Result**: The model successfully removes the confounding effect of carryover, maintaining power regardless of carryover strength.

---

## What Would Happen If We Used Hendrickson's Approach?

If we **removed** `carryover_effect` from our model (matching Hendrickson's analysis):

```r
# Hypothetical: NO carryover modeling (like Hendrickson)
model <- lmer(
  response ~ treatment * bm + week + (1 | participant_id)
)
```

**Expected results**:
- Power would **DECLINE** as `carryover_t1half` increases (0 → 1.0 → 2.0)
- Hybrid design would be MORE affected (more discontinuation periods)
- Standard errors would increase
- This would match the pattern I initially (incorrectly) predicted

---

## Key Insight: Model Misspecification

**Hendrickson's simulations contained an unmodeled confound:**

- **Data**: Includes carryover effects (exponential decay after discontinuation)
- **Model**: Does NOT account for carryover (Dbc is binary, not continuous)
- **Impact**: This adds noise and reduces power as carryover increases

**Our enhancement:**
- Explicitly models carryover in BOTH data generation AND analysis
- This is methodologically superior (no unmodeled confounds)
- But it changes the power dynamics compared to Hendrickson's approach

---

## Recommendation: Test Both Scenarios

To fully understand the impact of carryover, we should run TWO versions:

### Scenario 1: With Carryover Modeling (Current)
```r
# Model includes carryover control
model <- lmer(response ~ treatment * bm + week + carryover_effect + (1|participant_id))
```
**Expected**: Stable power across carryover levels

### Scenario 2: Without Carryover Modeling (Hendrickson-style)
```r
# Model does NOT include carryover (unmodeled confound)
model <- lmer(response ~ treatment * bm + week + (1|participant_id))
```
**Expected**: Declining power as carryover increases

This would allow us to:
1. Demonstrate the **cost of not modeling carryover** (Scenario 2)
2. Show the **benefit of our enhancement** (Scenario 1 vs. 2)
3. Understand both methodological approaches

---

## Additional Consideration: Biomarker-Dependent Carryover

There's a subtle issue in BOTH implementations:

**In data generation**, carryover magnitude depends on biomarker:
```r
# pm_functions.R:115-120
biomarker_interaction_effect <- interaction_strength * trial_data$tod * 2.0
bio_response_means <- base_bio_response + biomarker_interaction_effect

# Then carryover is applied to bio_response_means
# Which INCLUDES the biomarker interaction!
```

**In analysis models** (both Hendrickson's and ours):
- Carryover is modeled as **independent** of biomarker
- Hendrickson: `Dbc` (no interaction with bm for carryover)
- Ours: `carryover_effect` (no interaction with bm)

**Technically correct model** would be:
```r
response ~ treatment * bm + week + carryover_effect + carryover_effect * bm + (1|participant_id)
```

Or the carryover variable should capture the actual lagged response (which already includes biomarker effects).

---

## Conclusion

**Hendrickson did NOT model carryover effects in her statistical analysis**, even though she varied carryover parameters in data generation. This was an unmodeled confound that would reduce power as carryover increased.

**Our implementation DOES model carryover**, which is methodologically superior but creates different power dynamics. This explains why our power remains stable across carryover levels rather than declining.

This is not a bug - it's a **feature enhancement** over the original Hendrickson approach.

---

## Files Examined

1. `~/prj/c265/pmsimstats-master/vignettes/Produce_Publication_Results_1_generate_data.Rmd`
2. `~/prj/c265/pmsimstats-master/R/lme_analysis.R`
3. `~/prj/c265/pmsimstats-master/R/generateSimulatedResults.R`
4. `~/prj/c265/pmsimstats-master/R/generateData.R`

---

*Analysis completed: 2025-11-18*
