# Detailed Comparison: Hendrickson vs. Our Implementation

## Executive Summary

After deep analysis of both implementations, I've identified **critical differences** in:
1. **Biomarker interaction implementation** (data generation)
2. **Carryover modeling in analysis** (statistical model)
3. **Correlation structure computation**

These differences explain why power behaves differently in the two approaches.

---

## Part 1: Mean Construction (Data Generation)

### 1.1 Biological Response Means

#### Hendrickson's Approach (`generateData.R:84-100`)

```r
# Base response from Gompertz (function of time on drug)
brmeans <- modgompertz(d$tod, rp$max, rp$disp, rp$rate)

# Apply carryover sequentially through timepoints
if(nP>1){
  for(p in 2:nP){
    if(!d[p]$onDrug){
      if(d[p]$tsd>0){
        brmeans[p] <- brmeans[p] + brmeans[p-1]*(1/2)^(scalefactor * d$tsd[p]/modelparam$carryover_t1half)
      }
    }
  }
}
```

**Key points**:
- Gompertz function gives response as function of `tod` (time on drug)
- Carryover adds previous period's response × decay factor
- **NO biomarker interaction in means** - interaction is ONLY through correlations

---

#### Our Approach (`pm_functions.R:102-136`)

```r
# Base response from Gompertz
base_bio_response <- mod_gompertz(
  trial_data$tod,
  resp_param$max[resp_param$cat == "bio_response"],
  resp_param$disp[resp_param$cat == "bio_response"],
  resp_param$rate[resp_param$cat == "bio_response"]
)

# ADD POPULATION-LEVEL BIOMARKER INTERACTION TO MEANS
interaction_strength <- model_param$c.bm
biomarker_interaction_effect <- interaction_strength * trial_data$tod * 2.0
bio_response_means <- base_bio_response + biomarker_interaction_effect

# Apply carryover
bio_response_means <- apply_carryover_to_component(
  bio_response_means, trial_data, component_halflives,
  scale_factor, "br"
)
```

**Key points**:
- Same Gompertz base
- **ADDS population mean shift**: `c.bm * tod * 2.0`
- This is NOT participant-specific - it's a constant shift
- Same carryover logic as Hendrickson

---

### 1.2 Post-MVN Biomarker×Treatment Interaction

#### Hendrickson's Approach
**Does NOT apply post-MVN adjustments**. The biomarker effect manifests ONLY through the correlation between baseline biomarker and biological response in the MVN draw.

---

#### Our Approach (`pm_functions.R:643-668`)

```r
# AFTER MVN draw, add participant-specific biomarker interaction
for (br_col in br_columns) {
  timepoint_name <- sub("\\.br$", "", br_col)
  timepoint_idx <- which(trial_design$timepoint_name == timepoint_name)
  treatment_status <- trial_design$tod[timepoint_idx]

  # Add individual biomarker value × treatment × interaction strength
  participant_data[[br_col]] <- participant_data[[br_col]] +
    (participant_data$biomarker * treatment_status * interaction_strength)
}
```

**Key point**: This creates TRUE participant-specific biomarker×treatment interaction.

---

### 1.3 **CRITICAL ISSUE: Double-Counting Biomarker Effect**

Our implementation applies biomarker effect **THREE times**:

1. **Correlation matrix**: Biomarker correlated with bio_response (like Hendrickson) ✓
2. **Population mean shift**: `c.bm * tod * 2.0` added to means (NOT in Hendrickson) ❌
3. **Post-MVN adjustment**: `participant_biomarker * treatment * c.bm` (NOT in Hendrickson) ✓

**Issues**:
- Step 2 is redundant with step 3
- Step 2 inflates means without creating participant-specific effects
- Inflated means affect correlation scaling (see Part 2)
- This could destabilize the sigma matrix

**Recommendation**: Remove step 2 (population mean shift)

---

## Part 2: Correlation Matrix Construction

### 2.1 Biomarker Correlation Scaling

Both implementations use mean-ratio scaling for biomarker correlations when carryover is present.

#### Hendrickson (`generateData.R:148-152`)

```r
if (p > 1) {
  n0 <- paste(trialdesign$timeptname[p - 1], "br", sep = ".")
  mm1 <- means[which(n1 == labels)]  # Current mean (with carryover)
  mm0 <- means[which(n0 == labels)]  # Previous mean

  correlations["bm", n1] <- correlations[n1, "bm"] <-
    ifelse(brtest[p],
           ifelse(brmeans[p] == 0, 0, (mm1 / mm0) * modelparam$c.bm),
           modelparam$c.bm)
}
```

**Logic**:
- If off drug (`brtest[p]` TRUE):
  - If no carryover (`brmeans[p] == 0`): correlation = 0
  - Else: correlation = `(mm1/mm0) * c.bm`
- If on drug: correlation = `c.bm`

---

#### Our Approach (`pm_functions.R:414-429`)

```r
scaled_correlation <- ifelse(
  bio_response_test[timepoint_idx],
  ifelse(
    bio_response_means[timepoint_idx] == 0 ||
      abs(mean_value0) < 1e-10,
    0,
    (mean_value1 / max(mean_value0, 1e-10)) * model_param$c.bm
  ),
  model_param$c.bm
)

# CLAMP to valid range
scaled_correlation <- pmax(-0.99, pmin(0.99, scaled_correlation))
```

**Differences from Hendrickson**:
1. Added safety check for `mean_value0 < 1e-10`
2. **CLAMPING** to [-0.99, 0.99] (Hendrickson doesn't clamp)
3. Otherwise same logic

**Impact of inflated means from Part 1.3**:
- Our means include `c.bm * tod * 2.0` population shift
- This inflates `mm1` values
- Ratio `mm1/mm0` can become very large
- Even with clamping to 0.99, this affects correlation structure
- Could cause subtle differences in sigma matrix

---

### 2.2 Sigma Matrix Caching

#### Hendrickson
- Sigma built **fresh** for each batch of participants
- Called inside `generateData()` which is called once per parameter set per path
- No caching across Monte Carlo iterations

#### Our Approach
- Sigma built **once** per parameter combination
- **Cached** and reused for all Monte Carlo iterations
- Significant performance improvement
- But assumes sigma doesn't need to vary across iterations

**Both approaches are valid**, but caching makes our approach faster.

---

## Part 3: Statistical Model (LME Analysis)

### 3.1 Hendrickson's Model (`lme_analysis.R`)

#### Treatment Variable Construction (lines 107-109)

```r
# Create continuous drug variable Dbc
data.m2[,Dbc:=as.numeric(NA)]
data.m2[Db==TRUE,Dbc:=1]  # On drug: 1
data.m2[Db==FALSE,Dbc:=((1/2)^(op$carryover_scalefactor*tsd/op$carryover_t1half))]
# Off drug: exponential decay
```

**Critical**: This uses `op$carryover_t1half` from **analysis** options, NOT from data generation!

**Default** (lines 62-64):
```r
if(!("carryover_t1half"%in%names(op))){
  op$carryover_t1half=0  # DEFAULTS TO ZERO!
  op$carryover_scalefactor=1
}
```

**Result when `op$carryover_t1half = 0`**:
- `tsd/0 = Infinity`
- `(1/2)^Infinity = 0`
- Therefore: `Dbc = 0` when off drug

**Dbc becomes BINARY** (0 or 1), not continuous!

---

#### Model Formula (line 135)

```r
if(varInDb){
  modelbase=paste0(modelbase,"+Dbc+bm*Dbc")
}
# Final model: Sx ~ bm + t + Dbc + bm*Dbc + (1|ptID)
```

**Components**:
- `Sx` = symptom severity (outcome)
- `bm` = biomarker (continuous)
- `t` = time (week)
- `Dbc` = drug status (binary when op$carryover_t1half=0)
- `bm*Dbc` = biomarker × treatment interaction
- `(1|ptID)` = random intercept

**Key**: `Dbc` is BINARY (not continuous decay) in published simulations.

---

### 3.2 Our Model (`full_pmsim_analysis_hyb_versus_co.R`)

#### Carryover Variable Construction (lines 174-189)

```r
time_since_discontinuation = ifelse(
  treatment == 0,
  cumsum(treatment == 0),  # Counts weeks off treatment
  0
)

carryover_effect = if (params$carryover_t1half > 0) {
  ifelse(
    time_since_discontinuation > 0,
    (1/2)^(params$carryover_scale *
           time_since_discontinuation /
           params$carryover_t1half),
    0
  )
} else {
  0
}
```

**Key**: Uses `params$carryover_t1half` from **data generation**, NOT analysis options.

---

#### Model Formula (lines 204-213)

```r
if (params$carryover_t1half > 0) {
  model <- lmer(
    response ~ treatment * bm + week + carryover_effect + (1 | participant_id)
  )
} else {
  model <- lmer(
    response ~ treatment * bm + week + (1 | participant_id)
  )
}
```

**Components**:
- `response` = outcome
- `treatment` = binary (0/1)
- `bm` = biomarker (continuous)
- `week` = time
- `treatment * bm` = biomarker × treatment interaction
- `carryover_effect` = continuous exponential decay
- `(1|participant_id)` = random intercept

**Key differences from Hendrickson**:
1. `carryover_effect` is **continuous** (not binary)
2. `carryover_effect` matches data generation decay function
3. **Carryover IS modeled** (when t1half > 0)

---

## Part 4: Why Small Carryover Reduces Power in Hendrickson

### 4.1 The Mechanism

In Hendrickson's published simulations:

1. **Data generation**: Carryover t1half = 0, 0.1, or 0.2
   - Real carryover effects exist in data
   - Off-drug periods show decaying residual effects

2. **Analysis model**: `op$carryover_t1half` defaults to 0
   - `Dbc` becomes binary (0 when off drug, 1 when on)
   - Model does NOT capture continuous decay
   - Carryover is an **unmodeled confound**

3. **Impact on power**:
   - **With t1half = 0** (no carryover): Model is correctly specified ✓
   - **With t1half = 0.1** (small carryover):
     - Real data has decay patterns
     - Model treats all off-drug periods as identical (Dbc=0)
     - Adds unexplained variance
     - **Power drops significantly** ❌
   - **With t1half = 0.2** (more carryover):
     - Even more decay variation
     - Even more unexplained variance
     - **Power drops further** ❌

---

### 4.2 Why Even Small Carryover Has Large Impact

#### Model Misspecification

When carryover exists but isn't modeled:

```
True model: Y = β₀ + β₁*bm + β₂*treatment + β₃*bm*treatment + β₄*carryover + ε

Fitted model: Y = β₀ + β₁*bm + β₂*treatment + β₃*bm*treatment + ε'

Where: ε' = β₄*carryover + ε  (carryover absorbed into error term)
```

**Effect**:
- Residual variance ↑ (because carryover in error term)
- Standard errors ↑
- t-statistics ↓
- Power ↓

Even **small** carryover (t1half=0.1) can have **large** impact on power because:
1. Carryover affects MULTIPLE timepoints per participant
2. Creates systematic within-participant correlation in errors
3. Violates independence assumptions
4. Inflates variance estimates

---

### 4.3 Why Our Approach Shows Stable Power

In our simulations:

1. **Data generation**: Carryover t1half = 0, 1.0, or 2.0
   - Real carryover effects in data

2. **Analysis model**: Includes `carryover_effect` covariate
   - Uses **same** decay function as data generation
   - Perfect model specification
   - Carryover is **controlled**, not a confound

3. **Result**:
   - Model successfully removes carryover variance
   - Residual variance stays constant
   - Standard errors stay constant
   - **Power remains stable** ✓

---

## Part 5: Additional Differences

### 5.1 Scale Factor

- **Hendrickson**: Default `scalefactor = 2`
- **Ours**: Default `scale_factor = 2` in build_sigma, but =1 in analysis (from `carryover_scale` param)

**From base_params** (`full_pmsim_analysis_hyb_versus_co.R:493`):
```r
carryover_scale = 1.0
```

**Impact**: Our carryover decays slower (scale=1) vs Hendrickson (scale=2), meaning stronger/longer carryover effects.

---

### 5.2 Autocorrelation Structure

Both use same fixed autocorrelations:
- `c.tv = 0.8`
- `c.pb = 0.8`
- `c.br = 0.8`
- `c.cf1t = 0.2`
- `c.cfct = 0.1`

**No difference here** ✓

---

### 5.3 Design Structures

Both implement:
- 4-path hybrid N-of-1 design
- Traditional crossover (AB vs BA)
- Same phase durations

**No meaningful differences** ✓

---

## Part 6: Summary of Key Differences

| Aspect | Hendrickson | Our Implementation | Impact |
|--------|-------------|-------------------|--------|
| **Biomarker in means** | None | Population shift `c.bm*tod*2.0` | Inflates means ❌ |
| **Post-MVN interaction** | None | `participant_bm * treatment * c.bm` | Creates true interaction ✓ |
| **Carryover in analysis** | NOT modeled (Dbc binary) | MODELED (`carryover_effect`) | Major difference! |
| **Sigma caching** | No | Yes | Performance only ✓ |
| **Correlation clamping** | No | Yes (±0.99) | Stability ✓ |
| **Scale factor** | 2 | 1 (in analysis) | Stronger carryover in ours |
| **Model specification** | Dbc + bm*Dbc | treatment*bm + carryover_effect | Different parameterization |

---

## Part 7: Implications and Recommendations

### 7.1 The Biomarker Triple-Counting Issue

**Problem**: We apply biomarker effect via:
1. Correlation (correct)
2. Population mean shift (questionable)
3. Post-MVN adjustment (correct)

**Recommendation**: **Remove step 2** (population mean shift in `calculate_bio_response_with_interaction`)

**Modified code**:
```r
calculate_bio_response_with_interaction <- function(...) {
  base_bio_response <- mod_gompertz(trial_data$tod, ...)

  # REMOVE THIS:
  # biomarker_interaction_effect <- interaction_strength * trial_data$tod * 2.0
  # bio_response_means <- base_bio_response + biomarker_interaction_effect

  # JUST USE BASE:
  bio_response_means <- base_bio_response

  # Apply carryover
  bio_response_means <- apply_carryover_to_component(...)

  return(...)
}
```

The post-MVN adjustment (step 3) already creates the correct participant-specific interaction.

---

### 7.2 Understanding Power Differences

**Hendrickson's results** (power declines with carryover):
- Caused by **unmodeled confound**
- Not inherent to carryover itself
- Demonstrates cost of model misspecification

**Our results** (power stable):
- Caused by **correct model specification**
- Carryover properly controlled
- Demonstrates benefit of modeling carryover

**Both are informative**, but answer different questions:
- Hendrickson: "What happens if you ignore carryover?"
- Us: "What's the power if you properly model carryover?"

---

### 7.3 Validation Steps

To fully validate our implementation:

1. **Remove population mean shift** in biomarker interaction
2. **Run with carryover_scale = 2** (match Hendrickson)
3. **Compare results both WITH and WITHOUT carryover modeling**:
   - Scenario A: Include `carryover_effect` (current approach)
   - Scenario B: Remove `carryover_effect` (mimic Hendrickson)
4. **Verify Scenario B shows declining power** (like Hendrickson)
5. **Verify Scenario A shows stable power** (current results)

---

## Part 8: Theoretical Justification

### 8.1 Hendrickson's Approach (Correlation-Only)

**Biomarker effect through correlations**:
```
Cov(bm, br) = ρ * σ_bm * σ_br
```

When participant i has biomarker value `bm_i`:
- Draw `br_i` from MVN with correlation ρ
- High `bm_i` → high `br_i` (on average)
- Creates participant-specific treatment response

**Advantage**: Parsimonious, standard MVN approach
**Limitation**: Effect size limited by correlation strength

---

### 8.2 Our Approach (Correlation + Post-Adjustment)

**Step 1**: Draw from MVN with correlation
**Step 2**: Add `bm_i * treatment * c.bm`

**Advantage**: More flexible, can create stronger interactions
**Question**: Is this theoretically justified, or overcomplicated?

**Analysis**: The post-adjustment essentially changes the marginal distribution from MVN to a conditional structure:
```
br_i | bm_i, treatment ~ Normal(μ_br + c.bm * bm_i * treatment, σ²_br)
```

This is a **hierarchical/conditional model**, not pure MVN.

**Is this better?** Depends on biological plausibility:
- If biomarker×treatment interaction is ADDITIVE: Post-adjustment is correct
- If interaction emerges from CORRELATION: Hendrickson's approach is correct

**Recommendation**: Document this design choice explicitly. Consider making it optional.

---

## Conclusion

The main reasons for different power patterns:

1. **Carryover modeling**: We model it, Hendrickson doesn't (in analysis)
2. **Biomarker implementation**: We use post-MVN adjustment, Hendrickson uses correlations only
3. **Scale factor**: We use 1, Hendrickson uses 2

The most important difference is #1 - **carryover modeling in the statistical analysis**. This explains why our power is stable while Hendrickson's declines.

Both approaches are valid for different research questions. Ours demonstrates the benefit of properly modeling carryover; Hendrickson's demonstrates the cost of ignoring it.

---

*Analysis completed: 2025-11-18*
