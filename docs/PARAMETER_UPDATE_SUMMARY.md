# Parameter Update Summary

## Changes Made

Updated simulation parameters to improve power while maintaining positive-definiteness:

### Previous Parameters:
- **Autocorrelations**: c.tv = c.pb = c.br = 0.8 (Hendrickson values)
- **Biomarker correlation**: c.bm = 0.2, 0.34

### New Parameters:
- **Autocorrelations**: c.tv = c.pb = c.br = **0.6** (reduced from 0.8)
- **Biomarker correlation**: c.bm = **0, 0.3, 0.48** (following Hendrickson's 3-level design)

### Rationale:
1. **Lower autocorrelations (0.6)**: Less correlation between timepoints = more independent observations = higher effective sample size = better power
2. **Three c.bm levels**: Matches Hendrickson et al. (2020) design with 0 (no interaction), 0.3 (moderate), and maximum feasible (0.48)
3. **c.bm = 0.48**: Maximum value maintaining positive-definiteness with autocorr = 0.6
   - c.bm = 0.5 fails (non-PD)
   - c.bm = 0.6 fails even with very low autocorrelations (<0.3)

---

## Positive-Definiteness Boundary Testing

### Test Results:

**With autocorr = 0.6:**
- c.bm = 0 ✓ Valid
- c.bm = 0.3 ✓ Valid
- c.bm = 0.48 ✓ Valid (maximum)
- c.bm = 0.49 ✗ Failed
- c.bm = 0.5 ✗ Failed
- c.bm = 0.6 ✗ Failed

**With autocorr = 0.8 (Hendrickson):**
- c.bm = 0.34 ✓ Valid (maximum)
- c.bm = 0.35 ✗ Failed

**For c.bm = 0.6:**
- Would require autocorr < 0.3 (all tested values 0.3-0.6 failed)
- Such low autocorrelations may not be realistic

---

## Simulation Grid

**Total combinations**: 3 biomarker × 3 carryover × 2 designs × 2 approaches = 36 conditions
**Total simulations**: 36 × 20 iterations = **720 runs**

### Parameter Grid:
```
biomarker_correlation: 0, 0.3, 0.48
carryover_t1half: 0, 1.0, 2.0 weeks
design: hybrid, crossover
model_carryover: TRUE, FALSE
n_participants: 70
n_iterations: 20
```

### Fixed Parameters:
```
c.tv = 0.6
c.pb = 0.6
c.br = 0.6
c.cf1t = 0.2
c.cfct = 0.1
treatment_effect = 5.0
```

---

## Validation Status

✅ **All 18 design/parameter combinations produce valid PD matrices:**
- Hybrid + c.bm={0, 0.3, 0.48} + carryover={0, 1.0, 2.0}
- Crossover + c.bm={0, 0.3, 0.48} + carryover={0, 1.0, 2.0}

✅ **All 720 simulations completed successfully**

✅ **Visualizations generated:**
- `power_heatmaps_2x2.png/pdf` (2×2 grid showing all conditions)
- `figure4_equivalent_hendrickson_style.png/pdf` (original 2-panel format)

---

## Current Power Results

### Summary Statistics:

**Overall power range**: 0-0.20 (still low, but improved from previous 0-0.15)

**By biomarker correlation:**
- c.bm = 0: Power = 0% (expected - no interaction)
- c.bm = 0.3: Power = 10-15%
- c.bm = 0.48: Power = 10-20%

**By design:**
- Hybrid: Mean power = 8.3%, Max = 20%
- Crossover: Mean power = 8.3%, Max = 15%

**By analysis approach:**
- WITH carryover model: Mean = 8.1%
- WITHOUT carryover model: Mean = 8.6%

### Key Observations:

1. **c.bm = 0 shows 0% power** (correct - no interaction to detect)
2. **Power increases with c.bm** (0.3 → 0.48)
3. **Hybrid design shows slightly higher max power** (20% vs 15%)
4. **Minimal difference between WITH/WITHOUT carryover modeling** (-0.5% on average)
5. **Power still below conventional threshold** (80%)

---

## Why Power is Still Low

Despite increasing effect size (c.bm) and reducing autocorrelations, power remains low because:

1. **Sample size (N=70)** may still be too small for interaction detection
2. **Treatment effect (5.0)** may be modest relative to variability
3. **Interaction effect** is harder to detect than main effects
4. **20 iterations** may show sampling variation (need more for stable estimates)

---

## Recommendations to Increase Power

### Option 1: Increase Sample Size
```r
n_participants = c(100, 150, 200)
```
- Most direct way to increase power
- Test multiple values to create power curve

### Option 2: Increase Treatment Effect
```r
treatment_effect = c(7.0, 10.0)
```
- Larger main effect → larger interaction
- Consider clinical relevance

### Option 3: Increase Iterations
```r
n_iterations = 100
```
- More stable power estimates
- Takes longer to run (~30 minutes)

### Option 4: Reduce Variability
```r
within_subject_sd = 2.0  # (currently 2.8)
```
- Higher signal-to-noise ratio
- Consider realistic values

### Option 5: Combined Approach
- N = 100
- treatment_effect = 7.0
- n_iterations = 50
- Expected power: 40-60%+ (estimate)

---

## Files Updated

1. **`full_pmsim_analysis_hyb_versus_co.R`**
   - Lines 548-556: autocorrelations 0.8 → 0.6
   - Lines 524-528: biomarker_correlation to c(0, 0.3, 0.48)

2. **`visualize_heatmaps.R`**
   - Line 93: Updated caption to reflect new parameters

3. **`visualize_hendrickson_style.R`**
   - Lines 78-84: biomarker_correlation 0.34 → 0.48
   - Lines 113, 186: Updated titles and labels

---

## Next Steps

1. **Review current results** - Examine heatmaps to understand patterns
2. **Decide on power target** - What is acceptable power for this study?
3. **Choose enhancement strategy** - Sample size, effect size, or both?
4. **Re-run with new parameters** - Test higher N or larger effects
5. **Generate power curves** - Vary N to show relationship

---

*Last updated: 2025-11-18*
