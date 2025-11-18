# N-of-1 Trial Simulation Scripts

## Overview

This directory contains scripts for Monte Carlo simulation comparing Hybrid N-of-1 and Traditional Crossover trial designs, with a focus on statistical power when biomarker×treatment interactions exist and carryover effects are present.

---

## Quick Start

```bash
# Navigate to scripts directory
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts

# Step 1: Run simulations (~5-15 minutes)
Rscript full_pmsim_analysis_hyb_versus_co.R

# Step 2: Generate visualizations (~10 seconds)
Rscript visualize_hendrickson_style.R

# Step 3: View results
open ../output/figure4_equivalent_hendrickson_style.pdf
```

---

## Core Files

### Main Simulation

**`full_pmsim_analysis_hyb_versus_co.R`**
- Monte Carlo simulation comparing Hybrid vs. Crossover designs
- Two-scenario comparison: WITH vs. WITHOUT carryover modeling
- Generates 480 simulation runs (24 conditions × 20 iterations)
- **Output:** `../output/full_pmsim_analysis_hyb_versus_co.RData`

**Key features:**
- Aligned with Hendrickson et al. (2020) methodology
- Pure correlation-based biomarker interaction (no post-MVN adjustment)
- Fixed Hendrickson correlation parameters
- Exponential carryover decay with half-life
- Mixed-effects model analysis with random intercepts

### Core Functions

**`pm_functions.R`**
- `build_sigma_matrix()` - Constructs covariance matrices with PD validation
- `generate_data()` - MVN data generation for trial participants
- `calculate_bio_response_with_interaction()` - Computes BR means
- `validate_correlation_structure()` - Matrix validation utility

### Visualization

**`visualize_hendrickson_style.R`**
- Creates Figure 4-equivalent matching Hendrickson's format
- Shows both analysis scenarios side-by-side
- **Panel A:** No carryover, varying biomarker correlation
- **Panel B:** With carryover, varying half-life (fixed correlation = 0.4)
- **Output:**
  - `../output/figure4_equivalent_hendrickson_style.png` (300 DPI)
  - `../output/figure4_equivalent_hendrickson_style.pdf` (vector)

---

## Workflow

### 1. Run Main Simulation

```bash
Rscript full_pmsim_analysis_hyb_versus_co.R
```

**What it does:**
1. Builds sigma matrix cache (validates 12 parameter combinations)
2. Runs Monte Carlo simulation for each combination
3. Analyzes data with mixed-effects models (lmer)
4. Saves results to `../output/full_pmsim_analysis_hyb_versus_co.RData`

**Expected output:**
```
Building sigma matrix cache...
Sigma matrix positive definite: TRUE
  Eigenvalue range: 0.0234 to 2.156
✓ Valid sigma for hybrid design, biomarker_correlation = 0.2, carryover_t1half = 0
...
✓ Valid sigma for crossover design, biomarker_correlation = 0.4, carryover_t1half = 2.0

Sigma cache built. Valid combinations: 12

Running hybrid design with WITH_CARRYOVER_MODEL - params: 1 of 6
...

✓ Simulation complete!
Results saved to: ../output/full_pmsim_analysis_hyb_versus_co.RData
```

**Runtime:** ~5-15 minutes (480 simulations)

### 2. Generate Visualizations

```bash
Rscript visualize_hendrickson_style.R
```

**Prerequisites:** Must run Step 1 first!

**What it does:**
1. Loads simulation results from RData file
2. Creates two-panel figure (Panel A: no carryover, Panel B: with carryover)
3. Shows power heatmaps for both analysis approaches
4. Prints summary statistics and power loss calculations

**Expected output:**
```
✓ Figure 4-equivalent saved to:
  - ../output/figure4_equivalent_hendrickson_style.png
  - ../output/figure4_equivalent_hendrickson_style.pdf

================================================================================
SUMMARY STATISTICS
================================================================================
...

✓ Visualization complete!
```

**Runtime:** ~10 seconds

### 3. Examine Results

**View figures:**
```bash
# PDF (recommended)
open ../output/figure4_equivalent_hendrickson_style.pdf

# PNG
open ../output/figure4_equivalent_hendrickson_style.png
```

**Load in R:**
```r
load("../output/full_pmsim_analysis_hyb_versus_co.RData")
View(simulation_summary)  # Power by condition
View(simulation_results)  # Individual iterations
```

---

## Parameter Configuration

### Current Simulation Grid

Located in `full_pmsim_analysis_hyb_versus_co.R` (lines ~518-526):

```r
param_grid <- expand_grid(
  n_participants = c(70),
  biomarker_correlation = c(0.2, 0.4),
  carryover_t1half = c(0, 1.0, 2.0),
  treatment_effect = c(5.0)
)

n_iterations = 20
```

**Total simulations:**
- 1 sample size × 2 correlations × 3 carryover levels = 6 combinations
- × 2 designs (hybrid, crossover) = 12 combinations
- × 2 analysis approaches (WITH/WITHOUT carryover model) = 24 combinations
- × 20 iterations = **480 total simulations**

### Fixed Parameters

Located in `pm_functions.R` (lines ~38-90):

**Correlation structure (Hendrickson values):**
```r
c.tv = 0.8      # Time-variant autocorrelation
c.pb = 0.8      # Pharmacologic biomarker autocorrelation
c.br = 0.8      # Biological response autocorrelation
c.cf1t = 0.2    # Same-time cross-correlation
c.cfct = 0.1    # Different-time cross-correlation
```

**Response parameters:**
```r
variance.tv = 2.25
variance.pb = 2.25
variance.br = 2.25
treatment_effect = 5.0
baseline_mean = 15.0
```

**Carryover parameters:**
```r
scale_factor = 1.0    # BR-only carryover (affects biological response only)
n_weeks_on = 4        # Treatment period duration
n_weeks_off = 4       # Washout period duration
```

---

## Key Methodological Features

### 1. Two-Scenario Comparison

**Scenario 1: WITH Carryover Model**
- Analysis model includes `carryover_effect` term
- Controls for carryover as a covariate
- Our methodological enhancement over Hendrickson

**Scenario 2: WITHOUT Carryover Model**
- Analysis model does NOT include carryover term
- Replicates Hendrickson's original approach
- Carryover becomes unmodeled confound

**Purpose:** Quantify the cost of not modeling carryover effects

### 2. Pure Correlation-Based Biomarker Interaction

Following Option 1 alignment with Hendrickson:
- Biomarker×treatment interaction emerges from differential correlation
- NO population mean shift based on biomarker
- NO post-MVN adjustment
- Correlation `c.bm` varies by treatment status (on/off)

### 3. Quality Control

**Sigma matrix validation:**
- Rejects non-positive definite matrices (doesn't auto-fix)
- Biomarker correlation scaling with clamping to [-0.99, 0.99]
- Eigenvalue diagnostics printed during cache build

**Convergence tracking:**
- Mixed model convergence failures marked in results
- Error column in `simulation_results` tracks failed fits

---

## Output Files

### Generated by Simulation

**`../output/full_pmsim_analysis_hyb_versus_co.RData`**

Contains:
- `simulation_results` - Detailed results (480 rows)
  - Columns: `design`, `iteration`, `n_participants`, `biomarker_correlation`, `carryover_t1half`, `model_carryover`, `effect_size`, `std_error`, `t_value`, `p_value`, `significant`, `error`

- `simulation_summary` - Aggregated power (24 rows)
  - Columns: `design`, `n_participants`, `biomarker_correlation`, `carryover_t1half`, `model_carryover`, `power`, `mean_effect`, `mean_se`, `n_iterations`

- `param_grid` - Parameter combinations tested
- `model_params`, `resp_param`, `baseline_param` - Fixed parameters

### Generated by Visualization

**`../output/figure4_equivalent_hendrickson_style.png`**
- Raster format, 300 DPI
- 14 × 10 inches
- Good for presentations/web

**`../output/figure4_equivalent_hendrickson_style.pdf`**
- Vector format
- Publication-ready
- Recommended for manuscripts

---

## Modifying the Simulation

### Change Sample Sizes

Edit `full_pmsim_analysis_hyb_versus_co.R` line ~520:
```r
n_participants = c(50, 70, 90),  # Add multiple values
```

### Increase Iterations

Edit line ~630:
```r
n_iterations = 100  # More stable estimates (but slower)
```

### Add Biomarker Correlations

Edit line ~522:
```r
biomarker_correlation = c(0.2, 0.3, 0.4, 0.5),
```

**Caution:** Higher correlations (>0.5) may approach matrix stability limits

### Change Carryover Gradient

Edit line ~523:
```r
carryover_t1half = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
```

### Modify Fixed Correlations

Edit `pm_functions.R` lines ~60-66:
```r
c.tv = 0.8,
c.pb = 0.8,
c.br = 0.8,
c.cf1t = 0.2,
c.cfct = 0.1,
```

**Warning:** Changing these may affect matrix positive-definiteness. Check sigma cache output!

---

## Troubleshooting

### "Cannot find pm_functions.R"

**Solution:**
```bash
# Ensure you're in the correct directory
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
pwd  # Should show scripts directory
```

### "Cannot open RData file"

**Cause:** Simulation hasn't been run yet

**Solution:**
```bash
# Run simulation first
Rscript full_pmsim_analysis_hyb_versus_co.R
# Then run visualization
Rscript visualize_hendrickson_style.R
```

### "Non-positive definite matrix rejected"

**Cause:** Parameter combination creates invalid correlation matrix

**Expected:** This is quality control - invalid combinations are skipped

**Action:** Check console for which parameters failed. May need to reduce correlation or carryover values.

### Mixed model convergence warnings

**Expected:** Occasional warnings are normal with random data

**Monitor:** If >20% of iterations fail, review parameter choices

**Check:** `simulation_results$error` column shows which runs failed

---

## Performance

**Typical runtime (M1 Mac, 16GB RAM):**
- Sigma cache build: ~10 seconds
- Monte Carlo simulation: ~5-15 minutes (480 runs)
- Visualization: ~10 seconds
- **Total:** ~5-15 minutes

**Scaling:**
- Runtime scales linearly with iterations
- 100 iterations ≈ 25 minutes
- 1000 iterations ≈ 4 hours

---

## Related Documentation

### Comprehensive Guides

- **`../../docs/SIMULATION_WORKFLOW.md`** - Complete workflow documentation
- **`../../docs/TWO_SCENARIO_COMPARISON.md`** - Two-scenario design rationale
- **`../../docs/PSEUDOCODE_COMPARISON_HENDRICKSON_VS_OURS.md`** - Implementation comparison

### Technical References

- **`../../docs/technical_differences_scaling_and_pd.pdf`** - Matrix scaling and PD handling
- **`../../docs/WHY_POST_MVN_ADJUSTMENT_ANALYSIS.md`** - Rationale for removing post-MVN adjustment
- **`../../docs/HENDRICKSON_VS_OUR_DETAILED_COMPARISON.md`** - Detailed implementation differences

### Theory

- **`../../docs/carryover_correlation_theory.pdf`** - Mathematical foundations
- **`../../docs/correlation_structure_design.pdf`** - Correlation matrix construction

---

## Design Details

### Hybrid Design (4-Path Randomization)

Four balanced paths:
1. **Path 1:** A-B-A-B (alternating, starts A)
2. **Path 2:** B-A-B-A (alternating, starts B)
3. **Path 3:** A-A-B-B (blocks, starts A)
4. **Path 4:** B-B-A-A (blocks, starts B)

Each participant assigned to one path; 4 periods total.

### Crossover Design (2-Sequence)

Two sequences:
1. **Sequence 1 (AB):** A-B
2. **Sequence 2 (BA):** B-A

Participants alternately assigned; 2 periods total.

### Common Structure

- **Treatment A:** Active treatment (tod = 1)
- **Treatment B:** Control/placebo (tod = 0)
- **Each period:** 4 weeks on treatment + 4 weeks washout
- **Biomarker:** Measured at baseline (average of PB across timepoints)
- **Response:** Biological response (BR) measured each period

---

## Analysis Model

### Mixed-Effects Model Specification

**Scenario 1 (WITH carryover):**
```r
lmer(response ~ treatment * biomarker + week + carryover_effect +
     (1 | participant_id))
```

**Scenario 2 (WITHOUT carryover):**
```r
lmer(response ~ treatment * biomarker + week +
     (1 | participant_id))
```

**Key terms:**
- `treatment * biomarker` - Biomarker×treatment interaction (tested)
- `week` - Time effect (period/progression)
- `carryover_effect` - Exponential decay: (1/2)^(weeks_off/t1half)
- `(1 | participant_id)` - Random intercept

**Primary outcome:** Statistical power to detect `treatment:biomarker` interaction (p < 0.05)

---

## Citation

If using these scripts, please cite:

**Hendrickson, E., et al. (2020).** N-of-1 trials with multiple randomization structures for individualized treatment. *Statistics in Medicine*, 39(25), 3581-3599.

---

## Contact

For questions about this implementation:
- See documentation in `../../docs/`
- Check GitHub repository issues
- Review pseudocode comparison for implementation details

---

*Last updated: 2025-11-18*
