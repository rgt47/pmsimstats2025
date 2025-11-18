# Simulation Workflow Guide

## Overview

This document describes the complete workflow for running Monte Carlo simulations and generating visualizations for the N-of-1 trial comparison project.

---

## Quick Start (3 Steps)

```bash
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts

# Step 1: Run simulations
Rscript full_pmsim_analysis_hyb_versus_co.R

# Step 2: Generate visualizations
Rscript visualize_hendrickson_style.R

# Step 3: View results
open ../output/figure4_equivalent_hendrickson_style.pdf
```

---

## Detailed Workflow

### Prerequisites

**Required R packages:**
```r
tidyverse      # Data manipulation and ggplot2
lmerTest       # Mixed-effects models
viridis        # Color scales
MASS           # mvrnorm for MVN sampling
corpcor        # Matrix operations
patchwork      # Combining plots (for visualizations)
```

**Install if needed:**
```r
install.packages(c("tidyverse", "lmerTest", "viridis", "MASS", "corpcor", "patchwork"))
```

---

## Step 1: Run Main Simulation

### What It Does

The main simulation script (`full_pmsim_analysis_hyb_versus_co.R`) performs:

1. **Sigma matrix cache building** - Pre-validates all parameter combinations
2. **Monte Carlo simulation** - Generates trial data and analyzes it
3. **Two-scenario comparison** - Runs WITH and WITHOUT carryover modeling
4. **Results storage** - Saves detailed and summary results to RData file

### Current Parameter Grid

```r
param_grid <- expand_grid(
  n_participants = c(70),
  biomarker_correlation = c(0.2, 0.4),
  carryover_t1half = c(0, 1.0, 2.0),
  treatment_effect = c(5.0)
)

# Total combinations:
# 1 sample size × 2 biomarker correlations × 3 carryover levels = 6
# × 2 designs (hybrid, crossover) = 12
# × 2 analysis approaches (WITH/WITHOUT carryover) = 24
# × 20 iterations per combination = 480 total simulations
```

### Run the Simulation

**Option A: Command line**
```bash
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
Rscript full_pmsim_analysis_hyb_versus_co.R
```

**Option B: RStudio**
```r
setwd("/Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts")
source("full_pmsim_analysis_hyb_versus_co.R")
```

**Option C: Interactive R console**
```r
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
R
> source("full_pmsim_analysis_hyb_versus_co.R")
```

### Expected Runtime

- **Sigma cache building**: ~10 seconds (validates 12 combinations)
- **Monte Carlo simulation**: ~5-15 minutes (480 runs with mixed models)
- **Total**: ~5-15 minutes depending on system

### Output Files

**Primary output:**
```
analysis/output/full_pmsim_analysis_hyb_versus_co.RData
```

**Contents:**
- `simulation_results` - Detailed results (one row per iteration)
- `simulation_summary` - Aggregated power estimates
- `param_grid` - Parameter combinations tested
- `model_params`, `resp_param`, `baseline_param` - Fixed parameters

### Console Output

During execution, you'll see:
```
Building sigma matrix cache...
Sigma matrix positive definite: TRUE
  Eigenvalue range: 0.0234 to 2.156
✓ Valid sigma for hybrid design, biomarker_correlation = 0.2, carryover_t1half = 0
...
✓ Valid sigma for crossover design, biomarker_correlation = 0.4, carryover_t1half = 2.0

Sigma cache built. Valid combinations: 12

Running hybrid design with WITH_CARRYOVER_MODEL - params: 1 of 6
Running crossover design with WITH_CARRYOVER_MODEL - params: 1 of 6
Running hybrid design with NO_CARRYOVER_MODEL - params: 1 of 6
...

✓ Simulation complete!
Results saved to: analysis/output/full_pmsim_analysis_hyb_versus_co.RData
```

### Troubleshooting

**Problem: "Output directory not found"**
```bash
mkdir -p /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/output
```

**Problem: "pm_functions.R not found"**
```bash
# Make sure you're in the correct directory
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
```

**Problem: Package not installed**
```r
install.packages("package_name")
```

**Problem: Non-convergence warnings from lmer**
- Expected occasionally with random data
- Script continues and marks convergence failures
- Check `simulation_results$error` column

---

## Step 2: Generate Visualizations

### What It Does

The visualization script (`visualize_hendrickson_style.R`) creates:

1. **Panel A**: No carryover (t1/2 = 0), varying biomarker correlation
2. **Panel B**: With carryover, fixed biomarker correlation = 0.4
3. **Both panels** show: Hybrid vs. Crossover × WITH vs. WITHOUT carryover modeling

### Run Visualization

**Option A: Command line**
```bash
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
Rscript visualize_hendrickson_style.R
```

**Option B: RStudio/R console**
```r
setwd("/Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts")
source("visualize_hendrickson_style.R")
```

### Prerequisites

**Must run Step 1 first!** The script loads:
```r
load("analysis/output/full_pmsim_analysis_hyb_versus_co.RData")
```

**Additional package needed:**
```r
install.packages("patchwork")  # For combining plots
```

### Output Files

**Figures:**
```
analysis/output/figure4_equivalent_hendrickson_style.png  (300 DPI)
analysis/output/figure4_equivalent_hendrickson_style.pdf  (vector)
```

**Dimensions:**
- Width: 14 inches
- Height: 10 inches
- Format: Publication-ready

### Console Output

```
✓ Figure 4-equivalent saved to:
  - analysis/output/figure4_equivalent_hendrickson_style.png
  - analysis/output/figure4_equivalent_hendrickson_style.pdf

================================================================================
SUMMARY STATISTICS
================================================================================

Panel A (No Carryover, carryover_t1half = 0):
--------------------------------------------------------------------------------
# A tibble: 8 × 4
  design_label biomarker_correlation approach_label          power
  <chr>                        <dbl> <chr>                   <dbl>
1 Crossover                      0.2 WITHOUT Carryover Model 0.75
2 Crossover                      0.2 WITH Carryover Model    0.75
3 Crossover                      0.4 WITHOUT Carryover Model 0.85
...

Panel B (With Carryover, biomarker_correlation = 0.4):
--------------------------------------------------------------------------------
# A tibble: 8 × 4
  design_label carryover_t1half approach_label          power
  <chr>                   <dbl> <chr>                   <dbl>
1 Crossover                 1   WITHOUT Carryover Model 0.70
2 Crossover                 1   WITH Carryover Model    0.83
...

================================================================================
KEY FINDINGS
================================================================================

Power Lost by NOT Modeling Carryover:
--------------------------------------------------------------------------------
# A tibble: 8 × 6
  design    biomarker_correlation carryover_t1half `WITH Model` `WITHOUT Model` power_loss
  <chr>                     <dbl>            <dbl>        <dbl>           <dbl>      <dbl>
1 crossover                   0.2              1           0.72            0.65       0.07
2 crossover                   0.2              2           0.70            0.58       0.12
...

✓ Visualization complete!
```

---

## Step 3: View and Interpret Results

### View Figures

**PDF (recommended for quality):**
```bash
open /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/output/figure4_equivalent_hendrickson_style.pdf
```

**PNG (for presentations/web):**
```bash
open /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/output/figure4_equivalent_hendrickson_style.png
```

### Examine Results in R

```r
# Load results
load("/Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/output/full_pmsim_analysis_hyb_versus_co.RData")

# View summary statistics
View(simulation_summary)

# View detailed results
View(simulation_results)

# Check specific scenario
simulation_summary %>%
  filter(design == "hybrid",
         carryover_t1half == 1.0,
         biomarker_correlation == 0.4)

# Calculate power loss
power_comparison <- simulation_summary %>%
  filter(carryover_t1half > 0) %>%
  select(design, carryover_t1half, biomarker_correlation, model_carryover, power) %>%
  pivot_wider(names_from = model_carryover,
              values_from = power,
              names_prefix = "model_") %>%
  mutate(power_loss = model_TRUE - model_FALSE)

View(power_comparison)
```

### Key Metrics to Examine

**From `simulation_summary`:**
- `power` - Main outcome (proportion of significant results)
- `mean_effect` - Average biomarker×treatment interaction coefficient
- `mean_se` - Average standard error (higher = more noise)
- `n_iterations` - Number of successful runs (should be 20)

**From `simulation_results`:**
- `effect_size` - Individual iteration's interaction estimate
- `p_value` - Statistical significance
- `significant` - TRUE/FALSE for each iteration
- `error` - TRUE if convergence failed

---

## Interpretation Guide

### Panel A: No Carryover (Baseline)

**What to look for:**
- Both approaches (WITH/WITHOUT) should show SAME power
- No carryover in data → no difference in models
- Validates implementation alignment

**Expected pattern:**
```
biomarker_correlation = 0.2 → Power ~0.60-0.75
biomarker_correlation = 0.4 → Power ~0.80-0.90
```

Higher biomarker correlation → higher power (stronger interaction signal)

### Panel B: With Carryover

**What to look for:**
- WITHOUT modeling: Power DECLINES as carryover increases
- WITH modeling: Power STABLE across carryover levels
- Difference = cost of not modeling carryover

**Expected pattern (WITHOUT modeling):**
```
carryover_t1half = 0   → Power ~0.80
carryover_t1half = 1.0 → Power ~0.70  (decline!)
carryover_t1half = 2.0 → Power ~0.60  (larger decline!)
```

**Expected pattern (WITH modeling):**
```
carryover_t1half = 0   → Power ~0.80
carryover_t1half = 1.0 → Power ~0.78  (stable)
carryover_t1half = 2.0 → Power ~0.76  (slight decline from DF loss)
```

### Design Comparison

**Hybrid vs. Crossover:**
- Compare within each analysis approach
- Which design maintains power better?
- Which is more robust to carryover?

---

## Modifying the Simulation

### Change Sample Size

Edit `full_pmsim_analysis_hyb_versus_co.R` line ~520:
```r
param_grid <- expand_grid(
  n_participants = c(50, 70, 90),  # Add multiple values
  biomarker_correlation = c(0.2, 0.4),
  carryover_t1half = c(0, 1.0, 2.0),
  treatment_effect = c(5.0)
)
```

### Change Number of Iterations

Edit line ~630:
```r
n_iterations = 100  # Increase for more stable estimates (slower)
```

**Trade-off:**
- More iterations → more accurate power estimates
- More iterations → longer runtime
- Recommendation: 20 for development, 1000 for publication

### Add More Biomarker Correlations

Edit line ~520:
```r
biomarker_correlation = c(0.2, 0.3, 0.4, 0.5, 0.6),
```

**Note:** Higher correlations (>0.5) may approach matrix stability limits

### Change Carryover Half-lives

Edit line ~520:
```r
carryover_t1half = c(0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
```

Creates finer gradient to see power decline pattern

### Modify Fixed Parameters

Edit lines 38-90 in `pm_functions.R`:
```r
model_params <- list(
  c.tv = 0.8,      # Time-variant autocorrelation
  c.pb = 0.8,      # PB autocorrelation
  c.br = 0.8,      # BR autocorrelation
  c.bm = 0.4,      # Will be overridden by param_grid
  c.cf1t = 0.2,    # Same-time cross-correlation
  c.cfct = 0.1,    # Different-time cross-correlation
  N = 70           # Will be overridden
)

resp_param <- list(
  variance.tv = 2.25,
  variance.pb = 2.25,
  variance.br = 2.25,
  treatment_effect = 5.0,
  baseline_mean = 15.0
)

baseline_param <- list(
  variance.ue = 0.8,
  variance.ce = 0.0,
  t1half = 1.0,
  scale_factor = 1.0,  # BR-only carryover
  n_weeks_on = 4,
  n_weeks_off = 4
)
```

**Caution:** Changing these may affect matrix stability. Always check sigma cache output.

---

## Complete Example Session

```bash
# Start fresh
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts

# Ensure output directory exists
mkdir -p ../output

# Run simulation (5-15 minutes)
echo "Running simulation..."
Rscript full_pmsim_analysis_hyb_versus_co.R

# Generate visualizations (~10 seconds)
echo "Generating visualizations..."
Rscript visualize_hendrickson_style.R

# View results
echo "Opening results..."
open ../output/figure4_equivalent_hendrickson_style.pdf

# Examine in R
R --no-save <<EOF
load("../output/full_pmsim_analysis_hyb_versus_co.RData")
print("Summary statistics:")
print(simulation_summary)
print("\nParameter grid:")
print(param_grid)
EOF
```

---

## File Structure

```
pmsimstats2025/
├── analysis/
│   ├── scripts/
│   │   ├── pm_functions.R                        # Core functions
│   │   ├── full_pmsim_analysis_hyb_versus_co.R  # Main simulation
│   │   └── visualize_hendrickson_style.R        # Visualization
│   └── output/
│       ├── full_pmsim_analysis_hyb_versus_co.RData        # Results
│       ├── figure4_equivalent_hendrickson_style.png       # Figure (PNG)
│       └── figure4_equivalent_hendrickson_style.pdf       # Figure (PDF)
└── docs/
    ├── SIMULATION_WORKFLOW.md                    # This file
    ├── TWO_SCENARIO_COMPARISON.md                # Design documentation
    ├── PSEUDOCODE_COMPARISON_HENDRICKSON_VS_OURS.md
    └── technical_differences_scaling_and_pd.pdf
```

---

## Advanced: Parallel Processing

For large parameter grids, consider parallelization:

```r
# Edit full_pmsim_analysis_hyb_versus_co.R
library(parallel)
library(doParallel)

# Set up parallel backend
n_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Run iterations in parallel (modify run_monte_carlo)
# ... implementation details ...

stopCluster(cl)
```

**Note:** Parallel version not yet implemented. Current version is sequential.

---

## Troubleshooting Common Issues

### Issue: "Cannot find file pm_functions.R"

**Solution:**
```bash
# Check you're in the right directory
pwd
# Should show: /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts

# If not:
cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
```

### Issue: "Non-positive definite matrix rejected"

**Cause:** Parameter combination creates invalid correlation matrix

**Solution:**
- Check console output for which parameters failed
- Reduce biomarker correlation or carryover strength
- This is expected behavior (quality control)

### Issue: Convergence warnings from lmer

**Cause:** Some datasets have numerical issues with mixed model fitting

**Expected:** Occasional warnings are normal (marked in results)

**Concern:** If >20% of iterations fail, review parameters

### Issue: Power estimates are unstable

**Cause:** Too few iterations (only 20)

**Solution:** Increase `n_iterations` to 100 or 1000

### Issue: Visualization script fails

**Cause:** Simulation results not found

**Solution:** Run Step 1 first, ensure RData file exists

---

## Performance Benchmarks

**System:** M1 Mac, 16GB RAM

| Task | Time |
|------|------|
| Sigma cache build | ~10 sec |
| Single Monte Carlo run | ~0.5 sec |
| 480 total simulations | ~5 min |
| Visualization generation | ~10 sec |
| **Total workflow** | **~5-6 min** |

**Scaling:**
- Runtime scales linearly with iterations
- 1000 iterations ≈ 25 minutes
- Consider parallel processing for large grids

---

## Next Steps

After completing the workflow:

1. **Examine results** - Check if patterns match expectations
2. **Document findings** - Note power differences between scenarios
3. **Create additional visualizations** - Use ggplot2 with simulation_results
4. **Prepare manuscript** - Export publication-ready figures
5. **Archive results** - Commit RData files with git LFS

---

## References

- **Main documentation:** `docs/TWO_SCENARIO_COMPARISON.md`
- **Implementation details:** `docs/PSEUDOCODE_COMPARISON_HENDRICKSON_VS_OURS.md`
- **Technical notes:** `docs/technical_differences_scaling_and_pd.pdf`
- **Original paper:** Hendrickson et al. (2020), Statistics in Medicine

---

*Last updated: 2025-11-18*
