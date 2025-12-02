# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**pmsimstats2025** is an N-of-1 clinical trial simulation comparing multiple trial designs (Hybrid, Crossover, Open-Label, Parallel) with a focus on statistical power when biomarker×treatment interactions exist and carryover effects are present.

The work is aligned with **Hendrickson et al. (2020)** methodology for N-of-1 trials, with extensions for explicit biomarker moderation and carryover effect modeling.

---

## Quick Start Commands

### Run Simulations

```bash
# Clustered measurement schedule designs (Hybrid, OL+BDC)
Rscript analysis/scripts/simulation_clustered.R

# Evenly-spaced designs (OL, Crossover, Parallel)
Rscript analysis/scripts/simulation_evenly_spaced.R

# All designs combined
Rscript analysis/scripts/simulationplus.R
```

### Generate Visualizations

```bash
# Power heatmaps from simulation results
Rscript analysis/scripts/visualize_hendrickson_style.R
Rscript analysis/scripts/visualize_power_heatmaps_4panel.R
```

### Validate Package Dependencies

```bash
# Full validation with auto-fix (recommended before commits)
make check-renv

# Validation only (no auto-install)
make check-renv-no-fix
```

### Docker Workflow (No Local R Required)

```bash
# Start container with auto-detected profile
make r

# Or explicitly
make docker-run
```

---

## Architecture & Key Files

### Core Simulation Architecture

The simulation consists of a two-layer architecture:

#### Layer 1: Parameter Definition & Data Generation (`pm_functions.R`)

- **`build_sigma_matrix()`** - Constructs multivariate normal covariance matrices with PD validation
  - Uses fixed Hendrickson correlation values: `c.tv=0.8, c.pb=0.8, c.br=0.8, c.cf1t=0.2, c.cfct=0.1`
  - Clamps biomarker correlations to `[-0.99, 0.99]` to maintain positive definiteness
  - Returns `NULL` for invalid parameter combinations (rejects rather than auto-fixes)

- **`generate_data()`** - Generates multivariate normal trial data for one participant
  - Three factors: time-variant (TV), pharmacologic biomarker (PB), biological response (BR)
  - Each measured at multiple timepoints based on trial design

- **`calculate_bio_response_with_interaction()`** - Computes BR means with treatment & biomarker effects
  - BR varies by treatment status (on/off) and biomarker level
  - Implements carryover as exponential decay: `(1/2)^(weeks_off/t1half)`

- **`validate_correlation_structure()`** - Standalone validation utility
  - Checks eigenvalue range and condition number
  - Warns if condition number > 100 (ill-conditioned matrices)

#### Layer 2: Monte Carlo Simulation (`simulation_clustered.R`, `simulation_evenly_spaced.R`, `simulationplus.R`)

Three separate simulation programs, each with identical structure:

1. **Sigma Matrix Cache Building** - Pre-validates all parameter combinations
   - Tests each design × biomarker_correlation × carryover level combination
   - Prints diagnostic output with eigenvalue ranges
   - Only proceeds if all combinations produce valid (PD) matrices

2. **Monte Carlo Loop** - Runs multiple iterations per condition
   - Generates new data each iteration
   - Fits mixed-effects models with and without carryover
   - Tests `treatment × biomarker` interaction (α = 0.05)
   - Tracks convergence failures and effect sizes

3. **Power Summarization** - Aggregates results
   - Power = proportion of iterations with p < 0.05
   - Averages effect sizes and standard errors
   - Saves both individual results and summary tables

### Design-Specific Measurement Schedules

| Design | Schedule | Points | Structure | File |
|--------|----------|--------|-----------|------|
| Hybrid | `c(4, 8, 9, 10, 11, 12, 16, 20)` | 8 | Dense cluster at transition | `simulation_clustered.R` |
| OL+BDC | `c(4, 8, 12, 16, 17, 18, 19, 20)` | 8 | Dense cluster at discontinuation | `simulation_clustered.R` |
| OL | `c(2, 8, 14, 20)` | 4 | Evenly spaced (6-week gaps) | `simulation_evenly_spaced.R` |
| Crossover | `c(2, 8, 14, 20)` | 4 | Evenly spaced (6-week gaps) | `simulation_evenly_spaced.R` |
| Parallel | `c(2, 8, 14, 20)` | 4 | Evenly spaced (6-week gaps) | `simulation_evenly_spaced.R` |

**Key Design Principle**: Clustered measurement schedules (small gaps) enable 8 measurement points within a 20-week trial. Evenly-spaced designs with many points cause matrix ill-conditioning, so they use 4 points with 6-week gaps.

### Output Files

**Simulation Results** (RData format, loaded in R):
- `analysis/output/simulation_clustered_results.RData`
- `analysis/output/simulation_evenly_spaced_results.RData`
- `analysis/output/simulationplus_results.RData`

Each contains:
- `simulation_results` - 1 row per iteration with p-values, effect sizes, convergence status
- `simulation_summary` - Aggregated power by condition
- `param_grid`, `model_params`, `resp_param`, `baseline_param` - Parameter specifications

**Visualizations**:
- `analysis/output/power_heatmap_clustered.pdf` - Power heatmaps for Hybrid and OL+BDC
- `analysis/output/power_heatmap_evenly_spaced.pdf` - Power heatmaps for OL, Crossover, Parallel
- `analysis/output/power_heatmap_4panel.pdf` - Combined 4-panel heatmap

**Logs**:
- `analysis/output/simulation_clustered_log.txt`
- `analysis/output/simulation_evenly_spaced_log.txt`
- Timestamps and diagnostic output from matrix validation and analysis

---

## Critical Implementation Details

### 1. Positive Definiteness Handling

**Problem**: With certain parameter combinations (especially high biomarker correlations or many measurement points), covariance matrices become non-positive definite.

**Current Approach**: Rejects invalid combinations rather than auto-fixing
- `build_sigma_matrix()` returns `NULL` for non-PD matrices
- Simulation skips combinations that fail validation
- Diagnostic eigenvalue output helps identify problematic parameters

**Key Fix**: Biomarker correlation clamping to `[-0.99, 0.99]`
```r
# In pm_functions.R, around line 426
scaled_correlation <- pmax(-0.99, pmin(0.99, scaled_correlation))
```

### 2. Correlation Structure (Fixed Values)

Following Hendrickson et al. (2020), correlations are **fixed and do not depend on carryover parameters**:
- `c.tv = 0.8` - Time-variant autocorrelation (AR(1))
- `c.pb = 0.8` - Pharmacologic biomarker autocorrelation
- `c.br = 0.8` - Biological response autocorrelation
- `c.cf1t = 0.2` - Cross-correlation at same time point
- `c.cfct = 0.1` - Cross-correlation at different time points

**Principle**: Carryover affects MEANS only, not CORRELATIONS.

### 3. Carryover Modeling (Novel Extension)

Two analysis scenarios to quantify the cost of ignoring carryover:

**Scenario 1: WITH Carryover** - Model includes explicit carryover term
```r
lmer(response ~ treatment * biomarker + week + carryover_effect + (1 | participant_id))
```

**Scenario 2: WITHOUT Carryover** - Replicates Hendrickson's original approach
```r
lmer(response ~ treatment * biomarker + week + (1 | participant_id))
```

Carryover effect computed as exponential decay:
```r
carryover_effect = scale_factor * (0.5)^(weeks_off / t1half)
```
where `scale_factor = 1.0` (BR-only carryover) and `t1half` varies across conditions.

### 4. Measurement Schedule and Matrix Conditioning

Different designs require different measurement schedules due to **eigenvalue sensitivity**:

- **Clustered schedules** (Hybrid, OL+BDC): Measurement points cluster around transition times with sparse baseline/follow-up measurements. This reduces eigenvalue ratios despite 8 total points.

- **Evenly-spaced schedules** (OL, Crossover, Parallel): Must use only 4 points (6-week gaps) to maintain PD. Using 8 evenly-spaced points causes matrix ill-conditioning.

This is documented in `docs/technical_differences_scaling_and_pd.pdf`.

---

## Modifying Simulation Parameters

### Biomarker Correlation Grid

Located in simulation files (e.g., `simulation_clustered.R` lines ~500-520):

```r
param_grid <- expand_grid(
  biomarker_moderation = c(0.2, 0.4),        # Change here
  carryover_t1half = c(0, 1.0, 2.0),
  n_participants = 70
)
```

**Caution**: Higher correlations (>0.4) may cause PD failures. Check sigma cache output.

### Carryover Gradient

Edit the same section:
```r
carryover_t1half = c(0, 0.5, 1.0, 1.5, 2.0, 2.0, 3.0)  # Add values here
```

### Iteration Count

Around line 630 in simulation files:
```r
n_iterations = 20  # Change to 100 for more stable estimates (slower)
```

### Fixed Correlation Parameters

In `pm_functions.R` (lines ~60-66):
```r
c.tv = 0.8, c.pb = 0.8, c.br = 0.8, c.cf1t = 0.2, c.cfct = 0.1
```

**Warning**: Changing these may affect PD. Always check sigma cache output after changes.

---

## Documentation Structure

### Method & Theory
- **`analysis/scripts/README.md`** - Comprehensive workflow guide with parameter details
- **`docs/simulation_white_paper.md`** - Methodology alignment with Hendrickson, pseudocode algorithms
- **`docs/chat.Rmd`** - White paper comparing crossover vs. N-of-1 trial designs (renders to PDF)

### Technical
- **`docs/technical_differences_scaling_and_pd.pdf`** - Matrix scaling and PD handling rationale
- **`docs/correlation_structure_design.pdf`** - Covariance matrix construction guide
- **`analysis/scripts/CORRELATION_ALIGNMENT_CHANGES.md`** - Before/after correlation changes

### Theory
- **`analysis/scripts/carryover_correlation_theory.tex`** - Mathematical foundations

---

## Troubleshooting Common Issues

### "Non-positive definite matrix rejected"

**Normal behavior**: Invalid parameter combinations are skipped with diagnostic output.

**Example output**:
```
REJECTED: Non-positive definite matrix detected.
Eigenvalue range: 0.000234 to 3.421 (ratio: 14621)
Condition number: 14621 (ill-conditioned!)
Matrix size: 32×32
```

**Action**: Check which condition failed and either reduce correlation values or adjust measurement schedule.

### Mixed model convergence warnings

**Expected**: Occasional warnings on random data are normal.

**Monitor**: If >20% of iterations fail to converge, check the `error` column in `simulation_results` and review parameter choices.

### High eigenvalue ratios (>100 condition number)

**Cause**: Measurement schedules with many points (especially evenly-spaced) or high correlations.

**Solution**: Use clustered measurement schedules for more points, or reduce number of timepoints.

---

## Performance Notes

**Typical runtime** (M1 Mac, 16GB RAM, 20 iterations, 12-18 conditions):
- Sigma cache build: ~5-10 seconds
- Monte Carlo simulation: ~5-15 minutes
- Visualization: ~10 seconds
- **Total**: ~5-20 minutes

**Scaling**: Linear with iteration count
- 100 iterations ≈ 25 minutes
- 1000 iterations ≈ 4 hours

---

## Alignment with Hendrickson et al. (2020)

| Feature | Status | Implementation |
|---------|--------|-----------------|
| 4-path randomization | ✅ Aligned | Hybrid design in `simulation_clustered.R` |
| BR-only carryover | ✅ Aligned | `scale_factor = 1.0` in carryover computation |
| Fixed correlations | ✅ Aligned | `c.tv=0.8, c.pb=0.8, c.br=0.8, c.cf1t=0.2, c.cfct=0.1` |
| Time effect in model | ✅ Aligned | `week` term in all lmer formulas |
| Random intercept | ✅ Aligned | `(1 \| participant_id)` |
| **Carryover in model** | ✅ **Enhancement** | Novel: explicit `carryover_effect` term |
| **Biomarker moderation** | ✅ **Enhancement** | Novel: systematic variation of `biomarker_correlation` |

---

## References

**Primary**:
- Hendrickson, E., et al. (2020). N-of-1 trials with multiple randomization structures for individualized treatment. *Statistics in Medicine*, 39(25), 3581-3599.

**Crossover & Design**:
- Dwan, K., et al. (2019). CONSORT extension for reporting randomized controlled trials with non-pharmacological interventions: Guide and checklist. *BMJ*.
- Mills, E. J., et al. (2009). Design, analysis, and presentation of crossover trials. *Trials*, 10, 27.

**N-of-1 Methodology**:
- Zucker, D. R., et al. (2010). Evaluation and management of comorbid psychiatric disorders in chronic pain patients. *Current Pain & Headache Reports*, 14(1), 33-40.
- Lillie, E. O., et al. (2011). The n-of-1 clinical trial: The ultimate strategy for individualizing medicine? *Personalized Medicine*, 8(2), 161-173.
- Punja, S., et al. (2016). n-of-1 trials. *Journal of the American Medical Association*, 316(23), 2459.

---

## Git Repository

- **URL**: https://github.com/rgt47/pmsimstats2025
- **Status**: Public with collaborator access
- **Branch**: main
- **Last Updated**: 2025-11-26

---

## Key Contacts & Collaboration

- **Primary Collaborator**: rchendrickson (write access)
- **Package/Environment Manager**: renv (locked dependencies in `renv.lock`)
- **Docker Profile**: analysis (RStudio Server)

---

## Common Development Tasks

### Add a New Design Variant

1. Define measurement schedule in simulation file (around line 100)
2. Add to `param_grid` expansion (lines ~500-520)
3. Update heatmap labels (lines ~600-650)
4. Run simulation and check sigma cache output

### Test a New Biomarker Correlation

1. Edit `param_grid` in simulation file to include new value
2. Run sigma cache build first: check diagnostic output
3. If PD failure, examine eigenvalue diagnostics and adjust measurement schedule
4. Verify power heatmap makes logical sense

### Render White Paper to PDF

```r
# In R console
rmarkdown::render("docs/chat.Rmd")
```

Produces `/docs/chat.pdf` with:
- Table of contents
- LaTeX equations
- Code output and figures

### Validate Dependencies Before Commit

```bash
make check-renv  # Full validation with auto-fix
```

This checks that all packages used in code match `DESCRIPTION` and `renv.lock`, preventing dependency mismatches when others clone the repo.

---

*Last reviewed: 2025-12-02*
