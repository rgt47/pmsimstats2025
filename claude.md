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
  - Uses adjusted correlation values optimized for PD: `c.br=0.75, c.er=0.75, c.tr=0.75, c.cf1t=0.12, c.cfct=0.05`
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

- **`validate_parameter_grid()`** - Comprehensive pre-simulation validation function (NEW)
  - Tests all parameter combinations for positive definiteness
  - Computes condition numbers (κ) for numerical stability assessment
  - Returns structured list: `valid_combinations`, `invalid_combinations`, `n_valid`, `n_invalid`, `condition_numbers`, `invalid_reasons`
  - Enables filtering of problematic combinations before expensive Monte Carlo runs

- **`report_parameter_validation()`** - Generates human-readable validation summary (NEW)
  - Prints diagnostic output with PD status and condition number statistics
  - Reports problematic combinations and failure reasons
  - Provides recommendations for parameter adjustment

#### Layer 2: Monte Carlo Simulation (`simulation_clustered.R`, `simulation_evenly_spaced.R`, `simulationplus.R`)

Three separate simulation programs, each with identical structure:

1. **Pre-Simulation Sigma Validation** - Validates covariance matrices before simulation (NEW)
   - Tests 3 unique sigma structures corresponding to each design
   - Checks positive definiteness (eigenvalue criterion)
   - Computes condition numbers for numerical stability
   - Stops with diagnostic error if any sigma fails validation
   - Example output shows condition numbers (κ = 56-61) for well-conditioned matrices

2. **Sigma Matrix Cache Building** - Builds lookup table for data generation
   - Pre-computes sigma matrices for each parameter combination
   - Stores in-memory for fast repeated access during Monte Carlo iterations
   - Only proceeds if pre-simulation validation passed

3. **Monte Carlo Loop** - Runs multiple iterations per condition
   - Generates new data each iteration using pre-computed sigma matrices
   - Fits mixed-effects models with and without carryover
   - Tests `treatment × biomarker` interaction (α = 0.05)
   - Tracks convergence failures and effect sizes

4. **Power Summarization** - Aggregates results
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

### 1. Positive Definiteness Handling (Two-Level Validation)

**Problem**: Covariance matrices become non-positive definite with certain parameter combinations, especially high biomarker correlations or many measurement points.

**Current Approach**: Multi-level validation ensures PD before expensive Monte Carlo runs

**Level 1: Pre-Simulation Validation** (NEW)
- Tests 3 unique sigma structures before simulation starts
- `validate_parameter_grid()` checks all parameter combinations
- Stops with diagnostic error if any sigma fails PD test
- Avoids wasting computational resources on invalid combinations

**Level 2: Build-Time Validation**
- `build_sigma_matrix()` returns `NULL` for non-PD matrices during simulation
- Simulation logic rejects invalid combinations
- Diagnostic eigenvalue output helps identify problematic parameters

**Level 3: Numerical Stability**
- Condition number (κ) computed for all valid matrices
- Warnings issued if κ > 100 (ill-conditioned)
- Allows users to assess numerical stability before analysis

**Parameter Tuning**: Balanced reduction of Hendrickson correlation values
- Biomarker correlation clamping to `[-0.99, 0.99]`
```r
# In pm_functions.R, around line 426
scaled_correlation <- pmax(-0.99, pmin(0.99, scaled_correlation))
```
- Ensures all parameter combinations produce well-conditioned matrices (κ < 100)

### 2. Correlation Structure (Fixed Values)

Correlations are **fixed and do not depend on carryover parameters**. Values have been optimized to maintain positive definiteness across all parameter combinations:

- `c.tv = 0.8` - Time-variant autocorrelation (AR(1))
- `c.pb = 0.8` - Pharmacologic biomarker autocorrelation
- `c.br = 0.75` - Biological response autocorrelation (reduced from 0.8)
- `c.tr = 0.75` - Treatment response autocorrelation (reduced from 0.8)
- `c.er = 0.75` - Error/residual autocorrelation (reduced from 0.8)
- `c.cf1t = 0.12` - Cross-correlation at same time point (reduced from 0.2)
- `c.cfct = 0.05` - Cross-correlation at different time points (reduced from 0.1)
- `c.bm_baseline = 0.25` - Biomarker-baseline correlation (reduced from 0.3)
- `c.baseline_resp = 0.3` - Baseline-response correlation (reduced from 0.4)

**Principle**: Carryover affects MEANS only, not CORRELATIONS.

**Adjustment Rationale**: Original Hendrickson values created non-positive definite (non-PD) covariance matrices when combined with clustered measurement schedules (8 timepoints) and varying biomarker correlations. Balanced reduction across autocorrelation and cross-correlation parameters maintains correlation hierarchy (c.cfct < c.cf1t < c.autocorr) while ensuring PD for all parameter combinations. Condition numbers remain in well-conditioned range (κ < 100).

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

In `pm_functions.R` (lines ~38-46):
```r
c.tv = 0.8, c.pb = 0.8, c.br = 0.75, c.tr = 0.75, c.er = 0.75,
c.cf1t = 0.12, c.cfct = 0.05,
c.bm_baseline = 0.25, c.baseline_resp = 0.3
```

**Note**: These values have been optimized from Hendrickson originals to maintain positive definiteness across all parameter combinations while preserving correlation hierarchy. Original values: `c.br=0.8, c.cf1t=0.2, c.cfct=0.1, c.bm_baseline=0.3, c.baseline_resp=0.4`.

**Warning**: Changing these may affect PD. Always run pre-simulation validation via `validate_parameter_grid()` after changes.

---

## Documentation Structure

### Method & Theory
- **`analysis/scripts/README.md`** - Comprehensive workflow guide with parameter details
- **`docs/simulation_white_paper.md`** - Methodology alignment with Hendrickson, pseudocode algorithms
- **`docs/chat.Rmd`** - White paper comparing crossover vs. N-of-1 trial designs (renders to PDF)

### Mathematical Foundations (NEW - December 2025)
- **`docs/sigma_matrix_derivation.tex`** - Comprehensive mathematical derivation of Σ = D·R·D identity with block partitioning strategy, eigenvalue properties, and two-stage sampling algorithm
- **`docs/positive_definiteness_constraints.tex`** - Mathematical derivation of PD constraints using Sylvester's criterion, eigenvalue analysis, Gershgorin circles, and empirical guidelines for correlation parameters
- **`docs/biomarker_interaction_mechanism.tex`** - Detailed explanation of two-level biomarker-treatment interaction through covariance structure (c.bm) and mean-level modulation (biomarker_moderation)

### Parameter Validation (NEW - December 2025)
- **`analysis/scripts/PARAMETER_VALIDATION_GUIDE.md`** - Comprehensive user guide for parameter validation including function signatures, workflow examples, condition number interpretation, and best practices
- **`analysis/scripts/VALIDATION_OUTPUT_EXAMPLES.md`** - Detailed examples showing validation function outputs for three scenarios: all valid, some invalid, and ill-conditioned matrices
- **`analysis/scripts/example_parameter_validation.R`** - Runnable example demonstrating parameter validation workflow from parameter definition through simulation-ready grid

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
| Fixed correlations | ⚠️ Adapted | Hendrickson structure maintained; values optimized for PD: `c.br=0.75, c.cf1t=0.12, c.cfct=0.05` (see rationale below) |
| Time effect in model | ✅ Aligned | `week` term in all lmer formulas |
| Random intercept | ✅ Aligned | `(1 \| participant_id)` |
| **Carryover in model** | ✅ **Enhancement** | Novel: explicit `carryover_effect` term |
| **Biomarker moderation** | ✅ **Enhancement** | Novel: systematic variation of `biomarker_correlation` |
| **Parameter validation** | ✅ **Enhancement** | Novel: pre-simulation sigma validation via `validate_parameter_grid()` |

**Correlation Value Adjustment Rationale**: Original Hendrickson values were specified for different measurement schedules and did not account for the interaction between clustered measurement schedules (8 timepoints), block-partitioned covariance structure (26×26 decomposed), and varying biomarker correlations. Balanced reduction of autocorrelation and cross-correlation parameters maintains the theoretical structure while ensuring all parameter combinations produce positive definite matrices with good numerical conditioning (κ < 100).

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

---

## Recent Changes (December 2025)

### Pre-Simulation Sigma Validation Integration
- Added `validate_parameter_grid()` function to test all parameter combinations before simulation
- Integrated validation into both `simulation_clustered.R` and `simulation_evenly_spaced.R`
- Simulations now stop with diagnostic error if any sigma matrix fails PD validation
- Prevents wasting computational resources on invalid parameter combinations

### Mathematical Documentation
- Created `docs/sigma_matrix_derivation.tex` - Complete derivation of Σ = D·R·D with block partitioning strategy
- Created `docs/positive_definiteness_constraints.tex` - Mathematical analysis of PD constraints and correlation hierarchy
- Created `docs/biomarker_interaction_mechanism.tex` - Explanation of two-level biomarker-treatment interaction

### Parameter Validation Documentation
- Created `analysis/scripts/PARAMETER_VALIDATION_GUIDE.md` - User guide for validation workflow
- Created `analysis/scripts/VALIDATION_OUTPUT_EXAMPLES.md` - Detailed output examples and interpretation
- Created `analysis/scripts/example_parameter_validation.R` - Runnable validation example

### Correlation Parameter Optimization
- Adjusted correlation values from Hendrickson originals to optimize positive definiteness
- All parameter combinations now pass validation with well-conditioned matrices (κ < 100)
- Maintains correlation hierarchy and theoretical structure

### Code Quality Improvements
- Fixed string operator issues in `pm_functions.R` (string multiplication to `strrep()`)
- Simplified validation approach: tests 3 unique sigma structures before simulation
- Robust error handling: stops simulation immediately on validation failure

*Last reviewed: 2025-12-10*
