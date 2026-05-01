# pmsimstats2025

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/rgt47/pmsimstats2025/actions/workflows/r-package.yml/badge.svg)](https://github.com/rgt47/pmsimstats2025/actions/workflows/r-package.yml)
<!-- badges: end -->

**pmsimstats2025** is an N-of-1 clinical trial simulation study comparing
multiple trial designs (Hybrid, Crossover, Open-Label, Parallel) with a focus
on statistical power when biomarker-treatment interactions exist and carryover
effects are present.

The work is aligned with **Hendrickson et al. (2020)** methodology for N-of-1
trials, with extensions for explicit biomarker moderation and carryover effect
modeling.

## Features

- **Multiple Trial Designs**: Hybrid, Open-Label + Biomarker Discontinuation
  Confirmation (OL+BDC), Open-Label, Crossover, and Parallel designs
- **Biomarker-Treatment Interactions**: Systematic variation of biomarker
  correlation parameters to model treatment response heterogeneity
- **Carryover Effect Modeling**: Exponential decay carryover with configurable
  half-life parameters
- **Positive Definiteness Validation**: Pre-simulation sigma matrix validation
  with condition number monitoring
- **Mixed-Effects Analysis**: `lme` models (nlme) with continuous-time AR(1)
  correlation (`corCAR1`) and with/without explicit carryover terms

## Installation

Clone the repository and restore the R environment:

```bash
git clone https://github.com/rgt47/pmsimstats2025.git
cd pmsimstats2025
make r  # Enter Docker container with all dependencies
```

Or restore dependencies with renv:

```r
renv::restore()
```

## Quick Start

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

### Docker Workflow

```bash
# Start container with auto-detected profile
make r

# Or explicitly
make docker-run
```

## Trial Designs

| Design | Schedule | Points | Structure |
|:-------|:---------|:-------|:----------|
| Hybrid | c(4, 8, 9, 10, 11, 12, 16, 20) | 8 | Dense cluster at transition |
| OL+BDC | c(4, 8, 12, 16, 17, 18, 19, 20) | 8 | Dense cluster at discontinuation |
| OL | c(2, 8, 14, 20) | 4 | Evenly spaced (6-week gaps) |
| Crossover | c(2, 8, 14, 20) | 4 | Evenly spaced (6-week gaps) |
| Parallel | c(2, 8, 14, 20) | 4 | Evenly spaced (6-week gaps) |

## Core Architecture

### Layer 1: Parameter Definition (`pm_functions.R`)

- `build_sigma_matrix()` - Constructs multivariate normal covariance matrices
  with positive definiteness validation
- `generate_data()` - Generates multivariate normal trial data for one
  participant
- `calculate_bio_response_with_interaction()` - Computes biological response
  means with treatment and biomarker effects
- `validate_parameter_grid()` - Pre-simulation validation of all parameter
  combinations

### Layer 2: Monte Carlo Simulation

Three simulation programs with identical structure:

1. Pre-simulation sigma validation
2. Sigma matrix cache building
3. Monte Carlo loop with mixed-effects models
4. Power summarization

## Output Files

### Simulation Results

- `analysis/output/simulation_clustered_results.RData`
- `analysis/output/simulation_evenly_spaced_results.RData`
- `analysis/output/simulationplus_results.RData`

### Visualizations

- `analysis/output/power_heatmap_clustered.pdf`
- `analysis/output/power_heatmap_evenly_spaced.pdf`
- `analysis/output/power_heatmap_4panel.pdf`

## Documentation

### Method and Theory

- `analysis/scripts/README.md` - Workflow guide with parameter details
- `docs/simulation_white_paper.md` - Methodology alignment with Hendrickson

### Mathematical Foundations

- `docs/sigma_matrix_derivation.tex` - Mathematical derivation of covariance
  structure
- `docs/positive_definiteness_constraints.tex` - PD constraints and eigenvalue
  analysis
- `docs/biomarker_interaction_mechanism.tex` - Two-level biomarker-treatment
  interaction

### Parameter Validation

- `analysis/scripts/PARAMETER_VALIDATION_GUIDE.md` - User guide for validation
  workflow
- `analysis/scripts/VALIDATION_OUTPUT_EXAMPLES.md` - Output interpretation

## Carryover Modeling

Two analysis scenarios to quantify the cost of ignoring carryover:

**Scenario 1: WITH Carryover**

```r
lme(response ~ effective_drug_weeks * bm_centered + week,
    random = ~ 1 | participant_id,
    correlation = corCAR1(form = ~ week | participant_id),
    data = trial_data)
```

**Scenario 2: WITHOUT Carryover**

```r
lme(response ~ treatment * bm_centered + week,
    random = ~ 1 | participant_id,
    correlation = corCAR1(form = ~ week | participant_id),
    data = trial_data)
```

Effective drug weeks computed with exponential decay carryover:

```r
effective_drug_weeks <- weeks_on_drug * (0.5)^(time_since_discontinuation / t_half)
```

## Type I Error Control

Under null conditions (biomarker_moderation = 0), empirical Type I error rates
range from 3-9% across designs and carryover levels, consistent with the
nominal 5% level (100 iterations, 95% CI ≈ 0.01-0.11).

## Alignment with Hendrickson et al. (2020)

| Feature | Status | Implementation |
|:--------|:-------|:---------------|
| 4-path randomization | Aligned | Hybrid design |
| BR-only carryover | Aligned | scale_factor = 1.0 |
| Fixed correlations | Adapted | Reduced for PD stability |
| Time effect in model | Aligned | week term in lme |
| Random intercept | Aligned | (1 \| participant_id) |
| **AR(1) correlation** | **Enhancement** | corCAR1 in nlme |
| **Carryover in model** | **Enhancement** | Explicit carryover term |
| **Biomarker moderation** | **Enhancement** | Systematic variation |
| **Type I error eval** | **Enhancement** | Null condition included |
| **Parameter validation** | **Enhancement** | Pre-simulation validation |

## Dependencies

Key packages:

- `nlme` - Mixed-effects modeling with AR(1) correlation (corCAR1)
- `lme4`, `lmerTest` - Alternative mixed-effects modeling
- `dplyr`, `tidyverse` - Data manipulation
- `ggplot2`, `patchwork` - Visualization
- `furrr`, `future` - Parallel processing
- `MASS`, `corpcor` - Matrix operations

## Performance Notes

Typical runtime (M1 Mac, 16GB RAM, 20 iterations):

- Sigma cache build: 5-10 seconds
- Monte Carlo simulation: 5-15 minutes
- Visualization: 10 seconds
- **Total**: 5-20 minutes

Scaling: Linear with iteration count (100 iterations = 25 minutes)

## Reproducibility

This project is developed using the zzcollab framework for reproducible
research:

```bash
git clone https://github.com/rgt47/pmsimstats2025.git
cd pmsimstats2025
make r  # Enter Docker container with all dependencies
```

## References


## License

GPL-3

## Authors

- Rebecca C. Hendrickson (rebecca.hendrickson@va.gov)
- Ronald G. Thomas (rgthomas@ucsd.edu)
