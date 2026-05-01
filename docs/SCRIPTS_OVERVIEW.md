# R Scripts Overview

Summary of R scripts in the `pmsimstats2025/analysis/scripts/` directory.

## Setup & Infrastructure (00_*)

| Script | Scope |
|--------|-------|
| `00_database_setup.R` | Template for connecting to SQLite, PostgreSQL, MySQL, and ODBC databases |
| `00_setup_parallel.R` | Configures parallel processing (doParallel, future/furrr) for local, cluster, and cloud environments |
| `00_testing_guide.R` | Documentation/examples for testthat unit tests, integration tests, and data validation |
| `02_data_validation.R` | Data quality checks: missingness, types, outliers, duplicates, range validation (uses Palmer penguins as example) |

## Core Functions

| Script | Scope |
|--------|-------|
| `pm_functions.R` | Shared library: modified Gompertz growth curves, carryover decay models (exponential, linear, Weibull), bio-response calculations with biomarker×treatment interaction |

## Main Simulation Scripts

| Script | Scope |
|--------|-------|
| `simulation_clustered.R` | **Primary simulation**: Hybrid, OL+BDC, Crossover N-of-1 trials with clustered 8-point measurements; parallel processing; AR(1) correlation |
| `simulation_evenly_spaced.R` | OL, Crossover, Parallel designs with 4 evenly-spaced measurements; three-factor response model (BR, ER, TR) |
| `simulation_clustered_flexible_correlation.R` | Clustered designs with selectable correlation structures: AR1, exponential, power law, Matérn, rational quadratic |
| `simulation_evenly_spaced_flexible_correlation.R` | Evenly-spaced designs with the same flexible correlation structures |
| `simulation_clustered_exponential.R` | Clustered designs with exponential decay carryover model (half-life parameterization) |
| `simulation_clustered_twoweek.R` | Clustered designs variant (appears similar to flexible correlation) |
| `simulation_evenly_spaced_powerlaw.R` | Evenly-spaced designs with 8 timepoints using power law correlation |
| `simulationplus.R` | Simplified N-of-1 simulation with constant effect model for learning/prototyping |

## Sensitivity Analyses

| Script | Scope |
|--------|-------|
| `simulation_censoring.R` | Replicates Hendrickson Figure 4: power as function of dropout/censoring patterns (β₀ flat hazard, β₁ symptom-dependent) |
| `simulation_response_params.R` | Analogous to Hendrickson Figure 5: power sensitivity to response trajectory parameters (BR_rate, ER_rate, TR_rate) |

## Visualization

| Script | Scope |
|--------|-------|
| `plot_results.R` | Loads `simulation_clustered_results.RData`, generates power heatmaps by design, carryover, biomarker moderation |
| `plot_censoring.R` | Loads censoring simulation results, generates Hendrickson Figure 4A-style heatmaps |
| `visualize_heatmaps.R` | 4-panel heatmaps: carryover × biomarker correlation, faceted by design × carryover-model specification |
| `visualize_hendrickson_style.R` | Figure 4-equivalent: Panel A (no carryover, varying biomarker correlation), Panel B (varying carryover) |
| `visualize_power_heatmaps_4panel.R` | 4-panel power heatmaps: carryover vs effect size, comparing Hybrid vs Crossover with/without carryover model |

## Diagnostics & Debugging

| Script | Scope |
|--------|-------|
| `check_model_specification.R` | Verifies analysis model matches the data generating process; shows true coefficients |
| `check_model_specification2.R` | Part 2: fixes bm_centered scaling, verifies model recovery |
| `check_model_specification3.R` | Part 3: compares lmer (no AR1) vs lme (with AR1) |
| `debug_effect_recovery.R` | Checks if biomarker moderation effect is recovered unbiasedly across carryover conditions |
| `carryover_variance_analysis.R` | Explains how carryover reduces variance in `effective_drug_weeks`, thereby reducing power |
| `profile_simulation.R` | Performance profiling with profvis/bench to identify bottlenecks for potential Rcpp optimization |
| `example_parameter_validation.R` | Demonstrates pre-simulation validation of parameter combinations for positive definiteness |

## Positive Definiteness Testing

| Script | Scope |
|--------|-------|
| `test_pd_boundaries.R` | Systematically tests correlation parameters to find positive-definiteness boundaries for both designs |
| `test_new_params.R` | Quick test of autocorr=0.6, c.bm=0.6 |
| `test_boundary_autocorr_06.R` | Finds maximum c.bm with autocorr=0.6 |
| `test_cbm_06.R` | Tests c.bm=0, 0.3, 0.6 with autocorr=0.6 |
| `test_autocorr_for_cbm_06.R` | Finds minimum autocorrelation needed for c.bm=0.6 |

## Post-Simulation

| Script | Scope |
|--------|-------|
| `summarize_results.R` | Reads and summarizes `full_pmsim_analysis_hyb_versus_co.RData`; reports power statistics |
| `99_reproducibility_check.R` | Verifies reproducibility: R version, renv status, file integrity, data checksums |
