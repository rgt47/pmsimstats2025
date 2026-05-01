# Plan: Simplified Carryover-Focused Simulation

## Context

The main simulation (`simulation_clustered.R`, 1233 lines) sweeps 4
designs x 2 sample sizes x 4 biomarker moderation levels x 3 carryover
levels x 2 model specifications = ~144 conditions. This complexity
obscures the central finding: even 17-hour carryover half-lives destroy
power when unmodeled.

The goal is a focused, readable script that isolates this finding for
the Hybrid design, suitable for a reader who wants to understand the
carryover effect without navigating the full parameter grid.

## Key Discovery

`simulation_clustered.R` is self-contained -- it does NOT source
`pm_functions.R`. It defines its own sigma construction
(`build_sigma_guaranteed_pd`), data generation
(`generate_participant_twostage`), and design functions inline. The
simplified script must therefore carry its own function definitions
(extracted from simulation_clustered.R), not source a shared library.

## Structure

```
analysis/scripts/carryover_focus/
+-- PLAN.md                             # This file
+-- README.md                           # Purpose, how to run, what it shows
+-- carryover_power_simulation.R        # Single self-contained script
+-- output/                             # Generated results
```

## Script Design: `carryover_power_simulation.R`

One linear, self-contained R script (~350-450 lines). Organized in
clearly demarcated sections:

### Section 1: Setup (~20 lines)

- Packages: tidyverse, nlme, MASS, future, furrr, patchwork
- Output directory creation

### Section 2: Fixed Parameters (~40 lines)

- n_participants = 70 (fixed)
- n_iterations = 200 (default; env-var overridable for quick runs)
- biomarker_moderation = 0.4 (fixed -- moderate effect)
- biomarker_correlation = 0.3 (fixed)
- Response model parameters (BR/ER/TR rates, baseline, biomarker)
- Correlation parameters (c.br, c.er, c.tr, c.cf1t, c.cfct, etc.)

### Section 3: Helper Functions (~150 lines)

Extracted from simulation_clustered.R with minimal modification:

- `build_sigma_guaranteed_pd()` -- lines 208-330 of original
- `generate_participant_twostage()` -- lines 336-359 of original
- `create_hybrid_design()` -- lines 394-417 of original

### Section 4: Parameter Sweep (~15 lines)

Only two dimensions:

- carryover_halflife: c(0, 0.05, 0.1, 0.15, 0.2)
  (5 values -- finer than the original 3, to show the cliff)
- model_carryover: c(TRUE, FALSE)
- Remove redundant carryover=0 + model=FALSE row
- Total: 9 conditions (vs 144 in original)

### Section 5: Single Iteration Function (~120 lines)

Extracted from simulation_clustered.R `run_single_iteration()` but
simplified:

- Hardcoded to Hybrid design only (no design dispatch)
- Hardcoded to nlme::lme with corCAR1 (no engine selection)
- Same seed management: `set.seed(iter * 1000 + condition_id)`
- Same data generation, carryover calculation, model fitting logic
- Returns: iteration, carryover_halflife, model_carryover,
  effect_size, se, t_value, p_value, significant

### Section 6: Simulation Loop (~30 lines)

- Parallel execution via future_map_dfr (same pattern as original)
- Progress reporting per condition
- Bind results

### Section 7: Results Summary (~30 lines)

- Aggregate power, mean effect, SD, Wilson CI by condition
- Print formatted comparison table to console
- Compute power_difference and pct_change for modeled vs ignored

### Section 8: Visualization (~50 lines)

One figure with two panels (side by side):

- Left panel: Power vs carryover half-life, two lines
  (modeled=green, ignored=red)
- Right panel: Mean effect size vs carryover half-life, two lines
- Provenance annotation per CLAUDE.md convention
- Save as PDF and PNG to output/

### Section 9: Save Results (~10 lines)

- Save RData with results, summary, parameters
- Print completion message

## What Is Deliberately Omitted

- Three designs (Crossover, OL+BDC, OL) -- Hybrid only
- Sample size sweep (n=35) -- fixed at 70
- Biomarker moderation sweep (0, 0.2, 0.6) -- fixed at 0.4
- Type I error conditions (bm_mod=0)
- Sigma cache (only one sigma needed)
- Heatmap visualization (line plot is clearer for 1D sweep)
- Log file sink
- lmer alternative engine

## Source Lines Reference (simulation_clustered.R)

These are the exact line ranges to extract from the original:

- `build_sigma_guaranteed_pd()`: lines 208-330
- `generate_participant_twostage()`: lines 336-359
- `create_hybrid_design()`: lines 394-417
- `measurement_weeks_hybrid`: line 366
- `calculate_carryover` / `calculate_carryover_vec`: lines 104-118
- `run_single_iteration()` logic: lines 600-849
- Parameters: lines 28-81

## Files to Create

1. `analysis/scripts/carryover_focus/README.md`
2. `analysis/scripts/carryover_focus/carryover_power_simulation.R`
3. `analysis/scripts/carryover_focus/output/` (directory, via script)

## Verification

1. Run: `cd analysis/scripts/carryover_focus && Rscript carryover_power_simulation.R`
2. Confirm output/carryover_focus_results.RData is created
3. Confirm output/carryover_power_hybrid.pdf is created
4. Confirm console output shows the power cliff pattern:
   - carryover=0: modeled power ~0.80
   - carryover=0.1, modeled: ~0.85
   - carryover=0.1, ignored: ~0.30
5. Compare key values against simulation_clustered_results.RData
   to confirm consistency

## Implementation Progress

- [x] Read simulation_clustered.R (lines 1-849) -- functions extracted
- [ ] Create carryover_power_simulation.R
- [ ] Create README.md
- [ ] Run and verify
