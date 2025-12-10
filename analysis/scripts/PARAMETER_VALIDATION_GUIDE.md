# Parameter Validation Guide

## Overview

Pre-simulation parameter validation checks all parameter combinations for **positive definiteness** before running expensive Monte Carlo simulations. This prevents wasted computation on invalid parameter sets.

## Key Functions

### 1. `validate_parameter_grid()`

Validates all parameter combinations and returns a detailed report.

**Usage:**
```r
validation_result <- validate_parameter_grid(
  param_grid = param_grid,           # Your parameter combinations
  trial_design = trial_design,       # Trial design tibble
  model_params = model_params,       # Fixed correlation parameters
  resp_param = resp_param,           # Response parameters
  baseline_param = baseline_param,   # Baseline parameters
  verbose = TRUE                     # Print detailed output
)
```

**Returns a list with:**
- `valid_combinations`: Tibble of parameter combinations that pass validation
- `invalid_combinations`: Tibble of combinations that fail
- `n_valid`: Number of valid combinations
- `n_invalid`: Number of invalid combinations
- `condition_numbers`: Condition numbers (κ) for valid combinations
- `invalid_reasons`: Character vector explaining failures

### 2. `report_parameter_validation()`

Prints a detailed validation report with statistics and recommendations.

**Usage:**
```r
report_parameter_validation(validation_result, param_grid)
```

---

## Workflow: How to Use in Your Simulations

### Step 1: Define Parameters

```r
library(tidyverse)
source("pm_functions.R")

# Create trial design
trial_design <- expand_grid(
  participant_id = 1:70,
  week = c(4, 8, 9, 10, 11, 12, 16, 20)
)

# Fixed parameters (from Hendrickson et al.)
model_params <- list(
  c.tv = 0.8, c.pb = 0.8, c.br = 0.8,
  c.cf1t = 0.2, c.cfct = 0.1,
  c.bm = 0.3  # Will be overridden by param_grid
)

# Response and baseline parameters
resp_param <- tibble(...)
baseline_param <- tibble(...)
```

### Step 2: Create Parameter Grid

```r
param_grid <- expand_grid(
  n_participants = 70,
  biomarker_correlation = c(0.0, 0.2, 0.4, 0.6),  # Test range
  carryover_t1half = c(0, 1.0, 2.0)
)

# Result: 12 combinations (4 × 3)
```

### Step 3: Validate Before Simulation

```r
validation_result <- validate_parameter_grid(
  param_grid, trial_design, model_params, resp_param, baseline_param
)

# Output:
# ================================================================================
# PRE-SIMULATION PARAMETER VALIDATION
# ================================================================================
#
# Total combinations tested: 12
# ✓ Valid combinations:       10 (83.3%)
# ✗ Invalid combinations:     2 (16.7%)
#
# PROBLEMATIC COMBINATIONS:
# ────────────────────────────────────────────────────────────────────────────────
#   Row 11: Non-positive definite matrix (c.bm=0.60)
#   Row 12: Non-positive definite matrix (c.bm=0.60)
# ...
```

### Step 4: Filter to Valid Combinations (Optional)

```r
if (validation_result$n_invalid > 0) {
  # Use only valid combinations
  param_grid <- validation_result$valid_combinations

  cat(sprintf("Proceeding with %d valid combinations\n", nrow(param_grid)))
}
```

### Step 5: Run Simulation

```r
# Now run your simulation with the validated param_grid
# It will only include combinations that are guaranteed to be PD
source("simulation_clustered.R")
```

---

## Example Output

### Successful Validation (All Valid)

```
================================================================================
PRE-SIMULATION PARAMETER VALIDATION
================================================================================

Total combinations tested: 12
✓ Valid combinations:       12 (100.0%)
✗ Invalid combinations:     0 (0.0%)

Condition number (κ) statistics for valid combinations:
  Mean:   28.5
  Median: 27.3
  Min:    15.2 (best conditioned)
  Max:    42.8 (worst conditioned)

================================================================================
✓ ALL PARAMETER COMBINATIONS ARE VALID
  Ready to proceed with simulation!
================================================================================
```

### Failed Validation (Some Invalid)

```
================================================================================
PRE-SIMULATION PARAMETER VALIDATION
================================================================================

Total combinations tested: 12
✓ Valid combinations:       10 (83.3%)
✗ Invalid combinations:     2 (16.7%)

PROBLEMATIC COMBINATIONS:
────────────────────────────────────────────────────────────────────────────────
  Row 11: Non-positive definite matrix (c.bm=0.60)
  Row 12: Non-positive definite matrix (c.bm=0.60)

Details of invalid combinations:
# A tibble: 2 × 3
  n_participants biomarker_correlation carryover_t1half
           <dbl>                  <dbl>             <dbl>
1             70                    0.6                0
2             70                    0.6                2

RECOMMENDATIONS:
────────────────────────────────────────────────────────────────────────────────
  • All failures are due to biomarker correlation (c.bm) values being too high
  • Consider reducing the maximum biomarker_correlation in param_grid
  • Current constraint: c.bm ≤ 0.6
  • You may need to reduce to: c.bm ≤ 0.4 or lower

================================================================================
⚠ SOME PARAMETER COMBINATIONS ARE INVALID
  2 combinations excluded from simulation
================================================================================
```

---

## Interpreting Condition Numbers

The condition number (κ) measures numerical stability:

- **κ < 10**: **Excellent** - Well-conditioned, no stability concerns
- **κ 10-100**: **Good** - Acceptable numerical stability
- **κ 100-1000**: **Poor** - Ill-conditioned, monitor for numerical errors
- **κ > 1000**: **Very Poor** - High risk of numerical instability

**Example:**
```
Condition number (κ) statistics for valid combinations:
  Mean:   45.2
  Median: 42.8
  Min:    15.2 (best conditioned)
  Max:    98.7 (worst conditioned)
```

**Interpretation:** Most combinations are well-conditioned (κ < 100), but the worst case approaches the ill-conditioned threshold. Proceed with caution and monitor for convergence issues.

---

## Common Problems and Solutions

### Problem: "biomarker_correlation values too high"

**Symptom:** Multiple parameter combinations fail with non-PD errors

**Solution:**
```r
# OLD: Too many high values
biomarker_correlation = c(0.0, 0.3, 0.5, 0.7)  # 0.7 fails

# NEW: Reduced maximum
biomarker_correlation = c(0.0, 0.2, 0.4, 0.6)  # All pass
```

**Why:** The biomarker-response correlation is limited by the dimension of the problem. With 26 variables, c.bm cannot exceed ~0.6.

### Problem: "Some matrices ill-conditioned"

**Symptom:** All combinations pass PD test, but κ > 100 for some

**Solution Option 1:** Accept and monitor
```r
# Proceed with simulation but check for convergence warnings
# If lmer() has convergence issues, these parameters may be the cause
```

**Solution Option 2:** Exclude ill-conditioned combinations
```r
valid_well_conditioned <- validation_result$valid_combinations %>%
  filter(condition_numbers < 100)

param_grid <- valid_well_conditioned
```

**Solution Option 3:** Reduce correlation parameters
```r
# In model_params, reduce cross-correlations:
model_params$c.cf1t <- 0.15  # Was 0.2
model_params$c.cfct <- 0.08  # Was 0.1
```

### Problem: Validation is slow

**Cause:** Testing many parameter combinations with large designs

**Solution:** Reduce the parameter grid before full validation
```r
# Test a subset first
param_grid_small <- expand_grid(
  n_participants = 70,
  biomarker_correlation = c(0.2, 0.4, 0.6),  # 3 values
  carryover_t1half = c(0, 1.0)                # 2 values
)
# Much faster: 6 combinations instead of 12

# Once you find the boundary, expand strategically
param_grid_full <- expand_grid(
  n_participants = 70,
  biomarker_correlation = c(0.2, 0.3, 0.4, 0.5, 0.6),
  carryover_t1half = c(0, 0.5, 1.0, 1.5, 2.0)
)
```

---

## Integration with Existing Scripts

To add validation to `simulation_clustered.R`:

```r
# After line where param_grid is defined:
param_grid <- expand_grid(
  biomarker_moderation = c(0.2, 0.4),
  carryover_t1half = c(0, 1.0, 2.0),
  n_participants = 70
)

# ADD VALIDATION (new):
validation_result <- validate_parameter_grid(
  param_grid, trial_design_hybrid, model_params,
  resp_param, baseline_param
)

# Filter to valid combinations only
param_grid <- validation_result$valid_combinations

# Continue with simulation as before
# The param_grid now only contains valid combinations
```

---

## Best Practices

1. **Always validate before long simulations**
   - Takes minutes upfront, saves hours during simulation

2. **Check condition numbers**
   - κ < 100 ensures numerical stability
   - Monitor convergence warnings if κ > 100

3. **Document your validation**
   - Save the validation report with your results
   - Include with methods section of paper

4. **Iterative refinement**
   - Start with conservative parameter ranges
   - Use validation to identify the boundary
   - Expand grid strategically around the boundary

5. **Batch processing**
   - Validate once per design/trial configuration
   - Reuse results across multiple simulation runs
   - Cache validation results

---

## Example: Full Workflow

```r
#!/usr/bin/env Rscript

library(tidyverse)
source("pm_functions.R")

# 1. Setup (define trial_design, model_params, etc.)
# ... [omitted for brevity]

# 2. Create broad parameter grid
param_grid <- expand_grid(
  n_participants = 70,
  biomarker_correlation = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
  carryover_t1half = c(0, 1.0, 2.0)
)

# 3. VALIDATE
cat("Validating parameter combinations...\n")
validation <- validate_parameter_grid(
  param_grid, trial_design, model_params,
  resp_param, baseline_param
)

# 4. FILTER
param_grid <- validation$valid_combinations
cat(sprintf("Proceeding with %d valid combinations\n", nrow(param_grid)))

# 5. SIMULATE
cat("Running simulation...\n")
# [simulation code continues...]
```

---

## Reference: Function Signatures

```r
# Main validation function
validate_parameter_grid(
  param_grid,          # tibble: parameter combinations
  trial_design,        # tibble: trial design
  model_params,        # list: fixed correlation parameters
  resp_param,          # tibble: response parameters
  baseline_param,      # tibble: baseline parameters
  verbose = TRUE       # logical: print output?
) -> list(
  valid_combinations,  # tibble: passing combinations
  invalid_combinations,# tibble: failing combinations
  n_valid,            # integer: count of valid
  n_invalid,          # integer: count of invalid
  condition_numbers,  # numeric: κ values
  invalid_reasons     # character: failure reasons
)

# Reporting function
report_parameter_validation(
  validation_result,   # output from validate_parameter_grid()
  param_grid          # original parameter grid
) -> (invisible NULL, prints to console)
```

---

## See Also

- `pm_functions.R`: Source of validation functions
- `example_parameter_validation.R`: Complete worked example
- `positive_definiteness_constraints.tex`: Mathematical theory
