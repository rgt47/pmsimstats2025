# Parameter Validation Function Output

## Overview

The `validate_parameter_grid()` function produces **two types of output**:

1. **Console Output** (printed to screen during execution)
2. **Return Value** (R list object for programmatic use)

---

## EXAMPLE 1: All Combinations Valid

### Console Output

```
================================================================================
PRE-SIMULATION PARAMETER VALIDATION
================================================================================

Total combinations tested: 12
✓ Valid combinations:       12 (100.0%)
✗ Invalid combinations:     0 (0.0%)

Condition number (κ) statistics for valid combinations:
  Mean:   34.2
  Median: 31.8
  Min:    15.3 (best conditioned)
  Max:    58.9 (worst conditioned)

================================================================================
✓ ALL PARAMETER COMBINATIONS ARE VALID
  Ready to proceed with simulation!
================================================================================

```

### Return Value (R Object)

```r
validation_result <- validate_parameter_grid(...)

# validation_result is now a list containing:

$valid_combinations
# A tibble: 12 × 3
   n_participants biomarker_correlation carryover_t1half
            <dbl>                  <dbl>             <dbl>
 1             70                    0.2                 0
 2             70                    0.2               1
 3             70                    0.2               2
 4             70                    0.4                 0
 5             70                    0.4               1
 6             70                    0.4               2
 7             70                    0.6                 0
 8             70                    0.6               1
 9             70                    0.6               2
10             70                    0.8                 0
11             70                    0.8               1
12             70                    0.8               2

$invalid_combinations
# A tibble: 0 × 3  (empty - no failures)

$n_valid
[1] 12

$n_invalid
[1] 0

$condition_numbers
 [1] 15.3 18.2 20.1 28.5 31.8 34.2 45.6 48.3 52.1 55.2 56.8 58.9

$invalid_reasons
character(0)  # Empty - no failures
```

**Accessing components:**
```r
# Get valid combinations
valid_grid <- validation_result$valid_combinations

# Check how many failed
if (validation_result$n_invalid > 0) {
  cat("Some combinations failed\n")
}

# Get condition numbers for stability analysis
kappa_values <- validation_result$condition_numbers
mean_kappa <- mean(kappa_values)
```

---

## EXAMPLE 2: Some Combinations Invalid

### Console Output

```
================================================================================
PRE-SIMULATION PARAMETER VALIDATION
================================================================================

Total combinations tested: 15
✓ Valid combinations:       12 (80.0%)
✗ Invalid combinations:     3 (20.0%)

PROBLEMATIC COMBINATIONS:
────────────────────────────────────────────────────────────────────────────────
  Row 13: Non-positive definite matrix (c.bm=0.80)
  Row 14: Non-positive definite matrix (c.bm=0.80)
  Row 15: Non-positive definite matrix (c.bm=0.80)

Details of invalid combinations:
# A tibble: 3 × 3
  n_participants biomarker_correlation carryover_t1half
           <dbl>                  <dbl>             <dbl>
1             70                    0.8                 0
2             70                    0.8               1.0
3             70                    0.8               2.0

Condition number (κ) statistics for valid combinations:
  Mean:   38.5
  Median: 35.2
  Min:    16.4 (best conditioned)
  Max:    64.3 (worst conditioned)

RECOMMENDATIONS:
────────────────────────────────────────────────────────────────────────────────
  • All failures are due to biomarker correlation (c.bm) values being too high
  • Consider reducing the maximum biomarker_correlation in param_grid
  • Current constraint: c.bm ≤ 0.6
  • You may need to reduce to: c.bm ≤ 0.4 or lower

  • Alternatively, adjust correlation structure parameters:
    - Reduce c.cf1t (same-time cross-correlation, currently 0.2)
    - Reduce c.cfct (different-time cross-correlation, currently 0.1)
    - Increase c.autocorr (would require modifying Hendrickson parameters)

================================================================================
⚠ SOME PARAMETER COMBINATIONS ARE INVALID
  3 combinations excluded from simulation
================================================================================

```

### Return Value (R Object)

```r
validation_result <- validate_parameter_grid(...)

$valid_combinations
# A tibble: 12 × 3
   n_participants biomarker_correlation carryover_t1half
            <dbl>                  <dbl>             <dbl>
 1             70                    0.2                 0
 2             70                    0.2               1
 3             70                    0.2               2
 4             70                    0.4                 0
 5             70                    0.4               1
 6             70                    0.4               2
 7             70                    0.6                 0
 8             70                    0.6               1
 9             70                    0.6               2
10             70                   0.65                 0
11             70                   0.65               1
12             70                   0.65               2
# ... (12 valid rows)

$invalid_combinations
# A tibble: 3 × 3
  n_participants biomarker_correlation carryover_t1half
           <dbl>                  <dbl>             <dbl>
1             70                    0.8                 0
2             70                    0.8               1.0
3             70                    0.8               2.0

$n_valid
[1] 12

$n_invalid
[1] 3

$condition_numbers
 [1] 16.4 18.5 20.3 28.6 32.1 35.2 45.7 48.4 52.3 55.8 59.1 64.3

$invalid_reasons
[1] "Row 13: Non-positive definite matrix (c.bm=0.80)"
[2] "Row 14: Non-positive definite matrix (c.bm=0.80)"
[3] "Row 15: Non-positive definite matrix (c.bm=0.80)"
```

**Accessing components:**
```r
# Get the failed combinations
failed <- validation_result$invalid_combinations
cat("Failed parameters:\n")
print(failed)

# Get the reasons for failure
cat("Why they failed:\n")
for (reason in validation_result$invalid_reasons) {
  cat(sprintf("  • %s\n", reason))
}

# Filter original grid to valid combinations
valid_param_grid <- validation_result$valid_combinations
```

---

## EXAMPLE 3: Some Matrices Ill-Conditioned (Warning)

### Console Output

```
================================================================================
PRE-SIMULATION PARAMETER VALIDATION
================================================================================

Total combinations tested: 12
✓ Valid combinations:       12 (100.0%)
✗ Invalid combinations:     0 (0.0%)

⚠ WARNING [Row 5]: Ill-conditioned matrix (κ = 142.3, c.bm = 0.50)
⚠ WARNING [Row 9]: Ill-conditioned matrix (κ = 168.5, c.bm = 0.50)
⚠ WARNING [Row 11]: Ill-conditioned matrix (κ = 175.2, c.bm = 0.50)

Condition number (κ) statistics for valid combinations:
  Mean:   98.4
  Median: 87.3
  Min:    34.2 (best conditioned)
  Max:    175.2 (worst conditioned)

================================================================================
✓ ALL PARAMETER COMBINATIONS ARE VALID
  Ready to proceed with simulation!
================================================================================

```

**What this means:**
- All combinations pass the PD test ✓
- But some have poor numerical conditioning (κ > 100) ⚠
- You might see convergence warnings during simulation for those rows
- Consider filtering these out or monitoring more carefully

### Return Value

```r
$valid_combinations
# All 12 combinations (all pass PD test)

$n_valid
[1] 12

$condition_numbers
 [1]  34.2  52.1  68.3 142.3  95.6 115.2 168.5  78.9  87.3  71.4 175.2 102.8

# Note: Some κ values > 100 - these are the ill-conditioned ones
```

**Programmatic filtering:**
```r
# Keep only well-conditioned combinations
well_conditioned <- validation_result$valid_combinations %>%
  filter(condition_numbers < 100)

cat(sprintf("Filtering to well-conditioned: %d -> %d combinations\n",
            nrow(validation_result$valid_combinations),
            nrow(well_conditioned)))
```

---

## Return Value Structure (Complete Reference)

### As an R List

```r
validation_result <- list(

  # [1] Tibble of parameter combinations that passed validation
  valid_combinations = tibble(
    n_participants = numeric,
    biomarker_correlation = numeric,
    carryover_t1half = numeric
    # ... any other columns in original param_grid
  ),

  # [2] Tibble of parameter combinations that failed
  invalid_combinations = tibble(
    n_participants = numeric,
    biomarker_correlation = numeric,
    carryover_t1half = numeric
  ),

  # [3] Integer: count of valid combinations
  n_valid = integer,

  # [4] Integer: count of invalid combinations
  n_invalid = integer,

  # [5] Numeric vector: condition numbers (κ) for each valid combination
  #     Length = n_valid
  condition_numbers = numeric,

  # [6] Character vector: reason for each failure
  #     Length = n_invalid
  #     Format: "Row X: Non-positive definite matrix (c.bm=Y.YY)"
  invalid_reasons = character
)
```

### Accessing Values

```r
# Access each element
valid_grid <- validation_result$valid_combinations
invalid_grid <- validation_result$invalid_combinations
n_valid <- validation_result$n_valid
n_invalid <- validation_result$n_invalid
kappas <- validation_result$condition_numbers
reasons <- validation_result$invalid_reasons

# Check status
if (validation_result$n_invalid > 0) {
  cat("Some combinations failed!\n")
  print(validation_result$invalid_combinations)
}

# Analyze condition numbers
mean_kappa <- mean(validation_result$condition_numbers)
max_kappa <- max(validation_result$condition_numbers)

if (max_kappa > 100) {
  cat("Warning: Some matrices are ill-conditioned\n")
}
```

---

## Practical Example: Using the Output

```r
# Run validation
validation <- validate_parameter_grid(
  param_grid, trial_design, model_params,
  resp_param, baseline_param
)

# Case 1: No failures, proceed
if (validation$n_invalid == 0) {
  cat("✓ All combinations valid. Proceeding with simulation.\n")
  param_grid_final <- validation$valid_combinations
}

# Case 2: Some failures, filter and report
if (validation$n_invalid > 0) {
  cat(sprintf("⚠ %d combinations failed. Filtering...\n", validation$n_invalid))

  # Show what failed
  cat("\nFailed combinations:\n")
  print(validation$invalid_combinations)

  cat("\nReasons for failure:\n")
  for (reason in validation$invalid_reasons) {
    cat(sprintf("  • %s\n", reason))
  }

  # Use only valid ones
  param_grid_final <- validation$valid_combinations

  cat(sprintf("\nProceeding with %d valid combinations\n", nrow(param_grid_final)))
}

# Case 3: Check numerical stability
kappas <- validation$condition_numbers
if (max(kappas) > 100) {
  cat("\n⚠ Some matrices are ill-conditioned (κ > 100)\n")
  cat("  Monitor for convergence warnings during simulation\n")
}

# Now run simulation with final validated grid
source("simulation_clustered.R")
```

---

## Summary Table

| Output Component | Type | When Present | Use Case |
|------------------|------|--------------|----------|
| `valid_combinations` | Tibble | Always | Feed into simulation |
| `invalid_combinations` | Tibble | If n_invalid > 0 | Debug why parameters failed |
| `n_valid` | Integer | Always | Check how many passed |
| `n_invalid` | Integer | Always | Check if any failed |
| `condition_numbers` | Numeric vector | Always (for valid) | Assess numerical stability |
| `invalid_reasons` | Character vector | If n_invalid > 0 | Understand failure modes |

---

## Console Output vs. Return Value

| Aspect | Console Output | Return Value |
|--------|----------------|--------------|
| **When** | During execution | After execution |
| **How** | Printed to screen | R list object |
| **Use** | User feedback, logging | Programmatic filtering |
| **Persistence** | Lost unless redirected | Can be saved/reused |
| **Detail Level** | Human-readable summary | Complete raw data |

**Best practice:** Save the console output AND the return value for complete documentation.

```r
# Capture both
sink("validation_report.txt")  # Start logging to file
validation <- validate_parameter_grid(...)
sink()  # Stop logging

# Now have both:
# - validation_report.txt (console output)
# - validation (R object)
```
