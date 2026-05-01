# Two-Scenario Comparison: With vs. Without Carryover Modeling

## Overview

We now run simulations under **TWO analysis approaches** to answer different scientific questions:

### Scenario 1: WITH Carryover Modeling (Our Enhancement)
**Question**: *"What's the statistical power when you PROPERLY model carryover effects?"*

**Analysis Model**:
```r
response ~ treatment * biomarker + week + carryover_effect + (1|participant_id)
```

**Expected Result**: Power should remain **STABLE** across carryover levels because the model controls for the confound.

---

### Scenario 2: WITHOUT Carryover Modeling (Hendrickson Approach)
**Question**: *"What's the power when you IGNORE carryover effects?"*

**Analysis Model**:
```r
response ~ treatment * biomarker + week + (1|participant_id)
```

**Expected Result**: Power should **DECLINE** as carryover increases because it becomes an unmodeled confound.

---

## Implementation Details

### Data Generation (SAME for both scenarios)
- Pure Hendrickson approach
- Biomarker interaction emerges from correlation structure only
- No population mean shift
- No post-MVN adjustment
- Carryover effects added to bio_response means when `carryover_t1half > 0`

### Analysis Models (DIFFERENT)

**Scenario 1** (`model_carryover = TRUE`):
- Includes `carryover_effect` covariate
- Uses exponential decay: `(1/2)^(scale * weeks_off / t1half)`
- Matches data generation carryover function
- **Perfect model specification**

**Scenario 2** (`model_carryover = FALSE`):
- Does NOT include `carryover_effect`
- Carryover variance goes into residual error
- **Intentional model misspecification** (when carryover exists in data)

---

## Parameter Grid

```r
param_grid <- expand_grid(
  n_participants = c(70),
  biomarker_correlation = c(0.2, 0.4),
  carryover_t1half = c(0, 1.0, 2.0),  # None, moderate, strong
  treatment_effect = c(5.0)
)
```

**Total simulations per scenario**:
- 1 sample size × 2 biomarker correlations × 3 carryover levels = 6 combinations
- × 2 designs (hybrid, crossover) = 12 combinations
- × 20 iterations = 240 simulation runs
- × 2 analysis approaches = **480 total simulations**

---

## Expected Patterns

### When carryover_t1half = 0 (No Carryover)
**Both scenarios should show SAME power** because:
- No carryover in data
- No difference between models

### When carryover_t1half > 0 (Carryover Exists)

**Scenario 1 (WITH modeling)**:
- Power remains stable or slightly decreases
- Small decrease due to:
  - Extra parameter (carryover term) uses 1 degree of freedom
  - Slightly reduced power from DF loss
- But confound is controlled

**Scenario 2 (WITHOUT modeling)**:
- Power DECLINES significantly
- Larger decrease due to:
  - Carryover creates unexplained variance
  - Inflated standard errors
  - Reduced t-statistics
  - Unmodeled confound

**The DIFFERENCE between scenarios shows the COST of not modeling carryover.**

---

## Results Structure

### simulation_results
Each row represents one Monte Carlo iteration with columns:
- `design`: "hybrid" or "crossover"
- `n_participants`: 70
- `biomarker_correlation`: 0.2 or 0.4
- `carryover_t1half`: 0, 1.0, or 2.0
- `model_carryover`: TRUE or FALSE  ← **Key variable**
- `iteration`: 1 to 20
- `effect_size`: Estimated biomarker×treatment interaction
- `std_error`: Standard error of estimate
- `t_value`: t-statistic
- `p_value`: Two-sided p-value
- `significant`: TRUE if p < 0.05
- `error`: TRUE if model convergence failed

### simulation_summary
Aggregated power estimates:
- Groups by: design, n_participants, biomarker_correlation, carryover_t1half, **model_carryover**
- `power`: Proportion of significant results (main outcome)
- `mean_effect`: Average effect size
- `mean_se`: Average standard error
- `n_iterations`: Number of successful iterations

---

## How to Interpret Results

### Compare Within Design and Carryover Level

**Example**: Hybrid design, carryover_t1half = 1.0

| model_carryover | power | mean_se | interpretation |
|-----------------|-------|---------|----------------|
| TRUE            | 0.85  | 0.42    | Controlled carryover |
| FALSE           | 0.65  | 0.54    | Unmodeled confound |

**Interpretation**:
- **Power difference**: 0.85 - 0.65 = 0.20 (20 percentage points!)
- **SE inflation**: 0.54 / 0.42 = 1.29 (29% larger when not modeling)
- **Cost of ignoring carryover**: Substantial power loss

### Look for Hendrickson's Pattern

**Scenario 2** (WITHOUT modeling) should show power decline across carryover:

| carryover_t1half | power |
|------------------|-------|
| 0                | 0.80  |
| 1.0              | 0.65  |
| 2.0              | 0.55  |

This matches Hendrickson's published findings!

### Validate Our Enhancement

**Scenario 1** (WITH modeling) should show stable power:

| carryover_t1half | power |
|------------------|-------|
| 0                | 0.80  |
| 1.0              | 0.78  |
| 2.0              | 0.76  |

Small decline only from DF loss, not confounding.

---

## Scientific Contributions

### 1. Replicates Hendrickson's Finding
Scenario 2 demonstrates that **ignoring carryover reduces power**, validating her original results.

### 2. Shows Benefit of Our Enhancement
Scenario 1 demonstrates that **modeling carryover maintains power**, proving the value of our methodological improvement.

### 3. Quantifies the Cost
The difference between scenarios shows **exactly how much power you lose** by not accounting for carryover.

### 4. Design Comparison Across Approaches
We can compare hybrid vs. crossover under BOTH analysis strategies, showing which design is more robust to model misspecification.

---

## Key Implementation Changes

### Modified run_monte_carlo()
```r
run_monte_carlo <- function(design_name, params, sigma_cache = NULL,
                            model_carryover = TRUE) {
  # ...
  if (params$carryover_t1half > 0 && model_carryover) {
    model <- lmer(response ~ treatment * bm + week + carryover_effect +
                  (1|participant_id))
  } else {
    model <- lmer(response ~ treatment * bm + week + (1|participant_id))
  }
  # ...
}
```

### Main Simulation Loop
```r
for (design_name in c("hybrid", "crossover")) {
  for (model_carryover in c(TRUE, FALSE)) {
    design_results <- run_monte_carlo(
      design_name, current_params, sigma_cache,
      model_carryover = model_carryover
    )
    # ...
  }
}
```

### Results Tracking
All results include `model_carryover` column to distinguish scenarios.

---

## Visualization Recommendations

### 1. Power by Carryover (Faceted by Scenario)
```r
ggplot(simulation_summary,
       aes(x = carryover_t1half, y = power, color = design)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ model_carryover,
             labeller = labeller(model_carryover =
               c("TRUE" = "With Carryover Model",
                 "FALSE" = "Without Carryover Model (Hendrickson)"))) +
  labs(title = "Impact of Carryover Modeling on Statistical Power",
       x = "Carryover Half-life (weeks)",
       y = "Power to Detect Biomarker Interaction")
```

### 2. Power Loss from Not Modeling
```r
power_comparison <- simulation_summary %>%
  select(design, carryover_t1half, model_carryover, power) %>%
  pivot_wider(names_from = model_carryover,
              values_from = power,
              names_prefix = "model_") %>%
  mutate(power_loss = model_TRUE - model_FALSE)

ggplot(power_comparison,
       aes(x = carryover_t1half, y = power_loss, color = design)) +
  geom_line() +
  geom_point() +
  labs(title = "Power Lost by Not Modeling Carryover",
       y = "Power Difference (With - Without Modeling)")
```

### 3. Standard Error Inflation
```r
ggplot(simulation_summary,
       aes(x = carryover_t1half, y = mean_se,
           color = design, linetype = model_carryover)) +
  geom_line() +
  labs(title = "Standard Error Inflation from Unmodeled Carryover")
```

---

## Expected Conclusions

1. **Hendrickson was correct**: Not modeling carryover reduces power
2. **Our enhancement works**: Modeling carryover maintains power
3. **Quantified impact**: We can measure exactly how much power is lost
4. **Design comparison**: Hybrid vs. crossover performance under both approaches
5. **Methodological recommendation**: Always model carryover when it exists!

---

## Connection to Previous Documentation

- See `HENDRICKSON_VS_OUR_DETAILED_COMPARISON.md` for implementation differences
- See `WHY_POST_MVN_ADJUSTMENT_ANALYSIS.md` for biomarker interaction mechanism
- See `HENDRICKSON_CARRYOVER_ANALYSIS.md` for original findings

---

*Documentation created: 2025-11-18*
