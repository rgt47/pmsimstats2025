# Phase 2 Execution Guide: Crossover & N-of-1 Focus

**Objective**: Demonstrate that modifications reduce power dropoff as carryover increases

**Baseline metric**: Power curve stability = 1 - |Δpower from t1/2=0 to 0.2|
- Higher stability = flatter curve = more measured power reduction

---

## Baseline Pattern (from Phase 1)

### Crossover (CO), N=35, c.bm=0.3
```
t1/2 = 0.0:  Power = 62%
t1/2 = 0.1:  Power = 42%  (drop 20%)
t1/2 = 0.2:  Power = 54%  (net drop 8%)
Stability = 0.92
```
**Problem**: Sharp power loss at t1/2=0.1, only partial recovery at 0.2

### N-of-1 (Hybrid), N=35, c.bm=0.3
```
t1/2 = 0.0:  Power = 24%
t1/2 = 0.1:  Power = 30%  (gain 6%)
t1/2 = 0.2:  Power = 40%  (net gain 16%)
Stability = 0.84
```
**Observation**: More stable, but variability suggests sampling noise at low power

---

## Phase 2: Variant Testing Strategy

### Key Insight

The baseline shows CO is vulnerable to carryover (steep power loss), while N-of-1 is more resilient. Modifications should:

1. **Reduce power dropoff magnitude** in CO (goal: stability > 0.95)
2. **Stabilize power trajectory** by improving DGP-analysis alignment
3. **Show consistent direction** across both designs

### Execution Plan

For **each variant**:

1. Modify `analysis/scripts/pm_functions.R`
2. Run: `N_ITERATIONS=150 timeout 600 Rscript analysis/scripts/simulation_clustered_hendrickson.R`
3. Extract CO and N-of-1 results
4. Compute stability metric
5. Compare to baseline
6. Document power delta

---

## Variant Implementation Details

### Variant A: Scalefactor = 1 (vs baseline 2)

**File**: `pm_functions.R` line 55-56 in `simulation_clustered_hendrickson.R`

**Current**:
```r
dgp_scalefactor <- 2
analysis_scalefactor <- 2
```

**Change to**:
```r
dgp_scalefactor <- 1
analysis_scalefactor <- 1
```

**Rationale**: Standard PK parameterization

**Expected outcome**:
- Scalefactor=1 means carryover decay exponent is smaller
- Result: (1/2)^(1 × tsd / t1half) vs (1/2)^(2 × tsd / t1half)
- Smaller exponent = slower decay = MORE carryover strength remaining
- **Hypothesis**: More carryover in DGP but not in analysis model = WORSE power
- **Expected stability**: < 0.92 (worse than baseline)

### Variant B: Remove Ron Thomas Adjustment

**File**: `pm_functions.R` lines 1567-1574 in `build_path_sigma()`

**Current** (Ron Thomas):
```r
bm_br_corrs <- build_bm_br_correlations(
  brmeans, brtest, c.bm = c.bm, nP = n_timepoints
)
# This varies correlation based on BR means
```

**Change to** (constant correlation):
```r
bm_br_corrs <- rep(c.bm, n_timepoints)
# Use constant c.bm at all timepoints, don't scale by BR means
```

**Rationale**: Ron Thomas violates scale-invariance of correlation

**Expected outcome**:
- Removes complex mean-based scaling
- Correlation structure becomes more predictable
- **Hypothesis**: Cleaner covariance structure = steadier power estimates
- **Expected stability**: Similar to baseline or slightly improved

### Variant C: Correlation Structure (AR(1) → Compound Symmetry)

**File**: `pm_functions.R` lines 1490-1515 in `build_path_sigma()`

**Current** (AR(1)):
```r
# Time-based correlation: rho^|i-j|
# Newer timepoints more independent
```

**Change to** (Compound Symmetry):
```r
# Constant correlation between all pairs (except diagonal)
# All off-diagonal elements = same value
```

**Rationale**: Matches Hendrickson original

**Expected outcome**:
- AR(1) gives temporal structure; CS gives uniform correlations
- **Hypothesis**: CS removes temporal assumptions, may stabilize power
- **Expected stability**: Unclear (second-order effect)

### Variant D: Combine A + B

Run both Scalefactor=1 AND No Ron Thomas simultaneously

**Expected**: If A and B have independent effects, stability should improve (or worsen) additively

---

## Execution Workflow

### Step 1: Create variant versions of pm_functions.R

```bash
# Backup original
cp analysis/scripts/pm_functions.R analysis/scripts/pm_functions_ORIGINAL.R

# Create variants
cp analysis/scripts/pm_functions.R analysis/scripts/pm_functions_variant_A.R  # scalefactor=1
cp analysis/scripts/pm_functions.R analysis/scripts/pm_functions_variant_B.R  # no Ron Thomas
cp analysis/scripts/pm_functions.R analysis/scripts/pm_functions_variant_C.R  # CS instead AR1
```

### Step 2: Apply changes to each variant

Edit each `pm_functions_variant_X.R` with the modifications above.

### Step 3: Run simulations

```bash
# Baseline (already done in Phase 1)
# N_ITERATIONS=150 Rscript analysis/scripts/simulation_clustered_hendrickson.R

# Variant A: Scalefactor = 1
sed -i 's/dgp_scalefactor <- 2/dgp_scalefactor <- 1/' analysis/scripts/simulation_clustered_hendrickson.R
sed -i 's/analysis_scalefactor <- 2/analysis_scalefactor <- 1/' analysis/scripts/simulation_clustered_hendrickson.R
N_ITERATIONS=150 timeout 600 Rscript analysis/scripts/simulation_clustered_hendrickson.R
# Capture output as phase2a_variant_A_results.RData
mv analysis/output/simulation_clustered_hendrickson_results.RData analysis/output/phase2a_variant_A_results.RData

# Restore baseline scalefactors
git checkout analysis/scripts/simulation_clustered_hendrickson.R

# Repeat for Variant B, C, etc.
```

### Step 4: Compare variants

```r
# Load all results
baseline <- read_csv('analysis/output/phase1_baseline_co_nof1.csv')
load('analysis/output/phase2a_variant_A_results.RData')
variant_a <- summary_results %>%
  filter(scenario == 'S2_no_carryover', design %in% c('CO', 'N-of-1'))

# Compute stability for each variant
compute_stability <- function(data) {
  data %>%
    filter(c.bm == 0.3) %>%
    pivot_wider(
      id_cols = c(design, N),
      names_from = carryover_t1half,
      values_from = power
    ) %>%
    mutate(
      stability = 1 - abs(`0.2` - `0`),
      power_delta_at_0 = `0` - baseline_power_0,
      power_delta_at_02 = `0.2` - baseline_power_02
    )
}

# Create comparison table
comparison <- tibble(
  Variant = c('Baseline', 'A: sf=1', 'B: No RT', 'C: CS', 'D: A+B'),
  CO_N35_Stability = c(0.92, NA, NA, NA, NA),
  NOF1_N35_Stability = c(0.84, NA, NA, NA, NA),
  CO_N35_Power_at_0 = c(0.62, NA, NA, NA, NA),
  CO_N35_Power_at_02 = c(0.54, NA, NA, NA, NA)
)

# Fill in variant values (pseudocode)
for (variant_name in c('A', 'B', 'C', 'D')) {
  # Load variant results
  # Compute stability
  # Add to comparison table
}

# Visualize
comparison %>%
  ggplot(aes(x = Variant, y = CO_N35_Stability, fill = Variant)) +
  geom_col() +
  labs(title = 'Phase 2: Stability Comparison (CO, N=35, c.bm=0.3)',
       y = 'Stability (higher = flatter power curve)')
```

---

## Success Criteria

### For each variant:

1. **Power curves produced**: ✓ Results file created with same structure as baseline
2. **No convergence issues**: n_errors ≈ 0 (check summary_results)
3. **Type I error maintained**: c.bm=0 conditions show ~5% power
4. **Stability metric improved** (or documented if worsens):
   - Variant A: Expect stability < 0.92 (MORE carryover = WORSE)
   - Variant B: Expect stability ≈ 0.92-0.95 (modest change)
   - Variant C: Expect stability ≈ 0.92 (second-order effect)
   - Variant D: Additive effect of A + B

### Interpretation:

- **Increasing stability** = flatter power curve = better alignment of DGP and analysis
- **Decreasing stability** = steeper power curve = worse alignment (may still be useful if paired with Scenario 1 improvement)

---

## Phase 2 Deliverable

**File**: `analysis/reports/PHASE2_RESULTS_CO_NOF1.md`

```markdown
# Phase 2 Results: Crossover & N-of-1 Variants

## Summary Table

| Variant | CO (N=35) Stability | NOF1 (N=35) Stability | Interpretation |
|---------|-------------------|----------------------|-----------------|
| Baseline (sf=2, RT) | 0.92 | 0.84 | Starting point |
| A: sf=1 | ? | ? | MORE carryover strength |
| B: No Ron Thomas | ? | ? | Constant correlation |
| C: CS | ? | ? | Uniform temporal correlation |
| D: A+B | ? | ? | Combined effects |

## Detailed Findings

### Variant A: Scalefactor = 1

[Results and interpretation]

### Variant B: Remove Ron Thomas

[Results and interpretation]

### [Etc.]

## Phase 3 Connection

Based on Phase 2 results, the variants that:
- Increase stability = IMPROVE alignment
- These candidates for Phase 3 justification narrative
```

---

## Timeline

- **Phase 2A (sf=1)**: 10 min prep + 60 min simulation = 70 min
- **Phase 2B (No RT)**: 10 min prep + 60 min simulation = 70 min
- **Phase 2C (CS)**: 10 min prep + 60 min simulation = 70 min
- **Phase 2D (A+B)**: 5 min prep + 60 min simulation = 65 min
- **Phase 2 Analysis & Documentation**: 30 min

**Total Phase 2**: ~4.5 hours

---

## Next: Phase 3 Justification

After Phase 2 variants are tested, Phase 3 will:

1. Score each variant against 6 criteria
2. Rank which variants represent "improvements"
3. Build narrative explaining trade-offs

**Key question**: Does more stable power curve (less dropoff) come at cost of lower baseline power? If so, is it worth the trade-off?
