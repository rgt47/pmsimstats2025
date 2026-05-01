# Phase 2-4 Implementation Guide

**Purpose**: Document the systematic approach for testing design changes and
building improvement justifications.

**Status**: Ready for implementation
**Baseline**: Phase 1 complete (50 iterations, Scenario 2 baseline established)

---

## Phase 2: Systematic Change Impact Analysis

### Overview

Phase 2 tests five design modifications to measure how each affects power
curves compared to Phase 1 baseline (scalefactor=2, Ron Thomas enabled, etc.).

### Architecture

Rather than manually reimplement tests for each variant, use parametrized
simulation approach:

```
Create variants of simulation_clustered_hendrickson.R by:

1. Adding environment variables to control each parameter
2. Running each variant with identical seed/random stream as baseline
3. Capturing results to separate CSV files
4. Computing power deltas relative to baseline
```

### Change 1: Scalefactor = 1 (vs baseline 2)

**Location**: Line 55 and 226-227 of simulation_clustered_hendrickson.R

**Change required**:
```r
dgp_scalefactor <- as.integer(Sys.getenv('DGP_SCALEFACTOR', '2'))
analysis_scalefactor <- as.integer(Sys.getenv('ANALYSIS_SCALEFACTOR', '2'))
```

**Execution**:
```bash
# Baseline (sf=2)
N_ITERATIONS=100 timeout 600 Rscript analysis/scripts/simulation_clustered_hendrickson.R

# Variant (sf=1)
DGP_SCALEFACTOR=1 ANALYSIS_SCALEFACTOR=1 N_ITERATIONS=100 \
  timeout 600 Rscript analysis/scripts/simulation_clustered_hendrickson.R

# Compare power deltas
# Expected: Smaller scalefactor = stronger remaining carryover = lower power
```

**Expected pattern**: Power decrease at higher carryover levels (tsd > 0).

### Change 2: Remove Ron Thomas Adjustment

**Location**: Lines 1567-1574 in pm_functions.R (build_path_sigma)

**Change required**:
```r
# Current (Ron Thomas):
bm_br_corrs <- build_bm_br_correlations(
  brmeans, brtest, c.bm = c.bm, nP = n_timepoints
)

# Simplified (constant correlation):
bm_br_corrs <- rep(c.bm, n_timepoints)  # Or set to 0 if no correlation
```

**Execution**:
1. Create pm_functions_no_ronthomas.R (copy of pm_functions.R with Ron Thomas logic removed)
2. Source this variant in simulation
3. Run comparison test

**Expected pattern**: Modest power change (Ron Thomas has second-order effect).
Hypothesis: Removing Ron Thomas should increase power (less correlation
scaling noise), but effect magnitude unclear.

### Change 3: Correlation Structure (AR(1) → Compound Symmetry)

**Location**: Lines 1490-1515 in pm_functions.R (sigma matrix construction)

**Current**: AR(1) time-based correlation structure
**Target**: Compound symmetry (constant correlation across all timepoints)

**Change required**:
```r
# Replace AR(1) block with compound symmetry:
# sigma[i,j] = 1 if i==j, rho otherwise (for same factor type)
```

**Expected pattern**: Modest power change (temporal correlation structure
has second-order effect on test statistics).

### Change 4: BM-ER/BM-TR Correlations (0.5×c.bm → 0)

**Location**: Lines 1560-1576 in pm_functions.R

**Change required**:
```r
# Current: Creates additional correlations between BM and ER/TR
# Simplified: Set BM-ER and BM-TR correlations to 0
```

**Expected pattern**: Likely small power change (reduces within-subject
variance structure).

### Change 5: ER SD Scaling (constant → expectancy-weighted)

**Location**: Line 1358 in pm_functions.R (hendrickson_bl_params)

**Change required**:
```r
# Currently: ER has constant SD across all expectancy levels
# Target: Scale ER SD by expectancy factor (from Hendrickson orig)

er_sd <- bl_params$er$s * expectancies  # Scale by expectancy
```

**Expected pattern**: Larger change in power, particularly for designs with
varying expectancy (Hybrid, OL+BDC).

---

## Phase 2 Execution Workflow

### For each of 5 changes:

1. **Modify code**: Update pm_functions.R or simulation script
2. **Run baseline**: `N_ITERATIONS=100 Rscript simulation_clustered_hendrickson.R`
3. **Run variant**: `[PARAMETER] N_ITERATIONS=100 Rscript simulation_clustered_hendrickson.R`
4. **Compare results**:
   ```r
   baseline <- read_csv('simulation_clustered_hendrickson_results.csv')
   variant <- read_csv('simulation_clustered_hendrickson_results_variant.csv')

   comparison <- baseline %>%
     left_join(variant, by = c('design', 'N', 'c.bm', 'carryover_t1half'),
               suffix = c('_base', '_var')) %>%
     mutate(
       power_delta = power_var - power_base,
       pct_change = (power_delta / power_base) * 100
     )
   ```

5. **Document findings**:
   - Power delta at each condition
   - Overall impact magnitude
   - Direction of change (increased/decreased power)
   - Designs most affected

---

## Phase 3: Build Justification Narrative

### For each change, evaluate 6 criteria:

1. **Statistical validity**: Violates principle? (e.g., Ron Thomas violates
   scale-invariance of correlation)

2. **Standard practice**: Matches field conventions? (e.g., constant
   correlation is standard in econometrics)

3. **Clarity/simplicity**: How many new concepts to explain?

4. **DGP-analysis alignment**: Do DGP and analysis model agree?

5. **Realistic pharmacology**: Matches known drug kinetics?

6. **Power recovery benefit**: Does Scenario 1 (with modeling) show clear
   power recovery using this change?

### Template for each change:

```markdown
## Change [N]: [Name]

**Current Implementation**: [Details]
**Proposed**: [Details]

### Meets Criteria?

1. Statistical validity: [✓/✗] [Why]
2. Standard practice: [✓/✗] [Why]
3. Clarity/simplicity: [✓/✗] [Why]
4. DGP-analysis alignment: [✓/✗] [Why]
5. Realistic pharmacology: [✓/✗] [Why]
6. Power recovery benefit: [TBD] [Awaiting Phase 4]

### Phase 2 Impact on Power

[Insert comparison table from Phase 2]

**Interpretation**: [Explain delta in biological/statistical terms]

### Recommendation

- **Accept**: [Justification if ≥3 criteria met]
- **Reject**: [Justification if criteria not met]
- **Conditional**: [If pending additional evidence]
```

---

## Phase 4: Two-Scenario Power Recovery Demonstration

### Purpose

Demonstrate that Scenario 1 (with carryover modeling) provides power
recovery compared to Scenario 2 (ignoring carryover).

### Design

For each design, sample size, and carryover level:

**Scenario 1** (Hendrickson 2025 improved):
- DGP: Includes carryover (scalefactor, decay model)
- Analysis: `lmer(...~ bm * Dbc + ...)`  (models carryover as continuous)
- Power: Estimate via simulation

**Scenario 2** (Hendrickson orig):
- DGP: Same carryover as Scenario 1
- Analysis: `lmer(...~ bm * Db + ...)` (binary on-drug, no carryover decay)
- Power: From Phase 1 baseline

**Comparison**:
```r
power_recovery <- scenario1_results %>%
  left_join(scenario2_baseline,
            by = c('design', 'N', 'c.bm', 'carryover_t1half')) %>%
  mutate(
    power_gain = power_s1 - power_s2,
    pct_gain = (power_gain / power_s2) * 100
  ) %>%
  filter(carryover_t1half > 0)  # Only for conditions with carryover
```

### Success Criteria

- Scenario 1 power > Scenario 2 power at all carryover levels
- Power gain ≥ 10-15% at moderate carryover (t1/2 = 0.1-0.2)
- Gain is larger for designs that suffer most from carryover
  (Crossover > N-of-1 > OL+BDC)

### Execution

```bash
# Run Scenario 1 with final parameter choices
N_ITERATIONS=500 SCENARIO=1 Rscript analysis/scripts/simulation_clustered_hendrickson.R

# Compare to Phase 1 Scenario 2 baseline
Rscript -e "
  s1 <- read_csv('scenario1_results.csv')
  s2 <- read_csv('phase1_baseline_results.csv')
  # ... compute power recovery metrics
"
```

---

## Recommended Execution Order

**Phase 2** (Change impact, 1-2 hours per change × 5 = 5-10 hours total):
1. Change 1: Scalefactor
2. Change 2: Ron Thomas
3. Change 3: Correlation structure
4. Change 4: BM interactions
5. Change 5: ER scaling

**Phase 3** (Justification narratives, 30 min per change):
- Document findings and build improvement cases

**Phase 4** (Power recovery, 2-3 hours):
- Run Scenario 1 with final parameter choices
- Quantify power recovery benefits

---

## Critical Implementation Notes

### Reproducibility

- **Seed**: Keep set.seed(42) consistent across all variants
- **RNG Stream**: Use same parallelism approach for all runs
- **Iterations**: Use 100-200 for Phase 2 (faster than 500), scale to 500 for final publication

### Validation Checkpoints

After each Phase 2 change:
- [ ] Results file created
- [ ] Power estimates within reasonable range (0-100%)
- [ ] No models failed to converge (n_errors ≈ 0)
- [ ] Type I error maintained (~5% when c.bm=0)

### Documentation

Create summary table after Phase 2 completion:

| Change | Design | Power Delta | Type | Justification |
|--------|--------|-------------|------|---------------|
| SF1    | CO     | -12% (at CO=0.2) | Harm | Stronger carryover remains |
| SF1    | NO     | +2% | Neutral | Minimal effect |
| ... | ... | ... | ... | ... |

---

## Deliverables

1. **phase2_change1_scalefactor_results.csv** (5 files total)
2. **PHASE3_JUSTIFICATIONS.md** (narrative for each change)
3. **phase4_scenario1_comparison.csv** (power recovery metrics)
4. **validation_experiment_FINAL.Rmd** (populated with all results)

---

## Timeline Estimate

- Phase 2: 8-12 hours
- Phase 3: 2-3 hours
- Phase 4: 2-3 hours
- **Total**: 12-18 hours of computation + analysis

Recommend running Phases 2 and 4 in parallel (separate R sessions).
