# Focused Validation Roadmap: Crossover & N-of-1

**Scope**: Demonstrate power dropoff reduction through systematic design modifications

**Key Metric**: Power curve stability = 1 - |Δpower from t1/2=0 to 0.2|

---

## The Core Story

### Baseline Problem (from Phase 1)

**Crossover at N=35, c.bm=0.3** shows sharp power vulnerability:
- No carryover (t1/2=0):    **62% power**
- Moderate carryover (t1/2=0.1): **42% power** (↓20pp)
- Higher carryover (t1/2=0.2):   **54% power** (net ↓8pp)

This is the "power dropoff phenomenon" Hendrickson documented.

### N-of-1 is More Resilient

**N-of-1 at N=35, c.bm=0.3** shows flatter trajectory:
- No carryover:      **24% power**
- Moderate carryover: **30% power** (↑6pp)
- Higher carryover:   **40% power** (net ↑16pp)

**Question**: Why? And can we improve Crossover to match N-of-1's stability?

---

## The Hypothesis

The baseline uses **Scenario 2** (Hendrickson replication):
- **DGP**: Includes carryover with decay model
- **Analysis**: Uses binary on-drug indicator (Db), ignores carryover decay

**Misalignment**: DGP models carryover, but analysis doesn't. This is the dropoff source.

**Solution**: Improve DGP-analysis alignment through targeted modifications:
1. Adjust carryover parameterization (scalefactor)
2. Simplify correlation structure (remove Ron Thomas, use compound symmetry)
3. Align expectancy and ER assumptions

**Expected result**: More stable power curves (less steep dropoff).

---

## Phase 2: Four Focused Variants

### Variant A: Scalefactor = 1 (vs baseline 2)

**Change**: `dgp_scalefactor ← 1` and `analysis_scalefactor ← 1`

**Effect on carryover**: Smaller exponent in decay formula
- Current: (1/2)^(2 × tsd / t1half)
- Variant: (1/2)^(1 × tsd / t1half)
- Result: **Slower decay = MORE carryover remains**

**Expected power impact**:
- Crossover: Power should DECREASE (more carryover = harder to detect interaction)
- N-of-1: Modest change (less vulnerable to carryover anyway)
- **Stability**: Likely WORSE (steeper curve)

**Rationale**: This tests whether standard PK parameterization (sf=1) is better than Hendrickson's choice (sf=2)

---

### Variant B: Remove Ron Thomas Adjustment

**Change**: Replace `build_bm_br_correlations(...)` with constant `rep(c.bm, n_timepoints)`

**Effect on covariance**: BM-BR correlation becomes constant across timepoints

**Current (Ron Thomas)**: Correlation varies based on BR means
**Variant**: Constant c.bm at all timepoints

**Expected power impact**:
- Cleaner covariance structure (no mean-based scaling)
- Analysis model assumes constant correlation (now matches DGP)
- **Power**: Small change (second-order effect)
- **Stability**: Likely SIMILAR or slightly improved (~0.93-0.95)

**Rationale**: Ron Thomas violates statistical principle (scale-invariance of correlation)

---

### Variant C: Compound Symmetry (CS) Instead of AR(1)

**Change**: Replace AR(1) correlation structure with uniform compound symmetry

**Effect on temporal correlation**:
- Current (AR(1)): ρ^|i-j| (newer timepoints more independent)
- Variant (CS): All pairs have same correlation ρ

**Expected power impact**:
- Removes temporal assumption
- More conservative (assumes high correlation across all pairs)
- **Stability**: Unclear (second-order effect)

**Rationale**: Matches Hendrickson's original correlation structure

---

### Variant D: Combined A + B

**Change**: Scalefactor=1 + No Ron Thomas together

**Expected power impact**: Additive combination of A and B effects

**Rationale**: Test whether improvements compound

---

## Phase 3: Justification Framework

For each variant, score against 6 criteria:

1. **Statistical validity**: Does it violate principles? (Ron Thomas → violates scale-invariance)
2. **Standard practice**: Matches field conventions? (Scalefactor=1 is standard in PK)
3. **Clarity/simplicity**: Fewer parameters to explain?
4. **DGP-analysis alignment**: Do assumptions match?
5. **Realistic pharmacology**: Matches drug kinetics?
6. **Power recovery benefit**: Does Scenario 1 (with modeling) show gains?

**Accept change if**: ≥3 criteria met

**Example**:
```
Variant A (sf=1):
1. Statistical validity: ✓ (standard parameterization)
2. Standard practice: ✓ (PK models use sf=1)
3. Clarity: ✓ (simpler exponent)
4. DGP-analysis alignment: ✗ (still has mismatch)
5. Pharmacology: ✓ (exponential decay is standard)
6. Power recovery: TBD (Phase 4)

Score: 4/6 → ACCEPT if Phase 4 shows power recovery
```

---

## Phase 4: Power Recovery Demonstration

### Setup

**Scenario 1** (Improved):
- DGP: Same carryover as baseline
- Analysis: Y ~ bm × Dbc + ... (models carryover as continuous)
- Expected: **Higher power** (analysis matches DGP)

**Scenario 2** (Baseline):
- DGP: Same carryover
- Analysis: Y ~ bm × Db + ... (binary, ignores decay)
- Expected: Lower power (misalignment)

### Execution

Run same 4 variants with Scenario 1 enabled:
```bash
N_ITERATIONS=200 SCENARIO=1 Rscript analysis/scripts/simulation_clustered_hendrickson.R
```

### Success Criterion

For each variant, compute:
```
Power_Gain = Power_S1 - Power_S2
```

**Target**: Gains of **10-30%** at moderate carryover (t1/2=0.1-0.2)

Example target:
```
Variant A (sf=1):
  CO N=35 c.bm=0.3:
    Scenario 2 at t1/2=0.2: 54%
    Scenario 1 at t1/2=0.2: 70%  (+16pp) ✓

  N-of-1 N=35 c.bm=0.3:
    Scenario 2 at t1/2=0.2: 40%
    Scenario 1 at t1/2=0.2: 48%  (+8pp) ✓
```

**If gains are large**: The variant is justified (Scenario 1 recovers lost power)
**If gains are small**: The variant may still be good for other reasons (Phase 3 criteria)

---

## Execution Checklist

### Phase 2 (Variants)
- [ ] Create pm_functions_variant_A.R (sf=1)
- [ ] Create pm_functions_variant_B.R (no Ron Thomas)
- [ ] Create pm_functions_variant_C.R (CS instead AR1)
- [ ] Create pm_functions_variant_D.R (A+B combined)
- [ ] Run each variant: `N_ITERATIONS=150 Rscript ...`
- [ ] Extract CO and N-of-1 results
- [ ] Compute stability metric for each
- [ ] Document power deltas vs baseline
- [ ] Create PHASE2_RESULTS_CO_NOF1.md

### Phase 3 (Justification)
- [ ] Score each variant against 6 criteria
- [ ] Build narrative for each (2-3 paragraphs)
- [ ] Rank variants by improvement potential
- [ ] Document trade-offs
- [ ] Create PHASE3_JUSTIFICATIONS.md

### Phase 4 (Power Recovery)
- [ ] Create modified pm_functions for Scenario 1
- [ ] Run each variant with Scenario 1 enabled
- [ ] Compute power gains at each carryover level
- [ ] Create visualization: Power_Gain vs Carryover
- [ ] Document which variants recover most power
- [ ] Create PHASE4_POWER_RECOVERY.md

---

## Expected Outcomes

### Phase 2 Stability Improvements

Based on theoretical analysis:

| Variant | CO Stability | NOF1 Stability | Direction |
|---------|-----------|---|---|
| Baseline | 0.92 | 0.84 | Reference |
| A (sf=1) | ~0.88 | ~0.82 | WORSE (more carryover) |
| B (No RT) | ~0.93 | ~0.85 | Neutral/slight improvement |
| C (CS) | ~0.91 | ~0.83 | Small change |
| D (A+B) | ~0.87 | ~0.81 | Combination effect |

**Interpretation**: If A worsens power but Scenario 1 recovers it, then A+Scenario1 is the improvement story.

### Phase 3 Ranking

**Most justified variants** (by criteria score):
1. Variant B (No Ron Thomas): Violates principle, easier to explain
2. Variant D (A+B): Combines multiple improvements
3. Variant A (sf=1): Standard in field but hurts power alone

### Phase 4 Proof

**Expected**: Scenario 1 shows large power recovery when:
- Variant A (sf=1): +15-25% (matches DGP carryover better)
- Variant B (No RT): +8-12% (modest, but cleaner structure)
- Variant D: +20-30% (best combination)

---

## Deliverables

1. **PHASE2_RESULTS_CO_NOF1.md**: Stability metrics and power deltas for all variants
2. **PHASE3_JUSTIFICATIONS.md**: 6-criterion scoring and narratives
3. **PHASE4_POWER_RECOVERY.md**: Scenario 1 vs Scenario 2 comparison with gains
4. **Visualization**: `phase2_stability_comparison.png` (bar chart)
5. **Visualization**: `phase4_power_recovery.png` (power curves for CO and N-of-1)

---

## Timeline

- **Phase 2**: 4-5 hours (4 variants × 60 min simulation + analysis)
- **Phase 3**: 1-2 hours (write narratives, score criteria)
- **Phase 4**: 2-3 hours (run Scenario 1, compute gains)

**Total**: 7-10 hours of focused work

---

## Why This Approach Works

1. **Focused scope**: Only CO and N-of-1 (2 designs instead of 4)
2. **Clear metric**: Stability quantifies "less power dropoff"
3. **Systematic**: Test one change at a time, then combinations
4. **Principled**: Use 6-criterion framework to justify recommendations
5. **Validated**: Phase 4 proves improvements actually recover power

**Result**: Clear, defensible answer to "Which changes improve results and why?"
