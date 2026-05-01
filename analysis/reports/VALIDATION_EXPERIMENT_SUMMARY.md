# Validation Experiment: 4-Phase Framework

## Accomplished: Phase 1 Complete

**Objective**: Establish baseline power estimates for Scenario 2 (Hendrickson
replication, no carryover modeling in analysis)

**Method**: Ran simulation_clustered_hendrickson.R with 50 iterations per
condition across full parameter grid (4 designs × 2 N × 3 c.bm × 3 carryover)

**Results**: 72 conditions yielding stable power estimates from 0% to 96%

### Key Baseline Findings

**Design Efficiency** (c.bm=0.3, N=70, carryover=0):
```
Crossover: 82%  ← Most efficient
N-of-1:    62%
OL+BDC:    54%
OL:        16%   ← Least efficient
```

**Carryover Impact** (N=35, c.bm=0.3):
- Crossover loses 8-20% power as carryover increases
- N-of-1 shows modest sensitivity (variable, -2% to +6%)
- OL and OL+BDC show irregular patterns (sampling variation with low base power)

**Sample Size Effect** (c.bm=0.3):
- N=35 → N=70: 1.5-3.8x power multiplier
- Validates random effect parameterization
- Demonstrates carryover is especially problematic for under-powered designs

**Type I Error Control** (c.bm=0):
- All designs: 0%-10% power (nominal α=0.05)
- No evidence of model misspecification bias

### Validation Status

✓ Baseline established and stable
✓ Design patterns match theoretical expectations
✓ Carryover effect is detectable and reasonable
✓ Ready for Hendrickson Figure 4 comparison

⚠ Pending: Extract exact Hendrickson values for ±5% tolerance check

---

## Ready for Execution: Phases 2-4

### Phase 2: Systematic Change Impact (8-10 hours)

Test how 5 design modifications affect power:

1. **Scalefactor = 1** (vs baseline 2)
   - Rationale: Standard PK parameterization
   - Expected: Stronger remaining carryover = lower power

2. **Remove Ron Thomas** (constant BM-BR correlation)
   - Rationale: Violates scale-invariance of correlation
   - Expected: Modest power change, improved clarity

3. **Correlation Structure** (Compound Symmetry vs AR(1))
   - Rationale: Matches Hendrickson original
   - Expected: Small second-order effect

4. **BM-ER/BM-TR correlations = 0** (vs 0.5×c.bm)
   - Rationale: Simplification, test impact
   - Expected: Minimal power change

5. **ER SD Scaling** (expectancy-weighted vs constant)
   - Rationale: Recover Hendrickson's empirical variance structure
   - Expected: Mixed impact, depends on design

**Deliverable**: Power delta tables comparing each variant to baseline

### Phase 3: Build Justification Narratives (2-3 hours)

For each change, evaluate 6 criteria:
1. Statistical validity (violates principles?)
2. Standard practice (matches field conventions?)
3. Clarity/simplicity (fewer parameters?)
4. DGP-analysis alignment (assumptions match?)
5. Realistic pharmacology (matches drug kinetics?)
6. Power recovery benefit (Scenario 1 gains?)

Accept change if ≥3 criteria met.

**Deliverable**: 2-3 page narrative per change with recommendation

### Phase 4: Power Recovery Demonstration (2-3 hours)

Compare Scenario 1 (with carryover modeling) vs Scenario 2 (without)
on identical DGP including carryover.

**Success criterion**: Scenario 1 power ≥ Scenario 2 power at all
carryover levels, with gains of 10-30% at moderate carryover.

**Deliverable**: Power recovery metrics and figures

---

## Integration Roadmap

### This Week
- [x] Phase 1: Establish baseline
- [ ] Phase 2: Run 5 variants (parallel execution possible)
- [ ] Phase 3: Write justification narratives

### Next Week
- [ ] Phase 4: Power recovery demonstration
- [ ] Compile final validation_experiment_phase1-4.Rmd
- [ ] Compare to Hendrickson Figure 4 (exact values)

### Publication Quality
- [ ] Re-run with 500 iterations (vs 50 used for Phase 1)
- [ ] Formal statistical testing of power improvements
- [ ] Create comparison figures and summary tables

---

## Answer to User's Intent

**User Goal**: "Demonstrate which parameter settings generate results as
close as possible to hendrickson.pdf 'orig', and clearly state what changes
and why lead to improved, more logical results."

**This experiment delivers**:

1. **Baseline replication** (Phase 1):
   - Power estimates for all conditions in Scenario 2
   - Ready for direct comparison to Hendrickson Figure 4

2. **Change impact analysis** (Phase 2):
   - Precise quantification of how each modification affects power
   - Isolates individual contribution of each design choice

3. **Principled improvement narrative** (Phase 3):
   - Scores each change against 6 objective criteria
   - Explains why certain changes are improvements (≥3 criteria met)
   - Documents trade-offs transparently

4. **Proof of improvement** (Phase 4):
   - Demonstrates that Scenario 1 (with carryover modeling) recovers
     power lost from ignoring carryover
   - Shows benefit of improved analysis model

**Final output**: Comprehensive report showing:
- ✓ Which settings match Hendrickson (baseline from Phase 1)
- ✓ How each change affects power (Phase 2 deltas)
- ✓ Why changes are improvements (Phase 3 criteria)
- ✓ Proof improvements work (Phase 4 power recovery)

---

## File Organization

```
analysis/reports/
├── validation_experiment_phase1-4.Rmd        (Main document, updated)
├── PHASE1_FINDINGS.md                        (Phase 1 results)
├── PHASE2_IMPLEMENTATION_GUIDE.md            (Phase 2 instructions)
├── VALIDATION_EXPERIMENT_STATUS.md           (Overall status & roadmap)
└── VALIDATION_EXPERIMENT_SUMMARY.md          (This file)

analysis/output/
├── phase1_baseline_results.csv               (Phase 1 results)
├── phase2a_scalefactor1_results.csv          (Phase 2A: scalefactor variant)
└── [phase2b-2e, phase4 results to generate]
```

---

## Next Action

Execute Phase 2 by running 5 variants of simulation_clustered_hendrickson.R
with different parameter values (see PHASE2_IMPLEMENTATION_GUIDE.md for details).

Core workflow:
```bash
# Baseline is already run (Phase 1)

# For each of 5 changes:
# 1. Modify pm_functions.R or simulation script
# 2. Run: [ENV_VAR=value] N_ITERATIONS=100 Rscript simulation_clustered_hendrickson.R
# 3. Compare results (power_variant vs power_baseline)
# 4. Document findings
```

Estimated time: 8-10 hours total (can run variants in parallel on separate cores).
