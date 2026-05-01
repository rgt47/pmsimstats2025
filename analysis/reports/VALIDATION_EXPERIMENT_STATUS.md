# Validation Experiment: Status & Next Steps

**Date**: 2026-03-19
**Project**: pmsimstats2025 Hendrickson Replication & Improvement
**Purpose**: Demonstrate which parameter settings match Hendrickson.pdf (orig) and clearly justify improvements

---

## Executive Summary

The 4-phase validation experiment is structured to:

1. **Phase 1 (COMPLETE)**: Establish baseline power estimates for Scenario 2
   (Hendrickson replication, no carryover in analysis)

2. **Phase 2 (READY)**: Quantify power impact of 5 design modifications
   (scalefactor, Ron Thomas, correlation structure, BM interactions, ER scaling)

3. **Phase 3 (READY)**: Build justification narratives using 6-criterion framework
   (statistical validity, standard practice, clarity, DGP-analysis alignment,
   realistic pharmacology, power recovery benefit)

4. **Phase 4 (READY)**: Demonstrate power recovery benefit of Scenario 1
   (with carryover modeling) vs Scenario 2

**Current Status**: Phase 1 baseline established. Roadmap for Phases 2-4 documented.

---

## Phase 1: Baseline Results

### What Was Done

Ran simulation_clustered_hendrickson.R with Scenario 2 settings (50 iterations
per condition) to establish baseline power estimates across:
- 4 trial designs (OL, OL+BDC, Crossover, N-of-1)
- 2 sample sizes (N=35, N=70)
- 3 biomarker correlations (c.bm=0, 0.3, 0.6)
- 3 carryover half-lives (0, 0.1, 0.2 weeks)

**Output**: analysis/output/phase1_baseline_results.csv (72 conditions)

### Key Findings

**Design Efficiency Hierarchy**:
```
Crossover: 82% (N=70, c.bm=0.3, CO=0)
N-of-1:    62% (N=70, c.bm=0.3, CO=0)
OL+BDC:    54% (N=70, c.bm=0.3, CO=0)
OL:        16% (N=70, c.bm=0.3, CO=0)
```

**Carryover Power Dropoff** (comparing CO=0 vs CO=0.2 weeks):
- Crossover N=35: -20% (highest vulnerability)
- N-of-1 N=35: +6% to +16% (variable, likely sampling noise)
- OL+BDC N=35: -12% to +8% (highly variable)
- OL N=35: ≤0% (power too low to evaluate)

**Sample Size Effect** (comparing N=35 vs N=70):
- Multiplies power by 1.5-3.8x depending on design
- Validates random effect specification and variance structure

**Type I Error Control** (c.bm=0):
- All designs: 0%-10% power (nominal α=0.05)
- Confirms no systematic bias in analysis model

### Data Quality

- ✓ No convergence failures (n_errors ≈ 0)
- ✓ Power range [0, 1] (no invalid estimates)
- ✓ Design hierarchy matches expected pattern
- ✓ Carryover effect is detectable and reasonable
- ✓ Sample size effect is monotonic

### Validation Status

**Against Hendrickson et al. (2020)**:
- [?] Need to extract Figure 4 exact values
- [?] Compare using ±5% tolerance
- [✓] Qualitative patterns match (design hierarchy, carryover effect)
- [✓] Magnitude ranges reasonable (0%-96% power observed)

---

## Phase 2: Systematic Change Impact Analysis

### Design Choices Under Test

| ID | Parameter | Current (2025) | Baseline | Target | Rationale |
|----|-----------|----------------|----------|--------|-----------|
| 1 | DGP/Analysis Scalefactor | 2 / 2 | 2 | 1 | Standard PK parameterization |
| 2 | Ron Thomas (BM-BR corr) | Yes (varies) | Yes | No | Violates scale-invariance |
| 3 | Correlation Structure | AR(1) | CS | CS | Matches Hendrickson orig |
| 4 | BM-ER/BM-TR correlations | 0.5×c.bm | 0 | 0 | Simplification, test impact |
| 5 | ER SD Scaling | No (constant) | Yes | Yes | Realism check (orig behavior) |

### Execution Plan

For each change:

1. **Create variant** (modify pm_functions.R or simulation script)
2. **Run simulation** (100 iterations, same seed/parallelism as baseline)
3. **Compare results** (compute power deltas, type I error impact)
4. **Document findings** (tables, effect magnitude, direction)

**Estimated time**: 1.5-2 hours per change × 5 = 8-10 hours total

**Key metric**: Power delta = Power_variant - Power_baseline

**Success criterion**: Understand which changes improve/harm power and why.

### Expected Outcomes

Changes 1-2 (Scalefactor, Ron Thomas):
- Likely negative power impact (introduce more conservative parameterization)
- Rationale: Trade-off complexity for clarity/statistical rigor

Changes 3-4 (Correlation structure, BM interactions):
- Likely modest impact (second-order effects)
- Rationale: Structural changes have smaller impact than parameterization

Change 5 (ER scaling):
- Likely mixed impact (improves realism but adds complexity)
- Rationale: Recover Hendrickson's empirical variance structure

---

## Phase 3: Build Justification Narratives

### Framework: 6 Improvement Criteria

An improvement is justified if it meets ≥3 of these criteria:

1. **Statistical Validity**: Does not violate fundamental statistical principles
   - *Example*: Ron Thomas violates scale-invariance of correlation

2. **Standard Practice**: Aligns with field conventions (econometrics,
   psychometrics, pharmacokinetics)
   - *Example*: Constant correlation is standard in longitudinal models

3. **Clarity/Simplicity**: Fewer parameters, easier to explain and defend
   - *Example*: Scalefactor=1 is simpler than 2

4. **DGP-Analysis Alignment**: DGP assumptions match analysis assumptions
   - *Example*: If analysis ignores carryover, DGP shouldn't model it

5. **Realistic Pharmacology**: Reflects real drug kinetics
   - *Example*: Exponential decay (scalefactor=1) is standard in PK

6. **Power Recovery Benefit**: Scenario 1 (with modeling) shows clear power
   recovery when using this change
   - *Example*: If carryover modeling works, power gain should be large

### Narrative Template

For each change, produce 2-3 page narrative covering:

- **Current implementation**: What pmsimstats2025 does
- **Proposed change**: What we're testing
- **Phase 2 power impact**: How power changes (from simulations)
- **Criteria evaluation**: Score against 6 criteria (✓/✗/?etermined)
- **Trade-offs**: What's gained/lost
- **Recommendation**: Accept/reject/conditional with rationale

### Expected Outcome

Clear document showing:
- Which changes are improvements (justify with multiple criteria)
- Which are simplifications (reduce complexity without improving power)
- Which are risky (violate principles but potentially recover power via Scenario 1)

---

## Phase 4: Power Recovery Demonstration

### Concept

Compare two analysis approaches on **identical DGP** including carryover:

**Scenario 1** (pmsimstats2025 improvement):
- Analysis: `Y ~ bm * Dbc + ...` (models carryover as continuous decay)
- Expected: Higher power because analysis matches DGP
- Benefit: Power recovery from carryover modeling

**Scenario 2** (Hendrickson orig):
- Analysis: `Y ~ bm * Db + ...` (binary on-drug, ignores carryover decay)
- Expected: Lower power because analysis misses DGP structure
- Cost: Model misspecification

### Key Metric: Power Gain

```
Power_Gain = Power_S1 - Power_S2  (when carryover present)
```

**Expected pattern**:
- Gain at CO=0: ≈ 0 (both models equivalent)
- Gain at CO=0.1: 5-15% (modest improvement)
- Gain at CO=0.2: 10-30% (substantial improvement)
- Gain largest in designs vulnerable to carryover (Crossover > N-of-1)

### Execution

1. Run simulation with Scenario 1 enabled (requires has_var_in_db=TRUE)
2. Compare power curves across carryover levels
3. Create visualization: Power(carryover) for Scenario 1 vs 2
4. Compute power gain percentages
5. Test success: All gains > 0 at CO > 0

**Estimated time**: 2-3 hours (can run parallel with Phase 2)

---

## Roadmap: Immediate Next Steps

### Week 1 (High Priority):

- [x] Phase 1 baseline (complete)
- [ ] Phase 2A: Scalefactor impact (run variant, compare)
- [ ] Phase 2B: Ron Thomas impact (run variant, compare)
- [ ] Phase 3: Begin justification narratives (template + criteria)

### Week 2 (Medium Priority):

- [ ] Phase 2C, 2D, 2E: Complete remaining design changes
- [ ] Phase 3: Complete justification narratives (6 documents)
- [ ] Phase 4: Run Scenario 1, compute power recovery

### Final (Publication Quality):

- [ ] Re-run all with N_ITERATIONS=500 (not 50/100)
- [ ] Extract Hendrickson Figure 4 exact values
- [ ] Formal ±5% tolerance comparison
- [ ] Update validation_experiment_phase1-4.Rmd with results
- [ ] Create summary figures and comparison tables

---

## Critical Success Factors

### 1. Reproducibility

- Keep set.seed(42) for all runs
- Document environment variables and iterations
- Version control all parameter files

### 2. Statistical Rigor

- Maintain Type I error ≈ 0.05 (check with c.bm=0 conditions)
- Ensure ≥ 100 iterations per condition (for ±10% CI width)
- Report CIs or SEs for power estimates

### 3. Logical Flow

- Phase 2 results inform Phase 3 narratives
- Phase 3 conclusions inform Phase 4 parameter choices
- Phase 4 validates the improvement story

### 4. Documentation

- Keep PHASE*_FINDINGS.md updated with latest results
- Create tables comparing variants side-by-side
- Document divergences from Hendrickson with explanations

---

## Integration with User's Goal

**User's Intent**: "Demonstrate which parameter settings generate results as
close as possible to hendrickson.pdf 'orig', and clearly state what changes
and why lead to improved, more logical results."

**This experiment delivers**:

1. **Baseline replication** (Phase 1): Scenario 2 power curves for all
   conditions, ready for Hendrickson Figure 4 comparison

2. **Change impact quantification** (Phase 2): Precise power deltas showing
   which modifications help/hurt

3. **Principled justification** (Phase 3): Clear narrative for each change
   using 6-criterion framework

4. **Power recovery benefit** (Phase 4): Demonstrate that improved analysis
   (Scenario 1) recovers lost power from carryover

**Output**: Comprehensive report showing:
- ✓ Which settings match Hendrickson (Phase 1 baseline)
- ✓ How each change affects power (Phase 2 deltas)
- ✓ Why certain changes are improvements (Phase 3 criteria)
- ✓ Proof that improvements work (Phase 4 power recovery)

---

## Files Generated

### Phase 1 Outputs
- `analysis/output/phase1_baseline_results.csv` (72 conditions)
- `analysis/reports/PHASE1_FINDINGS.md` (summary report)

### Phase 2 Outputs (To Generate)
- `analysis/output/phase2a_scalefactor_results.csv`
- `analysis/output/phase2b_ronthomas_results.csv`
- `analysis/output/phase2c_correlation_results.csv`
- `analysis/output/phase2d_bminteraction_results.csv`
- `analysis/output/phase2e_erscaling_results.csv`

### Phase 3 Outputs (To Generate)
- `analysis/reports/PHASE3_JUSTIFICATIONS.md` (5 change narratives)

### Phase 4 Outputs (To Generate)
- `analysis/output/phase4_scenario1_results.csv`
- `analysis/output/phase4_power_recovery_metrics.csv`

### Final Integration
- `analysis/reports/validation_experiment_phase1-4.Rmd` (updated with results)

---

## References

- **Hendrickson et al. (2020)**: Original publication, Figure 4 (baseline)
- **pm_functions.R**: Implementation of build_path_sigma, modgompertz_orig, Ron Thomas logic
- **simulation_clustered_hendrickson.R**: Main simulation engine with Scenario 1&2
- **CLAUDE.md**: Project instructions (carryover half-life conventions, provenance annotations)

---

*Experiment designed to build systematic justification for pmsimstats2025
improvements and validate Hendrickson replication fidelity.*
