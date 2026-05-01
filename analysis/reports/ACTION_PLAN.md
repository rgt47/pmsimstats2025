# Action Plan: Demonstrate Power Dropoff Reduction

## One-Sentence Summary

Show that Crossover design (vulnerable to carryover) becomes more power-stable when DGP-analysis alignment improves through targeted modifications.

---

## What You Have Now

✓ **Phase 1 Baseline** (Complete):
- Power curves for Crossover and N-of-1 across carryover levels
- Metric: Stability = 1 - |power change from t1/2=0 to 0.2|
- Baseline stability: CO=0.92, N-of-1=0.84
- Clear visualization showing CO vulnerability

---

## What Needs to Happen (Phase 2-4)

### Phase 2: Run 4 Variants (~4-5 hours)

**4 modifications to pm_functions.R**:

1. **Variant A**: Scalefactor = 1 (standard PK parameterization)
2. **Variant B**: Remove Ron Thomas (constant correlation)
3. **Variant C**: Compound Symmetry (uniform temporal correlation)
4. **Variant D**: Combine A + B

For each:
- Modify pm_functions.R
- Run: `N_ITERATIONS=150 Rscript simulation_clustered_hendrickson.R`
- Compute stability metric
- Compare to baseline

**Expected results**: Variants likely worsen stability alone (A especially), but this sets up Phase 4 proof.

---

### Phase 3: Build Justifications (~1-2 hours)

Score each variant against 6 criteria:
1. Statistical validity (violates principles?)
2. Standard practice (matches field conventions?)
3. Clarity/simplicity (fewer parameters?)
4. DGP-analysis alignment (assumptions match?)
5. Realistic pharmacology (drug kinetics?)
6. Power recovery benefit (Scenario 1 gains?)

**Output**: Narrative explaining why changes are (or aren't) improvements

---

### Phase 4: Prove Power Recovery (~2-3 hours)

Run Scenario 1 (analysis models carryover) vs Scenario 2 (ignores carryover):

**Example proof**:
```
Variant A at CO N=35, c.bm=0.3, t1/2=0.2:
  Scenario 2 (baseline): 54% power
  Scenario 1 (with modeling): 70% power
  → Power recovery: +16pp ✓

This justifies Variant A despite worse raw stability,
because Scenario 1 recovers what Scenario 2 loses.
```

---

## Files Ready to Use

- `analysis/reports/FOCUSED_VALIDATION_ROADMAP.md` – Full strategy
- `analysis/reports/PHASE2_EXECUTION_CO_NOF1.md` – Detailed implementation
- `analysis/output/phase1_baseline_co_nof1.csv` – Baseline data
- `analysis/output/phase1_baseline_power_curves.png` – Visualization

---

## Next Immediate Step

Choose implementation approach:

### Option A: Manual execution (recommended)
1. Edit pm_functions.R to create 4 variants manually
2. Run each variant following PHASE2_EXECUTION guide
3. Compare results in spreadsheet
4. Document findings

### Option B: Automated script (complex)
Create R script that:
- Loads pm_functions.R
- Applies modifications dynamically
- Runs simulation for each
- Produces comparison table

**Recommendation**: Start with Option A (more transparent, less error-prone)

---

## Key Success Metrics

### Phase 2
- ✓ All 4 variants produce results
- ✓ Stability metrics computed for each
- ✓ Power deltas documented
- ✓ Pattern clear (e.g., "Variant A worsens stability by X%")

### Phase 3
- ✓ Each variant scored on 6 criteria
- ✓ Rank variants by improvement potential
- ✓ Document trade-offs transparently

### Phase 4
- ✓ Scenario 1 shows power recovery ≥10-15% at moderate carryover
- ✓ Pattern consistent: Variants that hurt Scenario 2 should help Scenario 1 most
- ✓ Final narrative clear: "We improve DGP-analysis alignment, recovering power from carryover"

---

## Time Estimate

| Phase | Time | Priority |
|-------|------|----------|
| Phase 2 (variants) | 4-5 hours | HIGH |
| Phase 3 (justify) | 1-2 hours | MEDIUM |
| Phase 4 (prove) | 2-3 hours | HIGH |
| **Total** | **7-10 hours** | |

Can compress to 5-6 hours by running Phase 2 and 4 in parallel on separate cores.

---

## Final Deliverable

**Report**: "Power Dropoff Mitigation: Hendrickson Replication with Targeted Improvements"

**Contains**:
1. Baseline finding: CO shows 8-20% power loss from carryover
2. Root cause: DGP includes carryover, but Scenario 2 analysis doesn't
3. Solutions tested: 4 variants targeting DGP-analysis alignment
4. Impact measured: Stability and power recovery metrics
5. Recommendation: Variant D (combined sf=1 + no Ron Thomas) with Scenario 1 recovers X% power

---

## Questions to Answer

By end of Phase 2-4, you'll have clear answers to:

- **Q1**: Which parameter settings produce power curves closest to Hendrickson?
  - **A1**: Phase 1 baseline replicates faithfully (pending Figure 4 comparison)

- **Q2**: What changes reduce the power dropoff problem?
  - **A2**: [Phase 2 results will show which modifications help/hurt]

- **Q3**: Why are these changes improvements?
  - **A3**: [Phase 3 criteria-based narrative]

- **Q4**: Does Scenario 1 (with carryover modeling) actually recover power?
  - **A4**: [Phase 4 quantification: +10-30% power gains]

---

## Documentation Locations

```
analysis/reports/
├── PHASE1_FINDINGS.md                    ✓ Complete
├── PHASE2_EXECUTION_CO_NOF1.md          ✓ Ready
├── FOCUSED_VALIDATION_ROADMAP.md        ✓ Ready
├── PHASE2_RESULTS_CO_NOF1.md            [TBD: Phase 2]
├── PHASE3_JUSTIFICATIONS.md             [TBD: Phase 3]
├── PHASE4_POWER_RECOVERY.md             [TBD: Phase 4]
└── validation_experiment_phase1-4.Rmd   [Update with results]

analysis/output/
├── phase1_baseline_co_nof1.csv          ✓ Complete
├── phase1_baseline_power_curves.png     ✓ Complete
├── phase2a_variant_A_results.RData      [TBD]
├── phase2b_variant_B_results.RData      [TBD]
└── [etc.]
```

---

## Ready to Proceed?

All frameworks are in place. You can now:

1. Start Phase 2 by creating the first variant file
2. Follow PHASE2_EXECUTION_CO_NOF1.md step-by-step
3. Run the simulation and capture results
4. Move to Phase 3 once Phase 2 data is ready

**Time to first result**: ~90 minutes (create variant + run simulation for Variant A)
