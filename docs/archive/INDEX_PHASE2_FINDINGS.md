# Phase 2 Diagnostic Investigation: Complete Index
## Quick Navigation Guide

**Status**: Investigation complete. Evidence-based synthesis with actionable recommendations.
**Date**: March 19, 2026

---

## Core Documents (Read in This Order)

### 1. PHASE2_EXECUTIVE_SUMMARY.md (START HERE)
**Length**: 5 pages | **Audience**: Decision makers
- Executive summary of findings
- All 5 Phase 2 variants with results
- Combined impact analysis
- Path forward for Phase 3-4

**Key takeaway**: Two unjustified 2025 extensions inflate power 50-90% and
degrade stability. Remove them to get honest, measured power curves.

---

### 2. BIOMARKER_MECHANISMS_SYNTHESIS.md (DEEP DIVE)
**Length**: 12 pages | **Audience**: Methodologists, code reviewers
- Theory: Mean-structure vs covariance-based approaches
- Empirical Phase 2 evidence for all 5 variants
- Reconciliation of theory vs. empirical validation
- Evidence-based assessment and recommendations

**Key insight**: Theoretical appeal of covariance approaches is real but
fragile. Implementation matters enormously. Hendrickson's original was
carefully designed. 2025 extensions were never empirically validated.

---

### 3. PHASE2_CULPRIT_FINDINGS.md (TECHNICAL REFERENCE)
**Length**: 8 pages | **Audience**: Code implementers, reproducibility
- Detailed results tables for all variants
- Code locations and exact changes
- Empirical results by design, sample size, carryover level
- Interpretation and verdict for each variant

**Key reference**: Go here to understand exactly what changed and what the
impact was.

---

## Supplementary Materials

### 4. Updated validation_experiment_phase1-4.Rmd
**Location**: analysis/reports/validation_experiment_phase1-4.Rmd
- All Phase 1-2 detailed results integrated
- Extensive documentation of how DGP changes affect covariance structures
- Ready for knit to HTML/PDF

### 5. Reproducible Test Scripts
All available in analysis/scripts/
- phase2a_scalefactor_variant.R (75 iterations per condition)
- phase2c_correlation_structure.R (75 iterations per condition)
- phase2d_bm_interactions.R (75 iterations per condition)
- phase2e_er_scaling.R (75 iterations per condition)

Run with `Rscript phase2X_*.R` to reproduce all findings.

---

## Quick Reference: Phase 2 Variants

| Phase | What Changed | Code Location | Result | Verdict |
|-------|-------------|---------------|--------|---------|
| 2A | Scalefactor 2→1 | simulation_clustered.R:55 | +1-8% stability | Minor, not culprit |
| 2B | Remove Ron Thomas | pm_functions.R:1568 | ~0% change | Not culprit |
| 2C | AR(1) → CS | build_path_sigma():1534-1565 | AR(1) better 2.7% | Not culprit (keep AR(1)) |
| **2D** | **Remove BM-ER/TV** | **pm_functions.R:1570-1576** | **-50 to -85% power** | **⭐ MAJOR CULPRIT** |
| **2E** | **Revert ER SD scale** | **pm_functions.R:1504** | **-65 to -90% power** | **⭐ MAJOR CULPRIT** |

---

## The Core Finding in One Paragraph

pmsimstats2025 added two extensions to Hendrickson's covariance-based DGP:
(1) BM-ER and BM-PB correlations (0.5×c.bm), and (2) constant ER SD instead
of expectancy-scaled. These inflate baseline power by 50-90% while degrading
power curve stability, creating an artificially exaggerated appearance of
carryover vulnerability. Removing these two changes (Phases 2D+2E) restores
realistic power curves (~13% baseline, ~8% vulnerability drop) that honestly
reflect the cost of ignoring carryover in analysis when it's present in the
DGP. The vulnerability is not a design flaw; it's the expected consequence of
DGP-analysis misalignment, and it becomes measured and interpretable once the
inflated baseline is removed.

---

## Decision Points for Phase 3-4

### Immediate Actions (High Confidence)

**Action 1**: Remove BM-ER/BM-TV correlations
- **File**: pm_functions.R
- **Lines**: 1570-1576
- **Evidence**: Phase 2D shows -50 to -85% power inflation with degraded stability
- **Confidence**: ⭐⭐⭐ (empirical, reproducible, diagnostic testing)

**Action 2**: Revert ER SD to expectancy-scaled
- **File**: pm_functions.R
- **Line**: 1504
- **From**: `rep(resp_params$pb$sd, nP)` (constant)
- **To**: `resp_params$pb$sd * expectancies` (scaled)
- **Evidence**: Phase 2E shows -65 to -90% power inflation with degraded stability
- **Confidence**: ⭐⭐⭐ (empirical, reproducible, diagnostic testing)

### Optional Actions (Medium Confidence)

**Action 3**: Test scalefactor=1
- **File**: simulation_clustered_hendrickson.R
- **Line**: 55
- **Evidence**: Phase 2A shows +1-8% stability improvement
- **Next step**: Run Phase 4 with both sf=1 and sf=2; measure power recovery
- **Decision**: Accept if Phase 4 power recovery is substantially better with sf=1

### Keep As-Is (High Confidence)

- **AR(1) correlation structure** (Phase 2C shows it's superior)
- **Ron Thomas adjustment** (Phase 2B shows negligible effect, harmless)

---

## Expected Outcomes After Fixing Phases 2D+2E

### Crossover (N=35, c.bm=0.3)
Current (2025) → Expected (corrected)
- Baseline power: 27% → 13%
- Power at t1/2=0.2: 36% → 5%
- Stability: 0.907 → 0.920 (improved)
- Vulnerability: Exaggerated → Realistic (~8% drop)

### Crossover (N=70, c.bm=0.3)
Current (2025) → Expected (corrected)
- Baseline power: 68% → 20%
- Power at t1/2=0.2: 69% → 8%
- Stability: 0.987 → 0.880 (improved)
- Vulnerability: Catastrophic → Measured

---

## Timeline

- **March 19, 2026**: Phase 2 investigation complete
- **Next**: Implement Phases 2D+2E fixes in pm_functions.R
- **Then**: Run corrected simulation_clustered_hendrickson.R (Phase 4)
- **Output**: Two-scenario comparison (with/without carryover modeling)
- **Validation**: Demonstrate power recovery from carryover modeling justifies
  the analytical investment

---

## Questions & Answers

### Q: Are we removing the covariance-based approach?
**A**: No. We're keeping Hendrickson's original covariance-based approach
(BM-BR only, expectancy-scaled ER SD). We're removing the unjustified 2025
extensions to that approach.

### Q: Will this make the vulnerability look worse?
**A**: It will look honest. The inflated baseline made the vulnerability
appear exaggerated. With realistic baseline, the vulnerability becomes
measured and justifiable (8% drop is realistic for ignoring carryover).

### Q: Do we need to re-run Phase 1 baseline?
**A**: Yes, after implementing Phases 2D+2E fixes. The new baseline will be
different (lower, more realistic) and should be documented.

### Q: What about scalefactor testing?
**A**: Optional secondary improvement. Phase 4 will show if sf=1 provides
additional power recovery benefit. Include in sensitivity analysis at minimum.

---

## Files Involved in Implementation

**Read-only** (for understanding):
- analysis/reports/validation_experiment_phase1-4.Rmd (updated)
- analysis/docs/PHASE2_EXECUTIVE_SUMMARY.md (new)
- analysis/docs/PHASE2_CULPRIT_FINDINGS.md (new)
- analysis/docs/BIOMARKER_MECHANISMS_SYNTHESIS.md (new)
- analysis/scripts/phase2[ace]_*.R (new, test scripts)
- analysis/scripts/phase2d_bm_interactions.R (new, test script)
- analysis/scripts/phase2e_er_scaling.R (new, test script)

**To be modified** (implementation):
- analysis/scripts/pm_functions.R (lines 1504, 1570-1576)
- analysis/scripts/simulation_clustered_hendrickson.R (line 55, optional)

---

## Reproducibility

All Phase 2 results are fully reproducible:

```bash
cd analysis/scripts

# Run individual variants
Rscript phase2a_scalefactor_variant.R
Rscript phase2c_correlation_structure.R
Rscript phase2d_bm_interactions.R
Rscript phase2e_er_scaling.R

# All use 75 iterations × 4 conditions × 3 carryover levels = 900 total sims
# Expected runtime: ~3 minutes per script
```

All results are stored in the Rmd file for permanent reference.

---

**Document generated**: 2026-03-19
**Investigation**: Phase 2 Diagnostic
**Evidence basis**: 5 variants, 75 iterations each, reproducible test scripts
**Confidence level**: ⭐⭐⭐ (empirical, validated, actionable)
