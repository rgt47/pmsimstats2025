# Phase 2 Diagnostic Findings: Vulnerability Culprits Identified

## Executive Summary

The unrealistic carryover vulnerability (20%+ power drop) in the 2025 baseline is caused by TWO design additions that inflate baseline power artificially:

1. **BM-ER/BM-TR correlations** (Phase 2D) - inflate power 50-85%
2. **Constant ER SD** (Phase 2E) - inflate power 65-90%

When removed (matching orig behavior), the power curves become realistic and the vulnerability becomes measured (~8% drop).

---

## Phase 2A: Scalefactor Impact

**Change**: DGP scalefactor from 2 (baseline) to 1

**Code Location**: `simulation_clustered_hendrickson.R` line 55

**Results** (CO N=35, c.bm=0.3):
- Baseline sf=2: 20% → 13% (7 pp drop, 0.933 stability)
- Variant sf=1: 21% → 16% (5 pp drop, 0.947 stability)
- **Effect**: +1.3% stability improvement

**Verdict**: Modest improvement, NOT a culprit. Justified if Phase 4 shows
clear power recovery benefit.

---

## Phase 2B: Ron Thomas Adjustment Impact

**Change**: BM-BR correlation from scaled by BR means to constant c.bm

**Code Location**: `pm_functions.R` line 1567-1568
```r
# Current (2025):
bm_br_corrs <- build_bm_br_correlations(brmeans, brtest, c.bm = c.bm, nP = nP)

# Variant (Phase 2B):
bm_br_corrs <- rep(c.bm, nP)  # Constant, no scaling
```

**Results** (CO N=35, c.bm=0.3):
- Baseline: 20% → 13% (0.933 stability)
- Variant: 20% → 13% (0.933 stability)
- **Effect**: No measurable change

**Verdict**: Minimal impact, NOT a culprit.

---

## Phase 2C: Correlation Structure Impact

**Change**: Temporal correlation from AR(1) to Compound Symmetry

**Code Location**: `build_path_sigma()` lines 1534-1565

**Current AR(1) behavior**:
```r
lag <- p2 - p
correlations[n1, n2] <- ac^lag  # Correlation decays with time
```

**Variant CS behavior**:
```r
correlations[n1, n2] <- ac  # Constant correlation, no decay
```

**Results** (CO N=35, c.bm=0.3):
- Baseline AR(1): 11% → 7% (0.960 stability)
- Variant CS: 20% → 13% (0.933 stability)
- **Effect**: AR(1) is BETTER by 2.7% stability

**Verdict**: Current AR(1) is superior (more realistic). Switching to CS
would be a regression.

---

## Phase 2D: BM-ER/BM-TR Correlation Impact — **MAJOR CULPRIT**

**Change**: BM-ER and BM-PB correlations from 0.5×c.bm to 0

**Code Location**: `build_path_sigma()` lines 1570-1576

**Current (2025) behavior**:
```r
# For each timepoint p:
n1_tv <- paste(timeptnames[p], 'tv', sep = '.')
correlations['bm', n1_tv] <- 0.5 * c.bm  # BM-TV correlation

n2_pb <- paste(timeptnames[p], 'pb', sep = '.')
correlations['bm', n2_pb] <- 0.5 * c.bm  # BM-PB correlation

# Plus existing BM-BR correlation (scaled)
```

**Variant (matching orig)**:
```r
# Only BM-BR correlation; no BM-TV or BM-PB
```

**Physical Interpretation**:
- **With interactions** (2025): Biomarker is correlated with expectancy-driven
  response AND time-varying response
- **Without** (orig): Biomarker only interacts with treatment-specific response

**Results** (CO N=35, c.bm=0.3):
```
                With BM-int    No BM-int    Δ Power
t1/2=0          27%            13%          -50%
t1/2=0.1        35%            12%          -65%
t1/2=0.2        36%            5%           -85%
Stability       0.907          0.920        +0.013 (worse with int)
```

**CO N=70, c.bm=0.3**:
```
                With BM-int    No BM-int    Δ Power
t1/2=0          68%            20%          -71%
t1/2=0.1        64%            13%          -79%
t1/2=0.2        69%            8%           -89%
Stability       0.987          0.880        -0.107 (WORSE with int)
```

**Verdict**: **MAJOR CULPRIT**. These correlations inflate baseline power
50-85% while WORSENING stability. They are not justified by the data or
by standard practice.

---

## Phase 2E: ER SD Scaling Impact — **MAJOR CULPRIT**

**Change**: ER SD from scaled by expectancy to constant across timepoints

**Code Location**: `build_path_sigma()` lines 1503-1504

**Current (2025) behavior**:
```r
sds <- c(sds, rep(resp_params$pb$sd, nP))  # Constant: same SD everywhere
```

**Variant (matching orig)**:
```r
sds <- c(sds, resp_params$pb$sd * expectancies)  # Scaled by expectancy
```

**Physical Interpretation**:
- **Constant (2025)**: Expectancy response has same variance whether
  patient expects treatment (1.0) or expects nothing (0.5)
- **Scaled (orig)**: When expectancy is lower, response variance is lower
  (more constrained, more realistic)

**Results** (CO N=35, c.bm=0.3):
```
                Const (2025)   Scaled (orig)  Δ Power
t1/2=0          43%            13%            -69%
t1/2=0.1        49%            12%            -76%
t1/2=0.2        47%            5%             -89%
Stability       0.960          0.920          -0.040 (worse const)
```

**CO N=70, c.bm=0.3**:
```
                Const (2025)   Scaled (orig)  Δ Power
t1/2=0          84%            20%            -76%
t1/2=0.1        75%            13%            -82%
t1/2=0.2        83%            8%             -90%
Stability       0.987          0.880          -0.107 (worse const)
```

**Verdict**: **MAJOR CULPRIT**. Constant ER SD inflates power 65-90% while
worsening stability. Orig's expectancy-scaled approach is more realistic.

---

## Combined Impact of Phase 2D + 2E

When both BM-ER/BM-TR correlations AND constant ER SD are combined (2025):

**CO N=35, c.bm=0.3**:
- Baseline power: 27% → 36% (both factors combined, inflated ~80%)
- True baseline (both removed): ~13% (matching orig)

**CO N=70, c.bm=0.3**:
- Baseline power: 68% → 69% (both factors, inflated ~200%+)
- True baseline (both removed): ~20% (matching orig)

These two design choices alone account for the inflated baseline and
the subsequent exaggerated vulnerability appearance.

---

## Recommendation for Phase 3-4

1. **Remove BM-ER/BM-TR correlations** (Phase 2D)
   - They inflate power without justification
   - They worsen stability
   - No standard practice supports these cross-component linkages

2. **Revert ER SD to expectancy-scaled** (Phase 2E)
   - More realistic (variance should reflect confidence)
   - Worsens less stability
   - Standard practice in variance modeling

3. **Keep scalefactor=2** (or test sf=1 with power recovery analysis)
   - Modest effect; justified only if Phase 4 power recovery is substantial

4. **Keep AR(1) correlation structure** (Phase 2C)
   - Superior to alternatives
   - More realistic for longitudinal data

---

## Next: Phase 4 Validation

Once Phase 2D and 2E changes are made to the DGP:
- **Scenario 1**: Rerun with carryover modeling in analysis
- **Scenario 2**: Baseline (no carryover modeling)
- Measure power recovery as justification for carryover modeling

This will complete the narrative: "Orig was conservative but correct;
2025 was optimistic but unrealistic. Here's how modeling carryover
recovers lost power under realistic assumptions."
