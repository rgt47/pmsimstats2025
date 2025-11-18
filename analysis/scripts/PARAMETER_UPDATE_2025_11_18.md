# Parameter Update: Carryover Half-Life Scale Correction

## Date: 2025-11-18

## Changes Made

### Carryover Half-Life Parameters

**File**: `full_pmsim_analysis_hyb_versus_co.R` (line 534)

**Previous values:**
```r
carryover_t1half = c(0, 1.0, 2.0)  # weeks
```

**Updated values:**
```r
carryover_t1half = c(0, 0.1, 0.2)  # weeks
```

**Scale change**: Reduced by 5-10x to match Hendrickson et al. (2020)

---

## Rationale

### Problem Identified

Through detailed code comparison with Hendrickson's original implementation (`pmsim-orig`), discovered that we were using carryover half-lives that were **5-10 times longer** than Hendrickson's:

| Condition | Our Previous | Hendrickson | Ratio |
|-----------|--------------|-------------|-------|
| None | 0 weeks | 0 weeks | 1x |
| Moderate | 1.0 week (7 days) | 0.1 weeks (0.7 days) | 10x |
| Strong | 2.0 weeks (14 days) | 0.2 weeks (1.4 days) | 10x |

### Impact of Long Half-Lives

With 1-2 week half-lives in a 1-week alternating design (A-B-A-B):

**Example with t½=2.0 weeks:**
```
Week 1 (ON):  BR = 10.0
Week 2 (OFF): BR = 0 + carryover(10.0 × 0.71) = 7.1  ← 71% persists!
Week 3 (ON):  BR = 10.0
Week 4 (OFF): BR = 0 + carryover(10.0 × 0.71) = 7.1
```
**Result**: Treatment contrast nearly eliminated

**Hendrickson with t½=0.2 weeks:**
```
Week 1 (ON):  BR = 10.0
Week 2 (OFF): BR = 0 + carryover(10.0 × 0.02) = 0.2  ← Only 2% persists
Week 3 (ON):  BR = 10.0
Week 4 (OFF): BR = 0 + carryover(10.0 × 0.02) = 0.2
```
**Result**: Clear treatment contrast maintained, carryover adds controlled noise

### Why This Caused Our Results to Differ

1. **Treatment contrast loss**: Long carryover blurred ON/OFF distinction
2. **Different mechanism**: Very long carryover creates a different pattern than short-term noise
3. **Time effect interaction**: Systematic long-term effects might be partially absorbed by `week` covariate
4. **No power decrease pattern**: Because the mechanism was fundamentally different from Hendrickson's

---

## Verification

### Code Comparison Confirmed

Reviewed Hendrickson's original code:
- `pmsim-orig/R/generateData.R` (data generation)
- `pmsim-orig/R/lme_analysis.R` (analysis model)
- `pmsim-orig/vignettes/Produce_Publication_Results_1_generate_data.Rmd` (parameter specification)
- `pmsim-orig/vignettes/visualize2.Rmd` (figure generation)

**Confirmed**:
- ✓ Implementation logic identical
- ✓ Analysis models equivalent
- ✓ Only difference is parameter scale

### Evidence

**From `Produce_Publication_Results_1_generate_data.Rmd:318`:**
```r
coremodelparams<-expand.grid(
  N=c(35,70),
  c.bm=c(0,.3,.6),
  carryover_t1half=c(0,.1,.2),  # ← Confirmed: 0.1 and 0.2 weeks
  c.tv=.8,c.pb=.8,c.br=.8,
  c.cf1t=.2,c.cfct=.1
)
```

**From `visualize2.Rmd:12`:**
> "different carryover half-life values (0, 0.1, 0.2 weeks)"

---

## Expected Impact of Change

### 1. Type I Error Rate

**Before** (20 iterations, c.bm=0):
- Observed: 0% (0/20 significant)
- Expected: ~5% (1/20)
- Issue: Too few iterations for stable estimate

**After** (100 iterations, c.bm=0):
- Expected: ~5% (5/100 significant)
- Should now match theoretical expectation

### 2. Power Pattern Across Carryover Levels

**Before** (t½ = 0, 1.0, 2.0 weeks):
- No systematic decrease
- Power varied randomly
- Pattern: 15% → 20% → 15%

**After** (t½ = 0, 0.1, 0.2 weeks):
- **Expected**: Power should DECREASE as carryover increases
- Should match Hendrickson's Figure 4 pattern
- Unmodeled carryover adds noise → reduces power

### 3. Hendrickson Alignment

All parameters now aligned with Hendrickson et al. (2020):

| Parameter | Our Value | Hendrickson | Status |
|-----------|-----------|-------------|--------|
| Autocorrelations | 0.6 | 0.8 | Different (for PD) |
| c.bm levels | 0, 0.3, 0.48 | 0, 0.3, 0.6 | Close (0.48 vs 0.6) |
| Carryover t½ | 0, 0.1, 0.2 | 0, 0.1, 0.2 | ✓ **Aligned** |
| N participants | 70 | 35, 70 | Partial |
| Time effect | Yes | Yes | ✓ Aligned |
| Analysis model | lmer | lmer | ✓ Aligned |

---

## Simulation Grid

**Updated grid** (line 523-534):
- 1 sample size: N = 70
- 3 effect sizes: c.bm = {0, 0.3, 0.48}
- 3 carryover levels: t½ = {0, 0.1, 0.2} weeks
- 2 designs: hybrid (4-path), crossover (2-sequence)
- 100 iterations per condition

**Total**: 3 × 3 × 2 = 18 conditions × 100 iterations = **1,800 simulation runs**

Each condition produces:
- WITH carryover model
- WITHOUT carryover model (testing unmodeled carryover impact)

**Total model fits**: 1,800 × 2 = **3,600 mixed models**

---

## Files Modified

1. **full_pmsim_analysis_hyb_versus_co.R** (line 534)
   - Changed carryover_t1half from c(0, 1.0, 2.0) to c(0, 0.1, 0.2)
   - Added detailed comment explaining the change

---

## Documentation Created

1. **HENDRICKSON_CODE_COMPARISON.md** - Detailed comparison of implementations
2. **CRITICAL_FINDING_CARRYOVER_HALFLIFE.md** - Analysis of scale discrepancy
3. **PARAMETER_UPDATE_2025_11_18.md** - This document

---

## Next Steps

1. **Run updated simulation**:
   ```bash
   cd /Users/zenn/Dropbox/prj/d08/pmsimstats2025/analysis/scripts
   Rscript full_pmsim_analysis_hyb_versus_co.R
   ```

2. **Verify expected patterns**:
   - Type I error ~5% at c.bm=0
   - Power decreases as carryover increases (when unmodeled)
   - Power higher for hybrid vs crossover design

3. **Generate visualizations** with corrected scale:
   ```bash
   Rscript visualize_heatmaps.R
   ```

4. **Update figures** showing carryover on correct scale (0, 0.1, 0.2 weeks)

---

*Parameter update completed: 2025-11-18*
