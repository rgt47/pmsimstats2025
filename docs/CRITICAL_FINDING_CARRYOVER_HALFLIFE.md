# CRITICAL FINDING: Carryover Half-Life Scale Difference

## Date: 2025-11-18

## Discovery

Found major discrepancy between Hendrickson's parameters and ours:

### Hendrickson's Carryover Half-Lives
**Source**: `pmsim-orig/vignettes/Produce_Publication_Results_1_generate_data.Rmd:318`

```r
carryover_t1half=c(0,.1,.2)
```

- **0 weeks** = No carryover
- **0.1 weeks** = 0.7 days (~17 hours)
- **0.2 weeks** = 1.4 days (~33 hours)

### Our Carryover Half-Lives
**Source**: `full_pmsim_analysis_hyb_versus_co.R:530`

```r
carryover_t1half = c(0, 1.0, 2.0)
```

- **0 weeks** = No carryover
- **1.0 week** = 7 days
- **2.0 weeks** = 14 days

## Impact

### Scale Difference: **5-10x longer**

Our carryover persists **5-10 times longer** than Hendrickson's!

### Why This Matters

1. **Different mechanism**: Very short half-lives (hours/days) vs long half-lives (weeks)

2. **Alternating design interaction**: In A-B-A-B design with 1-week periods:
   - Hendrickson's t½=0.1 weeks: Carryover decays to 0.1% after 1 week
   - Our t½=1.0 week: Carryover decays to 50% after 1 week
   - Our t½=2.0 weeks: Carryover decays to 70% after 1 week

3. **Cumulative effects**: With long half-lives, carryover accumulates over multiple periods, creating a completely different pattern

4. **Treatment contrast**: Short half-lives preserve ON/OFF contrast; long half-lives blur it

## Visualization of Difference

```
Hendrickson (t½=0.2 weeks = 1.4 days):
Week  | Treatment | BR Mean | Carryover | Total
------|-----------|---------|-----------|-------
  1   |    ON     |  10.0   |    0.0    | 10.0
  2   |    OFF    |   0.0   |    0.2    |  0.2  ← 2% of previous
  3   |    ON     |  10.0   |    0.0    | 10.0
  4   |    OFF    |   0.0   |    0.2    |  0.2

Our implementation (t½=2.0 weeks):
Week  | Treatment | BR Mean | Carryover | Total
------|-----------|---------|-----------|-------
  1   |    ON     |  10.0   |    0.0    | 10.0
  2   |    OFF    |   0.0   |    7.1    |  7.1  ← 71% of previous!
  3   |    ON     |  10.0   |    0.0    | 10.0
  4   |    OFF    |   0.0   |    7.1    |  7.1
```

## Why This Might Explain Our Results

### Hypothesis: Long Half-Lives Reduce Treatment Contrast

1. **With short carryover (Hendrickson)**:
   - Clear ON/OFF distinction maintained
   - Carryover adds small noise
   - Interaction still detectable

2. **With long carryover (us)**:
   - OFF periods have high residual BR
   - Treatment contrast reduced
   - May reduce power in a DIFFERENT way than Hendrickson observed

3. **Time effect interaction**:
   - Very long carryover might create patterns that the `week` covariate can partially absorb
   - Short carryover creates person-specific noise that can't be absorbed

## Hendrickson's Figure 4 Context

**Looking at the figure**: The carryover values shown are labeled as t½ in WEEKS:
- 0, 1, 2, 4 weeks

**BUT**: The vignette code uses 0, 0.1, 0.2 weeks

**Possibility 1**: The figure legend is wrong (should say "days" not "weeks")
**Possibility 2**: Different simulations used for the figure vs vignette
**Possibility 3**: The code in the vignette isn't what generated Figure 4

## Next Steps

1. **Check Hendrickson's visualize2.Rmd** - This might have the actual Figure 4 parameters

2. **Run our simulation with Hendrickson's scale**:
   ```r
   carryover_t1half = c(0, 0.1, 0.2)  # Match her vignette
   ```

3. **Also test intermediate values**:
   ```r
   carryover_t1half = c(0, 1.0, 2.0, 4.0)  # Match Figure 4 legend
   ```

4. **Compare power patterns** across all scales to understand the relationship

## Files to Check

1. `~/prj/c265/pmsimstats-master/pmsim-orig/vignettes/visualize2.Rmd`
2. `~/prj/c265/pmsimstats-master/pmsim-orig/vignettes/run_simulations.R`
3. Any saved simulation results that generated Figure 4

---

*Critical finding documented: 2025-11-18*
