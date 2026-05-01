# Visualization Guide

## Overview

Two visualization scripts are available to display the Monte Carlo simulation results:

1. **`visualize_hendrickson_style.R`** - Original 2-panel layout
2. **`visualize_heatmaps.R`** - NEW 2×2 heatmap layout (recommended)

---

## NEW: 2×2 Heatmap Visualization (Recommended)

### Script: `visualize_heatmaps.R`

Creates a 2×2 grid of heatmaps showing statistical power across all parameter combinations.

### Layout

```
┌─────────────────────────────────┬─────────────────────────────────┐
│  Hybrid Design                  │  Hybrid Design                  │
│  WITHOUT Carryover Model        │  WITH Carryover Model           │
├─────────────────────────────────┼─────────────────────────────────┤
│  Crossover Design               │  Crossover Design               │
│  WITHOUT Carryover Model        │  WITH Carryover Model           │
└─────────────────────────────────┴─────────────────────────────────┘
```

### Each Heatmap Shows:
- **X-axis**: Carryover half-life (0, 1.0, 2.0 weeks)
- **Y-axis**: Biomarker correlation (0.2, 0.34)
- **Color**: Statistical power (0-1 scale)
  - Red (low power) → Yellow (medium) → Green (high power)
- **Numbers**: Exact power values

### Output Files:
- `../output/power_heatmaps_2x2.png` (300 DPI, 16×12 inches)
- `../output/power_heatmaps_2x2.pdf` (vector, publication-ready)

### Usage:
```bash
Rscript visualize_heatmaps.R
```

### Key Features:
✅ Shows all 12 parameter combinations simultaneously
✅ Easy comparison across designs (Hybrid vs Crossover)
✅ Easy comparison across analysis approaches (WITH vs WITHOUT carryover model)
✅ Clear visualization of carryover effect gradient
✅ Clear visualization of biomarker correlation effect

---

## Original: 2-Panel Hendrickson-Style Visualization

### Script: `visualize_hendrickson_style.R`

Creates a 2-panel figure matching Hendrickson et al. (2020) layout.

### Layout

**Panel A**: No carryover (t½ = 0)
- Varies biomarker correlation (0.2, 0.34)
- Compares Hybrid vs Crossover designs
- Split by WITH/WITHOUT carryover model

**Panel B**: With carryover (fixed c.bm = 0.34)
- Varies carryover half-life (1.0, 2.0 weeks)
- Compares Hybrid vs Crossover designs
- Split by WITH/WITHOUT carryover model

### Output Files:
- `../output/figure4_equivalent_hendrickson_style.png` (300 DPI, 14×10 inches)
- `../output/figure4_equivalent_hendrickson_style.pdf` (vector, publication-ready)

### Usage:
```bash
Rscript visualize_hendrickson_style.R
```

---

## Comparison

| Feature | Heatmaps (NEW) | Hendrickson-Style |
|---------|----------------|-------------------|
| Shows all conditions | ✅ Yes (12) | ⚠️ Partial (8 of 12) |
| Design comparison | ✅ Easy | ✅ Easy |
| Carryover gradient | ✅ Clear | ⚠️ Split across panels |
| Biomarker gradient | ✅ Clear | ⚠️ Split across panels |
| Analysis approach | ✅ Side-by-side | ✅ Side-by-side |
| Publication match | ❌ Novel layout | ✅ Matches Hendrickson |
| Information density | ✅ High | ⚠️ Medium |

---

## Recommendations

### For Manuscript:
- **Use both visualizations**
  - Heatmaps (2×2) as main figure - shows complete picture
  - Hendrickson-style as supplementary - matches original paper format

### For Presentations:
- **Use heatmaps (2×2)** - clearer, more comprehensive

### For Quick Reference:
- **Use heatmaps (2×2)** - all results at a glance

---

## Simulation Parameters (Fixed)

Both visualizations display results from the same Monte Carlo simulation:

- **Sample size**: N = 70 participants
- **Autocorrelations**: c.tv = c.pb = c.br = 0.8 (Hendrickson values)
- **Cross-correlations**: c.cf1t = 0.2, c.cfct = 0.1 (Hendrickson values)
- **Treatment effect**: 5.0
- **Iterations**: 20 per combination

### Varied Parameters:
- **Biomarker correlation (c.bm)**: 0.2, 0.34
- **Carryover half-life (weeks)**: 0, 1.0, 2.0
- **Design**: Hybrid (4-path), Crossover (2-sequence)
- **Analysis approach**: WITH carryover model, WITHOUT carryover model

**Total combinations**: 2 × 3 × 2 × 2 = 24 conditions (480 simulations)

---

## Interpreting Results

### Power Scale:
- **0.00-0.20**: Very low power (red)
- **0.20-0.50**: Low power (orange)
- **0.50-0.80**: Medium power (yellow)
- **0.80-0.95**: Good power (light green)
- **0.95-1.00**: Excellent power (dark green)

### Key Patterns to Look For:

1. **Design Effect**:
   - Compare top row (Hybrid) vs bottom row (Crossover)
   - Hybrid generally shows higher or equal power

2. **Analysis Approach Effect**:
   - Compare left column (WITHOUT) vs right column (WITH)
   - Shows impact of modeling carryover

3. **Carryover Effect**:
   - Move along x-axis (0 → 1.0 → 2.0)
   - Shows how carryover strength affects power

4. **Biomarker Correlation Effect**:
   - Move along y-axis (0.2 → 0.34)
   - Higher correlation should increase interaction effect size

---

## Current Results Summary

### Main Findings:
1. **Overall power is low** (0-0.20 range) for current parameters
2. **Hybrid design** slightly outperforms Crossover (mean: 0.117 vs 0.079)
3. **Biomarker correlation c.bm = 0.2** shows higher power than c.bm = 0.34
   - This is unexpected and may warrant investigation
4. **Carryover modeling** shows minimal impact on power
   - Power difference (WITH - WITHOUT) ranges from -0.05 to 0.00
5. **Crossover design at c.bm = 0.34** shows 0% power across all conditions
   - May indicate issue with parameter selection or need for larger N

### PD Boundary:
- **Maximum c.bm = 0.34** for autocorr = 0.8 (Hendrickson values)
- c.bm = 0.35 produces non-positive definite matrices
- Boundary independent of carryover level
- Boundary identical for both designs

---

## Troubleshooting

### Visualization script fails with "cannot open RData file"
**Cause**: Simulation hasn't been run yet

**Solution**:
```bash
# Run simulation first
Rscript full_pmsim_analysis_hyb_versus_co.R
# Then run visualization
Rscript visualize_heatmaps.R
```

### Heatmap shows unexpected patterns
**Check**:
1. Verify simulation completed successfully (all 480 runs)
2. Check for convergence warnings in simulation output
3. Review parameter grid in `full_pmsim_analysis_hyb_versus_co.R`
4. Check for errors in `simulation_results$error` column

---

## Future Enhancements

### To Increase Power:
1. Increase sample size (N = 70 → 100 or 150)
2. Increase iterations (20 → 100 for more stable estimates)
3. Increase treatment effect (5.0 → 7.0 or 10.0)
4. Test different biomarker correlation values
5. Consider random slope for time effect

### Additional Visualizations:
1. Power curves as function of sample size
2. Effect size distributions
3. Convergence diagnostics
4. Sensitivity analysis plots

---

*Last updated: 2025-11-18*
