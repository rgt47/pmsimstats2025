# Hendrickson 4-Path Randomization Implementation

## Summary

This document describes the implementation of Hendrickson et al. (2020)'s 4-path randomization structure for the Hybrid N-of-1 trial design in `full_pmsim_analysis_hyb_versus_co.R`.

## What Was Changed

### 1. Hybrid Design Now Has 4 Randomization Paths

**Previous Implementation:**
- All participants followed the same deterministic sequence
- Weeks 1-8: On treatment (open label)
- Weeks 9-12: Off treatment (blinded discontinuation)
- Weeks 13-16: On treatment (crossover period 1)
- Weeks 17-20: Off treatment (crossover period 2)

**New Implementation (Hendrickson's Structure):**
- Participants are randomized to 1 of 4 paths representing 2×2 factorial design:

| Path | BD Phase (wk 9-12) | CO Sequence (wk 13-20) | Description |
|------|-------------------|----------------------|-------------|
| A    | Stay on drug      | Drug first (AB)      | Continuous + AB |
| B    | Stay on drug      | Placebo first (BA)   | Continuous + BA |
| C    | Discontinue       | Drug first (AB)      | Discontinuation + AB |
| D    | Discontinue       | Placebo first (BA)   | Discontinuation + BA |

All paths share the same open-label phase (weeks 1-8: on treatment).

### 2. Key Function Updates

#### `create_designs()` (lines 297-397)
- Randomly assigns participants to 1 of 4 paths using balanced randomization
- Implements path-specific treatment schedules
- Adds `expectancy` variable (1.0 for open-label, 0.5 for blinded phases)
- Adds `path` variable to track randomization assignment

#### `prepare_design_for_simulation()` (lines 407-443)
- Updated to use `expectancy` column from design if present
- Falls back to treatment-based expectancy if not provided
- Maintains compatibility with both old and new design structures

#### `run_monte_carlo()` (lines 33-285)
- **Major restructuring** to handle multiple paths
- Generates data separately for each path (following Hendrickson approach)
- Each path gets its own simulation with appropriate sample size
- Combines all path data before analysis
- Uses actual design information (via left_join with full_design) instead of hardcoded treatment sequences
- Simplified carryover calculation (uniform across designs)

#### Design Storage Section (lines 588-633)
- Creates separate design templates for each path
- Hybrid: 4 path templates (one per randomization arm)
- Crossover: 2 path templates (AB and BA sequences)
- Stores both `full_design` (with all participants/paths) and `design_paths` (list of templates)

#### Sigma Caching Section (lines 679-695)
- Updated to use first path template for sigma matrix building
- All paths share same timepoint structure, so one sigma matrix per design suffices

### 3. Bug Fixes

- **Removed `browser()` statement** at line 450 (now line 524)
- Treatment assignment now uses actual design data instead of hardcoded logic
- Carryover calculation simplified and unified across designs

## How It Works

### Data Generation Process

1. **Design Creation**: `create_designs()` randomly assigns each participant to a path
2. **Path Template Creation**: Each unique path gets its own design template
3. **Simulation Loop**: For each iteration:
   - Loop through all paths (1-4 for hybrid, 1-2 for crossover)
   - Find participants assigned to current path
   - Generate data for those participants using the path's design template
   - Combine data from all paths
4. **Analysis**: Mixed model analyzes combined dataset with all paths

### Randomization Balance

- Uses `rep(1:4, length.out = n_participants)` then `sample()` to ensure approximately equal allocation
- With N=70: ~17-18 participants per path
- With N=35: 8-9 participants per path

## Comparison to Hendrickson's Code

### Matches Hendrickson:
✓ 4-path randomization structure for hybrid design
✓ Separate data generation per path
✓ 2×2 factorial design (BD × CO sequence)
✓ Combined analysis across all paths
✓ Expectancy modeling (1.0 for open-label, 0.5 for blinded)

### Differences from Hendrickson:
- **Enhanced carryover**: Component-specific half-lives (br/pb/tv)
- **Dynamic correlations**: Adjust based on carryover strength
- **Sigma caching**: Pre-validation for numerical robustness
- **Simplified path handling**: Uses expand_grid + case_when instead of manual path construction

## Testing Recommendations

1. **Verify path balance**:
```r
original_designs$hybrid %>%
  group_by(path) %>%
  summarize(n = n_distinct(participant_id))
```

2. **Check treatment sequences by path**:
```r
original_designs$hybrid %>%
  filter(week %in% c(8, 12, 16, 20)) %>%
  select(path, week, treatment) %>%
  distinct() %>%
  arrange(path, week)
```

3. **Verify expectancy**:
```r
original_designs$hybrid %>%
  select(week, expectancy) %>%
  distinct() %>%
  arrange(week)
```

## Expected Impact on Results

### Power May Decrease Slightly

With the 4-path randomization:
- Paths B & D start crossover phase with placebo (less treatment exposure)
- Paths C & D discontinue during BD phase (less cumulative treatment)
- More heterogeneity in treatment patterns across participants

This is **more realistic** but may show **lower power** than the previous deterministic design where everyone got optimal treatment exposure.

### Design Comparison More Fair

Both hybrid and crossover now properly randomize participants, making the power comparison between designs more valid and aligned with Hendrickson's methodology.

## Next Steps

1. Run simulation with small n_iterations (5-10) to verify it works
2. Check that path assignments are balanced
3. Verify treatment sequences match expected patterns
4. Run full simulation (increase n_iterations to 300+)
5. Compare power results to Hendrickson's published findings

## References

Hendrickson, N. M., et al. (2020). "Comparing Four N-of-1 Trial Designs for Predictive Biomarker Validation." *Journal Name*, *Volume*(Issue), pages.
