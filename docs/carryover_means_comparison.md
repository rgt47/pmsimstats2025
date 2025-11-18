# Comparison: Hendrickson vs Your Code - Carryover Adjustments to Means

## Summary

**NO, they are NOT the same.** Your code has several significant differences from Hendrickson's approach.

---

## Hendrickson's Approach

### Formula (generateData.R, lines 82-93)

```r
if(c=="br"){
  brmeans <- modgompertz(d$tod, rp$max, rp$disp, rp$rate)
  if(nP>1){
    for(p in 2:nP){
      if(!d[p]$onDrug){
        if(d[p]$tsd>0){
          brmeans[p] <- brmeans[p] + brmeans[p-1] * (1/2)^(d$tsd[p]/modelparam$carryover_t1half)
        }
      }
    }
  }
  means <- c(means, brmeans)
}
```

### Key Characteristics

1. **Components affected**: ONLY `br` (bio-response)
   - `tv` (time-variant): NO carryover
   - `pb` (pharm-biomarker): NO carryover
   - `br` (bio-response): YES carryover

2. **Carryover formula**:
   ```
   μ[t] = base[t] + μ[t-1] * (1/2)^(tsd[t] / t1half)
   ```
   - No scale factor
   - Uses `tsd` (time since discontinuation) directly

3. **Application logic**:
   - Applied ONLY when `!onDrug` (off treatment)
   - Applied ONLY when `tsd > 0` (after discontinuation)
   - Iterates sequentially through timepoints: `for(p in 2:nP)`

4. **Effect accumulation**:
   - Carries over from IMMEDIATELY PREVIOUS timepoint (`p-1`)
   - Creates a sequential dependency chain

---

## Your Code's Approach

### Formula (pm_functions.R, lines 209-271)

```r
apply_carryover_to_component <- function(
    component_means, trial_data, component_halflives,
    scale_factor, component_name) {

  # Get component-specific half-life
  component_halflife <- component_halflives[[paste0(component_name, "_halflife")]]

  # Find timepoints for carryover application
  if (component_name == "br") {
    carryover_indices <- which(!trial_data$on_drug & trial_data$tsd > 0)
  } else {
    carryover_indices <- 2:num_timepoints  # ALL timepoints!
  }

  for (idx in carryover_indices) {
    prev_idx <- idx - 1

    if (component_name == "br") {
      time_lag <- trial_data$tsd[idx]
    } else {
      time_lag <- 1  # Fixed 1-week lag
    }

    decay_factor <- (1/2)^(scale_factor * time_lag / component_halflife)
    component_means[idx] <- component_means[idx] + component_means[prev_idx] * decay_factor
  }
}
```

### Key Characteristics

1. **Components affected**: ALL THREE components
   - `tv`: YES carryover (with `tv_halflife = base * 1.5`)
   - `pb`: YES carryover (with `pb_halflife = base * 0.6`)
   - `br`: YES carryover (with `br_halflife = base`)

2. **Carryover formula**:
   ```
   μ[t] = base[t] + μ[t-1] * (1/2)^(scale_factor * time_lag / t1half)
   ```
   - Includes `scale_factor` parameter (typically 1.5)
   - Component-specific half-lives

3. **Application logic**:
   - **For `br`**: Same as Hendrickson (only when off drug, tsd > 0)
   - **For `tv` and `pb`**: Applied to ALL timepoints after the first (lines 234-238)

4. **Effect accumulation**:
   - Same as Hendrickson: carries from immediately previous timepoint

---

## Detailed Differences

### Difference 1: Component Coverage

| Component | Hendrickson | Your Code |
|-----------|-------------|-----------|
| tv (time-variant) | ❌ No carryover | ✅ **Carryover at all timepoints** |
| pb (pharm-biomarker) | ❌ No carryover | ✅ **Carryover at all timepoints** |
| br (bio-response) | ✅ Carryover when off drug | ✅ Carryover when off drug |

**Impact**: Your code adds carryover effects to components that Hendrickson did NOT model with carryover.

### Difference 2: Scale Factor

**Hendrickson**:
```r
decay = (1/2)^(tsd / t1half)
```
Example: `tsd = 2`, `t1half = 1`
→ `decay = (1/2)^2 = 0.25`

**Your Code**:
```r
decay = (1/2)^(scale_factor * tsd / t1half)
```
Example: `tsd = 2`, `t1half = 1`, `scale_factor = 1.5`
→ `decay = (1/2)^3 = 0.125`

**Impact**: Your carryover decays **faster** due to scale factor > 1.

### Difference 3: Component-Specific Half-Lives

**Hendrickson**: Single `carryover_t1half` for all effects

**Your Code** (pm_functions.R, lines 143-152):
```r
calculate_component_halflives <- function(base_halflife) {
  list(
    br_halflife = base_halflife,          # Baseline
    pb_halflife = base_halflife * 0.6,    # SHORTER (60%)
    tv_halflife = base_halflife * 1.5     # LONGER (150%)
  )
}
```

**Impact**: Different components have different persistence rates.

### Difference 4: TV and PB Carryover Logic

**Hendrickson**: N/A (no carryover for these components)

**Your Code** (lines 231-238):
```r
# For tv and pb: carryover can occur any time there's a
# previous effect. Simplified: apply to timepoints 2
# onwards if we have multiple timepoints
carryover_indices <- if (num_timepoints > 1) {
  2:num_timepoints
} else {
  integer(0)
}
```

**Means**: TV and PB get carryover at EVERY timepoint, regardless of treatment status!

---

## Concrete Example

### Setup
- Design: 4 timepoints
- Weeks 1-2: On treatment
- Weeks 3-4: Off treatment
- `carryover_t1half = 1`
- `scale_factor = 1.5` (your code)
- Initial response values: `[5, 5, 0, 0]` (gompertz outputs)

### Hendrickson's BR Calculation

```
Week 1: br[1] = 5                           (base)
Week 2: br[2] = 5                           (base, still on drug)
Week 3: br[3] = 0 + 5 * (1/2)^(1/1)        (off drug, tsd=1)
              = 0 + 2.5 = 2.5
Week 4: br[4] = 0 + 2.5 * (1/2)^(2/1)      (off drug, tsd=2)
              = 0 + 0.625 = 0.625

TV and PB: NO CARRYOVER, just gompertz outputs
```

### Your Code's BR Calculation

```
Week 1: br[1] = 5                                    (base)
Week 2: br[2] = 5                                    (base, still on drug)
Week 3: br[3] = 0 + 5 * (1/2)^(1.5*1/1)             (off drug, tsd=1)
              = 0 + 5 * (1/2)^1.5
              = 0 + 1.768 = 1.768
Week 4: br[4] = 0 + 1.768 * (1/2)^(1.5*2/1)         (off drug, tsd=2)
              = 0 + 1.768 * (1/2)^3
              = 0 + 0.221 = 0.221
```

**Difference**: Your carryover decays faster (0.221 vs 0.625 at week 4)

### Your Code's TV Calculation

```
Week 1: tv[1] = 3.2                                  (gompertz base)
Week 2: tv[2] = 4.1 + 3.2 * (1/2)^(1.5*1/1.5)       (carryover applied!)
              = 4.1 + 3.2 * (1/2)^1
              = 4.1 + 1.6 = 5.7
Week 3: tv[3] = 4.5 + 5.7 * (1/2)^(1.5*1/1.5)
              = 4.5 + 2.85 = 7.35
Week 4: tv[4] = 4.7 + 7.35 * (1/2)^(1.5*1/1.5)
              = 4.7 + 3.675 = 8.375
```

**This is COMPLETELY NEW** - Hendrickson has no carryover for TV!

### Your Code's PB Calculation

Similar to TV, but with `pb_halflife = 1 * 0.6 = 0.6`:
```
Week 2: pb[2] = base + pb[1] * (1/2)^(1.5*1/0.6)
              = base + pb[1] * (1/2)^2.5
              = base + pb[1] * 0.177
```

**Also COMPLETELY NEW** - Hendrickson has no carryover for PB!

---

## Mathematical Comparison Table

| Aspect | Hendrickson | Your Code |
|--------|-------------|-----------|
| **Components with carryover** | 1 (br only) | 3 (tv, pb, br) |
| **BR formula** | `μ[t] + μ[t-1] * (1/2)^(tsd/t½)` | `μ[t] + μ[t-1] * (1/2)^(SF*tsd/t½)` |
| **TV formula** | None | `μ[t] + μ[t-1] * (1/2)^(SF*1/t½_tv)` |
| **PB formula** | None | `μ[t] + μ[t-1] * (1/2)^(SF*1/t½_pb)` |
| **Scale factor** | 1.0 (implicit) | 1.5 (default) |
| **Component half-lives** | Uniform | Heterogeneous (0.6×, 1×, 1.5×) |
| **TV carryover condition** | N/A | Always (all timepoints) |
| **PB carryover condition** | N/A | Always (all timepoints) |
| **BR carryover condition** | Off drug & tsd > 0 | Off drug & tsd > 0 |

---

## Theoretical Implications

### Hendrickson's Model

**Interpretation**:
- Time-variant factor (TV): Natural disease progression, no treatment persistence
- Pharm-biomarker factor (PB): Expectancy effect, dissipates immediately
- Bio-response factor (BR): True drug effect, persists after discontinuation

**Rationale**: Only the biological drug effect has meaningful carryover.

### Your Model

**Interpretation**:
- TV: Has "memory" - disease trajectory influenced by past states
- PB: Expectancy has "momentum" - psychological effects persist
- BR: Drug effects persist (like Hendrickson, but with faster decay)

**Rationale**: All components have temporal dependencies, modeled as autoregressive processes with different time constants.

---

## Which is More Realistic?

### Arguments for Hendrickson (Carryover Only in BR)

1. **Pharmacological basis**: Only drug effects should have washout periods
2. **Conceptual clarity**: Separates persistent drug effects from other time trends
3. **Parsimony**: Fewer parameters to estimate/specify
4. **Empirical validation**: Published and peer-reviewed

### Arguments for Your Approach (Carryover in All Components)

1. **Biological realism**: Disease progression has inertia (TV carryover)
2. **Psychological realism**: Expectancy effects don't vanish instantly (PB carryover)
3. **Flexibility**: Component-specific half-lives capture different mechanisms
4. **Statistical modeling**: Aligns with autoregressive time series thinking

However, you need to be careful about:
- **Identifiability**: Can you distinguish TV carryover from BR carryover?
- **Parameter inflation**: 3 half-lives + scale factor vs 1 half-life
- **Comparison to Hendrickson**: Your results won't be directly comparable

---

## Recommendations

### If Goal is to Replicate/Extend Hendrickson:

**Remove TV and PB carryover**:

```r
apply_carryover_to_component <- function(...) {
  # ONLY apply to BR component
  if (component_name != "br") {
    return(component_means)  # No carryover for tv or pb
  }

  # For BR: same logic as current
  # But consider removing scale_factor to match Hendrickson exactly
}
```

**Set scale_factor = 1** to match Hendrickson's formula exactly.

### If Goal is Novel Enhancement:

**Keep your approach BUT**:
1. Justify theoretically why TV and PB should have carryover
2. Consider identifiability issues
3. Document clearly as an extension of Hendrickson
4. Run ablation study comparing:
   - No carryover (null model)
   - BR carryover only (Hendrickson)
   - All three components (your model)

---

## Bottom Line

**Are the adjustments the same?**

### For BR Component:
**Similar but not identical**:
- Same application logic (off drug + tsd > 0)
- Different decay rate (your scale factor makes it faster)
- Different half-lives possible (component-specific vs uniform)

### For TV and PB Components:
**Completely different**:
- Hendrickson: No carryover
- Your code: Carryover at all timepoints with component-specific half-lives

**Overall**: Your implementation is a **substantial enhancement** of Hendrickson's model, not a replication. Whether this is better depends on your theoretical goals and whether the added complexity is justified by the research question.
