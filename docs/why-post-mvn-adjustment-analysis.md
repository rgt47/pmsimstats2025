# Why Are We Doing Post-MVN Adjustment? Deep Analysis

## The Core Question

Why does our implementation include this step:

```r
# Post-MVN adjustment
participant_data[[br_col]] <- participant_data[[br_col]] +
  (participant_data$biomarker * treatment_status * interaction_strength)
```

Is it necessary? Is it theoretically justified? Should we keep it?

---

## Two Fundamentally Different Mechanisms for Biomarker×Treatment Interaction

### Mechanism A: Correlation-Based (Hendrickson's Approach)

**The idea**: Biomarker and response are correlated because they share common underlying biology. The correlation varies by treatment status.

**Implementation**:
```
1. Define: Cor(biomarker, bio_response) = c.bm
2. But make this correlation DIFFERENT by treatment period:
   - On treatment: Cor(bm, br) = c.bm
   - Off treatment: Cor(bm, br) = 0 (or scaled down)
3. Draw from MVN
4. Done
```

**Mathematical Result**:
- Participants with high biomarker get higher responses ON treatment
- But this happens through CORRELATED DRAWS from MVN
- Effect size limited by correlation magnitude (can't exceed |1|)
- Conditional expectation: `E[br|bm] = μ_br + c.bm*(σ_br/σ_bm)*(bm - μ_bm)`
- **Slope of biomarker effect**: `c.bm*(σ_br/σ_bm)`

**Biological Interpretation**:
- Biomarker = inflammatory marker
- Treatment = anti-inflammatory drug
- Both reflect same underlying inflammatory state
- High biomarker → high baseline symptoms → more room to improve
- **Mechanism**: Shared biology/correlation, no direct causation

---

### Mechanism B: Direct Additive Effect (Our Post-MVN)

**The idea**: Biomarker DIRECTLY moderates treatment response through a mechanistic pathway (e.g., biomarker = enzyme activity that metabolizes drug).

**Implementation**:
```
1. Draw baseline biomarker and response from MVN (possibly uncorrelated)
2. THEN add: br_i = br_i + (bm_i × treatment × interaction_slope)
3. This is a direct, mechanistic effect
```

**Mathematical Result**:
- Participants with high biomarker get ADDITIONAL response boost
- Happens through DIRECT ADDITION, not correlation
- Effect size NOT limited by correlation
- **Slope of biomarker effect**: exactly `interaction_slope` (which we set to c.bm)

**Biological Interpretation**:
- Biomarker = drug-metabolizing enzyme level
- Treatment = drug that requires metabolism
- High enzyme → more active drug metabolite → better response
- **Mechanism**: Direct causal pathway

---

## The Problem: We're Doing BOTH

Currently our code:
1. **Uses c.bm as correlation** in sigma matrix → Mechanism A
2. **Uses c.bm as slope** in post-MVN adjustment → Mechanism B
3. **Same parameter (c.bm) serves double duty!**

### Mathematical Effect

**With correlation only** (Hendrickson):
- Slope of biomarker effect: `c.bm*(σ_br/σ_bm)`

**With post-MVN addition** (us):
- Correlation contributes: `c.bm*(σ_br/σ_bm)` to the slope
- Post-MVN adds: `c.bm` to the slope
- **Total slope**: `c.bm*(σ_br/σ_bm) + c.bm`

If σ_br/σ_bm ≈ 1, we're getting **DOUBLE the intended effect**!

---

## Why This Matters: Detection Power

Let's trace what happens in the statistical model:

### With Correlation Only (Hendrickson)
```r
# Data has differential correlation by treatment
# Model: response ~ treatment * biomarker
# Tests if slope of biomarker differs by treatment
# Effect size ≈ c.bm*(σ_br/σ_bm)
```

### With Post-MVN Addition (Us)
```r
# Data has EXPLICIT linear term: br += bm*treatment*c.bm
# Model: response ~ treatment * biomarker
# Tests same thing but with LARGER effect
# Effect size ≈ c.bm*(σ_br/σ_bm) + c.bm
```

**Result**: Our approach creates STRONGER interactions → EASIER to detect → HIGHER power

But is this the right comparison? Are we testing the same biological scenario as Hendrickson?

---

## What Did Hendrickson Intend c.bm to Mean?

From her documentation:
```r
#' \item{\code{c.bm}}{  Correlation between the biomarker and the biologic response}
```

It's explicitly a **CORRELATION**, not a regression slope!

### Correlation Properties
- Range: [-1, 1]
- Dimensionless
- Symmetric

### Our Post-MVN Treats It As a Slope
- Range: unbounded
- Has units (response units per biomarker unit)
- Directional

**We're mixing incompatible interpretations of the same parameter!**

---

## Why Was This Post-MVN Adjustment Added?

Looking at the initial commit description and code comments, there's no explicit justification for the post-MVN adjustment.

### Hypothesis 1: Thought correlation alone wouldn't create detectable interactions
- This is **FALSE** - correlation DOES create interactions (Hendrickson proves this)

### Hypothesis 2: Wanted stronger/more reliable interaction effects
- This is possible - makes power analysis less sensitive to correlation strength
- But inflates effect size beyond what the biological model (correlation) would produce

### Hypothesis 3: Implementing a different biological model than Hendrickson
- Direct pharmacogenetic effect rather than background correlation
- But then should use different parameters, not reuse c.bm

---

## The Biological Interpretation Question

What's the TRUE biological relationship between biomarker and treatment response?

### Scenario 1: Shared Biology (Correlation Model) ✓
- Biomarker and response both reflect underlying disease severity
- High biomarker → high baseline symptoms → more room to improve
- **Mechanism**: Correlation, no direct causation
- **Model**: Hendrickson's approach

### Scenario 2: Direct Moderation (Additive Model)
- Biomarker directly affects drug metabolism/efficacy
- High biomarker → mechanistically stronger drug effect
- **Mechanism**: Direct causal pathway
- **Model**: Post-MVN adjustment

**Both are biologically plausible, but they're DIFFERENT models!**

---

## Three Options Going Forward

### Option 1: Match Hendrickson (Pure Correlation) ⭐ RECOMMENDED

**Remove**:
- Population mean shift (`c.bm * tod * 2.0`)
- Post-MVN adjustment (`bm_i * treatment * c.bm`)

**Keep**:
- Correlation-based interaction only

**Advantages**:
- ✓ Conceptually clean
- ✓ Direct comparison with Hendrickson
- ✓ c.bm has clear meaning (correlation)
- ✓ Published, validated approach

**Disadvantages**:
- Effect size limited by correlation
- May need higher c.bm values for adequate power

**Rationale**:
1. **Scientific validity**: Hendrickson's correlation-based approach is theoretically sound and published
2. **Fair comparison**: We should use the same data-generating mechanism
3. **Parameter clarity**: c.bm as correlation has clear bounds and interpretation
4. **Sufficient effects**: Correlation DOES create detectable interactions (Hendrickson proves this)

---

### Option 2: Pure Additive Model

**Remove**:
- Population mean shift
- Biomarker-response correlation (or set to small background value)

**Keep**:
- Post-MVN adjustment
- **Rename c.bm** → `beta_interaction` or `interaction_slope`

**Advantages**:
- Explicit mechanistic model
- Effect size not limited
- Clear interpretation of parameter

**Disadvantages**:
- Different model than Hendrickson
- Not directly comparable
- Need to justify parameter choices

---

### Option 3: Hybrid (Both Mechanisms)

**Use TWO parameters**:
- `c.bm_correlation` = correlation component (e.g., 0.2)
- `beta_interaction` = direct effect component (e.g., 0.5)

**Remove**:
- Population mean shift (still redundant)

**Keep**:
- Both correlation and post-MVN adjustment
- Document that we're modeling BOTH mechanisms

**Advantages**:
- Most flexible
- Can model complex biology
- Can test relative contributions

**Disadvantages**:
- More complex
- Harder to interpret
- More parameters to justify

---

## Why Option 1 Is Recommended

### 1. The Post-MVN Adjustment Is Making the Problem EASIER Than It Should Be

By adding extra biomarker effect beyond what correlation alone would produce:
- We inflate the interaction signal
- Power estimates become **non-comparable** with Hendrickson's
- We're essentially "helping" the design by creating stronger signals
- Not testing the same biological scenario

### 2. Hendrickson's Approach Is Sufficient

Correlation-based interactions:
- ✓ Create detectable biomarker×treatment effects
- ✓ Have clear theoretical justification
- ✓ Match published methodology
- ✓ Don't require arbitrary post-hoc adjustments

### 3. Current Implementation Has Conceptual Issues

- Same parameter (c.bm) used for both correlation AND slope
- Incompatible interpretations mixed together
- Inflated means affect correlation scaling (mean-ratio logic)
- Potential contribution to matrix stability issues

### 4. Scientific Reproducibility

To claim alignment with Hendrickson et al. (2020):
- Must use same data-generating mechanism
- Correlation-based interaction is what they used
- Adding post-MVN adjustment changes the model fundamentally

---

## What Needs to Change

### Remove from `calculate_bio_response_with_interaction()`:

**Current code** (lines 113-120):
```r
# Create biomarker×treatment interaction effect
interaction_strength <- model_param$c.bm
biomarker_interaction_effect <- interaction_strength *
                                 trial_data$tod * 2.0

# Combine base response with biomarker interaction
bio_response_means <- base_bio_response +
                      biomarker_interaction_effect
```

**New code**:
```r
# No population mean shift needed
# Biomarker interaction emerges from correlation structure only
bio_response_means <- base_bio_response
```

---

### Remove from `generate_data()`:

**Current code** (lines 643-668):
```r
# Apply TRUE biomarker×treatment interaction to bio-response
# components. This creates individualized treatment responses
# based on each participant's biomarker level
interaction_strength <- model_param$c.bm

# Find bio-response columns
br_columns <- grep("\\.br$", names(participant_data), value = TRUE)

if (length(br_columns) > 0 && interaction_strength > 0) {
  # Get corresponding treatment status for each timepoint
  for (br_col in br_columns) {
    # Extract timepoint name (e.g., "W1" from "W1.br")
    timepoint_name <- sub("\\.br$", "", br_col)
    timepoint_idx <- which(trial_design$timepoint_name == timepoint_name)

    if (length(timepoint_idx) > 0) {
      # Get treatment status for this timepoint
      treatment_status <- trial_design$tod[timepoint_idx]

      # Apply biomarker×treatment interaction:
      # Individual biomarker value × treatment status × interaction strength
      participant_data[[br_col]] <- participant_data[[br_col]] +
        (participant_data$biomarker * treatment_status * interaction_strength)
    }
  }
}
```

**New code**:
```r
# No post-MVN adjustment needed
# Biomarker×treatment interaction emerges naturally from
# differential correlation structure in the MVN draw
# (Following Hendrickson et al. 2020 approach)
```

---

## Expected Impact of Changes

### On Data Generation
- Means will be lower (no `c.bm * tod * 2.0` inflation)
- Correlation scaling ratios will be smaller
- Sigma matrices may be more stable
- Biomarker effect emerges purely from correlation

### On Power
- May see **lower power** compared to current implementation
  - Because interaction signal is weaker (no double-counting)
  - This is CORRECT - current implementation artificially inflates power
- Power estimates will be **comparable to Hendrickson's**
- Will accurately reflect correlation-based interaction strength

### On Interpretation
- c.bm clearly means "correlation between biomarker and bio_response"
- No mixing of correlation and slope interpretations
- Clean conceptual model
- Direct comparison with published work

---

## Theoretical Justification for Correlation-Only Approach

### How Correlation Creates Interaction

When we have:
```
Cor(biomarker, bio_response_on_treatment) = c.bm
Cor(biomarker, bio_response_off_treatment) = 0
```

This creates a biomarker×treatment interaction because:

1. **On treatment**: High biomarker → high response (positive correlation)
2. **Off treatment**: Biomarker doesn't predict response (zero correlation)
3. **Interaction**: Treatment effect differs by biomarker level

### Mathematical Proof

From multivariate normal theory:
```
E[bio_response | biomarker, on_treatment] =
  μ_br + c.bm * (σ_br/σ_bm) * (biomarker - μ_bm)

E[bio_response | biomarker, off_treatment] = μ_br
```

The **difference** (treatment effect) is:
```
Treatment_effect(biomarker) = c.bm * (σ_br/σ_bm) * (biomarker - μ_bm)
```

This is a **linear function of biomarker** → a true interaction!

**No post-MVN adjustment needed** - the correlation structure ALREADY creates the interaction.

---

## Conclusion

The post-MVN adjustment is:
- ❌ **Unnecessary** for creating biomarker interactions
- ❌ **Inflates** effect size beyond correlation-based model
- ❌ **Mixes** incompatible interpretations of c.bm
- ❌ **Prevents** fair comparison with Hendrickson
- ❌ **Makes** the detection problem artificially easier

**Recommendation**: Remove both population mean shift AND post-MVN adjustment. Implement pure correlation-based approach (Option 1) to match Hendrickson's validated methodology.

---

*Analysis completed: 2025-11-18*
