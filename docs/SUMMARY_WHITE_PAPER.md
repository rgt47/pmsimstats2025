# Statistical Power for Detecting Biomarker-Treatment Interactions in N-of-1 Trials: A Monte Carlo Simulation Study

**pmsimstats2025 Project Summary**

---

## Abstract

This document summarizes the pmsimstats2025 simulation framework, which evaluates
statistical power for detecting biomarker-treatment interactions across multiple
N-of-1 trial designs. Building on Hendrickson et al. (2020), the framework
extends the methodology to explicitly model carryover effects and quantify the
statistical cost of ignoring carryover in analysis. Through Monte Carlo
simulation comparing five trial designs (Hybrid, OL+BDC, Open-Label, Crossover,
and Parallel), the project demonstrates that explicit carryover modeling
preserves statistical power across varying carryover conditions, while omitting
carryover terms leads to progressive power degradation.

---

## 1. Introduction

### 1.1 Background

Precision medicine seeks to optimize treatment selection based on individual
patient characteristics. Biomarker-guided treatment represents a promising
approach, but detecting biomarker-treatment interactions requires adequate
statistical power. N-of-1 trials offer a rigorous framework for evaluating
individual treatment effects, yet design choices significantly influence the
ability to detect moderating effects.

### 1.2 Research Questions

1. How does statistical power to detect biomarker-treatment interactions compare
   across hybrid and traditional trial designs?
2. What is the impact of carryover effects on power estimation?
3. What is the statistical cost of ignoring carryover in the analysis model?

### 1.3 Theoretical Foundation

This work builds directly on Hendrickson et al. (2020), "Optimizing Aggregated
N-Of-1 Trial Designs for Predictive Biomarker Validation: Statistical Methods
and Theoretical Findings" (Frontiers in Digital Health, 2:13), implementing
their three-factor response model (BR + ER + TR) and extending it with explicit
carryover modeling in the analysis phase.

The Hendrickson framework addresses a specific question: whether a baseline
biomarker (standing systolic blood pressure) predicts individual response to
prazosin for PTSD. The simulation approach generates correlated multivariate
normal data representing biological, expectancy, and time-variant response
components, then tests for significant biomarker × treatment interactions.

---

## 2. Methods

### 2.1 Trial Designs Evaluated

The framework compares five distinct trial designs:

| Design      | Measurement Schedule (weeks)   | Points | Structure                        |
|-------------|--------------------------------|--------|----------------------------------|
| Hybrid      | 4, 8, 9, 10, 11, 12, 16, 20    | 8      | Dense cluster at phase transition|
| OL+BDC      | 4, 8, 12, 16, 17, 18, 19, 20   | 8      | Dense at discontinuation         |
| Open-Label  | 2, 8, 14, 20                   | 4      | Evenly spaced (6-week intervals) |
| Crossover   | 2, 8, 14, 20                   | 4      | Traditional AB/BA sequence       |
| Parallel    | 2, 8, 14, 20                   | 4      | Between-subjects comparison      |

### 2.2 Three-Factor Response Model

The simulation generates participant responses as the sum of three components:

```
Response = Baseline + BR + ER + TR
```

Where:

- **Biological Response (BR)**: Accumulates during active treatment, moderated
  by biomarker level. Includes carryover with exponential decay.
- **Expectancy Response (ER)**: Based on treatment visibility, scaled by
  participant knowledge of treatment status (1.0 for open-label when participant
  knows they are on active drug, 0.5 for blinded conditions where participant
  believes there is 50% chance of active drug).
- **Time-variant Response (TR)**: Non-treatment-specific improvement over trial
  duration, including regression to mean and effects of study participation
  (contact with staff, symptom monitoring). Per Hendrickson's "tabula rasa"
  assumption, TR and ER are assumed to have equal magnitude and variance in the
  absence of empirical data to distinguish them.

### 2.3 Biomarker Moderation

The biomarker-treatment interaction is modeled as:

```
BR_rate × (1 + β_mod × biomarker_centered)
```

Where `β_mod` represents the moderation strength, varied systematically from 0
to 0.65 across simulation conditions.

### 2.4 Carryover Modeling

Carryover effects follow exponential decay, applied only to the biological
response (BR) component. The framework supports a unified parameterization that
accommodates both half-life and proportion-based specifications.

#### 2.4.1 Exponential Decay Model

The primary carryover model uses time-since-discontinuation (tsd):

```
carryover_retention = (1/2)^(tsd / t_half)
```

Where:
- `tsd`: weeks elapsed since participant discontinued active treatment
- `t_half`: half-life of carryover effect in weeks

#### 2.4.2 Unified Parameterization

The framework supports two equivalent parameterization modes:

| Mode        | Parameter      | Formula                        | Use Case                    |
|-------------|----------------|--------------------------------|-----------------------------|
| `halflife`  | t½ in weeks    | `(1/2)^(tsd / t_half)`         | Pharmacokinetic modeling    |
| `proportion`| p (0-1)        | Fixed retention proportion     | Simple sensitivity analysis |

**Conversion between modes** (at reference time t):
- Proportion → Half-life: `t_half = -t / log₂(p)`
- Half-life → Proportion: `p = (1/2)^(t / t_half)`

Example: 50% retention at 1 week ↔ t½ = 1 week

#### 2.4.3 Implementation Details

For each participant, the simulation:

1. Tracks when treatment is discontinued (`discontinuation_week`)
2. Computes `tsd` at each subsequent timepoint
3. Applies exponential decay to accumulated BR at discontinuation:
   ```
   BR_off_drug = BR_at_discontinuation × (1/2)^(tsd / t_half)
   ```
4. Generates continuous `carryover_effect` covariate for analysis model

#### 2.4.4 Scale Considerations

**Critical Note**: Hendrickson et al. (2020) used very short half-lives
(t½ = 0, 0.1, 0.2 weeks, i.e., 0.7-1.4 days), reflecting rapid pharmacokinetic
washout. The current implementation explores longer half-lives (0, 0.5, 1.0, 2.0
weeks) to capture potential neurobiological adaptation effects that persist
beyond plasma drug clearance.

| Scale       | t½ Range      | Interpretation                              |
|-------------|---------------|---------------------------------------------|
| Hendrickson | 0-0.2 weeks   | Plasma drug clearance                       |
| pmsimstats  | 0-2.0 weeks   | Clinical effect persistence                 |
| Prazosin PK | ~0.02 weeks   | Plasma half-life (2-3 hours)                |
| PTSD effect | Unknown       | To be estimated from trial data (Section 12)|

Parameters:
- `scale_factor`: 1.0 (biological response only, per Hendrickson)
- `t_half` (half-life): 0, 0.5, 1.0, or 2.0 weeks

### 2.5 Two-Scenario Comparison

The core methodological contribution compares analysis approaches. Hendrickson
et al. (2020) modeled carryover in data generation but did not include it in the
analysis model—they observed that even short carryover caused "precipitous
decline in power" for N-of-1 designs. This framework extends their work by
explicitly testing whether including carryover in the analysis model recovers
lost power.

**Scenario 1 (WITH Carryover Modeling)**:
```r
lmer(response ~ treatment * biomarker + week + carryover_effect +
     (1 | participant_id))
```

**Scenario 2 (WITHOUT Carryover Modeling)** — Hendrickson's original approach:
```r
lmer(response ~ treatment * biomarker + week +
     (1 | participant_id))
```

**Note on Design-Specific Models**: Following Hendrickson, open-label-only
designs use a time × biomarker interaction model rather than treatment ×
biomarker, since all participants experience the same treatment condition:
```r
lmer(response ~ biomarker + week + biomarker:week + (1 | participant_id))
```

### 2.6 Simulation Parameters

**Fixed Parameters**:

| Parameter              | Value | Description                         |
|------------------------|-------|-------------------------------------|
| n_participants         | 70    | Sample size per simulation          |
| n_iterations           | 100   | Monte Carlo replications            |
| BR_rate                | 0.5   | Biological response rate (pts/week) |
| ER_rate                | 0.2   | Expectancy rate (pts/week)          |
| TR_rate                | 0.1   | Time-variant rate (pts/week)        |
| baseline_mean          | 10.0  | Mean baseline score                 |
| within_subject_sd      | 1.8   | Within-subject variability          |
| between_subject_sd     | 2.0   | Between-subject variability         |

**Correlation Structure** (Hendrickson values):

| Parameter      | Value | Description                          |
|----------------|-------|--------------------------------------|
| c.br           | 0.75  | Biological response autocorrelation  |
| c.er           | 0.75  | Expectancy response autocorrelation  |
| c.tr           | 0.75  | Time-variant response autocorrelation|
| c.cf1t         | 0.12  | Same-time cross-correlation          |
| c.cfct         | 0.05  | Different-time cross-correlation     |
| c.bm_baseline  | 0.25  | Biomarker-baseline correlation       |

**Parameter Grid** (Two-Scenario Design):
```r
expand_grid(
  design = c("hybrid", "ol_bdc", "open_label", "crossover", "parallel"),
  biomarker_moderation = c(0, 0.25, 0.35, 0.45, 0.55, 0.65),
  biomarker_correlation = 0.3,
  carryover_halflife = c(0, 0.5, 1.0, 2.0),
  model_carryover = c(TRUE, FALSE)
)
```

The `model_carryover` parameter enables the two-scenario comparison:
- `TRUE`: Include `carryover_effect` term in analysis model
- `FALSE`: Omit carryover term (replicates Hendrickson approach)

When `carryover_halflife = 0`, the `model_carryover` parameter is redundant
(no carryover to model), so these duplicate rows are filtered from the grid.

### 2.7 Statistical Analysis

Power is computed as the proportion of simulations yielding a significant
treatment × biomarker interaction (p < 0.05) from linear mixed-effects models
with random intercepts for participants and CAR(1) correlation structure for
repeated measures.

---

## 3. Implementation

### 3.1 Software Architecture

The framework is implemented as an R package (`pmsimstats2025`) with:

- **Core Functions** (`pm_functions.R`, 1,301 lines): Data generation,
  covariance matrix construction, carryover calculation, and response modeling
- **Simulation Scripts**: Separate scripts for clustered designs
  (`simulation_clustered.R`) and evenly-spaced designs
  (`simulation_evenly_spaced.R`)
- **Visualization**: Publication-ready power heatmaps in Hendrickson style

### 3.2 Key Functions

| Function                              | Purpose                              |
|---------------------------------------|--------------------------------------|
| `build_sigma_matrix()`                | Construct positive-definite covariance matrix |
| `generate_data()`                     | Generate multivariate normal trial data |
| `calculate_carryover()`               | Model carryover decay (exponential, linear, Weibull) |
| `calculate_bio_response_with_interaction()` | Compute biomarker-moderated response |
| `validate_parameter_grid()`           | Pre-simulation parameter validation  |

### 3.3 Reproducibility Infrastructure

- **Docker**: Containerized environment via zzcollab framework
- **renv**: Locked dependencies (252 KB lock file)
- **Makefile**: 229 lines of automation targets
- **Testing**: Unit tests (testthat) and integration tests

### 3.4 Computational Requirements

| Configuration        | Runtime      |
|----------------------|--------------|
| 20 iterations        | 5-15 minutes |
| 100 iterations       | ~25 minutes  |
| 1,000 iterations     | ~4 hours     |

(M1 Mac, 16GB RAM baseline)

---

## 4. Results

### 4.1 Primary Outputs

The framework generates:

1. **Power Heatmaps**: Two-dimensional surfaces showing statistical power across
   biomarker correlation and carryover half-life parameters
2. **Coefficient Summaries**: Mixed model parameter estimates across conditions
3. **Publication Figures**: Hendrickson-style two-panel comparisons (WITH vs.
   WITHOUT carryover modeling)

### 4.2 Key Findings

The two-scenario comparison demonstrates:

- **With Carryover Modeling**: Statistical power remains stable across carryover
  conditions, as the confounding effect is explicitly controlled
- **Without Carryover Modeling**: Progressive power degradation occurs as
  carryover half-life increases, due to unmodeled confounding

### 4.3 Design Comparison

Preliminary results indicate that hybrid designs with dense measurement clusters
at phase transitions provide enhanced sensitivity for detecting biomarker-
treatment interactions compared to traditional evenly-spaced designs.

---

## 5. Project Structure

```
pmsimstats2025/
├── analysis/
│   ├── scripts/           # Simulation and visualization (47 files)
│   ├── output/            # Results and figures
│   ├── data/              # Raw and derived data
│   └── report/            # Manuscript template
├── docs/                  # Documentation (20 files)
├── vignettes/             # Tutorials (4 vignettes)
├── R/                     # Package functions
├── tests/                 # Test suite
├── renv/                  # Dependency management
├── Dockerfile             # Container specification
├── Makefile               # Build automation
└── README.md              # Project overview
```

---

## 6. Key Findings from Hendrickson et al. (2020)

The foundational paper established several critical findings that inform this
framework:

### 6.1 Design Comparison Results

| Design    | Power (no carryover) | Power (t½=0.1 wk) | Notes                    |
|-----------|---------------------|-------------------|--------------------------|
| Open-Label| Low                 | Low               | No within-subject contrast|
| OL+BDC    | Moderate            | Sharply reduced   | Vulnerable to carryover  |
| Crossover | Highest             | Moderately reduced| Most robust to carryover |
| N-of-1    | High                | Sharply reduced   | Best balance of power + practicality |

### 6.2 Critical Carryover Finding

> "The presence of even a short (0.1 weeks) carryover component resulted in a
> precipitous decline in power in both the N-of-1 and, to a lesser but still
> very significant degree, the open label + blinded discontinuation design."
> — Hendrickson et al. (2020), p. 7

This finding motivates the central question of pmsimstats2025: can explicit
carryover modeling in the analysis phase recover the lost power?

### 6.3 Practical Design Recommendation

Hendrickson recommended the N-of-1 design (8 weeks open-label → 4 weeks blinded
discontinuation → 8 weeks crossover) despite slightly lower power than
traditional crossover because it allows all participants to start on active
treatment—critical for recruiting symptomatic patients who would otherwise
decline randomization to placebo.

---

## 7. Extensions Beyond Hendrickson

| Feature                        | Hendrickson (2020)       | This Framework           |
|--------------------------------|--------------------------|--------------------------|
| Three-factor model (BR+ER+TR)  | Core                     | Implemented              |
| BR-only carryover in data gen  | Core                     | Implemented              |
| Expectancy scaling (1.0/0.5)   | Core                     | Implemented              |
| Tabula rasa assumption (TR=ER) | Core                     | Implemented              |
| Carryover half-life range      | 0, 0.1, 0.2 weeks        | 0, 0.5, 1.0, 2.0 weeks   |
| **Carryover in analysis model**| Not included             | Two-scenario comparison  |
| **Biomarker moderation grid**  | Fixed at observed values | Systematic variation     |
| **Parameter validation**       | Not described            | Pre-simulation PD checks |
| **Design-specific models**     | Implicit                 | Explicit documentation   |

**Key Methodological Extension**: Hendrickson observed that carryover caused
"precipitous decline in power" but did not test whether including carryover in
the analysis model could recover this power. This framework directly addresses
that gap through the two-scenario comparison.

---

## 8. Usage

### Quick Start

```bash
# Enter Docker environment
make r

# Run clustered design simulations
Rscript analysis/scripts/simulation_clustered.R

# Generate visualizations
Rscript analysis/scripts/visualize_hendrickson_style.R
```

### Vignettes

1. `01_generate_simulation_data.Rmd` - Data generation workflow
2. `02_analyze_and_visualize.Rmd` - Analysis and visualization
3. `03_apply_to_clinical_data.Rmd` - Clinical application
4. `quickstart.Rmd` - Rapid orientation

---

## 9. Dependencies

**Core Packages**:

- Statistical modeling: lme4, lmerTest, nlme
- Data manipulation: dplyr, tidyverse, data.table
- Visualization: ggplot2, patchwork, viridis
- Matrix operations: MASS, corpcor
- Parallel processing: foreach, doParallel, future, furrr

---

## 10. Citation

```
@software{pmsimstats2025,
  title = {pmsimstats2025: Monte Carlo Simulation for N-of-1 Trial Design},
  year = {2025},
  url = {https://github.com/[repository]},
  note = {R package version 1.0.0}
}
```

**Key Reference**:

Hendrickson RC, Thomas RG, Schork NJ, Raskind MA (2020). Optimizing Aggregated
N-Of-1 Trial Designs for Predictive Biomarker Validation: Statistical Methods
and Theoretical Findings. Frontiers in Digital Health, 2:13.
doi: 10.3389/fdgth.2020.00013

**Related Clinical Trial**: NCT03539614 (VA-funded N-of-1 trial of prazosin for
PTSD, currently enrolling)

---

## 11. License

GPL-3

---

## 12. Research Plan: Prazosin-PTSD Trial Analysis

This section outlines the methodological plan for applying the pmsimstats2025
framework to analyze actual prazosin trial data for PTSD treatment, with three
primary goals: (1) evaluating treatment efficacy, (2) characterizing carryover
effects, and (3) enhancing simulation modularity.

### 12.1 Goal 1: Did Prazosin Improve PTSD Symptoms?

#### 12.1.1 Analysis Challenge

The hybrid N-of-1 design presents unique analytical challenges compared to
traditional parallel-group RCTs:

- Participants serve as their own controls across treatment phases
- The 4-path randomization structure creates heterogeneous treatment sequences
- Carryover effects may confound treatment phase comparisons
- Biomarker moderation requires interaction testing, not simple group means

#### 12.1.2 Proposed Analysis Strategy

**Primary Analysis Model**:
```r
lmer(ptsd_score ~ treatment * biomarker_centered + week + carryover_effect +
     (1 | participant_id),
     data = prazosin_data)
```

This model estimates:

- `treatment`: Main effect of prazosin vs. placebo/off-drug
- `biomarker_centered`: Biomarker association with outcome
- `treatment:biomarker_centered`: Differential treatment effect by biomarker
  (primary hypothesis)
- `week`: Linear time trend (natural symptom trajectory)
- `carryover_effect`: Residual drug effect during washout periods
- `(1 | participant_id)`: Between-subject heterogeneity

**Interpretation Framework**:

| Coefficient               | Interpretation if Significant (p < 0.05)       |
|---------------------------|------------------------------------------------|
| treatment                 | Average prazosin effect across biomarker range |
| treatment:biomarker       | Biomarker moderates treatment response         |
| carryover_effect          | Residual effects persist after discontinuation |

#### 12.1.3 Secondary Analyses

1. **Path-Specific Effects**: Stratify by randomization path to assess treatment
   sequence effects
   ```r
   lmer(ptsd_score ~ treatment * biomarker * path + week + carryover +
        (1 | participant_id))
   ```

2. **Sensitivity to Carryover Specification**: Compare results with and without
   carryover term to quantify confounding risk (mirrors simulation two-scenario
   comparison)

3. **Individual Treatment Effects**: Extract participant-level BLUPs (Best
   Linear Unbiased Predictors) for personalized effect estimates
   ```r
   ranef(model)$participant_id
   ```

4. **Responder Analysis**: Classify participants by biomarker-predicted response
   magnitude

#### 12.1.4 Required Data Elements

| Variable            | Description                                    | Source          |
|---------------------|------------------------------------------------|-----------------|
| participant_id      | Unique identifier                              | Trial database  |
| path                | Randomization path (1-4)                       | Randomization   |
| week                | Assessment timepoint                           | Protocol        |
| treatment           | On prazosin (1) vs. off (0)                    | Protocol/path   |
| ptsd_score          | PTSD symptom measure (e.g., PCL-5, CAPS-5)     | Assessment      |
| biomarker           | Pre-specified moderator variable               | Baseline        |
| tsd                 | Time since discontinuation (weeks)             | Derived         |
| carryover_effect    | Computed decay function                        | Derived         |

#### 12.1.5 Implementation Steps

1. Import and validate clinical dataset structure
2. Compute treatment status indicators from path and week
3. Calculate time-since-discontinuation (tsd) for each observation
4. Apply carryover decay function (requires half-life specification)
5. Center biomarker variable at sample mean
6. Fit primary mixed-effects model
7. Extract coefficients, confidence intervals, and p-values
8. Conduct sensitivity analyses
9. Generate effect plots (treatment effect as function of biomarker)

---

### 12.2 Goal 2: Characterize Carryover in Actual Data

#### 12.2.1 The Carryover Quantification Problem

Carryover effects in prazosin trials reflect residual drug activity after
discontinuation. Characterizing carryover requires addressing:

1. **Temporal dynamics**: How quickly does the effect decay?
2. **Magnitude**: What proportion of on-drug effect persists?
3. **Individual variation**: Does carryover differ across participants?
4. **Outcome specificity**: Do different PTSD domains show different carryover?

#### 12.2.2 Proposed Carryover Quantification Structure

**A. Pharmacokinetic Anchoring**

Prazosin pharmacokinetics provide initial parameter bounds:

| Parameter                | Value           | Implication                        |
|--------------------------|-----------------|------------------------------------|
| Plasma half-life         | 2-3 hours       | Drug cleared within 24 hours       |
| CNS effect duration      | 8-12 hours      | Direct pharmacological effect      |
| Symptom response lag     | Days to weeks   | Clinical effect may persist longer |

The disconnect between plasma kinetics (hours) and symptom response (weeks)
suggests the relevant carryover operates at the neurobiological adaptation
level, not plasma drug concentration.

**B. Empirical Carryover Estimation**

Two complementary approaches:

**Approach 1: Model-Based Estimation**

Treat half-life as a parameter to estimate from data:

```r
# Grid search over candidate half-lives
half_lives <- c(0.5, 1, 2, 3, 4, 6, 8) # weeks

results <- map_dfr(half_lives, function(t_half) {
  data <- prazosin_data |>
    mutate(carryover = ifelse(tsd > 0,
                              (1/2)^(tsd / t_half),
                              0))

  model <- lmer(ptsd_score ~ treatment * biomarker + week + carryover +
                (1 | participant_id), data = data)

  tibble(
    t_half = t_half,
    AIC = AIC(model),
    BIC = BIC(model),
    carryover_coef = fixef(model)["carryover"],
    carryover_se = sqrt(vcov(model)["carryover", "carryover"])
  )
})
```

Select half-life minimizing AIC/BIC or maximizing carryover coefficient
precision.

**Approach 2: Non-Parametric Visualization**

Examine raw response patterns during washout periods:

```r
prazosin_data |>
  filter(tsd > 0) |>  # Post-discontinuation observations only
  ggplot(aes(x = tsd, y = ptsd_score, group = participant_id)) +
  geom_line(alpha = 0.3) +
  geom_smooth(aes(group = NULL), method = "loess") +
  labs(x = "Weeks Since Discontinuation",
       y = "PTSD Score",
       title = "Empirical Carryover Decay Pattern")
```

**C. Multi-Model Comparison**

Compare decay functional forms:

| Model       | Formula                              | Parameters     |
|-------------|--------------------------------------|----------------|
| Exponential | `(1/2)^(tsd / t_half)`               | t_half         |
| Linear      | `max(0, 1 - tsd / t_total)`          | t_total        |
| Weibull     | `exp(-(tsd / lambda)^k)`             | lambda, k      |
| Step        | `ifelse(tsd < threshold, 1, 0)`      | threshold      |

Use likelihood ratio tests or AIC comparison to select best-fitting model.

#### 12.2.3 Identification of Optimal Response Measure Sequence

**A. The Sequence Selection Problem**

PTSD assessments may include multiple instruments or subscales:

- PCL-5 total score
- PCL-5 cluster subscores (intrusion, avoidance, cognition/mood, arousal)
- CAPS-5 clinician rating
- Secondary measures (sleep quality, nightmare frequency)

Different measures may show different:

- Sensitivity to treatment
- Carryover persistence
- Biomarker moderation patterns

**B. Proposed Selection Strategy**

**Step 1: Univariate Screening**

For each candidate outcome measure:
```r
screen_outcome <- function(outcome_var, data) {
  formula <- as.formula(paste(outcome_var,
    "~ treatment * biomarker + week + (1 | participant_id)"))
  model <- lmer(formula, data = data)

  tibble(
    outcome = outcome_var,
    treatment_effect = fixef(model)["treatment"],
    interaction_effect = fixef(model)["treatment:biomarker"],
    interaction_p = summary(model)$coefficients["treatment:biomarker", "Pr(>|t|)"],
    icc = performance::icc(model)$ICC_adjusted,
    r2_marginal = performance::r2(model)$R2_marginal
  )
}
```

**Step 2: Carryover Sensitivity Comparison**

For each outcome, compare carryover model fit:
```r
compare_carryover_sensitivity <- function(outcome_var, data) {
  # Model without carryover
  m0 <- lmer(paste(outcome_var, "~ treatment * biomarker + week +
             (1 | participant_id)"), data = data)

  # Model with carryover
  m1 <- lmer(paste(outcome_var, "~ treatment * biomarker + week + carryover +
             (1 | participant_id)"), data = data)

  tibble(
    outcome = outcome_var,
    delta_AIC = AIC(m0) - AIC(m1),  # Positive = carryover helps
    carryover_coef = fixef(m1)["carryover"],
    lr_test_p = anova(m0, m1)$`Pr(>Chisq)`[2]
  )
}
```

**Step 3: Composite Ranking**

Rank outcomes by weighted criteria:

| Criterion                        | Weight | Rationale                          |
|----------------------------------|--------|------------------------------------|
| Interaction effect size          | 0.30   | Primary hypothesis sensitivity     |
| Carryover model improvement      | 0.25   | Carryover characterization utility |
| ICC (reliability)                | 0.20   | Measurement precision              |
| Variance explained (R²)          | 0.15   | Signal-to-noise ratio              |
| Clinical interpretability        | 0.10   | Practical relevance                |

**C. Recommended Sequence**

Based on typical PTSD trial structure, prioritize:

1. **Primary**: PCL-5 total score (validated, comparable across studies)
2. **Secondary**: CAPS-5 clinician rating (gold standard, if available)
3. **Exploratory**: PCL-5 arousal cluster (most likely to show prazosin-specific
   effects given mechanism of action)
4. **Mechanistic**: Nightmare frequency (direct target of prazosin)

#### 12.2.4 Deliverables for Goal 2

1. Estimated carryover half-life with confidence interval
2. Comparison of decay model functional forms
3. Ranked list of outcome measures by carryover sensitivity
4. Visualization of empirical carryover decay patterns
5. Participant-level carryover heterogeneity assessment

---

### 12.3 Goal 3: Enhance Simulation Modularity for Carryover Analysis

#### 12.3.1 Current Implementation Status

The existing `calculate_carryover()` function in `pm_functions.R` supports three
decay models (exponential, linear, Weibull) but has limited modularity:

**Current Limitations**:

- Carryover parameters embedded in simulation scripts, not configuration files
- Single half-life value applied uniformly (no individual variation)
- No mechanism to import empirically-estimated carryover from real data
- Limited integration with sensitivity analysis workflows

#### 12.3.2 Proposed Modular Architecture

**A. Carryover Configuration Object**

Create a dedicated configuration structure:

```r
carryover_config <- list(
  # Model specification
  model = "exponential",           # "exponential", "linear", "weibull", "step"

  # Core parameters
  half_life = 1.0,                 # weeks (can be vector for sensitivity)
  half_life_sd = 0.3,              # between-subject SD (optional)

  # Weibull-specific

  shape = 1.0,                     # shape parameter k

  # Component application
  apply_to = c("BR"),              # which response components
  scale_factors = c(BR = 1.0, ER = 0.0, TR = 0.0),


  # Estimation mode
  source = "specified",            # "specified", "estimated", "imported"
  estimated_from = NULL,           # path to estimation results

  # Sensitivity grid
  sensitivity_grid = c(0, 0.5, 1.0, 2.0)  # half-lives to test
)
```

**B. Modular Function Interface**

Refactor `calculate_carryover()` to accept configuration object:

```r
calculate_carryover <- function(tsd, config = carryover_config) {

  # Validate inputs

  stopifnot(config$model %in% c("exponential", "linear", "weibull", "step"))
  stopifnot(tsd >= 0)

  # Individual-level half-life (if heterogeneity specified)
  t_half <- if (!is.null(config$half_life_sd)) {
    rnorm(1, config$half_life, config$half_life_sd) |> max(0.01)
  } else {
    config$half_life
  }

  # Calculate decay
  decay <- switch(config$model,
    exponential = (1/2)^(tsd / t_half),
    linear = pmax(0, 1 - tsd / t_half),
    weibull = exp(-(tsd / t_half)^config$shape),
    step = as.numeric(tsd < t_half)
  )

  decay
}
```

**C. Integration with Empirical Estimates**

Create import pathway for real-data carryover estimates:

```r
import_carryover_estimate <- function(estimation_results_path) {
  # Load estimation results from Goal 2 analysis
  est <- readRDS(estimation_results_path)

  carryover_config <- list(
    model = est$best_model,
    half_life = est$half_life_estimate,
    half_life_sd = est$half_life_se,
    source = "estimated",
    estimated_from = estimation_results_path
  )

  carryover_config
}
```

**D. Sensitivity Analysis Wrapper**

Systematic sensitivity analysis across carryover specifications:

```r
run_carryover_sensitivity <- function(base_config, half_life_grid,
                                       n_iterations = 100) {

  results <- map_dfr(half_life_grid, function(t_half) {
    config <- modifyList(base_config, list(half_life = t_half))

    # Run simulation with this carryover specification
    sim_result <- run_simulation(
      carryover_config = config,
      n_iterations = n_iterations
    )

    sim_result |>
      mutate(carryover_half_life = t_half)
  })

  results
}
```

#### 12.3.3 Implementation Roadmap

**Phase 1: Configuration Externalization**

1. Define `carryover_config` schema in `R/carryover_config.R`
2. Create YAML configuration file support (`config/carryover.yaml`)
3. Refactor existing scripts to read from configuration
4. Add validation functions for configuration integrity

**Phase 2: Function Modularity**

1. Refactor `calculate_carryover()` to use configuration object
2. Add individual-level half-life heterogeneity option
3. Implement additional decay models (power law, two-compartment)
4. Create unit tests for each decay model

**Phase 3: Empirical Integration**

1. Build carryover estimation pipeline (Goal 2 deliverables)
2. Create `import_carryover_estimate()` function
3. Add validation for imported parameters
4. Document workflow for data-driven carryover specification

**Phase 4: Sensitivity Framework**

1. Implement `run_carryover_sensitivity()` wrapper
2. Create visualization functions for sensitivity results
3. Add parallel execution support for large grids
4. Generate sensitivity analysis vignette

#### 12.3.4 Deliverables for Goal 3

1. Externalized carryover configuration system (YAML + R interface)
2. Refactored `calculate_carryover()` with modular design
3. Import pathway for empirically-estimated carryover parameters
4. Sensitivity analysis framework with visualization
5. Updated vignettes demonstrating new modularity features
6. Unit test coverage for all carryover functions

---

### 12.4 Integration Timeline and Dependencies

```
Goal 1 (Prazosin Efficacy)
    │
    ├── Requires: Clinical data import
    ├── Requires: Carryover specification (from Goal 2 or assumed)
    └── Output: Treatment effect estimates
                    │
                    ▼
Goal 2 (Carryover Characterization)
    │
    ├── Requires: Clinical data with washout observations
    ├── Requires: Multiple outcome measures for comparison
    └── Output: Estimated carryover parameters
                    │
                    ▼
Goal 3 (Simulation Modularity)
    │
    ├── Input: Carryover estimates from Goal 2
    ├── Enables: Data-driven simulation parameterization
    └── Output: Enhanced simulation framework
                    │
                    ▼
    Validation: Re-simulate with empirical carryover
                Compare simulated vs. observed power
```

**Dependency Notes**:

- Goal 1 can proceed with assumed carryover values, then be re-run after Goal 2
- Goal 2 requires post-discontinuation data points to estimate decay
- Goal 3 benefits from Goal 2 but can proceed in parallel with placeholder
  values

---

### 12.5 Prazosin-Specific Considerations

#### 12.5.1 Mechanism of Action

Prazosin is an α1-adrenergic receptor antagonist. In PTSD:

- Reduces noradrenergic hyperactivity during sleep
- Primary effect on trauma-related nightmares
- Secondary effects on hyperarousal symptoms
- Mechanism suggests rapid offset (receptor binding ~hours) but potentially
  prolonged clinical effect (neuroadaptation ~weeks)

#### 12.5.2 Expected Carryover Profile

Based on pharmacology and clinical experience:

| Component           | Expected Carryover | Rationale                           |
|---------------------|-------------------|-------------------------------------|
| Nightmare reduction | Short (days)      | Direct receptor blockade            |
| Sleep quality       | Short-medium      | Secondary to nightmare reduction    |
| Daytime arousal     | Medium (1-2 weeks)| Sustained neuroadaptation           |
| Overall PTSD score  | Medium            | Composite of above                  |

#### 12.5.3 Biomarker Candidates

For treatment moderation analysis, consider:

- Baseline PTSD severity (dose-response relationship)
- Nightmare frequency at baseline (target symptom)
- Autonomic reactivity measures (physiological target)
- Trauma type/chronicity (moderates treatment response)
- Prior medication history (sensitization effects)

---

### 12.6 Summary: Three-Goal Integration

| Goal | Primary Question | Key Method | Output |
|------|------------------|------------|--------|
| 1 | Did prazosin work? | Mixed-effects model with carryover | Effect estimates, p-values |
| 2 | How long does effect persist? | Empirical carryover estimation | Half-life estimate, decay model |
| 3 | Can simulation use real parameters? | Modular configuration system | Enhanced simulation framework |

The three goals form a coherent pipeline: clinical data analysis (Goal 1)
informs carryover characterization (Goal 2), which parameterizes improved
simulations (Goal 3), enabling validation of the analytical approach against
known ground truth.

---

## Appendix A: Covariance Matrix Construction

The framework implements a three-stage partitioned approach for constructing
positive-definite covariance matrices:

1. **Stage 1**: Within-component temporal correlations (AR(1) structure)
2. **Stage 2**: Cross-component correlations (same-time and different-time)
3. **Stage 3**: Biomarker-baseline-response correlations

Eigenvalue validation ensures positive definiteness before simulation proceeds.

---

## Appendix B: Carryover Decay Models

Three decay models are supported:

1. **Exponential** (default): `(1/2)^(scale × t / t_half)`
2. **Linear**: `max(0, 1 - t / t_half)`
3. **Weibull**: `exp(-(t / t_half)^shape)`

The exponential model aligns with pharmacokinetic first-order elimination and is
recommended for most applications.

---

*Document generated: 2026-02-11*
*Project: pmsimstats2025 v1.0.0*
