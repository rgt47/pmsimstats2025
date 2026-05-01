# Modeling Carryover Effects in N-of-1 Trials: A Methodological Review

## Abstract

This white paper reviews current best practices for modeling carryover effects in
N-of-1 trials and crossover designs. We synthesize recommendations from the
statistical literature, compare analytical approaches, and provide guidance for
model specification when carryover effects are anticipated. Particular attention
is given to the parameterization of carryover decay functions and the separation
of treatment effects from residual carryover effects in mixed effects models.

## 1. Introduction

N-of-1 trials are multi-period crossover trials performed in a single individual,
designed to determine optimal treatment for that patient. A critical challenge in
these designs is the potential for carryover effects, where the influence of a
treatment persists after it has been discontinued, potentially biasing estimates
of subsequent treatment effects.

Carryover effects threaten the validity of treatment comparisons in crossover
studies. While washout periods are the most common design-based solution, they
may not always be feasible due to ethical concerns, patient burden, or practical
constraints. When washout periods are insufficient or absent, statistical methods
must account for carryover effects in the analysis model.

This review addresses three key questions:

1. When should carryover effects be modeled versus prevented by design?
2. What model specifications are recommended for carryover-adjusted analyses?
3. How should carryover decay be parameterized?

## 2. Design-Based Prevention versus Statistical Adjustment

### 2.1 The Case for Washout Periods

The statistical literature consistently emphasizes that preventing carryover
effects through study design is preferable to statistical adjustment. In a
methodological review of 74 N-of-1 trials, only 43.2% included a washout period,
with higher rates in studies of rare conditions (61.5%) compared to non-rare
conditions (39.3%) (Chen et al., 2024).

For pharmaceutical interventions, the standard recommendation is to set washout
periods at 3-4 times (or more) of the plasma elimination half-life (Lee & Lim,
2021). This duration ensures that treatment concentrations return to negligible
levels before the subsequent period begins.

### 2.2 Limitations of Design-Based Approaches

Washout periods are not always feasible:

- **Ethical concerns**: Taking patients off effective treatment may cause harm
- **Patient burden**: Extended trial duration reduces compliance
- **Practical constraints**: Some conditions require continuous management
- **Unknown pharmacodynamics**: The effect half-life may differ from the
  pharmacokinetic half-life

When washout periods cannot eliminate carryover effects, analytical approaches
become necessary.

## 3. Statistical Methods for Carryover-Adjusted Analysis

### 3.1 Comparison of Analytical Approaches

Chen et al. (2014) compared four methods for analyzing N-of-1 trials:

| Method | Carryover Handling | Performance |
|--------|-------------------|-------------|
| Paired t-test | Not modeled | Optimal Type I error, highest power when carryover absent |
| Mixed effects (difference) | Not modeled | Problematic Type I error |
| Mixed effects (full) | Explicitly modeled | Unaffected by carryover, lower power |
| Meta-analysis | Not modeled | Lower power, Type I error issues |

**Key finding**: The paired t-test is recommended for normally distributed data
when carryover effects are absent. When carryover effects are present, the mixed
effects model with explicit carryover terms provides unbiased estimates but with
reduced statistical power.

### 3.2 Mixed Effects Model Specification

The recommended mixed effects model for carryover-adjusted analysis includes:

```
Y_ijk = μ + α_i + π_j + τ_k + λ_{k(j-1)} + ε_ijk
```

Where:
- `Y_ijk` = response for subject i in period j receiving treatment k
- `μ` = overall mean
- `α_i` = random subject effect
- `π_j` = fixed period effect
- `τ_k` = fixed treatment effect (parameter of interest)
- `λ_{k(j-1)}` = carryover effect from treatment in period j-1
- `ε_ijk` = residual error

This formulation separates the direct treatment effect (`τ`) from the residual
carryover effect (`λ`), enabling unbiased estimation of both.

### 3.3 Evidence for Including Carryover Terms

Granholm et al. (2021) demonstrated through simulation that when carryover
effects are present:

- Models without carryover terms produce **biased treatment estimates**
- Models without carryover terms have **coverage below nominal 95%**
- Models with carryover terms provide **unbiased estimates** and
  **appropriate coverage probability**

The authors recommend using a mixed model with a carryover term and an
unstructured correlation matrix to obtain unbiased treatment effect estimates.

## 4. Parameterization of Carryover Decay

### 4.1 Decay Function Options

Three primary approaches exist for parameterizing how carryover effects decay
over time:

#### 4.1.1 Linear Decay

The simplest approach assumes carryover decreases linearly from 1 (at
discontinuation) to 0 over a predefined period:

```
Z(t) = max(0, 1 - t/T)
```

Where `t` is time since discontinuation and `T` is the assumed carryover
duration. This approach is computationally simple but may not reflect
pharmacological reality.

#### 4.1.2 Exponential (Geometric) Decay

More aligned with pharmacological principles, exponential decay follows:

```
Z(t) = (1/2)^(t/t_half)
```

Where `t_half` is the half-life of the carryover effect. This parameterization
reflects the common assumption that treatment effects decay at a constant
proportional rate.

#### 4.1.3 Distributed Lag Model

The most flexible approach, introduced by Chen & Ho (2023), uses lag coefficients
without imposing a specific decay function:

```
Y_t = β_0 X_t + β_1 X_{t-1} + β_2 X_{t-2} + ... + β_L X_{t-L} + ε_t
```

Where `β_0` captures the immediate effect and `β_1, ..., β_L` capture carryover
at successive lags. Bayesian shrinkage priors prevent overfitting by
constraining lag coefficients to decay smoothly.

### 4.2 Choosing a Decay Function

The choice of decay function should be guided by:

1. **Pharmacological knowledge**: If the drug's half-life is known, exponential
   decay with `t_half` matching the pharmacological half-life is appropriate.

2. **Uncertainty about decay shape**: If the decay pattern is unknown, the
   distributed lag model with Bayesian shrinkage allows the data to inform the
   decay shape while maintaining smoothness.

3. **Computational constraints**: Linear or exponential decay is simpler to
   implement and estimate than distributed lag models.

Chen & Ho (2023) found that when the true carryover follows exponential decay,
the Koyck distributed lag model (which assumes geometric decay) performs
optimally. However, when decay patterns deviate from exponential (e.g., slow
absorption curves), more flexible approaches are needed.

## 5. The Treatment-Carryover Separation Problem

### 5.1 The Aliasing Issue

A fundamental challenge in carryover modeling is that carryover effects can be
aliased (confounded) with other design factors. In the classic 2×2 crossover
design, carryover is confounded with sequence, making statistical separation
impossible without additional assumptions.

### 5.2 Implications for Model Parameterization

When modeling both treatment effects and carryover effects, care must be taken
to ensure these are identifiable. Two common approaches:

**Approach A: Combined Exposure Variable**

Combine treatment and carryover into a single "effective exposure" variable:

```
effective_exposure = treatment + carryover_decay × prior_treatment
```

This approach assumes the carryover effect is proportional to the treatment
effect. Power may be reduced when carryover is present because variance in the
combined predictor is reduced.

**Approach B: Separate Effect Terms**

Model treatment and carryover as separate additive terms:

```
Y = β_1 × treatment + β_2 × carryover_indicator + ...
```

This approach allows carryover to have a different magnitude than the treatment
effect but requires sufficient variation in both terms for identification.

### 5.3 Recommendation

For biomarker-moderated treatment effects, we recommend **Approach B** with
separate terms:

```
Y ~ treatment * biomarker + carryover * biomarker + covariates
```

This allows the biomarker moderation of active treatment to be estimated
independently from the moderation of decayed carryover effects.

## 6. Practical Recommendations

### 6.1 When to Model Carryover

Model carryover effects when:

1. Washout periods are absent or insufficient
2. The treatment is known to have effects persisting beyond administration
3. The pharmacodynamic half-life differs substantially from the
   pharmacokinetic half-life
4. Prior analyses suggest differential carryover effects between treatments

### 6.2 Model Specification Checklist

1. **Include explicit carryover terms** rather than ignoring carryover or
   relying on preliminary tests

2. **Separate treatment from carryover** when possible to avoid aliasing

3. **Use appropriate decay parameterization**:
   - Exponential decay if half-life is known
   - Linear decay for simplicity when decay shape is less critical
   - Distributed lag models when flexibility is needed

4. **Account for correlation structure** using appropriate random effects and
   autocorrelation specifications

5. **Consider power implications**: Including carryover terms typically reduces
   power; ensure sample size is adequate

### 6.3 Reporting Standards

When reporting carryover-adjusted analyses:

1. State the assumed carryover decay function and its parameters
2. Report both treatment and carryover effect estimates with confidence
   intervals
3. Conduct sensitivity analyses under alternative decay assumptions
4. Compare results with and without carryover adjustment

## 7. Application to Biomarker-Moderated Treatment Effects

When the research question involves biomarker moderation of treatment effects
in designs with potential carryover, the model should include:

```
response ~ treatment_term × biomarker +
           carryover_term × biomarker +
           time_covariates +
           (1 | subject)
```

Where:
- `treatment_term` indicates active treatment periods
- `carryover_term` indicates decayed residual effects during off-treatment
  periods
- The interaction with biomarker allows for differential moderation

This formulation maintains power to detect biomarker moderation of the active
treatment effect while properly accounting for carryover.

## 8. Conclusions

The statistical literature provides clear guidance on carryover effect modeling:

1. **Prevention is preferable to adjustment**: Design studies with adequate
   washout periods when feasible.

2. **Include carryover terms when effects are suspected**: Models without
   carryover terms produce biased estimates when carryover is present.

3. **Separate treatment from carryover**: Avoid combining these into a single
   exposure variable when possible, as this can reduce power and complicate
   interpretation.

4. **Choose decay functions based on pharmacological knowledge**: Exponential
   decay aligned with known half-lives is often appropriate; distributed lag
   models offer flexibility when decay patterns are uncertain.

5. **Accept power trade-offs**: Properly modeling carryover typically reduces
   statistical power, but this is preferable to obtaining biased estimates.

## References

Chen, X., & Ho, M. H. (2023). Analysis of N-of-1 trials using Bayesian
distributed lag model with autocorrelated errors. *Statistics in Medicine*,
42(11), 1718-1736. https://doi.org/10.1002/sim.9697

Chen, Y., Zhang, Z., Xu, P., Zhu, T., & Zhan, S. (2014). A comparison of four
methods for the analysis of N-of-1 trials. *PLOS ONE*, 9(2), e87752.
https://doi.org/10.1371/journal.pone.0087752

Chen, Y. J., et al. (2024). A methodological review of randomised n-of-1 trials.
*Trials*, 25, 100. https://doi.org/10.1186/s13063-024-08100-1

Granholm, A., Munch, M. W., Andersen-Ranberg, N. C., & Perner, A. (2021).
Statistical methods for testing carryover effects: A mixed effects model
approach. *Contemporary Clinical Trials Communications*, 22, 100764.
https://doi.org/10.1016/j.conctc.2021.100764

Lee, S., & Lim, H. S. (2021). Considerations for crossover design in clinical
study. *Korean Journal of Anesthesiology*, 74(4), 293-299.
https://doi.org/10.4097/kja.21165

Schmid, C. H., & Yang, X. (2014). Statistical design and analytic considerations
for N-of-1 trials. In *Design and implementation of N-of-1 trials: A user's
guide* (Chapter 4). Agency for Healthcare Research and Quality.
https://effectivehealthcare.ahrq.gov/products/n-1-trials/research-2014-1

Wienke, S., Weber, S., Guzek, V., & Wegscheider, K. (2023). Comparison of
Bayesian Networks, G-estimation and linear models to estimate causal treatment
effects in aggregated N-of-1 trials with carry-over effects. *BMC Medical
Research Methodology*, 23, 171.
https://doi.org/10.1186/s12874-023-02012-5

---

*White paper prepared for the pmsimstats2025 project, February 2026.*
