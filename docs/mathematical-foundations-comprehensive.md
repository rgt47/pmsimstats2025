# Mathematical Foundations of N-of-1 Trial Simulation: A Comprehensive Analysis

## Document Overview

This document provides an in-depth, comprehensive synthesis of all mathematical
elements underlying the `pmsimstats2025` N-of-1 clinical trial simulation
framework. The analysis draws from extensive documentation in the `docs/`
directory and represents a unified treatment of the project's theoretical
foundations.

---

## Part I: Covariance Matrix Architecture

### 1.1 The Core Identity: Σ = DRD

The fundamental mathematical identity governing the simulation is the
decomposition of the covariance matrix into correlation and variance components:

$$\Sigma = D \cdot R \cdot D$$

where:

- **Σ** is the full covariance matrix (26×26 for 8-timepoint designs)
- **R** is the correlation matrix (dimensionless, all diagonal entries = 1)
- **D** is a diagonal matrix of standard deviations

This decomposition separates two distinct concerns:

1. **Correlation structure (R)**: How variables relate to each other
   (independent of scale)
2. **Variance magnitude (D)**: The scale of each variable's variability

### 1.2 Block Partitioning Strategy

The full covariance matrix is partitioned into blocks for computational
efficiency:

$$\Sigma = \begin{bmatrix} \Sigma_{11} & \Sigma_{12} \\ \Sigma_{12}' & \Sigma_{22} \end{bmatrix}$$

where:

| Block | Dimension | Contents |
|-------|-----------|----------|
| Σ₁₁ | (3n_tp) × (3n_tp) | Response components (BR, ER, TR) at all timepoints |
| Σ₂₂ | 2 × 2 | Biomarker and baseline |
| Σ₁₂ | (3n_tp) × 2 | Cross-covariances |

For 8-timepoint designs: Σ₁₁ is 24×24, giving a full matrix of 26×26.

### 1.3 Intra-Block Structure

Within Σ₁₁, there is further structure:

$$\Sigma_{11} = \begin{bmatrix}
\Sigma_{BR} & C_{BR,ER} & C_{BR,TR} \\
C_{ER,BR} & \Sigma_{ER} & C_{ER,TR} \\
C_{TR,BR} & C_{TR,ER} & \Sigma_{TR}
\end{bmatrix}$$

Each diagonal block (e.g., Σ_BR) is an n_tp × n_tp matrix with:

- Diagonal entries = 1 (unit variance for correlation matrix)
- Off-diagonal entries = autocorrelation parameter (e.g., c.br = 0.75)

Cross-blocks (e.g., C_{BR,ER}) contain:

- Diagonal entries = c.cf1t (same-timepoint cross-correlation)
- Off-diagonal entries = c.cfct (different-timepoint cross-correlation)

---

## Part II: Conditional Multivariate Normal Sampling

### 2.1 The Two-Stage Algorithm

Rather than inverting the full 26×26 matrix, the simulation uses the Schur
complement to enable efficient two-stage sampling:

**Stage 1**: Generate (biomarker, baseline) from Σ₂₂

$$(x_2) \sim \mathcal{N}(\mu_2, \Sigma_{22})$$

**Stage 2**: Generate response components conditional on Stage 1

$$x_1 | x_2 \sim \mathcal{N}(\mu_{1|2}, \Sigma_{cond})$$

where:

$$\mu_{1|2} = \mu_1 + \Sigma_{12} \Sigma_{22}^{-1} (x_2 - \mu_2)$$

$$\Sigma_{cond} = \Sigma_{11} - \Sigma_{12} \Sigma_{22}^{-1} \Sigma_{12}'$$

### 2.2 Computational Advantages

This approach reduces computational complexity:

| Operation | Full Matrix | Two-Stage |
|-----------|-------------|-----------|
| Matrix inversion | 26×26 | 2×2 |
| Cholesky decomposition | O(26³) | O(2³) + O(24³) |
| Numerical stability | Lower | Higher |

The key insight is that Σ₂₂ is only 2×2, making its inversion trivial and
numerically stable.

### 2.3 Conditional Covariance (Σ_cond)

The conditional covariance Σ_cond is the central matrix for data generation.
It must be:

1. **Positive definite** (all eigenvalues > 0)
2. **Well-conditioned** (condition number κ < 100)
3. **Numerically stable** for Cholesky decomposition

---

## Part III: Positive Definiteness Constraints

### 3.1 Mathematical Definition

A symmetric matrix Σ is positive definite if and only if:

$$\mathbf{v}^T \Sigma \mathbf{v} > 0 \quad \forall \mathbf{v} \neq \mathbf{0}$$

Equivalently: All eigenvalues λᵢ > 0.

### 3.2 Verification Methods

Three methods are used to verify positive definiteness:

**Method 1: Sylvester's Criterion**

All leading principal minors must be positive:

$$\det(\Sigma_k) > 0 \quad \text{for } k = 1, 2, \ldots, n$$

where Σₖ is the k×k upper-left submatrix.

**Method 2: Eigenvalue Criterion**

$$\lambda_{min}(\Sigma) > 0$$

This is the primary method used in `build_sigma_matrix()`.

**Method 3: Cholesky Decomposition**

If Σ = LLᵀ exists (with L lower triangular, positive diagonal), then Σ is PD.

### 3.3 Gershgorin Circle Theorem

For eigenvalue bounds, the Gershgorin theorem states:

$$\lambda_i \in \bigcup_{j=1}^{n} \left\{ z \in \mathbb{C} : |z - a_{jj}| \leq R_j \right\}$$

where $R_j = \sum_{k \neq j} |a_{jk}|$ is the off-diagonal row sum.

For correlation matrices (diagonal = 1):

$$\lambda_i \in [1 - R_j, 1 + R_j]$$

This implies that if off-diagonal sums exceed 1, negative eigenvalues become
possible.

### 3.4 Correlation Hierarchy Constraint

The critical constraint for PD preservation is the correlation hierarchy:

$$\boxed{c_{cfct} < c_{cf1t} < \min(c_{tv}, c_{pb}, c_{br})}$$

**Rationale**: Measurements at different times cannot be more correlated than
measurements at the same time. Violating this creates impossible constraint
systems.

**Hendrickson's values satisfy this**:

- c.cfct = 0.05 < c.cf1t = 0.12 < c.br = 0.75 ✓

### 3.5 Biomarker Correlation Bound

Empirical analysis establishes:

$$c_{bm} \leq 0.6$$

for designs with 8 timepoints. Higher values risk non-PD matrices.

---

## Part IV: Condition Number Analysis

### 4.1 Definition and Interpretation

The condition number measures numerical stability:

$$\kappa = \frac{\lambda_{max}}{\lambda_{min}}$$

**Interpretation scale**:

| κ Value | Status | Interpretation |
|---------|--------|----------------|
| κ < 10 | Excellent | Near-perfect numerical stability |
| 10 < κ < 100 | Good | Acceptable for most applications |
| 100 < κ < 1000 | Poor | Numerical precision concerns |
| κ > 1000 | Severe | Cholesky decomposition unreliable |

### 4.2 The AR(1) Problem

AR(1) correlation structure: $\rho_{ij} = \rho^{|t_i - t_j|}$ with ρ = 0.75

For evenly-spaced designs (weeks 2, 8, 14, 20), time gaps are 6, 12, 18 weeks:

| Time Gap | AR(1) Correlation |
|----------|-------------------|
| 6 weeks | 0.75⁶ ≈ 0.178 |
| 12 weeks | 0.75¹² ≈ 0.032 |
| 18 weeks | 0.75¹⁸ ≈ 0.0056 |

**Correlation ratio**: 0.178 / 0.0056 ≈ 32:1

This extreme ratio creates condition numbers κ ≈ 600-14,621 for evenly-spaced
designs.

### 4.3 Why Clustering Solves This

Clustered designs (weeks 4, 8, 9, 10, 11, 12, 16, 20) have mixed gaps:

- Short gaps: 1-4 weeks → correlations 0.32-0.75
- Long gaps: 8-16 weeks → correlations 0.10-0.32

This creates **intermediate correlation values** spanning the correlation
space more evenly, producing balanced eigenvalue distributions with κ ≈ 56-61.

### 4.4 Alternative Correlation Structures

Five structures are implemented for flexible design:

| Structure | Formula | Best For |
|-----------|---------|----------|
| AR(1) | ρ^{\|t_i - t_j\|} | Clustered designs |
| Exponential | exp(-λ\|t_i - t_j\|) | Evenly-spaced (κ ≈ 18) |
| Power Law | (1 + α\|t_i - t_j\|)^{-β} | Evenly-spaced (κ ≈ 13) |
| Matérn | (1 + √3d/ρ)exp(-√3d/ρ) | Flexible smoothness |
| Rational Quadratic | (1 + t²/(2αλ²))^{-α} | Best stability (κ ≈ 10) |

---

## Part V: Three-Factor Response Model

### 5.1 Response Decomposition

Each participant's response is decomposed into three components:

$$\text{Response}_t = \text{Baseline} + \text{BR}_t + \text{ER}_t + \text{TR}_t$$

where:

| Component | Name | Meaning | Driver |
|-----------|------|---------|--------|
| BR | Biological Response | True pharmacological effect | Treatment status, biomarker |
| ER | Expectancy Response | Placebo/expectancy effect | Blinding status |
| TR | Time-variant Response | Natural disease progression | Time, individual variation |

### 5.2 Component Mean Formulas

**Biological Response (BR)**:

$$\text{BR}_{mean,t} = \begin{cases}
w_t \cdot r_{eff} & \text{if on treatment} \\
w_{accum} \cdot \gamma \cdot f(t_{off}) & \text{if off treatment (carryover)}
\end{cases}$$

where:

- $w_t$ = weeks on drug up to timepoint t
- $r_{eff} = r_0 (1 + \beta \cdot \text{BiomarkerCentered})$ = effective BR rate
- $\gamma$ = carryover proportion
- $f(t_{off})$ = decay function (exponential, power law, etc.)

**Expectancy Response (ER)**:

$$\text{ER}_{mean,t} = e_t \cdot r_{exp}$$

where:

- $e_t \in \{0, 0.5, 1\}$ = expectancy factor (design-dependent)
- $r_{exp}$ = expectancy response rate

**Time-variant Response (TR)**:

$$\text{TR}_{mean,t} = f_{time}(t)$$

representing natural disease trajectory independent of treatment.

### 5.3 Design-Specific Expectancy Patterns

| Design | Weeks with e=1.0 | Rationale |
|--------|------------------|-----------|
| OL+BDC | Weeks ≤ 16 | Open-label phase (know on treatment) |
| Hybrid | Weeks 4, 8 only | Only unblinded baseline period |
| Crossover | Varies by period | Blinding maintained throughout |

---

## Part VI: Biomarker×Treatment Interaction

### 6.1 Two-Level Mechanism

The biomarker-treatment interaction operates at two levels:

**Level 1: Covariance Structure (c.bm)**

The correlation between biomarker and BR components:

$$\text{Corr}(\text{biomarker}, \text{BR}_t) = c_{bm}$$

This creates covariance-level association: participants with high biomarker
values tend to have correlated BR values.

**Level 2: Mean Modulation (biomarker_moderation)**

The biomarker modifies treatment effect magnitude:

$$r_{eff} = r_0 (1 + \beta \cdot \text{BiomarkerCentered})$$

where β = `biomarker_moderation` parameter.

### 6.2 Mathematical Relationship

These two mechanisms are mathematically distinct:

| Aspect | c.bm | biomarker_moderation |
|--------|------|----------------------|
| Level | Covariance | Mean |
| Effect | Correlation in draws | Scaling of treatment effect |
| Parameter range | 0-0.6 (PD constraint) | 0-0.65 (power analysis) |
| Impact | Noise structure | Signal magnitude |

### 6.3 Detectability

The biomarker×treatment interaction is tested via:

$$H_0: \beta_{biomarker \times treatment} = 0$$

in the mixed model:

```r
lmer(response ~ biomarker * treatment + week + (1|participant_id))
```

Power depends on both c.bm and biomarker_moderation values.

---

## Part VII: Carryover Effect Modeling

### 7.1 The Fundamental Question

Should carryover effects modify correlation parameters in addition to means?

**Answer**: NO. Carryover affects means only.

### 7.2 Mathematical Argument

**Approach 1 (Incorrect)**: Carryover modifies correlations

$$\rho_{BR}(t_{on}, t_{off}) = \rho_0 \cdot g(\text{carryover})$$

This would mean participants' BR correlations change based on their treatment
history—biologically implausible for most drugs.

**Approach 2 (Correct)**: Carryover affects means only

$$\mu_{BR}(t_{off}) = \mu_{BR}(t_{on}) \cdot h(t_{since\_discontinuation})$$

The residual drug effect enters through the expected value, not the variance
structure.

### 7.3 Double-Counting Problem

If carryover modified both means and correlations:

1. Mean effect: BR_mean higher during off-period due to residual effect
2. Correlation effect: BR draws more correlated during off-period

This would count the same biological phenomenon twice, inflating carryover
impact and biasing estimates.

### 7.4 Five Carryover Models

| Model | Formula | Parameters | Best For |
|-------|---------|------------|----------|
| Fixed Proportion | γ (constant) | 1 | Simple drugs |
| Exponential | exp(-λt) | 2 | Most psychiatric drugs |
| Power Law | (1+αt)^{-β} | 3 | Tissue accumulation |
| Two-Compartment | Ae^{-αt} + Be^{-βt} | 4 | Known biphasic kinetics |
| Delayed Elimination | Lag + exponential | 3 | Prodrugs |

### 7.5 Half-Life Relationship

For exponential decay:

$$t_{1/2} = \frac{\ln(2)}{\lambda} \approx \frac{0.693}{\lambda}$$

Example values:

| λ (week⁻¹) | t₁/₂ (weeks) | Application |
|------------|--------------|-------------|
| 0.05 | 13.9 | Slow-elimination drugs |
| 0.10 | 6.9 | Typical psychiatric medications |
| 0.15 | 4.6 | Faster clearance |

---

## Part VIII: Mixed-Effects Model Specification

### 8.1 Core Model Formula

$$\text{Response}_{it} = \beta_0 + \beta_1 \cdot \text{Biomarker}_i + \beta_2 \cdot \text{Treatment}_{it} + \beta_3 \cdot \text{Biomarker}_i \times \text{Treatment}_{it} + \beta_4 \cdot \text{Week}_{it} + u_i + \epsilon_{it}$$

where:

- $\beta_3$ = biomarker×treatment interaction (primary endpoint)
- $u_i \sim N(0, \sigma_u^2)$ = random intercept
- $\epsilon_{it} \sim N(0, \sigma^2)$ = residual error

### 8.2 Model Variants

**With Carryover Term**:

```r
lmer(response ~ biomarker * treatment + week + carryover_effect +
     (1|participant_id))
```

**Without Carryover (Hendrickson Original)**:

```r
lmer(response ~ biomarker * treatment + week + (1|participant_id))
```

### 8.3 Random Effects Options

| Structure | Formula | Interpretation |
|-----------|---------|----------------|
| Random intercept | (1\|participant_id) | Individual baseline differences |
| Random slope for drug | (1+treatment\|participant_id) | Individual treatment response variation |
| Random slope for time | (1+week\|participant_id) | Individual disease trajectory variation |

### 8.4 Singularity Handling

When random effects structure is too complex for the data:

1. Model returns singular fit warning
2. Fallback to simpler structure (random intercept only)
3. Log singularity occurrence for diagnostic tracking

---

## Part IX: Design Comparison Framework

### 9.1 Five Trial Designs

| Design | Measurement Points | Schedule Type | Carryover Applicability |
|--------|-------------------|---------------|------------------------|
| Hybrid | 8 | Clustered | Moderate-High |
| OL+BDC | 8 | Clustered | High |
| Open-Label | 4 | Evenly-spaced | None |
| Crossover | 4 | Evenly-spaced | Critical |
| Parallel | 4 | Evenly-spaced | None |

### 9.2 Measurement Schedules

**Clustered (Hybrid)**:

Weeks: {4, 8, 9, 10, 11, 12, 16, 20}

- Dense cluster at transition (weeks 9-12)
- Captures carryover decay dynamics
- AR(1) works well (κ ≈ 56-61)

**Clustered (OL+BDC)**:

Weeks: {4, 8, 12, 16, 17, 18, 19, 20}

- Dense cluster at discontinuation (weeks 16-20)
- Extended washout observation
- AR(1) works well (κ ≈ 60)

**Evenly-Spaced**:

Weeks: {2, 8, 14, 20}

- Uniform 6-week gaps
- Requires alternative correlation structures
- Exponential or Power Law recommended (κ ≈ 13-18)

### 9.3 Expectancy Vulnerability

Based on analysis of expectancy allocation patterns:

**OL+BDC**: MORE vulnerable to expectancy confounding

- High expectancy (e=1.0) for 80% of trial (weeks 0-16)
- ER and BR effects coupled at randomization point
- Sudden drop in both components at week 16

**Hybrid**: LESS vulnerable

- High expectancy only for 20% of trial (weeks 4-8)
- Most measurements during blinded phase (e=0.5)
- Cleaner separation of BR and ER effects

---

## Part X: Validation and Quality Control

### 10.1 Pre-Simulation Sigma Validation

The `validate_parameter_grid()` function tests all parameter combinations
before expensive Monte Carlo runs:

1. Test each unique sigma structure for PD
2. Compute condition numbers
3. Flag ill-conditioned matrices (κ > 100)
4. Stop with diagnostic error if validation fails

### 10.2 Build-Time Validation

Within `build_sigma_matrix()`:

1. Compute eigenvalues
2. Check all eigenvalues > 0
3. Return NULL for non-PD combinations
4. Log diagnostic output

### 10.3 Correlation Clamping

To prevent invalid correlations from mean-ratio scaling:

```r
scaled_correlation <- pmax(-0.99, pmin(0.99, scaled_correlation))
```

This ensures correlation values remain in valid range (-1, 1) even when
mean ratios are extreme due to treatment and carryover effects.

### 10.4 Comparison: Reject vs. Auto-Fix

| Approach | Our Implementation | Hendrickson |
|----------|-------------------|-------------|
| Philosophy | Reject invalid | Auto-fix and continue |
| PD handling | Return NULL, skip | `make.positive.definite()` |
| Diagnostics | Explicit warnings | Silent |
| Result | Only valid matrices | Complete grid (may be perturbed) |

Our stricter approach ensures mathematical correctness at the cost of
potentially incomplete parameter grids.

---

## Part XI: Statistical Power Analysis

### 11.1 Power Definition

$$\text{Power} = P(\text{Reject } H_0 | H_1 \text{ true})$$

Estimated via Monte Carlo simulation:

$$\hat{\text{Power}} = \frac{\#(\text{iterations with } p < 0.05)}{N_{iterations}}$$

### 11.2 Power Determinants (Ranked by Impact)

Based on comprehensive parameter analysis:

| Rank | Parameter | Impact Mechanism |
|------|-----------|------------------|
| 1 | Design | Measurement schedule, carryover modeling |
| 2 | biomarker_moderation | Signal magnitude |
| 3 | carryover_t1half | Signal clarity during off-periods |
| 4 | c.bm | Biomarker-response correlation |
| 5 | n_participants | Statistical precision |
| 6 | Correlation structure | Numerical stability |

### 11.3 Carryover Cost Analysis

Comparing models with/without carryover term:

| Carryover Level | Without Term (Power) | With Term (Power) | Type I Error |
|-----------------|---------------------|-------------------|--------------|
| γ = 0 | ~0.50 | ~0.50 | 0.05 |
| γ = 0.5 | ~0.35 | ~0.45 | 0.07-0.10 (without) |
| γ = 1.0 | ~0.25 | ~0.35 | 0.10-0.12 (without) |

Ignoring carryover leads to:

- 10-20% power reduction
- Inflated Type I error
- Biased treatment effect estimates

---

## Part XII: Key Mathematical Insights

### 12.1 The Eigenvalue-Conditioning Nexus

The condition number problem is fundamentally about eigenvalue balance:

$$\kappa = \frac{\lambda_{max}}{\lambda_{min}}$$

**Good conditioning** requires eigenvalues to be similar in magnitude, which
requires correlation values to span the correlation space evenly rather than
clustering at extreme values.

### 12.2 The Clustering Insight

Clustered measurement schedules succeed because they create **intermediate
correlations**:

- AR(1) with 1-week gaps: ρ = 0.75
- AR(1) with 4-week gaps: ρ = 0.32
- AR(1) with 8-week gaps: ρ = 0.10

This natural variation from clustering balances eigenvalues.

Evenly-spaced designs with uniform large gaps (6 weeks) create extreme
correlation ratios (32:1), leading to extreme eigenvalue ratios.

### 12.3 The Carryover-Correlation Separation

The principle that carryover affects means only (not correlations) is
essential for:

1. **Avoiding double-counting** of carryover effects
2. **Maintaining interpretability** of correlation parameters
3. **Preserving positive definiteness** (carryover-dependent correlations
   could violate PD constraints for some treatment histories)

### 12.4 The Hierarchy Constraint

The correlation hierarchy ($c_{cfct} < c_{cf1t} < c_{autocorr}$) encodes
the temporal coherence principle:

> Measurements at different times cannot be more correlated than measurements
> at the same time.

Violating this creates mathematically impossible constraint systems.

---

## Part XIII: Practical Recommendations

### 13.1 Parameter Selection

**For Hendrickson comparison**:

```r
c.tv = 0.8, c.pb = 0.8, c.br = 0.75
c.cf1t = 0.12, c.cfct = 0.05
c.bm = 0.3 (vary: 0, 0.3, 0.6)
```

**For general use**:

```r
c.tv = 0.65, c.pb = 0.65, c.br = 0.65
c.cf1t = 0.18, c.cfct = 0.09
c.bm = 0.3
```

### 13.2 Correlation Structure Selection

| Schedule Type | Recommended Structure | Expected κ |
|---------------|----------------------|------------|
| Clustered (8pt) | AR(1) | 56-61 |
| Evenly-spaced (4pt) | Power Law or Rational Quadratic | 10-15 |
| Evenly-spaced (8pt) | Rational Quadratic | 15-25 |

### 13.3 Carryover Model Selection

| Drug Characteristics | Recommended Model |
|---------------------|-------------------|
| Unknown kinetics | Fixed Proportion (γ = 0.5) |
| Psychiatric medications | Exponential (λ = 0.10) |
| Tissue accumulation | Power Law |
| Known biphasic | Two-Compartment |

### 13.4 Validation Workflow

1. Define parameter grid
2. Run `validate_parameter_grid()` before simulation
3. Check condition numbers (all κ < 100)
4. Investigate any rejected combinations
5. Proceed only with validated grid

---

## Appendix A: Key Formulas Reference

### Covariance Matrix

$$\Sigma = D \cdot R \cdot D$$

### Conditional Mean

$$\mu_{1|2} = \mu_1 + \Sigma_{12} \Sigma_{22}^{-1} (x_2 - \mu_2)$$

### Conditional Covariance

$$\Sigma_{cond} = \Sigma_{11} - \Sigma_{12} \Sigma_{22}^{-1} \Sigma_{12}'$$

### Condition Number

$$\kappa = \frac{\lambda_{max}}{\lambda_{min}}$$

### AR(1) Correlation

$$\rho_{ij} = \rho^{|t_i - t_j|}$$

### Exponential Correlation

$$\rho_{ij} = \exp(-\lambda |t_i - t_j|)$$

### Carryover Exponential Decay

$$\text{Carryover}_t = \text{Effect}_0 \cdot \exp\left(-\frac{\ln(2)}{t_{1/2}} \cdot t_{off}\right)$$

### Half-Life Relationship

$$t_{1/2} = \frac{\ln(2)}{\lambda}$$

### Hierarchy Constraint

$$c_{cfct} < c_{cf1t} < \min(c_{tv}, c_{pb}, c_{br})$$

---

## Appendix B: Document Sources

This synthesis draws from the following documentation files:

1. `sigma_matrix_derivation.tex` - Core Σ = DRD derivation
2. `positive_definiteness_constraints.tex` - PD mathematical analysis
3. `correlation_structure_design.tex` - Block structure design
4. `carryover_correlation_theory.tex` - Carryover-correlation separation
5. `correlation_structure_white_paper.Rmd` - AR(1) vs alternatives analysis
6. `carryover_modeling_white_paper.Rmd` - Five carryover models
7. `mixed_model_comparison_and_best_practices.md` - Model specification guidance
8. `technical_differences_scaling_and_pd.tex` - Implementation differences
9. `correlation_parameters_guide.tex` - Parameter interpretation guide
10. `whitepaper.tex` - Crossover vs N-of-1 comparison

---

*Document generated: December 17, 2025*
*Repository: pmsimstats2025*
*Total mathematical documentation reviewed: 10+ documents, ~15,000 lines*
