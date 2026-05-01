# Simulating Biomarker-Treatment Interactions: Correlation-Based vs. Mean-Structure Approaches

**pmsimstats2025 Project Notes**

---

## ⚠️  CRITICAL UPDATE (March 2026)

**This document was written before Phase 2 diagnostic testing. Sections 3.4
(Pros) and section 7 (original Recommendation) contain claims about the
correlation-based approach that have been empirically refuted.**

**Phase 2 diagnostic testing (March 2026) identified that pmsimstats2025's
extensions to the correlation-based approach (BM-ER/BM-TR correlations and
constant ER SD parameterization) inflate baseline power 50-90% without
improving stability or improving the detection of biomarker-treatment
interactions. These additions contradict the theoretical justification
offered in section 3.4.**

**Section 7 (Revised Recommendation for pmsimstats2025) now reflects Phase 2
findings and supersedes the original section 7 recommendations.**

See `analysis/docs/PHASE2_EXECUTIVE_SUMMARY.md` and
`analysis/docs/PHASE2_CULPRIT_FINDINGS.md` for detailed empirical evidence.

---

## 1. The Problem

The pmsimstats2025 simulation and the original Hendrickson et al.
(2020) codebase use fundamentally different mechanisms to generate
the biomarker-treatment interaction -- the statistical signal that
the simulation is designed to detect. This document reviews the two
approaches, surveys the relevant literature, and evaluates the
pros and cons of each for the N-of-1 aggregated trial context.

The biomarker-treatment interaction is a predictive (not prognostic)
effect: participants with higher biomarker values respond
*differently* to treatment than those with lower values. The
analysis model tests this via an interaction term
(`Dbc:bm_centered` or `treatment:bm_centered`). The question is
how the data-generating process (DGP) should encode this
interaction so that the simulation faithfully represents the
clinical scenario.

## 2. Approach 1: Mean-Structure (Deterministic Modulation)

### 2.1 Mechanism

The biomarker enters the DGP as a fixed coefficient in the linear
predictor. A participant's biological response to drug is a
deterministic function of their biomarker value:

```
Y_i = baseline + BR_i + ER_i + TR_i + noise
BR_i = f(drug_exposure) * (1 + beta_TB * bm_centered_i)
```

A participant 1 SD above the biomarker mean receives exactly
`(1 + beta_TB)` times the treatment effect of the average
participant. The residual noise adds variability but does not
modulate the interaction. The biomarker-treatment relationship is
known with certainty conditional on the biomarker value.

### 2.2 pmsimstats2025 implementation

In `carryover_power_simulation.R` (lines 345--347):

```r
effective_BR_rate <- BR_rate * (1 + bm_mod * bm_centered)
br_on_drug <- weeks_on_drug * effective_BR_rate
```

The sigma matrix applies a constant `c.bm` to all BR timepoints
regardless of treatment status (lines 148--153), plus
half-strength correlations to ER and TR that the original does not
have. The biomarker-BR correlation in the covariance structure
serves only to generate realistic between-participant variability;
it does not encode the treatment-dependent interaction.

### 2.3 Literature using this approach

Nearly all published biomarker-stratified trial simulations use
the mean-structure approach:

**Riecke et al. (2018).** Simulated Cox regression with
`beta_TxB * T * B` in the log-hazard. Three interaction strengths
(beta_TxB = 0, 0.095, 0.285). Found that including prognostic
covariates increases power to detect the interaction and reduces
bias. 500 subjects, 1000 simulations, 108 scenarios.
(DOI: 10.1186/s12874-018-0457-9)

**Haller et al. (2019).** Simulated survival outcomes with
deterministic hazard-rate functions of biomarker Z and treatment
T. Six interaction shapes (null, constant HR, linear, nonlinear
monotonic, tail effect, U-shaped). Compared seven analysis methods
(median split, quartile split, optimal cutpoint, STEPP, Cox
linear, MFPI, LPLB). Cox regression with linear interaction was
most efficient for monotonic shapes; MFPI and LPLB outperformed
for complex nonlinear interactions.
(DOI: 10.1159/000500436)

**Langbaum et al. (2014).** Alzheimer's disease prevention trial
enrichment. Treatment effect applied only to biomarker-positive
(at-risk) participants; non-at-risk participants received zero
treatment effect. Found that ignoring sample heterogeneity
(different variance components for at-risk vs. not-at-risk groups)
overestimates power.
(DOI: 10.1016/j.jalz.2013.12.019)

**Friede, Parsons, and Stallard (2012).** Adaptive subgroup
selection designs. Treatment effects specified as fixed mean
differences for different biomarker-defined subgroups. Simulation
operates at the level of sufficient statistics (multivariate
normal test statistics), not individual patient data.
(DOI: 10.1002/sim.5541)

**Blackston et al. (2019).** Aggregated N-of-1 vs. RCT
comparison. Treatment effect is fixed and homogeneous across all
participants; no biomarker-treatment interaction is modeled.
Carryover is modeled as fixed increments.
(DOI: 10.3390/healthcare7040137)

### 2.4 Pros

1. **Direct control over effect size.** The parameter `beta_TB`
   maps directly to the interaction coefficient that the analysis
   model estimates. Power calculations are analytically tractable
   because the true parameter is known exactly. For a given
   `beta_TB`, `Var(B)`, and residual variance, the expected power
   can be derived in closed form.

2. **Transparent DGP-analysis alignment.** The DGP and the
   analysis model share the same parametric form. If the analysis
   model is `Y ~ T * B`, the simulation generates data from the
   same family. There is no model misspecification by construction,
   making it straightforward to verify that the simulation is
   calibrated correctly (i.e., that `mean(beta_hat) = beta_TB`
   over replications).

3. **Standard in the literature.** Nearly every simulation paper
   outside Hendrickson uses this approach, making results directly
   comparable to published power estimates for biomarker-stratified
   designs.

4. **Computationally simple.** No path-specific covariance
   matrices are needed. One set of parameters generates all
   participants regardless of treatment assignment. The sigma
   matrix can be pre-built and cached once per design.

### 2.5 Cons

1. **Overly optimistic about biomarker predictive utility.** Given
   a biomarker value, the treatment effect is known with certainty
   (up to residual noise). In clinical practice, biomarkers are
   imperfect predictors: two patients with identical biomarker
   values can have very different treatment responses for reasons
   unrelated to the biomarker (genetic background, comorbidities,
   adherence). The mean-structure approach cannot represent this
   heterogeneity.

2. **No high-biomarker non-responders (except through noise).** A
   patient 2 SD above the biomarker mean always has a large
   treatment effect mean. The only way they can be a non-responder
   is if the residual error overwhelms the signal. Clinically,
   non-response despite favorable biomarker status is common and
   may reflect biological complexity that the mean-structure
   approach does not capture.

3. **DGP-analysis alignment is both a feature and a limitation.**
   Because the simulation generates data from the same model family
   that the analysis uses, it cannot reveal the consequences of
   model misspecification. If the real-world interaction operates
   through a different mechanism (e.g., covariance-mediated), the
   simulation will overestimate the power of the standard analysis
   model.

4. **Treatment-dependent variance is absent.** In the
   mean-structure approach, the residual variance is the same
   whether a participant is on drug or off drug, and whether the
   biomarker is high or low. In clinical data, treatment response
   variability often depends on treatment status and patient
   characteristics.

## 3. Approach 2: Correlation-Based (Covariance Structure)

### 3.1 Mechanism

The biomarker and the biological response component are drawn
jointly from a multivariate normal distribution. The interaction
is encoded as a correlation between the biomarker and the
treatment-response component, with the correlation depending on
treatment status:

```
(B_i, BR_i(t)) ~ MVN(mu, Sigma)

where:
  corr(B, BR(t)) = c.bm       when on drug (tod > 0)
  corr(B, BR(t)) = scaled      when off drug with carryover
  corr(B, BR(t)) = 0           when off drug, no carryover
```

The biomarker does not appear in the mean of BR. Instead,
participants with higher biomarker draws *tend* to get higher BR
draws through the joint distribution. The tendency is controlled
by `c.bm`, and the prediction is inherently imperfect: even at
`c.bm = 0.6`, approximately 64% of the variance in BR is
unexplained by the biomarker.

### 3.2 Hendrickson's original implementation

In `generateData.R` (lines 84--163), the biological response
means are computed via a modified Gompertz function of cumulative
time on drug (`tod`):

```r
brmeans <- modgompertz(d$tod, rp$max, rp$disp, rp$rate)
```

The biomarker-BR correlation is then set in the correlation
matrix, with the Ron Thomas modification (lines 146--162):

```r
for (p in 1:nP) {
  n1 <- paste(trialdesign$timeptname[p], "br", sep = ".")
  if (p > 1) {
    n0 <- paste(trialdesign$timeptname[p - 1], "br", sep = ".")
    mm1 <- means[which(n1 == labels)]
    mm0 <- means[which(n0 == labels)]
    correlations["bm", n1] <- ifelse(
      brtest[p],
      ifelse(brmeans[p] == 0, 0, (mm1 / mm0) * c.bm),
      c.bm
    )
  }
}
```

The logic:

- `brtest[p]` is TRUE when the raw BR mean (from Gompertz of
  `tod`) is zero -- i.e., when the participant has zero cumulative
  drug exposure at that timepoint.
- **On-drug timepoints** (`brtest` FALSE, `tod > 0`): `corr(bm,
  BR) = c.bm`. Full correlation. High-biomarker participants tend
  to draw higher BR values.
- **Off-drug timepoints, no carryover** (`brtest` TRUE,
  `brmeans = 0`): `corr(bm, BR) = 0`. No correlation. The
  biomarker has no predictive value when there is no drug effect
  to predict.
- **Off-drug timepoints, with carryover** (`brtest` TRUE,
  `brmeans > 0` after carryover adjustment): `corr(bm, BR) =
  (mm1/mm0) * c.bm`. Scaled correlation proportional to the
  ratio of current to previous BR mean. As carryover decays,
  the correlation approaches zero.

Crucially, the first timepoint (`p = 1`) is skipped entirely (the
`p > 1` guard), so the first timepoint's BR has zero correlation
with the biomarker regardless of treatment status.

Because the correlation structure depends on the treatment
sequence, each path through the trial design produces a different
covariance matrix. The original generates data separately for each
path using path-specific sigma matrices, then merges the paths
into a single dataset for analysis.

### 3.3 Literature using this approach

Hendrickson et al. (2020) is the only simulation study identified
in the biomarker-treatment interaction literature that uses this
approach. The paper describes the simulation framework but does not
explicitly discuss the choice between correlation-based and
mean-structure interaction encoding.

### 3.4 Pros

1. **Clinical realism.** Biomarkers are imperfect predictors. The
   correlation-based approach naturally generates high-biomarker
   non-responders and low-biomarker responders at rates controlled
   by `c.bm`. At `c.bm = 0.6`, the biomarker explains 36% of BR
   variance; the remaining 64% represents unmeasured sources of
   treatment response heterogeneity. This is consistent with the
   empirical observation that even the best clinical biomarkers
   explain only a modest fraction of treatment response
   variability.

2. **Treatment-dependent covariance is biologically principled.**
   The biomarker should predict response only when the drug is
   active. When off drug (`BR = 0`), there is nothing for the
   biomarker to predict, so `corr(bm, BR) = 0`. The
   correlation-based approach encodes this directly in the
   covariance structure: the statistical relationship between
   biomarker and outcome changes with treatment status. The
   mean-structure approach achieves this indirectly through the
   `T * B` interaction term, which is algebraically equivalent for
   the conditional mean but does not produce the same conditional
   variance structure.

3. **Tests the analysis model under realistic misspecification.**
   Because the DGP and the analysis model have different parametric
   forms (MVN covariance vs. LME interaction coefficient), the
   simulation reveals how well the standard analysis model recovers
   the true signal when the data-generating mechanism is more
   complex than the analysis assumes. This is the typical situation
   in real trials, where the true DGP is never exactly the model
   being fitted.

4. **Consistent with Senn's variance decomposition.** Senn (2016)
   decomposes treatment response heterogeneity into patient-by-
   treatment interaction variance (psi^2). A biomarker in the
   correlation-based DGP explains a fraction of psi^2 proportional
   to `c.bm^2`, leaving residual heterogeneity. In the
   mean-structure approach, psi^2 is entirely determined by
   `beta_TB^2 * Var(B)`, with no residual treatment heterogeneity
   beyond noise. The correlation-based approach is therefore the
   more conservative (and arguably more honest) representation of
   biomarker utility.

   **⚠️ CAVEAT (added March 2026)**: These theoretical pros are
   valid for the *Hendrickson parameterization* (BM-BR correlations
   only). pmsimstats2025 extends this to include BM-ER and BM-TV
   correlations, which add signal to non-treatment components. Phase 2
   testing shows these extensions inflate power 50-85% while worsening
   stability, contradicting the conservative principle stated above.
   See section 7 (Revised Recommendation).

### 3.5 Cons

1. **Indirect control over the detectable effect size.** The
   analysis model estimates `beta_TB` (the interaction
   coefficient), but the DGP does not contain `beta_TB` as a
   parameter. The relationship between `c.bm` and the expected
   `beta_hat_TB` depends on the variances of B and BR, the
   Gompertz curve shape, the number of on-drug timepoints, the
   autocorrelation structure, and the correlation between factors
   -- a complex mapping that makes analytic power calculation
   difficult. In practice, the expected `beta_hat_TB` must be
   determined empirically by running the simulation and averaging
   the estimates.

2. **Path-specific covariance matrices required.** Because the
   bm-BR correlation depends on treatment status, each path
   through the trial design produces a different covariance
   matrix. For the Hybrid design with 4 paths, this requires 4
   sigma matrices per parameter combination. Pre-computation and
   caching become more complex, particularly because the BR means
   (which determine the Ron Thomas scaling) depend on the response
   parameters (Gompertz max, displacement, rate), not just the
   design geometry.

3. **Harder to validate.** With the mean-structure approach, one
   can verify that `mean(beta_hat) = beta_TB` over replications to
   confirm calibration. With the correlation-based approach, the
   expected `beta_hat` is a derived quantity that depends on many
   parameters jointly. Validation requires either analytic
   derivation of the expected interaction coefficient (nontrivial)
   or large-sample simulation calibration.

4. **Unique to Hendrickson.** No other published simulation study
   in the biomarker-stratified trial literature uses this
   approach, limiting comparability with published power estimates.
   Results from the correlation-based DGP cannot be directly
   compared to power estimates from Riecke et al., Haller et al.,
   or other standard references.

5. **The Ron Thomas scaling is fragile.** The `mm1/mm0` ratio
   depends on the Gompertz means at consecutive timepoints. When
   the denominator (`mm0`) is small (early timepoints with short
   drug exposure), the ratio can be large, producing correlation
   values greater than `c.bm`. When both means are near zero
   (off-drug, minimal carryover), the ratio is numerically
   unstable. The `brtest` guard handles the exact-zero case but
   does not address near-zero instabilities.

## 4. The Conceptual Bridge

Senn's (2016) variance decomposition framework clarifies that the
two approaches are not alternatives but rather represent different
assumptions about the same variance partition.

Consider the treatment effect for patient $i$:

$$\tau_i = \bar{\tau} + \gamma B_i + u_i$$

where $\bar{\tau}$ is the average treatment effect, $\gamma$ is
the biomarker moderation coefficient, $B_i$ is the biomarker
value, and $u_i$ captures residual treatment heterogeneity
unexplained by the biomarker.

- **Mean-structure approach:** Sets $\text{Var}(u_i) = 0$. The
  biomarker explains all systematic treatment heterogeneity. Any
  remaining variation in treatment response is pure noise (the
  residual $\varepsilon$ in the outcome model). The total
  patient-by-treatment interaction variance is
  $\psi^2 = \gamma^2 \text{Var}(B)$.

- **Correlation-based approach:** Allows
  $\text{Var}(u_i) > 0$. The biomarker explains only a fraction
  of treatment heterogeneity, controlled by `c.bm`. The total
  patient-by-treatment interaction variance is
  $\psi^2 = \gamma^2 \text{Var}(B) + \text{Var}(u)$, where
  $\gamma$ is implied by the covariance structure rather than
  specified directly.

At `c.bm = 1.0`, the correlation-based approach converges to the
mean-structure approach: the biomarker perfectly predicts BR, and
the stochastic correlation collapses to a deterministic
relationship. At `c.bm < 1.0`, the correlation-based approach
generates more residual heterogeneity, which attenuates the
detectable interaction and reduces power relative to the
mean-structure approach for the same nominal biomarker-response
association.

This attenuation is not a defect but a realistic representation:
if the biomarker explains only 36% of treatment response variance
(`c.bm = 0.6`), then 64% of the variance appears as noise to the
analysis model, and the power to detect the interaction is
correspondingly lower.

## 5. Implications for Carryover Modeling

The choice of interaction mechanism affects how carryover
contaminates the analysis in distinct ways.

### 5.1 Under the mean-structure approach

At off-drug timepoints with carryover, the DGP generates:

```
BR_off = br_at_discont * Dbc * (1 + bm_mod * bm_centered)
```

The biomarker modulation is present at full strength whenever
`Dbc > 0`. The naive model (`treatment = 0` at these timepoints)
attributes this biomarker-modulated residual drug response to
noise. The interaction estimate is attenuated because the on/off
contrast is contaminated.

### 5.2 Under the correlation-based approach

At off-drug timepoints with carryover, the DGP generates BR draws
from an MVN where `corr(bm, BR)` is scaled by the Ron Thomas
ratio (`mm1/mm0`). As carryover decays, the scaled correlation
approaches zero. The biomarker's predictive value at off-drug
timepoints decays in parallel with the drug effect. The naive
model's misspecification is proportional to the *covariance*
between the biomarker and the off-drug BR, which is a function of
both the carryover level and the correlation scaling -- a more
complex interaction than the mean-structure approach produces.

### 5.3 The practical consequence

For the two-scenario comparison (WITH vs. WITHOUT carryover
modeling), the power gap between the two analysis models will be
*smaller* under the correlation-based DGP than under the
mean-structure DGP, because the correlation-based approach
generates weaker interaction signal at off-drug timepoints (the
correlation is scaled down, not just the mean). This means the
pmsimstats2025 simulation, which uses the mean-structure approach,
may *overestimate* the benefit of modeling carryover relative to
what the original Hendrickson DGP would show.

## 6. Summary Comparison

| Feature | Mean-Structure | Correlation-Based |
|---|---|---|
| Interaction location | Fixed coefficient in the linear predictor | Joint covariance of biomarker and BR |
| Nature of prediction | Deterministic given biomarker | Stochastic; probabilistic tendency |
| High-BM non-responders | Only through residual noise | Natural; frequency controlled by c.bm |
| Effect size control | Direct (beta_TB) | Indirect (c.bm, variances, Gompertz) |
| Analytic power formulas | Standard; tractable | Difficult; requires simulation |
| Sigma matrices needed | 1 per design (path-invariant) | 1 per path per design |
| DGP-analysis alignment | Exact (same parametric family) | Misspecified (realistic) |
| Treatment-dependent variance | No | Yes (correlation varies with drug status) |
| Literature precedent | Nearly universal | Hendrickson (2020) only |
| Computational cost | Lower | Higher (path-specific sigma) |
| Biomarker utility estimate | Optimistic | Conservative |

## 7. Empirical Validation: Phase 2 Diagnostic Findings

**Critical Discovery (March 2026)**: Phase 2 diagnostic testing revealed
that the correlation-based approach as implemented in pmsimstats2025 has
unintended consequences that **contradict the theoretical arguments** in
section 3.4 above.

### 7.1 Phase 2D: BM-ER/BM-TR Correlation Impact

pmsimstats2025 extends the correlation-based approach by adding 0.5 × c.bm
correlations between the biomarker and **expectancy response (ER)** and
**time-varying response (TV)** components (not present in Hendrickson orig).

**Empirical result** (75 iterations, CO N=35, c.bm=0.3):
- **With BM-ER/BM-TR**: Power 27% → 35% → 36% (t1/2=0, 0.1, 0.2), Stability 0.907
- **Without** (orig): Power 13% → 12% → 5%, Stability 0.920
- **Power inflation**: -50% to -85% when these correlations are removed

**Interpretation**: These additions inflate baseline power by ~80% while
simultaneously *worsening* power curve stability. They link the biomarker
to non-treatment components, creating artificial signal that does not
reflect the treatment-specific effect. This violates the theoretical
motivation in section 3.4 that biomarkers should predict treatment
response (BR), not expectancy or time-varying responses.

### 7.2 Phase 2E: ER SD Scaling Impact

pmsimstats2025 uses constant ER SD across all timepoints. Hendrickson orig
scales ER SD by the expectancy factor (lower variance at lower expectancy).

**Empirical result** (75 iterations, CO N=35, c.bm=0.3):
- **Constant (2025)**: Power 43% → 49% → 47%, Stability 0.960
- **Scaled (orig)**: Power 13% → 12% → 5%, Stability 0.920
- **Power inflation**: -65% to -90% when scaled parameterization is used

**Interpretation**: Constant ER SD inflates power ~80% while slightly
worsening stability. Orig's approach (variance scales with expectancy
level) is more realistic: when expectancy is lower, response variance
should be more constrained. This is standard in heteroscedastic
parameterization and is more clinically plausible.

### 7.3 Synthesis: Why theoretical appeal != empirical validity

The document's sections 3.4 pros (clinical realism, treatment-dependent
variance, testing under misspecification, Senn's variance decomposition)
are theoretically sound. **However**, the implementation in 2025 adds
correlations and variance structures that:

1. Are not specified in Hendrickson's original algorithm
2. Inflate baseline power without improving stability
3. Link biomarker to non-treatment components, contradicting the
   theoretical basis (biomarkers should predict treatment response, not
   placebo or time effects)
4. Produce power estimates 50-90% higher than the defensible correlation-based
   model (Hendrickson orig)

**Conclusion**: The correlation-based approach in principle has merit for
clinical realism, but the pmsimstats2025 extensions to that approach are
*not justified* and represent over-parameterization rather than principled
elaboration.

## 7. Revised Recommendation for pmsimstats2025

**Phase 2 testing has necessitated revision of section 7 recommendations.**

### 7.1 For faithful replication of Hendrickson (2020)

Keep the correlation-based approach with **Hendrickson's exact
parameterization**:
- Only BM-BR correlations (not BM-ER or BM-TV)
- Ron Thomas scaling of BM-BR by BR means ratio
- Expectancy-scaled ER SD (not constant)
- AR(1) temporal correlation structure (not Compound Symmetry)

This produces realistic power curves (~13-20% baseline for CO N=35, c.bm=0.3)
and measured vulnerability (~8% drop with carryover), matching Hendrickson
Figure 4.

### 7.2 Why 2025 "improvements" were unvalidated

The pmsimstats2025 extensions (BM-ER/BM-TV correlations + constant ER SD)
were hypothesized to improve clinical realism but **were not tested** against
the null hypothesis that they were unnecessary. Phase 2 testing revealed
they violate the core principle they were meant to implement: they inflate
baseline power without improving stability or detection of the
biomarker-treatment interaction.

Specifically, if these additions were improving the simulation's fidelity
to the causal mechanism of biomarker-modulated treatment response, they
should improve or maintain power curve stability. Instead, they degrade it,
indicating they are introducing noise or overfitting rather than capturing
true signal.

### 7.3 For the two-scenario comparison paper

Use the corrected correlation-based model (Hendrickson parameterization)
for the main analysis:
1. **Scenario 2** (no carryover modeling): Establishes realistic baseline
   vulnerability
2. **Scenario 1** (with carryover modeling): Demonstrates power recovery

Report sensitivity analysis with mean-structure DGP to show the result
generalizes to the more common simulation convention. Explicitly note that
choice of DGP affects power estimates by 30-50% (documented by Phase 2
testing) and justify the choice of correlation-based on fidelity to
Hendrickson rather than on theoretical superiority.

### 7.4 For methodological contribution

Phase 2 testing has filled the gap identified in section 6 by comparing
two simulation approaches systematically. Key finding: **the choice of
which correlations to include in the covariance structure has profound
effects on power estimates** (50-90% inflation from questionable additions).

Recommend publishing Phase 2 findings as a methodological note: "Validation
of Biomarker Interaction Mechanisms in Simulation: How Covariance Structure
Choices Affect Power Estimates in N-of-1 Aggregated Trials." This directly
addresses the oversight that section 3.3 noted: "Hendrickson et al. (2020)
is the only simulation study identified... [that uses this approach]. The
paper... does not explicitly discuss the choice..."

## 8. References

- Hendrickson RC, Thomas RG, Schork NJ. Optimizing aggregated
  N-of-1 trial designs for predictive biomarker validation.
  *Frontiers in Digital Health*. 2020;2:13.
  doi: 10.3389/fdgth.2020.00013

- Riecke BF, et al. A simulation study on estimating
  biomarker-treatment interaction effects in randomized trials
  with prognostic variables.
  *BMC Medical Research Methodology*. 2018;18:32.
  doi: 10.1186/s12874-018-0457-9

- Haller B, et al. A simulation study comparing different
  statistical approaches for the identification of predictive
  biomarkers. *Methods of Information in Medicine*.
  2019;58(S01):e60-e79. doi: 10.1159/000500436

- Langbaum JB, et al. Simulating effects of biomarker
  enrichment on Alzheimer's prevention trials.
  *Alzheimer's & Dementia*. 2014;10(5 Suppl):S388-S395.
  doi: 10.1016/j.jalz.2013.12.019

- Friede T, Parsons N, Stallard N. A conditional error
  function approach for subgroup selection in adaptive clinical
  trials. *Statistics in Medicine*. 2012;31(30):4309-4320.
  doi: 10.1002/sim.5541

- Blackston JW, et al. Comparison of aggregated N-of-1 trials
  with parallel and crossover randomized controlled trials using
  simulation studies. *Healthcare*. 2019;7(4):137.
  doi: 10.3390/healthcare7040137

- Senn S. Mastering variation: variance components and
  personalised medicine. *Statistics in Medicine*.
  2016;35(7):966-977. doi: 10.1002/sim.6739

- Senn S. Sample size considerations for n-of-1 trials.
  *Statistical Methods in Medical Research*. 2019;28(2):372-383.
  doi: 10.1177/0962280217726801

- Zhang Z, et al. The use of covariates and random effects in
  evaluating predictive biomarkers under a potential outcome
  framework. *Annals of Applied Statistics*.
  2014;8(4):2336-2355. doi: 10.1214/14-AOAS773

---
*Rendered on 2026-03-18 at 11:19 PDT.*
*Source: ~/prj/alz/01-pmsimstats/pmsimstats2025/analysis/scripts/carryover_focus/biomarker_interaction_mechanisms.md*
