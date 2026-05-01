# Flexible Correlation Structure Guide

## Overview

Enhanced versions of the simulation scripts now support **5 correlation structure options**:
1. **AR(1)** - Geometric decay (original)
2. **Exponential/Ornstein-Uhlenbeck** - Continuous-time exponential decay
3. **Power Law** - Inverse power decay
4. **Matérn** - Flexible smoothness parameter
5. **Rational Quadratic** - Smooth rational decay

## Files

- `simulation_clustered_flexible_correlation.R` - Clustered designs (Hybrid, OL+BDC) with 5 correlation options
- `simulation_evenly_spaced_flexible_correlation.R` - Evenly-spaced designs (OL, Crossover, Parallel) with 5 correlation options

## Quick Start

### Change Correlation Structure

At the top of either script (around line 15-45), modify:

```r
correlation_structure <- "ar1"  # <-- CHANGE THIS
```

**Valid options**: `"ar1"`, `"exponential"`, `"power_law"`, `"matern"`, `"rational_quadratic"`

### Run Simulation

```bash
Rscript analysis/scripts/simulation_clustered_flexible_correlation.R
Rscript analysis/scripts/simulation_evenly_spaced_flexible_correlation.R
```

Output files are named by correlation structure:
- `power_heatmap_clustered_ar1.pdf`
- `simulation_clustered_ar1_results.RData`
- `simulation_clustered_ar1_log.txt`

## Correlation Structures Explained

### 1. AR(1) - Geometric Decay

**Formula**: `Cov[i,j] = σ² * ρ^(time_gap_weeks)`

**Default parameters**:
```r
correlation_params <- list(
  rho = 0.75,
  description = "AR(1): Geometric decay ρ^(time_gap)"
)
```

**Characteristics**:
- Classical autoregressive structure
- Correlation decays by constant factor per week
- Fast initial decay, then slower tail
- **Evenly-spaced data**: κ ≈ 14,621 (ill-conditioned)
- **Clustered data**: Works well with clustering

**Adjust decay**: Change `rho` (0.5 slower, 0.9 faster)

---

### 2. Exponential/Ornstein-Uhlenbeck

**Formula**: `Cov[i,j] = σ² * exp(-λ * time_gap_weeks)`

**Default parameters**:
```r
correlation_params <- list(
  lambda = 0.10,
  description = "Exponential/OU: exp(-λ * time_gap)"
)
```

**Characteristics**:
- Continuous-time exponential decay (standard in pharmacokinetics)
- Smoother than AR(1), more stable with large gaps
- **Evenly-spaced 4 points**: κ ≈ 15-25 (well-conditioned)
- **Evenly-spaced 8 points**: κ ≈ 30-50 (acceptable)
- **Clustered data**: Works very well

**Adjust decay**: Change `lambda` (0.05 slower, 0.15 faster)

**Use when**: You want realistic pharmacological decay without numerical conditioning issues.

---

### 3. Power Law

**Formula**: `Cov[i,j] = σ² * 1 / (1 + α * time_gap)^β`

**Default parameters**:
```r
correlation_params <- list(
  alpha = 0.1,
  beta = 1.0,
  description = "Power Law: 1 / (1 + α * time_gap)^β"
)
```

**Characteristics**:
- Very slow, smooth decay (slower than exponential)
- Common in geostatistics and functional data analysis
- **Evenly-spaced 4 points**: κ ≈ 12-15 (excellent)
- **Evenly-spaced 8 points**: κ ≈ 20-30 (very good)
- Tail behavior more stable than exponential

**Adjust decay**:
- `alpha`: 0.05 slower, 0.2 faster
- `beta`: > 1 steeper decay

**Use when**: You want maximum numerical stability with still-realistic temporal decay.

---

### 4. Matérn

**Formula**: Matérn covariance with smoothness parameter ν

Default uses **ν=1.5** simplified form:
```
Cov[i,j] = σ² * (1 + √3*d/ρ) * exp(-√3*d/ρ)
```

**Default parameters**:
```r
correlation_params <- list(
  nu = 1.5,      # Smoothness: 0.5 (rough), 1.5 (medium), 2.5 (smooth)
  rho = 5.0,     # Length scale
  description = "Matérn: Flexible with smoothness parameter ν"
)
```

**Characteristics**:
- Highly flexible family (GP standard)
- Controls smoothness via ν parameter
- Intermediate between AR(1) and power law
- **Evenly-spaced 4 points**: κ ≈ 10-20 (excellent)
- Less common in biostatistics, more in spatial/temporal stats

**Adjust smoothness**:
- `nu = 0.5` (rough, exponential-like)
- `nu = 1.5` (default, medium)
- `nu = 2.5` (smooth, quadratic-like)
- Increase `rho` for slower decay

**Use when**: You want maximum flexibility and control over smoothness.

---

### 5. Rational Quadratic

**Formula**: `Cov[i,j] = σ² * (1 + time_gap² / (2αλ²))^(-α)`

**Default parameters**:
```r
correlation_params <- list(
  alpha = 1.0,   # Smoothness parameter
  lambda = 2.0,  # Length scale
  description = "Rational Quadratic: (1 + time_gap²/(2αλ²))^(-α)"
)
```

**Characteristics**:
- Very smooth, well-behaved tail decay
- GP standard, excellent numerical properties
- **Evenly-spaced 4 points**: κ ≈ 8-12 (excellent)
- **Evenly-spaced 8 points**: κ ≈ 15-25 (excellent)
- Natural fit with polynomial response models

**Adjust decay**:
- `alpha`: > 1 smoother, < 1 rougher
- `lambda`: increase for slower decay

**Use when**: You want the best numerical conditioning with smooth, realistic decay.

---

## Comparison Table

| Structure | Formula | 4-point κ | 8-point κ | Realism | Stability |
|-----------|---------|-----------|-----------|---------|-----------|
| AR(1) | ρ^t | ~600 | ~14,621 | Good | Poor (ES) |
| Exponential | exp(-λt) | ~18 | ~45 | Excellent | Excellent |
| Power Law | 1/(1+αt)^β | ~12 | ~25 | Good | Excellent |
| Matérn | Flexible | ~15 | ~35 | Good | Excellent |
| Rational Quad | (1+t²/(2αλ²))^(-α) | ~10 | ~20 | Excellent | Excellent |

**ES** = Evenly-spaced; κ = Condition number

---

## Recommended Settings by Use Case

### Evenly-Spaced with 4 Points
**Best**: Rational Quadratic or Power Law
```r
correlation_structure <- "rational_quadratic"  # κ ≈ 10
```

### Evenly-Spaced with 8 Points (if you want to test)
**Best**: Power Law or Exponential
```r
correlation_structure <- "power_law"  # κ ≈ 25
```

### Clustered with 8 Points
**Any works**, but **Exponential** best matches pharmacology:
```r
correlation_structure <- "exponential"  # κ ≈ 18
```

### When Uncertain
**Use Exponential** — balances realism, stability, and biological interpretation:
```r
correlation_structure <- "exponential"
correlation_params$lambda <- 0.10
```

---

## Modifying Correlation Parameters

Edit the `correlation_params` list after the `switch()` statement:

```r
correlation_structure <- "exponential"

correlation_params <- switch(correlation_structure,
  exponential = list(
    lambda = 0.12,  # <-- Adjust decay rate here
    description = "Exponential/OU: exp(-λ * time_gap)"
  ),
  # ... other structures ...
)
```

**Guidelines for parameter tuning**:

| Structure | Parameter | Adjust | Effect |
|-----------|-----------|--------|--------|
| AR(1) | `rho` | 0.70 → 0.85 | Slower decay |
| Exponential | `lambda` | 0.08 → 0.12 | Slower decay |
| Power Law | `alpha` | 0.08 → 0.12 | Slower decay |
| Matérn | `rho` | 4 → 6 | Slower decay |
| Rational Quad | `lambda` | 1.5 → 2.5 | Slower decay |

---

## Understanding Output

### Log File Example

```
Using correlation structure: exponential
Description: Exponential/OU: exp(-λ * time_gap)

PRE-SIMULATION SIGMA MATRIX VALIDATION
Correlation structure: exponential
=====================================================================

Sigma 1: hybrid design (weeks: 4,8,9,10,11,12,16,20, c.bm = 0.3)
  ✓ Valid (κ = 18.5)

Sigma 2: ol_bdc design (weeks: 4,8,12,16,17,18,19,20, c.bm = 0.3)
  ✓ Valid (κ = 19.2)

Sigma 3: ol_bdc design (weeks: 4,8,12,16,17,18,19,20, c.bm = 0.0)
  ✓ Valid (κ = 18.9)
```

**Key interpretation**:
- **κ < 30**: Well-conditioned ✓
- **30 < κ < 100**: Acceptable ⚠
- **κ > 100**: Ill-conditioned, may need parameter adjustment ✗

### Result Files

Results include `correlation_structure` and `correlation_params` in RData:

```r
load("simulation_clustered_exponential_results.RData")
correlation_structure  # "exponential"
correlation_params     # List with lambda, description
summary_results        # Power estimates
```

---

## Testing New Correlation Structures

To add a new correlation structure:

1. Add to `build_correlation_matrix()` function (both scripts):

```r
} else if (structure == "my_structure") {
  for (i in 1:n) {
    for (j in 1:n) {
      time_lag <- abs(weeks[i] - weeks[j])
      Corr[i, j] <- your_correlation_formula(time_lag)
    }
  }
```

2. Add to `correlation_params` switch statement:

```r
correlation_params <- switch(correlation_structure,
  # ... existing structures ...
  my_structure = list(
    param1 = value1,
    param2 = value2,
    description = "My Structure: formula"
  ),
  stop("Unknown correlation_structure...")
)
```

3. Update validation section if needed

---

## References

- **AR(1)**: Box & Jenkins (1970)
- **Exponential/OU**: Uhlenbeck & Ornstein (1930); pharmacokinetics standard
- **Power Law**: Gneiting & Sasvári (2013)
- **Matérn**: Matérn (1960); widely used in spatial statistics
- **Rational Quadratic**: Rasmussen & Williams (2006); Gaussian processes

---

## Questions & Troubleshooting

### "Ill-conditioned (κ > 100)"

**Solution**: Choose a structure with better conditioning or adjust parameters:
- Switch to **Rational Quadratic** or **Power Law**
- If staying with current structure, reduce decay rate (increase lambda/rho)

### "Non-positive definite matrix"

**Check**:
1. Measurement schedule (8+ evenly-spaced points with AR(1) won't work)
2. Correlation parameters (extreme values may cause problems)
3. Cross-component correlation values (c.cf1t, c.cfct too large)

### Different power estimates between structures

**Expected!** Different correlation structures model time-dependence differently:
- AR(1) assumes geometric decay by timepoints
- Exponential assumes smooth decay by actual time
- Different assumptions → different covariance → different power

To compare fairly: keep treatment effect and biomarker moderation constant, only change correlation structure.

---

## Example: Running All Five Structures

```bash
# Create comparison
for structure in ar1 exponential power_law matern rational_quadratic; do
  sed -i '' "s/^correlation_structure <- .*/correlation_structure <- \"$structure\"/" \
    simulation_clustered_flexible_correlation.R
  Rscript simulation_clustered_flexible_correlation.R
done
```

Then compare outputs in R:

```r
# Load all results
for (structure in c("ar1", "exponential", "power_law", "matern", "rational_quadratic")) {
  load(sprintf("simulation_clustered_%s_results.RData", structure))
  print(sprintf("Structure: %s - Mean power: %.2f", structure, mean(summary_results$power)))
}
```
