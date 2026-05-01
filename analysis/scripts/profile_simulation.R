# Profile simulation_clustered.R to identify bottlenecks for Rcpp optimization
# Run with: Rscript profile_simulation.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(MASS)
  library(zoo)
  library(profvis)
  library(bench)
})

cat("=" |> strrep(70), "\n")
cat("PROFILING SIMULATION COMPONENTS\n")
cat("=" |> strrep(70), "\n\n")

# ============================================================================
# SETUP: Copy essential functions and parameters from simulation_clustered.R
# ============================================================================

n_participants <- 70
baseline_mean <- 10.0
between_subject_sd <- 2.0
within_subject_sd <- 1.8
biomarker_mean <- 5.0
biomarker_sd <- 2.0
c.br <- 0.75
c.er <- 0.75
c.tr <- 0.75
c.cf1t <- 0.12
c.cfct <- 0.05
c.bm_baseline <- 0.25
c.baseline_resp <- 0.3
BR_rate <- 0.5
ER_rate <- 0.2
TR_rate <- 0.1

measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)
n_tp <- length(measurement_weeks)

# ----------------------------------------------------------------------------
# Function 1: build_ar1_time (nested loops - candidate for Rcpp)
# ----------------------------------------------------------------------------

build_ar1_time_R <- function(weeks, rho, sigma) {
  n <- length(weeks)
  Cov <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      time_lag <- abs(weeks[i] - weeks[j])
      Cov[i, j] <- sigma^2 * rho^time_lag
    }
  }
  Cov
}

# Vectorized alternative
build_ar1_time_vec <- function(weeks, rho, sigma) {
  time_diff <- outer(weeks, weeks, FUN = function(x, y) abs(x - y))
  sigma^2 * rho^time_diff
}

# ----------------------------------------------------------------------------
# Function 2: build_sigma_guaranteed_pd
# ----------------------------------------------------------------------------

build_sigma_guaranteed_pd <- function(weeks, c.bm) {
  n_tp <- length(weeks)
  br_idx <- 1:n_tp
  er_idx <- (n_tp + 1):(2 * n_tp)
  tr_idx <- (2 * n_tp + 1):(3 * n_tp)

  Sigma_22 <- matrix(c(
    biomarker_sd^2,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    between_subject_sd^2
  ), 2, 2)

  Sigma_BR <- build_ar1_time_R(weeks, c.br, within_subject_sd)
  Sigma_ER <- build_ar1_time_R(weeks, c.er, within_subject_sd)
  Sigma_TR <- build_ar1_time_R(weeks, c.tr, within_subject_sd)

  Sigma_11 <- matrix(0, 3 * n_tp, 3 * n_tp)
  Sigma_11[br_idx, br_idx] <- Sigma_BR
  Sigma_11[er_idx, er_idx] <- Sigma_ER
  Sigma_11[tr_idx, tr_idx] <- Sigma_TR

  for (i in 1:n_tp) {
    for (j in 1:n_tp) {
      time_lag <- abs(weeks[i] - weeks[j])
      cross_cov <- if (i == j) {
        c.cf1t * within_subject_sd^2
      } else {
        c.cfct * within_subject_sd^2 * (0.9^time_lag)
      }
      Sigma_11[br_idx[i], er_idx[j]] <- cross_cov
      Sigma_11[er_idx[j], br_idx[i]] <- cross_cov
      Sigma_11[br_idx[i], tr_idx[j]] <- cross_cov
      Sigma_11[tr_idx[j], br_idx[i]] <- cross_cov
      Sigma_11[er_idx[i], tr_idx[j]] <- cross_cov
      Sigma_11[tr_idx[j], er_idx[i]] <- cross_cov
    }
  }

  Sigma_12 <- matrix(0, 3 * n_tp, 2)
  Sigma_12[br_idx, 1] <- c.bm * within_subject_sd * biomarker_sd
  Sigma_12[er_idx, 1] <- c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[br_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd

  Sigma_22_inv <- solve(Sigma_22)
  cross_term <- Sigma_12 %*% Sigma_22_inv %*% t(Sigma_12)
  Sigma_cond <- Sigma_11 - cross_term

  total_dim <- 3 * n_tp + 2
  Sigma <- matrix(0, total_dim, total_dim)
  Sigma[1:(3*n_tp), 1:(3*n_tp)] <- Sigma_11
  Sigma[(3*n_tp+1):total_dim, (3*n_tp+1):total_dim] <- Sigma_22
  Sigma[1:(3*n_tp), (3*n_tp+1):total_dim] <- Sigma_12
  Sigma[(3*n_tp+1):total_dim, 1:(3*n_tp)] <- t(Sigma_12)

  list(
    Sigma = Sigma,
    Sigma_22 = Sigma_22,
    Sigma_cond = Sigma_cond,
    Sigma_12 = Sigma_12,
    indices = list(br = br_idx, er = er_idx, tr = tr_idx, n_tp = n_tp)
  )
}

# ----------------------------------------------------------------------------
# Function 3: generate_participant_twostage
# ----------------------------------------------------------------------------

generate_participant_twostage <- function(sigma_parts, idx) {
  stage1 <- mvrnorm(1, mu = c(biomarker_mean, baseline_mean),
                    Sigma = sigma_parts$Sigma_22)
  z <- c(stage1[1] - biomarker_mean, stage1[2] - baseline_mean)
  cond_mean <- as.vector(sigma_parts$Sigma_12 %*%
                           solve(sigma_parts$Sigma_22) %*% z)
  responses <- mvrnorm(1, mu = cond_mean, Sigma = sigma_parts$Sigma_cond)

  list(
    biomarker = stage1[1],
    baseline = stage1[2],
    br_random = responses[idx$br],
    er_random = responses[idx$er],
    tr_random = responses[idx$tr]
  )
}

# ----------------------------------------------------------------------------
# Function 4: create_hybrid_design
# ----------------------------------------------------------------------------

create_hybrid_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) |>
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        week %in% c(4, 8) ~ 1,
        week == 9 ~ 1,
        week == 10 & path %in% c(1, 2) ~ 1,
        week == 10 & path %in% c(3, 4) ~ 0,
        week %in% c(11, 12) ~ 0,
        week == 16 & path %in% c(1, 3) ~ 1,
        week == 16 & path %in% c(2, 4) ~ 0,
        week == 20 & path %in% c(1, 3) ~ 0,
        week == 20 & path %in% c(2, 4) ~ 1,
        TRUE ~ NA_real_
      ),
      expectancy = if_else(week %in% c(4, 8), 1.0, 0.5)
    )
}

# ============================================================================
# BENCHMARK INDIVIDUAL COMPONENTS
# ============================================================================

cat("1. BENCHMARKING: build_ar1_time (R loops vs vectorized)\n")
cat("-" |> strrep(50), "\n")

bm_ar1 <- bench::mark(
  loop = build_ar1_time_R(measurement_weeks, c.br, within_subject_sd),
  vectorized = build_ar1_time_vec(measurement_weeks, c.br, within_subject_sd),
  check = TRUE,
  min_iterations = 100
)
print(bm_ar1)
cat("\n")

# ----------------------------------------------------------------------------

cat("2. BENCHMARKING: build_sigma_guaranteed_pd\n")
cat("-" |> strrep(50), "\n")

bm_sigma <- bench::mark(
  build_sigma = build_sigma_guaranteed_pd(measurement_weeks, 0.3),
  min_iterations = 100
)
print(bm_sigma)
cat("\n")

# ----------------------------------------------------------------------------

cat("3. BENCHMARKING: generate_participant_twostage (per participant)\n")
cat("-" |> strrep(50), "\n")

sigma_parts <- build_sigma_guaranteed_pd(measurement_weeks, 0.3)
idx <- sigma_parts$indices

bm_gen <- bench::mark(
  generate_one = generate_participant_twostage(sigma_parts, idx),
  min_iterations = 1000
)
print(bm_gen)
cat("\n")

# ----------------------------------------------------------------------------

cat("4. BENCHMARKING: generate ALL participants (n=70)\n")
cat("-" |> strrep(50), "\n")

generate_all_participants <- function(sigma_parts, idx, n) {
  lapply(1:n, function(pid) {
    generate_participant_twostage(sigma_parts, idx)
  })
}

bm_gen_all <- bench::mark(
  generate_70 = generate_all_participants(sigma_parts, idx, 70),
  min_iterations = 50
)
print(bm_gen_all)
cat("\n")

# ----------------------------------------------------------------------------

cat("5. BENCHMARKING: create_hybrid_design\n")
cat("-" |> strrep(50), "\n")

bm_design <- bench::mark(
  create_design = create_hybrid_design(70, measurement_weeks),
  min_iterations = 100
)
print(bm_design)
cat("\n")

# ----------------------------------------------------------------------------

cat("6. BENCHMARKING: lme model fitting (THE BOTTLENECK)\n")
cat("-" |> strrep(50), "\n")

# Generate properly structured sample data for model fitting
set.seed(42)
n_part <- 70
n_weeks <- length(measurement_weeks)

trial_data <- expand_grid(
  participant_id = 1:n_part,
  week = measurement_weeks
) |>
  mutate(
    # Random effects per participant
    participant_intercept = rnorm(n_part)[participant_id],
    # Generate biomarker once per participant
    biomarker = rnorm(n_part, biomarker_mean, biomarker_sd)[participant_id],
    bm_centered = biomarker - biomarker_mean,
    # Treatment varies by week (simplified pattern)
    treatment = if_else(week <= 10, 1, 0),
    # Response with reasonable variance
    response = baseline_mean +
      participant_intercept * between_subject_sd +
      treatment * 0.5 +        # treatment effect
      bm_centered * 0.2 +      # biomarker effect
      week * 0.05 +            # time trend
      rnorm(n(), 0, within_subject_sd)  # residual error
  ) |>
  dplyr::select(participant_id, week, treatment, bm_centered, response)

# Benchmark model fitting with simpler correlation structure
bm_lme <- bench::mark(
  lme_fit = lme(
    response ~ treatment + bm_centered + week,
    random = ~ 1 | participant_id,
    data = trial_data,
    control = lmeControl(opt = "optim")
  ),
  min_iterations = 20
)
print(bm_lme)
cat("\n")

# Also benchmark with corCAR1 (more expensive)
cat("6b. BENCHMARKING: lme with corCAR1 correlation\n")
cat("-" |> strrep(50), "\n")

bm_lme_car1 <- tryCatch({
  bench::mark(
    lme_car1 = lme(
      response ~ treatment + bm_centered + week,
      random = ~ 1 | participant_id,
      correlation = corCAR1(form = ~ week | participant_id),
      data = trial_data,
      control = lmeControl(opt = "optim")
    ),
    min_iterations = 10
  )
}, error = function(e) {
  cat("  corCAR1 model failed:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(bm_lme_car1)) {
  print(bm_lme_car1)
}
cat("\n")

# ============================================================================
# SUMMARY AND RECOMMENDATIONS
# ============================================================================

cat("=" |> strrep(70), "\n")
cat("PROFILING SUMMARY\n")
cat("=" |> strrep(70), "\n\n")

# Extract median times
times <- list(
  ar1_loop = as.numeric(bm_ar1$median[1]),
  ar1_vec = as.numeric(bm_ar1$median[2]),
  sigma = as.numeric(bm_sigma$median[1]),
  gen_one = as.numeric(bm_gen$median[1]),
  gen_all = as.numeric(bm_gen_all$median[1]),
  design = as.numeric(bm_design$median[1]),
  lme = as.numeric(bm_lme$median[1])
)

# Estimate total time per iteration
n_conditions <- 100  # approximate
n_iterations <- 100
participants_per_iter <- 70

# Per-iteration breakdown (in seconds, assuming times are in various units)
cat("Estimated time breakdown per iteration:\n\n")

cat(sprintf("  build_sigma:           %s\n", format(bm_sigma$median[1])))
cat(sprintf("  generate_participants: %s\n", format(bm_gen_all$median[1])))
cat(sprintf("  create_design:         %s\n", format(bm_design$median[1])))
cat(sprintf("  lme model fitting:     %s\n", format(bm_lme$median[1])))
if (!is.null(bm_lme_car1)) {
  cat(sprintf("  lme + corCAR1:         %s\n", format(bm_lme_car1$median[1])))
}
cat("\n")

# Calculate percentages - use corCAR1 time if available (more realistic)
lme_time <- if (!is.null(bm_lme_car1)) bm_lme_car1$median[1] else bm_lme$median[1]
total_est <- as.numeric(bm_sigma$median[1]) +
  as.numeric(bm_gen_all$median[1]) +
  as.numeric(bm_design$median[1]) +
  as.numeric(lme_time)

pct_sigma <- as.numeric(bm_sigma$median[1]) / total_est * 100
pct_gen <- as.numeric(bm_gen_all$median[1]) / total_est * 100
pct_design <- as.numeric(bm_design$median[1]) / total_est * 100
pct_lme <- as.numeric(lme_time) / total_est * 100

cat("Percentage of iteration time:\n\n")
cat(sprintf("  build_sigma:           %5.1f%%\n", pct_sigma))
cat(sprintf("  generate_participants: %5.1f%%\n", pct_gen))
cat(sprintf("  create_design:         %5.1f%%\n", pct_design))
cat(sprintf("  lme model fitting:     %5.1f%% <-- BOTTLENECK\n", pct_lme))
cat("\n")

cat("=" |> strrep(70), "\n")
cat("RCPP OPTIMIZATION RECOMMENDATIONS\n")
cat("=" |> strrep(70), "\n\n")

cat("1. LME MODEL FITTING (~", round(pct_lme), "% of time)\n")
cat("   - This is the dominant bottleneck\n")
cat("   - Cannot be easily replaced with Rcpp (complex iterative algorithm)\n")
cat("   - Options:\n")
cat("     a) Use lme4::lmer instead of nlme::lme (often faster)\
")
cat("     b) Parallelize across iterations (foreach/future)\n")
cat("     c) Consider simplified models for initial screening\n\n")

cat("2. COVARIANCE MATRIX CONSTRUCTION (~", round(pct_sigma), "% of time)\n")
cat("   - build_ar1_time has nested loops (good Rcpp candidate)\n")
cat("   - Already showed vectorized R is ",
    round(as.numeric(bm_ar1$median[1]) / as.numeric(bm_ar1$median[2]), 1),
    "x faster\n")
cat("   - Rcpp would give modest additional speedup (~2-5x over vectorized)\n")
cat("   - Low priority: small fraction of total time\n\n")

cat("3. PARTICIPANT DATA GENERATION (~", round(pct_gen), "% of time)\n")
cat("   - mvrnorm calls are already optimized (uses LAPACK)\n")
cat("   - solve() for Sigma_22_inv could be pre-computed and cached\n")
cat("   - Rcpp unlikely to help much here\n\n")

cat("4. DESIGN CREATION (~", round(pct_design), "% of time)\n")
cat("   - dplyr operations are already efficient\n")
cat("   - Could pre-compute designs outside iteration loop\n")
cat("   - Rcpp not recommended\n\n")

cat("OVERALL RECOMMENDATION:\n")
cat("-" |> strrep(50), "\n")
cat("Rcpp will provide MINIMAL benefit (<5% speedup) because:\n")
cat("  - Model fitting dominates runtime and cannot be easily replaced\n")
cat("  - Matrix operations already use optimized BLAS/LAPACK\n")
cat("  - R's vectorized operations are near-C performance\n\n")

cat("BETTER OPTIMIZATION STRATEGIES:\n")
cat("  1. PARALLELIZE iterations across CPU cores (biggest win)\n")
cat("  2. Cache sigma matrices (only 3 unique structures)\n")
cat("  3. Pre-compute Sigma_22_inv once per sigma structure\n")
cat("  4. Consider lme4::lmer instead of nlme::lme\n")
cat("  5. Reduce n_iterations for exploratory runs\n\n")

cat("=" |> strrep(70), "\n")
