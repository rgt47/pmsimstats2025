# Simplified N-of-1 Trial Simulation - Constant Effect Model
# Maximum simplification for learning

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(MASS)
  library(conflicted)
})

# Resolve conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(lmerTest::lmer)

# =============================================================================
# PARAMETERS
# =============================================================================

n_participants <- 70
n_iterations <- 20

# Three-factor response model - RATE-BASED (points per week)
BR_rate <- 0.5                 # Biological Response: drug improvement rate
ER_rate <- 0.2                 # Expectancy Response: placebo improvement rate
TR_rate <- 0.1                 # Time-variant Response: natural improvement rate
treatment_effect <- BR_rate    # Alias for display purposes

# Biomarker moderation of treatment effect
# Higher biomarker → stronger treatment response
# This creates the treatment × biomarker interaction
# NOTE: biomarker_moderation is now set in the param_grid, not here

baseline_mean <- 10.0          # Mean baseline response
between_subject_sd <- 2.0      # SD of participant random effects
within_subject_sd <- 1.8       # SD of measurement noise
biomarker_mean <- 5.0          # Mean biomarker value
biomarker_sd <- 2.0            # SD of biomarker

# Correlation parameters (from Hendrickson)
# Autocorrelations within each response type
c.br <- 0.8                    # BR autocorrelation across timepoints
c.er <- 0.8                    # ER autocorrelation across timepoints
c.tr <- 0.8                    # TR autocorrelation across timepoints

# Cross-correlations between response types
c.cf1t <- 0.2                  # Same-time cross-correlation (BR-ER, BR-TR, ER-TR)
c.cfct <- 0.1                  # Different-time cross-correlation

# Biomarker and baseline correlations
c.bm_baseline <- 0.3           # Biomarker-baseline correlation
c.baseline_resp <- 0.4         # Baseline-response correlation

# Allowed correlation values (finite grid for consistent results)
allowed_correlations <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

# Parameter grid - what we're testing
# Note: carryover only applies to hybrid design (parallel has no treatment switches)
# Include biomarker_moderation = 0 to evaluate Type I error (size of test)
param_grid <- bind_rows(
  # Hybrid design with carryover variations
  expand_grid(
    design = "hybrid",
    biomarker_moderation = c(0, 0.25, 0.35, 0.45),
    biomarker_correlation = c(0.3),
    carryover_decay_rate = c(0, 0.5)
  ),
  # Parallel design (carryover = 0, not applicable)
  expand_grid(
    design = "parallel",
    biomarker_moderation = c(0, 0.25, 0.35, 0.45),
    biomarker_correlation = c(0.3),
    carryover_decay_rate = c(0)
  )
)

# =============================================================================
# BUILD SIGMA WITH GUARANTEED PD (Time-Based AR(1))
# =============================================================================

# Build sigma using independent construction with automatic scaling
# Structure: responses (24) conditioned on participant vars (2)
build_sigma_guaranteed_pd <- function(weeks, c.bm, params) {
  n_tp <- length(weeks)

  # Index ranges for final 26x26 (for compatibility)
  br_idx <- 1:n_tp
  er_idx <- (n_tp + 1):(2 * n_tp)
  tr_idx <- (2 * n_tp + 1):(3 * n_tp)
  bm_idx <- 3 * n_tp + 1
  bl_idx <- 3 * n_tp + 2

  # =========================================================================
  # STAGE 1: Build Σ₂₂ (2×2) - Biomarker & Baseline
  # =========================================================================
  # Always PD for |correlation| < 1

  Sigma_22 <- matrix(c(
    biomarker_sd^2,                                    c.bm_baseline * biomarker_sd * between_subject_sd,
    c.bm_baseline * biomarker_sd * between_subject_sd, between_subject_sd^2
  ), 2, 2)

  Sigma_22_inv <- solve(Sigma_22)

  # =========================================================================
  # STAGE 2: Build Σ₁₁ (24×24) - Responses with Time-Based AR(1)
  # =========================================================================
  # Guaranteed PD by construction

  # Helper: Build AR(1) covariance block based on actual time lags
  build_ar1_time <- function(weeks, rho, sigma) {
    n <- length(weeks)
    Cov <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        time_lag <- abs(weeks[i] - weeks[j])
        Cov[i, j] <- sigma^2 * rho^time_lag
      }
    }
    return(Cov)
  }

  # Build diagonal blocks (within each response type)
  Sigma_BR <- build_ar1_time(weeks, c.br, within_subject_sd)
  Sigma_ER <- build_ar1_time(weeks, c.er, within_subject_sd)
  Sigma_TR <- build_ar1_time(weeks, c.tr, within_subject_sd)

  # Initialize 24x24 with diagonal blocks
  Sigma_11 <- matrix(0, 3 * n_tp, 3 * n_tp)
  Sigma_11[br_idx, br_idx] <- Sigma_BR
  Sigma_11[er_idx, er_idx] <- Sigma_ER
  Sigma_11[tr_idx, tr_idx] <- Sigma_TR

  # Add cross-correlations between response types
  # Use time-based decay for cross-correlations too
  for (i in 1:n_tp) {
    for (j in 1:n_tp) {
      time_lag <- abs(weeks[i] - weeks[j])

      if (i == j) {
        # Same timepoint: use c.cf1t
        cross_cov <- c.cf1t * within_subject_sd^2
      } else {
        # Different timepoint: use c.cfct with time decay
        cross_cov <- c.cfct * within_subject_sd^2 * (0.9^time_lag)
      }

      # BR-ER
      Sigma_11[br_idx[i], er_idx[j]] <- cross_cov
      Sigma_11[er_idx[j], br_idx[i]] <- cross_cov

      # BR-TR
      Sigma_11[br_idx[i], tr_idx[j]] <- cross_cov
      Sigma_11[tr_idx[j], br_idx[i]] <- cross_cov

      # ER-TR
      Sigma_11[er_idx[i], tr_idx[j]] <- cross_cov
      Sigma_11[tr_idx[j], er_idx[i]] <- cross_cov
    }
  }

  # =========================================================================
  # STAGE 3: Build Σ₁₂ (24×2) - Cross-covariance (regression coefficients)
  # =========================================================================

  Sigma_12 <- matrix(0, 3 * n_tp, 2)

  # Column 1: Biomarker effects on responses
  # BR strongly correlated with biomarker (this is what we're testing)
  Sigma_12[br_idx, 1] <- c.bm * within_subject_sd * biomarker_sd
  # ER and TR weakly correlated with biomarker
  Sigma_12[er_idx, 1] <- c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- c.bm * 0.5 * within_subject_sd * biomarker_sd

  # Column 2: Baseline effects on responses
  Sigma_12[br_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- c.baseline_resp * within_subject_sd * between_subject_sd

  # =========================================================================
  # STAGE 4: Check Schur Complement & Snap to Grid if Needed
  # =========================================================================

  # Function to find largest allowed correlation that keeps matrix PD
  find_valid_correlation <- function(Sigma_11, Sigma_22_inv, c.bm_requested,
                                     allowed_vals, sigma_resp, sigma_bm, sigma_bl,
                                     br_idx, er_idx, tr_idx) {
    # Try correlations from requested down to 0
    valid_vals <- allowed_vals[allowed_vals <= c.bm_requested]
    valid_vals <- sort(valid_vals, decreasing = TRUE)

    for (c.bm_try in valid_vals) {
      # Rebuild Sigma_12 with this correlation
      Sigma_12_try <- matrix(0, nrow(Sigma_11), 2)
      Sigma_12_try[br_idx, 1] <- c.bm_try * sigma_resp * sigma_bm
      Sigma_12_try[er_idx, 1] <- c.bm_try * 0.5 * sigma_resp * sigma_bm
      Sigma_12_try[tr_idx, 1] <- c.bm_try * 0.5 * sigma_resp * sigma_bm
      Sigma_12_try[br_idx, 2] <- c.baseline_resp * sigma_resp * sigma_bl
      Sigma_12_try[er_idx, 2] <- c.baseline_resp * sigma_resp * sigma_bl
      Sigma_12_try[tr_idx, 2] <- c.baseline_resp * sigma_resp * sigma_bl

      # Check if PD
      cross_term <- Sigma_12_try %*% Sigma_22_inv %*% t(Sigma_12_try)
      Sigma_cond_try <- Sigma_11 - cross_term
      min_eig <- min(eigen(Sigma_cond_try, only.values = TRUE)$values)

      if (min_eig > 1e-6) {
        return(list(c.bm = c.bm_try, Sigma_12 = Sigma_12_try, Sigma_cond = Sigma_cond_try))
      }
    }

    # Fallback to 0 correlation
    Sigma_12_try <- matrix(0, nrow(Sigma_11), 2)
    Sigma_12_try[br_idx, 2] <- c.baseline_resp * sigma_resp * sigma_bl
    Sigma_12_try[er_idx, 2] <- c.baseline_resp * sigma_resp * sigma_bl
    Sigma_12_try[tr_idx, 2] <- c.baseline_resp * sigma_resp * sigma_bl
    cross_term <- Sigma_12_try %*% Sigma_22_inv %*% t(Sigma_12_try)
    Sigma_cond_try <- Sigma_11 - cross_term

    return(list(c.bm = 0, Sigma_12 = Sigma_12_try, Sigma_cond = Sigma_cond_try))
  }

  # Find valid correlation on the grid
  result <- find_valid_correlation(
    Sigma_11, Sigma_22_inv, c.bm,
    allowed_correlations, within_subject_sd, biomarker_sd, between_subject_sd,
    br_idx, er_idx, tr_idx
  )

  Sigma_12 <- result$Sigma_12
  Sigma_cond <- result$Sigma_cond
  effective_c.bm <- result$c.bm

  if (effective_c.bm < c.bm) {
    cat(sprintf("  Snapped biomarker correlation: %.2f → %.2f (to ensure PD)\n",
                c.bm, effective_c.bm))
  }

  # =========================================================================
  # Return partitioned result
  # =========================================================================

  return(list(
    Sigma_11 = Sigma_11,
    Sigma_22 = Sigma_22,
    Sigma_12 = Sigma_12,
    Sigma_cond = Sigma_cond,
    Sigma_22_inv = Sigma_22_inv,
    requested_c.bm = c.bm,
    effective_c.bm = effective_c.bm,
    indices = list(br = br_idx, er = er_idx, tr = tr_idx, bm = bm_idx, bl = bl_idx)
  ))
}

# =============================================================================
# TWO-STAGE DATA GENERATION (2x2 + 24x24)
# =============================================================================

# Generate one participant's data using two-stage approach
# Uses pre-computed partitioned sigma from build_sigma_guaranteed_pd
generate_participant_twostage <- function(sigma_parts, idx) {
  # Stage 1: Generate participant variables (biomarker, baseline) from 2x2
  x2 <- mvrnorm(1, mu = c(0, 0), Sigma = sigma_parts$Sigma_22)

  # Stage 2: Compute conditional mean given participant vars
  # Conditional mean: μ₁|₂ = Σ₁₂ Σ₂₂⁻¹ x₂
  mu_cond <- sigma_parts$Sigma_12 %*% sigma_parts$Sigma_22_inv %*% x2

  # Generate responses from conditional distribution
  # Sigma_cond is pre-computed and guaranteed PD
  x1 <- mvrnorm(1, mu = as.vector(mu_cond), Sigma = sigma_parts$Sigma_cond)

  # Extract components
  n_tp <- length(idx$br)
  list(
    biomarker = x2[1] + biomarker_mean,
    baseline = x2[2] + baseline_mean,
    br_random = x1[1:n_tp],
    er_random = x1[(n_tp+1):(2*n_tp)],
    tr_random = x1[(2*n_tp+1):(3*n_tp)]
  )
}

# =============================================================================
# VERIFICATION: Compare two-stage vs direct 26x26
# =============================================================================

verify_twostage_equivalence <- function(Sigma, idx, n_samples = 10000) {
  cat("Verifying two-stage equivalence with", n_samples, "samples...\n")

  # Partition sigma
  partitioned <- partition_sigma(Sigma, idx)

  # Generate samples using direct method
  set.seed(999)
  samples_direct <- mvrnorm(n_samples, mu = rep(0, nrow(Sigma)), Sigma = Sigma)

  # Generate samples using two-stage method
  set.seed(888)
  samples_twostage <- matrix(0, n_samples, nrow(Sigma))
  resp_idx <- c(idx$br, idx$er, idx$tr)
  part_idx <- c(idx$bm, idx$bl)

  for (i in 1:n_samples) {
    result <- generate_participant_twostage(partitioned, idx)
    samples_twostage[i, resp_idx] <- c(result$br_random, result$er_random, result$tr_random)
    samples_twostage[i, part_idx] <- c(result$biomarker - biomarker_mean,
                                        result$baseline - baseline_mean)
  }

  # Compare covariance matrices
  cov_direct <- cov(samples_direct)
  cov_twostage <- cov(samples_twostage)

  max_cov_diff <- max(abs(cov_direct - cov_twostage))
  mean_cov_diff <- mean(abs(cov_direct - cov_twostage))

  # Compare means (should both be ~0)
  max_mean_diff <- max(abs(colMeans(samples_direct) - colMeans(samples_twostage)))

  cat(sprintf("  Max covariance difference: %.6f\n", max_cov_diff))
  cat(sprintf("  Mean covariance difference: %.6f\n", mean_cov_diff))
  cat(sprintf("  Max mean difference: %.6f\n", max_mean_diff))

  # Check if differences are within sampling error
  # Expected sampling error for covariance ~= 2*sigma^2/sqrt(n)
  expected_error <- 2 * max(diag(Sigma))^2 / sqrt(n_samples)

  if (max_cov_diff < expected_error * 3) {
    cat("  ✓ PASSED: Differences within expected sampling error\n")
    return(TRUE)
  } else {
    cat("  ✗ FAILED: Differences exceed expected sampling error\n")
    return(FALSE)
  }
}

# =============================================================================
# DESIGN STRUCTURES
# =============================================================================

# Measurement schedule (same for all designs)
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

# -----------------------------------------------------------------------------
# HYBRID DESIGN: Open-label run-in + blinded crossover
# -----------------------------------------------------------------------------
create_hybrid_design <- function(n_participants, measurement_weeks) {
  # Randomize participants to 4 paths (balanced)
  path_assignment <- sample(rep(1:4, length.out = n_participants))

  # Create design matrix
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],

      # Treatment assignment based on path and week
      # Weeks 4, 8: Open-label, all on active
      # Week 9: Blinded, all on active
      # Week 10: Blinded, randomized (paths 1,2 = active; paths 3,4 = placebo)
      # Weeks 11, 12: Blinded, all on placebo
      # Week 16: Blinded crossover
      # Week 20: Blinded crossover
      treatment = case_when(
        week %in% c(4, 8) ~ 1,                              # Open-label: all active
        week == 9 ~ 1,                                      # Blinded: all active
        week == 10 & path %in% c(1, 2) ~ 1,                 # Randomized: paths 1,2 active
        week == 10 & path %in% c(3, 4) ~ 0,                 # Randomized: paths 3,4 placebo
        week %in% c(11, 12) ~ 0,                            # All on placebo
        week == 16 & path %in% c(1, 3) ~ 1,                 # Crossover period
        week == 16 & path %in% c(2, 4) ~ 0,
        week == 20 & path %in% c(1, 3) ~ 0,                 # Crossover period
        week == 20 & path %in% c(2, 4) ~ 1,
        TRUE ~ NA_real_
      ),

      # Expectancy: 1.0 = open-label (know they're on drug), 0.5 = blinded
      expectancy = if_else(week %in% c(4, 8), 1.0, 0.5)
    )
}

# -----------------------------------------------------------------------------
# PARALLEL DESIGN: Randomized to treatment or placebo for entire study
# -----------------------------------------------------------------------------
create_parallel_design <- function(n_participants, measurement_weeks) {
  # Randomize participants to treatment (1) or placebo (0) - balanced
  treatment_assignment <- sample(rep(0:1, length.out = n_participants))

  # Create design matrix
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = treatment_assignment[participant_id] + 1,  # Path 1 = placebo, Path 2 = treatment

      # Treatment: same throughout study based on randomization
      treatment = treatment_assignment[participant_id],

      # Expectancy: 0.5 throughout (all blinded)
      expectancy = 0.5
    )
}

# =============================================================================
# RUN SIMULATION
# =============================================================================

cat("Running simulation with constant effect model...\n")
cat("Treatment effect =", treatment_effect, "\n")
cat("Participants =", n_participants, "\n")
cat("Iterations per condition =", n_iterations, "\n\n")

results <- tibble()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf("Condition %d: design=%s, bm_mod=%.2f, carryover=%.2f\n",
              i, params$design, params$biomarker_moderation, params$carryover_decay_rate))

  # Run iterations
  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    # Create design with fresh path randomization for this iteration
    if (params$design == "hybrid") {
      trial_design <- create_hybrid_design(n_participants, measurement_weeks)
    } else if (params$design == "parallel") {
      trial_design <- create_parallel_design(n_participants, measurement_weeks)
    }

    n_timepoints <- length(measurement_weeks)

    # Build sigma with guaranteed PD using time-based AR(1)
    sigma_parts <- build_sigma_guaranteed_pd(measurement_weeks, params$biomarker_correlation, params)
    idx <- sigma_parts$indices
    effective_bm_corr <- sigma_parts$effective_c.bm

    # Generate correlated data for each participant using TWO-STAGE approach
    all_participant_data <- list()

    for (pid in 1:n_participants) {
      # Two-stage generation: (biomarker, baseline) then (BR, ER, TR | biomarker, baseline)
      result <- generate_participant_twostage(sigma_parts, idx)

      all_participant_data[[pid]] <- tibble(
        participant_id = pid,
        biomarker = result$biomarker,
        baseline = result$baseline,
        br_random = result$br_random,
        er_random = result$er_random,
        tr_random = result$tr_random
      )
    }

    participant_data <- bind_rows(all_participant_data) %>%
      group_by(participant_id) %>%
      mutate(timepoint_idx = row_number()) %>%
      ungroup()

    # Get biomarker moderation for this condition
    bm_mod <- params$biomarker_moderation

    # Generate trial data
    trial_data <- trial_design %>%
      group_by(participant_id) %>%
      mutate(timepoint_idx = row_number()) %>%
      ungroup() %>%
      left_join(participant_data, by = c("participant_id", "timepoint_idx")) %>%
      group_by(participant_id) %>%
      arrange(week) %>%
      mutate(
        # ===========================================
        # CUMULATIVE TIME TRACKING
        # ===========================================

        # Cumulative weeks on drug
        weeks_on_drug = cumsum(treatment),

        # Cumulative weeks with expectancy (weighted by expectancy level)
        weeks_with_expectancy = cumsum(expectancy),

        # Cumulative weeks in trial (for time trend)
        weeks_in_trial = week - min(week),

        # Track when drug stops for carryover calculation
        time_off = cumsum(treatment == 0),

        # Carryover effect indicator (first week off drug after being on)
        carryover_effect = as.numeric(treatment == 0 & lag(treatment, default = 1) == 1),

        # ===========================================
        # THREE-FACTOR RESPONSE MODEL (RATE-BASED + RANDOM VARIATION)
        # ===========================================

        # Center biomarker for moderation calculation
        bm_centered = (biomarker - biomarker_mean) / biomarker_sd,

        # 1. BR (Biological Response) - accumulates while on drug + random variation
        #    - Mean: BR_rate points per week on drug
        #    - BIOMARKER MODERATION: effect scales with biomarker level
        #    - Carryover: partial at first off-drug timepoint, then 0
        #    - Random: from 26x26 sigma matrix (correlated across time)
        BR_mean = {
          # Treatment effect moderated by biomarker
          # Higher biomarker → stronger BR effect
          effective_BR_rate <- BR_rate * (1 + bm_mod * bm_centered)

          br_accumulated <- lag(weeks_on_drug, default = 0) * effective_BR_rate
          first_off <- treatment == 0 & lag(treatment, default = 1) == 1

          ifelse(treatment == 1,
                 weeks_on_drug * effective_BR_rate,
                 ifelse(first_off,
                        br_accumulated * params$carryover_decay_rate,
                        0))
        },
        BR = BR_mean + br_random,

        # 2. ER (Expectancy Response) - accumulates based on expectancy + random variation
        #    - Mean: ER_rate points per week × expectancy level
        #    - Random: from 26x26 sigma matrix
        ER_mean = weeks_with_expectancy * ER_rate,
        ER = ER_mean + er_random,

        # 3. TR (Time-variant Response) - linear time trend + random variation
        #    - Mean: TR_rate points per week
        #    - Random: from 26x26 sigma matrix
        TR_mean = weeks_in_trial * TR_rate,
        TR = TR_mean + tr_random,

        # ===========================================
        # TOTAL RESPONSE
        # ===========================================
        # Response = baseline + BR + ER + TR
        # Note: baseline already includes between-subject variability from sigma
        # Note: biomarker moderation is applied explicitly to BR_mean above
        response = baseline + BR + ER + TR
      ) %>%
      ungroup()

    # Fit mixed model
    model_result <- tryCatch({
      # Re-center biomarker using sample mean (for model interpretation)
      # Note: data generation uses population parameters for standardization
      trial_data <- trial_data %>%
        mutate(bm_centered = biomarker - mean(biomarker))

      if (params$carryover_decay_rate > 0) {
        model <- lmer(response ~ treatment * bm_centered + week + carryover_effect +
                      (1 | participant_id),
                      data = trial_data)
      } else {
        model <- lmer(response ~ treatment * bm_centered + week +
                      (1 | participant_id),
                      data = trial_data)
      }

      # Extract treatment × biomarker interaction
      coefs <- summary(model)$coefficients
      idx <- which(rownames(coefs) == "treatment:bm_centered")

      if (length(idx) == 0) {
        stop("Interaction term not found")
      }

      # Calculate p-value
      t_val <- coefs[idx, "t value"]
      df_approx <- nrow(trial_data) - nrow(coefs)
      p_value <- 2 * pt(-abs(t_val), df = df_approx)

      tibble(
        iteration = iter,
        design = params$design,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,  # Use effective (grid) value
        carryover_decay_rate = params$carryover_decay_rate,
        effect_size = coefs[idx, "Estimate"],
        se = coefs[idx, "Std. Error"],
        t_value = t_val,
        p_value = p_value,
        significant = p_value < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(iteration = iter, design = params$design,
             biomarker_moderation = params$biomarker_moderation,
             biomarker_correlation = effective_bm_corr,
             carryover_decay_rate = params$carryover_decay_rate, effect_size = NA,
             se = NA, t_value = NA, p_value = NA, significant = NA)
    })

    results <- bind_rows(results, model_result)
  }
}

# =============================================================================
# SUMMARIZE RESULTS
# =============================================================================

cat("\n", strrep("=", 50), "\n")
cat("RESULTS\n")
cat(strrep("=", 50), "\n\n")

summary_results <- results %>%
  group_by(design, biomarker_moderation, biomarker_correlation, carryover_decay_rate) %>%
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(summary_results)

# =============================================================================
# VISUALIZATION
# =============================================================================

library(viridis)

# Adaptive visualization based on data structure
plot_power_results <- function(data) {
  has_moderation <- n_distinct(data$biomarker_moderation) > 1
  has_carryover <- n_distinct(data$carryover_decay_rate) > 1
  has_design <- n_distinct(data$design) > 1

  # Convert to factors
  data <- data %>%
    mutate(
      biomarker_moderation = factor(biomarker_moderation),
      carryover_decay_rate = factor(carryover_decay_rate),
      design = factor(design)
    )

  if (has_moderation) {
    # Primary: biomarker moderation on x-axis
    p <- ggplot(data, aes(x = biomarker_moderation, y = power, fill = design)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100), group = design),
                position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Biomarker Moderation (effect size)", y = "Power")

    # Facet by carryover if multiple values
    if (has_carryover) {
      p <- p + facet_wrap(~ carryover_decay_rate,
                         labeller = labeller(carryover_decay_rate = function(x) paste("Carryover:", x)))
    }

  } else if (has_carryover) {
    # Bar chart by carryover
    p <- ggplot(data, aes(x = carryover_decay_rate, y = power, fill = design)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100), group = design),
                position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Carryover Decay Rate (points/week)", y = "Power")

  } else {
    # Single bar per design
    p <- ggplot(data, aes(x = design, y = power, fill = design)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Design", y = "Power")
  }

  # Common styling
  p <- p +
    scale_fill_viridis_d(name = "Design") +
    labs(title = "Statistical Power by Biomarker Moderation",
         subtitle = sprintf("N=%d, %d iterations per condition",
                           n_participants, n_iterations)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  return(p)
}

# Generate and save plot
p <- plot_power_results(summary_results)
print(p)

# Save outputs
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "power_results.pdf"), p, width = 8, height = 6)
save(results, summary_results, file = file.path(output_dir, "simulation_results.RData"))

cat("\nDone! Results saved to", output_dir, "\n")
cat("- power_results.pdf\n")
cat("- simulation_results.RData\n")
