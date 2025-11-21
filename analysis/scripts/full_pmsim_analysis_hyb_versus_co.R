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
n_iterations <- 5

# Three-factor response model - RATE-BASED (points per week)
BR_rate <- 0.5                 # Biological Response: drug improvement rate
ER_rate <- 0.2                 # Expectancy Response: placebo improvement rate
TR_rate <- 0.1                 # Time-variant Response: natural improvement rate

baseline_mean <- 10.0          # Mean baseline response
between_subject_sd <- 3.0      # SD of participant random effects
within_subject_sd <- 2.8       # SD of measurement noise
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

# Parameter grid - what we're testing
param_grid <- expand_grid(
  biomarker_correlation = c(0.3),
  carryover_decay_rate = c(0, 0.25)    # Linear decay: points lost per week off drug
)

# =============================================================================
# BUILD 26x26 HENDRICKSON SIGMA MATRIX
# =============================================================================

# Structure: [BR(8) | ER(8) | TR(8) | Biomarker(1) | Baseline(1)]
build_hendrickson_sigma <- function(n_timepoints, c.bm, params) {
  n <- 3 * n_timepoints + 2  # 26 for 8 timepoints
  R <- matrix(0, n, n)

  # Index ranges
  br_idx <- 1:n_timepoints
  er_idx <- (n_timepoints + 1):(2 * n_timepoints)
  tr_idx <- (2 * n_timepoints + 1):(3 * n_timepoints)
  bm_idx <- 3 * n_timepoints + 1
  bl_idx <- 3 * n_timepoints + 2

  # Helper: fill autocorrelation block
  fill_auto <- function(idx, rho) {
    for (i in seq_along(idx)) {
      for (j in seq_along(idx)) {
        if (i == j) {
          R[idx[i], idx[j]] <<- 1
        } else {
          R[idx[i], idx[j]] <<- rho
        }
      }
    }
  }

  # Helper: fill cross-correlation block
  fill_cross <- function(idx1, idx2, rho_same, rho_diff) {
    for (i in seq_along(idx1)) {
      for (j in seq_along(idx2)) {
        if (i == j) {
          R[idx1[i], idx2[j]] <<- rho_same
          R[idx2[j], idx1[i]] <<- rho_same
        } else {
          R[idx1[i], idx2[j]] <<- rho_diff
          R[idx2[j], idx1[i]] <<- rho_diff
        }
      }
    }
  }

  # 1. Autocorrelations within each response type
  fill_auto(br_idx, c.br)
  fill_auto(er_idx, c.er)
  fill_auto(tr_idx, c.tr)

  # 2. Cross-correlations between response types
  fill_cross(br_idx, er_idx, c.cf1t, c.cfct)
  fill_cross(br_idx, tr_idx, c.cf1t, c.cfct)
  fill_cross(er_idx, tr_idx, c.cf1t, c.cfct)

  # 3. Biomarker correlations
  R[bm_idx, bm_idx] <- 1
  # Biomarker-BR correlation (the key parameter we're testing)
  for (i in br_idx) {
    R[bm_idx, i] <- c.bm
    R[i, bm_idx] <- c.bm
  }
  # Biomarker-ER and Biomarker-TR (weaker)
  for (i in c(er_idx, tr_idx)) {
    R[bm_idx, i] <- c.bm * 0.5
    R[i, bm_idx] <- c.bm * 0.5
  }

  # 4. Baseline correlations
  R[bl_idx, bl_idx] <- 1
  R[bm_idx, bl_idx] <- c.bm_baseline
  R[bl_idx, bm_idx] <- c.bm_baseline
  # Baseline-response correlations
  for (i in c(br_idx, er_idx, tr_idx)) {
    R[bl_idx, i] <- c.baseline_resp
    R[i, bl_idx] <- c.baseline_resp
  }

  # 5. Convert correlation to covariance matrix
  # SD vector: within_subject_sd for responses, biomarker_sd, between_subject_sd for baseline
  sd_vec <- c(
    rep(within_subject_sd, 3 * n_timepoints),  # BR, ER, TR
    biomarker_sd,                               # Biomarker
    between_subject_sd                          # Baseline
  )

  # Sigma = diag(SD) %*% R %*% diag(SD)
  Sigma <- diag(sd_vec) %*% R %*% diag(sd_vec)

  # Check positive definiteness
  eigenvalues <- eigen(Sigma, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    warning("Sigma matrix is not positive definite!")
    return(NULL)
  }

  return(list(
    Sigma = Sigma,
    R = R,
    indices = list(br = br_idx, er = er_idx, tr = tr_idx, bm = bm_idx, bl = bl_idx)
  ))
}

# =============================================================================
# HYBRID DESIGN STRUCTURE
# =============================================================================

# Measurement schedule
measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)

# Function to create design with randomized path assignment
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
  cat(sprintf("Condition %d: biomarker=%.2f, carryover_decay=%.2f\n",
              i, params$biomarker_correlation, params$carryover_decay_rate))

  # Run iterations
  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    # Create design with fresh path randomization for this iteration
    hybrid_design <- create_hybrid_design(n_participants, measurement_weeks)

    n_timepoints <- length(measurement_weeks)

    # Build 26x26 Hendrickson sigma matrix
    sigma_result <- build_hendrickson_sigma(n_timepoints, params$biomarker_correlation, params)
    if (is.null(sigma_result)) {
      cat("  Skipping iteration", iter, "- non-positive definite sigma\n")
      next
    }
    Sigma <- sigma_result$Sigma
    idx <- sigma_result$indices

    # Mean vector (all zeros - we add means separately)
    mu <- rep(0, nrow(Sigma))

    # Generate correlated data for each participant from 26-dim MVN
    all_participant_data <- list()

    for (pid in 1:n_participants) {
      # Generate from multivariate normal
      mvn_draw <- mvrnorm(1, mu = mu, Sigma = Sigma)

      # Extract components
      br_values <- mvn_draw[idx$br]         # BR residuals at each timepoint
      er_values <- mvn_draw[idx$er]         # ER residuals at each timepoint
      tr_values <- mvn_draw[idx$tr]         # TR residuals at each timepoint
      biomarker <- mvn_draw[idx$bm] + biomarker_mean   # Biomarker (add mean)
      baseline <- mvn_draw[idx$bl] + baseline_mean     # Baseline (add mean)

      all_participant_data[[pid]] <- tibble(
        participant_id = pid,
        biomarker = biomarker,
        baseline = baseline,
        br_random = br_values,
        er_random = er_values,
        tr_random = tr_values
      )
    }

    participant_data <- bind_rows(all_participant_data) %>%
      group_by(participant_id) %>%
      mutate(timepoint_idx = row_number()) %>%
      ungroup()

    # Generate trial data
    trial_data <- hybrid_design %>%
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

        # ===========================================
        # THREE-FACTOR RESPONSE MODEL (RATE-BASED + RANDOM VARIATION)
        # ===========================================

        # 1. BR (Biological Response) - accumulates while on drug + random variation
        #    - Mean: BR_rate points per week on drug
        #    - Carryover: partial at first off-drug timepoint, then 0
        #    - Random: from 26x26 sigma matrix (correlated across time)
        BR_mean = {
          br_accumulated <- lag(weeks_on_drug, default = 0) * BR_rate
          first_off <- treatment == 0 & lag(treatment, default = 1) == 1

          ifelse(treatment == 1,
                 weeks_on_drug * BR_rate,
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
        # Note: biomarker moderation is captured in sigma correlations
        response = baseline + BR + ER + TR
      ) %>%
      ungroup()

    # Fit mixed model
    model_result <- tryCatch({
      # Center biomarker for better interpretation
      trial_data <- trial_data %>%
        mutate(bm_centered = biomarker - mean(biomarker))

      if (params$carryover_t1half > 0) {
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
        biomarker_correlation = params$biomarker_correlation,
        carryover_decay_rate = params$carryover_decay_rate,
        effect_size = coefs[idx, "Estimate"],
        se = coefs[idx, "Std. Error"],
        t_value = t_val,
        p_value = p_value,
        significant = p_value < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(iteration = iter, biomarker_correlation = params$biomarker_correlation,
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
  group_by(biomarker_correlation, carryover_decay_rate) %>%
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
  has_carryover <- n_distinct(data$carryover_decay_rate) > 1
  has_biomarker <- n_distinct(data$biomarker_correlation) > 1

  # Convert to factors
  data <- data %>%
    mutate(
      carryover_decay_rate = factor(carryover_decay_rate),
      biomarker_correlation = factor(biomarker_correlation)
    )

  if (has_carryover && has_biomarker) {
    # 2D heatmap
    p <- ggplot(data, aes(x = carryover_decay_rate, y = biomarker_correlation, fill = power)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), color = "black", size = 5) +
      labs(x = "Carryover Decay Rate (points/week)", y = "Biomarker Correlation")

  } else if (has_carryover) {
    # Bar chart by carryover
    p <- ggplot(data, aes(x = carryover_decay_rate, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Carryover Decay Rate (points/week)", y = "Power")

  } else if (has_biomarker) {
    # Bar chart by biomarker
    p <- ggplot(data, aes(x = biomarker_correlation, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Biomarker Correlation", y = "Power")

  } else {
    # Single bar
    p <- ggplot(data, aes(x = 1, y = power, fill = power)) +
      geom_col(width = 0.4) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 6) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "", y = "Power") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  # Common styling
  p <- p +
    scale_fill_viridis_c(name = "Power", labels = scales::percent, limits = c(0, 1)) +
    labs(title = "Statistical Power: Constant Effect Model",
         subtitle = sprintf("N=%d, treatment effect=%.1f, %d iterations",
                           n_participants, treatment_effect, n_iterations)) +
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
