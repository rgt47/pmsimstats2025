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
autocorr <- 0.8                # Autocorrelation between adjacent timepoints
autocorr_decay <- 0.9          # How fast autocorrelation decays with lag

# Parameter grid - what we're testing
param_grid <- expand_grid(
  biomarker_correlation = c(0.3),
  carryover_t1half = c(0, 0.5)
)

# =============================================================================
# BUILD CORRELATION MATRIX
# =============================================================================

# Function to build AR(1)-like correlation matrix for timepoints
build_corr_matrix <- function(n_timepoints, rho, decay = 0.9) {
  # Create correlation matrix where corr decreases with lag
  # corr(t1, t2) = rho * decay^|t1 - t2|
  mat <- matrix(1, n_timepoints, n_timepoints)
  for (i in 1:n_timepoints) {
    for (j in 1:n_timepoints) {
      if (i != j) {
        lag <- abs(i - j)
        mat[i, j] <- rho * (decay ^ (lag - 1))
      }
    }
  }
  return(mat)
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
  cat(sprintf("Condition %d: biomarker=%.1f, carryover=%.1f\n",
              i, params$biomarker_correlation, params$carryover_t1half))

  # Run iterations
  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    # Create design with fresh path randomization for this iteration
    hybrid_design <- create_hybrid_design(n_participants, measurement_weeks)

    n_timepoints <- length(measurement_weeks)

    # Build correlation matrix for residuals
    R <- build_corr_matrix(n_timepoints, autocorr, autocorr_decay)

    # Convert to covariance matrix
    Sigma <- R * within_subject_sd^2

    # Generate correlated data for each participant
    all_participant_data <- list()

    for (pid in 1:n_participants) {
      # Generate biomarker (correlated with response via shared factor)
      biomarker <- rnorm(1, biomarker_mean, biomarker_sd)

      # Random intercept (between-subject variability)
      random_intercept <- rnorm(1, 0, between_subject_sd)

      # Generate correlated residuals across timepoints
      correlated_residuals <- mvrnorm(1, mu = rep(0, n_timepoints), Sigma = Sigma)

      # Biomarker-correlated component of residuals
      # Higher biomarker -> residuals shifted by biomarker_correlation
      bm_effect_on_residuals <- (biomarker - biomarker_mean) / biomarker_sd *
                                 params$biomarker_correlation * within_subject_sd

      all_participant_data[[pid]] <- tibble(
        participant_id = pid,
        biomarker = biomarker,
        random_intercept = random_intercept,
        residual = correlated_residuals + bm_effect_on_residuals
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
        drug_stopped = treatment == 0 & lag(treatment, default = 1) == 1,
        time_off = cumsum(treatment == 0),

        # BR level at last on-drug timepoint (for carryover decay)
        br_at_stop = weeks_on_drug * BR_rate,

        # ===========================================
        # THREE-FACTOR RESPONSE MODEL (RATE-BASED)
        # ===========================================

        # 1. BR (Biological Response) - accumulates while on drug
        #    - Rate: BR_rate points per week on drug
        #    - When off: decays from accumulated level
        BR = if (params$carryover_t1half > 0) {
          ifelse(treatment == 1,
                 weeks_on_drug * BR_rate,                    # Accumulating while on
                 lag(weeks_on_drug, default = 0) * BR_rate * # Decay from last level
                   (0.5)^(time_off / params$carryover_t1half))
        } else {
          ifelse(treatment == 1,
                 weeks_on_drug * BR_rate,                    # Accumulating while on
                 0)                                          # Instant drop when off
        },

        # 2. ER (Expectancy Response) - accumulates based on expectancy
        #    - Rate: ER_rate points per week × expectancy level
        #    - Open-label: full rate (1.0 × ER_rate)
        #    - Blinded: half rate (0.5 × ER_rate)
        ER = weeks_with_expectancy * ER_rate,

        # 3. TR (Time-variant Response) - linear time trend
        #    - Rate: TR_rate points per week
        #    - Accumulates regardless of treatment
        TR = weeks_in_trial * TR_rate,

        # Biomarker moderation: biomarker amplifies BR (not ER or TR)
        # Higher biomarker = stronger biological response to treatment
        biomarker_moderation = (biomarker - biomarker_mean) / biomarker_sd *
                               params$biomarker_correlation * BR,

        # ===========================================
        # TOTAL RESPONSE
        # ===========================================
        # Response = baseline + random effect + BR + ER + TR + biomarker moderation + noise
        response = baseline_mean +
                   random_intercept +
                   BR +                      # Biological (drug) effect
                   ER +                      # Expectancy (placebo) effect
                   TR +                      # Time-variant effect
                   biomarker_moderation +    # Biomarker amplifies BR
                   residual
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
        carryover_t1half = params$carryover_t1half,
        effect_size = coefs[idx, "Estimate"],
        se = coefs[idx, "Std. Error"],
        t_value = t_val,
        p_value = p_value,
        significant = p_value < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(iteration = iter, biomarker_correlation = params$biomarker_correlation,
             carryover_t1half = params$carryover_t1half, effect_size = NA,
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
  group_by(biomarker_correlation, carryover_t1half) %>%
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
  has_carryover <- n_distinct(data$carryover_t1half) > 1
  has_biomarker <- n_distinct(data$biomarker_correlation) > 1

  # Convert to factors
  data <- data %>%
    mutate(
      carryover_t1half = factor(carryover_t1half),
      biomarker_correlation = factor(biomarker_correlation)
    )

  if (has_carryover && has_biomarker) {
    # 2D heatmap
    p <- ggplot(data, aes(x = carryover_t1half, y = biomarker_correlation, fill = power)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), color = "black", size = 5) +
      labs(x = "Carryover Half-life (weeks)", y = "Biomarker Correlation")

  } else if (has_carryover) {
    # Bar chart by carryover
    p <- ggplot(data, aes(x = carryover_t1half, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Carryover Half-life (weeks)", y = "Power")

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
