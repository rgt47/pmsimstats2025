# N-of-1 Trial Simulation - Evenly-Spaced Designs with Power Law Correlation
# OL, Crossover, Parallel designs with 8 evenly-spaced measurement points
# Uses power law correlation structure: 1 / (1 + α * time_gap)^β
# Tests whether faster-decaying correlation enables 8 timepoints for evenly-spaced designs

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(MASS)
  library(zoo)
  library(conflicted)
})

# Resolve conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(nlme::lme)

# ============================================================================
# PARAMETERS
# ============================================================================

n_participants <- 70
n_iterations <- 20

# Three-factor response model - RATE-BASED (points per week)
BR_rate <- 0.5  # Biological Response: drug improvement rate
ER_rate <- 0.2  # Expectancy Response: placebo improvement rate
TR_rate <- 0.1  # Time-variant Response: natural improvement rate
treatment_effect <- BR_rate  # Alias for display purposes

baseline_mean <- 10.0       # Mean baseline response
between_subject_sd <- 2.0   # SD of participant random effects
within_subject_sd <- 1.8    # SD of measurement noise
biomarker_mean <- 5.0       # Mean biomarker value
biomarker_sd <- 2.0         # SD of biomarker

# Correlation parameters (adapted from Hendrickson)
c.br <- 0.75  # BR autocorrelation
c.er <- 0.75  # ER autocorrelation
c.tr <- 0.75  # TR autocorrelation
c.cf1t <- 0.12  # Same-time cross-correlation
c.cfct <- 0.05  # Different-time cross-correlation
c.bm_baseline <- 0.25  # Biomarker-baseline correlation
c.baseline_resp <- 0.3 # Baseline-response correlation

# Power Law correlation parameters
# Faster decay than AR(1) to maintain PD with 8 evenly-spaced points
power_law_alpha <- 0.15  # Inverse power coefficient (higher = faster decay)
power_law_beta <- 1.2    # Power exponent (higher = faster decay)

# Parameter grid - OL, Crossover, Parallel designs with exponential carryover
param_grid <- bind_rows(
  # Open-label design (power conditions)
  expand_grid(
    design = "ol",
    biomarker_moderation = c(0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    t1half = c(0, 0.1, 0.2)  # Hendrickson values (weeks)
  ),
  # Type I error condition for OL
  expand_grid(
    design = "ol",
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    t1half = c(0, 0.1, 0.2)  # Hendrickson values (weeks)
  ),
  # Crossover design
  expand_grid(
    design = "crossover",
    biomarker_moderation = c(0, 0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    t1half = c(0, 0.1, 0.2)  # Hendrickson values (weeks)
  ),
  # Parallel design
  expand_grid(
    design = "parallel",
    biomarker_moderation = c(0, 0.25, 0.35, 0.45, 0.55, 0.65),
    biomarker_correlation = c(0.3),
    t1half = c(0, 0.1, 0.2)  # Hendrickson values (weeks)
  )
)

# ============================================================================
# POWER LAW CORRELATION MATRIX
# ============================================================================

build_power_law_correlation <- function(weeks, alpha, beta) {
  n <- length(weeks)
  Corr <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      time_lag <- abs(weeks[i] - weeks[j])
      Corr[i, j] <- 1 / (1 + alpha * time_lag)^beta
    }
  }
  return(Corr)
}

# ============================================================================
# BUILD SIGMA WITH GUARANTEED PD (Power Law Correlation)
# ============================================================================

build_sigma_guaranteed_pd <- function(weeks, c.bm, params) {
  n_tp <- length(weeks)

  br_idx <- 1:n_tp
  er_idx <- (n_tp + 1):(2 * n_tp)
  tr_idx <- (2 * n_tp + 1):(3 * n_tp)
  bm_idx <- 3 * n_tp + 1
  bl_idx <- 3 * n_tp + 2

  # STAGE 1: Build Sigma_22 (2x2) - Biomarker & Baseline
  Sigma_22 <- matrix(c(
    biomarker_sd^2,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    c.bm_baseline * biomarker_sd * between_subject_sd,
    between_subject_sd^2
  ), 2, 2)

  # STAGE 2: Build Sigma_11 (3*n_tp x 3*n_tp) - Response components
  # Build power law correlation matrices
  Corr_BR <- build_power_law_correlation(weeks, power_law_alpha, power_law_beta)
  Corr_ER <- build_power_law_correlation(weeks, power_law_alpha, power_law_beta)
  Corr_TR <- build_power_law_correlation(weeks, power_law_alpha, power_law_beta)

  # Convert to covariance
  Sigma_BR <- within_subject_sd^2 * Corr_BR
  Sigma_ER <- within_subject_sd^2 * Corr_ER
  Sigma_TR <- within_subject_sd^2 * Corr_TR

  Sigma_11 <- matrix(0, 3 * n_tp, 3 * n_tp)
  Sigma_11[br_idx, br_idx] <- Sigma_BR
  Sigma_11[er_idx, er_idx] <- Sigma_ER
  Sigma_11[tr_idx, tr_idx] <- Sigma_TR

  # Cross-component covariances
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

  # STAGE 3: Build Sigma_12
  effective_c.bm <- c.bm
  effective_c.baseline <- c.baseline_resp

  Sigma_12 <- matrix(0, 3 * n_tp, 2)
  Sigma_12[br_idx, 1] <- effective_c.bm * within_subject_sd * biomarker_sd
  Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[br_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd

  # Verify positive definiteness
  Sigma_22_inv <- solve(Sigma_22)
  cross_term <- Sigma_12 %*% Sigma_22_inv %*% t(Sigma_12)
  Sigma_cond <- Sigma_11 - cross_term
  min_eig <- min(eigen(Sigma_cond, only.values = TRUE)$values)

  if (min_eig <= 1e-6) {
    stop(sprintf(
      paste(
        "\nUNEXPECTED: Parameter combination failed positive definiteness check!",
        "This should NOT happen - parameters were pre-validated.",
        "c.bm = %.2f",
        "Min eigenvalue = %.2e",
        "",
        "Power law parameters: α = %.2f, β = %.2f",
        "Consider adjusting power_law_alpha or power_law_beta.",
        sep = "\n"
      ),
      effective_c.bm, min_eig, power_law_alpha, power_law_beta))
  }

  # Build full sigma
  total_dim <- 3 * n_tp + 2
  Sigma <- matrix(0, total_dim, total_dim)
  Sigma[1:(3*n_tp), 1:(3*n_tp)] <- Sigma_11
  Sigma[(3*n_tp+1):total_dim, (3*n_tp+1):total_dim] <- Sigma_22
  Sigma[1:(3*n_tp), (3*n_tp+1):total_dim] <- Sigma_12
  Sigma[(3*n_tp+1):total_dim, 1:(3*n_tp)] <- t(Sigma_12)

  indices <- list(
    br = br_idx, er = er_idx, tr = tr_idx,
    bm = bm_idx, bl = bl_idx, n_tp = n_tp
  )

  return(list(
    Sigma = Sigma,
    Sigma_22 = Sigma_22,
    Sigma_cond = Sigma_cond,
    Sigma_12 = Sigma_12,
    indices = indices,
    effective_c.bm = effective_c.bm
  ))
}

# ============================================================================
# TWO-STAGE GENERATION
# ============================================================================

generate_participant_twostage <- function(sigma_parts, idx) {
  Sigma_22 <- sigma_parts$Sigma_22
  Sigma_cond <- sigma_parts$Sigma_cond
  Sigma_12 <- sigma_parts$Sigma_12
  n_tp <- idx$n_tp

  # Stage 1: Generate (biomarker, baseline)
  stage1 <- mvrnorm(1, mu = c(biomarker_mean, baseline_mean), Sigma = Sigma_22)
  biomarker <- stage1[1]
  baseline <- stage1[2]

  # Stage 2: Generate responses conditional on stage 1
  z <- c(biomarker - biomarker_mean, baseline - baseline_mean)
  cond_mean <- as.vector(Sigma_12 %*% solve(Sigma_22) %*% z)
  responses <- mvrnorm(1, mu = cond_mean, Sigma = Sigma_cond)

  list(
    biomarker = biomarker,
    baseline = baseline,
    br_random = responses[idx$br],
    er_random = responses[idx$er],
    tr_random = responses[idx$tr]
  )
}

# ============================================================================
# DESIGN STRUCTURES - EVENLY-SPACED WITH 8 POINTS
# ============================================================================

# 8 evenly-spaced points (6-week gaps)
measurement_weeks_evenly <- c(2, 8, 14, 20, 26, 32, 38, 44)

create_ol_design <- function(n_participants, measurement_weeks) {
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      treatment = 1,  # Always on treatment (open-label)
      expectancy = 1.0
    )
}

create_crossover_design <- function(n_participants, measurement_weeks) {
  # Two-period crossover: on for weeks 2-20, off for weeks 26-44
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      treatment = if_else(week <= 20, 1, 0),
      expectancy = 1.0
    )
}

create_parallel_design <- function(n_participants, measurement_weeks) {
  # Parallel: half on treatment all 44 weeks, half on placebo
  path_assignment <- sample(rep(1:2, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],
      treatment = if_else(path == 1, 1, 0),
      expectancy = 1.0
    )
}

# ============================================================================
# PRE-SIMULATION SIGMA VALIDATION (4 unique structures)
# ============================================================================

cat("\n")
cat(paste0("=" , strrep("=", 79) , "\n"))
cat("PRE-SIMULATION SIGMA MATRIX VALIDATION\n")
cat("Power Law Correlation Model (8 Evenly-Spaced Timepoints)\n")
cat(sprintf("Power Law Parameters: α = %.2f, β = %.2f\n", power_law_alpha, power_law_beta))
cat(paste0("=" , strrep("=", 79) , "\n\n"))

# Get unique biomarker correlations in param_grid
unique_bm_corr <- unique(param_grid$biomarker_correlation)

# Define sigma structures to validate
sigma_structures <- list(
  list(design = "ol", weeks = measurement_weeks_evenly, c_bm = 0.3),
  list(design = "crossover", weeks = measurement_weeks_evenly, c_bm = 0.3),
  list(design = "parallel", weeks = measurement_weeks_evenly, c_bm = 0.3),
  list(design = "ol", weeks = measurement_weeks_evenly, c_bm = 0)
)

# Filter to only those that exist in param_grid
sigma_structures <- Filter(function(x) x$c_bm %in% unique_bm_corr, sigma_structures)

# Test each sigma structure
all_valid <- TRUE

for (i in seq_along(sigma_structures)) {
  sigma_def <- sigma_structures[[i]]

  cat(sprintf("Sigma %d: %s design (weeks: %s, c.bm = %.1f)\n",
              i, sigma_def$design, paste(sigma_def$weeks, collapse = ","), sigma_def$c_bm))

  # Build sigma matrix
  sigma_result <- tryCatch({
    build_sigma_guaranteed_pd(
      weeks = sigma_def$weeks,
      c.bm = sigma_def$c_bm,
      params = list()
    )
  }, error = function(e) {
    cat(sprintf("  ✗ FAILED: %s\n", conditionMessage(e)))
    return(NULL)
  })

  if (is.null(sigma_result)) {
    cat("  ✗ Non-positive definite\n\n")
    all_valid <- FALSE
  } else {
    # Check eigenvalues
    eigs <- eigen(sigma_result$Sigma, only.values = TRUE)$values
    min_eig <- min(eigs)
    cond_number <- max(eigs) / max(min_eig, 1e-10)

    cat(sprintf("  ✓ Valid (κ = %.1f)\n", cond_number))

    if (cond_number > 100) {
      cat(sprintf("    ⚠ WARNING: Ill-conditioned (κ > 100)\n"))
    }
    cat("\n")
  }
}

if (!all_valid) {
  stop("\n✗ VALIDATION FAILED: One or more sigma matrices are not positive definite.\n",
       "   Adjust power_law_alpha or power_law_beta.")
}

cat(paste0("=" , strrep("=", 79) , "\n"))
cat("✓ All sigma matrices valid - proceeding with simulation\n")
cat(paste0("=" , strrep("=", 79) , "\n\n"))

# ============================================================================
# RUN SIMULATION
# ============================================================================

# Set up logging
log_file <- "../output/simulation_evenly_spaced_powerlaw_log.txt"
if (!dir.exists("../output")) dir.create("../output", recursive = TRUE)
sink(log_file, split = TRUE)

cat("Running EVENLY-SPACED designs simulation...\n")
cat("Correlation Model: POWER LAW (Alternative for 8 evenly-spaced points)\n")
cat(sprintf("Power Law Parameters: α = %.2f, β = %.2f\n", power_law_alpha, power_law_beta))
cat("Designs: OL, Crossover, Parallel (8 evenly-spaced measurement points)\n")
cat("Treatment effect =", treatment_effect, "\n")
cat("Participants =", n_participants, "\n")
cat("Iterations per condition =", n_iterations, "\n")
cat("Half-life values (weeks):", paste(unique(param_grid$t1half), collapse = ", "), "\n\n")

results <- tibble()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf(
    "Condition %d: design=%s, bm_mod=%.2f, t1half=%.1f\n",
    i, params$design, params$biomarker_moderation, params$t1half
  ))

  for (iter in 1:n_iterations) {
    set.seed(iter * 1000 + i)

    # Create design based on type
    if (params$design == "ol") {
      trial_design <- create_ol_design(n_participants, measurement_weeks_evenly)
    } else if (params$design == "crossover") {
      trial_design <- create_crossover_design(n_participants, measurement_weeks_evenly)
    } else if (params$design == "parallel") {
      trial_design <- create_parallel_design(n_participants, measurement_weeks_evenly)
    }

    n_timepoints <- length(measurement_weeks_evenly)

    sigma_parts <- build_sigma_guaranteed_pd(
      measurement_weeks_evenly,
      params$biomarker_correlation,
      params
    )
    idx <- sigma_parts$indices
    effective_bm_corr <- sigma_parts$effective_c.bm

    all_participant_data <- list()
    for (pid in 1:n_participants) {
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

    bm_mod <- params$biomarker_moderation
    t1half <- params$t1half

    trial_data <- trial_design %>%
      group_by(participant_id) %>%
      mutate(timepoint_idx = row_number()) %>%
      ungroup() %>%
      left_join(participant_data, by = c("participant_id", "timepoint_idx")) %>%
      group_by(participant_id) %>%
      arrange(week) %>%
      mutate(
        weeks_on_drug = cumsum(treatment),
        weeks_with_expectancy = cumsum(expectancy),
        weeks_in_trial = week - min(week),
        # Identify the discontinuation timepoint
        disc_week = if_else(
          treatment == 0 & lag(treatment, default = 1) == 1,
          week,
          NA_real_
        )
      ) %>%
      fill(disc_week, .direction = "down") %>%
      mutate(
        weeks_off_treatment = if_else(
          treatment == 0 & !is.na(disc_week),
          week - disc_week,
          0
        ),
        bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
        BR_mean = {
          effective_BR_rate <- BR_rate * (1 + bm_mod * bm_centered)
          br_accumulated <- lag(weeks_on_drug, default = 0) * effective_BR_rate

          if (t1half == 0) {
            # No carryover: immediate washout
            ifelse(treatment == 1,
              weeks_on_drug * effective_BR_rate,
              0
            )
          } else {
            # Exponential decay: accum * 0.5^(weeks_off / t1half)
            ifelse(treatment == 1,
              weeks_on_drug * effective_BR_rate,
              br_accumulated * (0.5^(weeks_off_treatment / t1half))
            )
          }
        },
        BR = BR_mean + br_random,
        ER_mean = weeks_with_expectancy * ER_rate,
        ER = ER_mean + er_random,
        TR_mean = weeks_in_trial * TR_rate,
        TR = TR_mean + tr_random,
        response = baseline + BR + ER + TR
      ) %>%
      ungroup()

    model_result <- tryCatch({
      trial_data <- trial_data %>%
        mutate(bm_centered = biomarker - mean(biomarker))

      # Include carryover effect only if t1half > 0
      if (params$t1half > 0) {
        trial_data <- trial_data %>%
          mutate(
            carryover_effect = as.numeric(
              treatment == 0 & lag(treatment, default = 1) == 1
            )
          )

        model <- lme(
          response ~ treatment * bm_centered + week + carryover_effect,
          random = ~ 1 | participant_id,
          correlation = corCAR1(form = ~ week | participant_id),
          data = trial_data,
          control = lmeControl(opt = "optim")
        )
      } else {
        model <- lme(
          response ~ treatment * bm_centered + week,
          random = ~ 1 | participant_id,
          correlation = corCAR1(form = ~ week | participant_id),
          data = trial_data,
          control = lmeControl(opt = "optim")
        )
      }

      coefs <- summary(model)$tTable
      interaction_term <- "treatment:bm_centered"
      idx_coef <- which(rownames(coefs) == interaction_term)

      tibble(
        iteration = iter,
        design = params$design,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,
        t1half = params$t1half,
        effect_size = coefs[idx_coef, "Value"],
        se = coefs[idx_coef, "Std.Error"],
        t_value = coefs[idx_coef, "t-value"],
        p_value = coefs[idx_coef, "p-value"],
        significant = coefs[idx_coef, "p-value"] < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(
        iteration = iter,
        design = params$design,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,
        t1half = params$t1half,
        effect_size = NA, se = NA, t_value = NA, p_value = NA, significant = NA
      )
    })

    results <- bind_rows(results, model_result)
  }
}

# ============================================================================
# SUMMARIZE AND SAVE
# ============================================================================

cat("\n", strrep("=", 50), "\n")
cat("RESULTS - EVENLY-SPACED DESIGNS (POWER LAW)\n")
cat(strrep("=", 50), "\n\n")

summary_results <- results %>%
  group_by(design, biomarker_moderation, biomarker_correlation, t1half) %>%
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

print(summary_results)

# ============================================================================
# VISUALIZATION
# ============================================================================

library(viridis)

plot_power_heatmap <- function(data) {
  design_labels <- c(
    "ol" = "Open-Label\n(OL)",
    "crossover" = "Crossover\n(CO)",
    "parallel" = "Parallel\n(PA)"
  )

  data <- data %>%
    mutate(
      t1half = factor(t1half),
      biomarker_moderation = factor(
        biomarker_moderation,
        levels = rev(sort(unique(biomarker_moderation)))
      ),
      design = factor(
        design_labels[as.character(design)],
        levels = design_labels
      )
    )

  p <- ggplot(
    data,
    aes(x = t1half, y = biomarker_moderation, fill = power)
  ) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.0f%%", power * 100)),
      color = "black",
      size = 4,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      name = "Power",
      low = "#d73027",
      mid = "#fee08b",
      high = "#1a9850",
      midpoint = 0.5,
      limits = c(0, 1),
      labels = scales::percent
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    facet_wrap(~design, ncol = 3) +
    labs(
      title = "Statistical Power - Evenly-Spaced Designs with Power Law Correlation",
      subtitle = sprintf(
        "N=%d participants, %d iterations per condition\nPower Law: α=%.2f, β=%.2f",
        n_participants, n_iterations, power_law_alpha, power_law_beta
      ),
      x = "Half-Life Parameter (weeks)",
      y = "Biomarker Moderation\n(effect size)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey90", color = NA)
    )

  return(p)
}

p_heatmap <- plot_power_heatmap(summary_results)
print(p_heatmap)

# Save outputs
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
ggsave(
  file.path(output_dir, sprintf("power_heatmap_evenly_spaced_powerlaw_%s.pdf", timestamp)),
  p_heatmap,
  width = 12,
  height = 7
)

save(
  results, summary_results,
  file = file.path(output_dir, "simulation_evenly_spaced_powerlaw_results.RData")
)

cat("\nDone! Results saved to", output_dir, "\n")
cat(sprintf("- power_heatmap_evenly_spaced_powerlaw_%s.pdf\n", timestamp))
cat("- simulation_evenly_spaced_powerlaw_results.RData\n")

# Close log file
sink()
