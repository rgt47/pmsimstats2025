# Simplified N-of-1 Trial Simulation - Hybrid Design Only
# Stripped down for learning

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(MASS)
  library(corpcor)
  library(conflicted)
})

# Resolve conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(lmerTest::lmer)

source("pm_functions.R")

# =============================================================================
# PARAMETERS
# =============================================================================

n_participants <- 70
n_iterations <- 5
treatment_effect <- 5.0
carryover_scale <- 1
between_subject_sd <- 3.0
within_subject_sd <- 2.8

# Parameter grid - what we're testing
param_grid <- expand_grid(
  biomarker_correlation = c(0.3),
  carryover_t1half = c(0, 0.5)
)

# Model parameters for sigma matrix
model_params <- list(
  N = n_participants,
  c.bm = 0.3,
  carryover_t1half = 0,
  c.tv = 0.6, c.pb = 0.6, c.br = 0.6,
  c.cf1t = 0.2, c.cfct = 0.1
)

# Response parameters
resp_param <- tibble(
  cat = c("time_variant", "pharm_biomarker", "bio_response"),
  max = c(1.0, 1.0, treatment_effect),
  disp = c(2.0, 2.0, 2.0),
  rate = c(0.3, 0.3, 0.3),
  sd = c(within_subject_sd, within_subject_sd, within_subject_sd)
)

# Baseline parameters
baseline_param <- tibble(
  cat = c("biomarker", "baseline"),
  m = c(5.0, 10.0),
  sd = c(2.0, between_subject_sd)
)

# =============================================================================
# CREATE HYBRID DESIGN
# =============================================================================

measurement_weeks <- c(4, 8, 9, 10, 11, 12, 16, 20)
path_assignment <- sample(rep(1:4, length.out = n_participants))

hybrid_design <- expand_grid(
  participant_id = 1:n_participants,
  week = measurement_weeks
) %>%
  mutate(
    path = path_assignment[participant_id],
    treatment = case_when(
      week %in% c(4, 8) ~ 1,
      week == 9 ~ 1,
      week %in% c(10, 11, 12) & path %in% c(1, 2) ~ 1,
      week %in% c(10, 11, 12) & path %in% c(3, 4) ~ 0,
      week == 16 & path %in% c(1, 3) ~ 1,
      week == 16 & path %in% c(2, 4) ~ 0,
      week == 20 & path %in% c(1, 3) ~ 0,
      week == 20 & path %in% c(2, 4) ~ 1,
      TRUE ~ NA_real_
    ),
    expectancy = if_else(week %in% c(4, 8), 1, 0.5)
  )

# Create path templates for simulation
prepare_design <- function(design_df) {
  design_df %>%
    arrange(participant_id, week) %>%
    mutate(
      timepoint_name = paste0("W", week),
      tod = treatment,
      e = expectancy,
      t_wk = 1
    ) %>%
    group_by(participant_id) %>%
    mutate(
      tpb = cumsum(e),
      tsd = if_else(tod == 0, cumsum(tod == 0), 0)
    ) %>%
    ungroup() %>%
    select(timepoint_name, t_wk, e, tod, tsd, tpb)
}

# Get one template per path
design_paths <- map(1:4, function(path_id) {
  hybrid_design %>%
    filter(path == path_id) %>%
    filter(participant_id == min(participant_id)) %>%
    prepare_design()
})

# =============================================================================
# RUN SIMULATION
# =============================================================================

cat("Running simulation...\n")
results <- tibble()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf("\nCondition %d: biomarker=%.1f, carryover=%.1f\n",
              i, params$biomarker_correlation, params$carryover_t1half))

  # Update model params
  model_params$c.bm <- params$biomarker_correlation
  model_params$carryover_t1half <- params$carryover_t1half

  # Build sigma matrix
  sigma <- build_sigma_matrix(
    model_params, resp_param, baseline_param, design_paths[[1]],
    c("time_variant", "pharm_biomarker", "bio_response"),
    c("tv", "pb", "br"), verbose = FALSE
  )

  if (is.null(sigma)) {
    cat("  Skipping - non-positive definite sigma\n")
    next
  }

  # Run iterations
  for (iter in 1:n_iterations) {
    # Generate data for each path
    all_data <- map_dfr(1:4, function(path_id) {
      path_participants <- hybrid_design %>%
        filter(path == path_id) %>%
        pull(participant_id) %>%
        unique()

      n_in_path <- length(path_participants)
      if (n_in_path == 0) return(tibble())

      path_model_params <- model_params
      path_model_params$N <- n_in_path

      sim_data <- generate_data(
        model_param = path_model_params,
        resp_param = resp_param,
        baseline_param = baseline_param,
        trial_design = design_paths[[path_id]],
        empirical = FALSE,
        make_positive_definite = FALSE,
        seed = iter * 1000 + path_id,
        scale_factor = carryover_scale,
        cached_sigma = sigma
      )

      sim_data %>%
        mutate(path = path_id, participant_id = path_participants[participant_id])
    })

    # Reshape for analysis
    timepoint_cols <- grep("^W\\d+$", names(all_data), value = TRUE)

    analysis_data <- all_data %>%
      select(participant_id, biomarker, path, all_of(timepoint_cols)) %>%
      rename(bm = biomarker) %>%
      pivot_longer(cols = all_of(timepoint_cols), names_to = "timepoint", values_to = "response") %>%
      mutate(week = as.integer(str_replace(timepoint, "W", ""))) %>%
      left_join(hybrid_design %>% select(participant_id, week, treatment),
                by = c("participant_id", "week")) %>%
      group_by(participant_id) %>%
      arrange(week) %>%
      mutate(
        carryover_effect = if (params$carryover_t1half > 0) {
          tsd <- cumsum(treatment == 0 & lag(treatment, default = 0) == 1)
          ifelse(treatment == 0 & tsd > 0, (0.5)^(tsd / params$carryover_t1half), 0)
        } else 0
      ) %>%
      ungroup()

    # Fit model
    model_result <- tryCatch({
      if (params$carryover_t1half > 0) {
        model <- lmer(response ~ treatment * bm + week + carryover_effect + (1 | participant_id),
                      data = analysis_data)
      } else {
        model <- lmer(response ~ treatment * bm + week + (1 | participant_id),
                      data = analysis_data)
      }

      coefs <- summary(model)$coefficients
      idx <- which(rownames(coefs) == "treatment:bm")
      p_value <- 2 * pt(-abs(coefs[idx, "t value"]), df = nrow(analysis_data) - nrow(coefs))

      tibble(
        iteration = iter,
        biomarker_correlation = params$biomarker_correlation,
        carryover_t1half = params$carryover_t1half,
        effect_size = coefs[idx, "Estimate"],
        p_value = p_value,
        significant = p_value < 0.05
      )
    }, error = function(e) {
      cat("  Error in iteration", iter, ":", conditionMessage(e), "\n")
      tibble(iteration = iter, biomarker_correlation = params$biomarker_correlation,
             carryover_t1half = params$carryover_t1half, effect_size = NA,
             p_value = NA, significant = NA)
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
    n = n(),
    .groups = "drop"
  )

print(summary_results)

# =============================================================================
# VISUALIZATION
# =============================================================================

library(viridis)

# General heatmap function - adapts to available variables
plot_power_heatmap <- function(data) {
  # Detect which variables have multiple values
  has_design <- "design" %in% names(data) && n_distinct(data$design) > 1
  has_carryover <- "carryover_t1half" %in% names(data) && n_distinct(data$carryover_t1half) > 1
  has_biomarker <- "biomarker_correlation" %in% names(data) && n_distinct(data$biomarker_correlation) > 1

  # Ensure factors for proper ordering
  if ("carryover_t1half" %in% names(data)) {
    data <- data %>% mutate(carryover_t1half = factor(carryover_t1half))
  }
  if ("biomarker_correlation" %in% names(data)) {
    data <- data %>% mutate(biomarker_correlation = factor(biomarker_correlation))
  }
  if ("design" %in% names(data)) {
    data <- data %>% mutate(design = factor(design))
  }

  # Determine x and y axes based on what's varying
  if (has_carryover && has_biomarker) {
    # 2D heatmap: carryover x biomarker
    p <- ggplot(data, aes(x = carryover_t1half, y = biomarker_correlation, fill = power)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), color = "black", size = 4) +
      labs(x = "Carryover Half-life (weeks)", y = "Biomarker Correlation")

    if (has_design) {
      p <- p + facet_wrap(~ design)
    }

  } else if (has_carryover && has_design) {
    # 2D heatmap: carryover x design
    p <- ggplot(data, aes(x = carryover_t1half, y = design, fill = power)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), color = "black", size = 4) +
      labs(x = "Carryover Half-life (weeks)", y = "Design")

    if (has_biomarker) {
      p <- p + facet_wrap(~ biomarker_correlation, labeller = label_both)
    }

  } else if (has_biomarker && has_design) {
    # 2D heatmap: biomarker x design
    p <- ggplot(data, aes(x = biomarker_correlation, y = design, fill = power)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), color = "black", size = 4) +
      labs(x = "Biomarker Correlation", y = "Design")

    if (has_carryover) {
      p <- p + facet_wrap(~ carryover_t1half, labeller = label_both)
    }

  } else if (has_carryover) {
    # 1D bar chart: carryover only
    p <- ggplot(data, aes(x = carryover_t1half, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Carryover Half-life (weeks)", y = "Power")

  } else if (has_biomarker) {
    # 1D bar chart: biomarker only
    p <- ggplot(data, aes(x = biomarker_correlation, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Biomarker Correlation", y = "Power")

  } else if (has_design) {
    # 1D bar chart: design only
    p <- ggplot(data, aes(x = design, y = power, fill = power)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%.0f%%", power * 100)), vjust = -0.5, size = 5) +
      scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
      labs(x = "Design", y = "Power")

  } else {
    # Single value - just show it
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
    labs(title = "Statistical Power Analysis") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  return(p)
}

# Generate plot
p <- plot_power_heatmap(summary_results)
print(p)

# Save plot
ggsave("../output/power_heatmap.pdf", p, width = 8, height = 6)

# Save results
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
save(results, summary_results, file = file.path(output_dir, "full_pmsim_analysis_hyb_versus_co.RData"))

cat("\nDone! Results and plot saved.\n")
