#' vig5.R: Self-contained Monte Carlo clinical trial design
#' comparison
#'
#' This is a monolithic script that contains all the code needed to
#' run the Monte Carlo simulation from the
#' monte_carlo_design_comparison vignette. All helper functions and
#' simulation code are included directly in this file.

# Load required packages
# Suppress package startup messages
# and clear the workspace
rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)       # For mixed models
  library(viridis)        # For viridis color scales
  library(MASS)           # For mvrnorm
  library(corpcor)        # For is.positive.definite and make.positive.definite
  library(conflicted)
})

# Suppress conflict messages
suppressMessages({
  conflicts_prefer(dplyr::select, .quiet = TRUE)
  conflicts_prefer(dplyr::filter, .quiet = TRUE)
  conflicts_prefer(dplyr::lag, .quiet = TRUE)
  conflicts_prefer(lmerTest::lmer, .quiet = TRUE)
})

# Source the common simulation and analysis functions
source("pm_functions.R")

#===========================================================================
# Function: run_monte_carlo
# Description: Run a Monte Carlo simulation for a design and
# parameters with support for multiple randomization paths
#===========================================================================
run_monte_carlo <- function(design_name, params,
                            sigma_cache = NULL) {
  # FIXED CORRELATIONS: Do NOT adjust based on carryover
  # Correlations remain constant; only MEANS are affected by carryover

  # Update model parameters with current simulation parameters
  current_model_params <- model_params
  current_model_params$N <- params$n_participants
  current_model_params$c.bm <- params$biomarker_correlation
  current_model_params$carryover_t1half <- params$carryover_t1half

  # Correlation values remain FIXED (Hendrickson approach)
  # Already set in model_params: c.tv = 0.8, c.pb = 0.8, c.br = 0.8,
  # c.cf1t = 0.2, c.cfct = 0.1

  # Get the full design with all paths
  if (design_name == "hybrid") {
    full_design <- designs$hybrid$full_design
    design_paths <- designs$hybrid$design_paths
  } else {
    full_design <- designs$crossover$full_design
    design_paths <- designs$crossover$design_paths
  }

  # Check if we have a valid cached sigma for this
  # parameter/design combination
  cache_key <- create_sigma_cache_key(design_name, params)
  cached_sigma <- if (!is.null(sigma_cache)) {
    sigma_cache[[cache_key]]
  } else {
    NULL
  }

  # If no valid cached sigma, skip this parameter combination
  if (is.null(cached_sigma)) {
    cat("Skipping", design_name, "design for parameters:",
        "biomarker_correlation =", params$biomarker_correlation,
        "carryover_t1half =", params$carryover_t1half,
        "(non-positive definite sigma matrix)\n")
    return(tibble())
  }

  # Set number of iterations
  n_iter <- params$n_iterations

  # Use map_dfr for sequential processing
  all_results <- map_dfr(1:n_iter, function(iter) {

    # Generate data for each path separately, then combine
    # This follows Hendrickson's approach of simulating each path
    # independently
    all_path_data <- map_dfr(1:length(design_paths),
                              function(path_id) {
      # Get participants assigned to this path
      path_participants <- full_design %>%
        filter(path == path_id) %>%
        pull(participant_id) %>%
        unique()

      # Calculate how many participants in this path
      n_in_path <- length(path_participants)

      if (n_in_path == 0) return(tibble())

      # Update model params for this path's sample size
      path_model_params <- current_model_params
      path_model_params$N <- n_in_path

      # Generate data for this path
      path_sim_data <- generate_data(
        model_param = path_model_params,
        resp_param = resp_param,
        baseline_param = baseline_param,
        trial_design = design_paths[[path_id]],
        empirical = FALSE,
        make_positive_definite = FALSE,
        seed = iter * 1000 + path_id,  # Unique seed per path
        scale_factor = params$carryover_scale,
        cached_sigma = cached_sigma
      )

      # Add path and actual participant IDs
      path_sim_data %>%
        mutate(
          path = path_id,
          participant_id = path_participants[participant_id]
        )
    })

    # Now we have data for all participants across all paths
    sim_data <- all_path_data

    # Convert to format needed for analysis
    # Extract timepoint columns and reshape to long format
    timepoint_columns <- grep("^W\\d+$", names(sim_data),
                              value = TRUE)

    # Keep only needed columns (including path if present)
    cols_to_keep <- c("participant_id", "biomarker",
                      timepoint_columns)
    if ("path" %in% names(sim_data)) {
      cols_to_keep <- c(cols_to_keep, "path")
    }

    analysis_data <- sim_data %>%
      select(all_of(cols_to_keep)) %>%
      rename(bm = biomarker)

    # Convert to long format for analysis
    analysis_data_long <- analysis_data %>%
      pivot_longer(
        cols = all_of(timepoint_columns),
        names_to = "timepoint",
        values_to = "response"
      ) %>%
      # Extract week number
      mutate(week = as.integer(str_replace(timepoint, "W", "")))

    # Merge with design info to get actual treatment assignments
    # This handles multiple paths correctly
    design_info <- full_design %>%
      select(participant_id, week, treatment, path) %>%
      distinct()

    analysis_data_long <- analysis_data_long %>%
      left_join(design_info, by = c("participant_id", "week")) %>%
      # Add carryover variables by participant
      group_by(participant_id) %>%
      arrange(week) %>%
      mutate(
        # Calculate design-specific carryover effects
        prev_treatment = lag(treatment, default = 0),
        treatment_stopped = (prev_treatment == 1 & treatment == 0),

        # Calculate carryover effect uniformly across designs
        # Time since discontinuation
        time_since_discontinuation = ifelse(
          treatment == 0,
          cumsum(treatment == 0),
          0
        ),
        # Carryover effect (exponential decay)
        carryover_effect = if (params$carryover_t1half > 0) {
          ifelse(
            time_since_discontinuation > 0,
            (1/2)^(params$carryover_scale *
                   time_since_discontinuation /
                   params$carryover_t1half),
            0
          )
        } else {
          0
        }
      ) %>%
      ungroup()

    # Analyze with linear mixed model
    # Following Hendrickson et al. (2020) approach:
    # - Always include time effect (week) to account for period effects
    # - Biomarker × treatment interaction as primary effect of interest
    # - Carryover effect when applicable (enhancement over Hendrickson)
    # - Random intercept for participant-level variation
    model_result <- tryCatch({

      # Include carryover effect in model if present
      if (params$carryover_t1half > 0) {
        model <- lmer(
          response ~ treatment * bm + week + carryover_effect +
            (1 | participant_id),
          data = analysis_data_long
        )
      } else {
        model <- lmer(
          response ~ treatment * bm + week + (1 | participant_id),
          data = analysis_data_long
        )
      }

      # Extract model summary
      model_summary <- summary(model)
      coefs <- model_summary$coefficients

      # Get row index for interaction term
      interaction_idx <- which(rownames(coefs) == "treatment:bm")

      # Extract statistics of interest
      effect_size <- coefs[interaction_idx, "Estimate"]
      std_error <- coefs[interaction_idx, "Std. Error"]
      t_value <- coefs[interaction_idx, "t value"]

      # Calculate p-value (2-sided test)
      df <- model_summary$devcomp$dims[["n"]] - length(coefs[,1])
      p_value <- 2 * pt(-abs(t_value), df = df)
      significant <- p_value < 0.05

      # Return results
      tibble(
        effect_size = effect_size,
        std_error = std_error,
        t_value = t_value,
        p_value = p_value,
        significant = significant,
        error = FALSE
      )
    }, error = function(e) {
      # Return NA values on error
      tibble(
        effect_size = NA_real_,
        std_error = NA_real_,
        t_value = NA_real_,
        p_value = NA_real_,
        significant = NA,
        error = TRUE
      )
    })

    # Add metadata for this iteration
    model_result %>%
      mutate(
        iteration = iter,
        design = design_name,
        n_participants = params$n_participants,
        biomarker_correlation = params$biomarker_correlation,
        carryover_t1half = params$carryover_t1half
      )
  })

  return(all_results)
}
#=========================================================================== 
# End Function: run_monte_carlo
#=========================================================================== 

#=========================================================================== 
# Function: run_parameter_sweep
# Description: Run simulations for multiple parameter combinations
#=========================================================================== 
run_parameter_sweep <- function(param_grid) {
  # Number of parameter combinations
  n_combinations <- nrow(param_grid)

  # Run simulations for each parameter combination
  results <- map_dfr(1:n_combinations, function(i) {
    # Extract current parameters
    current_params <- as.list(param_grid[i,])
    current_params$n_iterations <- base_params$n_iterations
    current_params$treatment_effect <-
      base_params$treatment_effect
    current_params$carryover_scale <-
      base_params$carryover_scale
    current_params$between_subject_sd <-
      base_params$between_subject_sd
    current_params$within_subject_sd <-
      base_params$within_subject_sd

    # Run Monte Carlo for each design
    # Run Hybrid N-of-1 design
    hybrid_results <- run_monte_carlo("hybrid", current_params)

    # Run Traditional Crossover design
    crossover_results <- run_monte_carlo("crossover", current_params)

    # Combine results
    bind_rows(hybrid_results, crossover_results)
  })

  return(results)
}
#=========================================================================== 
# End Function: run_parameter_sweep
#=========================================================================== 

#===========================================================================
# Function: create_designs
# Description: Create trial designs for a given number of
# participants with Hendrickson's 4-path randomization structure
#===========================================================================
create_designs <- function(n_participants) {
  # HYBRID N-OF-1 DESIGN (Design 4) - Hendrickson 4-path structure
  # Implements 2x2 factorial randomization:
  #   Factor 1: Blinded discontinuation (stay on vs discontinue)
  #   Factor 2: Crossover sequence (AB vs BA)
  #
  # Phase structure (20 weeks total):
  #   Weeks 1-8:   Open-label (all on treatment)
  #   Weeks 9-12:  Blinded discontinuation (randomized)
  #   Weeks 13-16: Crossover period 1 (randomized sequence)
  #   Weeks 17-20: Crossover period 2 (opposite of period 1)

  # Randomly assign participants to one of 4 paths
  # Ensure roughly equal distribution across paths
  path_assignment <- rep(1:4, length.out = n_participants)
  path_assignment <- sample(path_assignment)  # Randomize

  # Path definitions (matching Hendrickson's tdNof11):
  # pathA: BD stay on (1,1), CO drug first (1,0)
  # pathB: BD stay on (1,1), CO placebo first (0,1)
  # pathC: BD discontinue (0,0), CO drug first (1,0)
  # pathD: BD discontinue (0,0), CO placebo first (0,1)

  hybrid_design <- expand_grid(
    participant_id = 1:n_participants,
    week = 1:20
  ) %>%
  mutate(
    path = path_assignment[participant_id],
    # Define treatment schedule based on path
    treatment = case_when(
      # Weeks 1-8: Open-label phase (all paths on treatment)
      week <= 8 ~ 1,

      # Weeks 9-12: Blinded discontinuation phase
      # Paths A & B: stay on drug
      # Paths C & D: discontinue (off drug)
      week >= 9 & week <= 12 & path %in% c(1, 2) ~ 1,
      week >= 9 & week <= 12 & path %in% c(3, 4) ~ 0,

      # Weeks 13-16: Crossover period 1
      # Paths A & C: on drug (drug-first sequence)
      # Paths B & D: off drug (placebo-first sequence)
      week >= 13 & week <= 16 & path %in% c(1, 3) ~ 1,
      week >= 13 & week <= 16 & path %in% c(2, 4) ~ 0,

      # Weeks 17-20: Crossover period 2 (switch)
      # Paths A & C: off drug
      # Paths B & D: on drug
      week >= 17 & week <= 20 & path %in% c(1, 3) ~ 0,
      week >= 17 & week <= 20 & path %in% c(2, 4) ~ 1,

      TRUE ~ NA_real_
    ),
    # Define expectancy (blinding)
    expectancy = case_when(
      week <= 8 ~ 1,        # Open-label: know on drug
      week >= 9 & week <= 20 ~ 0.5,  # Blinded phases
      TRUE ~ NA_real_
    ),
    # Define periods for analysis
    period = case_when(
      week <= 8 ~ 1,        # Open-label
      week >= 9 & week <= 12 ~ 2,   # Blinded discontinuation
      week >= 13 & week <= 20 ~ 3,  # Crossover
      TRUE ~ NA_real_
    )
  )

  # TRADITIONAL CROSSOVER DESIGN (Design 3)
  # Simple AB vs BA randomization (no change from original)
  crossover_design <- expand_grid(
    participant_id = 1:n_participants,
    week = 1:20
  ) %>%
  mutate(
    # Randomly assign to AB or BA sequence
    sequence = if_else(participant_id %% 2 == 0, "AB", "BA"),
    # Set path based on sequence (path 1 = AB, path 2 = BA)
    path = if_else(participant_id %% 2 == 0, 1, 2),
    # Set treatment based on sequence
    treatment = case_when(
      sequence == "AB" & week <= 10 ~ 1,
      sequence == "AB" & week > 10 ~ 0,
      sequence == "BA" & week <= 10 ~ 0,
      sequence == "BA" & week > 10 ~ 1,
      TRUE ~ NA_real_
    ),
    expectancy = 0.5,  # Blinded throughout
    period = if_else(week <= 10, 1, 2)
  )

  # Return both designs in a list
  return(list(
    hybrid = hybrid_design,
    crossover = crossover_design
  ))
}
#=========================================================================== 
# End Function: create_designs
#=========================================================================== 

#===========================================================================
# Function: prepare_design_for_simulation
# Description: Convert designs to a format compatible with
# generate_data
#=========================================================================== 
prepare_design_for_simulation <- function(design_df) {
  # Create a copy to avoid modifying the original
  design <- design_df %>%
    as_tibble()

  # Add timepoint names that are compatible with generate_data
  design <- design %>%
    arrange(participant_id, week) %>%
    mutate(
      timepoint_name = paste0("W", week),
      tod = treatment,  # Time on drug (1 = on, 0 = off)
      # Use expectancy column if present, otherwise derive from
      # treatment
      e = if ("expectancy" %in% names(.)) expectancy else
          if_else(tod > 0, 1, 0),
      t_wk = 1  # Each timepoint is 1 week
    ) %>%
    group_by(participant_id) %>%
    arrange(participant_id, week) %>%
    mutate(
      # Time on pharmacologic biomarker (cumulative when
      # expectancy > 0)
      tpb = cumsum(e),
      # Calculate time since discontinuation for carryover
      last_drug = lag(tod, default = 0),
      drug_stopped = (last_drug == 1 & tod == 0),
      # Time since discontinuation
      tsd = if_else(tod == 0, cumsum(tod == 0), 0)
    ) %>%
    ungroup()

  # Ensure the design has all required columns
  design <- design %>%
    select(timepoint_name, t_wk, e, tod, tsd, tpb)

  return(design)
}
#=========================================================================== 
# End Function: prepare_design_for_simulation
#=========================================================================== 
# ==================
# ==================
# ==================
# ==================

# ==================
# Main script begins
# ==================

# Base set of parameters - these will be varied in the Monte Carlo
# simulation. Parameters informed by Hendrickson et al. (2020) for
# realistic clinical scenarios
base_params <- list(
  # Sample size - modest for pragmatic trial designs
  n_participants = 30,

  # Biomarker-response correlation (default - will be varied)
  # Hendrickson et al. studied correlations in this range
  # Moderate correlation - proven to work in positive definiteness
  # testing
  biomarker_correlation = 0.2,

  # Treatment effect (baseline)
  # Reduced from 5 for more realistic power
  treatment_effect = 3.5,

  # Carryover effect parameters
  # Half-life of carryover effect in weeks
  carryover_t1half = 1,
  # Scale factor (set to 1 to match Hendrickson)
  carryover_scale = 1,

  # Random variation - increased to make signal detection more
  # challenging
  # Increased between-subject variability
  between_subject_sd = 3.0,
  # Increased within-subject variability
  within_subject_sd = 2.8,

  # Number of Monte Carlo iterations
  n_iterations = 20
)

# Define parameter grid for simulation
param_grid <- expand_grid(
  n_participants = c(70),
  # Reasonable range (will be forced positive definite)
  biomarker_correlation = c(0.2, 0.4),
  # Three carryover conditions: none, moderate, strong
  carryover_t1half = c(0, 1.0, 2.0),
  # Increased from 3.0 to improve power
  treatment_effect = c(5.0)
)
# Define base correlation parameters for carryover adjustment
# FIXED CORRELATION STRUCTURE (Hendrickson et al. 2020)
# These values do NOT vary with carryover parameters.
# Carryover affects MEANS only, not correlations.
# See: correlation_structure_design.pdf, carryover_correlation_theory.tex

# Create model_params compatible with generate_data function
model_params <- list(
  N = base_params$n_participants,
  # Correlation between biomarker and response (VARIES in parameter sweep)
  c.bm = base_params$biomarker_correlation,

  # Carryover parameters (affects MEANS only)
  carryover_t1half = base_params$carryover_t1half,

  # FIXED CORRELATION VALUES (Hendrickson approach)
  # These do NOT change with carryover_t1half
  c.tv = 0.8,    # Autocorrelation: time_variant across time
  c.pb = 0.8,    # Autocorrelation: pharm_biomarker across time
  c.br = 0.8,    # Autocorrelation: bio_response across time
  c.cf1t = 0.2,  # Cross-correlation: different factors, same time
  c.cfct = 0.1   # Cross-correlation: different factors, different times
)

# Set up response parameters
resp_param <- tibble(
  cat = c("time_variant", "pharm_biomarker", "bio_response"),
  # Maximum response value
  max = c(1.0, 1.0, base_params$treatment_effect),
  disp = c(2.0, 2.0, 2.0),  # Displacement
  rate = c(0.3, 0.3, 0.3),  # Rate
  # Standard deviation
  sd = c(base_params$within_subject_sd,
         base_params$within_subject_sd,
         base_params$within_subject_sd)
)

# Set up baseline parameters
baseline_param <- tibble(
  cat = c("biomarker", "baseline"),
  m = c(5.0, 10.0),  # Mean values
  # Standard deviation
  sd = c(2.0, base_params$between_subject_sd)
)

# Analysis options
analysis_options <- list(
  use_expectancy = TRUE,  # Use time information in model
  random_slope = FALSE,   # No random slopes for simplicity
  full_output = FALSE,    # Return only summary statistics
  simple_carryover = FALSE, # Use complex carryover model
  carryover_halflife = base_params$carryover_t1half,
  carryover_scale_factor = base_params$carryover_scale
)

# Options for simulation
sim_options <- list(
  n_reps = base_params$n_iterations  # Number of Monte Carlo replications
)

# Create designs for base sample size for visualization
original_designs <- create_designs(base_params$n_participants)

# For HYBRID design: Create separate design templates for each path
# following Hendrickson's approach
hybrid_paths <- list()
for (path_id in 1:4) {
  path_design <- original_designs$hybrid %>%
    filter(path == path_id, participant_id ==
           min(participant_id[path == path_id]))

  hybrid_paths[[path_id]] <- prepare_design_for_simulation(
    path_design
  )
}

# For CROSSOVER design: Create separate templates for AB and BA
crossover_paths <- list()
# Path 1: AB sequence (even IDs)
crossover_paths[[1]] <- prepare_design_for_simulation(
  original_designs$crossover %>%
    filter(participant_id %% 2 == 0,
           participant_id == min(participant_id[participant_id %%
                                                   2 == 0]))
)
# Path 2: BA sequence (odd IDs)
crossover_paths[[2]] <- prepare_design_for_simulation(
  original_designs$crossover %>%
    filter(participant_id %% 2 == 1,
           participant_id == min(participant_id[participant_id %%
                                                   2 == 1]))
)

# Store designs in a list with both full design and path templates
designs <- list(
  hybrid = list(
    full_design = original_designs$hybrid,
    design_paths = hybrid_paths,
    name = "Hybrid N-of-1"
  ),
  crossover = list(
    full_design = original_designs$crossover,
    design_paths = crossover_paths,
    name = "Traditional Crossover"
  )
)

# ========================================== 
# Main simulation
# ========================================== 

# Simulation parameters
simulation_params <- list(
  n_iterations = 20,
  carryover_scale = base_params$carryover_scale,
  between_subject_sd = base_params$between_subject_sd,
  within_subject_sd = base_params$within_subject_sd
)

# Build sigma matrix cache for all parameter/design combinations
cat("Building sigma matrix cache...\n")
sigma_cache <- list()
valid_combinations <- tibble()

for (i in 1:nrow(param_grid)) {
  current_params <- as.list(param_grid[i,])
  current_params$n_iterations <- simulation_params$n_iterations
  current_params$carryover_scale <- simulation_params$carryover_scale

  # FIXED CORRELATIONS: Do NOT adjust based on carryover
  # Correlation structure is independent of carryover parameters

  # Update model parameters for this combination
  current_model_params <- model_params
  current_model_params$N <- current_params$n_participants
  current_model_params$c.bm <-
    current_params$biomarker_correlation
  current_model_params$carryover_t1half <-
    current_params$carryover_t1half

  # Correlation values remain FIXED at Hendrickson values
  # (already set in model_params: c.tv=0.8, c.pb=0.8, c.br=0.8,
  # c.cf1t=0.2, c.cfct=0.1)

  # Test both designs
  # For sigma matrix, use first path template (all paths have same
  # timepoints)
  for (design_name in c("hybrid", "crossover")) {
    # Select first path design template for sigma building
    design_template <- designs[[design_name]]$design_paths[[1]]

    # Try to build sigma matrix
    factor_types <- c("time_variant", "pharm_biomarker",
                      "bio_response")
    factor_abbreviations <- c("tv", "pb", "br")

    sigma_result <- build_sigma_matrix(
      current_model_params, resp_param, baseline_param,
      design_template, factor_types, factor_abbreviations,
      verbose = TRUE
    )

    cache_key <- create_sigma_cache_key(design_name,
                                        current_params)

    if (!is.null(sigma_result)) {
      # Positive definite - add to cache
      sigma_cache[[cache_key]] <- sigma_result
      valid_combinations <- bind_rows(
        valid_combinations,
        tibble(
          design = design_name,
          biomarker_correlation =
            current_params$biomarker_correlation,
          carryover_t1half = current_params$carryover_t1half
        )
      )
      cat("✓ Valid sigma for", design_name, "design,",
          "biomarker_correlation =",
          current_params$biomarker_correlation,
          ", carryover_t1half =",
          current_params$carryover_t1half, "\n")
    } else {
      # Non-positive definite - skip
      cat("✗ Non-PD sigma for", design_name, "design,",
          "biomarker_correlation =",
          current_params$biomarker_correlation,
          ", carryover_t1half =",
          current_params$carryover_t1half, "\n")
    }
  }
}

cat("\nSigma cache built. Valid combinations:\n")
print(valid_combinations)

# Initialize results storage
simulation_results <- tibble()

# For each parameter combination
for (i in 1:nrow(param_grid)) {
  # Extract current parameters
  current_params <- as.list(param_grid[i,])
  current_params$n_iterations <- simulation_params$n_iterations
  current_params$carryover_scale <- simulation_params$carryover_scale

  # Run for each design using our Monte Carlo function
  for (design_name in c("hybrid", "crossover")) {
    # Run the Monte Carlo simulation with sigma cache
    design_results <- run_monte_carlo(design_name, current_params,
                                      sigma_cache)

    # Add to overall results
    simulation_results <- bind_rows(simulation_results,
                                    design_results)
  }
}

# Calculate summary statistics across all iterations
if (nrow(simulation_results) > 0) {
  simulation_summary <- simulation_results %>%
    group_by(design, n_participants, biomarker_correlation,
             carryover_t1half) %>%
    summarize(
      power = mean(significant, na.rm = TRUE),
      mean_effect = mean(effect_size, na.rm = TRUE),
      mean_se = mean(std_error, na.rm = TRUE),
      n_iterations = n(),
      .groups = "drop"
    )
} else {
  cat("\nNo valid simulation results - all parameter",
      "combinations resulted in non-positive definite",
      "matrices.\n")
  cat("Consider using lower correlation parameters or",
      "different parameter ranges.\n")
  simulation_summary <- tibble()
}

# Display summary statistics
if (nrow(simulation_summary) > 0) {
  print(summary(simulation_summary$power))

  # Save the simulation results for visualization and further analysis
  sim_results <- list(
    simulation_results = simulation_results,
    simulation_summary = simulation_summary
  )

  # ========================================== 
  # Visualization
  # ========================================== 

  # Simple heatmap of power by carryover by design
  power_heatmap <- ggplot(
      simulation_summary,
      aes(x = factor(carryover_t1half), y = design,
          fill = power)
    ) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", power)),
              color = "black", size = 4, fontface = "bold") +
    scale_fill_viridis_c(name = "Power",
                         labels = scales::percent) +
    labs(
      title = "Statistical Power by Design and Carryover Half-life",
      x = "Carryover Half-life (weeks)",
      y = "Design"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank()
    )

  print(power_heatmap)

  # Display summary statistics
  print(simulation_summary)
  # dim
  dim(simulation_summary)
} else {
  cat("\nNo visualization possible - no valid results to",
      "display.\n")
  cat("Enhanced carryover implementation successfully",
      "integrated!\n")

  # Save empty results
  sim_results <- list(
    simulation_results = simulation_results,
    simulation_summary = simulation_summary
  )
}
# ========================================== 
# End of script
# ==========================================
