# N-of-1 Trial Simulation - Clustered/8-Point Measurement Designs
# Hybrid, OL+BDC, and Crossover designs with 8 measurement points each
# Hybrid/OL+BDC have clustered measurements; Crossover is evenly spaced

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(lme4)        # Strategy 4: faster lmer alternative
  library(MASS)
  library(zoo)
  library(patchwork)
  library(conflicted)
  library(future)      # Strategy 1: parallel processing
  library(furrr)       # Strategy 1: parallel purrr
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
n_iterations <- 500

# ============================================================================
# OPTIMIZATION SETTINGS (Strategy 5: exploratory mode)
# ============================================================================
# Set exploratory_mode = TRUE for faster runs during development
# Set exploratory_mode = FALSE for final publication-quality results
# Can also be set via environment variable: EXPLORATORY_MODE=TRUE
# Can also override iterations directly: N_ITERATIONS=100

exploratory_mode <- as.logical(Sys.getenv("EXPLORATORY_MODE", "FALSE"))
exploratory_iterations <- 20  # Reduced iterations for exploratory runs

# Check for N_ITERATIONS environment variable override
env_iterations <- Sys.getenv("N_ITERATIONS", "")
if (nchar(env_iterations) > 0) {
  n_iterations <- as.integer(env_iterations)
  cat("*** CUSTOM ITERATIONS:", n_iterations, "***\n\n")
} else if (exploratory_mode) {
  n_iterations <- exploratory_iterations
  cat("*** EXPLORATORY MODE: Running with", n_iterations, "iterations ***\n\n")
}

# Strategy 1: Parallel processing setup
# Detect available cores (leave 1 for system)
n_cores <- max(1, parallel::detectCores() - 1)
use_parallel <- TRUE  # Set to FALSE to disable parallel processing

# Strategy 4: Model fitting engine
# Options: "lme" (nlme::lme) or "lmer" (lme4::lmer)
# Use lme with corCAR1 for correct AR(1) correlation structure matching the DGP
model_engine <- "lme"

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
# Reduced to improve matrix conditioning and positive definiteness
c.br <- 0.75  # BR autocorrelation (reduced from 0.8)
c.er <- 0.75  # ER autocorrelation (reduced from 0.8)
c.tr <- 0.75  # TR autocorrelation (reduced from 0.8)
c.cf1t <- 0.12  # Same-time cross-correlation (reduced from 0.2)
c.cfct <- 0.05  # Different-time cross-correlation (reduced from 0.1)
c.bm_baseline <- 0.25  # Biomarker-baseline correlation (reduced from 0.3)
c.baseline_resp <- 0.3 # Baseline-response correlation (reduced from 0.4)

# Allowed correlation values (finite grid for consistent results)
allowed_correlations <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

# ============================================================================
# CARRYOVER PARAMETERIZATION
# ============================================================================
# Unified carryover specification supporting both:
#   - "halflife": exponential decay with t_half in weeks
#   - "proportion": fixed proportion retained at first off-drug timepoint
#
# The exponential decay model: retention = (1/2)^(tsd / t_half)
# where tsd = time since discontinuation (weeks)
#
# Conversion between parameterizations:
#   Given proportion p retained at time t: t_half = -t / log2(p)
#   Given halflife h, proportion at time t: p = (1/2)^(t/h)

carryover_mode <- "halflife"  # "halflife" or "proportion"

# Carryover function: computes retention given time since discontinuation
calculate_carryover <- function(tsd, carryover_param, mode = carryover_mode) {
 if (carryover_param == 0 || is.na(tsd) || tsd <= 0) {

    0
  } else if (mode == "halflife") {
    (1/2)^(tsd / carryover_param)
  } else if (mode == "proportion") {
    carryover_param
  } else {
    stop("Unknown carryover mode: ", mode)
  }
}

# Vectorized version for use in mutate
calculate_carryover_vec <- Vectorize(calculate_carryover)

# Parameter grid - Hybrid, OL+BDC, and Crossover designs
# Now includes model_carryover for two-scenario comparison:
#   - model_carryover = TRUE: include carryover term in analysis model
#   - model_carryover = FALSE: omit carryover term (Hendrickson approach)
param_grid <- bind_rows(
  # Open-label + Blinded Discontinuation design (power conditions)
  expand_grid(
    design = "ol_bdc",
    n_participants = c(35, 70),
    biomarker_moderation = c(0.2, 0.4, 0.6),
    biomarker_correlation = c(0.3),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE, FALSE)
  ),
  # Type I error condition for OL+BDC (null: no biomarker moderation)
  expand_grid(
    design = "ol_bdc",
    n_participants = c(35, 70),
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE, FALSE)
  ),
  # Hybrid design (power conditions)
  expand_grid(
    design = "hybrid",
    n_participants = c(35, 70),
    biomarker_moderation = c(0.2, 0.4, 0.6),
    biomarker_correlation = c(0.3),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE, FALSE)
  ),
  # Type I error condition for Hybrid (null: no biomarker moderation)
  expand_grid(
    design = "hybrid",
    n_participants = c(35, 70),
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE, FALSE)
  ),
  # Crossover design (power conditions) - Hendrickson Figure 2C
  expand_grid(
    design = "crossover",
    n_participants = c(35, 70),
    biomarker_moderation = c(0.2, 0.4, 0.6),
    biomarker_correlation = c(0.3),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE, FALSE)
  ),
  # Type I error condition for Crossover (null: no biomarker moderation)
  expand_grid(
    design = "crossover",
    n_participants = c(35, 70),
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE, FALSE)
  ),
  # Open-Label design (power conditions) - Hendrickson Figure 2A
  # NOTE: OL uses only model_carryover=TRUE because treatment=1 for all obs
  # (no treatment variance, so treatment*bm_centered interaction is unestimable)
  expand_grid(
    design = "ol",
    n_participants = c(35, 70),
    biomarker_moderation = c(0.2, 0.4, 0.6),
    biomarker_correlation = c(0.3),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE)
  ),
  # Type I error condition for OL (null: no biomarker moderation)
  expand_grid(
    design = "ol",
    n_participants = c(35, 70),
    biomarker_moderation = c(0),
    biomarker_correlation = c(0),
    carryover_param = c(0, 0.1, 0.2),  # Hendrickson values
    model_carryover = c(TRUE)
  )
) |>
  # Remove redundant rows: when carryover_param = 0, model_carryover doesn't matter
  filter(!(carryover_param == 0 & model_carryover == FALSE)) |>
  rename(carryover_halflife = carryover_param)

# ============================================================================
# BUILD SIGMA WITH GUARANTEED PD (Time-Based AR(1))
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

  Sigma_BR <- build_ar1_time(weeks, c.br, within_subject_sd)
  Sigma_ER <- build_ar1_time(weeks, c.er, within_subject_sd)
  Sigma_TR <- build_ar1_time(weeks, c.tr, within_subject_sd)

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

  # STAGE 3: Build Sigma_12 (no on-the-fly adjustment - parameters are pre-validated)
  effective_c.bm <- c.bm
  effective_c.baseline <- c.baseline_resp

  Sigma_12 <- matrix(0, 3 * n_tp, 2)
  Sigma_12[br_idx, 1] <- effective_c.bm * within_subject_sd * biomarker_sd
  Sigma_12[er_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[tr_idx, 1] <- effective_c.bm * 0.5 * within_subject_sd * biomarker_sd
  Sigma_12[br_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[er_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd
  Sigma_12[tr_idx, 2] <- effective_c.baseline * within_subject_sd * between_subject_sd

  # Verify positive definiteness (should always pass if parameters were validated correctly)
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
        "This indicates either:",
        "  1. A bug in the validate_parameter_grid() function, or",
        "  2. Inconsistent covariance matrix construction between validation and simulation",
        "",
        "Please investigate the mismatch between validate_parameter_grid()",
        "and build_sigma_guaranteed_pd().",
        sep = "\n"
      ),
      effective_c.bm, min_eig))
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

  # Strategy 3: Pre-compute Sigma_22_inv for use in generate_participant_twostage
  Sigma_22_inv <- solve(Sigma_22)

  # Also pre-compute the conditional mean transformation matrix
  cond_mean_transform <- Sigma_12 %*% Sigma_22_inv

  list(
    Sigma = Sigma,
    Sigma_22 = Sigma_22,
    Sigma_22_inv = Sigma_22_inv,           # Strategy 3: cached inverse
    Sigma_cond = Sigma_cond,
    Sigma_12 = Sigma_12,
    cond_mean_transform = cond_mean_transform,  # Strategy 3: cached transform
    indices = indices,
    effective_c.bm = effective_c.bm
  )
}

# ============================================================================
# TWO-STAGE GENERATION
# ============================================================================

generate_participant_twostage <- function(sigma_parts, idx) {
  Sigma_22 <- sigma_parts$Sigma_22
  Sigma_cond <- sigma_parts$Sigma_cond
  n_tp <- idx$n_tp

  # Stage 1: Generate (biomarker, baseline)
  stage1 <- mvrnorm(1, mu = c(biomarker_mean, baseline_mean), Sigma = Sigma_22)
  biomarker <- stage1[1]
  baseline <- stage1[2]

  # Stage 2: Generate responses conditional on stage 1
  # Strategy 3: Use pre-computed transform instead of solve() each time
  z <- c(biomarker - biomarker_mean, baseline - baseline_mean)
  cond_mean <- as.vector(sigma_parts$cond_mean_transform %*% z)
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
# DESIGN STRUCTURES - CLUSTERED SCHEDULES
# ============================================================================

# 8 points with clusters for Hybrid and OL+BDC
measurement_weeks_hybrid <- c(4, 8, 9, 10, 11, 12, 16, 20)
measurement_weeks_ol_bdc <- c(4, 8, 12, 16, 17, 18, 19, 20)
# 8 points evenly spaced for Crossover (Hendrickson Figure 2C)
measurement_weeks_crossover <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)
# 8 points evenly spaced for Open-Label (Hendrickson Figure 2A)
measurement_weeks_ol <- c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)

create_ol_bdc_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        week <= 16 ~ 1,
        week == 17 ~ 1,
        week == 18 & path %in% c(1, 2) ~ 1,
        week == 18 & path %in% c(3, 4) ~ 0,
        week %in% c(19, 20) ~ 0,
        TRUE ~ 1
      ),
      expectancy = if_else(week <= 16, 1.0, 0.5)
    )
}

create_hybrid_design <- function(n_participants, measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
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

create_crossover_design <- function(n_participants, measurement_weeks) {
  # Traditional crossover: 2 sequences (AB and BA)
  # Period 1 (weeks <= 10): sequence AB = drug, sequence BA = placebo
  # Period 2 (weeks > 10): sequence AB = placebo, sequence BA = drug
  path_assignment <- sample(rep(1:2, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        # Period 1: weeks 2.5, 5, 7.5, 10
        week <= 10 & path == 1 ~ 1,  # AB sequence: drug first
        week <= 10 & path == 2 ~ 0,  # BA sequence: placebo first
        # Period 2: weeks 12.5, 15, 17.5, 20
        week > 10 & path == 1 ~ 0,   # AB sequence: placebo second
        week > 10 & path == 2 ~ 1,   # BA sequence: drug second
        TRUE ~ NA_real_
      ),
      expectancy = 0.5  # Blinded throughout
    )
}

create_ol_design <- function(n_participants, measurement_weeks) {
  # Open-Label design (Hendrickson Figure 2A): all participants on drug
  # No randomization needed - everyone gets active treatment
  # expectancy = 1.0 (open-label, participants know they are on drug)
  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = 1,  # Single path for all
      treatment = 1,  # Always on drug
      expectancy = 1.0  # Open-label: full expectancy effect
    )
}

# ============================================================================
# PRE-SIMULATION SIGMA VALIDATION (3 unique structures)
# ============================================================================

cat("\n")
cat(paste0("=", strrep("=", 79), "\n"))
cat("PRE-SIMULATION SIGMA MATRIX VALIDATION\n")
cat(paste0("=", strrep("=", 79), "\n\n"))

# Get unique biomarker correlations in param_grid
unique_bm_corr <- unique(param_grid$biomarker_correlation)

# Define sigma structures to validate: (design, measurement_weeks, c.bm)
sigma_structures <- list(
  list(design = "hybrid", weeks = measurement_weeks_hybrid, c_bm = 0.3),
  list(design = "hybrid", weeks = measurement_weeks_hybrid, c_bm = 0),
  list(design = "ol_bdc", weeks = measurement_weeks_ol_bdc, c_bm = 0.3),
  list(design = "ol_bdc", weeks = measurement_weeks_ol_bdc, c_bm = 0),
  list(design = "crossover", weeks = measurement_weeks_crossover, c_bm = 0.3),
  list(design = "crossover", weeks = measurement_weeks_crossover, c_bm = 0),
  list(design = "ol", weeks = measurement_weeks_ol, c_bm = 0.3),
  list(design = "ol", weeks = measurement_weeks_ol, c_bm = 0)
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
      params = list()  # params passed but not used in PD check
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
       "   Check parameter grid or measurement schedule.")
}

cat(paste0("=", strrep("=", 79), "\n"))
cat("✓ All sigma matrices valid - proceeding with simulation\n")
cat(paste0("=", strrep("=", 79), "\n\n"))

# ============================================================================
# STRATEGY 2: BUILD AND CACHE SIGMA MATRICES
# ============================================================================
# Pre-build all unique sigma structures to avoid redundant computation

cat("Building and caching sigma matrices...\n")

# Create cache with named list for fast lookup
sigma_cache <- list()

for (i in seq_along(sigma_structures)) {
  sigma_def <- sigma_structures[[i]]
  cache_key <- paste(sigma_def$design, sigma_def$c_bm, sep = "_")

  sigma_cache[[cache_key]] <- list(
    sigma_parts = build_sigma_guaranteed_pd(
      weeks = sigma_def$weeks,
      c.bm = sigma_def$c_bm,
      params = list()
    ),
    weeks = sigma_def$weeks
  )

  cat(sprintf("  Cached: %s (c.bm=%.1f)\n", sigma_def$design, sigma_def$c_bm))
}

# Helper function to retrieve cached sigma (takes cache as parameter for parallel)
get_cached_sigma <- function(design, c_bm, cache) {
  cache_key <- paste(design, c_bm, sep = "_")
  if (!cache_key %in% names(cache)) {
    stop("Sigma not found in cache: ", cache_key)
  }
  cache[[cache_key]]
}

# Pre-generate designs (they can be reused with different seeds for path_assignment)
# Actually, designs need fresh random path assignments per iteration, so we keep
# the create_*_design functions as-is

cat("Sigma cache ready with", length(sigma_cache), "structures\n\n")

# ============================================================================
# RUN SIMULATION (Optimized with Strategies 1-5)
# ============================================================================

# Set up logging
log_file <- "../output/simulation_clustered_log.txt"
if (!dir.exists("../output")) dir.create("../output", recursive = TRUE)
sink(log_file, split = TRUE)

cat("Running CLUSTERED designs simulation...\n")
cat("Designs: Open-Label, OL+BDC, Hybrid, Crossover (8 measurement points each)\n")
cat("Carryover mode:", carryover_mode, "\n")
cat("Treatment effect =", treatment_effect, "\n")
cat("Participants = 35, 70 (from param_grid)\n")
cat("Iterations per condition =", n_iterations, "\n")
cat("Total conditions =", nrow(param_grid), "\n")
cat("Model engine:", model_engine, "\n")
cat("Parallel processing:", if (use_parallel) paste0("YES (", n_cores, " cores)") else "NO", "\n")
if (exploratory_mode) cat("*** EXPLORATORY MODE ***\n")
cat("\n")

# ----------------------------------------------------------------------------
# Define single iteration function (for parallel execution)
# ----------------------------------------------------------------------------

run_single_iteration <- function(iter, params, sigma_cache, model_engine) {
  # Set seed for reproducibility
  set.seed(iter * 1000 + params$condition_id)

  # Get cached sigma (Strategy 2)
  cached <- get_cached_sigma(params$design, params$biomarker_correlation, sigma_cache)
  sigma_parts <- cached$sigma_parts
  measurement_weeks <- cached$weeks
  idx <- sigma_parts$indices
  effective_bm_corr <- sigma_parts$effective_c.bm

  # Get n_participants from params
  local_n <- params$n_participants

  # Create design
  if (params$design == "hybrid") {
    trial_design <- create_hybrid_design(local_n, measurement_weeks)
  } else if (params$design == "crossover") {
    trial_design <- create_crossover_design(local_n, measurement_weeks)
  } else if (params$design == "ol") {
    trial_design <- create_ol_design(local_n, measurement_weeks)
  } else {
    trial_design <- create_ol_bdc_design(local_n, measurement_weeks)
  }

  # Generate participant data (Strategy 3: uses pre-computed Sigma_22_inv)
  all_participant_data <- lapply(1:local_n, function(pid) {
    result <- generate_participant_twostage(sigma_parts, idx)
    tibble(
      participant_id = pid,
      biomarker = result$biomarker,
      baseline = result$baseline,
      br_random = result$br_random,
      er_random = result$er_random,
      tr_random = result$tr_random
    )
  })

  participant_data <- bind_rows(all_participant_data) |>
    group_by(participant_id) |>
    mutate(timepoint_idx = row_number()) |>
    ungroup()

  bm_mod <- params$biomarker_moderation
  carryover_t_half <- params$carryover_halflife

  # Build trial data with proper time-since-discontinuation carryover
  trial_data <- trial_design |>
    group_by(participant_id) |>
    mutate(timepoint_idx = row_number()) |>
    ungroup() |>
    left_join(participant_data, by = c("participant_id", "timepoint_idx")) |>
    group_by(participant_id) |>
    arrange(week) |>
    mutate(
      weeks_on_drug = cumsum(treatment),
      weeks_with_expectancy = cumsum(expectancy),
      weeks_in_trial = week - min(week),
      just_discontinued = treatment == 0 & lag(treatment, default = 0) == 1,
      discontinuation_week = if_else(just_discontinued, week, NA_real_),
      discontinuation_week = zoo::na.locf(discontinuation_week, na.rm = FALSE),
      tsd = if_else(
        treatment == 0 & !is.na(discontinuation_week),
        week - discontinuation_week,
        0
      ),
      bm_centered = (biomarker - biomarker_mean) / biomarker_sd,
      effective_BR_rate = BR_rate * (1 + bm_mod * bm_centered),
      # Cumulative BR model with exponential decay carryover
      br_on_drug = weeks_on_drug * effective_BR_rate,
      # Track BR at time of discontinuation
      br_at_discont = if_else(just_discontinued, lag(br_on_drug, default = 0), NA_real_),
      br_at_discont = zoo::na.locf(br_at_discont, na.rm = FALSE),
      # carryover_effect: decay factor (0-1) - used in both data gen and analysis
      # t½=0 means NO carryover: effect drops to 0 immediately when off drug
      # t½>0 means carryover: effect decays exponentially
      carryover_effect = case_when(
        treatment == 1 ~ 1,  # On drug: full effect
        treatment == 0 & carryover_t_half == 0 ~ 0,  # Off drug, no carryover
        treatment == 0 & carryover_t_half > 0 ~ (1/2)^(tsd / carryover_t_half),  # Off drug, decay
        TRUE ~ 0
      ),
      # BR when off drug = decayed accumulated effect (0 if no carryover)
      br_off_drug = if_else(
        treatment == 0 & !is.na(br_at_discont),
        br_at_discont * carryover_effect,
        0
      ),
      BR_mean = if_else(treatment == 1, br_on_drug, br_off_drug),
      # effective_drug_weeks: cumulative drug exposure accounting for carryover
      # On drug: weeks_on_drug (cumulative)
      # Off drug: weeks_at_discont * carryover_effect (decayed cumulative)
      weeks_at_discont = if_else(just_discontinued, lag(weeks_on_drug, default = 0), NA_real_),
      weeks_at_discont = zoo::na.locf(weeks_at_discont, na.rm = FALSE),
      # effective_drug_weeks: analysis predictor for carryover model
      # On drug: weeks_on_drug (cumulative)
      # Off drug: weeks_at_discont * carryover_effect (decayed; 0 if t½=0)
      effective_drug_weeks = case_when(
        treatment == 1 ~ weeks_on_drug,
        treatment == 0 & !is.na(weeks_at_discont) ~ weeks_at_discont * carryover_effect,
        TRUE ~ 0
      ),
      BR = BR_mean + br_random,
      ER_mean = weeks_with_expectancy * ER_rate,
      ER = ER_mean + er_random,
      TR_mean = weeks_in_trial * TR_rate,
      TR = TR_mean + tr_random,
      response = baseline + BR + ER + TR
    ) |>
    ungroup()

  # Standardize biomarker using sample mean and SD for analysis
  # This matches the standardization used in the DGP (line 556)
  trial_data <- trial_data |>
    mutate(bm_centered = (biomarker - mean(biomarker)) / sd(biomarker))

  # Fit analysis model - Strategy 4: choice of lme or lmer
  # Two scenarios:
  # 1. model_carryover=TRUE: use effective_drug_weeks which accounts for cumulative
  #    exposure including carryover decay
  # 2. model_carryover=FALSE: use weeks_on_drug (ignores carryover, treats off-drug as 0)
  model_result <- tryCatch({
    if (model_engine == "lmer") {
      # Strategy 4: Use lme4::lmer (faster)
      if (params$design == "ol") {
        # OL design: everyone always on drug, so effective_drug_weeks is collinear
        # with week. Use week * bm_centered directly (no separate + week term)
        model <- lmer(
          response ~ week * bm_centered + (1 | participant_id),
          data = trial_data,
          REML = TRUE
        )
        interaction_term <- "week:bm_centered"
      } else if (params$model_carryover) {
        # Use effective_drug_weeks: cumulative exposure with carryover decay
        # At t½=0, effective_drug_weeks = weeks_on_drug when on drug, 0 when off
        # At t½>0, effective_drug_weeks = weeks_on_drug when on drug, decayed when off
        model <- lmer(
          response ~ effective_drug_weeks * bm_centered + week +
            (1 | participant_id),
          data = trial_data,
          REML = TRUE
        )
        interaction_term <- "effective_drug_weeks:bm_centered"
      } else {
        # Model carryover = FALSE: use weeks_on_drug (doesn't decay when off drug)
        # This misspecifies the model when there IS carryover in the data
        model <- lmer(
          response ~ weeks_on_drug * bm_centered + week + (1 | participant_id),
          data = trial_data,
          REML = TRUE
        )
        interaction_term <- "weeks_on_drug:bm_centered"
      }

      # Extract coefficients from lmer
      coefs <- summary(model)$coefficients
      idx_coef <- which(rownames(coefs) == interaction_term)

      # lmer doesn't provide p-values directly; use Satterthwaite via lmerTest
      # or normal approximation for large samples
      est <- coefs[idx_coef, "Estimate"]
      se <- coefs[idx_coef, "Std. Error"]
      t_val <- coefs[idx_coef, "t value"]
      # Normal approximation for p-value (valid for n=70 participants)
      p_val <- 2 * pnorm(-abs(t_val))

      tibble(
        iteration = iter,
        design = params$design,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,
        carryover_halflife = params$carryover_halflife,
        model_carryover = params$model_carryover,
        effect_size = est,
        se = se,
        t_value = t_val,
        p_value = p_val,
        significant = p_val < 0.05
      )
    } else {
      # Original nlme::lme approach
      if (params$design == "ol") {
        # OL design: everyone always on drug, so effective_drug_weeks is collinear
        # with week. Use week * bm_centered directly (no separate + week term)
        model <- lme(
          response ~ week * bm_centered,
          random = ~ 1 | participant_id,
          correlation = corCAR1(form = ~ week | participant_id),
          data = trial_data,
          control = lmeControl(opt = "optim")
        )
        interaction_term <- "week:bm_centered"
      } else if (params$model_carryover) {
        model <- lme(
          response ~ effective_drug_weeks * bm_centered + week,
          random = ~ 1 | participant_id,
          correlation = corCAR1(form = ~ week | participant_id),
          data = trial_data,
          control = lmeControl(opt = "optim")
        )
        interaction_term <- "effective_drug_weeks:bm_centered"
      } else {
        model <- lme(
          response ~ weeks_on_drug * bm_centered + week,
          random = ~ 1 | participant_id,
          correlation = corCAR1(form = ~ week | participant_id),
          data = trial_data,
          control = lmeControl(opt = "optim")
        )
        interaction_term <- "weeks_on_drug:bm_centered"
      }

      coefs <- summary(model)$tTable
      idx_coef <- which(rownames(coefs) == interaction_term)

      tibble(
        iteration = iter,
        design = params$design,
        n_participants = local_n,
        biomarker_moderation = params$biomarker_moderation,
        biomarker_correlation = effective_bm_corr,
        carryover_halflife = params$carryover_halflife,
        model_carryover = params$model_carryover,
        effect_size = coefs[idx_coef, "Value"],
        se = coefs[idx_coef, "Std.Error"],
        t_value = coefs[idx_coef, "t-value"],
        p_value = coefs[idx_coef, "p-value"],
        significant = coefs[idx_coef, "p-value"] < 0.05
      )
    }
  }, error = function(e) {
    tibble(
      iteration = iter,
      design = params$design,
      n_participants = local_n,
      biomarker_moderation = params$biomarker_moderation,
      biomarker_correlation = effective_bm_corr,
      carryover_halflife = params$carryover_halflife,
      model_carryover = params$model_carryover,
      effect_size = NA_real_,
      se = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      significant = NA
    )
  })

  model_result
}

# ----------------------------------------------------------------------------
# Strategy 1: Parallel execution across iterations
# ----------------------------------------------------------------------------

# Add condition_id to param_grid for seed management
param_grid <- param_grid |>
  mutate(condition_id = row_number())

# Set up parallel backend
if (use_parallel && n_cores > 1) {
  plan(multisession, workers = n_cores)
  cat(sprintf("Parallel backend initialized with %d workers\n\n", n_cores))
} else {
  plan(sequential)
  cat("Running in sequential mode\n\n")
}

# Progress tracking
start_time <- Sys.time()
total_iterations <- nrow(param_grid) * n_iterations

# Run simulation - parallelize across conditions
results_list <- list()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf(
    "Condition %d/%d: design=%s, bm_mod=%.2f, carryover_t½=%.1f wks, model_carryover=%s\n",
    i, nrow(param_grid), params$design, params$biomarker_moderation,
    params$carryover_halflife, params$model_carryover
  ))

  # Parallel execution across iterations (Strategy 1)
  # Capture params in local scope for lambda
  local_params <- params
  local_sigma_cache <- sigma_cache
  local_model_engine <- model_engine

  condition_results <- future_map_dfr(
    1:n_iterations,
    function(iter) {
      run_single_iteration(iter, local_params, local_sigma_cache, local_model_engine)
    },
    .options = furrr_options(
      seed = TRUE,
      globals = list(
        run_single_iteration = run_single_iteration,
        get_cached_sigma = get_cached_sigma,
        generate_participant_twostage = generate_participant_twostage,
        create_hybrid_design = create_hybrid_design,
        create_ol_bdc_design = create_ol_bdc_design,
        create_crossover_design = create_crossover_design,
        create_ol_design = create_ol_design,
        local_params = local_params,
        local_sigma_cache = local_sigma_cache,
        local_model_engine = local_model_engine,
        n_participants = n_participants,
        biomarker_mean = biomarker_mean,
        baseline_mean = baseline_mean,
        biomarker_sd = biomarker_sd,
        BR_rate = BR_rate,
        ER_rate = ER_rate,
        TR_rate = TR_rate
      ),
      packages = c("tidyverse", "nlme", "lme4", "MASS", "zoo")
    )
  )

  results_list[[i]] <- condition_results

  # Progress update
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  completed <- i * n_iterations
  rate <- completed / elapsed
  remaining <- (total_iterations - completed) / rate
  cat(sprintf("  Completed in %.1f sec (ETA: %.1f sec remaining)\n",
              elapsed, remaining))
}

results <- bind_rows(results_list)

# Reset parallel backend
plan(sequential)

# ============================================================================
# SUMMARIZE AND SAVE
# ============================================================================

cat("\n", strrep("=", 70), "\n")
cat("RESULTS - CLUSTERED DESIGNS (Two-Scenario Comparison)\n")
cat(strrep("=", 70), "\n\n")

summary_results <- results |>
  group_by(design, n_participants, biomarker_moderation, biomarker_correlation,
           carryover_halflife, model_carryover) |>
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect_size, na.rm = TRUE),
    sd_effect = sd(effect_size, na.rm = TRUE),
    n = n(),
    n_errors = sum(is.na(significant)),
    .groups = "drop"
  ) |>
  mutate(
    # 95% CI for binomial proportion (Wilson score interval)
    n_success = round(power * n),
    ci_lower = (power + 1.96^2 / (2 * n) -
      1.96 * sqrt((power * (1 - power) + 1.96^2 / (4 * n)) / n)) /
      (1 + 1.96^2 / n),
    ci_upper = (power + 1.96^2 / (2 * n) +
      1.96 * sqrt((power * (1 - power) + 1.96^2 / (4 * n)) / n)) /
      (1 + 1.96^2 / n),
    ci_lower = pmax(0, ci_lower),
    ci_upper = pmin(1, ci_upper)
  )

print(summary_results, n = 50)

# Two-scenario comparison summary
cat("\n", strrep("=", 70), "\n")
cat("TWO-SCENARIO COMPARISON: Impact of Carryover Modeling\n")
cat(strrep("=", 70), "\n\n")

scenario_comparison <- summary_results |>
  filter(
    carryover_halflife > 0,
    biomarker_moderation > 0,
    design != "ol"  # OL only has model_carryover=TRUE
  ) |>
  select(design, n_participants, biomarker_moderation, carryover_halflife,
         model_carryover, power) |>
  pivot_wider(
    names_from = model_carryover,
    values_from = power,
    names_prefix = "model_carryover_"
  ) |>
  mutate(
    power_difference = model_carryover_TRUE - model_carryover_FALSE,
    pct_change = (power_difference / model_carryover_FALSE) * 100
  ) |>
  arrange(design, n_participants, carryover_halflife, biomarker_moderation)

print(scenario_comparison, n = 50)

# ============================================================================
# VISUALIZATION - Two Separate Figures
# ============================================================================
# Figure 1: Main designs (Crossover, Hybrid, OL+BDC) - WITH vs WITHOUT carryover
# Figure 2: Open-Label only (separate plot)

library(viridis)
library(patchwork)

# Prepare plot data with labels
plot_data <- summary_results |>
  filter(!is.na(power)) |>
  mutate(
    bm_mod_label = factor(
      sprintf("%.1f", biomarker_moderation),
      levels = rev(sprintf("%.1f", sort(unique(biomarker_moderation))))
    ),
    carryover_label = factor(
      paste0("t½=", carryover_halflife, " wk"),
      levels = paste0("t½=", sort(unique(carryover_halflife)), " wk")
    ),
    n_label = factor(
      paste0("N=", n_participants),
      levels = c("N=35", "N=70")
    ),
    design_label = factor(
      case_when(
        design == "crossover" ~ "Crossover",
        design == "hybrid" ~ "Hybrid/N-of-1",
        design == "ol_bdc" ~ "OL+BDC",
        design == "ol" ~ "Open-Label"
      ),
      levels = c("Crossover", "Hybrid/N-of-1", "OL+BDC", "Open-Label")
    ),
    power_pct = round(power * 100),
    ci_lower_pct = round(ci_lower * 100),
    ci_upper_pct = round(ci_upper * 100)
  )

# ==========================================================================
# FIGURE 1: Main designs (excluding Open-Label)
# ==========================================================================
# For t½=0 (no carryover), duplicate results for both WITH and WITHOUT panels
# since the models are equivalent when there's no carryover in the data
main_data_base <- plot_data |>
  filter(design != "ol")

# Create WITHOUT carryover data (includes t½=0 from model_carryover=TRUE)
without_data <- main_data_base |>
  filter(model_carryover == FALSE | (model_carryover == TRUE & carryover_halflife == 0)) |>
  mutate(
    scenario_label = factor(
      paste0("A. WITHOUT Carryover\n(N=", n_participants, ")"),
      levels = c("A. WITHOUT Carryover\n(N=35)", "A. WITHOUT Carryover\n(N=70)")
    )
  )

# Create WITH carryover data
with_data <- main_data_base |>
  filter(model_carryover == TRUE) |>
  mutate(
    scenario_label = factor(
      paste0("B. WITH Carryover\n(N=", n_participants, ")"),
      levels = c("B. WITH Carryover\n(N=35)", "B. WITH Carryover\n(N=70)")
    )
  )

# Combine
main_data <- bind_rows(without_data, with_data) |>
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = c(
        "A. WITHOUT Carryover\n(N=35)", "A. WITHOUT Carryover\n(N=70)",
        "B. WITH Carryover\n(N=35)", "B. WITH Carryover\n(N=70)"
      )
    )
  )

p_main <- ggplot(main_data, aes(x = carryover_label, y = bm_mod_label, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d", power_pct)),
    color = "black", size = 3.0, fontface = "bold",
    vjust = -0.3
  ) +
  geom_text(
    aes(label = sprintf("(%d, %d)", ci_lower_pct, ci_upper_pct)),
    color = "black", size = 1.8,
    vjust = 1.2
  ) +
  scale_fill_gradient2(
    name = "Power (%)",
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1),
    labels = function(x) sprintf("%.0f", x * 100)
  ) +
  facet_grid(
    scenario_label ~ design_label,
    labeller = label_value
  ) +
  labs(
    title = "Statistical Power for Biomarker-Treatment Interaction",
    subtitle = sprintf(
      "N=35 and N=70 | %d iterations | Clustered designs (8 timepoints)",
      n_iterations
    ),
    x = "Carryover Half-Life",
    y = "Biomarker Moderation",
    caption = sprintf("Source: simulation_clustered.R | %s", format(Sys.time(), "%Y-%m-%d"))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 9),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    plot.caption = element_text(size = 7, hjust = 1, color = "gray50"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ==========================================================================
# FIGURE 2: Open-Label only (separate figure)
# ==========================================================================
ol_data <- plot_data |>
  filter(design == "ol") |>
  mutate(
    scenario_label = factor(
      paste0("N=", n_participants),
      levels = c("N=35", "N=70")
    )
  )

p_ol <- ggplot(ol_data, aes(x = carryover_label, y = bm_mod_label, fill = power)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%d", power_pct)),
    color = "black", size = 3.5, fontface = "bold",
    vjust = -0.3
  ) +
  geom_text(
    aes(label = sprintf("(%d, %d)", ci_lower_pct, ci_upper_pct)),
    color = "black", size = 2.0,
    vjust = 1.2
  ) +
  scale_fill_gradient2(
    name = "Power (%)",
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.5, limits = c(0, 1),
    labels = function(x) sprintf("%.0f", x * 100)
  ) +
  facet_grid(
    scenario_label ~ .,
    labeller = label_value
  ) +
  labs(
    title = "Statistical Power: Open-Label Design",
    subtitle = sprintf(
      "%d iterations | 8 timepoints | WITH carryover in analysis model",
      n_iterations
    ),
    x = "Carryover Half-Life",
    y = "Biomarker Moderation",
    caption = sprintf("Source: simulation_clustered.R | %s", format(Sys.time(), "%Y-%m-%d"))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 10),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    plot.caption = element_text(size = 7, hjust = 1, color = "gray50"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Display plots
print(p_main)
print(p_ol)

# Save outputs
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Save main figure (3 designs)
ggsave(
  file.path(output_dir, sprintf("power_main_%s.pdf", timestamp)),
  p_main,
  width = 12, height = 10
)
ggsave(
  file.path(output_dir, sprintf("power_main_%s.png", timestamp)),
  p_main,
  width = 12, height = 10, dpi = 300
)

# Save Open-Label figure (separate)
ggsave(
  file.path(output_dir, sprintf("power_openlabel_%s.pdf", timestamp)),
  p_ol,
  width = 6, height = 5
)
ggsave(
  file.path(output_dir, sprintf("power_openlabel_%s.png", timestamp)),
  p_ol,
  width = 6, height = 5, dpi = 300
)

# Save results
save(
  results, summary_results, scenario_comparison,
  carryover_mode, n_iterations,
  file = file.path(output_dir, "simulation_clustered_results.RData")
)

cat("\n", strrep("=", 70), "\n")
cat("SIMULATION COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("Results saved to:", output_dir, "\n\n")
cat("Output files:\n")
cat(sprintf("  - power_main_%s.pdf/png (3 designs: Crossover, Hybrid, OL+BDC)\n", timestamp))
cat(sprintf("  - power_openlabel_%s.pdf/png (Open-Label only)\n", timestamp))
cat("  - simulation_clustered_results.RData\n")
cat("\nKey objects in RData:\n")
cat("  - results: raw iteration-level results\n")
cat("  - summary_results: aggregated power by condition\n")
cat("  - scenario_comparison: WITH vs WITHOUT carryover comparison\n")

# Close log file
sink()
