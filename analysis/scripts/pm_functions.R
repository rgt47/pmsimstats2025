#' pm_functions.R: Common functions for precision medicine trial simulation
#'
#' This script contains shared functions extracted from vig5.R for use in
#' multiple vignettes and analysis scripts.

#===========================================================================
# Function: mod_gompertz
# Description: Modified Gompertz function for modeling growth curves
#===========================================================================
mod_gompertz <- function(time, max_value, displacement, rate) {
  max_value * (1 - exp(-displacement * exp(-rate * time)))
}
#===========================================================================
# End Function: mod_gompertz
#===========================================================================



#===========================================================================
# Function: calculate_carryover
# Description: Calculate carryover effect using different decay
# models
#===========================================================================
calculate_carryover <- function(time_since_discontinuation,
                                previous_effect,
                                model = "exponential", params) {
  switch(model,
    "exponential" = {
      # Exponential decay: previous_effect *
      # (1/2)^(scale_factor * time/halflife)
      scale_factor <- if(is.null(params$scale_factor)) {
        2
      } else {
        params$scale_factor
      }
      halflife <- if(is.null(params$carryover_t1half)) {
        1
      } else {
        params$carryover_t1half
      }
      previous_effect * (1/2)^(scale_factor *
                               time_since_discontinuation /
                               halflife)
    },

    "linear" = {
      # Linear decay to zero over total_time weeks
      halflife_for_total <- if(is.null(params$carryover_t1half)) {
        1
      } else {
        params$carryover_t1half
      }
      total_time <- if(is.null(params$total_time)) {
        (2 * halflife_for_total)
      } else {
        params$total_time
      }
      decay_factor <- pmax(0, 1 - time_since_discontinuation /
                                   total_time)
      previous_effect * decay_factor
    },

    "weibull" = {
      # Weibull decay: exp(-(time/lambda)^k)
      # k=1 gives exponential, k>1 gives accelerating decay
      k <- if(is.null(params$k)) 1 else params$k
      halflife <- if(is.null(params$carryover_t1half)) {
        1
      } else {
        params$carryover_t1half
      }
      lambda <- if(is.null(params$lambda)) {
        (halflife / (log(2)^(1/k)))
      } else {
        params$lambda
      }
      previous_effect * exp(-(time_since_discontinuation /
                              lambda)^k)
    },

    stop("Unknown carryover model: ", model)
  )
}
#===========================================================================
# End Function: calculate_carryover
#===========================================================================

#===========================================================================
# Enhanced Carryover Functions
#===========================================================================

#' Calculate bio-response with biomarker×treatment interaction
#'
#' @param trial_data Trial design data with tod column
#' @param model_param Model parameters including c.bm interaction
#'   strength
#' @param resp_param Response parameters for bio_response category
#' @param component_halflives Component-specific half-lives for
#'   carryover
#' @param scale_factor Carryover scale factor
#' @return List with bio_response_means and bio_response_test
calculate_bio_response_with_interaction <- function(
    trial_data, model_param, resp_param,
    component_halflives, scale_factor) {
  # Base response follows mod_gompertz for treatment periods
  base_bio_response <- mod_gompertz(
    trial_data$tod,
    resp_param$max[resp_param$cat == "bio_response"],
    resp_param$disp[resp_param$cat == "bio_response"],
    resp_param$rate[resp_param$cat == "bio_response"]
  )

  # NOTE: Biomarker×treatment interaction emerges from correlation structure
  # in the MVN draw, following Hendrickson et al. (2020) approach.
  # No population mean shift needed - interaction is created by differential
  # correlation between biomarker and bio_response by treatment status.
  # See: WHY_POST_MVN_ADJUSTMENT_ANALYSIS.md
  bio_response_means <- base_bio_response

  # Check for zero values (for correlation matrix logic)
  bio_response_test <- base_bio_response == 0

  # Apply enhanced carryover to bio-response component
  bio_response_means <- apply_carryover_to_component(
    bio_response_means, trial_data, component_halflives,
    scale_factor, "br"
  )

  return(list(
    bio_response_means = bio_response_means,
    bio_response_test = bio_response_test,
    raw_bio_response_means = base_bio_response
  ))
}

#' Calculate component-specific carryover half-lives
#'
#' @param base_halflife Base carryover half-life from model
#'   parameters
#' @return List of component-specific half-lives
#'
#' @details Following Hendrickson et al. (2020), all components use
#' the same half-life. Only BR component actually receives carryover
#' adjustment (tv and pb are returned unchanged by
#' apply_carryover_to_component).
calculate_component_halflives <- function(base_halflife) {
  list(
    # Uniform half-life for all components (Hendrickson approach)
    br_halflife = base_halflife,
    pb_halflife = base_halflife,
    tv_halflife = base_halflife
  )
}

#' DEPRECATED: Calculate carryover-adjusted correlation parameters
#'
#' This function has been DEPRECATED and is no longer used.
#'
#' REASON FOR DEPRECATION:
#' Based on theoretical analysis (see carryover_correlation_theory.tex):
#' - Carryover affects MEANS (systematic component), not CORRELATIONS
#' - Dynamic adjustment violates positive definite constraints
#' - Contradicts standard practice in longitudinal modeling
#'
#' CURRENT APPROACH:
#' Use fixed correlation values (Hendrickson et al. 2020) that do not
#' vary with carryover parameters. Carryover is modeled in the mean
#' structure only.
#'
#' @param base_correlations Base correlation values (no longer used)
#' @param carryover_halflife Carryover half-life (no longer used)
#' @return Fixed correlation parameters (Hendrickson values)
calculate_carryover_adjusted_correlations <- function(
    base_correlations, carryover_halflife) {

  # DEPRECATED: This function now returns FIXED values
  # regardless of input parameters

  warning("calculate_carryover_adjusted_correlations() is deprecated. ",
          "Use fixed Hendrickson correlation values instead. ",
          "See correlation_structure_design.pdf for details.")

  # Return fixed Hendrickson values (no dynamic adjustment)
  return(list(
    c.tv = 0.8,    # Fixed autocorrelation (Hendrickson)
    c.pb = 0.8,    # Fixed autocorrelation (Hendrickson)
    c.br = 0.8,    # Fixed autocorrelation (Hendrickson)
    c.cf1t = 0.2,  # Fixed same-time cross-correlation (Hendrickson)
    c.cfct = 0.1   # Fixed diff-time cross-correlation (Hendrickson)
  ))

  # OLD CODE (commented out for reference):
  # carryover_strength <- carryover_halflife / (1 + carryover_halflife)
  # base_autocorr <- base_correlations$base_autocorr
  # base_cross_same <- base_correlations$base_cross_same
  # base_cross_diff <- base_correlations$base_cross_diff
  # autocorr_boost <- 0.3 * carryover_strength
  # adjusted_tv <- min(0.95, base_autocorr + autocorr_boost)
  # adjusted_pb <- min(0.95, base_autocorr + autocorr_boost)
  # adjusted_br <- min(0.95, base_autocorr + autocorr_boost * 1.2)
  # cross_time_boost <- 0.4 * carryover_strength
  # adjusted_cross_diff <- min(0.8, base_cross_diff + cross_time_boost)
  # cross_same_reduction <- 0.1 * carryover_strength
  # adjusted_cross_same <- max(0.05, base_cross_same - cross_same_reduction)
}

#' Carryover calculation for bio-response component (Hendrickson method)
#'
#' @param component_means Current component means
#' @param trial_data Trial design data
#' @param component_halflives Component-specific half-lives
#' @param scale_factor Carryover scale factor (kept for compatibility,
#'   but set to 1)
#' @param component_name Component abbreviation ("tv", "pb", "br")
#' @return Adjusted means with carryover effects
#'
#' @details Following Hendrickson et al. (2020), carryover is applied
#' ONLY to the bio-response (br) component when participants are off
#' drug. The formula is: mu[t] = base[t] + mu[t-1] * (1/2)^(tsd/t1half)
#' where tsd is time since discontinuation.
apply_carryover_to_component <- function(
    component_means, trial_data, component_halflives,
    scale_factor, component_name) {

  num_timepoints <- length(component_means)

  # ONLY apply carryover to BR component (Hendrickson approach)
  if (component_name != "br") {
    return(component_means)  # No carryover for tv or pb
  }

  if (num_timepoints > 1) {
    # Get component-specific half-life
    halflife_key <- paste0(component_name, "_halflife")
    component_halflife <- component_halflives[[halflife_key]]

    if (is.null(component_halflife) || component_halflife == 0) {
      return(component_means)  # No carryover if halflife is 0
    }

    # Bio-response: carryover when off drug and time since
    # discontinuation > 0
    carryover_indices <- which(!trial_data$on_drug &
                                trial_data$tsd > 0)

    if (length(carryover_indices) > 0) {
      for (idx in carryover_indices) {
        prev_idx <- idx - 1

        # Safety check: ensure prev_idx is valid
        if (prev_idx >= 1 &&
            prev_idx <= length(component_means) &&
            idx >= 1 && idx <= length(component_means)) {

          # Use time since discontinuation for bio-response
          time_lag <- trial_data$tsd[idx]

          # Calculate decay factor
          # NOTE: scale_factor should be 1 to match Hendrickson
          decay_factor <- (1/2)^(scale_factor * time_lag /
                                 component_halflife)

          # Apply carryover effect
          component_means[idx] <- component_means[idx] +
            component_means[prev_idx] * decay_factor
        }
      }
    }
  }

  return(component_means)
}

#===========================================================================
# Helper Functions for Data Generation
#===========================================================================

prepare_trial_data <- function(trial_design) {
  # Ensure trial_design is a data frame
  if (!is.data.frame(trial_design)) {
    stop("Trial design is not a data frame. Class: ", class(trial_design),
         ", Type: ", typeof(trial_design))
  }

  trial_data <- as_tibble(trial_design)

  # Check if t_wk exists and calculate cumulative if it does
  if ("t_wk" %in% names(trial_data)) {
    trial_data$t_wk_cumulative <- cumsum(trial_data$t_wk)
  } else {
    # If t_wk doesn't exist, calculate it or use week column directly
    if ("week" %in% names(trial_data)) {
      trial_data$t_wk <- c(trial_data$week[1], diff(trial_data$week))
      trial_data$t_wk_cumulative <- trial_data$week
    } else {
      stop("Neither 't_wk' nor 'week' column found in trial design data")
    }
  }

  trial_data <- trial_data %>% mutate(on_drug = (tod > 0))
  return(trial_data)
}

build_correlation_matrix <- function(
    labels, trial_design, model_param, num_timepoints,
    factor_types, factor_abbreviations,
    bio_response_test = NULL, bio_response_means = NULL,
    means = NULL) {
  # Build a correlation matrix
  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  # Apply correlations between factors
  for (factor_idx in seq_along(factor_types)) {
    current_factor <- factor_abbreviations[factor_idx]

    # Build autocorrelations across time - VECTORIZED
    if (num_timepoints > 1) {
      autocorrelation <- model_param[[paste("c", current_factor, sep = ".")]]

      # Create all combinations of time points at once
      point_indices <- expand.grid(
        p1 = 1:(num_timepoints-1),
        p2 = (2:num_timepoints)
      )
      # Filter valid combinations where p2 > p1
      point_indices <- point_indices[
        point_indices$p2 > point_indices$p1,
      ]

      # Create name vectors for efficient indexing
      name1 <- paste(trial_design$timepoint_name[
                       point_indices$p1],
                     current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name[
                       point_indices$p2],
                     current_factor, sep = ".")

      # Set all correlations at once using matrix indexing
      for (idx in 1:nrow(point_indices)) {
        correlations[name1[idx], name2[idx]] <- autocorrelation
        correlations[name2[idx], name1[idx]] <- autocorrelation
      }
    }

    # Build autocorrelations across factors - VECTORIZED
    for (factor2_idx in setdiff(seq_along(factor_types), factor_idx)) {
      other_factor <- factor_abbreviations[factor2_idx]

      # Same timepoint cross-factor correlations
      name1 <- paste(trial_design$timepoint_name,
                     current_factor, sep = ".")
      name2 <- paste(trial_design$timepoint_name,
                     other_factor, sep = ".")

      # Set all same-timepoint correlations at once
      for (idx in 1:length(name1)) {
        correlations[name1[idx], name2[idx]] <- model_param$c.cf1t
        correlations[name2[idx], name1[idx]] <- model_param$c.cf1t
      }

      # Different timepoint cross-factor correlations
      if (num_timepoints > 1) {
        # Create all combinations of time points at once
        point_indices <- expand.grid(
          p1 = 1:(num_timepoints-1),
          p2 = (2:num_timepoints)
        )
        # Filter valid combinations where p2 > p1
        point_indices <- point_indices[
          point_indices$p2 > point_indices$p1,
        ]

        # Create name vectors for efficient indexing
        name1 <- paste(trial_design$timepoint_name[
                         point_indices$p1],
                       current_factor, sep = ".")
        name2 <- paste(trial_design$timepoint_name[
                         point_indices$p2],
                       other_factor, sep = ".")

        # Set all correlations at once
        for (idx in 1:nrow(point_indices)) {
          correlations[name1[idx], name2[idx]] <- model_param$c.cfct
          correlations[name2[idx], name1[idx]] <- model_param$c.cfct
        }
      }
    }

    # Special handling for biomarker correlation
    if (current_factor == "br") {
      for (timepoint_idx in 1:num_timepoints) {
        name1 <- paste(trial_design$timepoint_name[
                         timepoint_idx], "br", sep = ".")

        if (timepoint_idx > 1) {
          # Use advanced biomarker correlation logic if data
          # is provided
          if (!is.null(bio_response_test) &&
              !is.null(bio_response_means) &&
              !is.null(means)) {
            name0 <- paste(
              trial_design$timepoint_name[timepoint_idx - 1],
              "br", sep = "."
            )
            mean_value1 <- means[which(name1 == labels)]
            mean_value0 <- means[which(name0 == labels)]

            # Handle special cases with careful checks
            # Calculate correlation with mean-value scaling
            scaled_correlation <- ifelse(
              bio_response_test[timepoint_idx],
              ifelse(
                bio_response_means[timepoint_idx] == 0 ||
                  abs(mean_value0) < 1e-10,
                0,
                (mean_value1 / max(mean_value0, 1e-10)) *
                  model_param$c.bm
              ),
              model_param$c.bm
            )

            # CRITICAL: Clamp correlation to valid range [-0.99, 0.99]
            # Prevent correlations from exceeding 1 or going below -1
            # Leave small margin (0.99) to ensure positive definiteness
            scaled_correlation <- pmax(-0.99, pmin(0.99, scaled_correlation))

            correlations["biomarker", name1] <-
              correlations[name1, "biomarker"] <-
              scaled_correlation
          } else {
            # Simple biomarker correlation
            correlations["biomarker", name1] <-
              correlations[name1, "biomarker"] <-
              model_param$c.bm
          }
        }
      }
    }
  }

  return(correlations)
}

#===========================================================================
# Function: generate_data
# Description: Primary data generation function for trial
# simulation
#===========================================================================
generate_data <- function(
    model_param, resp_param, baseline_param, trial_design,
    empirical, make_positive_definite, seed = NA,
    scale_factor = 2, verbose = FALSE, track_pd_stats = TRUE,
    cached_sigma = NULL) {
  # Initialize variables for tracking sigma matrix statistics
  sigma_count <- 0
  non_positive_definite_count <- 0

  # Use cached sigma if provided, otherwise build from scratch
  if (!is.null(cached_sigma)) {
    # Use pre-computed sigma matrix
    sigma <- cached_sigma$sigma
    labels <- cached_sigma$labels
    standard_deviations <- cached_sigma$standard_deviations
    correlations <- cached_sigma$correlations

    # Extract needed variables for data generation
    trial_data <- prepare_trial_data(trial_design)
    num_timepoints <- dim(trial_design)[1]
    factor_types <- c("time_variant", "pharm_biomarker", "bio_response")
    factor_abbreviations <- c("tv", "pb", "br")

    # Skip sigma construction, means will be calculated below

  } else {
    # Original sigma construction path
    # I. Turn the trial design information into something easier to use
    trial_data <- prepare_trial_data(trial_design)
    num_timepoints <- dim(trial_design)[1]

    # II. Set up variables to track - baseline parameters for
    # each participant ("biomarker","baseline"), and the three
    # modeled factors for each stage of the trial.

    # Set up the variable names
    factor_types <- c("time_variant", "pharm_biomarker", "bio_response")
    factor_abbreviations <- c("tv", "pb", "br")

    labels <- c(
      c("biomarker", "baseline"),
      paste(trial_design$timepoint_name, factor_abbreviations[1], sep = "."),
      paste(trial_design$timepoint_name, factor_abbreviations[2], sep = "."),
      paste(trial_design$timepoint_name, factor_abbreviations[3], sep = ".")
    )

    # Set up vectors with the standard deviations and means
    standard_deviations <- c(
      baseline_param$sd[baseline_param$cat == "biomarker"],
      baseline_param$sd[baseline_param$cat == "baseline"]
    )
    standard_deviations <- c(
      standard_deviations,
      rep(resp_param$sd[resp_param$cat == "time_variant"], num_timepoints)
    )
    standard_deviations <- c(
      standard_deviations,
      rep(resp_param$sd[resp_param$cat == "pharm_biomarker"], num_timepoints) * trial_design$e
    )
    standard_deviations <- c(
      standard_deviations,
      rep(resp_param$sd[resp_param$cat == "bio_response"], num_timepoints)
    )

    # Note: bio_response variables will be calculated below, so use NULL for now
    # The correlation matrix will be rebuilt after means calculation
    correlations <- NULL

    # Sigma matrix will be built after means calculation

  } # End of else block for sigma construction

  # Calculate component-specific half-lives for enhanced carryover
  component_halflives <- calculate_component_halflives(model_param$carryover_t1half)

  # Calculate means with enhanced carryover (works for both cached and non-cached)
  means <- c(
    baseline_param$m[baseline_param$cat == "biomarker"],
    baseline_param$m[baseline_param$cat == "baseline"]
  )

  for (factor_idx in seq_along(factor_types)) {
    current_factor <- factor_types[factor_idx]
    factor_abbrev <- factor_abbreviations[factor_idx]

    if (factor_abbrev == "tv") {
      time_variant_means <- mod_gompertz(
        trial_data$t_wk_cumulative,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      )

      # Apply enhanced carryover to time-variant component
      time_variant_means <- apply_carryover_to_component(
        time_variant_means, trial_data, component_halflives, scale_factor, "tv"
      )

      means <- c(means, time_variant_means)
    }

    if (factor_abbrev == "pb") {
      pharm_biomarker_means <- mod_gompertz(
        trial_data$tpb,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      ) * trial_design$e

      # Apply enhanced carryover to pharmacological-biomarker component
      pharm_biomarker_means <- apply_carryover_to_component(
        pharm_biomarker_means, trial_data, component_halflives, scale_factor, "pb"
      )

      means <- c(means, pharm_biomarker_means)
    }

    if (factor_abbrev == "br") {
      # Calculate bio-response with biomarker×treatment interaction using helper function
      br_result <- calculate_bio_response_with_interaction(
        trial_data, model_param, resp_param, component_halflives, scale_factor
      )
      
      bio_response_means <- br_result$bio_response_means
      bio_response_test <- br_result$bio_response_test

      # Find indices for bio_response (br) columns dynamically
      br_indices <- grep("\\.br$", labels)

      # Set names for debugging - using dynamic indices
      if (length(br_indices) > 0) {
        names(bio_response_test) <- labels[br_indices]
      }

      means <- c(means, bio_response_means)
    }
  }

  # Build correlation matrix for non-cached path (if needed)
  if (is.null(cached_sigma)) {
    correlations <- build_correlation_matrix(
      labels, trial_design, model_param, num_timepoints,
      factor_types, factor_abbreviations,
      bio_response_test, bio_response_means, means
    )

    # Track statistics about sigma matrix if requested
    if (track_pd_stats) {
      sigma_count <- sigma_count + 1
    }

    # Fast path for positive definiteness handling
    is_positive_definite <- TRUE  # Assume matrix is positive definite until proven otherwise
    need_pd_check <- make_positive_definite || track_pd_stats

    # Turn correlation matrix into covariance matrix using efficient outer product
    sigma <- outer(standard_deviations, standard_deviations) * correlations

    # Check/fix positive definiteness if required
    if (need_pd_check) {
      is_positive_definite <- corpcor::is.positive.definite(sigma)

      # Update statistics if tracking
      if (track_pd_stats && !is_positive_definite) {
        non_positive_definite_count <- non_positive_definite_count + 1
      }
    }

    # Fix if turned on and not positive definite
    if (make_positive_definite && !is_positive_definite) {
      # For more robust conversion to positive definite
      sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
    }
  }

  # Set the seed if provided for reproducibility
  if (!is.na(seed)) {
    set.seed(seed)
  }

  # Generate multivariate normal data
  participant_data <- MASS::mvrnorm(
    n = model_param$N,
    mu = means,
    Sigma = sigma,
    empirical = empirical
  )
  participant_data <- as_tibble(participant_data)
  colnames(participant_data) <- labels

  # NOTE: No post-MVN biomarker adjustment needed.
  # Biomarker×treatment interaction emerges naturally from differential
  # correlation structure in the MVN draw (Hendrickson et al. 2020 approach).
  # The correlation matrix specifies:
  #   - Cor(biomarker, bio_response) = c.bm when on treatment
  #   - Cor(biomarker, bio_response) = 0 (or scaled) when off treatment
  # This differential correlation IS the biomarker×treatment interaction.
  # See: WHY_POST_MVN_ADJUSTMENT_ANALYSIS.md for detailed explanation.

  # Process participant data using helper function
  participant_data <- process_participant_data(
    participant_data, trial_design$timepoint_name,
    factor_abbreviations, model_param$N
  )

  # Attach matrix tracking statistics as attributes if requested
  if (track_pd_stats) {
    non_positive_definite_rate <- non_positive_definite_count / sigma_count
    attr(participant_data, "sigma_count") <- sigma_count
    attr(participant_data, "non_positive_definite_count") <- non_positive_definite_count
    attr(participant_data, "non_positive_definite_rate") <- non_positive_definite_rate
  }

  return(participant_data)
}
#===========================================================================
# End Function: generate_data
#===========================================================================

#===========================================================================
# Helper function for building and caching sigma matrices
#===========================================================================
build_sigma_matrix <- function(model_param, resp_param, baseline_param, trial_design,
                              factor_types, factor_abbreviations, verbose = FALSE) {
  # This function builds just the sigma matrix without generating data

  trial_data <- prepare_trial_data(trial_design)
  num_timepoints <- dim(trial_design)[1]

  # Set up variable names (same as generate_data)
  labels <- c(
    c("biomarker", "baseline"),
    paste(trial_design$timepoint_name, factor_abbreviations[1], sep = "."),
    paste(trial_design$timepoint_name, factor_abbreviations[2], sep = "."),
    paste(trial_design$timepoint_name, factor_abbreviations[3], sep = ".")
  )

  # Set up standard deviations (same as generate_data)
  standard_deviations <- c(
    baseline_param$sd[baseline_param$cat == "biomarker"],
    baseline_param$sd[baseline_param$cat == "baseline"]
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "time_variant"], num_timepoints)
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "pharm_biomarker"], num_timepoints) * trial_design$e
  )
  standard_deviations <- c(
    standard_deviations,
    rep(resp_param$sd[resp_param$cat == "bio_response"], num_timepoints)
  )

  # Calculate component-specific half-lives for enhanced carryover (same as generate_data)
  component_halflives <- calculate_component_halflives(model_param$carryover_t1half)

  # Calculate means for biomarker correlation logic with enhanced carryover
  means <- c(
    baseline_param$m[baseline_param$cat == "biomarker"],
    baseline_param$m[baseline_param$cat == "baseline"]
  )

  # Calculate means for each factor with enhanced carryover
  for (factor_idx in seq_along(factor_types)) {
    current_factor <- factor_types[factor_idx]
    factor_abbrev <- factor_abbreviations[factor_idx]

    if (factor_abbrev == "tv") {
      time_variant_means <- mod_gompertz(
        trial_data$t_wk_cumulative,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      )

      # Apply enhanced carryover to time-variant component
      time_variant_means <- apply_carryover_to_component(
        time_variant_means, trial_data, component_halflives, 2, "tv"
      )

      means <- c(means, time_variant_means)
    }

    if (factor_abbrev == "pb") {
      pharm_biomarker_means <- mod_gompertz(
        trial_data$tpb,
        resp_param$max[resp_param$cat == current_factor],
        resp_param$disp[resp_param$cat == current_factor],
        resp_param$rate[resp_param$cat == current_factor]
      ) * trial_design$e

      # Apply enhanced carryover to pharmacological-biomarker component
      pharm_biomarker_means <- apply_carryover_to_component(
        pharm_biomarker_means, trial_data, component_halflives, 2, "pb"
      )

      means <- c(means, pharm_biomarker_means)
    }

    if (factor_abbrev == "br") {
      # Calculate bio-response with biomarker×treatment interaction using helper function
      br_result <- calculate_bio_response_with_interaction(
        trial_data, model_param, resp_param, component_halflives, 2
      )
      
      bio_response_means <- br_result$bio_response_means
      bio_response_test <- br_result$bio_response_test

      means <- c(means, bio_response_means)
    }
  }

  # Build correlation matrix with full biomarker logic
  correlations <- build_correlation_matrix(labels, trial_design, model_param, num_timepoints,
                                          factor_types, factor_abbreviations,
                                          bio_response_test, bio_response_means, means)

  # Convert to covariance matrix
  sigma <- outer(standard_deviations, standard_deviations) * correlations

  # Check positive definiteness - REJECT non-PD matrices
  is_pd <- corpcor::is.positive.definite(sigma)

  if (verbose) {
    cat("Sigma matrix positive definite:", is_pd, "\n")
  }

  if (!is_pd) {
    if (verbose) {
      cat("REJECTED: Non-positive definite matrix detected.\n")
      cat("Choose better correlation parameters to fix this issue.\n")

      # Compute eigenvalues for diagnostics
      eigenvalues <- eigen(sigma, only.values = TRUE)$values
      cat("Eigenvalue range: [", min(eigenvalues), ", ",
          max(eigenvalues), "]\n")
      cat("Number of negative eigenvalues:",
          sum(eigenvalues < 0), "\n")
    }
    return(NULL)  # Return NULL to signal invalid parameter combination
  }

  return(list(
    sigma = sigma,
    labels = labels,
    standard_deviations = standard_deviations,
    correlations = correlations
  ))
}

#===========================================================================
# Function: validate_correlation_structure
# Description: Validate that correlation parameters produce a positive
# definite covariance matrix with good conditioning
#
# Returns: TRUE if valid, FALSE if non-positive definite
# Side effects: Prints diagnostic information about eigenvalues and
#               condition number
#===========================================================================
validate_correlation_structure <- function(model_params,
                                          resp_param,
                                          baseline_param,
                                          trial_design) {
  # Build sigma matrix
  sigma_result <- build_sigma_matrix(
    model_params, resp_param, baseline_param,
    trial_design,
    factor_types = c("time_variant", "pharm_biomarker",
                     "bio_response"),
    factor_abbreviations = c("tv", "pb", "br"),
    verbose = TRUE
  )

  if (is.null(sigma_result)) {
    cat("FAILED: Non-positive definite\n")
    return(FALSE)
  }

  # Check condition number
  sigma <- sigma_result$sigma
  eigenvalues <- eigen(sigma, only.values = TRUE)$values
  condition_number <- max(eigenvalues) / min(eigenvalues)

  cat("Eigenvalue range: [", min(eigenvalues), ", ",
      max(eigenvalues), "]\n")
  cat("Condition number: ", condition_number, "\n")

  # Well-conditioned if condition number < 100
  if (condition_number > 100) {
    warning("Matrix is ill-conditioned (condition number = ",
            condition_number, ")")
  }

  return(TRUE)
}

create_sigma_cache_key <- function(design_name, params) {
  # Create unique key for each parameter/design combination
  paste(design_name, 
        params$n_participants,
        params$biomarker_correlation,
        params$carryover_t1half,
        sep = "_")
}

#===========================================================================
# Helper function for processing participant data
#===========================================================================
process_participant_data <- function(participant_data, timepoint_names, factor_abbreviations, N) {
  # Add participant ID column
  participant_data <- participant_data %>%
    mutate(participant_id = 1:N)

  # Pre-allocate all result columns to avoid growing the data frame
  new_cols <- c(paste0("D_", timepoint_names), timepoint_names)
  zeros_df <- as_tibble(matrix(0, nrow = N, ncol = length(new_cols)))
  colnames(zeros_df) <- new_cols
  participant_data <- bind_cols(participant_data, zeros_df)

  # Calculate deltas (sums of factors) for all timepoints
  for (timepoint_idx in 1:length(timepoint_names)) {
    timepoint_name <- timepoint_names[timepoint_idx]
    delta_col <- paste0("D_", timepoint_name)
    components <- paste(timepoint_name, factor_abbreviations, sep = ".")

    # Calculate delta for all participants at once
    participant_data <- participant_data %>%
      mutate(!!delta_col := rowSums(select(., all_of(components))))
  }

  # Calculate timepoint scores from baseline and factors
  for (timepoint_idx in 1:length(timepoint_names)) {
    timepoint_name <- timepoint_names[timepoint_idx]
    components <- paste(timepoint_name, factor_abbreviations, sep = ".")

    # Calculate timepoint value from baseline and factors
    participant_data <- participant_data %>%
      mutate(!!timepoint_name := baseline - rowSums(select(., all_of(components))))
  }

  return(participant_data)
}

#===========================================================================
# Helper function for data preparation
#===========================================================================
prepare_analysis_data <- function(trial_design_set, data) {
  num_groups <- length(trial_design_set)

  # Pre-process all trial designs at once
  trial_designs_with_baseline <- map(1:num_groups, function(group_idx) {
    trial_design <- trial_design_set[[group_idx]]
    bind_rows(
      tibble(
        timepoint_names = "baseline",
        t_wk = 0,
        e = 0,
        tod = 0,
        tsd = 0,
        tpb = 0
      ),
      as_tibble(trial_design)
    ) %>%
    mutate(
      t = cumsum(t_wk),
      group = group_idx
    )
  })

  # Add path column if it doesn't exist
  if (!"path" %in% names(data)) {
    data$path <- 1
    message("Note: 'path' column not found in data. Adding path=1 to all rows.")
  }

  # Process data for all groups
  all_data <- tibble()
  last_participant_id <- 0

  for (group_idx in 1:num_groups) {
    trial_design <- trial_designs_with_baseline[[group_idx]]
    group_data <- data %>%
      filter(path == group_idx) %>%
      as_tibble()

    timepoint_names <- c("baseline", unique(trial_design$timepoint_names))
    valid_timepoint_names <- timepoint_names[timepoint_names %in% names(group_data)]

    if (length(valid_timepoint_names) < length(timepoint_names)) {
      missing_timepoints <- setdiff(timepoint_names, valid_timepoint_names)
      message("Warning: Missing timepoint columns: ", paste(missing_timepoints, collapse=", "))
    }


    group_processed <- group_data %>%
      select(participant_id, biomarker, all_of(valid_timepoint_names)) %>%
      mutate(participant_id = participant_id + last_participant_id)

    last_participant_id <- max(group_processed$participant_id)

    group_long <- group_processed %>%
      pivot_longer(
        cols = all_of(valid_timepoint_names),
        names_to = "timepoint_names",
        values_to = "symptoms",
        values_drop_na = FALSE
      )

    group_merged <- group_long %>%
      left_join(trial_design, by = "timepoint_names")

    all_data <- bind_rows(all_data, group_merged)
  }

  return(all_data)
}

#===========================================================================
# Function: lme_analysis
# Description: Linear mixed effects analysis of trial data
#===========================================================================
lme_analysis <- function(trial_design_set, data, options = list()) {
  # Set default options if not provided
  default_options <- list(
    use_expectancy = TRUE,
    random_slope = FALSE,
    full_output = FALSE,
    simple_carryover = FALSE,
    carryover_halflife = 0,
    carryover_scale_factor = 1
  )

  # Merge provided options with defaults
  options <- modifyList(default_options, options)

  # Validate incompatible options
  if ((options$carryover_halflife > 0) && (options$simple_carryover == TRUE)) {
    stop("Cannot use both simple_carryover and carryover_halflife models simultaneously")
  }

  # Prepare analysis data using helper function
  all_data <- prepare_analysis_data(trial_design_set, data)

  # Calculate derived variables for model
  data_for_model <- all_data %>%
    mutate(
      drug_binary = ifelse(tod > 0, 1, 0)
    )

  # Apply carryover effects if specified - optimized calculation
  if (options$carryover_halflife > 0) {
    # Pre-calculate the half-life factor for efficiency
    half_life_factor <- options$carryover_scale_factor / options$carryover_halflife

    # Calculate carryover effects using tidyverse approach
    data_for_model <- data_for_model %>%
      group_by(participant_id) %>%
      mutate(
        # Calculate exponential falloff for carryover
        # Handle NA values safely
        drug_binary = case_when(
          is.na(tod) ~ NA_real_,
          tod > 0 ~ 1,
          is.na(tsd) ~ 0,
          tsd <= 0 ~ 0,
          TRUE ~ dplyr::lag(drug_binary, default = 0) * (1/2)^(half_life_factor * tsd)
        )
      ) %>%
      ungroup()
  } else if (options$simple_carryover) {
    # Simple carryover is already efficient
    data_for_model <- data_for_model %>%
      mutate(tsd_effect = tsd)
  }

  # Calculate group average at each timepoint after baseline
  data_for_model <- data_for_model %>%
    group_by(timepoint_names) %>%
    mutate(mean_symptoms = mean(symptoms, na.rm = TRUE)) %>%
    ungroup()

  # Ensure consistent column order
  data_for_model <- data_for_model %>%
    select(participant_id, biomarker, timepoint_names, t, symptoms, drug_binary, everything())

  # Build formula based on options
  formula_str <- "symptoms ~ biomarker"

  # Add drug effect (drug_binary)
  formula_str <- paste(formula_str, "+ drug_binary")

  # Add interaction term (biomarker * drug effect)
  formula_str <- paste(formula_str, "+ biomarker:drug_binary")

  # Add time effect if using it
  if (options$use_expectancy) {
    formula_str <- paste(formula_str, "+ t")
  }

  # Add simple carryover effect if using it
  if (options$simple_carryover) {
    formula_str <- paste(formula_str, "+ tsd_effect")
  }

  # Add random effects
  if (options$random_slope) {
    formula_str <- paste(formula_str, "+ (1 + drug_binary|participant_id)")
  } else {
    formula_str <- paste(formula_str, "+ (1|participant_id)")
  }

  # Convert to formula
  formula <- as.formula(formula_str)

  # Fit the model and capture any warnings
  fit_messages <- character()
  withCallingHandlers(
    model <- lmer(formula, data = data_for_model),
    warning = function(w) {
      fit_messages <<- c(fit_messages, w$message)
    }
  )

  # Extract relevant summary information
  model_summary <- summary(model)

  # Find the relevant coefficient (biomarker:drug_binary interaction)
  beta_row <- which(rownames(model_summary$coefficients) == "biomarker:drug_binary")

  # Create output
  output <- tibble(
    beta = model_summary$coefficients[beta_row, 1],
    beta_se = model_summary$coefficients[beta_row, 2],
    p_value = model_summary$coefficients[beta_row, 5],
    is_singular = is_singular(model),
    warning = paste(fit_messages, collapse = "; ")
  )

  # Return appropriate output based on options
  if (options$full_output) {
    return(list(
      formula = formula,
      model = model,
      data = data_for_model,
      summary = output
    ))
  } else {
    return(output)
  }
}
#===========================================================================
# End Function: lme_analysis
#===========================================================================
