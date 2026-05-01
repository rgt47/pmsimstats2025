# simulation_clustered_hendrickson.R
# Hendrickson-aligned N-of-1 trial simulation
#
# This script replicates the Hendrickson et al. (2020) DGP exactly,
# then compares two analysis scenarios:
#   Scenario 1: Model carryover in analysis (pmsimstats2025 contribution)
#   Scenario 2: Ignore carryover in analysis (Hendrickson replication)
#
# Key differences from simulation_clustered.R:
# 1. Uses modgompertz_orig (offset-rescaled, f(0)=0)
# 2. Per-path sigma construction via build_path_sigma
# 3. Interval-based tod/tsd (not counting-based)
# 4. Ron Thomas BM-BR correlation scaling
# 5. DGP scalefactor=2, analysis scalefactor=1 (matching orig)
# 6. Hendrickson's extracted parameter values
# 7. Orig's lmer formula: Sx ~ bm + t + Dbc + bm*Dbc + (1|ptID)

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(MASS)
  library(patchwork)
  library(conflicted)
  library(future)
  library(furrr)
  library(viridis)
})

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)

source('analysis/scripts/pm_functions.R')

# ============================================================================
# PARAMETERS
# ============================================================================

n_iterations <- 1000

exploratory_mode <- as.logical(Sys.getenv('EXPLORATORY_MODE', 'FALSE'))
env_iterations <- Sys.getenv('N_ITERATIONS', '')
if (nchar(env_iterations) > 0) {
  n_iterations <- as.integer(env_iterations)
  cat('*** CUSTOM ITERATIONS:', n_iterations, '***\n\n')
} else if (exploratory_mode) {
  n_iterations <- 20
  cat('*** EXPLORATORY MODE:', n_iterations, 'iterations ***\n\n')
}

n_cores <- max(1, parallel::detectCores() - 1)
use_parallel <- TRUE

dgp_scalefactor <- 2
analysis_scalefactor <- 2

# ============================================================================
# TRIAL DESIGN DEFINITIONS (matching Hendrickson vignette exactly)
# ============================================================================

trial_designs <- list(
  `N-of-1` = list(
    timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    timeptnames = c('OL1', 'OL2', 'BD1', 'BD2',
                    'BD3', 'BD4', 'COd', 'COp'),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    paths = list(
      pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
      pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
      pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
      pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
    )
  ),
  OL = list(
    timepoints = cumsum(rep(2.5, 8)),
    timeptnames = paste0('OL', 1:8),
    expectancies = rep(1, 8),
    paths = list(
      pathA = rep(1, 8)
    )
  ),
  `OL+BDC` = list(
    timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
    timeptnames = c('OL1', 'OL2', 'OL3', 'OL4',
                    'BD1', 'BD2', 'BD3', 'BD4'),
    expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
    paths = list(
      pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
      pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
    )
  ),
  CO = list(
    timepoints = cumsum(rep(2.5, 8)),
    timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    expectancies = rep(0.5, 8),
    paths = list(
      pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
      pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
    )
  )
)

# ============================================================================
# PARAMETER GRID (matching Hendrickson's coremodelparams)
# ============================================================================

param_grid <- expand_grid(
  design = names(trial_designs),
  N = c(35, 70),
  c.bm = c(0, 0.3, 0.6),
  carryover_t1half = c(0, 0.1, 0.2, 1.0, 2.0)
) |>
  mutate(condition_id = row_number())

cat(sprintf(
  'Parameter grid: %d conditions x %d iterations = %d total fits\n',
  nrow(param_grid), n_iterations,
  nrow(param_grid) * n_iterations * 2
))

# ============================================================================
# PRE-BUILD AND CACHE SIGMA MATRICES
# ============================================================================

cat('\nBuilding sigma cache...\n')

sigma_cache <- list()

for (i in 1:nrow(param_grid)) {
  row <- param_grid[i, ]
  td <- trial_designs[[row$design]]

  for (path_name in names(td$paths)) {
    cache_key <- paste(row$design, path_name,
      row$c.bm, row$carryover_t1half, sep = '_')

    if (!cache_key %in% names(sigma_cache)) {
      sigma_cache[[cache_key]] <- build_path_sigma(
        timepoints = td$timepoints,
        ondrug = td$paths[[path_name]],
        expectancies = td$expectancies,
        c.bm = row$c.bm,
        carryover_t1half = row$carryover_t1half,
        scalefactor = dgp_scalefactor
      )
    }
  }
}

cat(sprintf('  Cached %d unique sigma matrices\n', length(sigma_cache)))

n_non_pd <- sum(
  !sapply(sigma_cache, function(x) x$diagnostics$was_pd)
)
if (n_non_pd > 0) {
  cat(sprintf(
    '  WARNING: %d matrices required PD enforcement\n', n_non_pd
  ))
}
cat('\n')

# ============================================================================
# SIMULATION FUNCTION
# ============================================================================

run_single_iteration <- function(iter, params, sigma_cache,
                                 trial_designs,
                                 analysis_scalefactor) {
  set.seed(iter * 1000 + params$condition_id)

  td <- trial_designs[[params$design]]
  n_paths <- length(td$paths)
  participants_per_path <- params$N %/% n_paths
  remainder <- params$N %% n_paths

  all_analysis_data <- list()
  global_pid <- 0

  for (path_idx in seq_along(td$paths)) {
    path_name <- names(td$paths)[path_idx]
    ondrug <- td$paths[[path_idx]]

    cache_key <- paste(params$design, path_name,
      params$c.bm, params$carryover_t1half, sep = '_')
    ps <- sigma_cache[[cache_key]]

    n_this_path <- participants_per_path +
      (if (path_idx <= remainder) 1 else 0)

    dat <- mvrnorm(
      n = n_this_path,
      mu = ps$means,
      Sigma = ps$sigma
    )
    if (n_this_path == 1) dat <- matrix(dat, nrow = 1)
    colnames(dat) <- ps$labels

    nP <- ps$indices$nP
    tpnames <- ps$path_info$timeptnames
    t_cumul <- ps$path_info$t_cumulative
    tsd_vec <- ps$path_info$tsd

    for (pid in 1:n_this_path) {
      global_pid <- global_pid + 1

      for (tp in 1:nP) {
        tv_col <- paste(tpnames[tp], 'tv', sep = '.')
        pb_col <- paste(tpnames[tp], 'pb', sep = '.')
        br_col <- paste(tpnames[tp], 'br', sep = '.')

        Sx <- dat[pid, 'BL'] - (
          dat[pid, tv_col] +
          dat[pid, pb_col] +
          dat[pid, br_col]
        )

        Db <- as.numeric(ondrug[tp] == 1)

        Dbc <- if (ondrug[tp] == 1) {
          1
        } else if (params$carryover_t1half == 0 ||
                   tsd_vec[tp] == 0) {
          0
        } else {
          (1/2)^(analysis_scalefactor * tsd_vec[tp] /
                 params$carryover_t1half)
        }

        all_analysis_data[[length(all_analysis_data) + 1]] <- list(
          ptID = global_pid,
          bm = dat[pid, 'bm'],
          Sx = Sx,
          t = t_cumul[tp],
          Db = Db,
          Dbc = Dbc,
          tsd = tsd_vec[tp],
          path = path_idx
        )
      }
    }
  }

  merged <- bind_rows(lapply(all_analysis_data, as_tibble))

  has_var_in_db <- length(unique(merged$Db)) > 1

  results <- list()

  # Scenario 2: Hendrickson replication (no carryover in analysis)
  fit_s2 <- tryCatch({
    if (has_var_in_db) {
      lmer(Sx ~ bm + t + Db + bm * Db + (1 | ptID),
        data = merged)
    } else {
      merged_ondrug <- merged |>
        group_by(ptID) |>
        filter(any(Db == 1)) |>
        ungroup()
      lmer(Sx ~ bm + t + bm * t + (1 | ptID),
        data = merged_ondrug)
    }
  }, error = function(e) NULL)

  if (!is.null(fit_s2)) {
    coefs_s2 <- summary(fit_s2)$coefficients
    target_s2 <- if (has_var_in_db) 'bm:Db' else 'bm:t'
    if (target_s2 %in% rownames(coefs_s2)) {
      t_val <- coefs_s2[target_s2, 't value']
      p_val <- 2 * pnorm(-abs(t_val))
      results$s2 <- tibble(
        iteration = iter,
        scenario = 'S2_no_carryover',
        design = params$design,
        N = params$N,
        c.bm = params$c.bm,
        carryover_t1half = params$carryover_t1half,
        effect = coefs_s2[target_s2, 'Estimate'],
        se = coefs_s2[target_s2, 'Std. Error'],
        t_value = t_val,
        p_value = p_val,
        significant = p_val < 0.05
      )
    }
  }

  # Scenario 1: With carryover modeling (pmsimstats2025 contribution)
  if (has_var_in_db && params$carryover_t1half > 0) {
    fit_s1 <- tryCatch({
      lmer(Sx ~ bm + t + Dbc + bm * Dbc + (1 | ptID),
        data = merged)
    }, error = function(e) NULL)

    if (!is.null(fit_s1)) {
      coefs_s1 <- summary(fit_s1)$coefficients
      if ('bm:Dbc' %in% rownames(coefs_s1)) {
        t_val <- coefs_s1['bm:Dbc', 't value']
        p_val <- 2 * pnorm(-abs(t_val))
        results$s1 <- tibble(
          iteration = iter,
          scenario = 'S1_with_carryover',
          design = params$design,
          N = params$N,
          c.bm = params$c.bm,
          carryover_t1half = params$carryover_t1half,
          effect = coefs_s1['bm:Dbc', 'Estimate'],
          se = coefs_s1['bm:Dbc', 'Std. Error'],
          t_value = t_val,
          p_value = p_val,
          significant = p_val < 0.05
        )
      }
    }
  }

  if (length(results) == 0) {
    tibble(
      iteration = iter,
      scenario = NA_character_,
      design = params$design,
      N = params$N,
      c.bm = params$c.bm,
      carryover_t1half = params$carryover_t1half,
      effect = NA_real_,
      se = NA_real_,
      t_value = NA_real_,
      p_value = NA_real_,
      significant = NA
    )
  } else {
    bind_rows(results)
  }
}

# ============================================================================
# RUN SIMULATION
# ============================================================================

output_dir <- 'analysis/output'
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

log_file <- file.path(output_dir,
  'simulation_clustered_hendrickson_log.txt')
sink(log_file, split = TRUE)

cat('Running HENDRICKSON-ALIGNED simulation...\n')
cat('Designs:', paste(names(trial_designs), collapse = ', '), '\n')
cat(sprintf('DGP scalefactor: %d\n', dgp_scalefactor))
cat(sprintf('Analysis scalefactor: %d\n', analysis_scalefactor))
cat(sprintf('Iterations: %d\n', n_iterations))
cat(sprintf('Conditions: %d\n', nrow(param_grid)))
cat(sprintf('Parallel: %s (%d cores)\n\n',
  if (use_parallel) 'YES' else 'NO', n_cores))

if (use_parallel && n_cores > 1) {
  plan(multisession, workers = n_cores)
  cat(sprintf('Parallel backend: %d workers\n\n', n_cores))
} else {
  plan(sequential)
  cat('Sequential mode\n\n')
}

start_time <- Sys.time()
results_list <- list()

for (i in 1:nrow(param_grid)) {
  params <- as.list(param_grid[i, ])
  cat(sprintf(
    'Condition %d/%d: %s N=%d c.bm=%.1f t1half=%.1f\n',
    i, nrow(param_grid), params$design, params$N,
    params$c.bm, params$carryover_t1half
  ))

  local_params <- params
  local_cache <- sigma_cache
  local_designs <- trial_designs
  local_sf <- analysis_scalefactor

  condition_results <- future_map_dfr(
    1:n_iterations,
    function(iter) {
      run_single_iteration(
        iter, local_params, local_cache,
        local_designs, local_sf
      )
    },
    .options = furrr_options(
      seed = TRUE,
      globals = list(
        run_single_iteration = run_single_iteration,
        local_params = local_params,
        local_cache = local_cache,
        local_designs = local_designs,
        local_sf = local_sf,
        build_path_sigma = build_path_sigma,
        modgompertz_orig = modgompertz_orig,
        compute_tod_orig = compute_tod_orig,
        compute_tsd_orig = compute_tsd_orig,
        compute_tpb_orig = compute_tpb_orig,
        build_bm_br_correlations = build_bm_br_correlations,
        hendrickson_resp_params = hendrickson_resp_params,
        hendrickson_bl_params = hendrickson_bl_params,
        hendrickson_corr_params = hendrickson_corr_params
      ),
      packages = c('tidyverse', 'lme4', 'MASS', 'corpcor')
    )
  )

  results_list[[i]] <- condition_results

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = 'secs'))
  completed <- i
  rate <- completed / elapsed
  remaining <- (nrow(param_grid) - completed) / rate
  cat(sprintf('  %.1f sec (ETA: %.0f sec)\n', elapsed, remaining))
}

results <- bind_rows(results_list)
plan(sequential)

# ============================================================================
# SUMMARIZE
# ============================================================================

cat('\n', strrep('=', 70), '\n')
cat('RESULTS SUMMARY\n')
cat(strrep('=', 70), '\n\n')

summary_results <- results |>
  filter(!is.na(scenario)) |>
  group_by(scenario, design, N, c.bm, carryover_t1half) |>
  summarize(
    power = mean(significant, na.rm = TRUE),
    mean_effect = mean(effect, na.rm = TRUE),
    sd_effect = sd(effect, na.rm = TRUE),
    n = n(),
    n_errors = sum(is.na(significant)),
    .groups = 'drop'
  ) |>
  mutate(
    n_success = round(power * n),
    ci_lower = pmax(0,
      (power + 1.96^2 / (2 * n) -
       1.96 * sqrt((power * (1 - power) +
         1.96^2 / (4 * n)) / n)) /
      (1 + 1.96^2 / n)),
    ci_upper = pmin(1,
      (power + 1.96^2 / (2 * n) +
       1.96 * sqrt((power * (1 - power) +
         1.96^2 / (4 * n)) / n)) /
      (1 + 1.96^2 / n))
  )

print(summary_results, n = 100)

# Scenario comparison (S1 vs S2)
cat('\n', strrep('=', 70), '\n')
cat('SCENARIO COMPARISON: S1 (with carryover) vs S2 (without)\n')
cat(strrep('=', 70), '\n\n')

scenario_compare <- summary_results |>
  filter(carryover_t1half > 0, c.bm > 0) |>
  select(scenario, design, N, c.bm, carryover_t1half, power) |>
  pivot_wider(
    names_from = scenario,
    values_from = power
  ) |>
  mutate(
    power_gain = S1_with_carryover - S2_no_carryover,
    pct_gain = round(
      (power_gain / pmax(S2_no_carryover, 0.01)) * 100, 1
    )
  ) |>
  arrange(design, N, carryover_t1half, c.bm)

print(scenario_compare, n = 100)

# ============================================================================
# VISUALIZATION
# ============================================================================

plot_data <- summary_results |>
  filter(!is.na(power)) |>
  mutate(
    c_bm_label = factor(
      sprintf('c.bm=%.1f', c.bm),
      levels = rev(sprintf('c.bm=%.1f',
        sort(unique(summary_results$c.bm))))
    ),
    carryover_label = factor(
      sprintf('t1/2=%.1f wk', carryover_t1half),
      levels = sprintf('t1/2=%.1f wk',
        sort(unique(carryover_t1half)))
    ),
    n_label = factor(paste0('N=', N),
      levels = c('N=35', 'N=70')),
    scenario_label = factor(
      case_when(
        scenario == 'S1_with_carryover' ~
          'S1: With Carryover',
        scenario == 'S2_no_carryover' ~
          'S2: Without Carryover'
      ),
      levels = c('S2: Without Carryover',
                 'S1: With Carryover')
    ),
    power_pct = round(power * 100),
    ci_lower_pct = round(ci_lower * 100),
    ci_upper_pct = round(ci_upper * 100)
  )

for (design_name in unique(plot_data$design)) {
  d_data <- plot_data |> filter(design == design_name)

  if (nrow(d_data) == 0) next

  p <- ggplot(d_data,
    aes(x = carryover_label, y = c_bm_label, fill = power)) +
    geom_tile(color = 'white', linewidth = 0.5) +
    geom_text(
      aes(label = sprintf('%d', power_pct)),
      color = 'black', size = 3.0, fontface = 'bold',
      vjust = -0.3
    ) +
    geom_text(
      aes(label = sprintf('(%d,%d)',
        ci_lower_pct, ci_upper_pct)),
      color = 'black', size = 1.8, vjust = 1.2
    ) +
    scale_fill_gradient2(
      name = 'Power (%)',
      low = '#d73027', mid = '#fee08b', high = '#1a9850',
      midpoint = 0.5, limits = c(0, 1),
      labels = function(x) sprintf('%.0f', x * 100)
    ) +
    facet_grid(
      n_label ~ scenario_label,
      labeller = label_value
    ) +
    labs(
      title = sprintf('Power: %s Design (Hendrickson-Aligned)',
        design_name),
      subtitle = sprintf(
        '%d iterations | DGP sf=%d, analysis sf=%d',
        n_iterations, dgp_scalefactor, analysis_scalefactor
      ),
      x = 'Carryover Half-Life',
      y = 'Biomarker-BR Correlation',
      caption = sprintf(
        'Source: simulation_clustered_hendrickson.R | %s',
        format(Sys.time(), '%Y-%m-%d'))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'right',
      strip.text = element_text(face = 'bold', size = 9),
      strip.text.y = element_text(angle = 0),
      axis.text.x = element_text(
        angle = 45, hjust = 1, size = 8),
      plot.title = element_text(
        face = 'bold', size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.caption = element_text(
        size = 7, hjust = 1, color = 'gray50')
    )

  print(p)

  timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')
  ggsave(
    file.path(output_dir, sprintf(
      'power_hendrickson_%s_%s.pdf',
      tolower(gsub('[^a-zA-Z0-9]', '', design_name)),
      timestamp)),
    p, width = 10, height = 6
  )
  ggsave(
    file.path(output_dir, sprintf(
      'power_hendrickson_%s_%s.png',
      tolower(gsub('[^a-zA-Z0-9]', '', design_name)),
      timestamp)),
    p, width = 10, height = 6, dpi = 300
  )
}

# Type I error check
cat('\n', strrep('=', 70), '\n')
cat('TYPE I ERROR RATES (c.bm = 0)\n')
cat(strrep('=', 70), '\n\n')

type1 <- summary_results |>
  filter(c.bm == 0) |>
  select(scenario, design, N, carryover_t1half, power,
    ci_lower, ci_upper)

print(type1, n = 50)

# ============================================================================
# SAVE
# ============================================================================

timestamp <- format(Sys.time(), '%Y%m%d_%H%M%S')

save(
  results, summary_results, scenario_compare,
  param_grid, n_iterations,
  dgp_scalefactor, analysis_scalefactor,
  file = file.path(output_dir,
    'simulation_clustered_hendrickson_results.RData')
)

cat('\n', strrep('=', 70), '\n')
cat('SIMULATION COMPLETE\n')
cat(strrep('=', 70), '\n')
cat('Results saved to:', output_dir, '\n')
cat('Key output files:\n')
cat('  - simulation_clustered_hendrickson_results.RData\n')
cat('  - power_hendrickson_*_*.pdf/png (per design)\n')
cat('  - simulation_clustered_hendrickson_log.txt\n')

sink()
