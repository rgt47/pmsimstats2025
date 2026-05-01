# Phase 2A: Scalefactor Variant Test
# Compare scalefactor=1 vs scalefactor=2 (baseline)
# Focus: CO and N-of-1 designs only

library(tidyverse)
library(lme4)
library(mvtnorm)

set.seed(42)
source('analysis/scripts/pm_functions.R')

cat('\n==============================================\n')
cat('PHASE 2A: Scalefactor Variant (sf=1 vs sf=2)\n')
cat('==============================================\n')
cat('Sample: CO, N-of-1 | N: 35, 70 | c.bm: 0.3\n')
cat('Iterations: 75 per condition\n\n')

# Trial designs
designs <- list(
  CO = list(
    timepoints = cumsum(c(4, 4, 1, 1, 1, 1, 4, 4)),
    ondrug = c(1, 1, 1, 1, 0, 0, 1, 0),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    name = 'CO'
  ),
  NOF1 = list(
    timepoints = cumsum(c(4, 4, 1, 1, 1, 1, 4, 4)),
    ondrug = c(1, 1, 1, 1, 0, 0, 1, 0),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    name = 'N-of-1'
  )
)

# Test grid: CO and N-of-1 only, both sf variants
test_grid <- expand_grid(
  design_key = c('CO', 'NOF1'),
  N = c(35, 70),
  c.bm = 0.3,
  carryover_t1half = c(0, 0.1, 0.2),
  scalefactor = c(1, 2)
)

all_results <- list()

for (i in 1:nrow(test_grid)) {
  row <- test_grid[i, ]
  design <- designs[[row$design_key]]

  power_count <- 0

  for (iter in 1:75) {
    if (iter %% 25 == 1) cat('.')

    # Build sigma with specified scalefactor
    path_result <- build_path_sigma(
      timepoints = design$timepoints,
      ondrug = design$ondrug,
      expectancies = design$expectancies,
      c.bm = row$c.bm,
      carryover_t1half = row$carryover_t1half,
      scalefactor = row$scalefactor
    )

    # Generate data
    all_data <- list()
    for (pid in 1:row$N) {
      latent <- mvtnorm::rmvnorm(1, mean = path_result$means,
                                 sigma = path_result$sigma)[1, ]
      bl_bm <- latent[1]

      for (tp in 1:8) {
        idx <- (tp - 1) * 3 + 1
        all_data[[length(all_data) + 1]] <- tibble(
          participant = pid,
          timepoint = tp,
          BM = latent[1 + idx],
          Y = latent[1 + idx] + latent[2 + idx] + latent[3 + idx],
          baseline_BM = bl_bm,
          Drug = design$ondrug[tp]
        )
      }
    }

    analysis_data <- bind_rows(all_data) %>%
      mutate(BM_c = BM - baseline_BM)

    tryCatch({
      mod <- lmer(Y ~ BM_c * Drug + (1 | participant),
                 data = analysis_data, REML = FALSE)
      coefs <- fixef(mod)
      if ('BM_c:Drug' %in% names(coefs)) {
        se <- sqrt(vcov(mod)['BM_c:Drug', 'BM_c:Drug'])
        t_stat <- abs(coefs['BM_c:Drug'] / se)
        p_val <- 2 * pnorm(-t_stat)
        if (p_val < 0.05) power_count <- power_count + 1
      }
    }, error = function(e) {})
  }

  all_results[[length(all_results) + 1]] <- tibble(
    design = row$design_key,
    N = row$N,
    c.bm = row$c.bm,
    carryover_t1half = row$carryover_t1half,
    scalefactor = row$scalefactor,
    power = power_count / 75
  )

  if (i %% 6 == 0) cat('\n')
}

cat('\n\n=== RESULTS ===\n\n')
results_df <- bind_rows(all_results)

# Reshape for direct comparison
comparison <- results_df %>%
  pivot_wider(
    id_cols = c(design, N, c.bm, carryover_t1half),
    names_from = scalefactor,
    values_from = power,
    names_prefix = 'sf'
  ) %>%
  mutate(
    delta = sf1 - sf2,
    pct_diff = if_else(sf2 > 0, (delta / sf2) * 100, NA_real_)
  )

print(comparison %>%
  mutate(
    sf1_str = sprintf('%.0f%%', sf1 * 100),
    sf2_str = sprintf('%.0f%%', sf2 * 100),
    delta_str = sprintf('%+.0f%%', delta * 100),
    pct_diff_str = sprintf('%+.1f%%', pct_diff)
  ) %>%
  select(design, N, carryover_t1half, sf1_str, sf2_str,
         delta_str, pct_diff_str) %>%
  rename('t1/2' = carryover_t1half, 'sf=1' = sf1_str, 'sf=2' = sf2_str,
         'Δ' = delta_str, 'Δ%' = pct_diff_str))

cat('\n\n=== STABILITY ANALYSIS ===\n')
cat('Stability = 1 - |power drop from t1/2=0 to t1/2=0.2|\n\n')

stability <- results_df %>%
  pivot_wider(
    id_cols = c(design, N, scalefactor),
    names_from = carryover_t1half,
    values_from = power
  ) %>%
  mutate(
    stability = 1 - abs(`0.2` - `0`)
  ) %>%
  select(design, N, scalefactor, stability) %>%
  pivot_wider(
    id_cols = c(design, N),
    names_from = scalefactor,
    values_from = stability,
    names_prefix = 'sf'
  ) %>%
  mutate(stability_delta = sf1 - sf2)

print(stability %>%
  mutate(
    sf1 = sprintf('%.3f', sf1),
    sf2 = sprintf('%.3f', sf2),
    stability_delta = sprintf('%+.3f', stability_delta)
  ) %>%
  rename('sf=1' = sf1, 'sf=2' = sf2, 'Δ' = stability_delta))

cat('\n=== KEY FINDING ===\n')
cat('Does scalefactor=1 reduce vulnerability compared to sf=2?\n')
cat('If Δstability > 0 for a design, then sf=1 improves power curve stability.\n\n')

for_reporting <- stability %>%
  select(design, N, stability_delta)

print(for_reporting %>%
  mutate(
    stability_delta_str = sprintf('%+.3f', stability_delta),
    interpretation = if_else(stability_delta > 0,
                             'IMPROVED (sf=1 better)',
                             'WORSE (sf=2 better)')
  ) %>%
  select(design, N, stability_delta_str, interpretation) %>%
  rename('Design' = design, 'Δ Stability' = stability_delta_str))

cat('\n')
