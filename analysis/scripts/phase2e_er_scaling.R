# Phase 2E: Expectancy Response SD Scaling Impact
# Test ER SD scaling (orig) vs. constant (current 2025)
# In orig: expectancy SD = resp_params$pb$sd * expectancy factor
# In 2025: expectancy SD = constant (no scaling)
# Focus: CO and N-of-1, c.bm=0.3 only

library(tidyverse)
library(lme4)
library(mvtnorm)

set.seed(42)
source('analysis/scripts/pm_functions.R')

cat('\n============================================\n')
cat('PHASE 2E: ER SD Scaling Impact\n')
cat('Scaled by expectancy vs. constant\n')
cat('============================================\n\n')

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

# Modified build_path_sigma with ER scaling toggle
build_path_sigma_er_scaling <- function(
    timepoints,
    ondrug,
    expectancies,
    c.bm,
    carryover_t1half,
    scalefactor = 2,
    scale_er_sd = FALSE,
    resp_params = hendrickson_resp_params,
    bl_params = hendrickson_bl_params,
    corr_params = hendrickson_corr_params,
    make_pd = TRUE) {

  nP <- length(ondrug)
  t_wk <- c(timepoints[1], diff(timepoints))
  tod <- compute_tod_orig(t_wk, ondrug)
  tsd <- compute_tsd_orig(t_wk, ondrug)
  tpb <- compute_tpb_orig(t_wk, expectancies)
  t_cumulative <- cumsum(t_wk)

  timeptnames <- paste0('T', seq_len(nP))
  cl <- c('tv', 'pb', 'br')

  labels <- c(
    'bm', 'BL',
    paste(timeptnames, cl[1], sep = '.'),
    paste(timeptnames, cl[2], sep = '.'),
    paste(timeptnames, cl[3], sep = '.')
  )

  sds <- c(bl_params$bm$sd, bl_params$BL$sd)
  sds <- c(sds, rep(resp_params$tv$sd, nP))

  if (scale_er_sd) {
    # Scale PB SD by expectancy (orig behavior)
    sds <- c(sds, resp_params$pb$sd * expectancies)
  } else {
    # Constant PB SD (current 2025)
    sds <- c(sds, rep(resp_params$pb$sd, nP))
  }

  sds <- c(sds, rep(resp_params$br$sd, nP))

  tv_means <- modgompertz_orig(t_cumulative,
    resp_params$tv$max, resp_params$tv$disp, resp_params$tv$rate)
  pb_means <- modgompertz_orig(tpb,
    resp_params$pb$max, resp_params$pb$disp, resp_params$pb$rate) *
    expectancies

  brmeans <- modgompertz_orig(tod,
    resp_params$br$max, resp_params$br$disp, resp_params$br$rate)
  brtest <- brmeans == 0
  rawbrmeans <- brmeans

  if (nP > 1 && carryover_t1half > 0) {
    for (p in 2:nP) {
      if (ondrug[p] == 0 && tsd[p] > 0) {
        brmeans[p] <- brmeans[p] +
          brmeans[p - 1] * (1/2)^(scalefactor * tsd[p] / carryover_t1half)
      }
    }
  }

  means <- c(bl_params$bm$m, bl_params$BL$m,
             tv_means, pb_means, brmeans)

  correlations <- diag(length(labels))
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  for (cc in cl) {
    ac <- corr_params[[paste0('c.', cc)]]
    if (nP > 1) {
      for (p in 1:(nP - 1)) {
        for (p2 in (p + 1):nP) {
          n1 <- paste(timeptnames[p], cc, sep = '.')
          n2 <- paste(timeptnames[p2], cc, sep = '.')
          lag <- p2 - p
          correlations[n1, n2] <- ac^lag
          correlations[n2, n1] <- ac^lag
        }
      }
    }

    for (c2 in setdiff(cl, cc)) {
      for (p in 1:nP) {
        n1 <- paste(timeptnames[p], cc, sep = '.')
        n2 <- paste(timeptnames[p], c2, sep = '.')
        correlations[n1, n2] <- corr_params$c.cf1t
        correlations[n2, n1] <- corr_params$c.cf1t
      }
      if (nP > 1) {
        for (p in 1:(nP - 1)) {
          for (p2 in (p + 1):nP) {
            n1 <- paste(timeptnames[p], cc, sep = '.')
            n2 <- paste(timeptnames[p2], c2, sep = '.')
            correlations[n1, n2] <- corr_params$c.cfct
            correlations[n2, n1] <- corr_params$c.cfct
          }
        }
      }
    }
  }

  bm_br_corrs <- build_bm_br_correlations(
    brmeans, brtest, c.bm = c.bm, nP = nP
  )
  for (p in 1:nP) {
    n1 <- paste(timeptnames[p], 'br', sep = '.')
    if (p > 1) {
      correlations['bm', n1] <- bm_br_corrs[p]
      correlations[n1, 'bm'] <- bm_br_corrs[p]
    }
  }

  sigma <- outer(sds, sds) * correlations
  was_pd <- corpcor::is.positive.definite(sigma)
  if (make_pd && !was_pd) {
    sigma <- corpcor::make.positive.definite(sigma, tol = 1e-3)
  }

  list(
    sigma = sigma,
    means = means,
    labels = labels,
    sds = sds,
    correlations = correlations,
    indices = list(bm = 1, bl = 2, nP = nP),
    diagnostics = list(
      brmeans = brmeans,
      bm_br_corrs = bm_br_corrs,
      was_pd = was_pd,
      eigenvalues = eigen(sigma, only.values = TRUE)$values
    )
  )
}

# Test grid
test_grid <- expand_grid(
  design_key = c('CO', 'NOF1'),
  N = c(35, 70),
  c.bm = 0.3,
  carryover_t1half = c(0, 0.1, 0.2),
  scale_er = c('no', 'yes')
)

all_results <- list()

for (i in 1:nrow(test_grid)) {
  row <- test_grid[i, ]
  design <- designs[[row$design_key]]
  scale_er_sd <- row$scale_er == 'yes'

  power_count <- 0

  for (iter in 1:75) {
    if (iter %% 25 == 1) cat('.')

    path_result <- build_path_sigma_er_scaling(
      timepoints = design$timepoints,
      ondrug = design$ondrug,
      expectancies = design$expectancies,
      c.bm = row$c.bm,
      carryover_t1half = row$carryover_t1half,
      scalefactor = 2,
      scale_er_sd = scale_er_sd
    )

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
    scale_er = row$scale_er,
    power = power_count / 75
  )

  if (i %% 6 == 0) cat('\n')
}

cat('\n\n=== RESULTS ===\n\n')
results_df <- bind_rows(all_results)

# Reshape for comparison
comparison <- results_df %>%
  pivot_wider(
    id_cols = c(design, N, c.bm, carryover_t1half),
    names_from = scale_er,
    values_from = power
  ) %>%
  mutate(
    delta = yes - no,
    pct_diff = if_else(no > 0, (delta / no) * 100, NA_real_)
  )

print(comparison %>%
  mutate(
    no_str = sprintf('%.0f%%', no * 100),
    yes_str = sprintf('%.0f%%', yes * 100),
    delta_str = sprintf('%+.0f%%', delta * 100),
    pct_diff_str = sprintf('%+.1f%%', pct_diff)
  ) %>%
  select(design, N, carryover_t1half, no_str, yes_str,
         delta_str, pct_diff_str) %>%
  rename('t1/2' = carryover_t1half, 'No Scale' = no_str,
         'With Scale' = yes_str, 'Δ' = delta_str, 'Δ%' = pct_diff_str))

cat('\n\n=== STABILITY ANALYSIS ===\n')
cat('Stability = 1 - |power drop from t1/2=0 to t1/2=0.2|\n\n')

stability <- results_df %>%
  pivot_wider(
    id_cols = c(design, N, scale_er),
    names_from = carryover_t1half,
    values_from = power
  ) %>%
  mutate(
    stability = 1 - abs(`0.2` - `0`)
  ) %>%
  select(design, N, scale_er, stability) %>%
  pivot_wider(
    id_cols = c(design, N),
    names_from = scale_er,
    values_from = stability
  ) %>%
  mutate(stability_delta = yes - no)

print(stability %>%
  mutate(
    no = sprintf('%.3f', no),
    yes = sprintf('%.3f', yes),
    stability_delta = sprintf('%+.3f', stability_delta)
  ) %>%
  rename('No Scale' = no, 'With Scale' = yes))

cat('\n=== KEY FINDING ===\n')
cat('Does ER SD scaling improve stability?\n')
cat('If Δstability > 0, scaling is BENEFICIAL.\n\n')

print(stability %>%
  select(design, N, stability_delta) %>%
  mutate(
    direction = if_else(stability_delta > 0, 'Scaling helps', 'Constant is better'),
    magnitude = abs(stability_delta)
  ) %>%
  rename('Δ Stability' = stability_delta, 'Direction' = direction))

cat('\n')
