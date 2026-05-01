# Phase 2C: Correlation Structure Impact
# Test AR(1) vs Compound Symmetry correlation models
# Focus: CO and N-of-1, c.bm=0.3 only

library(tidyverse)
library(lme4)
library(mvtnorm)

set.seed(42)
source('analysis/scripts/pm_functions.R')

cat('\n============================================\n')
cat('PHASE 2C: Correlation Structure Impact\n')
cat('AR(1) vs Compound Symmetry (CS)\n')
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

# Helper: Build correlation matrix with Compound Symmetry structure
build_cs_correlations <- function(labels, c.tv, c.pb, c.br, c.cf1t, c.cfct) {
  # Compound symmetry: all same-type correlations are constant,
  # all cross-type correlations are constant (ignoring time)

  n <- length(labels)
  correlations <- diag(n)
  rownames(correlations) <- labels
  colnames(correlations) <- labels

  cl <- c('tv', 'pb', 'br')
  timeptnames <- gsub('\\.(tv|pb|br)$', '', labels[grep('\\.(tv|pb|br)$', labels)])
  timeptnames <- unique(timeptnames)

  for (cc in cl) {
    # Within-type correlations (constant across all timepoints)
    pattern <- paste0('\\.', cc, '$')
    indices <- grep(pattern, labels)
    if (length(indices) > 1) {
      for (i in seq_len(length(indices) - 1)) {
        for (j in (i + 1):length(indices)) {
          correlations[indices[i], indices[j]] <- c.tv  # Same const for all types
          correlations[indices[j], indices[i]] <- c.tv
        }
      }
    }

    # Cross-type correlations within same timepoint
    for (tp in timeptnames) {
      for (c2 in setdiff(cl, cc)) {
        n1 <- paste0(tp, '.', cc)
        n2 <- paste0(tp, '.', c2)
        if (n1 %in% labels && n2 %in% labels) {
          correlations[n1, n2] <- c.cf1t
          correlations[n2, n1] <- c.cf1t
        }
      }

      # Cross-type, cross-timepoint
      for (c2 in setdiff(cl, cc)) {
        for (tp2 in timeptnames) {
          if (tp != tp2) {
            n1 <- paste0(tp, '.', cc)
            n2 <- paste0(tp2, '.', c2)
            if (n1 %in% labels && n2 %in% labels) {
              correlations[n1, n2] <- c.cfct
              correlations[n2, n1] <- c.cfct
            }
          }
        }
      }
    }
  }

  correlations
}

# Modified version of build_path_sigma with CS option
build_path_sigma_cs <- function(
    timepoints,
    ondrug,
    expectancies,
    c.bm,
    carryover_t1half,
    scalefactor = 2,
    use_cs = FALSE,
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
  sds <- c(sds, rep(resp_params$pb$sd, nP) * expectancies)
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

  if (use_cs) {
    correlations <- build_cs_correlations(
      labels,
      c.tv = corr_params$c.tv,
      c.pb = corr_params$c.pb,
      c.br = corr_params$c.br,
      c.cf1t = corr_params$c.cf1t,
      c.cfct = corr_params$c.cfct
    )
  } else {
    # Use existing AR(1) structure (baseline)
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
            # AR(1): correlation = ac^(lag)
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
  corr_type = c('AR1', 'CS')
)

all_results <- list()

for (i in 1:nrow(test_grid)) {
  row <- test_grid[i, ]
  design <- designs[[row$design_key]]
  use_cs <- row$corr_type == 'CS'

  power_count <- 0

  for (iter in 1:75) {
    if (iter %% 25 == 1) cat('.')

    path_result <- build_path_sigma_cs(
      timepoints = design$timepoints,
      ondrug = design$ondrug,
      expectancies = design$expectancies,
      c.bm = row$c.bm,
      carryover_t1half = row$carryover_t1half,
      scalefactor = 2,
      use_cs = use_cs
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
    corr_type = row$corr_type,
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
    names_from = corr_type,
    values_from = power
  ) %>%
  mutate(
    delta = CS - AR1,
    pct_diff = if_else(AR1 > 0, (delta / AR1) * 100, NA_real_)
  )

print(comparison %>%
  mutate(
    AR1_str = sprintf('%.0f%%', AR1 * 100),
    CS_str = sprintf('%.0f%%', CS * 100),
    delta_str = sprintf('%+.0f%%', delta * 100),
    pct_diff_str = sprintf('%+.1f%%', pct_diff)
  ) %>%
  select(design, N, carryover_t1half, AR1_str, CS_str,
         delta_str, pct_diff_str) %>%
  rename('t1/2' = carryover_t1half, 'AR(1)' = AR1_str, 'CS' = CS_str,
         'Δ' = delta_str, 'Δ%' = pct_diff_str))

cat('\n\n=== STABILITY ANALYSIS ===\n')
cat('Stability = 1 - |power drop from t1/2=0 to t1/2=0.2|\n\n')

stability <- results_df %>%
  pivot_wider(
    id_cols = c(design, N, corr_type),
    names_from = carryover_t1half,
    values_from = power
  ) %>%
  mutate(
    stability = 1 - abs(`0.2` - `0`)
  ) %>%
  select(design, N, corr_type, stability) %>%
  pivot_wider(
    id_cols = c(design, N),
    names_from = corr_type,
    values_from = stability
  ) %>%
  mutate(stability_delta = CS - AR1)

print(stability %>%
  mutate(
    AR1 = sprintf('%.3f', AR1),
    CS = sprintf('%.3f', CS),
    stability_delta = sprintf('%+.3f', stability_delta)
  ) %>%
  rename('AR(1)' = AR1))

cat('\n=== KEY FINDING ===\n')
cat('Does CS improve stability compared to AR(1)?\n')
cat('If Δstability < 0, CS is WORSE (AR(1) is better).\n\n')

print(stability %>%
  select(design, N, stability_delta) %>%
  mutate(
    direction = if_else(stability_delta > 0, 'CS better', 'AR(1) better'),
    magnitude = abs(stability_delta)
  ) %>%
  rename('Δ Stability' = stability_delta, 'Direction' = direction))

cat('\n')
