# diagnose_carryover_power.R
# Phase 3C: Diagnostic script for carryover power dropoff
#
# For each carryover level, inspects sigma matrix properties to
# understand how carryover affects the BM-BR correlation structure
# and thus statistical power.

suppressPackageStartupMessages({
  library(tidyverse)
  library(MASS)
})

source('analysis/scripts/pm_functions.R')

cat(strrep('=', 70), '\n')
cat('CARRYOVER POWER DIAGNOSTICS\n')
cat(strrep('=', 70), '\n\n')

carryover_levels <- c(0, 0.1, 0.2)
hybrid_path_a <- c(1, 1, 1, 1, 0, 0, 1, 0)
timepoints <- c(4, 8, 9, 10, 11, 12, 16, 20)
expectancies <- c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

for (t1half in carryover_levels) {
  cat(strrep('-', 70), '\n')
  cat(sprintf('Carryover t1half = %.1f weeks\n', t1half))
  cat(strrep('-', 70), '\n\n')

  result <- build_path_sigma(
    timepoints = timepoints,
    ondrug = hybrid_path_a,
    expectancies = expectancies,
    c.bm = 0.3,
    carryover_t1half = t1half
  )

  cat('BR means (raw vs adjusted):\n')
  br_df <- tibble(
    tp = 1:8,
    ondrug = hybrid_path_a,
    raw = round(result$diagnostics$rawbrmeans, 4),
    adjusted = round(result$diagnostics$brmeans, 4),
    diff = round(
      result$diagnostics$brmeans - result$diagnostics$rawbrmeans,
      4
    )
  )
  print(br_df, n = 8)

  cat('\nBM-BR correlations:\n')
  cat(sprintf('  [%s]\n',
    paste(round(result$diagnostics$bm_br_corrs, 4),
      collapse = ', ')))

  eigs <- result$diagnostics$eigenvalues
  cat(sprintf('\nEigenvalues (top 5): %s\n',
    paste(round(sort(eigs, decreasing = TRUE)[1:5], 4),
      collapse = ', ')))
  cat(sprintf('Eigenvalues (bottom 5): %s\n',
    paste(round(sort(eigs)[1:5], 4), collapse = ', ')))
  cat(sprintf('Condition number: %.1f\n', max(eigs) / min(eigs)))
  cat(sprintf('Was PD before enforcement: %s\n',
    result$diagnostics$was_pd))

  if (t1half > 0) {
    ref <- build_path_sigma(
      timepoints = timepoints,
      ondrug = hybrid_path_a,
      expectancies = expectancies,
      c.bm = 0.3,
      carryover_t1half = 0
    )
    sigma_diff <- result$sigma - ref$sigma
    cat(sprintf(
      '\nSigma diff from t1half=0: max_abs=%.4e, mean_abs=%.4e\n',
      max(abs(sigma_diff)),
      mean(abs(sigma_diff))
    ))
  }
  cat('\n')
}

cat(strrep('=', 70), '\n')
cat('MINI-SIMULATION: Power at each carryover level\n')
cat(strrep('=', 70), '\n\n')

n_iter <- 50
N <- 70

for (t1half in carryover_levels) {
  cat(sprintf('t1half=%.1f: ', t1half))

  paths <- list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )

  path_sigmas <- lapply(paths, function(od) {
    build_path_sigma(
      timepoints = timepoints,
      ondrug = od,
      expectancies = expectancies,
      c.bm = 0.3,
      carryover_t1half = t1half
    )
  })

  pvals <- numeric(n_iter)
  for (iter in 1:n_iter) {
    set.seed(iter * 1000)
    n_paths <- length(paths)
    base_per_path <- N %/% n_paths
    remainder <- N %% n_paths

    all_data <- list()
    for (path_idx in seq_along(paths)) {
      n_this <- base_per_path + (if (path_idx <= remainder) 1 else 0)
      ps <- path_sigmas[[path_idx]]
      dat <- mvrnorm(
        n = n_this,
        mu = ps$means,
        Sigma = ps$sigma
      )
      if (n_this == 1) dat <- matrix(dat, nrow = 1)
      dat <- as_tibble(dat, .name_repair = 'minimal')
      colnames(dat) <- ps$labels

      nP <- ps$indices$nP
      timeptnames <- ps$path_info$timeptnames
      ondrug <- ps$path_info$ondrug
      t_cumulative <- ps$path_info$t_cumulative
      tsd_vec <- ps$path_info$tsd

      for (pid in 1:n_this) {
        global_pid <- sum(
          base_per_path * (path_idx - 1) +
          min(path_idx - 1, remainder)
        ) + pid

        for (tp in 1:nP) {
          br_cols <- paste(timeptnames[tp], 'br', sep = '.')
          tv_cols <- paste(timeptnames[tp], 'tv', sep = '.')
          pb_cols <- paste(timeptnames[tp], 'pb', sep = '.')
          Sx <- dat$BL[pid] - (
            dat[[tv_cols]][pid] +
            dat[[pb_cols]][pid] +
            dat[[br_cols]][pid]
          )
          Dbc <- if (ondrug[tp] == 1) {
            1
          } else if (t1half == 0 || tsd_vec[tp] == 0) {
            0
          } else {
            (1/2)^(1 * tsd_vec[tp] / t1half)
          }
          all_data[[length(all_data) + 1]] <- tibble(
            ptID = global_pid,
            bm = dat$bm[pid],
            Sx = Sx,
            t = t_cumulative[tp],
            Dbc = Dbc,
            tsd = tsd_vec[tp]
          )
        }
      }
    }
    merged <- bind_rows(all_data)

    fit <- tryCatch({
      lme4::lmer(Sx ~ bm + t + Dbc + bm * Dbc + (1 | ptID),
        data = merged)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      coefs <- summary(fit)$coefficients
      if ('bm:Dbc' %in% rownames(coefs)) {
        t_val <- coefs['bm:Dbc', 't value']
        pvals[iter] <- 2 * pnorm(-abs(t_val))
      } else {
        pvals[iter] <- NA
      }
    } else {
      pvals[iter] <- NA
    }
  }
  power <- mean(pvals < 0.05, na.rm = TRUE)
  n_fail <- sum(is.na(pvals))
  cat(sprintf('power=%.0f%% (%d failures)\n',
    power * 100, n_fail))
}

cat('\n')
cat(strrep('=', 70), '\n')
cat('DIAGNOSTICS COMPLETE\n')
cat(strrep('=', 70), '\n')
