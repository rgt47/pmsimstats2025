# validate_hendrickson_alignment.R
# Phase 3A: Variable chain validation against Hendrickson et al. (2020)
#
# For each design and path, compute and compare key intermediate
# variables (tod, tsd, brmeans, BM-BR correlations, Dbc) against
# hand-computed reference values from the original vignette.

suppressPackageStartupMessages({
  library(tidyverse)
})

source('analysis/scripts/pm_functions.R')

cat(strrep('=', 70), '\n')
cat('HENDRICKSON ALIGNMENT VALIDATION\n')
cat(strrep('=', 70), '\n\n')

designs <- list(
  `N-of-1 (Hybrid)` = list(
    timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    timeptnames = c('OL1', 'OL2', 'BD1', 'BD2',
                    'BD3', 'BD4', 'COd', 'COp'),
    paths = list(
      pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
      pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
      pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
      pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
    )
  ),
  `Open-Label` = list(
    timepoints = cumsum(rep(2.5, 8)),
    expectancies = rep(1, 8),
    timeptnames = paste0('OL', 1:8),
    paths = list(
      pathA = rep(1, 8)
    )
  ),
  `OL+BDC` = list(
    timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
    expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
    timeptnames = c('OL1', 'OL2', 'OL3', 'OL4',
                    'BD1', 'BD2', 'BD3', 'BD4'),
    paths = list(
      pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
      pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
    )
  ),
  `Crossover` = list(
    timepoints = cumsum(rep(2.5, 8)),
    expectancies = rep(0.5, 8),
    timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
    paths = list(
      pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
      pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
    )
  )
)

carryover_levels <- c(0, 0.1, 0.2)
c_bm_levels <- c(0, 0.3)

for (design_name in names(designs)) {
  design <- designs[[design_name]]
  cat(strrep('-', 70), '\n')
  cat(sprintf('Design: %s\n', design_name))
  cat(strrep('-', 70), '\n')
  cat(sprintf('Timepoints: %s\n',
    paste(design$timepoints, collapse = ', ')))
  cat(sprintf('Expectancies: %s\n',
    paste(design$expectancies, collapse = ', ')))

  t_wk <- c(design$timepoints[1], diff(design$timepoints))
  cat(sprintf('Intervals (t_wk): %s\n',
    paste(t_wk, collapse = ', ')))

  for (path_name in names(design$paths)) {
    ondrug <- design$paths[[path_name]]
    cat(sprintf('\n  %s: ondrug = [%s]\n', path_name,
      paste(ondrug, collapse = ', ')))

    tod <- compute_tod_orig(t_wk, ondrug)
    tsd <- compute_tsd_orig(t_wk, ondrug)
    tpb <- compute_tpb_orig(t_wk, design$expectancies)

    cat(sprintf('    tod = [%s]\n', paste(tod, collapse = ', ')))
    cat(sprintf('    tsd = [%s]\n', paste(tsd, collapse = ', ')))
    cat(sprintf('    tpb = [%s]\n', paste(tpb, collapse = ', ')))

    br_raw <- modgompertz_orig(tod,
      hendrickson_resp_params$br$max,
      hendrickson_resp_params$br$disp,
      hendrickson_resp_params$br$rate)
    cat(sprintf('    brmeans_raw = [%s]\n',
      paste(round(br_raw, 3), collapse = ', ')))

    for (t1half in carryover_levels) {
      result <- build_path_sigma(
        timepoints = design$timepoints,
        ondrug = ondrug,
        expectancies = design$expectancies,
        c.bm = 0.3,
        carryover_t1half = t1half
      )
      cat(sprintf('    t1half=%.1f: brmeans_adj = [%s]\n',
        t1half,
        paste(round(result$diagnostics$brmeans, 3),
          collapse = ', ')))
      cat(sprintf('             bm_br_corrs = [%s]\n',
        paste(round(result$diagnostics$bm_br_corrs, 4),
          collapse = ', ')))

      if (t1half > 0) {
        dbc_analysis <- ifelse(ondrug == 1, 1,
          (1/2)^(1 * tsd / t1half))
        cat(sprintf('             Dbc(sf=1)   = [%s]\n',
          paste(round(dbc_analysis, 4), collapse = ', ')))
      }

      cat(sprintf('             PD: %s | min_eig: %.4e\n',
        result$diagnostics$was_pd,
        min(result$diagnostics$eigenvalues)))
    }
  }
  cat('\n')
}

cat(strrep('=', 70), '\n')
cat('GOMPERTZ COMPARISON: mod_gompertz vs modgompertz_orig\n')
cat(strrep('=', 70), '\n\n')

t_test <- seq(0, 20, by = 2.5)
cat(sprintf('%-6s  %-12s  %-12s  %-10s\n',
  't', 'orig', 'pmsim2025', 'diff'))
for (t_val in t_test) {
  v_orig <- modgompertz_orig(t_val, 10.98604, 5, 0.42)
  v_new <- mod_gompertz(t_val, 10.98604, 5, 0.42)
  cat(sprintf('%5.1f   %11.5f   %11.5f   %9.5f\n',
    t_val, v_orig, v_new, v_orig - v_new))
}

cat('\nKey insight: mod_gompertz(0, ...) =',
  round(mod_gompertz(0, 10.98604, 5, 0.42), 5),
  '(should be 0)\n')
cat('         modgompertz_orig(0, ...) =',
  round(modgompertz_orig(0, 10.98604, 5, 0.42), 5),
  '\n')

cat('\n')
cat(strrep('=', 70), '\n')
cat('VALIDATION COMPLETE\n')
cat(strrep('=', 70), '\n')
