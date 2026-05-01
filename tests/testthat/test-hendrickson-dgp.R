# test-hendrickson-dgp.R
# Unit tests for Hendrickson-aligned DGP functions (Phase 3B)

source(file.path(
  rprojroot::find_root(rprojroot::is_r_package),
  'analysis', 'scripts', 'pm_functions.R'
))

test_that('modgompertz_orig returns 0 at t=0', {
  expect_equal(modgompertz_orig(0, 10, 5, 0.42), 0)
  expect_equal(modgompertz_orig(0, 6.50647, 5, 0.35), 0)
})

test_that('modgompertz_orig asymptotes at maxr', {
  expect_equal(
    modgompertz_orig(100, 10, 5, 0.42), 10,
    tolerance = 0.01
  )
  expect_equal(
    modgompertz_orig(100, 10.98604, 5, 0.42), 10.98604,
    tolerance = 0.01
  )
})

test_that('modgompertz_orig is monotonically increasing', {
  t_seq <- seq(0, 20, by = 0.5)
  vals <- modgompertz_orig(t_seq, 10.98604, 5, 0.42)
  expect_true(all(diff(vals) >= 0))
})

test_that('compute_tsd_orig matches Hybrid Path A', {
  t_wk <- c(4, 4, 1, 1, 1, 1, 4, 4)
  od <- c(1, 1, 1, 1, 0, 0, 1, 0)
  result <- compute_tsd_orig(t_wk, od)
  expect_equal(result, c(0, 0, 0, 0, 1, 2, 0, 4))
})

test_that('compute_tsd_orig zeroes pre-drug periods', {
  t_wk <- c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5)
  od <- c(0, 0, 0, 0, 1, 1, 1, 1)
  result <- compute_tsd_orig(t_wk, od)
  # Before ever on drug: tsd should be 0
  expect_equal(result[1:4], c(0, 0, 0, 0))
  # After on drug: on drug so tsd = 0

  expect_equal(result[5:8], c(0, 0, 0, 0))
})

test_that('compute_tsd_orig for CO Path B (off then on)', {
  t_wk <- c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5)
  od <- c(0, 0, 0, 0, 1, 1, 1, 1)
  result <- compute_tsd_orig(t_wk, od)
  expect_equal(result, c(0, 0, 0, 0, 0, 0, 0, 0))
})

test_that('compute_tod_orig matches Hybrid Path A', {
  t_wk <- c(4, 4, 1, 1, 1, 1, 4, 4)
  od <- c(1, 1, 1, 1, 0, 0, 1, 0)
  result <- compute_tod_orig(t_wk, od)
  expect_equal(result, c(4, 8, 9, 10, 0, 0, 4, 0))
})

test_that('compute_tod_orig matches CO Path A', {
  t_wk <- rep(2.5, 8)
  od <- c(1, 1, 1, 1, 0, 0, 0, 0)
  result <- compute_tod_orig(t_wk, od)
  expect_equal(result, c(2.5, 5, 7.5, 10, 0, 0, 0, 0))
})

test_that('compute_tpb_orig accumulates correctly', {
  t_wk <- c(4, 4, 1, 1, 1, 1, 4, 4)
  exp <- c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  result <- compute_tpb_orig(t_wk, exp)
  expect_equal(result, c(4, 8, 9, 10, 11, 12, 16, 20))
})

test_that('Ron Thomas corr: on-drug periods get c.bm', {
  brmeans <- c(4, 8, 9, 10, 0.01, 0, 4, 0)
  brtest <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)
  corrs <- build_bm_br_correlations(brmeans, brtest, c.bm = 0.3, nP = 8)
  expect_equal(corrs[1], 0)
  expect_equal(corrs[2], 0.3)
  expect_equal(corrs[3], 0.3)
  expect_equal(corrs[4], 0.3)
  expect_equal(corrs[7], 0.3)
})

test_that('Ron Thomas corr: off-drug with zero brmeans gets 0', {
  brmeans <- c(4, 8, 9, 10, 0.01, 0, 4, 0)
  brtest <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)
  corrs <- build_bm_br_correlations(brmeans, brtest, c.bm = 0.3, nP = 8)
  expect_equal(corrs[6], 0)
  expect_equal(corrs[8], 0)
})

test_that('Ron Thomas corr: off-drug with carryover gets scaled', {
  brmeans <- c(4, 8, 9, 10, 0.5, 0, 4, 0)
  brtest <- c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE)
  corrs <- build_bm_br_correlations(brmeans, brtest, c.bm = 0.3, nP = 8)
  # p=5: brtest=TRUE, brmeans=0.5 != 0, so (0.5/10)*0.3 = 0.015
  expect_equal(corrs[5], (0.5 / 10) * 0.3)
})

test_that('build_path_sigma returns valid sigma for Hybrid Path A', {
  result <- build_path_sigma(
    timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    ondrug = c(1, 1, 1, 1, 0, 0, 1, 0),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    c.bm = 0.3,
    carryover_t1half = 0.1
  )
  expect_true(corpcor::is.positive.definite(result$sigma))
  expect_equal(length(result$means), 26)  # 2 + 3*8
  expect_equal(result$means[1], hendrickson_bl_params$bm$m)
})

test_that('build_path_sigma with no carryover has zero off-drug brmeans', {
  result <- build_path_sigma(
    timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    ondrug = c(1, 1, 1, 1, 0, 0, 1, 0),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    c.bm = 0.3,
    carryover_t1half = 0
  )
  # Off-drug timepoints (5,6,8) should have brmeans = 0
  expect_equal(result$diagnostics$brmeans[5], 0)
  expect_equal(result$diagnostics$brmeans[6], 0)
  expect_equal(result$diagnostics$brmeans[8], 0)
  # On-drug timepoints should have brmeans > 0
  expect_true(result$diagnostics$brmeans[1] > 0)
})

test_that('build_path_sigma: c.bm=0 produces zero BM-BR correlations', {
  result <- build_path_sigma(
    timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
    ondrug = c(1, 1, 1, 1, 0, 0, 1, 0),
    expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    c.bm = 0,
    carryover_t1half = 0.1
  )
  expect_true(all(result$diagnostics$bm_br_corrs == 0))
})

test_that('build_path_sigma dimensions are correct for OL', {
  result <- build_path_sigma(
    timepoints = cumsum(rep(2.5, 8)),
    ondrug = rep(1, 8),
    expectancies = rep(1, 8),
    c.bm = 0.3,
    carryover_t1half = 0
  )
  # 2 baseline + 3*8 factors = 26
  expect_equal(nrow(result$sigma), 26)
  expect_equal(ncol(result$sigma), 26)
  expect_equal(length(result$labels), 26)
})

test_that('OL design: all tod > 0, all tsd = 0', {
  result <- build_path_sigma(
    timepoints = cumsum(rep(2.5, 8)),
    ondrug = rep(1, 8),
    expectancies = rep(1, 8),
    c.bm = 0.3,
    carryover_t1half = 0
  )
  expect_true(all(result$path_info$tod > 0))
  expect_true(all(result$path_info$tsd == 0))
})
