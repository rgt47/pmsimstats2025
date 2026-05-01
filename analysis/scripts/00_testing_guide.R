# Testing Guide for Research Projects
# This script provides examples and templates for adding tests to your research repository

cat("=== RESEARCH PROJECT TESTING GUIDE ===\\n")
cat("This guide covers testing at multiple levels:\\n")
cat("1. Unit tests for R package functions\\n")
cat("2. Integration tests for analysis scripts\\n") 
cat("3. Data validation tests\\n")
cat("4. Reproducibility tests\\n\\n")

# Load required packages for testing
library(testthat)
library(here)

# 1. UNIT TESTS FOR PACKAGE FUNCTIONS ====
cat("=== 1. UNIT TESTS FOR PACKAGE FUNCTIONS ===\\n")
cat("Location: tests/testthat/\\n")
cat("Framework: testthat package\\n\\n")

cat("Example test file structure:\\n")
cat("tests/testthat/test-my-function.R:\\n")
cat("\\n")
cat("test_that(\\"my_function works correctly\\", {\\n")
cat("  # Test basic functionality\\n")
cat("  result <- my_function(input_data)\\n")
cat("  expect_equal(nrow(result), expected_rows)\\n")
cat("  expect_true(all(result$column > 0))\\n")
cat("})\\n\\n")

cat("test_that(\\"my_function handles edge cases\\", {\\n")
cat("  # Test with empty input\\n")
cat("  expect_error(my_function(data.frame()))\\n")
cat("  \\n")
cat("  # Test with NA values\\n")
cat("  data_with_na <- data.frame(x = c(1, NA, 3))\\n")
cat("  result <- my_function(data_with_na)\\n")
cat("  expect_false(any(is.na(result)))\\n")
cat("})\\n\\n")

cat("To run package tests:\\n")
cat("- devtools::test()     # Run all tests\\n")
cat("- testthat::test_file(\\"tests/testthat/test-my-function.R\\")  # Run specific test\\n")
cat("- make test           # Using Makefile\\n\\n")

# 2. INTEGRATION TESTS FOR ANALYSIS SCRIPTS ====
cat("=== 2. INTEGRATION TESTS FOR ANALYSIS SCRIPTS ===\\n")
cat("Location: tests/integration/\\n")
cat("Purpose: Test that analysis scripts run without errors\\n\\n")

# Create integration test directory and example
dir.create(here("tests", "integration"), recursive = TRUE, showWarnings = FALSE)

integration_test_example <- "# Integration Test Example
# tests/integration/test-analysis-pipeline.R

# Test that main analysis scripts can run without errors
test_that(\\"data validation script runs without errors\\", {
  expect_no_error({
    source(here(\\"scripts\\", \\"02_data_validation.R\\"))
  })
})

test_that(\\"parallel setup configures correctly\\", {
  expect_no_error({
    source(here(\\"scripts\\", \\"00_setup_parallel.R\\"))
  })
})

# Test with mock data
test_that(\\"analysis works with test data\\", {
  # Create mock dataset
  test_data <- data.frame(
    id = 1:100,
    treatment = sample(c(\\"A\\", \\"B\\"), 100, replace = TRUE),
    outcome = rnorm(100, mean = 50, sd = 10)
  )
  
  # Save test data
  temp_file <- tempfile(fileext = \\".csv\\")
  write.csv(test_data, temp_file, row.names = FALSE)
  
  # Test analysis functions
  # Add your specific analysis tests here
  expect_true(nrow(test_data) == 100)
  
  # Cleanup
  unlink(temp_file)
})"

cat("\\nExample integration test:\\n")
cat(integration_test_example)
cat("\\n\\n")

# 3. DATA VALIDATION TESTS ====
cat("=== 3. DATA VALIDATION TESTS ===\\n")
cat("Location: tests/data/\\n")
cat("Purpose: Ensure data quality and consistency\\n\\n")

data_test_example <- "# Data Validation Test Example
# tests/data/test-data-quality.R

test_that(\\"raw data meets quality standards\\", {
  # Skip if no data files exist
  data_dir <- here(\\"data\\", \\"raw_data\\")
  skip_if_not(dir.exists(data_dir), \\"No raw data directory found\\")
  
  data_files <- list.files(data_dir, pattern = \\"\\\\.csv$\\", full.names = TRUE)
  skip_if(length(data_files) == 0, \\"No CSV files found in raw data\\")
  
  for (file in data_files) {
    data <- read.csv(file)
    
    # Test basic properties
    expect_true(nrow(data) > 0, info = paste(\\"Empty data file:\\", basename(file)))
    expect_true(ncol(data) > 0, info = paste(\\"No columns in:\\", basename(file)))
    
    # Test for required columns (customize for your data)
    # expect_true(\\"id\\" %in% names(data), info = \\"Missing id column\\")
    
    # Test data types and ranges (customize for your data)
    # if(\\"age\\" %in% names(data)) {
    #   expect_true(all(data$age >= 0 & data$age <= 120, na.rm = TRUE))
    # }
  }
})"

cat("Example data validation test:\\n")
cat(data_test_example)
cat("\\n\\n")

# 4. REPRODUCIBILITY TESTS ====
cat("=== 4. REPRODUCIBILITY TESTS ===\\n")
cat("Location: tests/reproducibility/\\n")
cat("Purpose: Ensure analysis can be reproduced\\n\\n")

repro_test_example <- "# Reproducibility Test Example
# tests/reproducibility/test-reproducibility.R

test_that(\\"analysis environment is reproducible\\", {
  # Test renv lockfile exists
  expect_true(file.exists(\\"renv.lock\\"), \\"renv.lock file missing\\")
  
  # Test package versions are locked
  if (require(\\"renv\\", quietly = TRUE)) {
    status <- renv::status()
    expect_true(length(status) == 0, \\"Package environment not synchronized\\")
  }
})

test_that(\\"analysis scripts have consistent output\\", {
  # Test that scripts produce consistent results
  # This is especially important for analyses with random components
  
  set.seed(12345)  # Set seed for reproducibility
  
  # Run analysis script
  # result1 <- source(here(\\"scripts\\", \\"my_analysis.R\\"))
  
  set.seed(12345)  # Reset seed
  
  # Run again
  # result2 <- source(here(\\"scripts\\", \\"my_analysis.R\\"))
  
  # Compare results
  # expect_equal(result1, result2, info = \\"Analysis results not reproducible\\")
})"

cat("Example reproducibility test:\\n")
cat(repro_test_example)
cat("\\n\\n")

# 5. TESTING WORKFLOW ====
cat("=== 5. TESTING WORKFLOW ===\\n")
cat("Recommended testing workflow:\\n\\n")

cat("1. DEVELOPMENT CYCLE:\\n")
cat("   - Write function in R/\\n")
cat("   - Write test in tests/testthat/\\n")
cat("   - Run devtools::test() to verify\\n")
cat("   - Iterate until tests pass\\n\\n")

cat("2. ANALYSIS TESTING:\\n")
cat("   - Create integration tests for scripts\\n")
cat("   - Test with sample/mock data\\n")
cat("   - Validate data quality regularly\\n")
cat("   - Check reproducibility periodically\\n\\n")

cat("3. CI/CD TESTING:\\n")
cat("   - GitHub Actions runs tests automatically\\n")
cat("   - Tests run on multiple R versions\\n")
cat("   - Tests run in clean environment\\n")
cat("   - Failures block merges to main branch\\n\\n")

cat("4. USEFUL TESTING COMMANDS:\\n")
cat("   devtools::test()                    # Run all package tests\\n")
cat("   devtools::check()                   # Full package check\\n")
cat("   testthat::test_dir(\\"tests/data\\")    # Run data tests\\n")
cat("   source(\\"analysis/scripts/99_reproducibility_check.R\\")  # Check reproducibility\\n")
cat("   make test                          # Run via Makefile\\n\\n")

cat("5. TEST COVERAGE:\\n")
cat("   # Install covr package for coverage analysis\\n")
cat("   library(covr)\\n")
cat("   coverage <- package_coverage()\\n")
cat("   report(coverage)  # Generate HTML coverage report\\n\\n")

cat("=== TESTING SETUP COMPLETE ===\\n")
cat("Next steps:\\n")
cat("1. Add specific tests for your functions in tests/testthat/\\n")
cat("2. Create integration tests for your analysis scripts\\n")
cat("3. Set up data validation tests for your datasets\\n")
cat("4. Run tests regularly during development\\n")
cat("5. Check test coverage with covr package\\n")
