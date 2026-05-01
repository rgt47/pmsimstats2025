# Contributing to pmsimstats2025

Thank you for your interest in contributing to this N-of-1 clinical trial
simulation project.

## Development Setup

### Using Docker (Recommended)

```bash
git clone https://github.com/rgt47/pmsimstats2025.git
cd pmsimstats2025
make r  # Starts RStudio Server at localhost:8787
```

### Local R Installation

```bash
git clone https://github.com/rgt47/pmsimstats2025.git
cd pmsimstats2025
```

Then in R:

```r
renv::restore()
```

## Code Style

This project follows the tidyverse style guide with these specifics:

- **Pipe operator**: Use native R pipe `|>` (not magrittr `%>%`)
- **Assignment**: Use `<-` for assignment (not `=`)
- **Naming**: Use `snake_case` for variables and functions
- **Returns**: Use implicit returns (omit `return()` unless for early exit)
- **Documentation**: Use roxygen2 style for function documentation
- **Comments**: Only add "why" comments, not "what" comments

### Example

```r
calculate_power <- function(results, alpha = 0.05) {
  results |>
    filter(!is.na(p_value)) |>
    summarize(
      power = mean(p_value < alpha),
      n_valid = n()
    )
}
```

## Testing

Run tests before submitting changes:

```bash
# All unit tests
Rscript -e 'devtools::test()'

# Single test file
Rscript -e 'testthat::test_file("tests/testthat/test-utils.R")'

# In Docker
make docker-test
```

## Dependency Validation

Before committing, validate that all dependencies are properly tracked:

```bash
make check-renv
```

This ensures packages used in code are listed in both `DESCRIPTION` and
`renv.lock`.

## Pull Request Process

1. **Fork** the repository and create a feature branch from `main`
2. **Make changes** with clear, descriptive commit messages
3. **Run tests** and dependency validation
4. **Submit PR** with a description of changes and motivation

### Commit Message Format

Use descriptive commit messages:

```
Add parameter validation for biomarker correlations

- Implement validate_parameter_grid() function
- Add condition number computation for numerical stability
- Update simulation scripts to use pre-validation
```

## Simulation-Specific Guidelines

### Modifying Correlation Parameters

Changes to correlation values in `pm_functions.R` can affect positive
definiteness. Always:

1. Run `validate_parameter_grid()` after changes
2. Check condition numbers are below 100
3. Test with existing simulation scripts

### Adding New Trial Designs

1. Define measurement schedule in simulation file
2. Add design to `param_grid` expansion
3. Update visualization labels
4. Document the design rationale

## Questions

For questions about the methodology or implementation, please open an issue on
GitHub.

## License

By contributing, you agree that your contributions will be licensed under the
GPL-3 license.
