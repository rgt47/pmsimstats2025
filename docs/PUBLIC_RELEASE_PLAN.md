# Public Repository Release Plan

This document outlines the cleanup and preparation steps required before making
the `pmsimstats2025` repository public on GitHub.

## Executive Summary

The repository is a research compendium for N-of-1 clinical trial simulations.
Before public release, it requires metadata corrections, file cleanup, test
fixes, and documentation consolidation.

---

## Priority 1: Critical Fixes (Must Do)

### 1.1 Fix DESCRIPTION File Metadata

**Current state**: DONE. DESCRIPTION metadata has been corrected (package name,
authors, description).

**File**: `DESCRIPTION`

**Required changes**:

```
Package: pmsimstats2025
Title: N-of-1 Clinical Trial Simulation with Biomarker Interactions
Version: 0.1.0
Authors@R: c(
    person("Rebecca C.", "Hendrickson", email = "rebecca.hendrickson@va.gov",
           role = "aut"),
    person("Ronald G.", "Thomas", email = "rgthomas@ucsd.edu",
           role = c("aut", "cre"))
    )
Description: Monte Carlo simulation framework for comparing N-of-1 clinical
    trial designs (Hybrid, Crossover, Open-Label, Parallel) with focus on
    statistical power when biomarker-treatment interactions and carryover
    effects are present. Based on Hendrickson et al. (2020) methodology.
```

### 1.2 Fix testthat.R Package Reference

**Current state**: DONE. `tests/testthat.R` now correctly references
`pmsimstats2025`.

**File**: `tests/testthat.R`

**Current content**:

```r
library(testthat)
library(d08)
test_check("d08")
```

**Required change**:

```r
library(testthat)
library(pmsimstats2025)
test_check("pmsimstats2025")
```

### 1.3 Remove CLAUDE.md from .gitignore

**Current state**: `.gitignore` contains `CLAUDE.md` on line 69, preventing it
from being committed.

**File**: `.gitignore`

**Action**: Remove the line `CLAUDE.md` from `.gitignore` if you want the file
to be part of the public repository. Alternatively, keep it ignored if you
prefer CLAUDE.md to remain local-only.

---

## Priority 2: File Cleanup (Should Do)

### 2.1 Remove Temporary and Test Files

**Files to delete**:

```
.zz0098.R                                    # Root directory temp file
analysis/scripts/.zzvim_r_temp.R             # Vim temp file
analysis/scripts/.zz8571.R                   # Temp file
analysis/report/.zzvim_r_temp.R              # Vim temp file (if exists)
```

**Command**:

```bash
rm -f .zz0098.R
rm -f analysis/scripts/.zzvim_r_temp.R
rm -f analysis/scripts/.zz8571.R
rm -f analysis/report/.zzvim_r_temp.R
```

### 2.2 Add Temp Files to .gitignore

Add these patterns to `.gitignore`:

```
# Vim temporary files
.zzvim_r_temp.R
.zz*.R

# Claude local files (optional - remove if you want CLAUDE.md public)
# CLAUDE.md
```

### 2.3 Review and Consolidate Test Scripts

**Current state**: Multiple `test_*.R` files in `analysis/scripts/` appear to
be development/exploratory scripts rather than formal tests:

- `test_pd_boundaries.R`
- `test_new_params.R`
- `test_boundary_autocorr_06.R`
- `test_cbm_06.R`
- `test_autocorr_for_cbm_06.R`

**Recommendation**: Either:

1. Move to `tests/exploratory/` directory (create it)
2. Delete if no longer needed
3. Rename to `explore_*.R` to clarify they are not formal tests

### 2.4 Remove Redundant Symlink

**File**: `ZZCOLLAB_USER_GUIDE.md` (symlink to `docs/zzcollab-user-guide.md`)

**Action**: Delete the symlink; users can find the guide in `docs/`.

```bash
rm ZZCOLLAB_USER_GUIDE.md
```

### 2.5 Remove Duplicate claude.md

**File**: `claude.md` (lowercase, 20KB)

**Current state**: Both `CLAUDE.md` and `claude.md` exist. The lowercase version
appears to be older/duplicate.

**Action**: Delete `claude.md` after confirming `CLAUDE.md` is the authoritative
version.

```bash
rm claude.md
```

---

## Priority 3: Documentation Cleanup (Good Practice)

### 3.1 Consolidate Markdown Documentation

**Current state**: 30+ markdown files spread across `docs/` and
`analysis/scripts/`. Many appear to be development notes.

**Recommendation**: Create a clear hierarchy:

```
docs/
├── README.md                    # Index of documentation
├── methodology/
│   ├── simulation_white_paper.md
│   ├── hendrickson-alignment-changes.md
│   └── correlation_structure_discussion.md
├── technical/
│   ├── positive_definiteness_constraints.tex
│   ├── sigma_matrix_derivation.tex
│   └── biomarker_interaction_mechanism.tex
├── guides/
│   ├── data-workflow-guide.md
│   └── zzcollab-user-guide.md
└── archive/                     # Development notes (optional)
    └── ... (historical notes)
```

### 3.2 Review SIMPLIFICATION_TUTORIAL.md

**Location**: Root directory

**Action**: Determine if this should be:

- Moved to `docs/guides/`
- Deleted if outdated
- Kept in root if it's a primary onboarding document

### 3.3 Add CONTRIBUTING.md

Create a standard `CONTRIBUTING.md` file for external contributors:

```markdown
# Contributing to pmsimstats2025

## Development Setup

1. Clone the repository
2. Run `make r` to start the Docker environment
3. Run `make check-renv` to validate dependencies

## Code Style

- Use native R pipe `|>` (not `%>%`)
- Use `<-` for assignment
- Follow snake_case naming
- Use roxygen2 for function documentation

## Testing

Run tests before submitting PRs:

```bash
make docker-test
```

## Pull Request Process

1. Create a feature branch
2. Make changes with clear commit messages
3. Ensure tests pass
4. Submit PR with description of changes
```

### 3.4 Add CODE_OF_CONDUCT.md

Consider adding a standard code of conduct (e.g., Contributor Covenant).

---

## Priority 4: Dependency and Environment Cleanup

### 4.1 Audit DESCRIPTION Imports

**Audit completed 2025-01-28**. Packages fall into two categories:

**Core simulation packages (KEEP)**:

- `tidyverse`, `dplyr`, `ggplot2`, `readr` - Data manipulation and visualization
- `MASS`, `corpcor` - Matrix operations and positive definiteness
- `lme4`, `lmerTest`, `nlme` - Mixed-effects modeling
- `conflicted` - Namespace conflict resolution
- `viridis`, `scales`, `patchwork` - Visualization
- `here` - Path management
- `foreach`, `doParallel`, `future`, `furrr` - Parallel processing
- `testthat`, `devtools`, `covr` - Testing and development
- `rmarkdown`, `knitr`, `bookdown` - Report generation
- `renv` - Dependency management

**Template/example packages (CONSIDER REMOVING)**:

These are from the zzcollab scaffold and only used in template scripts:

- `DBI`, `RSQLite`, `RPostgres`, `RMySQL`, `odbc` - Database setup template
  (`00_database_setup.R`)
- `palmerpenguins` - Example dataset (`02_data_validation.R`,
  `test-data-pipeline.R`)
- `visdat`, `naniar`, `skimr`, `janitor` - Data validation template
  (`02_data_validation.R`)
- `emmeans` - Not used in any script
- `digest`, `sessioninfo` - Utility packages, optional

**Recommendation**: If you want a minimal installation footprint, remove the
template packages and their associated scripts. Otherwise, keep them for
users who want the full research compendium scaffold.

**To remove template packages**:

1. Delete template scripts:
   - `analysis/scripts/00_database_setup.R`
   - `analysis/scripts/02_data_validation.R`
   - `tests/integration/test-data-pipeline.R` (uses palmerpenguins)

2. Remove from DESCRIPTION Imports:
   - `DBI`, `RSQLite`, `RPostgres`, `RMySQL`, `odbc`
   - `palmerpenguins`, `visdat`, `naniar`, `skimr`, `janitor`
   - `emmeans`, `digest`, `sessioninfo`

3. Run `renv::snapshot()` to update renv.lock

### 4.2 Verify renv.lock Consistency

```bash
make check-renv
```

Ensure all packages in DESCRIPTION are in renv.lock and vice versa.

### 4.3 Test Docker Build

```bash
make docker-build
make docker-test
```

Verify the Docker image builds successfully and tests pass.

---

## Priority 5: GitHub Configuration

### 5.1 Add Repository Topics

When making the repository public, add topics for discoverability:

- `n-of-1-trials`
- `clinical-trials`
- `biomarker`
- `mixed-effects-models`
- `monte-carlo-simulation`
- `r-package`
- `reproducible-research`

### 5.2 Configure GitHub Actions

**Current workflows**:

- `.github/workflows/r-package.yml`
- `.github/workflows/render-report.yml`

**Action**: Verify workflows run successfully. Consider adding:

- Dependency caching
- Status badges to README.md

### 5.3 Add License Badge and Citation

Update README.md to include:

```markdown
<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/rgt47/pmsimstats2025/workflows/R-CMD-check/badge.svg)](https://github.com/rgt47/pmsimstats2025/actions)
<!-- badges: end -->
```

### 5.4 Create CITATION.cff

Add a `CITATION.cff` file for proper academic citation:

```yaml
cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
  - family-names: "Thomas"
    given-names: "Ronald G."
    orcid: "https://orcid.org/YOUR-ORCID"
title: "pmsimstats2025: N-of-1 Clinical Trial Simulation"
version: 0.1.0
date-released: 2025-01-28
url: "https://github.com/rgt47/pmsimstats2025"
```

---

## Priority 6: Security Review

### 6.1 Scan for Sensitive Information

Before making public, verify no sensitive data exists:

```bash
# Search for potential credentials
grep -r "password" --include="*.R" --include="*.Rmd"
grep -r "api_key" --include="*.R" --include="*.Rmd"
grep -r "secret" --include="*.R" --include="*.Rmd"

# Search for email addresses (review if appropriate)
grep -rE "[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}" --include="*.R"
```

### 6.2 Review .Rprofile

**File**: `.Rprofile`

Ensure no personal paths or credentials are included.

### 6.3 Check Git History

If sensitive information was ever committed, consider using `git filter-branch`
or BFG Repo-Cleaner before making public.

---

## Implementation Checklist

### Critical (Before Public Release)

- [x] Fix DESCRIPTION package name and author info *(completed 2025-01-28)*
- [x] Fix tests/testthat.R package reference *(completed 2025-01-28)*
- [x] Decide on CLAUDE.md inclusion (update .gitignore) *(completed 2025-01-28)*
- [x] Delete `claude.md` (lowercase duplicate) *(completed 2025-01-28)*
- [ ] Run `make check-renv` and fix any issues
- [ ] Run `make docker-test` and fix any failures
- [x] Security scan for credentials/sensitive data *(completed 2025-01-28 - PASSED)*

### Recommended (Before or Shortly After)

- [x] Delete temporary files (.zz*.R) *(completed 2025-01-28)*
- [x] Delete ZZCOLLAB_USER_GUIDE.md symlink *(completed 2025-01-28)*
- [x] Add CONTRIBUTING.md *(completed 2025-01-28)*
- [x] Add CITATION.cff *(completed 2025-01-28)*
- [x] Update README.md with badges *(completed 2025-01-28)*
- [ ] Verify GitHub Actions work

### Nice to Have (Can Be Done Later)

- [ ] Consolidate documentation structure
- [ ] Rename test_*.R exploratory scripts
- [x] Audit DESCRIPTION dependencies *(completed 2025-01-28 - see Priority 4)*
- [ ] Trim unused template packages (optional)
- [ ] Add repository topics on GitHub
- [ ] Add CODE_OF_CONDUCT.md

---

## Estimated Effort

| Priority | Tasks | Estimated Time |
|----------|-------|----------------|
| P1: Critical | 4 tasks | 30 minutes |
| P2: Cleanup | 5 tasks | 30 minutes |
| P3: Documentation | 4 tasks | 1-2 hours |
| P4: Dependencies | 3 tasks | 30 minutes |
| P5: GitHub Config | 4 tasks | 30 minutes |
| P6: Security | 3 tasks | 15 minutes |

**Total estimated time**: 3-4 hours for comprehensive cleanup

---

## Quick Start Script

Run these commands to perform the critical fixes:

```bash
# 1. Delete temporary files
rm -f .zz0098.R
rm -f analysis/scripts/.zzvim_r_temp.R
rm -f analysis/scripts/.zz8571.R
rm -f claude.md
rm -f ZZCOLLAB_USER_GUIDE.md

# 2. Update .gitignore (remove CLAUDE.md line, add temp patterns)
# Manual edit required

# 3. Fix DESCRIPTION (manual edit required)

# 4. Fix tests/testthat.R (manual edit required)

# 5. Validate
make check-renv
make docker-test

# 6. Security scan
grep -r "password\|api_key\|secret" --include="*.R" --include="*.Rmd"
```

---

*Document created: 2025-01-28*
*Repository: pmsimstats2025*
