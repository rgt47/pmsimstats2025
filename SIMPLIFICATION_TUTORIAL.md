# Git-Based Simplification Tutorial
## Creating a Learning Version of Complex Codebase

---

## Quick Start

```bash
cd ~/prj/d08/pmsimstats2025

# Create and switch to new branch
git checkout -b simplified-learning

# Verify you're on the new branch
git branch
```

You should see `* simplified-learning` indicating you're on the new branch.

---

## Simplification Plan

### Phase 1: Minimal Working Simulation (Start Here)

#### Goal
Single script that runs one simulation and prints results.

#### Keep Only:
1. **Single design type**: Hybrid only (simpler than crossover)
2. **Single parameter set**:
   - 70 participants
   - biomarker_correlation = 0.2
   - carryover_t1half = 0 (no carryover initially)
   - 5 iterations (instead of 20)
3. **Core functions only**:
   - `build_sigma_matrix()`
   - `generate_trial_data()`
   - `run_monte_carlo()`

#### Remove Temporarily:
- Crossover design logic
- Parameter grid expansion
- Sigma cache optimization
- Advanced diagnostics
- Validation functions
- Multiple carryover scenarios

#### Learning Goals:
- How does sigma matrix work?
- How is trial data generated?
- What does the mixed model do?

---

### Phase 2: Add Complexity Gradually

#### 2a. Add Parameter Variation
```r
# Test different biomarker correlations
biomarker_correlations <- c(0.2, 0.4)
```
**Learn**: How correlations affect power

#### 2b. Add Carryover
```r
# Introduce carryover effect
carryover_t1half <- 1.0
```
**Learn**: Mean vs correlation effects

#### 2c. Add Crossover Design
```r
# Implement second design
designs <- c("hybrid", "crossover")
```
**Learn**: Design comparison

#### 2d. Add Optimization
- Sigma cache
- Validation functions
- Diagnostic output

**Learn**: Edge cases and performance optimization

---

### Phase 3: Restore Full Complexity

Merge insights back to main branch (see Merging section below).

---

## Git Workflow

### Working on Simplified Branch

```bash
# Make changes to simplify code
# Edit files, remove complexity, add learning comments

# Commit your changes
git add -A
git commit -m "Phase 1: Minimal working simulation"

# Continue iterating
git add file.R
git commit -m "Added clearer comments to sigma matrix function"
```

### Comparing with Original

```bash
# See all differences between branches
git diff main simplified-learning

# See differences for specific file
git diff main simplified-learning -- analysis/scripts/pm_functions.R

# See just the file names that differ
git diff --name-only main simplified-learning
```

### Switching Between Branches

```bash
# Switch back to see original complex version
git checkout main

# Return to simplified version
git checkout simplified-learning

# Quick check: which branch am I on?
git branch
```

---

## Merging Insights Back to Main

### Option 1: Cherry-Pick Specific Commits (Recommended)

Use this when you make an improvement in `simplified-learning` that you want in `main`:

```bash
# 1. Make your improvement on simplified-learning branch
git checkout simplified-learning
# ... make changes ...
git add your_file.R
git commit -m "Improved function X with better clarity"

# 2. Note the commit hash
git log -1
# Copy the commit hash (looks like: abc1234def5678)

# 3. Switch to main
git checkout main

# 4. Cherry-pick just that commit
git cherry-pick abc1234def5678

# 5. Return to simplified branch
git checkout simplified-learning
```

**When to use**: You found a bug fix, better comment, or clearer implementation.

---

### Option 2: Selective File Merge

Use this when you want specific files but not the whole commit:

```bash
# 1. Switch to main branch
git checkout main

# 2. Get specific file from simplified-learning
git checkout simplified-learning -- analysis/scripts/pm_functions.R

# 3. Review what changed
git diff

# 4. Commit the changes
git commit -m "Merged improved comments from simplified branch"
```

**When to use**: You rewrote a single file to be clearer and want that version in main.

---

### Option 3: Full Merge (Use Carefully)

Use this only if you want ALL changes from simplified branch:

```bash
# Switch to main
git checkout main

# Merge everything from simplified-learning
git merge simplified-learning
```

**Warning**: This will try to merge all your simplifications into main. You'll likely need to resolve conflicts. Only use if:
- Simplified branch became the "better" version
- You're ready to replace complex version entirely

**More common approach**: Keep branches separate and only cherry-pick improvements.

---

### Option 4: Interactive Rebase (Advanced)

Use this to selectively apply multiple commits:

```bash
# 1. Create a temporary branch from main
git checkout main
git checkout -b temp-merge

# 2. Cherry-pick a range of commits
git log simplified-learning  # Find commit range
git cherry-pick abc123..def456

# 3. Review and merge to main
git checkout main
git merge temp-merge

# 4. Clean up
git branch -d temp-merge
```

**When to use**: You have multiple related commits to bring over.

---

## Typical Workflow Example

### Day 1: Setup
```bash
git checkout -b simplified-learning
# Remove 90% of complexity from main script
git commit -m "Phase 1: Stripped down to single simulation"
```

### Day 2: Understanding
```bash
# Add lots of learning comments
git commit -m "Added detailed comments explaining sigma matrix"
```

### Day 3: Discovery
```bash
# Found a bug in correlation calculation
git commit -m "Fixed: correlation clamping bug in edge case"

# This bug fix should go to main!
git log -1  # Copy commit hash
git checkout main
git cherry-pick <hash>
git checkout simplified-learning
```

### Day 4: Learning More
```bash
# Add parameter variation
git commit -m "Phase 2a: Added biomarker correlation variation"
```

### Day 5: Insight
```bash
# Realized a better way to structure validation
git commit -m "Improved validation function clarity"

# Want this in main too
git log -1
git checkout main
git cherry-pick <hash>
git checkout simplified-learning
```

---

## Best Practices

### 1. Commit Often on Simplified Branch
```bash
# Small, focused commits make cherry-picking easier
git commit -m "Simplified sigma matrix function"
git commit -m "Added learning comments to data generation"
git commit -m "Removed carryover complexity for now"
```

### 2. Keep Branches Independent
- **simplified-learning**: Learning version, lots of comments, reduced complexity
- **main**: Production version, full complexity, minimal comments

### 3. Cherry-Pick Improvements, Not Simplifications
**DO cherry-pick to main**:
- Bug fixes
- Performance improvements
- Better algorithms
- Useful comments that don't clutter

**DON'T cherry-pick to main**:
- Removed features
- Learning-only comments
- Simplified parameter grids

### 4. Document Your Learning
```bash
# Add notes as you go
echo "## Day 1 Learning" >> LEARNING_NOTES.md
echo "- Sigma matrix is built from correlation parameters" >> LEARNING_NOTES.md
git add LEARNING_NOTES.md
git commit -m "Added learning notes"
```

---

## Example Simplifications for This Project

### 1. Simplify Main Script

**Original** (line 628):
```r
expand_grid(
  design = c("hybrid", "crossover"),
  biomarker_correlation = c(0.2, 0.4),
  carryover_t1half = c(0, 1.0, 2.0)
)
```

**Simplified**:
```r
# Just one scenario to understand the basics
params <- list(
  design = "hybrid",
  biomarker_correlation = 0.2,
  carryover_t1half = 0
)
```

### 2. Simplify Functions

**Original** `build_sigma_matrix()`: 200+ lines with validation

**Simplified**:
```r
# Remove all validation for now, just build the matrix
# Add comments explaining each correlation parameter
build_sigma_matrix_simple <- function(...) {
  # Focus on understanding the structure
  # Add print statements to see intermediate values
}
```

### 3. Add Learning Comments

```r
# LEARNING: The sigma matrix represents correlations between:
# - Time-variant measurements (rows 1-4)
# - Biomarker measurements (rows 5-8)
# - Response measurements (rows 9-12)
# Each set of 4 represents the 4 time points in the trial
```

---

## Troubleshooting

### I Want to Start Over
```bash
# Reset simplified branch to match main
git checkout simplified-learning
git reset --hard main

# Or delete and recreate
git checkout main
git branch -D simplified-learning
git checkout -b simplified-learning
```

### I Accidentally Committed to Main
```bash
# Move commit to simplified-learning
git log -1  # Get commit hash
git checkout simplified-learning
git cherry-pick <hash>
git checkout main
git reset --hard HEAD~1  # Remove last commit from main
```

### Merge Conflict During Cherry-Pick
```bash
# Git will pause and show conflicts
# Edit the files to resolve conflicts
git add <resolved-files>
git cherry-pick --continue

# Or abort
git cherry-pick --abort
```

---

## Next Steps

1. **Create branch**: `git checkout -b simplified-learning`
2. **Simplify**: Remove complexity, add comments
3. **Commit**: `git commit -m "Phase 1 complete"`
4. **Learn**: Run simplified version, understand it
5. **Improve**: Find better ways to do things
6. **Cherry-pick**: Move improvements back to main
7. **Repeat**: Gradually add complexity back

---

## Resources

### Git Commands Quick Reference
```bash
# Branch management
git branch                          # List branches
git checkout <branch>               # Switch branch
git checkout -b <branch>            # Create and switch

# Comparing
git diff main simplified-learning   # See all differences
git log main..simplified-learning   # See commits only in simplified

# Merging
git cherry-pick <hash>              # Copy one commit
git merge <branch>                  # Merge all changes

# Undoing
git reset --hard HEAD~1             # Remove last commit
git checkout -- <file>              # Discard changes to file
```

### Visualization
```bash
# Install gitk (if not already installed)
# On macOS: brew install git-gui

# Visualize branch history
gitk --all
```

---

## Project-Specific Tips

### Start With These Files:
1. `analysis/scripts/pm_functions.R` - Core functions
2. `analysis/scripts/full_pmsim_analysis_hyb_versus_co.R` - Main script

### Key Concepts to Understand:
1. **Sigma matrix structure**: How correlation parameters create the covariance matrix
2. **Carryover mechanism**: How it affects means but not correlations
3. **Mixed model**: What the `lmer()` formula means
4. **Design differences**: Hybrid vs Crossover randomization

### Simplification Targets:
- Lines 511-649 in main script: Parameter grid → single params
- `build_sigma_matrix()`: Remove validation, add print statements
- `run_monte_carlo()`: Reduce iterations from 20 to 5
- Remove crossover design entirely at first

---

*Created: 2025-11-19*
*For project: pmsimstats2025*
