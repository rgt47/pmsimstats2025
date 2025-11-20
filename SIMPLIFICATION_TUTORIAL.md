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

---
---

# PART 2: Step-by-Step Walkthrough with Narrative and Intuition

---

## Step 1: Understanding Why We Use Branches for Learning

### The Narrative

Imagine you have a beautifully engineered Swiss watch (your complex codebase). It works perfectly, but you can't understand how all the gears interact. You want to learn, but you're afraid to take it apart because:

1. You might not be able to put it back together
2. You might lose pieces
3. You need it to keep working while you learn

**Git branches are like having a perfect photograph of the watch that you can restore anytime.**

### The Mental Model

Think of your repository as a tree:

```
main branch (complex, working code)
     |
     |---- simplified-learning branch (your learning playground)
```

- **Main branch**: The "production" version. Complex, complete, proven to work.
- **Simplified-learning branch**: Your personal laboratory. You can break things, simplify drastically, add tons of explanatory comments, and experiment freely.

### The Key Insight

**The magic**: At any moment, you can:
- Look back at the original (switch to `main`)
- Compare your changes side-by-side (`git diff`)
- Copy specific improvements back to the original (cherry-pick)
- Throw away your experiments and start over (reset)

**The safety**: Your complex working code on `main` never changes unless you explicitly choose to update it.

### Why This Beats Copying to a New Repo

Copying creates two completely separate universes with no easy way to:
- Compare what you changed
- Move improvements between versions
- See the history of what you learned
- Switch back and forth instantly

Branches keep them connected but separate - like parallel dimensions you can hop between.

---

## Step 2: Creating Your Learning Branch

### The Narrative

You're about to create a "parallel universe" of your code. Think of this moment as Neo in The Matrix being offered the red pill - you're about to see how your code really works, but you can always come back.

### The Command Breakdown

```bash
cd ~/prj/d08/pmsimstats2025
git checkout -b simplified-learning
```

Let's decode what's actually happening:

#### `git checkout -b simplified-learning`

This is actually TWO commands combined:

```bash
# What it really does:
git branch simplified-learning    # Creates new branch
git checkout simplified-learning  # Switches to it
```

**The `-b` flag** means "create **b**efore checking out"

### The Intuition: What Just Happened?

#### Before the command:
```
main branch (you are here)
  ├─ file1.R
  ├─ file2.R
  └─ file3.R
```

#### After the command:
```
main branch
  ├─ file1.R (original)
  ├─ file2.R (original)
  └─ file3.R (original)

simplified-learning branch (you are here now)
  ├─ file1.R (identical copy for now)
  ├─ file2.R (identical copy for now)
  └─ file3.R (identical copy for now)
```

**Key insight**: Right now, both branches are IDENTICAL. They're pointing to the exact same commit in history. They only diverge when you make changes.

### The Mental Model: Branches Are Pointers

Think of branches as **bookmarks** pointing to specific moments in your code's history:

```
History Timeline:
[commit1] → [commit2] → [commit3] ← main branch (bookmark)
                                  ← simplified-learning (bookmark)
```

Both bookmarks point to the same place right now. When you make a change on `simplified-learning`, it moves forward while `main` stays put:

```
After you make changes:
[commit1] → [commit2] → [commit3] ← main (unchanged)
                            ↓
                        [commit4] ← simplified-learning (moved forward)
```

### Verifying You're on the New Branch

```bash
git branch
```

You should see:
```
  main
* simplified-learning
```

**The asterisk (`*`)** is like a "You Are Here" marker on a map. It shows which branch you're currently on.

### What You Can Do Now

At this moment:
- ✅ All your files look identical to `main`
- ✅ You can start editing freely
- ✅ Your changes won't affect `main`
- ✅ You can switch back to `main` anytime with `git checkout main`

### The Safety Net

**Important realization**: Nothing you do on this branch can hurt `main` unless you explicitly merge or cherry-pick changes. You could:
- Delete entire files
- Completely rewrite functions
- Break everything
- Experiment wildly

...and `main` remains untouched, like it's in a time capsule.

### Try This Right Now

Let's prove the branches are independent:

```bash
# See that files are currently identical
ls -la

# This is your playground now. You're standing on simplified-learning.
# main is still there, frozen in time, waiting for you.
```

---

## Step 5: Comparing Branches with git diff

### The Narrative

Imagine you're a detective with two photographs of the same room taken at different times. Your job is to spot every difference - what was moved, what was added, what was removed.

`git diff` is your magnifying glass. It shows you **exactly** what changed between any two versions of your code.

### The Mental Model: Two Snapshots Side by Side

```
main branch                    simplified-learning branch
┌─────────────────┐            ┌─────────────────┐
│ c(0, 0.3, 0.48) │  ───────►  │ c(0.3)          │
│ c(0, 0.5, 1.0)  │  changed   │ c(0)            │
│ n_iterations=20 │  ───────►  │ n_iterations=5  │
│ hybrid+crossover│            │ hybrid only     │
└─────────────────┘            └─────────────────┘
```

### Essential Diff Commands

#### See summary of all differences:
```bash
git diff main simplified-learning --stat
```

#### See actual code changes:
```bash
git diff main simplified-learning -- path/to/file.R
```

#### See just file names that changed:
```bash
git diff --name-only main simplified-learning
```

### Reading Diff Output

| Symbol | Meaning |
|--------|---------|
| `-` (red) | Line **removed** from main |
| `+` (green) | Line **added** in simplified-learning |
| `@@` | Location marker (line numbers) |

### Example Output Explained

```diff
-  biomarker_correlation = c(0, 0.3, 0.48),    ← REMOVED (3 values)
+  biomarker_correlation = c(0.3),              ← ADDED (1 value)
```

This shows you replaced 3 parameter values with 1 - exactly what you intended for simplification.

### The Key Insight

**You can always see exactly what you changed.** This means:
- No fear of breaking things (you can see what to undo)
- Clear documentation of your learning journey
- Easy to explain changes to others

---
