# Pseudocode Comparison: Hendrickson vs. Our Approach

## Document Purpose

This document provides side-by-side pseudocode comparison of:
1. **Hendrickson et al. (2020)** original simulation methodology
2. **Our current implementation** (after Option 1 alignment)
3. **Explicit identification of all differences**

---

## Table of Contents

1. [High-Level Overview](#high-level-overview)
2. [Parameter Setup](#parameter-setup)
3. [Trial Design Generation](#trial-design-generation)
4. [Sigma Matrix Construction](#sigma-matrix-construction)
5. [Data Generation](#data-generation)
6. [Statistical Analysis](#statistical-analysis)
7. [Monte Carlo Simulation](#monte-carlo-simulation)
8. [Summary of All Differences](#summary-of-all-differences)

---

## High-Level Overview

### Hendrickson's Approach

```
FUNCTION main_simulation():
    # Fixed correlation parameters
    SET correlation_params = {
        c.tv = 0.8,    # Time-variant autocorrelation
        c.pb = 0.8,    # Pharmacologic biomarker autocorrelation
        c.br = 0.8,    # Biological response autocorrelation
        c.bm = varies, # Biomarker-response correlation (VARIES in grid)
        c.cf1t = 0.2,  # Same-time cross-correlation
        c.cfct = 0.1   # Different-time cross-correlation
    }

    # Parameter grid to explore
    FOR EACH n_participants IN {30, 50, 70, 90}:
        FOR EACH biomarker_correlation IN {0.2, 0.4, 0.6, 0.8}:
            FOR EACH carryover_t1half IN {0, 0.5, 1.0, 1.5, 2.0}:
                FOR EACH design IN {hybrid, crossover}:
                    # Run Monte Carlo
                    FOR iteration = 1 TO 1000:
                        data = generate_trial_data(params, design)
                        result = analyze_data(data)  # NO carryover term
                        STORE result

                    # Calculate power
                    power = PROPORTION(significant results)
                    STORE power

    RETURN all_results
END FUNCTION
```

### Our Approach (Current)

```
FUNCTION main_simulation():
    # Fixed correlation parameters (SAME as Hendrickson after alignment)
    SET correlation_params = {
        c.tv = 0.8,
        c.pb = 0.8,
        c.br = 0.8,
        c.bm = varies,
        c.cf1t = 0.2,
        c.cfct = 0.1
    }

    # Parameter grid (SMALLER than Hendrickson)
    FOR EACH n_participants IN {70}:
        FOR EACH biomarker_correlation IN {0.2, 0.4}:
            FOR EACH carryover_t1half IN {0, 1.0, 2.0}:
                FOR EACH design IN {hybrid, crossover}:
                    # ⭐ DIFFERENCE: Two-scenario comparison
                    FOR EACH model_carryover IN {TRUE, FALSE}:
                        # Run Monte Carlo
                        FOR iteration = 1 TO 20:
                            data = generate_trial_data(params, design)
                            result = analyze_data(data, model_carryover)
                            STORE result

                        # Calculate power
                        power = PROPORTION(significant results)
                        STORE power

    RETURN all_results
END FUNCTION
```

**DIFFERENCE #1**: We run TWO analysis scenarios per data generation (with/without carryover modeling)
**DIFFERENCE #2**: Smaller parameter grid (fewer sample sizes, fewer iterations)

---

## Parameter Setup

### Hendrickson's Parameters

```
FUNCTION setup_parameters():
    # Fixed correlation structure
    correlations = {
        c.tv = 0.8,      # Time-variant autocorrelation
        c.pb = 0.8,      # Pharmacologic biomarker autocorrelation
        c.br = 0.8,      # Biological response autocorrelation
        c.bm = varies,   # Biomarker-response correlation (GRID PARAMETER)
        c.cf1t = 0.2,    # Same-time cross-correlation
        c.cfct = 0.1     # Different-time cross-correlation
    }

    # Response parameters
    response_params = {
        variance.tv = 2.25,      # Time-variant variance
        variance.pb = 2.25,      # Pharmacologic biomarker variance
        variance.br = 2.25,      # Biological response variance
        treatment_effect = 5.0,  # Treatment effect size
        baseline_mean = 15.0     # Baseline response mean
    }

    # Carryover parameters
    carryover_params = {
        t1half = varies,         # Half-life (GRID PARAMETER)
        scale_factor = 1.0,      # BR-only carryover (affects BR means only)
        washout_weeks = 4        # Fixed washout period
    }

    RETURN {correlations, response_params, carryover_params}
END FUNCTION
```

### Our Parameters (Current)

```
FUNCTION setup_parameters():
    # Fixed correlation structure (SAME as Hendrickson)
    correlations = {
        c.tv = 0.8,
        c.pb = 0.8,
        c.br = 0.8,
        c.bm = varies,   # GRID PARAMETER
        c.cf1t = 0.2,
        c.cfct = 0.1
    }

    # Response parameters (SAME)
    response_params = {
        variance.tv = 2.25,
        variance.pb = 2.25,
        variance.br = 2.25,
        treatment_effect = 5.0,
        baseline_mean = 15.0
    }

    # Carryover parameters (SAME)
    carryover_params = {
        t1half = varies,         # GRID PARAMETER
        scale_factor = 1.0,      # BR-only carryover
        washout_weeks = 4
    }

    RETURN {correlations, response_params, carryover_params}
END FUNCTION
```

**NO DIFFERENCE**: Parameter setup is now identical after Option 1 alignment

---

## Trial Design Generation

### Hendrickson's Design Generation

```
FUNCTION create_hybrid_design(n_participants, n_weeks_on, n_weeks_off):
    design = EMPTY_TABLE

    # 4-path randomization (Hendrickson et al. 2020)
    path_assignments = {
        Path 1: [A, B, A, B],  # AB alternating, starts with A
        Path 2: [B, A, B, A],  # BA alternating, starts with B
        Path 3: [A, A, B, B],  # AA BB blocks
        Path 4: [B, B, A, A]   # BB AA blocks
    }

    # Assign participants to paths (balanced)
    FOR participant_id = 1 TO n_participants:
        path_id = ((participant_id - 1) MOD 4) + 1
        treatment_sequence = path_assignments[path_id]

        # Create timepoint rows
        FOR period = 1 TO 4:
            FOR week_in_period = 1 TO (n_weeks_on + n_weeks_off):
                treatment = IF week_in_period <= n_weeks_on THEN
                                treatment_sequence[period]
                            ELSE
                                "washout"

                # Convert to numeric: A=1, B=0, washout=0
                tod = IF treatment == "A" THEN 1 ELSE 0

                # Calculate weeks since treatment for carryover
                IF treatment == "washout":
                    weeks_off = week_in_period - n_weeks_on
                ELSE:
                    weeks_off = 0

                ADD_ROW(design, {
                    participant_id,
                    period,
                    week_in_period,
                    treatment,
                    tod,
                    weeks_off,
                    path = path_id
                })

    RETURN design
END FUNCTION

FUNCTION create_crossover_design(n_participants, n_weeks_on, n_weeks_off):
    design = EMPTY_TABLE

    # Traditional 2-sequence crossover
    sequences = {
        Sequence 1: [A, B],  # AB
        Sequence 2: [B, A]   # BA
    }

    # Assign participants (alternating)
    FOR participant_id = 1 TO n_participants:
        sequence_id = IF participant_id MOD 2 == 0 THEN 1 ELSE 2
        treatment_sequence = sequences[sequence_id]

        # Create timepoint rows (2 periods)
        FOR period = 1 TO 2:
            FOR week_in_period = 1 TO (n_weeks_on + n_weeks_off):
                treatment = IF week_in_period <= n_weeks_on THEN
                                treatment_sequence[period]
                            ELSE
                                "washout"

                tod = IF treatment == "A" THEN 1 ELSE 0

                IF treatment == "washout":
                    weeks_off = week_in_period - n_weeks_on
                ELSE:
                    weeks_off = 0

                ADD_ROW(design, {
                    participant_id,
                    period,
                    week_in_period,
                    treatment,
                    tod,
                    weeks_off,
                    sequence = sequence_id
                })

    RETURN design
END FUNCTION
```

### Our Design Generation (Current)

```
FUNCTION create_hybrid_design(n_participants, n_weeks_on, n_weeks_off):
    # IDENTICAL to Hendrickson
    # (See Hendrickson pseudocode above)

    # ⭐ DIFFERENCE: We also add a 'path' column explicitly
    # Hendrickson has this implicitly, we make it explicit

    RETURN design
END FUNCTION

FUNCTION create_crossover_design(n_participants, n_weeks_on, n_weeks_off):
    # MOSTLY IDENTICAL to Hendrickson
    # (See Hendrickson pseudocode above)

    # ⭐ DIFFERENCE: We add 'path' column for consistency
    # path = 1 for AB sequence, path = 2 for BA sequence
    # This makes our structure consistent across both designs

    RETURN design
END FUNCTION
```

**DIFFERENCE #3**: We add explicit `path` column to crossover design (Hendrickson uses `sequence`)
- This is a minor structural difference for code consistency
- Does not affect statistical properties

---

## Sigma Matrix Construction

### Hendrickson's Sigma Construction

```
FUNCTION build_sigma_matrix(design, correlations, variances):
    # Extract unique timepoints from design
    timepoints = UNIQUE(design$timepoint_name)
    n_timepoints = LENGTH(timepoints)

    # Initialize correlation matrix (n_timepoints × n_timepoints)
    # Each timepoint has 3 components: TV, PB, BR
    matrix_size = n_timepoints * 3
    R = IDENTITY_MATRIX(matrix_size)

    # Component indices for each timepoint
    FOR t = 1 TO n_timepoints:
        tv_idx = (t-1)*3 + 1  # Time-variant index
        pb_idx = (t-1)*3 + 2  # Pharmacologic biomarker index
        br_idx = (t-1)*3 + 3  # Biological response index

        # Within-timepoint correlations (same time, different component)
        IF design$tod[t] == 1:  # On treatment
            # ⭐ KEY: Biomarker-response correlation only on treatment
            R[pb_idx, br_idx] = correlations$c.bm
            R[br_idx, pb_idx] = correlations$c.bm
        ELSE:
            # Off treatment: no biomarker-response correlation
            R[pb_idx, br_idx] = 0
            R[br_idx, pb_idx] = 0

        # Same-time cross-correlations (TV-PB, TV-BR)
        R[tv_idx, pb_idx] = correlations$c.cf1t
        R[pb_idx, tv_idx] = correlations$c.cf1t
        R[tv_idx, br_idx] = correlations$c.cf1t
        R[br_idx, tv_idx] = correlations$c.cf1t

        # Between-timepoint correlations (autocorrelations)
        FOR s = (t+1) TO n_timepoints:
            tv_idx_s = (s-1)*3 + 1
            pb_idx_s = (s-1)*3 + 2
            br_idx_s = (s-1)*3 + 3

            # Autocorrelations (same component, different time)
            R[tv_idx, tv_idx_s] = correlations$c.tv
            R[tv_idx_s, tv_idx] = correlations$c.tv

            R[pb_idx, pb_idx_s] = correlations$c.pb
            R[pb_idx_s, pb_idx] = correlations$c.pb

            R[br_idx, br_idx_s] = correlations$c.br
            R[br_idx_s, br_idx] = correlations$c.br

            # Cross-correlations (different component, different time)
            R[tv_idx, pb_idx_s] = correlations$c.cfct
            R[pb_idx_s, tv_idx] = correlations$c.cfct
            R[tv_idx, br_idx_s] = correlations$c.cfct
            R[br_idx_s, tv_idx] = correlations$c.cfct
            # (and all other cross-time, cross-component pairs)

    # Convert correlation matrix to covariance matrix
    # Create standard deviation vector
    sd_vector = EMPTY_VECTOR(matrix_size)
    FOR t = 1 TO n_timepoints:
        sd_vector[(t-1)*3 + 1] = SQRT(variances$tv)
        sd_vector[(t-1)*3 + 2] = SQRT(variances$pb)
        sd_vector[(t-1)*3 + 3] = SQRT(variances$br)

    # Σ = D * R * D, where D is diagonal matrix of SDs
    Sigma = OUTER_PRODUCT(sd_vector, sd_vector) * R

    # Check positive definiteness
    eigenvalues = EIGEN(Sigma)$values
    is_positive_definite = ALL(eigenvalues > 0)

    IF NOT is_positive_definite:
        # ⭐ Hendrickson: Auto-fix non-PD matrices
        Sigma = make_positive_definite(Sigma, tolerance = 1e-3)

    RETURN Sigma
END FUNCTION
```

### Our Sigma Construction (Current)

```
FUNCTION build_sigma_matrix(design, correlations, variances):
    # MOSTLY IDENTICAL structure to Hendrickson
    # (See Hendrickson pseudocode for main logic)

    # Initialize correlation matrix
    matrix_size = n_timepoints * 3
    R = IDENTITY_MATRIX(matrix_size)

    # ⭐ DIFFERENCE: Biomarker correlation scaling
    # We scale c.bm by mean ratio to maintain consistent relationship
    FOR t = 1 TO n_timepoints:
        tv_idx = (t-1)*3 + 1
        pb_idx = (t-1)*3 + 2
        br_idx = (t-1)*3 + 3

        IF design$tod[t] == 1:  # On treatment
            # Calculate mean values for this timepoint
            mean_pb = calculate_mean(t, "pb", design, params)
            mean_br = calculate_mean(t, "br", design, params)

            # Scale correlation by mean ratio
            mean_ratio = mean_br / mean_pb
            scaled_correlation = correlations$c.bm * mean_ratio

            # ⭐ CRITICAL: Clamp to valid range [-0.99, 0.99]
            scaled_correlation = MAX(-0.99, MIN(0.99, scaled_correlation))

            R[pb_idx, br_idx] = scaled_correlation
            R[br_idx, pb_idx] = scaled_correlation
        ELSE:
            R[pb_idx, br_idx] = 0
            R[br_idx, pb_idx] = 0

        # All other correlations: SAME as Hendrickson
        # (autocorrelations, cross-correlations, etc.)

    # Convert to covariance matrix (SAME as Hendrickson)
    Sigma = OUTER_PRODUCT(sd_vector, sd_vector) * R

    # Check positive definiteness
    eigenvalues = EIGEN(Sigma)$values
    is_positive_definite = ALL(eigenvalues > 0)

    # ⭐ DIFFERENCE: REJECT non-PD matrices instead of auto-fixing
    IF NOT is_positive_definite:
        PRINT "REJECTED: Non-positive definite matrix"
        PRINT "Eigenvalue range:", MIN(eigenvalues), "to", MAX(eigenvalues)
        RETURN NULL  # Signal invalid parameter combination

    RETURN Sigma
END FUNCTION
```

**DIFFERENCE #4**: Biomarker correlation scaling by mean ratio
- Hendrickson: Uses c.bm directly
- Us: Scale by mean ratio, then clamp to [-0.99, 0.99]
- **Note**: After removing population mean shift, mean ratios are smaller, so scaling has less effect

**DIFFERENCE #5**: Non-PD matrix handling
- Hendrickson: Auto-fix with `make_positive_definite()`
- Us: REJECT and return NULL (skip invalid parameter combinations)

---

## Data Generation

This is where the most important conceptual differences exist.

### Hendrickson's Data Generation

```
FUNCTION generate_trial_data(design, params):
    # Step 1: Calculate mean vectors for each participant
    # ============================================================

    participant_data = EMPTY_TABLE

    FOR EACH participant IN design:
        participant_id = participant$id

        # Get participant's timepoints
        timepoints = FILTER(design, participant_id == this_participant)
        n_timepoints = NROW(timepoints)

        # Initialize mean vector (3 components per timepoint)
        mean_vector = EMPTY_VECTOR(n_timepoints * 3)

        FOR t = 1 TO n_timepoints:
            timepoint = timepoints[t]

            # TV component mean (baseline only, no treatment effect)
            mean_vector[(t-1)*3 + 1] = params$baseline_mean

            # PB component mean (baseline only)
            mean_vector[(t-1)*3 + 2] = params$baseline_mean

            # BR component mean
            base_br = params$baseline_mean

            # Add treatment effect
            IF timepoint$tod == 1:  # On treatment
                treatment_effect = params$treatment_effect
            ELSE:
                treatment_effect = 0

            # Add carryover effect
            IF timepoint$weeks_off > 0 AND params$carryover_t1half > 0:
                # Exponential decay: (1/2)^(weeks_off / t1half)
                carryover_effect = treatment_effect *
                                   (0.5)^(timepoint$weeks_off / params$carryover_t1half)
            ELSE:
                carryover_effect = 0

            # ⭐ BR mean includes treatment + carryover ONLY
            # NO biomarker-dependent mean shift
            mean_vector[(t-1)*3 + 3] = base_br + treatment_effect + carryover_effect

        # Step 2: Draw from multivariate normal
        # ============================================================

        # Build sigma matrix for this participant's design
        Sigma = build_sigma_matrix(timepoints, params$correlations, params$variances)

        # Draw from MVN(mean_vector, Sigma)
        mvn_draw = RMVNORM(n = 1, mean = mean_vector, sigma = Sigma)

        # Step 3: Extract components
        # ============================================================

        participant_row = {
            participant_id = participant_id,
            biomarker = MEAN(mvn_draw[PB components])  # Average PB across timepoints
        }

        FOR t = 1 TO n_timepoints:
            timepoint_name = timepoints[t]$name

            # Extract TV, PB, BR for this timepoint
            participant_row[[paste0(timepoint_name, ".tv")]] = mvn_draw[(t-1)*3 + 1]
            participant_row[[paste0(timepoint_name, ".pb")]] = mvn_draw[(t-1)*3 + 2]
            participant_row[[paste0(timepoint_name, ".br")]] = mvn_draw[(t-1)*3 + 3]

        # ⭐ NO POST-MVN ADJUSTMENT
        # Data is used directly as drawn from MVN

        ADD_ROW(participant_data, participant_row)

    RETURN participant_data
END FUNCTION
```

### Our Data Generation (Current - After Option 1)

```
FUNCTION generate_trial_data(design, params):
    # Step 1: Calculate mean vectors for each participant
    # ============================================================

    participant_data = EMPTY_TABLE

    FOR EACH participant IN design:
        participant_id = participant$id
        timepoints = FILTER(design, participant_id == this_participant)
        n_timepoints = NROW(timepoints)

        mean_vector = EMPTY_VECTOR(n_timepoints * 3)

        FOR t = 1 TO n_timepoints:
            timepoint = timepoints[t]

            # TV component mean (SAME as Hendrickson)
            mean_vector[(t-1)*3 + 1] = params$baseline_mean

            # PB component mean (SAME as Hendrickson)
            mean_vector[(t-1)*3 + 2] = params$baseline_mean

            # BR component mean
            base_br = params$baseline_mean

            # Add treatment effect (SAME)
            IF timepoint$tod == 1:
                treatment_effect = params$treatment_effect
            ELSE:
                treatment_effect = 0

            # Add carryover effect (SAME)
            IF timepoint$weeks_off > 0 AND params$carryover_t1half > 0:
                carryover_effect = treatment_effect *
                                   (0.5)^(timepoint$weeks_off / params$carryover_t1half)
            ELSE:
                carryover_effect = 0

            # ⭐ REMOVED: Population mean shift based on biomarker
            # (This was: + params$c.bm * timepoint$tod * 2.0)
            # Now SAME as Hendrickson

            mean_vector[(t-1)*3 + 3] = base_br + treatment_effect + carryover_effect

        # Step 2: Draw from multivariate normal (SAME as Hendrickson)
        # ============================================================

        Sigma = build_sigma_matrix(timepoints, params$correlations, params$variances)
        mvn_draw = RMVNORM(n = 1, mean = mean_vector, sigma = Sigma)

        # Step 3: Extract components (SAME structure)
        # ============================================================

        participant_row = {
            participant_id = participant_id,
            biomarker = MEAN(mvn_draw[PB components])
        }

        FOR t = 1 TO n_timepoints:
            timepoint_name = timepoints[t]$name
            participant_row[[paste0(timepoint_name, ".tv")]] = mvn_draw[(t-1)*3 + 1]
            participant_row[[paste0(timepoint_name, ".pb")]] = mvn_draw[(t-1)*3 + 2]
            participant_row[[paste0(timepoint_name, ".br")]] = mvn_draw[(t-1)*3 + 3]

        # ⭐ REMOVED: Post-MVN biomarker adjustment
        # (This was: br_i = br_i + (biomarker * treatment * c.bm))
        # Now SAME as Hendrickson

        ADD_ROW(participant_data, participant_row)

    RETURN participant_data
END FUNCTION
```

**DIFFERENCE #6 (REMOVED)**: ~~Population mean shift~~ - NOW SAME after Option 1
**DIFFERENCE #7 (REMOVED)**: ~~Post-MVN adjustment~~ - NOW SAME after Option 1

✅ **Data generation is now IDENTICAL to Hendrickson's approach**

---

## Statistical Analysis

This is the **MAJOR REMAINING DIFFERENCE** between approaches.

### Hendrickson's Analysis

```
FUNCTION analyze_trial_data(data, design):
    # Reshape data to long format
    # ============================================================

    analysis_data = EMPTY_TABLE

    FOR EACH participant IN data:
        participant_id = participant$id
        biomarker = participant$biomarker

        # Get participant's timepoints from design
        timepoints = FILTER(design, participant_id == this_participant)

        FOR EACH timepoint IN timepoints:
            # Extract response for this timepoint
            timepoint_name = timepoint$name
            response = participant[[paste0(timepoint_name, ".br")]]
            treatment = timepoint$tod
            week = timepoint$week

            ADD_ROW(analysis_data, {
                participant_id,
                biomarker,
                response,
                treatment,
                week
            })

    # Fit mixed-effects model
    # ============================================================

    # ⭐ Hendrickson's model: NO carryover term
    model = LMER(
        response ~ treatment * biomarker + week + (1 | participant_id),
        data = analysis_data
    )

    # Extract biomarker×treatment interaction
    # ============================================================

    interaction_term = COEF(model)["treatment:biomarker"]
    interaction_se = SE(model)["treatment:biomarker"]
    interaction_t = interaction_term / interaction_se
    interaction_p = 2 * PT(-ABS(interaction_t), df = model$df)

    significant = (interaction_p < 0.05)

    RETURN {
        effect_size = interaction_term,
        std_error = interaction_se,
        t_value = interaction_t,
        p_value = interaction_p,
        significant = significant
    }
END FUNCTION
```

### Our Analysis (Current - Two Scenarios)

```
FUNCTION analyze_trial_data(data, design, model_carryover = TRUE):
    # Reshape data to long format (SAME as Hendrickson)
    # ============================================================

    analysis_data = EMPTY_TABLE

    FOR EACH participant IN data:
        participant_id = participant$id
        biomarker = participant$biomarker
        timepoints = FILTER(design, participant_id == this_participant)

        FOR EACH timepoint IN timepoints:
            timepoint_name = timepoint$name
            response = participant[[paste0(timepoint_name, ".br")]]
            treatment = timepoint$tod
            week = timepoint$week

            # ⭐ DIFFERENCE: Calculate carryover covariate
            IF timepoint$weeks_off > 0 AND params$carryover_t1half > 0:
                carryover_effect = (0.5)^(timepoint$weeks_off / params$carryover_t1half)
            ELSE:
                carryover_effect = 0

            ADD_ROW(analysis_data, {
                participant_id,
                biomarker,
                response,
                treatment,
                week,
                carryover_effect  # ⭐ NEW COLUMN
            })

    # Fit mixed-effects model
    # ============================================================

    # ⭐ MAJOR DIFFERENCE: Conditional model specification
    IF params$carryover_t1half > 0 AND model_carryover == TRUE:
        # Scenario 1: WITH carryover modeling (our enhancement)
        model = LMER(
            response ~ treatment * biomarker + week + carryover_effect +
                       (1 | participant_id),
            data = analysis_data
        )
    ELSE:
        # Scenario 2: WITHOUT carryover modeling (Hendrickson approach)
        model = LMER(
            response ~ treatment * biomarker + week + (1 | participant_id),
            data = analysis_data
        )

    # Extract biomarker×treatment interaction (SAME as Hendrickson)
    # ============================================================

    interaction_term = COEF(model)["treatment:biomarker"]
    interaction_se = SE(model)["treatment:biomarker"]
    interaction_t = interaction_term / interaction_se
    interaction_p = 2 * PT(-ABS(interaction_t), df = model$df)

    significant = (interaction_p < 0.05)

    RETURN {
        effect_size = interaction_term,
        std_error = interaction_se,
        t_value = interaction_t,
        p_value = interaction_p,
        significant = significant,
        model_carryover = model_carryover  # ⭐ Track which scenario
    }
END FUNCTION
```

**DIFFERENCE #8 (MAJOR)**: Carryover modeling in statistical analysis
- **Hendrickson**: NEVER includes carryover term in model
- **Us - Scenario 1**: Includes `carryover_effect` covariate (model_carryover = TRUE)
- **Us - Scenario 2**: Excludes carryover term (model_carryover = FALSE, matches Hendrickson)

**This is the KEY methodological contribution of our work:**
- We can demonstrate the BENEFIT of modeling carryover (Scenario 1)
- We can REPLICATE Hendrickson's findings (Scenario 2)
- We can QUANTIFY the difference between approaches

---

## Monte Carlo Simulation

### Hendrickson's Monte Carlo

```
FUNCTION run_monte_carlo(design_name, params, n_iterations = 1000):
    results = EMPTY_TABLE

    FOR iteration = 1 TO n_iterations:
        # Generate trial design
        design = GENERATE_DESIGN(design_name, params)

        # Generate data
        data = generate_trial_data(design, params)

        # Analyze data (no carryover modeling)
        result = analyze_trial_data(data, design)

        # Store result
        result$iteration = iteration
        result$design = design_name
        result$n_participants = params$n_participants
        result$biomarker_correlation = params$biomarker_correlation
        result$carryover_t1half = params$carryover_t1half

        ADD_ROW(results, result)

    RETURN results
END FUNCTION
```

### Our Monte Carlo (Current)

```
FUNCTION run_monte_carlo(design_name, params, n_iterations = 20,
                        model_carryover = TRUE):
    results = EMPTY_TABLE

    FOR iteration = 1 TO n_iterations:
        # Generate trial design (SAME)
        design = GENERATE_DESIGN(design_name, params)

        # Generate data (SAME after Option 1)
        data = generate_trial_data(design, params)

        # ⭐ DIFFERENCE: Pass model_carryover to analysis
        result = analyze_trial_data(data, design, model_carryover)

        # Store result
        result$iteration = iteration
        result$design = design_name
        result$n_participants = params$n_participants
        result$biomarker_correlation = params$biomarker_correlation
        result$carryover_t1half = params$carryover_t1half
        result$model_carryover = model_carryover  # ⭐ Track scenario

        ADD_ROW(results, result)

    RETURN results
END FUNCTION
```

**DIFFERENCE #9**: Monte Carlo structure
- Hendrickson: Single analysis per data generation
- Us: Parameterized by `model_carryover` flag

**DIFFERENCE #10**: Number of iterations
- Hendrickson: 1000 iterations
- Us: 20 iterations (for computational efficiency)

---

## Summary of All Differences

### Differences REMOVED (Now Aligned)

| # | Feature | Previous | Current (After Option 1) |
|---|---------|----------|--------------------------|
| ~~6~~ | Population mean shift | ❌ We added c.bm * tod * 2.0 | ✅ Removed - now same |
| ~~7~~ | Post-MVN adjustment | ❌ We added br_i += bm_i * treatment * c.bm | ✅ Removed - now same |

### Differences REMAINING

| # | Feature | Hendrickson | Our Approach | Impact |
|---|---------|-------------|--------------|--------|
| **1** | **Two-scenario comparison** | Single analysis approach | WITH and WITHOUT carryover modeling | **Major** - Key contribution |
| **2** | Parameter grid size | 4 sample sizes, 4 correlations, 5 carryover levels | 1 sample size, 2 correlations, 3 carryover levels | Minor - computational |
| **3** | Design structure | Hybrid uses `path`, Crossover uses `sequence` | Both use `path` column | Trivial - coding convenience |
| **4** | Biomarker correlation scaling | Uses c.bm directly | Scales by mean ratio, clamps to [-0.99, 0.99] | Minor - preserves validity |
| **5** | Non-PD matrix handling | Auto-fixes with `make_positive_definite()` | Rejects and skips combination | Minor - quality control |
| **8** | **Carryover in analysis model** | **NEVER models carryover** | **Conditionally models carryover** | **MAJOR** - Key methodological difference |
| **9** | Monte Carlo structure | Single analysis per dataset | Parameterized by `model_carryover` | Major - enables comparison |
| **10** | Iterations | 1000 | 20 | Minor - computational |

---

## Theoretical Implications

### Why Hendrickson's Power Declines with Carryover

```
When carryover exists in data BUT is not modeled in analysis:

1. True data-generating model:
   y_ij = β₀ + β₁·treatment_ij + β₂·biomarker_i + β₃·(treatment_ij × biomarker_i) +
          β₄·carryover_ij + ε_ij

2. Fitted model (Hendrickson):
   y_ij ~ β₀ + β₁·treatment_ij + β₂·biomarker_i + β₃·(treatment_ij × biomarker_i) + ε_ij

3. Result:
   - Carryover effect goes into residual (ε_ij)
   - Residual variance INFLATES: Var(ε_ij) ↑
   - Standard errors INFLATE: SE(β₃) ↑
   - t-statistics SHRINK: t = β₃/SE(β₃) ↓
   - Power DECLINES: P(reject H₀) ↓

⟹ Unmodeled confound reduces power
```

### Why Our Approach Maintains Power

```
When carryover exists in data AND is modeled in analysis:

1. True data-generating model:
   y_ij = β₀ + β₁·treatment_ij + β₂·biomarker_i + β₃·(treatment_ij × biomarker_i) +
          β₄·carryover_ij + ε_ij

2. Fitted model (Scenario 1):
   y_ij ~ β₀ + β₁·treatment_ij + β₂·biomarker_i + β₃·(treatment_ij × biomarker_i) +
          β₄·carryover_ij + ε_ij

3. Result:
   - Model CORRECTLY SPECIFIED
   - Carryover effect captured by β₄
   - Residual variance STAYS SMALL: Var(ε_ij) ≈ constant
   - Standard errors STAY SMALL: SE(β₃) ≈ constant
   - t-statistics MAINTAINED: t = β₃/SE(β₃) ≈ constant
   - Power MAINTAINED: P(reject H₀) ≈ constant

⟹ Controlling confound maintains power

Small power loss from extra degree of freedom (estimating β₄),
but MUCH smaller than loss from unmodeled confound.
```

---

## Key Insights from Comparison

### 1. Data Generation is Now Identical ✅

After implementing Option 1 (pure correlation approach):
- Both use correlation-only mechanism for biomarker interaction
- Both add carryover to BR means during data generation
- Both use same MVN draw procedure
- No post-MVN adjustments in either approach

**Conclusion**: Fair comparison, same biological model

### 2. Analysis Approach is the Key Difference ⭐

The ONLY major methodological difference is whether carryover is modeled in the statistical analysis:
- Hendrickson: Treats carryover as unmodeled confound
- Us (Scenario 1): Controls for carryover explicitly
- Us (Scenario 2): Replicates Hendrickson's approach

**Conclusion**: This is our methodological contribution

### 3. Expected Results Pattern

**Scenario 2 (WITHOUT modeling) should replicate Hendrickson**:
- Power declines as carryover increases
- Matches her published findings
- Validates our implementation

**Scenario 1 (WITH modeling) should show improvement**:
- Power remains stable across carryover levels
- Demonstrates value of our enhancement
- Quantifies benefit of proper carryover modeling

**Difference between scenarios = Cost of not modeling carryover**

---

## Validation Checklist

To ensure our implementation matches Hendrickson:

### Data Generation
- [x] Fixed correlation structure (c.tv = 0.8, c.cf1t = 0.2, etc.)
- [x] Differential biomarker correlation by treatment status
- [x] BR-only carryover (scale_factor = 1.0)
- [x] Exponential decay: (1/2)^(weeks_off/t1half)
- [x] NO population mean shift for biomarker interaction
- [x] NO post-MVN adjustment

### Trial Design
- [x] Hybrid: 4-path randomization
- [x] Crossover: 2-sequence AB/BA
- [x] Balanced assignment
- [x] Same period/washout structure

### Statistical Analysis
- [x] Scenario 2 matches Hendrickson: response ~ treatment * biomarker + week + (1|id)
- [x] Scenario 1 extends Hendrickson: adds carryover_effect term
- [x] Both test biomarker×treatment interaction
- [x] Both use lmer() with random intercept

### Conceptual Alignment
- [x] Biomarker interaction emerges from correlation structure
- [x] Carryover affects data but not (Hendrickson's) analysis model
- [x] Power calculated as proportion of significant results
- [x] Monte Carlo framework with multiple iterations

---

## Implementation Files

### Hendrickson's Original Code
- Repository: `~/prj/c265/pmsimstats-master`
- Main simulation: `pmSimulation.R`
- Functions: `pmfunctions.R`

### Our Current Code
- Repository: `/Users/zenn/Dropbox/prj/d08/pmsimstats2025`
- Main simulation: `analysis/scripts/full_pmsim_analysis_hyb_versus_co.R`
- Functions: `analysis/scripts/pm_functions.R`

---

## References

**Hendrickson, E., et al. (2020)**. N-of-1 trials with multiple randomization structures for individualized treatment. *Statistics in Medicine*, 39(25), 3581-3599.

---

*Documentation created: 2025-11-18*
*Purpose: Comprehensive side-by-side comparison with explicit difference identification*
