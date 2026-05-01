# Power Dropoff Analysis: The On/Off-Drug Variable Chain

**Technical analysis of the mechanism by which carryover
effects reduce statistical power in N-of-1 trial designs**

---

## 1. Overview

The power to detect a biomarker-treatment interaction in
N-of-1 and crossover trial designs depends entirely on the
contrast between on-drug and off-drug observations. Carryover
effects erode this contrast by allowing residual drug effect
to persist into off-drug periods. This document traces the
exact variable chain -- from trial design specification
through derived quantities to the analysis model interaction
term -- that mediates the power dropoff.

Five derived variables form this chain:

1. `ondrug` / `treatment` -- binary path assignment
2. `tod` / `weeks_on_drug` -- cumulative time on drug
3. `tsd` / `time_since_discontinuation` -- time off drug
4. `Dbc` / `effective_drug_weeks` -- continuous treatment
   indicator (the analysis predictor)
5. The interaction term (`bm * treatment_indicator`) --
   the coefficient whose significance defines power

Each variable transforms the preceding one. Carryover enters
at step 3 and propagates through steps 4 and 5 to reduce
the variance of the interaction term, inflating its standard
error and thus reducing power.

---

## 2. Variable 1: Path Assignment (`ondrug` / `treatment`)

The binary indicator of whether a participant receives active
drug at each measurement timepoint. This is the only variable
that differs across paths within a design.

### 2.1 orig implementation

Paths are defined as explicit binary vectors passed to
`buildtrialdesign()`:

```r
# N-of-1 (Hybrid) design, 8 timepoints
# Weeks: 4, 8, 9, 10, 11, 12, 16, 20
tdNof11 <- buildtrialdesign(
  ondrug = list(
    pathA = c(1,1,1,1,0,0,1,0),
    pathB = c(1,1,1,1,0,0,0,1),
    pathC = c(1,1,1,0,0,0,1,0),
    pathD = c(1,1,1,0,0,0,0,1)
  )
)
```

Each path vector has length equal to the number of
timepoints. Participants are randomly assigned to one path.
The function `buildtrialdesign()` then derives `tod`, `tsd`,
and `tpb` from these vectors.

### 2.2 main implementation

Paths are constructed procedurally in design-specific
functions. The Hybrid design:

```r
create_hybrid_design <- function(n_participants,
                                 measurement_weeks) {
  path_assignment <- sample(rep(1:4, length.out = n_participants))

  expand_grid(
    participant_id = 1:n_participants,
    week = measurement_weeks
  ) %>%
    mutate(
      path = path_assignment[participant_id],
      treatment = case_when(
        week %in% c(4, 8) ~ 1,           # OL phase: all on drug
        week == 9 ~ 1,                    # BD phase start
        week == 10 & path %in% c(1,2) ~ 1,
        week == 10 & path %in% c(3,4) ~ 0,
        week %in% c(11, 12) ~ 0,         # All off drug
        week == 16 & path %in% c(1,3) ~ 1,
        week == 16 & path %in% c(2,4) ~ 0,
        week == 20 & path %in% c(1,3) ~ 0,
        week == 20 & path %in% c(2,4) ~ 1,
        TRUE ~ NA_real_
      ),
      expectancy = if_else(week %in% c(4, 8), 1.0, 0.5)
    )
}
```

### 2.3 simple implementation

Uses a switch-based `randomize_assignment()` with block
randomization per design segment:

```r
randomize_assignment <- function(design) {
  on_drug <- design$on_drug_fixed
  if (design$name == "CO") {
    if (runif(1) < 0.5)
      on_drug <- c(rep(TRUE, 4), rep(FALSE, 4))
    else
      on_drug <- c(rep(FALSE, 4), rep(TRUE, 4))
  }
  # ... similar for N-of-1 and OL+BDC
  on_drug
}
```

### 2.4 Path structure comparison

The treatment indicator across the 8 timepoints for each
design's paths:

```
HYBRID (N-of-1) -- 4 paths
  Week:     4   8   9  10  11  12  16  20
  Path A:   1   1   1   1   0   0   1   0
  Path B:   1   1   1   1   0   0   0   1
  Path C:   1   1   1   0   0   0   1   0
  Path D:   1   1   1   0   0   0   0   1
  -------   OL phase  BD phase    CO phase
            (all on)  (staggered) (crossed)

OL+BDC -- 2 paths (orig) / 4 paths (main)
  Week:     4   8  12  16  17  18  19  20
  Path A:   1   1   1   1   1   1   0   0
  Path B:   1   1   1   1   1   0   0   0
  -------   OL phase (all on)    BD (off)

  main extends to 4 paths with finer staggering at
  weeks 17-18

CROSSOVER -- 2 paths
  Week:    2.5   5  7.5  10  12.5  15  17.5  20
  Path A:   1    1   1    1   0     0    0     0
  Path B:   0    0   0    0   1     1    1     1

OPEN-LABEL -- 1 path
  Week:    2.5   5  7.5  10  12.5  15  17.5  20
  Path A:   1    1   1    1   1     1    1     1
```

The number of off-drug timepoints per design determines the
upper bound on treatment contrast:

| Design   | On-drug TPs | Off-drug TPs | Treatment variance |
|----------|:-----------:|:------------:|:------------------:|
| OL       |     8       |      0       | 0 (no contrast)    |
| CO       |     4       |      4       | Maximum (balanced)  |
| OL+BDC   |    ~6       |     ~2       | Low (late only)     |
| Hybrid   |    ~5       |     ~3       | Moderate (spread)   |

---

## 3. Variable 2: Cumulative Time on Drug (`tod` / `weeks_on_drug`)

This is the cumulative sum of treatment-weighted time
intervals. It increases while on drug and plateaus (freezes)
when off drug. It feeds the Gompertz growth curve to
determine the biological response (BR) mean at each
timepoint.

### 3.1 orig implementation

In `buildtrialdesign()` (lines 90-96):

```r
od_intervals <- t_wk * od       # interval duration * ondrug
tod <- od_intervals
for (iTP in 2:length(timepoints)) {
  if (od_intervals[iTP] > 0) {
    tod[iTP] <- tod[iTP] + tod[iTP-1]  # on-drug: accumulate
  }
  # off-drug: tod stays at previous accumulated value
}
```

### 3.2 main implementation

In `simulation_clustered.R` (line 655):

```r
weeks_on_drug = cumsum(treatment)
```

Note: the main uses a counting-based cumulative sum (number
of on-drug timepoints) rather than a time-weighted sum.
For evenly spaced designs this is proportional; for
clustered designs, the distinction matters.

### 3.3 Worked example: Hybrid Path A

```
Week:            4    8    9   10   11   12   16   20
treatment:       1    1    1    1    0    0    1    0
weeks_on_drug:   1    2    3    4    4    4    5    5
                 ↑ accumulating ↑    ↑frozen↑  ↑+1  ↑frozen
```

The BR mean is computed from `tod`:

```r
BR_mean = modgompertz(tod, max, disp, rate)
```

Off-drug, `tod` freezes at 4, so the Gompertz-derived BR
mean plateaus. Without carryover, the BR mean drops to zero
immediately when treatment stops. With carryover, the
plateau value decays exponentially (see Variable 3).

---

## 4. Variable 3: Time Since Discontinuation (`tsd`)

This is the elapsed time since the participant was last on
drug. It is the direct input to the carryover decay function.

### 4.1 orig implementation

In `buildtrialdesign()` (lines 92-105):

```r
tsd <- t_wk * (od != 1)         # interval * off-drug
for (iTP in 2:length(timepoints)) {
  if (od_intervals[iTP] == 0) {  # off-drug
    tsd[iTP] <- tsd[iTP] + tsd[iTP-1]  # accumulate
  }
  # on-drug: tsd resets (stays at interval contribution)
}
everondrug <- (cumulative(od) > 0)
tsd <- everondrug * tsd          # zero out pre-drug periods
```

### 4.2 main implementation

In `simulation_clustered.R` (lines 658-665):

```r
just_discontinued = treatment == 0 &
                    lag(treatment, default = 0) == 1,
discontinuation_week = if_else(
  just_discontinued, week, NA_real_),
discontinuation_week = zoo::na.locf(
  discontinuation_week, na.rm = FALSE),
tsd = if_else(
  treatment == 0 & !is.na(discontinuation_week),
  week - discontinuation_week,
  0
)
```

The main uses absolute week differences rather than
accumulated intervals. This is more precise for irregularly
spaced timepoints.

### 4.3 simple implementation

In `generate_participant()`:

```r
if (on_drug[w]) {
  last_drug_week <- w
  carryover_factors[w] <- 0
} else if (last_drug_week > 0) {
  weeks_since <- w - last_drug_week
  carryover_factors[w] <- calc_carryover(
    weeks_since, carryover_halflife)
}
```

### 4.4 Worked example: Hybrid Path A

```
Week:            4    8    9   10   11   12   16   20
treatment:       1    1    1    1    0    0    1    0
tsd:             0    0    0    0    1    2    0    4

Carryover at t½ = 0.1 weeks (Hendrickson):
  week 11: (1/2)^(1/0.1) = (1/2)^10 = 0.001 (0.1%)
  week 12: (1/2)^(2/0.1) = (1/2)^20 = 0.000001 (~0%)
  week 20: (1/2)^(4/0.1) = (1/2)^40 ~ 0 (~0%)

Carryover at t½ = 0.2 weeks (Hendrickson):
  week 11: (1/2)^(1/0.2) = (1/2)^5 = 0.031 (3.1%)
  week 12: (1/2)^(2/0.2) = (1/2)^10 = 0.001 (0.1%)
  week 20: (1/2)^(4/0.2) = (1/2)^20 ~ 0 (~0%)

Carryover at t½ = 1 week (simple's values):
  week 11: (1/2)^(1/1) = 0.5 (50%)
  week 12: (1/2)^(2/1) = 0.25 (25%)
  week 20: (1/2)^(4/1) = 0.063 (6.3%)
```

This example reveals why the half-life parameterization
matters so much. At Hendrickson's values (0.1-0.2 weeks),
carryover is negligible after one week of discontinuation.
At the simple's values (0.5-2 weeks), substantial drug
effect persists for weeks, heavily contaminating the
off-drug contrast.

### 4.5 Scale factor complication

The orig and main apply a scale factor (default 2) inside
the carryover formula:

```r
# orig (generateData.R, line 96):
brmeans[p] <- brmeans[p] +
  brmeans[p-1] * (1/2)^(scalefactor * tsd / t1half)

# main (pm_functions.R, line 41):
previous_effect * (1/2)^(scale_factor * tsd / halflife)
```

With `scale_factor = 2`, the effective half-life is halved.
A nominal `t_half = 0.2` weeks becomes an effective
`t_half = 0.1` weeks. The simple omits this scale factor,
so its nominal values are also the effective values.

Effective half-lives across codebases:

| Codebase | Nominal t½ (wk) | Scale factor | Effective t½ (wk) |
|----------|:---------------:|:------------:|:------------------:|
| orig     | 0.1, 0.2        | 2            | 0.05, 0.10         |
| main*    | 0.1, 0.2        | 1**          | 0.10, 0.20         |
| simple   | 0.5, 1, 2       | 1            | 0.50, 1.00, 2.00   |

*main's `simulation_clustered.R` uses its own
`calculate_carryover()` without scale factor;
`pm_functions.R` retains scale_factor = 2 for the shared
function.

**This means the main's clustered simulation and the orig
produce different carryover decay rates even for the same
nominal half-life value.**

---

## 5. Variable 4: Continuous Treatment Indicator

This is the variable that actually enters the analysis model
as a predictor. It transforms the binary treatment state
into a continuous measure of drug exposure that accounts for
(or ignores) carryover.

### 5.1 orig: `Dbc`

In `lme_analysis.R` (lines 107-109):

```r
data.m2[Db == TRUE, Dbc := 1]
data.m2[Db == FALSE,
  Dbc := (1/2)^(op$carryover_scalefactor * tsd /
                op$carryover_t1half)]
```

`Dbc` is bounded [0, 1]:
- On drug: 1 (full exposure)
- Off drug: exponentially decaying from 1 toward 0

### 5.2 main: `effective_drug_weeks`

In `simulation_clustered.R` (lines 697-701):

```r
effective_drug_weeks = case_when(
  treatment == 1 ~ weeks_on_drug,
  treatment == 0 & !is.na(weeks_at_discont) ~
    weeks_at_discont * carryover_effect,
  TRUE ~ 0
)
```

`effective_drug_weeks` is bounded [0, max(weeks_on_drug)]:
- On drug: cumulative count (1, 2, 3, ...)
- Off drug: accumulated count at discontinuation *
  exponential decay

The misspecified alternative (`model_carryover = FALSE`) uses
`weeks_on_drug` directly, which plateaus off-drug rather than
decaying.

### 5.3 simple: carryover adjustment (optional)

In `analyze_trial()`:

```r
# With adjustment:
adjusted_response = response - carryover * br_rate

# Without adjustment:
# raw response used directly
```

Rather than modifying the predictor, the simple subtracts
the estimated carryover effect from the response variable
before computing phase means.

### 5.4 Worked example: Hybrid Path A

```
Week:                 4    8    9   10   11   12   16   20
treatment:            1    1    1    1    0    0    1    0

orig Dbc (t½=0.2, sf=2):
                      1    1    1    1  .001   ~0    1   ~0

main effective_drug_weeks (t½=0.2):
  weeks_on_drug:      1    2    3    4    4    4    5    5
  carryover_effect:   1    1    1    1  .031 .001    1   ~0
  effective:          1    2    3    4  .124 .004    5   ~0

main weeks_on_drug (no carryover model):
                      1    2    3    4    4    4    5    5
```

Note the difference between the two main scenarios:
- `effective_drug_weeks` drops sharply off-drug (good
  contrast with on-drug)
- `weeks_on_drug` stays at 4 off-drug (poor contrast:
  the model sees accumulated exposure = 4 both on and
  off drug)

---

## 6. Variable 5: The Interaction Term

The coefficient whose p-value defines statistical power.

### 6.1 Analysis model formulas

**orig** (`lme_analysis.R`):

When treatment varies across participants (CO, Hybrid,
OL+BDC):
```
Sx ~ bm + t + Dbc + bm * Dbc + (1 | ptID)
```
Interaction of interest: `bm:Dbc`

When treatment is constant (OL):
```
Sx ~ bm + t + bm * t + (1 | ptID)
```
Interaction of interest: `bm:t`

**main** (`simulation_clustered.R`):

With carryover modeled:
```
response ~ effective_drug_weeks * bm_centered + week
         + (1 | participant_id)
```
or with `nlme::lme`:
```
response ~ effective_drug_weeks * bm_centered + week,
  random = ~ 1 | participant_id,
  correlation = corCAR1(form = ~ week | participant_id)
```
Interaction of interest: `effective_drug_weeks:bm_centered`

Without carryover modeled:
```
response ~ weeks_on_drug * bm_centered + week
         + (1 | participant_id)
```
Interaction of interest: `weeks_on_drug:bm_centered`

**simple** (`analyze_trial()`):

For CO and N-of-1:
```
drug_effect ~ biomarker
```
where `drug_effect = mean(on_drug) - mean(off_drug)`
per participant.

For OL+BDC:
```
(blinded - open) ~ biomarker * group
```
Interaction of interest: `biomarker:groupplacebo`

### 6.2 How carryover degrades the interaction

The interaction coefficient estimates the degree to which
the biomarker predicts differential response to treatment.
Its standard error is inversely proportional to the variance
of the treatment indicator:

```
SE(beta_interaction) ~ 1 / sqrt(Var(treatment_indicator))
```

Carryover reduces `Var(treatment_indicator)` by making
off-drug observations resemble on-drug observations:

```
Without carryover:
  On-drug predictor values:  [1, 2, 3, 4, 5]
  Off-drug predictor values: [0, 0, 0]
  Contrast: large

With carryover (t½ = 1 week):
  On-drug predictor values:  [1, 2, 3, 4, 5]
  Off-drug predictor values: [2.0, 1.0, 0.3]
  Contrast: reduced

With carryover (t½ = 2 weeks):
  On-drug predictor values:  [1, 2, 3, 4, 5]
  Off-drug predictor values: [2.8, 2.0, 1.2]
  Contrast: severely reduced
```

When the treatment indicator has low variance, the
regression cannot distinguish the biomarker-moderated
drug effect from noise, and power drops.

---

## 7. Design-Specific Vulnerability

### 7.1 Open-Label (OL)

No off-drug timepoints exist. The interaction uses `bm * t`
(biomarker x elapsed time) rather than
`bm * treatment_indicator`. Carryover is irrelevant because
there is no discontinuation.

Power depends on: biomarker-time interaction strength.
Carryover sensitivity: **none**.

### 7.2 Traditional Crossover (CO)

Four on-drug and four off-drug timepoints provide maximum
treatment contrast. However, the entire second period is
vulnerable to carryover from the first period. For Path A
(drug first), all off-drug measurements (weeks 12.5-20)
carry residual drug effect.

```
Path A (drug -> placebo):
  Week:   2.5   5  7.5   10  12.5  15   17.5  20
  Tx:      1    1    1    1    0     0     0    0
  tsd:     0    0    0    0   2.5   5.0   7.5  10
  CO at t½=0.2: --  --   --   --  .0001  ~0   ~0   ~0
  CO at t½=2:   --  --   --   --  .42   .18   .07  .03
```

At short half-lives (Hendrickson), carryover is negligible
by the first off-drug measurement. At long half-lives
(simple), substantial contamination persists.

Power depends on: on/off contrast in second half.
Carryover sensitivity: **high** (long tsd values help
with Hendrickson values; problematic with longer t½).

### 7.3 Open-Label + Blinded Discontinuation (OL+BDC)

Long open-label phase (weeks 4-16 on drug) followed by
short blinded discontinuation (weeks 17-20, staggered
off-drug). Only 2-4 off-drug timepoints exist, all at the
end of the trial.

```
OL+BDC Path A (orig):
  Week:    4    8   12   16   17   18   19   20
  Tx:      1    1    1    1    1    1    0    0
  tsd:     0    0    0    0    0    0    1    2
```

The tsd values are short (1-2 weeks from last on-drug
timepoint to off-drug measurements), which means even
modest carryover half-lives can produce substantial
contamination:

```
  CO at t½=0.1: week 19: 0.001, week 20: ~0
  CO at t½=0.2: week 19: 0.031, week 20: 0.001
  CO at t½=1.0: week 19: 0.500, week 20: 0.250
```

Power depends on: late on/off contrast with few off-drug
timepoints.
Carryover sensitivity: **moderate to high** (few off-drug
timepoints; short tsd means even moderate t½ contaminates).

### 7.4 Hybrid (N-of-1)

Combines open-label, blinded discontinuation, and crossover
phases. Four paths create varied discontinuation patterns:

```
Hybrid Path A:
  Week:    4    8    9   10   11   12   16   20
  Tx:      1    1    1    1    0    0    1    0
  tsd:     0    0    0    0    1    2    0    4

Hybrid Path D:
  Week:    4    8    9   10   11   12   16   20
  Tx:      1    1    1    0    0    0    0    1
  tsd:     0    0    0    0    1    2    3    0
```

Multiple discontinuation events create multiple contrast
windows with different tsd values, providing the analysis
model with richer information about the decay trajectory.

Power depends on: multiple contrast windows across phases.
Carryover sensitivity: **moderate** (varied tsd values mean
some contrasts are clean even when others are contaminated).

---

## 8. Summary Table

```
+-----------+----------+----------+----------+-----------+
| Variable  | orig     | main     | simple   | Role      |
+-----------+----------+----------+----------+-----------+
| treatment | ondrug   | treatment| on_drug  | Binary    |
| (V1)      | (list of | (case_   | (block   | path      |
|           | vectors) | when)    | random)  | assignment|
+-----------+----------+----------+----------+-----------+
| cum. drug | tod      | weeks_on | (not     | Gompertz  |
| time (V2) | (accum   | _drug    | cumul.)  | input for |
|           | loop)    | (cumsum) |          | BR mean   |
+-----------+----------+----------+----------+-----------+
| time off  | tsd      | tsd      | weeks_   | Carryover |
| drug (V3) | (accum   | (week -  | since    | decay     |
|           | loop)    | discont) |          | input     |
+-----------+----------+----------+----------+-----------+
| continuous| Dbc      | effective| carryover| Analysis  |
| treatment | [0,1]    | _drug_wks| _factors | model     |
| (V4)      |          | [0,max]  | [0,1]    | predictor |
+-----------+----------+----------+----------+-----------+
| interact. | bm:Dbc   | eff_drug | biomarker| Coeff     |
| term (V5) | or bm:t  | :bm_cent | (on diff)| tested    |
+-----------+----------+----------+----------+-----------+
```

---

## 9. Variable Equivalence: orig vs. main

This section traces each derived variable through both
codebases for every design, using the exact inputs from the
orig vignette and the main's `simulation_clustered.R`. The
goal is to determine whether the two produce identical
treatment-related variables.

### 9.1 Hybrid (N-of-1) Design

**Inputs (orig vignette, lines 48-61):**

```r
timepoints = c(4, 8, 9, 10, 11, 12, 16, 20)
ondrug = list(
  pathA = c(1,1,1,1,0,0,1,0),
  pathB = c(1,1,1,1,0,0,0,1),
  pathC = c(1,1,1,0,0,0,1,0),
  pathD = c(1,1,1,0,0,0,0,1)
)
expectancies = c(1, 1, .5, .5, .5, .5, .5, .5)
```

**Inputs (main, lines 366, 394-417):**

```r
measurement_weeks_hybrid = c(4, 8, 9, 10, 11, 12, 16, 20)
treatment = case_when(
  week %in% c(4, 8) ~ 1,
  week == 9 ~ 1,
  week == 10 & path %in% c(1, 2) ~ 1,
  week == 10 & path %in% c(3, 4) ~ 0,
  week %in% c(11, 12) ~ 0,
  week == 16 & path %in% c(1, 3) ~ 1,
  week == 16 & path %in% c(2, 4) ~ 0,
  week == 20 & path %in% c(1, 3) ~ 0,
  week == 20 & path %in% c(2, 4) ~ 1
)
```

#### 9.1.1 treatment (V1): EQUAL

```
Week:       4  8  9 10 11 12 16 20
orig A:     1  1  1  1  0  0  1  0
main 1:     1  1  1  1  0  0  1  0  ✓

orig B:     1  1  1  1  0  0  0  1
main 2:     1  1  1  1  0  0  0  1  ✓

orig C:     1  1  1  0  0  0  1  0
main 3:     1  1  1  0  0  0  1  0  ✓

orig D:     1  1  1  0  0  0  0  1
main 4:     1  1  1  0  0  0  0  1  ✓
```

The binary treatment vectors are identical across all four
paths.

#### 9.1.2 tod / weeks_on_drug (V2): NOT EQUAL

The orig computes `tod` as cumulative *time-weighted*
exposure (in weeks). The main computes `weeks_on_drug` as
a simple *count* of on-drug timepoints via `cumsum(treatment)`.

**orig `buildtrialdesign()` logic:**

```r
t_wk = c(4, 4, 1, 1, 1, 1, 4, 4)  # interval durations
od_intervals = t_wk * od            # weighted by ondrug
```

For Path A (`od = c(1,1,1,1,0,0,1,0)`):

```r
od_intervals = c(4, 4, 1, 1, 0, 0, 4, 0)
```

Accumulation loop: `tod[i] = tod[i] + tod[i-1]` when
`od_intervals[i] > 0`:

```
i=1: tod = 4
i=2: od_intervals[2]=4 > 0: tod = 4 + 4 = 8
i=3: od_intervals[3]=1 > 0: tod = 1 + 8 = 9
i=4: od_intervals[4]=1 > 0: tod = 1 + 9 = 10
i=5: od_intervals[5]=0:     tod = 0 (no accumulation)
i=6: od_intervals[6]=0:     tod = 0
i=7: od_intervals[7]=4 > 0: tod = 4 + 0 = 4
     ^^^ RESETS! Does NOT carry forward from i=4.
i=8: od_intervals[8]=0:     tod = 0
```

Wait -- that cannot be right. Let me re-read the loop
more carefully. The accumulation is:

```r
tod <- od_intervals          # initialize
for (iTP in 2:length) {
  if (od_intervals[iTP] > 0) {
    tod[iTP] <- tod[iTP] + tod[iTP-1]
  }
}
```

So `tod[iTP]` carries forward from `tod[iTP-1]`, not from
the last on-drug value. At i=7, `tod[iTP-1] = tod[6] = 0`,
so `tod[7] = 4 + 0 = 4`.

**orig tod for Path A:**

```
Week:    4    8    9   10   11   12   16   20
tod:     4    8    9   10    0    0    4    0
```

**main weeks_on_drug for Path 1:**

```r
weeks_on_drug = cumsum(treatment)
treatment:       1    1    1    1    0    0    1    0
cumsum:          1    2    3    4    4    4    5    5
```

```
Week:    4    8    9   10   11   12   16   20
wod:     1    2    3    4    4    4    5    5
```

**These are NOT equal.** They differ in two fundamental
ways:

1. **Units.** The orig's `tod` is in *weeks* of drug
   exposure (weighted by interval duration). The main's
   `weeks_on_drug` is a *count* of on-drug measurement
   occasions.

2. **Continuity across gaps.** The orig's `tod` resets
   after an off-drug gap because it only carries forward
   from `tod[iTP-1]`, which is 0 during off-drug periods.
   After resuming drug at week 16, `tod = 4` (the 4-week
   interval), not `10 + 4 = 14`. The main's `cumsum`
   never resets -- it monotonically increases, reaching
   5 at week 16.

**Consequence for BR mean computation:**

In the orig, BR mean is computed via `modgompertz(tod, ...)`.
After an off-drug gap, `tod` restarts from the current
interval duration. This means the biological response curve
*restarts* from a low value after re-initiation, as though
the previous drug exposure was partially forgotten.

In the main, BR mean is computed as
`weeks_on_drug * effective_BR_rate`. Since `weeks_on_drug`
monotonically increases, the BR mean after re-initiation
continues from the accumulated count (5, not restarting
from 1). The growth trajectory does not reset.

This is a **substantive methodological difference** in how
drug re-initiation is modeled. The orig treats each drug
episode as partially independent (Gompertz restarts from
a reduced `tod`). The main treats total exposure as
cumulative regardless of gaps.

#### 9.1.3 tsd (V3): NOT EQUAL

**orig logic:**

```r
tsd <- t_wk * (od != 1)     # interval * off-drug
for (iTP in 2:length) {
  if (od_intervals[iTP] == 0) {  # off-drug
    tsd[iTP] <- tsd[iTP] + tsd[iTP-1]
  }
}
everondrug <- (cumulative(od) > 0)
tsd <- everondrug * tsd
```

For Path A:

```r
t_wk = c(4, 4, 1, 1, 1, 1, 4, 4)
od   = c(1, 1, 1, 1, 0, 0, 1, 0)

tsd_init = t_wk * (od != 1) = c(0, 0, 0, 0, 1, 1, 0, 4)

Accumulation (only when od_intervals[i] == 0):
i=2: od_int[2]=4 > 0, skip
i=3: od_int[3]=1 > 0, skip
i=4: od_int[4]=1 > 0, skip
i=5: od_int[5]=0: tsd[5] = 1 + tsd[4] = 1 + 0 = 1
i=6: od_int[6]=0: tsd[6] = 1 + tsd[5] = 1 + 1 = 2
i=7: od_int[7]=4 > 0, skip
i=8: od_int[8]=0: tsd[8] = 4 + tsd[7] = 4 + 0 = 4

everondrug = cumulative(c(1,1,1,1,0,0,1,0))
           = c(1, 2, 3, 4, 4, 4, 5, 5) > 0
           = c(T, T, T, T, T, T, T, T)
tsd = everondrug * tsd (all TRUE, no effect)
```

**orig tsd for Path A:**

```
Week:    4    8    9   10   11   12   16   20
tsd:     0    0    0    0    1    2    0    4
```

**main logic (lines 658-665):**

```r
just_discontinued = treatment == 0 &
                    lag(treatment, default = 0) == 1
discontinuation_week = if_else(just_discontinued,
                               week, NA_real_)
discontinuation_week = zoo::na.locf(discontinuation_week,
                                    na.rm = FALSE)
tsd = if_else(
  treatment == 0 & !is.na(discontinuation_week),
  week - discontinuation_week, 0)
```

For Path 1:

```
Week:           4   8   9  10  11  12  16  20
treatment:      1   1   1   1   0   0   1   0
lag(treatment): 0   1   1   1   1   0   0   1
just_discont:   F   F   F   F   T   F   F   T

discontinuation_week before locf:
                NA  NA  NA  NA  11  NA  NA  20

discontinuation_week after locf:
                NA  NA  NA  NA  11  11  11  20

tsd:
  wk 4:  treatment=1 -> 0
  wk 8:  treatment=1 -> 0
  wk 9:  treatment=1 -> 0
  wk 10: treatment=1 -> 0
  wk 11: treatment=0, discont_wk=11: 11-11 = 0
  wk 12: treatment=0, discont_wk=11: 12-11 = 1
  wk 16: treatment=1 -> 0
  wk 20: treatment=0, discont_wk=20: 20-20 = 0
```

**main tsd for Path 1:**

```
Week:    4    8    9   10   11   12   16   20
tsd:     0    0    0    0    0    1    0    0
```

**orig vs. main tsd for Path A/1:**

```
Week:    4    8    9   10   11   12   16   20
orig:    0    0    0    0    1    2    0    4
main:    0    0    0    0    0    1    0    0
                            ^^^  ^^^       ^^^
```

**These are NOT equal.** They differ at weeks 11, 12,
and 20:

- **Week 11:** orig = 1, main = 0. The orig counts the
  interval from week 10 to week 11 (1 week) as time
  since drug. The main computes `week - discontinuation_week
  = 11 - 11 = 0` because the discontinuation happens
  *at* week 11 (the first off-drug measurement), so
  tsd = 0 at the moment of discontinuation.

- **Week 12:** orig = 2, main = 1. Difference of 1 week,
  consistent with the offset at week 11.

- **Week 20:** orig = 4, main = 0. The orig accumulates
  4 weeks of off-drug time (the interval from week 16 to
  week 20). The main sets `discontinuation_week = 20`
  (week 20 is the first off-drug timepoint after resuming
  at week 16), so `tsd = 20 - 20 = 0`.

**Root cause:** The two codebases define "time since
discontinuation" differently:

- **orig:** `tsd` at a given timepoint is the cumulative
  off-drug *interval duration* leading up to and including
  that timepoint. It measures "how much off-drug time has
  elapsed since the last on-drug interval ended."

- **main:** `tsd` at a given timepoint is
  `week - discontinuation_week`, where
  `discontinuation_week` is the *week number of the first
  off-drug observation*. At the first off-drug timepoint
  itself, `tsd = 0`.

The orig's convention is that tsd represents the time
*already elapsed* off drug by the time a measurement is
taken. The main's convention is that tsd represents the
time elapsed since the *transition point*.

This produces a systematic offset: **the orig's tsd is
approximately one measurement interval larger than the
main's at every off-drug timepoint.** For the Hybrid
design's BD phase (1-week intervals), the offset is
1 week. For the CO phase (4-week interval at week 20),
the offset is 4 weeks.

#### 9.1.4 Db (derived on-drug indicator): NOT EQUAL

The orig derives `Db` in `lme_analysis.R` (line 102) as:

```r
Db = (tod > 0)
```

For Path A:

```
Week:    4    8    9   10   11   12   16   20
tod:     4    8    9   10    0    0    4    0
Db:      T    T    T    T    F    F    T    F
```

The main uses `treatment` directly (which equals `ondrug`):

```
Week:    4    8    9   10   11   12   16   20
tx:      1    1    1    1    0    0    1    0
```

These happen to be equal for Path A because `tod > 0`
coincides with `ondrug == 1` (tod is zero exactly when
ondrug is zero). **Equal for this design.**

#### 9.1.5 Carryover effect: NOT EQUAL

The carryover decay applied to the BR mean uses the
tsd values computed above. Since tsd differs, the
carryover effect differs.

**orig carryover (generateData.R, line 96):**

```r
brmeans[p] <- brmeans[p] +
  brmeans[p-1] * (1/2)^(scalefactor * tsd[p] / t1half)
```

With `scalefactor = 2`, `t1half = 0.2`:

```
Week 11: tsd=1: (1/2)^(2*1/0.2) = (1/2)^10 = 0.00098
Week 12: tsd=2: (1/2)^(2*2/0.2) = (1/2)^20 ≈ 0
Week 20: tsd=4: (1/2)^(2*4/0.2) = (1/2)^40 ≈ 0
```

**main carryover (simulation_clustered.R, line 679):**

```r
carryover_effect = (1/2)^(tsd / carryover_t_half)
```

With `t1half = 0.2`:

```
Week 11: tsd=0: carryover_effect = 1 (on-drug formula
         applies since treatment=0 but tsd=0 means
         the case_when hits: (1/2)^(0/0.2) = 1)
         But wait -- treatment=0, so this goes to the
         off-drug branch. (1/2)^(0/0.2) = (1/2)^0 = 1
Week 12: tsd=1: (1/2)^(1/0.2) = (1/2)^5 = 0.03125
Week 20: tsd=0: (1/2)^(0/0.2) = 1
```

**Comparison at t½ = 0.2 weeks:**

```
Week:         11       12       20
orig CO:      0.00098  ~0       ~0
main CO:      1.000    0.031    1.000
              ^^^^^    ^^^^^    ^^^^^
```

**These are dramatically different.** The main computes
a carryover_effect of 1.0 (full retention) at weeks 11
and 20 because tsd = 0 at those timepoints. The orig
computes near-zero retention because its tsd is 1 and 4
weeks respectively.

The consequence is that the main's BR mean at week 11
(first off-drug measurement) retains the full accumulated
drug effect, while the orig's BR mean decays it almost
completely. At week 20 (off-drug after re-initiation),
the main retains the full effect while the orig decays
it to zero.

This is a **critical divergence** that directly affects
power calculations. The main's model allows much more
residual drug effect during off-drug periods than the
orig's model, particularly at the first off-drug
measurement after each discontinuation.

#### 9.1.6 expectancy: EQUAL

```
Week:       4    8    9   10   11   12   16   20
orig:       1    1   .5   .5   .5   .5   .5   .5
main:       1    1   .5   .5   .5   .5   .5   .5  ✓
```

### 9.2 OL+BDC Design

**Inputs (orig vignette, lines 75-85):**

```r
timepoints = c(4, 8, 12, 16, 17, 18, 19, 20)
ondrug = list(
  pathA = c(1,1,1,1,1,1,0,0),
  pathB = c(1,1,1,1,1,0,0,0)
)
expectancies = c(1, 1, 1, 1, .5, .5, .5, .5)
```

**Inputs (main, lines 373-392):**

```r
measurement_weeks_ol_bdc = c(4, 8, 12, 16, 17, 18, 19, 20)
treatment = case_when(
  week <= 16 ~ 1,
  week == 17 ~ 1,
  week == 18 & path %in% c(1, 2) ~ 1,
  week == 18 & path %in% c(3, 4) ~ 0,
  week %in% c(19, 20) ~ 0,
  TRUE ~ 1
)
```

#### 9.2.1 treatment (V1): NOT EQUAL (different paths)

The orig defines **2 paths**:

```
Week:       4  8 12 16 17 18 19 20
orig A:     1  1  1  1  1  1  0  0
orig B:     1  1  1  1  1  0  0  0
```

The main defines **4 paths** with a different staggering:

```
Week:       4  8 12 16 17 18 19 20
main 1:     1  1  1  1  1  1  0  0
main 2:     1  1  1  1  1  1  0  0
main 3:     1  1  1  1  1  0  0  0
main 4:     1  1  1  1  1  0  0  0
```

Main paths 1-2 are identical to orig path A; main paths
3-4 are identical to orig path B. The main uses 4 paths
(equal assignment) where the orig uses 2 paths. **The
effective treatment vectors are the same two patterns,
but with 4-path assignment in the main vs. 2-path in
the orig.** The statistical content is equivalent --
both produce a 50/50 split between the two discontinuation
schedules.

However, there is a subtle difference at week 17:

```
orig A: week 17, ondrug = 1
orig B: week 17, ondrug = 1
main 1: week 17, treatment = 1 (week <= 16 fails,
         week == 17 ~ 1)
main 3: week 17, treatment = 1 (same)
```

Both agree: all paths are on drug at week 17. The first
off-drug measurement is at week 18 (paths B/3-4) or
week 19 (paths A/1-2).

**The 2-vs-4 path distinction is cosmetic for OL+BDC.**
Both produce identical treatment patterns at 50/50 ratio.

#### 9.2.2 tod / weeks_on_drug (V2): NOT EQUAL

**orig tod for Path A:**

```r
t_wk = c(4, 4, 4, 4, 1, 1, 1, 1)
od   = c(1, 1, 1, 1, 1, 1, 0, 0)

od_intervals = c(4, 4, 4, 4, 1, 1, 0, 0)
tod accumulation:
  i=1: 4
  i=2: 4+4  = 8
  i=3: 4+8  = 12
  i=4: 4+12 = 16
  i=5: 1+16 = 17
  i=6: 1+17 = 18
  i=7: od_int=0, no accum: tod = 0
  i=8: od_int=0, no accum: tod = 0
```

```
Week:    4    8   12   16   17   18   19   20
tod:     4    8   12   16   17   18    0    0
```

**main weeks_on_drug for Path 1:**

```
Week:    4    8   12   16   17   18   19   20
tx:      1    1    1    1    1    1    0    0
cumsum:  1    2    3    4    5    6    6    6
```

**Differences:** Same two issues as Hybrid:
1. Units (weeks vs. count)
2. Off-drug behavior: orig drops to 0, main plateaus at 6

#### 9.2.3 tsd (V3): NOT EQUAL

**orig tsd for Path A:**

```r
tsd_init = t_wk * (od != 1) = c(0,0,0,0,0,0,1,1)
Accumulation:
  i=7: od_int=0: tsd = 1 + tsd[6] = 1 + 0 = 1
  i=8: od_int=0: tsd = 1 + tsd[7] = 1 + 1 = 2
```

```
Week:    4    8   12   16   17   18   19   20
tsd:     0    0    0    0    0    0    1    2
```

**main tsd for Path 1:**

```
Week:           4   8  12  16  17  18  19  20
treatment:      1   1   1   1   1   1   0   0
lag(treatment): 0   1   1   1   1   1   1   0
just_discont:   F   F   F   F   F   F   T   F
discont_week:  NA  NA  NA  NA  NA  NA  19  19
tsd:            0   0   0   0   0   0   0   1
```

```
Week:    4    8   12   16   17   18   19   20
orig:    0    0    0    0    0    0    1    2
main:    0    0    0    0    0    0    0    1
                                      ^^   ^^
```

Same systematic offset: the orig's tsd is one interval
ahead of the main's. At the first off-drug timepoint
(week 19), orig says tsd = 1 (one week has elapsed off
drug) while main says tsd = 0 (just discontinued).

### 9.3 Traditional Crossover (CO) Design

**Inputs (orig vignette, lines 87-98):**

```r
timepoints = cumulative(rep(2.5, 8))
# = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)
ondrug = list(
  pathA = c(1,1,1,1,0,0,0,0),
  pathB = c(0,0,0,0,1,1,1,1)
)
expectancies = rep(.5, 8)
```

**Inputs (main, lines 419-442):**

```r
measurement_weeks_crossover =
  c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20)
treatment = case_when(
  week <= 10 & path == 1 ~ 1,
  week <= 10 & path == 2 ~ 0,
  week >  10 & path == 1 ~ 0,
  week >  10 & path == 2 ~ 1
)
```

#### 9.3.1 treatment (V1): EQUAL

```
Week:      2.5  5  7.5  10  12.5  15  17.5  20
orig A:     1   1   1    1    0    0    0     0
main 1:     1   1   1    1    0    0    0     0  ✓

orig B:     0   0   0    0    1    1    1     1
main 2:     0   0   0    0    1    1    1     1  ✓
```

#### 9.3.2 tod / weeks_on_drug (V2): NOT EQUAL

**orig tod for Path A:**

```r
t_wk = rep(2.5, 8)
od_intervals = c(2.5, 2.5, 2.5, 2.5, 0, 0, 0, 0)
tod:  2.5  5.0  7.5  10.0  0  0  0  0
```

**main weeks_on_drug for Path 1:**

```
cumsum(1,1,1,1,0,0,0,0) = 1  2  3  4  4  4  4  4
```

Different units (weeks vs. count) and different off-drug
behavior (0 vs. plateau at 4).

#### 9.3.3 tsd (V3): NOT EQUAL

**orig tsd for Path A:**

```r
tsd_init = c(0, 0, 0, 0, 2.5, 2.5, 2.5, 2.5)
Accumulation:
  i=5: 2.5 + 0   = 2.5
  i=6: 2.5 + 2.5 = 5.0
  i=7: 2.5 + 5.0 = 7.5
  i=8: 2.5 + 7.5 = 10.0
```

```
Week:   2.5    5   7.5   10   12.5   15   17.5   20
tsd:     0     0    0     0    2.5   5.0   7.5   10.0
```

**main tsd for Path 1:**

```
just_discont at week 12.5: TRUE (tx=0, lag=1)
discont_week:
  NA  NA  NA  NA  12.5  12.5  12.5  12.5
tsd:
  0   0   0   0   0     2.5   5.0   7.5
```

```
Week:   2.5    5   7.5   10   12.5   15   17.5   20
orig:    0     0    0     0    2.5   5.0   7.5   10.0
main:    0     0    0     0    0     2.5   5.0    7.5
                                ^^^   ^^^   ^^^    ^^^
```

Systematic offset of 2.5 weeks (one measurement interval)
at every off-drug timepoint.

#### 9.3.4 tsd for Path B (placebo-first): DIFFERENT ISSUE

**orig tsd for Path B:**

```r
od = c(0,0,0,0,1,1,1,1)
tsd_init = c(2.5, 2.5, 2.5, 2.5, 0, 0, 0, 0)

Before accumulation:
  i=2: od_int=0: tsd = 2.5 + 2.5 = 5.0
  i=3: od_int=0: tsd = 2.5 + 5.0 = 7.5
  i=4: od_int=0: tsd = 2.5 + 7.5 = 10.0
  i=5-8: od_int > 0: skip

everondrug = cumulative(0,0,0,0,1,1,1,1) > 0
           = c(0,0,0,0,1,2,3,4) > 0
           = c(F, F, F, F, T, T, T, T)
tsd = everondrug * tsd
    = c(0, 0, 0, 0, 0, 0, 0, 0)
```

```
Week:   2.5    5   7.5   10   12.5   15   17.5   20
tsd:     0     0    0     0     0     0     0      0
```

Path B never has off-drug time *after* being on drug, so
tsd is zero everywhere. This is correct: carryover only
applies after drug exposure.

**main tsd for Path 2:**

Path 2 has `treatment = c(0,0,0,0,1,1,1,1)`. There is
no discontinuation (no transition from 1 to 0), so
`just_discontinued` is never TRUE. `tsd = 0` everywhere.

**Path B/2: EQUAL** (both zero everywhere).

### 9.4 Open-Label (OL) Design

**Inputs (orig, lines 64-73):**

```r
timepoints = cumulative(rep(2.5, 8))
ondrug = list(pathA = rep(1, 8))
expectancies = rep(1, 8)
```

**Inputs (main, lines 444-457):**

```r
treatment = 1  # Always on drug
expectancy = 1.0
```

All variables trivially equal: treatment = 1 everywhere,
tod accumulates continuously, tsd = 0 everywhere.
**EQUAL.**

### 9.5 Summary of Divergences

```
+-------------+-----------+-------------+-----------+
| Variable    | OL        | CO          | Hybrid    |
|             |           |             | & OL+BDC  |
+-------------+-----------+-------------+-----------+
| treatment   | Equal     | Equal       | Equal     |
+-------------+-----------+-------------+-----------+
| tod /       | Not equal | Not equal   | Not equal |
| weeks_on_   | (units    | (units +    | (units +  |
| drug        | only)     | off-drug    | off-drug  |
|             |           | behavior)   | behavior) |
+-------------+-----------+-------------+-----------+
| tsd         | Equal     | Not equal   | Not equal |
|             | (both 0)  | (offset by  | (offset   |
|             |           | 1 interval) | + reset)  |
+-------------+-----------+-------------+-----------+
| Db / tx     | Equal     | Equal       | Equal     |
| (binary)    |           |             |           |
+-------------+-----------+-------------+-----------+
| carryover   | Equal     | Not equal   | Not equal |
| effect      | (both 0)  | (offset     | (offset   |
|             |           | propagates) | + reset)  |
+-------------+-----------+-------------+-----------+
```

### 9.6 Root Causes

**Divergence 1: `tod` units and accumulation.**
The orig accumulates actual time-on-drug in *weeks*
(weighted by interval duration). The main counts
on-drug *measurement occasions* via `cumsum(treatment)`.
For evenly spaced designs (CO, OL) these differ by a
constant factor (2.5x for CO). For irregularly spaced
designs (Hybrid, OL+BDC), the relationship is
non-linear because intervals range from 1 to 4 weeks.

This affects the BR mean computation. The orig feeds
`tod` (in weeks) into a Gompertz curve. The main feeds
`weeks_on_drug` (count) into a linear rate model
(`weeks_on_drug * BR_rate`). These produce different
dose-response shapes independent of the unit difference.

**Divergence 2: `tod` off-drug behavior.**
The orig's `tod` carries forward only from the immediately
preceding timepoint. During off-drug periods, `tod[i-1]`
is 0 (since `tod` was not accumulated), so when drug
resumes, `tod` restarts from the new interval's duration
rather than from the pre-gap accumulated value. The main's
`cumsum(treatment)` is monotonically non-decreasing and
plateaus at the last on-drug count.

This means the orig models drug re-initiation as a
partial restart (lower Gompertz input), while the main
models it as continuation of cumulative exposure.

**Divergence 3: `tsd` measurement convention.**
The orig's `tsd` at a timepoint measures the total
off-drug interval duration from the last on-drug interval
through the current timepoint. At the *first* off-drug
measurement, this equals the interval duration (typically
1-2.5 weeks). The main's `tsd` measures elapsed time
since the discontinuation week itself, so at the first
off-drug measurement `tsd = 0`.

This produces a systematic offset of approximately one
measurement interval. The consequence is that the main
applies *more* carryover retention at the first off-drug
timepoint than the orig (because `(1/2)^0 = 1` vs.
`(1/2)^(sf * interval / t_half)`).

**Divergence 4: Scale factor.**
The orig's carryover formula includes `scalefactor = 2`
by default. The main's `simulation_clustered.R` uses a
direct `(1/2)^(tsd / t_half)` without scale factor.
This doubles the effective decay rate in the orig relative
to the main for the same nominal half-life.

### 9.7 Combined Impact on Power Estimates

These divergences compound. For the Hybrid design at
t½ = 0.2 weeks, the carryover retention at the first
off-drug timepoint (week 11) is:

- **orig:** `(1/2)^(2 * 1 / 0.2) = (1/2)^10 ≈ 0.001`
  (0.1% retention)
- **main:** `(1/2)^(0 / 0.2) = (1/2)^0 = 1.0`
  (100% retention)

This is a three-order-of-magnitude difference in
carryover effect at the most informative contrast point.
The main effectively treats the first off-drug measurement
as having full drug effect, while the orig treats it as
having nearly zero residual effect.

The implications for power estimation are:

1. The **main will show larger power dropoff** from
   carryover than the orig, because off-drug observations
   retain more drug effect and thus reduce the on/off
   contrast more severely.

2. The **orig's power estimates are less sensitive to
   carryover** because the first off-drug measurement
   already has minimal residual effect, preserving the
   contrast.

3. The two codebases will produce **non-comparable power
   curves** even for identical nominal parameters
   (same t½, same design, same sample size).

Whether the orig or main convention is "correct" depends
on the intended interpretation:

- If `tsd` should measure "how long has the participant
  been off drug by the time this measurement is taken,"
  the **orig's convention** is more natural -- the
  participant has been off drug for the full interval
  duration by measurement time.

- If `tsd` should measure "how long since the treatment
  state changed," the **main's convention** is more
  natural -- the transition happens at the week boundary
  and time accumulates from there.

The Hendrickson et al. (2020) paper should be consulted
to determine which convention matches the published
methodology.

---

## 10. Recommendations for pmsimstats2025

The goal of pmsimstats2025 is to replicate the orig's
simulation methodology as closely as possible using modern
R practices, and then to diagnose why the orig shows power
falling too rapidly as carryover increases. The
recommendations below are ordered to achieve faithful
replication first (items 1-4), then enable diagnosis of
the power dropoff (items 5-6).

### 10.1 Replication: Fix tsd convention (Critical)

The main's `tsd = week - discontinuation_week` produces
`tsd = 0` at the first off-drug timepoint. The orig's
interval-accumulation approach produces `tsd = interval
duration` (1-2.5 weeks) at the first off-drug timepoint.
This is the single largest source of divergence.

Recommended fix in `simulation_clustered.R`:

```r
# Replace current tsd computation (lines 658-665)
# with last-on-drug-week anchor:
last_on_drug_week = if_else(
  treatment == 1, week, NA_real_),
last_on_drug_week = zoo::na.locf(
  last_on_drug_week, na.rm = FALSE),
tsd = if_else(
  treatment == 0 & !is.na(last_on_drug_week),
  week - last_on_drug_week,
  0
)
```

This anchors tsd to the last *on-drug* measurement week
rather than the first *off-drug* measurement week,
producing values that match the orig:

```
Hybrid Path A:
  Week:       4   8   9  10  11  12  16  20
  Current:    0   0   0   0   0   1   0   0
  Fixed:      0   0   0   0   1   2   0   4
  Orig:       0   0   0   0   1   2   0   4  ✓
```

### 10.2 Replication: Time-weighted tod with gap reset (Critical)

Replace `cumsum(treatment)` with time-weighted accumulation
that resets after off-drug gaps, matching the orig's `tod`
behavior.

The orig's `tod` at timepoint i equals:
- On drug: sum of interval durations for all consecutive
  on-drug timepoints ending at i
- Off drug: 0

This means drug re-initiation restarts the exposure clock.
At week 16 of Hybrid Path A, `tod = 4` (the 4-week
interval since week 12), not `14` (total lifetime exposure).

Recommended approach: compute `tod` with the same loop
logic as the orig, or equivalently:

```r
interval = week - lag(week, default = 0),
drug_interval = treatment * interval,
# Reset-aware cumulative sum: restarts when
# drug_interval is preceded by 0
tod = ave(drug_interval, cumsum(drug_interval == 0),
          FUN = cumsum)
```

This also changes the BR mean computation. The orig uses
`modgompertz(tod, ...)` while the main uses
`weeks_on_drug * BR_rate`. To replicate the orig, the
main should either:
- (a) Implement the modified Gompertz and feed it `tod`,
  matching the orig's DGP exactly, or
- (b) Use `tod * BR_rate` as a linear approximation,
  accepting that the dose-response shape differs

Option (a) is more faithful; option (b) is simpler and
may suffice if the Gompertz is approximately linear over
the trial duration.

### 10.3 Replication: Restore scalefactor = 2 (Critical)

The orig's carryover formula is:

```r
brmeans[p] + brmeans[p-1] *
  (1/2)^(scalefactor * tsd / t1half)
```

with `scalefactor = 2`. The main's
`simulation_clustered.R` omits the scale factor:

```r
(1/2)^(tsd / carryover_t_half)
```

This halves the effective decay rate in the main relative
to the orig. To replicate the orig, restore
`scalefactor = 2`:

```r
carryover_effect = case_when(
  treatment == 1 ~ 1,
  treatment == 0 & carryover_t_half == 0 ~ 0,
  treatment == 0 & carryover_t_half > 0 ~
    (1/2)^(2 * tsd / carryover_t_half),
  TRUE ~ 0
)
```

Note: `pm_functions.R` already has `scale_factor = 2` as
default. The divergence is only in `simulation_clustered.R`
which implements its own inline carryover computation.

### 10.4 Replication: Implement Ron Thomas correlation scaling (Critical)

The orig's `generateData.R` (lines 148-153) scales the
biomarker-BR correlation at off-drug timepoints using a
ratio of BR means:

```r
correlations["bm", n1] <- ifelse(
  brtest[p],
  ifelse(brmeans[p] == 0, 0,
         (mm1 / mm0) * modelparam$c.bm),
  modelparam$c.bm)
```

where `brtest[p]` is TRUE when the raw (pre-carryover)
BR mean was zero, `mm1` is the carryover-adjusted BR
mean at the current timepoint, and `mm0` is the BR mean
at the previous timepoint.

This scaling is the mechanism by which the biomarker-
treatment interaction varies across on-drug and off-drug
timepoints in the MVN covariance structure. The main's
`simulation_clustered.R` does not implement this scaling
-- it uses a fixed `c.bm` for all BR-biomarker
correlations regardless of treatment status. This is a
fundamental methodological difference that must be
resolved for faithful replication.

The interaction between this correlation scaling and
the carryover parameters is the most likely source of
the "power falls too rapidly" problem (see Section 10.5).

### 10.5 Diagnosis: Correlation matrix perturbation analysis (High)

Once the main replicates the orig's variable chain, the
power dropoff can be diagnosed by examining how the
biomarker-BR correlation varies across carryover levels.

With the orig's aggressive decay (`scalefactor = 2`,
`t_half = 0.2`, `tsd >= 1` at first off-drug timepoint):

```
Carryover retention at week 11: 0.001
Correlation scaling: mm1/mm0 ≈ 0.001
Effective c.bm at week 11: 0.001 * 0.3 ≈ 0.0003
```

This near-zero correlation at off-drug timepoints means
the off-drug observations contribute almost no information
about the biomarker-treatment interaction. The result is
effectively the same as having no off-drug data at all.

But with zero carryover (`t_half = 0`), the off-drug BR
mean is exactly zero and the correlation is set to zero.
The off-drug observations also contribute no interaction
information.

So the question becomes: **why does power differ between
`t_half = 0` and `t_half = 0.2` if both produce near-zero
off-drug correlations?**

Possible mechanisms to investigate:

1. **PD enforcement distortion.** The
   `make.positive.definite(sigma, tol = 1e-3)` step may
   redistribute variance differently when small nonzero
   correlations are present (carryover > 0) versus exactly
   zero (carryover = 0). Even small perturbations to the
   on-drug block of the covariance matrix could degrade
   the detectable interaction.

2. **BR mean contamination.** Even 0.1% residual BR at
   off-drug timepoints changes the Gompertz trajectory
   through the accumulation logic. At week 12 of the
   Hybrid design, `brmeans[6] = brmeans[5] * 0.001`, but
   `brmeans[5]` itself was `brmeans[4] * 0.001`. These
   cascading small values may interact with the
   correlation scaling in unexpected ways (e.g.,
   `mm1/mm0 = 0.001/0.001` producing unstable ratios).

3. **Analysis model sensitivity.** The analysis model uses
   `Dbc` (the continuous treatment indicator), which in
   the orig incorporates the same carryover formula. When
   `Dbc` has small nonzero values off-drug instead of
   exactly zero, the `bm * Dbc` interaction has slightly
   different properties. The regression coefficient may
   be biased by the near-zero off-drug `Dbc` values that
   still carry some variance.

Recommended diagnostic script:

```r
# For each carryover level (0, 0.1, 0.2):
# 1. Print the full 26x26 correlation matrix
# 2. Print eigenvalues before and after PD enforcement
# 3. Print the biomarker-BR correlations at each TP
# 4. Print the BR means at each TP per path
# 5. Compare sigma matrices element-wise between
#    t_half = 0 and t_half > 0
```

### 10.6 Diagnosis: Targeted power comparison (High)

Run a small simulation (50-100 iterations) with the
replicated orig logic at three carryover levels and
compare power curves. If the main's replicated results
match the orig's published power curves, the replication
is confirmed and the power dropoff is inherent to the
orig's methodology. If they still differ, there is a
remaining implementation gap to close.

### 10.7 Supporting changes

| # | Action                                  | Priority |
|---|-----------------------------------------|:--------:|
| 7 | Validation script: print all 5 vars     | High     |
|   | for each design/path against hand-       |          |
|   | computed reference values                |          |
| 8 | Unit tests for tsd, tod, carryover at    | Medium   |
|   | boundary conditions                      |          |
| 9 | Document all convention choices           | Medium   |
|   | (tsd anchor, tod reset, scale factor,    |          |
|   | correlation scaling) in a single         |          |
|   | reference document                       |          |
| 10| Trim DESCRIPTION: remove unused          | Low      |
|   | database drivers, palmerpenguins, etc.   |          |

---

*Rendered on 2026-03-18 at 09:25 PDT.*
*Source: ~/prj/alz/pmsim-power-dropoff-analysis.md*
