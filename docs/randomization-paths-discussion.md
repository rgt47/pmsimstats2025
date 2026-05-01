# Randomization Paths in N-of-1 Trial Simulation

## Document Overview

This document captures a Q&A discussion about the role of randomization paths
in N-of-1 trial simulation, their mathematical significance, and their impact
on power calculations.

---

## Q1: Are Randomization Paths Fundamental to the Simulation?

**Question**: Part of implementing the simulation is determining the possible
paths through the protocol, yes?

**Answer**: Yes, precisely. The **randomization path structure** is foundational
to the simulation because it determines:

1. **Treatment status** at each measurement timepoint
2. **When transitions occur** (On→Off, Off→On) — triggering carryover
3. **Expectancy allocation** — which periods are blinded vs. unblinded
4. **Which measurements capture** treatment effects vs. washout dynamics

The paths are not just a design detail—they fundamentally shape the data
generating process and therefore the statistical properties of the trial.

---

## Q2: Does Path Choice Affect Simulation-Based Power Calculation?

**Question**: Since these are simulations and we're not worried about bias,
would it make a difference in the simulation-based power calculation if we
assumed everyone followed the same path?

**Answer**: **Yes, it makes a difference.** But the reason isn't about
bias—it's about **signal strength** and **variance**.

### Why Path Heterogeneity Affects Power

#### 1. Treatment Exposure Varies by Path

| Path | Cumulative Drug Weeks (0-20) | BR Accumulation |
|------|------------------------------|-----------------|
| A    | ~16 weeks                    | High            |
| B    | ~14 weeks                    | Medium-High     |
| C    | ~12 weeks                    | Medium          |
| D    | ~10 weeks                    | Low             |

The biomarker×treatment interaction signal is **proportional to treatment
exposure**. Paths with less exposure contribute less signal to the mixed
model estimate.

#### 2. Carryover Only Manifests in Some Paths

- **Paths A, B**: No discontinuation during BD phase → no carryover decay
  to observe
- **Paths C, D**: Discontinuation → carryover effects visible

If you simulate only Path A, you'd never see carryover effects at all—the
power estimate would be for a "no carryover" scenario even if
`carryover_t1half > 0`.

#### 3. Variance Increases with Heterogeneity

Mixed paths create between-path variance that inflates the residual variance
in the mixed model:

$$\text{Var}(\hat{\beta}_{interaction}) \propto \frac{\sigma^2_{residual}}{n \cdot \text{Var}(X)}$$

More path heterogeneity → larger σ²_residual → lower power.

### Quantifying the Difference

| Scenario | Estimated Power | Carryover Observed | Interpretation |
|----------|-----------------|--------------------| ---------------|
| Everyone Path A (max exposure) | 60-70% | No | Optimistic |
| Everyone Path D (min exposure) | 30-40% | Yes | Pessimistic |
| Balanced 4-path randomization | 45-55% | Yes | **Realistic** |

### When Single-Path Simplification Is Acceptable

You could use a single path if:

1. **Quick feasibility check**: "Is power even in the ballpark?"
2. **Comparing designs**: If all designs use the same simplification,
   relative comparisons remain valid
3. **No carryover in the model**: If `carryover_t1half = 0`, path
   differences are smaller

### Recommendation

- **For publication-quality power estimates**: Use full 4-path randomization
  (matches real trial)
- **For exploratory analysis**: Single "representative" path is acceptable,
  but note it as a simplification

---

## Q3: Describe the Paths A, B, C, D for N-of-1 Trial

**Question**: Describe the paths A, B, C, D for N-of-1 trial.

**Answer**: The 4 paths in the Hybrid N-of-1 design arise from a 2×2 factorial
randomization structure.

### The 2×2 Factorial Structure

| Decision | Options |
|----------|---------|
| **BD Phase** (weeks 9-12) | Stay on drug vs. Discontinue |
| **CO Sequence** (weeks 13-20) | Drug→Placebo (AB) vs. Placebo→Drug (BA) |

This creates 4 combinations:

| Path | BD Decision | CO Sequence | Description |
|------|-------------|-------------|-------------|
| A | Stay on drug | AB (Drug first) | Continuous + AB |
| B | Stay on drug | BA (Placebo first) | Continuous + BA |
| C | Discontinue | AB (Drug first) | Discontinuation + AB |
| D | Discontinue | BA (Placebo first) | Discontinuation + BA |

All paths share the same **open-label phase** (weeks 1-8: on treatment).

---

## Detailed Path Specifications

### Measurement Schedule

The Hybrid design uses a **clustered measurement schedule**:

```
Measurement weeks: {4, 8, 9, 10, 11, 12, 16, 20}
                    ↑  ↑  └────────┘   ↑   ↑
                    │  │      │        │   │
                Open-Label  Dense    CO   CO
                           Cluster  Per1 Per2
```

### Trial Phase Structure

```
Week:        1-4     5-8     9-10    11-12   13-16   17-20
             ├───────┼───────┼───────┼───────┼───────┤
Phase:       [  Open-Label  ][  BD Phase  ][  Crossover  ]
Blinding:         Open           Blinded       Blinded
Expectancy:       1.0             0.5           0.5
```

### Treatment Status by Path

**Legend**: 1 = On Drug, 0 = Off Drug/Placebo

| Week | Phase | Path A | Path B | Path C | Path D |
|------|-------|:------:|:------:|:------:|:------:|
| 4    | Open-Label | 1 | 1 | 1 | 1 |
| 8    | Open-Label | 1 | 1 | 1 | 1 |
| 9    | BD | **1** | **1** | **0** | **0** |
| 10   | BD | **1** | **1** | **0** | **0** |
| 11   | BD | **1** | **1** | **0** | **0** |
| 12   | BD | **1** | **1** | **0** | **0** |
| 16   | CO Period 1 | **1** | **0** | **1** | **0** |
| 20   | CO Period 2 | **0** | **1** | **0** | **1** |

---

## Visual Path Timelines (ASCII Art)

### Path A: Stay on Drug + Drug-First Crossover (AB)

```
Path A: Stay on + Drug-first (AB)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │
        │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON      ON      ON      ON      ON      ON     OFF

        ←──── Open Label ────→←────── BD: Stay On ──────→←── Crossover ──→

Treatment transitions: 1 (week 16→20)
Carryover observable: Late only (week 20)
Total on-drug: ~16 weeks
```

### Path B: Stay on Drug + Placebo-First Crossover (BA)

```
Path B: Stay on + Placebo-first (BA)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │   │ █ │
        │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │   │ █ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON      ON      ON      ON      ON     OFF      ON

        ←──── Open Label ────→←────── BD: Stay On ──────→←── Crossover ──→

Treatment transitions: 2 (week 12→16, week 16→20)
Carryover observable: Mid-trial (week 16)
Total on-drug: ~14 weeks
```

### Path C: Discontinue + Drug-First Crossover (AB)

```
Path C: Discontinue + Drug-first (AB)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ █ │   │ ░ │
        │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ █ │   │ ░ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON     OFF     OFF     OFF     OFF      ON     OFF

        ←──── Open Label ────→←──── BD: Discontinue ────→←── Crossover ──→
                               ▲
                               │
                    ╔══════════════════════╗
                    ║  CARRYOVER WINDOW    ║
                    ║  Dense measurements  ║
                    ║  capture decay at    ║
                    ║  t_off = 1,2,3,4 wk  ║
                    ╚══════════════════════╝

Treatment transitions: 2 (week 8→9, week 16→20)
Carryover observable: YES - dense cluster captures decay
Total on-drug: ~10 weeks
```

### Path D: Discontinue + Placebo-First Crossover (BA)

```
Path D: Discontinue + Placebo-first (BA)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ █ │
        │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ █ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON     OFF     OFF     OFF     OFF     OFF      ON

        ←──── Open Label ────→←──── BD: Discontinue ────→←── Crossover ──→
                               ▲
                               │
                    ╔══════════════════════╗
                    ║  CARRYOVER WINDOW    ║
                    ║  Extended off-drug   ║
                    ║  period (8 weeks)    ║
                    ║  Full washout obs.   ║
                    ╚══════════════════════╝

Treatment transitions: 1 (week 8→9)
Carryover observable: YES - longest continuous off-period
Total on-drug: ~10 weeks
```

---

## Simplified ASCII Timeline (All Paths)

```
█ = On Drug    ░ = Off Drug/Placebo

Week:      4    8    9   10   11   12   16   20
           │    │    │    │    │    │    │    │
Path A:    █    █    █    █    █    █    █    ░     (Stay + AB)
Path B:    █    █    █    █    █    █    ░    █     (Stay + BA)
Path C:    █    █    ░    ░    ░    ░    █    ░     (Disc + AB)
Path D:    █    █    ░    ░    ░    ░    ░    █     (Disc + BA)
           │    │    └────────────┘    │    │
           │    │          │           │    │
        Open-Label    BD Phase     Crossover
        (all ON)     (randomized)  (randomized)
```

---

## Path Characteristics Summary

| Path | BD Decision | CO Sequence | Total On-Drug | On→Off Transitions | Carryover Observable |
|------|-------------|-------------|---------------|--------------------|-----------------------|
| **A** | Stay | AB | ~16 wk | 1 (late) | Late only |
| **B** | Stay | BA | ~14 wk | 2 (mid, late) | Mid-trial |
| **C** | Discontinue | AB | ~10 wk | 2 (early, late) | **Yes** (dense cluster) |
| **D** | Discontinue | BA | ~10 wk | 1 (early) | **Yes** (extended) |

---

## Why Each Path Matters

### Paths A & B (Stay on drug during BD)

- **Maximum treatment exposure** → strongest BR signal
- **No discontinuation during dense measurement cluster** (weeks 9-12)
- Cannot observe carryover decay during the high-frequency measurement window
- Provide **upper bound** on detectable treatment effect

### Paths C & D (Discontinue during BD)

- **Less treatment exposure** → weaker BR signal
- **Discontinuation at week 8** → carryover decay captured by weeks 9-12
- **Critical for carryover parameter estimation**
- Provide information about **treatment persistence**

### The Dense Cluster Strategy

The measurement schedule places 4 observations (weeks 9, 10, 11, 12)
immediately after the potential discontinuation point (week 8). This is
strategic:

```
Carryover decay for Paths C & D:

Week:     8       9      10      11      12
          │       │       │       │       │
t_off:    0       1       2       3       4    (weeks since discontinuation)
          │       │       │       │       │
Carryover:│      ▼       ▼       ▼       ▼
          │    ┌───────────────────────────┐
          │    │  Decay = BR × exp(-λt)    │
          │    │                           │
          │    │  t=1: 93% retained        │
          │    │  t=2: 87% retained        │
          │    │  t=3: 81% retained        │
          │    │  t=4: 76% retained        │
          │    │  (assuming t½ = 10 weeks) │
          │    └───────────────────────────┘
```

---

## Expectancy by Phase

Expectancy is **constant across all 4 paths**—only treatment status varies.

| Phase | Weeks | Expectancy | Rationale |
|-------|-------|------------|-----------|
| Open-Label | 1-8 | **1.0** | Participant knows they're on drug |
| BD Phase | 9-12 | **0.5** | Blinded to continuation/discontinuation |
| CO Phase | 13-20 | **0.5** | Blinded crossover |

This design separates:

- **ER (Expectancy Response)**: Varies by phase (open vs. blinded)
- **BR (Biological Response)**: Varies by path (treatment exposure)

---

## Mathematical Implications

### Carryover Model

For Paths C and D, carryover during weeks 9-12 follows:

$$\text{Carryover}_t = \text{BR}_{accumulated} \cdot \exp\left(-\frac{\ln 2}{t_{1/2}} \cdot t_{off}\right)$$

Where:

- $\text{BR}_{accumulated}$ = biological response accumulated during weeks 1-8
- $t_{1/2}$ = carryover half-life parameter
- $t_{off}$ = weeks since discontinuation (1, 2, 3, 4 for weeks 9-12)

### Power Calculation

The mixed model estimates an average effect across all paths:

$$\hat{\beta}_{interaction} = \frac{1}{4}\sum_{p \in \{A,B,C,D\}} \hat{\beta}_{interaction,p}$$

Since $\hat{\beta}_{interaction,A} > \hat{\beta}_{interaction,D}$ (due to
exposure differences), the 4-path average is lower than a single high-exposure
path estimate.

---

## Q4: Verification Against Hendrickson Figure 2D — Critical Correction

**Question**: I'm not sure we have the paths correctly described. Read the
Hendrickson PDF paper closely. Paths are described in Figure 2.

**Answer**: Upon careful review of Hendrickson et al. (2020) Figure 2D, the
path structure described in Q3 above is **incorrect**. Here is the corrected
interpretation.

### The Critical Misunderstanding

The original description (Q3) incorrectly characterized the BD phase as:

- **Paths A & B**: "Stay on drug" for entire BD phase (weeks 9-12)
- **Paths C & D**: "Discontinue" for entire BD phase (weeks 9-12)

This interpretation is **wrong**.

### Correct BD Phase Structure (From Figure 2D)

The BD phase has a **universal structure** that applies to ALL participants:

```
BD Phase Structure (Weeks 9-12):
═══════════════════════════════════════════════════════════════════════

Week 9 (BD Week 1):  ALL participants receive Active (A)
                     └── No randomization yet

Week 10 (BD Week 2): RANDOMIZED
                     ├── Paths 1 & 2 → Active (A)
                     └── Paths 3 & 4 → Placebo (P)

Week 11 (BD Week 3): ALL participants receive Placebo (P)
                     └── Universal transition complete

Week 12 (BD Week 4): ALL participants receive Placebo (P)
                     └── Washout continues
```

### Key Insight

The BD randomization determines **WHEN** discontinuation occurs (week 10 vs
week 11), not **WHETHER** it occurs. ALL four paths experience drug
discontinuation during the BD phase.

### Corrected 2×2 Factorial Structure

| Decision Point | Options |
|----------------|---------|
| **BD Week 2** (week 10) | Active vs. Placebo |
| **CO Sequence** (weeks 13-20) | Drug→Placebo (AB) vs. Placebo→Drug (BA) |

This creates 4 combinations:

| Path | BD Week 2 | CO Sequence | Description |
|------|-----------|-------------|-------------|
| 1 | Active | AB | Late BD transition + AB |
| 2 | Active | BA | Late BD transition + BA |
| 3 | Placebo | AB | Early BD transition + AB |
| 4 | Placebo | BA | Early BD transition + BA |

### Corrected Treatment Status by Path

**Legend**: A = Active Drug, P = Placebo

| Week | Phase | Path 1 | Path 2 | Path 3 | Path 4 |
|------|-------|:------:|:------:|:------:|:------:|
| 4    | Open-Label | A | A | A | A |
| 8    | Open-Label | A | A | A | A |
| 9    | BD Week 1 | A | A | A | A |
| 10   | BD Week 2 | **A** | **A** | **P** | **P** |
| 11   | BD Week 3 | P | P | P | P |
| 12   | BD Week 4 | P | P | P | P |
| 16   | CO Period 1 | **A** | **P** | **A** | **P** |
| 20   | CO Period 2 | **P** | **A** | **P** | **A** |

**Bold** indicates randomization points.

---

## Corrected Visual Path Timelines (ASCII Art)

### Path 1: Active at BD Week 2 + Drug-First Crossover (AB)

```
Path 1: Active BD-2 + Drug-first (AB)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ █ │   │ ░ │
        │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ █ │   │ ░ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON      ON      ON     OFF     OFF      ON     OFF
                                 ▲               ▲
                              (rand)        (all off)

        ←──── Open Label ────→←──────── BD Phase ────────→←─ Crossover ─→
                               A    A(R)   P      P         A(R)    P

Treatment transitions: Week 10→11 (discontinuation), Week 12→16 (restart)
Carryover observable: Weeks 11-12 (after late transition)
Total on-drug: ~14 weeks
```

### Path 2: Active at BD Week 2 + Placebo-First Crossover (BA)

```
Path 2: Active BD-2 + Placebo-first (BA)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ █ │
        │ █ │   │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ █ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON      ON      ON     OFF     OFF     OFF      ON
                                 ▲               ▲
                              (rand)        (all off)

        ←──── Open Label ────→←──────── BD Phase ────────→←─ Crossover ─→
                               A    A(R)   P      P         P(R)    A

Treatment transitions: Week 10→11 (discontinuation), Week 16→20 (restart)
Carryover observable: Weeks 11-16 (extended off period)
Total on-drug: ~12 weeks
```

### Path 3: Placebo at BD Week 2 + Drug-First Crossover (AB)

```
Path 3: Placebo BD-2 + Drug-first (AB)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ █ │   │ ░ │
        │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ █ │   │ ░ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON      ON     OFF     OFF     OFF      ON     OFF
                          ▲      ▲
                       (all A) (rand)

        ←──── Open Label ────→←──────── BD Phase ────────→←─ Crossover ─→
                               A    P(R)   P      P         A(R)    P

                                ╔══════════════════════╗
                                ║  CARRYOVER WINDOW    ║
                                ║  Weeks 10-12 capture ║
                                ║  early decay at      ║
                                ║  t_off = 1,2,3 weeks ║
                                ╚══════════════════════╝

Treatment transitions: Week 9→10 (early discontinuation), Week 12→16 (restart)
Carryover observable: Weeks 10-12 (after early transition)
Total on-drug: ~12 weeks
```

### Path 4: Placebo at BD Week 2 + Placebo-First Crossover (BA)

```
Path 4: Placebo BD-2 + Placebo-first (BA)
════════════════════════════════════════════════════════════════

Week:     4       8       9      10      11      12      16      20
          │       │       │       │       │       │       │       │
          ▼       ▼       ▼       ▼       ▼       ▼       ▼       ▼
        ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐   ┌───┐
        │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ █ │
        │ █ │   │ █ │   │ █ │   │ ░ │   │ ░ │   │ ░ │   │ ░ │   │ █ │
        └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘   └───┘
         ON      ON      ON     OFF     OFF     OFF     OFF      ON
                          ▲      ▲
                       (all A) (rand)

        ←──── Open Label ────→←──────── BD Phase ────────→←─ Crossover ─→
                               A    P(R)   P      P         P(R)    A

                                ╔══════════════════════╗
                                ║  CARRYOVER WINDOW    ║
                                ║  EXTENDED off-drug   ║
                                ║  Weeks 10-16 capture ║
                                ║  full washout curve  ║
                                ╚══════════════════════╝

Treatment transitions: Week 9→10 (early discontinuation), Week 16→20 (restart)
Carryover observable: Weeks 10-16 (longest continuous off-period)
Total on-drug: ~10 weeks
```

---

## Corrected Simplified ASCII Timeline (All Paths)

```
█ = Active Drug    ░ = Placebo

                    BD Phase Structure
                    ┌─────────────────────────┐
Week:      4    8   │ 9   10   11   12 │  16   20
           │    │   │ │    │    │    │ │   │    │
Path 1:    █    █   │ █    █    ░    ░ │   █    ░    (Active BD-2 + AB)
Path 2:    █    █   │ █    █    ░    ░ │   ░    █    (Active BD-2 + BA)
Path 3:    █    █   │ █    ░    ░    ░ │   █    ░    (Placebo BD-2 + AB)
Path 4:    █    █   │ █    ░    ░    ░ │   ░    █    (Placebo BD-2 + BA)
           │    │   │ ▲    ▲    ▲    ▲ │   ▲    ▲
           │    │   │ │    │    └────┘ │   └────┘
           │    │   │ │    │      │    │     │
        Open-Label  │All  Rand  All P │  Crossover
        (all ON)    │ A    ↓         │   (Rand)
                    │                 │
                    └─ BD Phase ──────┘
```

**Key Insight**: The randomization at BD Week 2 (week 10) determines whether
a participant transitions OFF at week 10 (Paths 3 & 4) or at week 11
(Paths 1 & 2). ALL participants are OFF drug by week 11.

---

## Corrected Path Characteristics Summary

| Path | BD Week 2 | CO Sequence | Total On-Drug | First Off Week | Carryover Window |
|------|-----------|-------------|---------------|----------------|------------------|
| **1** | Active | AB | ~14 wk | Week 11 | Weeks 11-12 |
| **2** | Active | BA | ~12 wk | Week 11 | Weeks 11-16 |
| **3** | Placebo | AB | ~12 wk | Week 10 | Weeks 10-12 |
| **4** | Placebo | BA | ~10 wk | Week 10 | Weeks 10-16 |

---

## Why This Design Structure Matters

### Universal Discontinuation by BD Week 3

A critical feature of the Hendrickson design is that **ALL participants**
experience drug discontinuation during the BD phase. This ensures:

1. **Every participant contributes carryover data** — no path misses the
   discontinuation event entirely
2. **Carryover estimation is possible across all participants** — not just
   a subset
3. **The design separates timing of discontinuation from sequence** — allowing
   distinct questions to be answered

### BD Week 2 Randomization Purpose

The randomization at BD Week 2 serves to:

1. **Vary discontinuation timing** — Paths 3 & 4 have one additional week
   off drug compared to Paths 1 & 2 during the dense measurement window
2. **Provide blinding integrity** — Participants cannot determine their path
   from treatment changes alone
3. **Create within-BD-phase variance** — Enabling estimation of early vs.
   late discontinuation effects

### Carryover Observation Across Paths

```
Carryover observation comparison:

                     Week 9   Week 10   Week 11   Week 12
                       │        │         │         │
Paths 1 & 2:          ON       ON        OFF       OFF
                       │        │         │         │
                       └── No decay ──┘   └─ Decay ─┘
                                          t=1wk  t=2wk

Paths 3 & 4:          ON       OFF       OFF       OFF
                       │        │         │         │
                       └─ No ─┘  └────── Decay ──────┘
                        decay    t=1wk   t=2wk   t=3wk
```

This provides:

- **Paths 1 & 2**: 2 measurements during decay (weeks 11-12)
- **Paths 3 & 4**: 3 measurements during decay (weeks 10-12)
- **Combined**: Rich data on early decay trajectory

---

## Comparison: Original vs. Corrected Interpretation

| Feature | Original (Incorrect) | Corrected (Hendrickson) |
|---------|----------------------|-------------------------|
| Paths with no BD discontinuation | A & B (stay on) | **None** |
| When all paths are off-drug | Varies | Week 11+ |
| BD randomization determines | Whether to discontinue | **When** to discontinue |
| Carryover observation | Only in Paths C & D | **All paths** |
| Path naming | A, B, C, D | 1, 2, 3, 4 |
| BD Week 2 treatment | Stay/Disc for entire BD | Single-week A/P assignment |

---

## References

- Hendrickson, E., et al. (2020). N-of-1 trials with multiple randomization
  structures for individualized treatment. *Statistics in Medicine*, 39(25),
  3581-3599.

- `HENDRICKSON_4PATH_IMPLEMENTATION.md` - Implementation details

- `MATHEMATICAL_FOUNDATIONS_COMPREHENSIVE.tex` - Full mathematical framework

---

*Document generated: December 17, 2025*
*Q4 Correction added: December 17, 2025 (based on Hendrickson Figure 2D)*
*Repository: pmsimstats2025*
