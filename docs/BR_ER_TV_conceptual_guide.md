# Understanding the BR-ER-TV Decomposition
## An Intuitive Guide to Treatment Response in N-of-1 Trials

*March 19, 2026 at 17:54*

---

## The Big Picture: Why Split Responses?

Imagine you're a PTSD patient trying prazosin (a blood pressure medication that also reduces nightmares). After starting treatment, your PTSD symptoms improve. But **why** did they improve?

Total improvement = Drug effect + Expectation effect + Natural healing

$Y = BR + ER + TV + noise$

This is the fundamental question: your overall response is a mixture of multiple causes. If we lump them together, we lose information about **what's actually happening**.

### The Problem with One-Component Models

A naive model might say: `improvement ~ treatment + error`

But this commits a fundamental mistake: the `error` term becomes a wastebasket for everything we didn't measure.

```
improvement = 10 (treatment) + 15 (expectation) + 3 (natural healing) + 2 (noise) = 30 total

Naive model sees: improvement = 10 (treatment) + 20 (error)
```

The model gets the right total but **completely misses the structure**. You think the drug is more powerful than it is, because you've conflated drug effect with placebo effect.

---

## Component 1: Biological Response (BR)

### What It Is

**BR is the physiological response to the drug itself.** It's what happens because of the medication's chemical effect on your body, independent of what you believe or expect.

### Real Example: PTSD and Prazosin

Prazosin is an alpha-blocker originally developed for high blood pressure. It also reduces nightmares in PTSD patients—but only when it's actually in their system.

- **Week 1 on drug**: Brain prazosin concentration increases → nightmares decrease by 5 points
- **Week 2 on drug**: Continued exposure → effect plateaus at 8 points (brain is saturated)
- **Week 3 on drug**: Effect stays at 8 points (steady-state)
- **Week 4 OFF drug**: Within days, prazosin washes out → nightmare symptoms return to baseline

This is BR: the pharmacological effect that depends only on drug exposure, not on belief.

### Mathematical Pattern

BR follows a **saturating curve** (modified Gompertz):

$$BR(t) = \text{max\_BR} \times (1 - \exp(-\text{rate} \times (1 + \text{disp})^{-\text{time\_on\_drug}}))$$

**In English**:
- Starts at zero (no drug = no effect)
- Increases rapidly at first (each dose matters a lot)
- Levels off (more drug doesn't help more—law of diminishing returns)
- Never exceeds `max_BR` (biological ceiling)

### Why Not Just Use a Linear Trend?

Linear would assume: *each week on drug adds the same amount of improvement*

$$BR = 2 \text{ points/week} \times \text{weeks\_on\_drug}$$

After 1 week: 2 points
After 2 weeks: 4 points
After 10 weeks: 20 points ← Biologically implausible!

Real drugs don't work linearly. Your body reaches saturation. The Gompertz curve captures this.

### Key Parameters

| Parameter | Meaning | Intuition |
|-----------|---------|-----------|
| `max_BR` | Ceiling effect | "Best this drug can possibly do" (e.g., 15 points max improvement) |
| `rate` | Speed of onset | How quickly you reach half-maximum (fast = 1 week, slow = 4 weeks) |
| `disp` | Shape of curve | Whether improvement is early (steep) or gradual (flat) |

**Concrete example**: Prazosin for PTSD
- `max_BR = 10` (drug can reduce nightmares by at most 10 points on a 0-100 scale)
- `rate = 0.3` (takes about 3-4 weeks to reach half-maximum effect)
- `disp = 1.0` (smooth, classical saturation curve)

---

## Component 2: Expectancy Response (ER)

### What It Is

**ER is the placebo effect—how much you improve just because you believe you're being treated.**

This is NOT "all in your head" (a dismissive phrase). The placebo effect is neurobiologically real: belief activates dopamine, reduces inflammation, and decreases pain perception.

### Real Example: The Blinded Discontinuation Design

Here's where it gets clever. In the Hendrickson design, patients don't know when they're switched to placebo.

- **Weeks 1-8 (Open-label)**: Everyone knows they're on prazosin
  - BR effect: 8 points (biological)
  - ER effect: 4 points (expectation: "I'm on the drug, so I should improve")
  - **Total: 12 points**

- **Weeks 9-12 (Blinded discontinuation)**: Half switch to placebo *without knowing*
  - Drug group: BR = 8 (still on drug), ER = 2 (lower expectation: "might be placebo now")
  - Placebo group: BR = 0 (not on drug), ER = 3 (expectation: "might still be on drug")
  - **Total for drug: 10 points**
  - **Total for placebo: 3 points**

The 7-point gap is where the real drug effect reveals itself.

### Why ER Needs Its Own Component

If ER was just absorbed into "error", you'd miss this pattern entirely.

```
Naive model sees:
- Open label: improvement = 12
- Blinded: improvement = 10 vs 3
- Conclusion: "drug is good, people are a bit noisy"

Smart model sees:
- Open label: BR=8 + ER=4
- Blinded on-drug: BR=8 + ER=2
- Blinded placebo: BR=0 + ER=3
- Conclusion: "drug effect is 8, expectancy boosts it by 4, test is sensitive"
```

### The Expectancy Parameter

ER also follows a Gompertz curve, but **scaled by expectancy**:

$$ER(t) = \text{max\_ER} \times \text{expectancy\_factor} \times (1 - \exp(-\ldots))$$

where `expectancy_factor` varies by phase:

| Phase | Expectancy | Reasoning |
|-------|-----------|-----------|
| Open-label "on drug" | 1.0 | Full belief in treatment |
| Blinded discontinuation | 0.5 | Uncertain if still on drug |
| Open-label placebo | 0.5 | Told it's placebo or washout |

**Why 0.5 not 0?** Because even when told "you might be on placebo," patients retain partial expectation: "maybe I'm in the lucky group" or "the previous dose still helps."

### Why Not Constant Variance?

A naive model might assume: ER variance is the same regardless of phase.

But that's wrong! When expectancy is high (open-label, 1.0), outcomes are more variable:
- Some patients respond hugely to placebo (get 8 points)
- Others don't respond at all (get 0 points)
- Variance = 4-5 points

When expectancy is low (blinded, 0.5), outcomes are more uniform:
- Weaker effect, less room for individual variation
- Most patients get small improvements
- Variance = 1-2 points

**This is why ER variance is expectancy-scaled**: $SD_{ER} = \text{baseline\_SD} \times \text{expectancy}$

---

## Component 3: Time-Varying Response (TV)

### What It Is

**TV is natural disease progression or decline, independent of any treatment.** It's what would happen if the patient did absolutely nothing.

### Real Example: PTSD Symptom Trajectory

Some PTSD patients naturally improve over time (weeks and months of experience, talking to friends, life circumstances change). Some naturally worsen (stress accumulates, trauma anniversaries, loss).

- **Patient A**: PTSD symptoms naturally decline by 2 points per week (therapy is helping, life is improving)
- **Patient B**: PTSD symptoms naturally increase by 0.5 points per week (anniversary of trauma, new stress)

This has **nothing to do with the drug**. It's just the patient's natural trajectory.

### Why This Matters (And Why You Can't Absorb It Into Error)

Imagine we didn't include TV:

```
Observed response = BR + ER + TV + error

What we'd measure:
Patient A (improving naturally): response = 8 + 4 + 3 + noise = 15
Patient B (worsening naturally): response = 8 + 4 - 1 + noise = 11

Naive model: response ~ treatment
- Patient A looks like treatment is very effective
- Patient B looks like treatment is less effective
- Conclusion: maybe the drug works better for some people? (FALSE!)
```

The issue: **TV is confounded with treatment effect**. If you don't account for it, you'll think the drug is more powerful (or less powerful) than it actually is.

Throwing it in "error" doesn't help—the error term becomes biased and heteroscedastic (non-constant variance).

### Mathematical Pattern

TV also follows a Gompertz curve:

$$TV(t) = \text{max\_TV} \times (1 - \exp(-\text{rate\_TV} \times (1 + \text{disp\_TV})^{-\text{time\_elapsed}}))$$

But unlike BR and ER, TV is **independent of treatment assignment**. Both drug and placebo groups should have the same TV trajectory.

### Why Can't We Just Control for Time?

You might think: "Just add `week` as a covariate and call it a day."

That's close, but not quite right:

- Linear time (`week`) assumes constant improvement/decline
- Reality: symptoms might decline quickly at first, then plateau
- A Gompertz curve captures S-shaped or saturating patterns
- TV is individual-specific (Patient A's natural trajectory ≠ Patient B's)

Ignoring the curve shape means you'll misestimate "when to measure" and whether the effect is real.

---

## Putting It Together: Why Three Components?

### Example: A Patient's Actual Trajectory

Meet Sarah, a PTSD patient in the prazosin study.

**Week 1-8 (Open-label on drug)**:
- BR = 4 (prazosin is working: -4 nightmare severity)
- ER = 6 (placebo: she believes she's on the drug: -6 severity)
- TV = 1 (natural improvement from time: -1 severity)
- **Total = -11 points, she reports feeling great**

**Week 9-12 (Blinded discontinuation, she gets placebo)**:
- BR = 0 (not really on drug)
- ER = 2 (lower expectation, but still thinks might be on drug: -2 severity)
- TV = 1 (natural healing continues: -1 severity)
- **Total = -3 points, she's disappointed**

**Week 13-16 (Crossover: back on real drug)**:
- BR = 4 (prazosin returns: -4 severity)
- ER = 3 (skeptical now after the "placebo"phase: -3 severity)
- TV = 0.5 (natural healing has plateaued: -0.5 severity)
- **Total = -7.5 points**

**What a one-component model would conclude**: "Response dropped from 11 to 3 to 7.5—very noisy, hard to interpret"

**What the three-component model shows**: "BR is consistent (4, 0, 4). ER varies with belief (6, 2, 3). TV plateaus (1, 1, 0.5). The drug works; expectations modulate the effect; natural healing helps initially."

### Information Gain

By decomposing into BR, ER, TV, we can now ask:

1. **How well does the drug work?** → Look at BR only = 4 points
2. **How much of the improvement is placebo?** → Look at ER = 2-6 points (bigger than BR!)
3. **Would she improve without treatment?** → Look at TV = 0-1 points (modest natural healing)
4. **Is the drug effect stable?** → BR should be consistent across phases (yes, 4, 0, 4)
5. **Does belief amplify the drug?** → High ER when expectation high (yes)

**With a one-component model, you can't ask any of these questions.**

---

## The Math Behind It (Intuitive Version)

### The Core Model

$$Y_{ij} = BM_i + BR_j(t) + ER_j(t) + TV_j(t) + \text{noise}_{ij}$$

where:
- $i$ = participant
- $j$ = timepoint
- $BM_i$ = baseline biomarker (doesn't change)
- $BR, ER, TV$ = components that change over time

### Why Gompertz Curves?

All three components follow similar curves:

$$\text{Component}(t) = \max \times (1 - \exp(-\text{rate} \times (1 + \text{disp})^{-t}))$$

This captures real biological/psychological patterns:
- **Slow start**: Takes time for drug to accumulate or for belief to form
- **Fast middle**: Once started, rapid improvement
- **Leveling off**: Eventually you hit diminishing returns

### Why Separate Covariance?

Each component has its own variance because they operate on different timescales:

- **BR variance**: Depends on individual pharmacokinetics (metabolism varies)
- **ER variance**: Depends on suggestibility and expectation level
- **TV variance**: Depends on life circumstances and disease trajectory

$$\text{Var}(Y) = \text{Var}(BR) + \text{Var}(ER) + \text{Var}(TV) + \text{Var}(\text{noise})$$

NOT $\text{Var}(Y) = \text{Var}(BR + ER + TV + \text{noise})$

Because the components are partially independent (correlated, not additive), you need a full covariance matrix to capture the relationships.

---

## Common Questions

### Q: Why not just use a random intercept and slope for each patient?

**A**: You could, but you'd lose information. A random slope model assumes linear time effects:

$$Y_{ij} = \alpha_i + \beta_i \times t_j + \text{error}$$

Better: $$Y_{ij} = \alpha_i + BR_j + ER_j + TV_j + \text{error}$$

The three-component model is more specific about *why* the trajectory is non-linear.

### Q: Isn't TV just disease natural history? Can't you measure it separately?

**A**: In principle, yes. But:
1. **Ethical**: Can't give some patients placebo to measure TV (some need treatment)
2. **Practical**: Takes years to observe without intervention
3. **Statistical**: The model extracts TV from data more efficiently

The three-component model **infers** TV from patterns in the data.

### Q: Why is biomarker only correlated with BR, not ER?

**A**: Because the biomarker (e.g., baseline blood pressure) is a **biological marker**. It predicts how your body responds to the drug, not how much you believe the treatment will work.

If baseline blood pressure predicted placebo response, that would be suspicious—it would suggest the biomarker is confounded with psychological factors.

### Q: What if TV goes negative (symptoms worsen over time)?

**A**: That's okay! TV can be negative for patients whose condition naturally worsens. The model handles this. In fact, detecting TV < 0 is informative: "This condition naturally deteriorates; treatment effect is even larger than it appears."

---

## Why This Matters for Biomarker Validation

The whole goal is to ask: **"Does the biomarker predict treatment response?"**

### The Wrong Way

```
response = biomarker + error
```

But if response = BR + ER + TV, and you only model one thing, you're:
- Confounding biomarker × BR with biomarker × ER
- Attributing TV changes to biomarker
- Losing statistical power

### The Right Way

```
BR = β_0 + β_biomarker × biomarker + error
```

Now you're asking: "Does biomarker predict the *biological* response?" That's the precision medicine question.

The ER and TV components are nuisances—they reduce power but aren't your target. Modeling them separately means you can control for their effects and isolate the biomarker signal.

---

## Summary Table

| Component | Source | Timescale | Measured By | Why Important |
|-----------|--------|-----------|-------------|---------------|
| **BR** | Drug chemistry | Weeks | Phase contrast (on vs off drug) | Answers "Does the drug work?" |
| **ER** | Patient belief | Hours-weeks | Expectancy manipulation (open vs blinded) | Quantifies placebo effect |
| **TV** | Disease natural history | Weeks-months | Long-term trends | Controls for confounding |
| **Error** | Measurement noise | All | Residuals | Model fit |

---

## Key Takeaway

The BR-ER-TV decomposition isn't just technical bookkeeping. It's **principled measurement of different causes** of outcome changes.

When you lump everything together into a single outcome, you can't distinguish drug effect from placebo from natural healing. The decomposition lets you see the structure underneath.

And that structure is what lets you validate biomarkers as precision medicine tools.

---

## Further Reading

- **Hendrickson et al. (2020)** — Original paper introducing this framework for N-of-1 trials
- **Senn (2016)** — Variance decomposition in personalized medicine
- **Psychological mechanisms**: Placebo research shows ER is neurobiologically real, not "fake"
- **Pharmacokinetics**: Why BR follows Gompertz (saturation kinetics)

---

*Document created: March 2026*

*Audience: Undergraduate statistics students learning about mixed effects models and simulation-based power analysis*
