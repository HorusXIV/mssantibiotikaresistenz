# Micro Time-Scale Calibration

This document defines how the micro simulation maps discrete steps to real time.
The production definition is `steps_per_day: 12` because the micro layer is evaluated during a 12-hour overnight active window. Under this definition, one micro step corresponds to one real hour within that nightly window, while the macro layer still advances one full day at a time.

## Literature and Data Anchors

Use these sources as calibration anchors, not as direct YAML values. The model is
a 14-trait effective within-host model, so biological rates must be mapped to
model-scale rates.

| Area | Model use | Plausible anchor |
|---|---|---|
| Nasal carriage and colonization ecology | Carrier prevalence, host-risk stratification, persistence assumptions | S. aureus colonizes multiple body sites, with the anterior nares as a key reservoir; reviews summarize persistent/intermittent carriage and risk factors. Sources: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6186810/>, <https://www.nature.com/articles/nrmicro.2017.104> |
| Mutation and replication time | Upper/lower plausibility for `base_mutation_rate` and step length | In human nasal colonization, S. aureus generation-time estimates are about 73-81 minutes for two MRSA sequence types, with mutation-accumulation estimates around `2.0e-10` to `2.8e-10` mutations per nucleotide per generation. Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6425579/> |
| MRSA carriage duration / clearance | Target for `p_clearance` and clearance-threshold behavior | One hospital-discharge cohort estimated median MRSA clearance time at 8.5 months. Source: <https://pubmed.ncbi.nlm.nih.gov/?term=11317238> |
| Community and hospital MRSA prevalence | Initial carrier fraction, replacement carrier fraction, MRSA surveillance targets | ECDC Surveillance Atlas and CDC AR & Patient Safety Portal provide public MRSA/AMR surveillance data. Sources: <https://atlas.ecdc.europa.eu/>, <https://www.cdc.gov/healthcare-associated-infections/php/data/ar-patient-safety-portal.html> |
| Antibiotic stress / SOS response | Direction and magnitude checks for `stress_mutation_boost` | Ciprofloxacin induces a SOS-mediated S. aureus response; beta-lactams can induce SOS response and horizontal transfer of virulence factors. Sources: <https://pubmed.ncbi.nlm.nih.gov/17085555/>, <https://doi.org/10.1128/jb.188.7.2726-2729.2006> |
| Horizontal gene transfer | Plausibility for `base_hgt_rate` and transfer-gene block | Staphylococcal plasmid conjugation and mobilization contribute to antimicrobial-resistance gene spread. Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4993578/> |
| Persisters, dormancy, damage | Directional support for dormancy, repair, tolerance, and stress blocks | Persister/tolerance literature supports low-metabolic-state survival under antibiotics while cautioning that dormancy alone is insufficient. Sources: <https://pubmed.ncbi.nlm.nih.gov/27398229/>, <https://pmc.ncbi.nlm.nih.gov/articles/PMC10398667/> |

## Time-Scale Rules

Let `S_ref` be the step count at which the current parameters were calibrated and
`S_new` the target step count within the active micro window. For the project default, the active window is 12 hours and `S_new = 12`, so each micro step is one hour.

### Event probabilities

Use this for per-step Bernoulli hazards such as basal death and HGT event
probabilities:

```text
p_new = 1 - (1 - p_ref)^(opportunities_ref / opportunities_new)
```

For ordinary per-step hazards, `opportunities = steps_per_day` within the active micro window. For HGT in the
current engine, attempts happen every third step, so use `steps_per_day / 3`.

### Poisson event intensities

Use this for mutation intensity:

```text
lambda_new = lambda_ref * S_ref / S_new
```

This preserves expected events per gene per day.

### Compound fractional rates

Use this for growth-like fractional rates where repeated application compounds:

```text
r_new = exp(log(1 + r_ref) * S_ref / S_new) - 1
```

### Additive rates

Use this for damage and repair increments accumulated per step:

```text
x_new = x_ref * S_ref / S_new
```

### Step-count scales

Use inverse scaling for parameters expressed in steps:

```text
half_life_steps_new = half_life_steps_ref * S_new / S_ref
```

## Implemented Workflow

The reusable implementation is in `src/mss/simulation/micro/time_calibration.py`.
It provides:

- `rescale_micro_config_for_step_duration()`: returns a calibrated config copy.
- `describe_time_scaling()`: reports every changed parameter and factor.
- `run_micro_time_scale_ensemble()`: runs controlled within-host episodes.
- `summarize_ensemble()`: aggregates diagnostics across seeds.

The CLI is:

```bash
uv run mss-micro-time-calibrate --config config/simulation_realistic.yml \
  --target-steps-per-day 12 --active-window-hours 12 --n-days 30 --n-seeds 5
```

It writes:

- `candidate_micro_12_steps.yml`: candidate config for the 12-step overnight micro window.
- `micro_time_scaling.parquet`: parameter-by-parameter rescaling table.
- `micro_time_ensemble.parquet`: seed/day diagnostic rows.
- `micro_time_summary.parquet`: daily ensemble mean/std diagnostics.

## Calibration Framework

1. Define the real-time semantics: `steps_per_day = 12` means twelve one-hour micro updates during the overnight active window, while one macro step remains one day.
2. Convert all step-sensitive parameters with the rules above.
3. Run paired reference and target diagnostics with identical seeds.
4. Compare daily summary metrics, not single trajectories: `total_population`,
   `resistant_fraction`, `p_clearance`, `n_strains`, transmissibility, and
   severity.
5. Tune only residual deviations after rescaling. Prefer grouped calibration:
   population/growth first, clearance second, mutation/HGT third, dormancy and
   damage last.
6. Freeze the active-window candidate into a scenario YAML only after the diagnostics
   are stable across multiple seeds and at least one antibiotic-stress scenario.

## Validation Strategy

Use three validation layers:

1. **Time-window tests:** the 12-step production run should report `step_duration_hours = 1.0` for `active_window_hours = 12`; if another resolution is tested, paired reference and rescaled target runs should show similar ensemble means over 30-90 days.
2. **Empirical plausibility:** clearance should remain compatible with MRSA
   carriage-duration targets; mutation/HGT outputs should not imply implausibly
   rapid resistant-dominant sweeps without antibiotic pressure.
3. **Coupled-model checks:** run the full macro/micro simulation and compare
   prevalence, acquisition, mean resistance fraction, and LOS/mortality outputs
   against the existing calibrated scenario and surveillance targets.

Treat online literature as constraints on observable outcomes and orders of
magnitude. Parameters such as `base_mutation_rate`, `mutation_std`,
`base_hgt_rate`, damage scales, and dormancy synergies are effective model
parameters and remain identifiable only through simulation-output calibration.
