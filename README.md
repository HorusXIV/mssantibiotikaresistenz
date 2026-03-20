# MSS - Antibiotic Resistance Simulation

A multi-scale simulation framework for modeling antibiotic resistance dynamics in hospital environments. The system operates on two interconnected levels: a macro-level graph-based hospital network simulation and a micro-level evolutionary algorithm simulating bacterial strain dynamics within individual patients.

## System Architecture

```
+------------------+       +------------------+       +------------------+
|                  |       |                  |       |                  |
|  Macro Layer     | <---> |  Patient         | <---> |  Micro Layer     |
|  (Hospital Net)  |       |  (Interface)     |       |  (Evolution Sim) |
|                  |       |                  |       |                  |
+------------------+       +------------------+       +------------------+
     |                           |                           |
     | - Hospital graph          | - Base stats              | - Genetic algorithm
     | - Patient transfers       | - Modifiers               | - Strain evolution
     | - Admissions/discharges   | - State tracking          | - Resistance dynamics
     | - Department management   | - Data exchange           | - Fitness landscapes
     +---------------------------+---------------------------+
```

### Design Philosophy: Base Stats and Modifiers

The simulation employs a base stats and modifiers approach for tractability:

- **Base stats** represent intrinsic, stable patient characteristics (age, immune strength, vulnerability)
- **Modifiers** are dynamic multipliers applied by environmental factors (hospital hygiene, antibiotic regimens) or strain-driven effects (lethality, transmissibility)

This separation allows clear attribution of effects and simplifies parameter tuning.

## Simulation Levels

### Macro Level (Hospital Network)

The macro simulation models a graph-based hospital network with daily time steps. It manages:

| Component | Description |
|-----------|-------------|
| **Hospitals** | Nodes with diagnostic speed, patient capacity |
| **Departments** | WARD and ICU with distinct hygiene levels |
| **Patient Flow** | Admissions, discharges, inter-hospital transfers |
| **Isolation** | Per-department isolation effectiveness |
| **Antibiotic Policy** | Drug class selection, dosage levels |
| **Hygiene Standards** | Per-department hygiene levels (0-1 scale) |

The macro layer provides daily context to each patient via `PatientDailyContext`.

### Micro Level (Evolutionary Simulation)

The micro simulation runs 12 steps per day for carrier patients, using evolutionary algorithms to model:

| Aspect | Description |
|--------|-------------|
| **Bacterial Strains** | Multiple genotypes with distinct resistance profiles |
| **Selection Pressure** | Antibiotic exposure drives resistance evolution |
| **Fitness Landscapes** | Trade-offs between resistance and growth rate |
| **Within-Host Dynamics** | Resistant fraction, dominant genotype tracking |

Strain-driven outputs returned to macro:
- `lethality_modifier`: Increases patient baseline death risk
- `severity_modifier`: Affects disease severity
- `relative_transmissibility`: Modifies transmission probability
- `p_clearance`: Daily probability of carriage clearance

## Macro-Simulation Implementation

The macro-simulation module (`macro_simulation/`) implements a hospital-network carrier model with daily transmission and clearance events.

### SimulationConfig

All tunable parameters are collected in `SimulationConfig` (`macro_simulation/simulation.py`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `base_hygiene` | 0.7 | Hospital-wide hygiene level (0–1) |
| `base_isolation_effectiveness` | 0.8 | Fraction by which isolation reduces transmission |
| `base_diagnostic_speed` | 0.5 | Carrier detection speed (passed to PatientDailyContext) |
| `base_transmission_rate` | 0.05 | Per-carrier-per-susceptible daily beta |
| `daily_contact_attempts` | 80.0 | Effective daily contact opportunities used in the macro hazard model |
| `icu_abx_probability` | 0.60 | Daily P(ABX=on) in ICU |
| `ward_abx_probability` | 0.15 | Daily P(ABX=on) in WARD |
| `carrier_isolation_probability` | 0.30 | Daily detection probability for non-isolated carriers |

### Transmission Model

Colonisation probability per susceptible patient per day:

```
p_colonize = min(1, carrier_force × susceptibility_multiplier_for_macro())
```

**Carrier Force / Hazard:**
- Summed over all carriers in the same hospital node
- Each carrier contributes `base_transmission_rate × transmission_multiplier_for_macro()`
- `transmission_multiplier_for_macro()` = `sociability × relative_transmissibility`
- Isolated carriers are reduced by `(1 − isolation_effectiveness)`
- Total infectiousness is normalized by hospital occupancy
- Colonisation uses a hazard model with `daily_contact_attempts`, which prevents immediate saturation in large hospitals

**Susceptibility:**
- `susceptibility_multiplier_for_macro()` = `vulnerability / immune_strength`
- Prior-infection history flag halves susceptibility

### Antibiotic Profiles

| Department | P(ABX=on) | Classes available | Dose levels |
|------------|-----------|-------------------|-------------|
| ICU | 60 % | beta_lactam, fluoroquinolone, glycopeptide, macrolide, aminoglycoside | low / std / high |
| WARD | 15 % | beta_lactam, fluoroquinolone, glycopeptide, macrolide, aminoglycoside | low / std / high |

### Daily Simulation Loop

Each call to `step()` processes all hospitals in sequence:

```
for each hospital:
    1. CLEARANCE:    patient.should_clear_today(rng) → C→S via patient.clear_carriage()
    2. CONTEXT:      patient.update_context(ctx) for every patient
    3. MICRO (opt):  for each carrier, make_micro_request() -> micro.process_batch() -> apply_micro_response()
                     (one micro day per macro day; default micro day = 12 internal steps)
    4. TRANSMISSION: force-of-infection model → S→C for susceptible patients
```

**Deterministic RNG:**
- Single seeded `random.Random` instance for all stochastic decisions
- Context building consumes exactly 4 draws per patient per step

### State Classification

| State | Value | Description |
|-------|-------|-------------|
| SUSCEPTIBLE | `"S"` | No active carriage; `episode_id=None`, `resistant_fraction=0.0` |
| CARRIER | `"C"` | Active carriage; micro-derived fields populated |

### Context Provided to Patient

After each step, every patient holds an updated `PatientDailyContext`:

```python
PatientDailyContext(
    hospital_id="hospital_001",
    department=Department.WARD,          # or Department.ICU
    hygiene_level=0.7,                   # base_hygiene
    isolation_effectiveness=0.8,         # base_isolation_effectiveness
    diagnostic_speed=0.5,               # base_diagnostic_speed
    is_isolated=False,                   # updated daily by detection model
    regimen=AntibioticRegimen(on=True, abx_class="beta_lactam", dose_level="std")
)
```

### Usage

```python
from macro_simulation.simulator import MacroSimulator
from macro_simulation.simulation import SimulationConfig

config = SimulationConfig(base_hygiene=0.8, icu_abx_probability=0.7)
sim = MacroSimulator(config=config, n_hospitals=10, seed=42)

sim.admit(patient, hospital_id="hospital_001", department=Department.ICU)
sim.transfer(patient, to_hospital_id="hospital_002")
sim.discharge(patient)

# Advance one day
sim.step()
```

### Coupled Runner (Macro + Micro)

Use the dedicated project runner to execute a full coupled simulation:

```bash
uv run python run_coupled_simulation.py
```

All tunable runner, population, macro, and micro parameters are read from:

- [shared/config.yml](/home/lukas/src/Sem4/MSS/shared/config.yml)

That file includes:

- `run`: simulation duration, seed, run id, quiet mode
- `population`: hospital count, initial department, susceptible/carrier counts, patient templates
- `macro`: hospital/network transmission parameters
- `micro`: within-host evolution parameters and worker count

Quick smoke run after editing the YAML:

```bash
uv run python run_coupled_simulation.py
```

## Micro-Simulation Implementation

The micro-simulation module (`micro_simulation/`) implements a within-host evolutionary algorithm that models bacterial population dynamics under antibiotic pressure.

### Bacterial Genome

Each bacterial strain is represented by a 14-gene genome with normalized float values (0.0-1.0):

| Gene | Index | Category | Function |
|------|-------|----------|----------|
| `GROWTH_BASE` | 0 | Metabolism | Baseline replication rate |
| `METABOLIC_OPTIMIZATION` | 1 | Metabolism | Compensates resistance fitness costs |
| `EFFLUX_PUMPS` | 2 | Resistance | Active drug efflux mechanism |
| `TARGET_MODIFICATION` | 3 | Resistance | Altered antibiotic target sites |
| `PERMEABILITY_REDUCTION` | 4 | Resistance | Reduced membrane permeability |
| `VIRULENCE` | 5 | Virulence | Disease severity, lethality |
| `STEALTH` | 6 | Survival | Immune evasion capability |
| `ADHESION` | 7 | Survival | Host colonization, transmissibility |
| `MUTATION_RATE_MODIFIER` | 8 | Evolvability | Intrinsic mutation rate |
| `HGT_COMPETENCE` | 9 | Evolvability | Horizontal gene transfer receptivity |
| `DNA_REPAIR` | 10 | Lifecycle | Repairs lineage damage, slows senescence |
| `DORMANCY_PROPENSITY` | 11 | Lifecycle | Enables low-growth persistence under stress |
| `STRESS_RESPONSE` | 12 | Lifecycle | Activates protective stress programs |
| `DAMAGE_TOLERANCE` | 13 | Lifecycle | Buffers mortality from accumulated damage |

### Fitness Calculation

Fitness is computed dynamically based on the current environment:

```
Fitness = (growth_base - net_resistance_costs) * ABX_survival * Immune_survival
```

**Epistasis and Compensation:**
- Each resistance gene carries a base fitness cost (efflux: 0.15, target modification: 0.12, permeability: 0.08)
- `metabolic_optimization` reduces total resistance costs by up to 80% (compensatory mutations)
- Net cost = raw_cost * (1.0 - optimization * 0.8)

**Antibiotic Survival:**
- Each antibiotic class has a profile defining efficacy of each resistance mechanism
- Protection = weighted sum of (resistance_gene * mechanism_efficacy)
- Survival = 1 - effective_kill_rate * (1 - protection)

**Immune Survival:**
- Base clearance rate modulated by immune_strength and immune_status
- `stealth` gene provides up to 70% immune evasion

### Antibiotic Profiles

| Class | Efflux Efficacy | Target Mod Efficacy | Permeability Efficacy | Base Kill Rate |
|-------|-----------------|---------------------|----------------------|----------------|
| beta_lactam | 0.3 | 0.8 | 0.4 | 0.75 |
| fluoroquinolone | 0.6 | 0.7 | 0.3 | 0.80 |
| aminoglycoside | 0.4 | 0.5 | 0.6 | 0.70 |
| macrolide | 0.7 | 0.4 | 0.3 | 0.65 |
| tetracycline | 0.8 | 0.3 | 0.2 | 0.60 |
| glycopeptide | 0.2 | 0.9 | 0.5 | 0.85 |

### Population Model

The simulation tracks bacterial populations as **strains** rather than individual agents for computational efficiency:

- `StrainPopulation`: Contains genome array (n_strains x 10) and population counts
- Maximum 50 strains tracked simultaneously
- Strains below threshold population are pruned
- Carrying capacity: 10^9 bacteria

### Daily Simulation Loop (12 Steps)

Each day is divided into 12 discrete time steps:

```
for step in 1..12:
    1. SELECTION: Grow/shrink populations based on relative fitness
    2. MUTATION:  Apply Gaussian noise to genes, create new strains
    3. HGT:       Transfer resistance genes between strains (every 3rd step)
    4. CONSOLIDATE: Remove extinct strains, enforce max strain limit
```

**Mutation Dynamics:**
- Base mutation rate: 0.01 per gene per step
- Stress-induced boost: 3x under antibiotic pressure
- Genome's `mutation_rate_modifier` further modulates rate
- Mutations apply Gaussian noise (std=0.05) to selected genes

**Horizontal Gene Transfer:**
- Probability based on `hgt_competence` gene
- Primarily transfers resistance genes (efflux, target modification, permeability, metabolic optimization)
- Creates recombinant strains with blended gene values

### Genotype Classification

Strains are classified based on average resistance score:

| Genotype | Resistance Score | Description |
|----------|------------------|-------------|
| S | < 0.2 | Susceptible |
| R1 | 0.2 - 0.4 | Low resistance |
| R2 | 0.4 - 0.7 | Medium resistance |
| R3 | >= 0.7 | High resistance |

### Output to Macro Layer

After 12 steps, the micro-simulation returns:

```python
{
    "updated_state": {
        "resistant_fraction": 0.15,      # Population-weighted resistance
        "dominant_genotype": "R1"        # Most populous strain's classification
    },
    "derived_effects": {
        "relative_transmissibility": 1.2,  # Based on adhesion + virulence
        "lethality_modifier": 0.95,        # Based on virulence
        "severity_modifier": 1.1,          # Based on virulence + adhesion
        "p_clearance": 0.02                # Based on population size + stealth
    }
}
```

### Usage

```python
from micro_simulation import MicroSimulator, SimulationConfig

# Configure simulation
config = SimulationConfig(
    steps_per_day=12,
    base_mutation_rate=0.01,
    stress_mutation_boost=3.0
)

# Create simulator (maintains episode state)
simulator = MicroSimulator(config=config)

# Process patient request
request = patient.make_micro_request(run_id="sim1", day=5, dt_days=1, seed=42)
response = simulator.process_request(request)

# Apply to patient
patient.apply_micro_response(response)
```

### Batch Processing

For HPC scenarios with thousands of patients:

```python
# Process multiple patients in parallel
requests = [p.make_micro_request(...) for p in carrier_patients]
responses = simulator.process_batch(requests, parallel=True)

for patient, response in zip(carrier_patients, responses):
    patient.apply_micro_response(response)
```

## Patient Interface

The `Patient` class (in `exchange/patient.py`) serves as the data exchange interface between macro and micro layers.

### Patient Base Stats

| Attribute | Type | Description |
|-----------|------|-------------|
| `age_years` | int | Patient age, affects multiple modifiers |
| `compliance` | float (0-1) | Adherence to treatment, isolation protocols |
| `vulnerability` | float | Susceptibility multiplier |
| `immune_strength` | float | Immune system effectiveness multiplier |
| `immune_status` | str | "normal" or "suppressed" |
| `sociability` | float | Contact rate multiplier for transmission |

### Patient State

| Attribute | Type | Description |
|-----------|------|-------------|
| `state` | HealthState | SUSCEPTIBLE (S) or CARRIER (C) |
| `treatment_phase` | TreatmentPhase | NONE, ON_ABX, POST_ABX |
| `history_flags` | Set[str] | Medical history markers |
| `is_isolated` | bool | Current isolation status |

### Micro-Driven Modifiers

| Attribute | Type | Source | Effect |
|-----------|------|--------|--------|
| `resistant_fraction` | float (0-1) | Micro | Proportion of resistant bacteria |
| `dominant_genotype` | str | Micro | Current dominant strain identifier |
| `relative_transmissibility` | float | Micro | Transmission rate multiplier |
| `lethality_modifier` | float | Micro | Death risk multiplier |
| `severity_modifier` | float | Micro | Disease severity multiplier |
| `p_clearance` | float | Micro | Daily clearance probability |

## Data Flow

### Daily Simulation Cycle

1. **Macro provides context** via `patient.update_context(ctx)`:
   - Hospital/department assignment
   - Hygiene and isolation parameters
   - Antibiotic regimen

2. **Patient generates micro request** via `patient.make_micro_request()`:
   - Host factors (age, immune status, vulnerability)
   - Current antibiotic exposure
   - Initial bacterial state

3. **Micro returns response** applied via `patient.apply_micro_response()`:
   - Updated resistant fraction and dominant genotype
   - Derived modifiers (transmissibility, lethality, clearance)

4. **Macro uses patient methods** for transmission decisions:
   - `transmission_multiplier_for_macro()`: Combined sociability and strain transmissibility
   - `susceptibility_multiplier_for_macro()`: Vulnerability adjusted by immune status
   - `daily_death_risk_multiplier()`: Vulnerability times lethality modifier
   - `should_clear_today(rng)`: Stochastic clearance check

## Project Structure

```
MSS/
├── exchange/
│   └── patient.py              # Patient interface class
├── macro_simulation/
│   ├── __init__.py
│   ├── simulation.py           # SimulationConfig dataclass
│   └── simulator.py            # MacroSimulator (admit/discharge/transfer/step)
├── micro_simulation/
│   ├── __init__.py             # Module exports
│   ├── genome.py               # BacterialGenome, fitness functions
│   ├── simulation.py           # StrainPopulation, simulate_day()
│   └── simulator.py            # MicroSimulator batch interface
├── tests/
│   └── integration_tests/
│       ├── test_macro_patient_integration.py
│       └── test_micro_patient_integration.py
├── System_Overview/
│   ├── build_gephi_graphs.py   # Network visualization export
│   ├── amr_system_map.*        # System-level network files
│   └── amr_transfer_network.*  # Transfer network files
└── .gitlab-ci.yml              # CI: lint (ruff, black) + integration tests
```

## Requirements

- Python 3.10+
- Dependencies listed in `pyproject.toml`

## References

This simulation framework is designed for research into antibiotic resistance spread dynamics. The model abstracts complex biological and epidemiological processes while maintaining sufficient fidelity for policy analysis and intervention testing.
