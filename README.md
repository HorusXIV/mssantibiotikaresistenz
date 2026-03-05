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
│   └── patient.py          # Patient interface class
├── macro_simulation/       # Hospital network simulation
├── micro_simulation/       # Evolutionary bacterial simulation
├── build_gephi_graphs.py   # Network visualization export
├── amr_system_map.*        # System-level network files
└── amr_transfer_network.*  # Transfer network files
```

## Requirements

- Python 3.10+
- Dependencies listed in `pyproject.toml`

## References

This simulation framework is designed for research into antibiotic resistance spread dynamics. The model abstracts complex biological and epidemiological processes while maintaining sufficient fidelity for policy analysis and intervention testing.
