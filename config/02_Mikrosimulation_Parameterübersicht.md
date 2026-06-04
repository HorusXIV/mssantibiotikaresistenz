# Mikrosimulation: Parameterübersicht

Diese Übersicht beschreibt die Mikrosimulation so, wie sie aktuell im Projekt ausgeführt wird. Referenz ist die produktive Konfiguration `config/simulation_realistic.yml`; die Code-Defaults in `src/mss/simulation/micro/engine.py` weichen an mehreren Stellen davon ab und sind deshalb nicht die primäre Quelle für die Projektläufe.

Jeder Parameter ist einem Typ zugeordnet:

- **geschätzt** - Literatur oder öffentliche Surveillance-Daten liefern eine brauchbare Zielgröße oder Größenordnung.
- **kalibrieren** - der Wert ist ein effektiver Modellparameter; Literatur kann höchstens eine Richtung oder Größenordnung liefern.
- **technisch** - numerischer/architektonischer Parameter, der Rechenkosten, Logging oder Stabilität steuert.
- **Modellannahme** - bewusst gesetzte Strukturannahme, weil der Code ein vereinfachtes Trait-Modell statt ein biophysikalisches Detailmodell verwendet.

---

## Was die Mikrosimulation modelliert

Die Mikrosimulation ist ein zustandsbehaftetes Within-Host-Modell für Carrier-Patienten. Sie simuliert nicht einzelne Bakterienzellen und keine echten DNA-Sequenzen, sondern eine Population aus Stämmen. Jeder Stamm hat:

- einen kontinuierlichen Genomvektor mit `NUM_GENES = 14` Trait-Werten in `[0, 1]`
- eine absolute Populationsgröße
- ein Linienalter (`lineage_ages`)
- eine akkumulierte Schadenslast (`damage_loads`)
- Strain-, Parent-, Donor- und Founder-IDs für Rohlogs

Ein Makro-Tag wird in `steps_per_day` Mikro-Schritte zerlegt. In jedem Schritt passiert:

1. **Selektion und Populationsdynamik**: Wachstum und Tod werden aus Fitness, Antibiotikadruck, Immundruck, Dormanz, Schaden und Alterung berechnet.
2. **Mutation**: Genwerte einzelner Stämme werden verrauscht; daraus entstehen neue Teilstämme.
3. **Horizontaler Gentransfer (HGT)**: alle drei Schritte können Resistenz-/Persistenzgene zwischen Stämmen gemischt werden.
4. **Konsolidierung**: sehr kleine Stämme werden entfernt; bei zu vielen Stämmen bleiben die größten erhalten.

Der Episodenzustand wird über Tage gespeichert. Eine neue Tagesberechnung beginnt also mit der Population des Vortages, nicht mit einem frischen Seed.

---

## Kopplung mit der Makrosimulation

### Makro -> Patient -> Mikro

`MacroSimulator._build_context()` erzeugt pro Tag den Patientenkontext:

- Station (`ward`, `icu`, `isolation`)
- Isolation und Hygiene
- Antibiotika-Regime (`on`, `abx_class`, `dose_level`)

`Patient.update_context()` übernimmt diesen Kontext und leitet `adherence` aus `compliance` ab. ICU-Patienten erhalten dabei einen Bonus von `+0.1`, gedeckelt auf `1.0`.

`Patient.make_micro_request()` baut nur für Carrier mit `episode_id` einen Mikro-Request. Aktiv genutzt werden in der Engine:

| Request-Feld | Wirkung in Mikro |
|---|---|
| `episode_id` | Schlüssel für persistenten Episodenzustand |
| `patient_id`, `t_day`, `seed` | Identität, Logging, Reproduzierbarkeit |
| `abx.on`, `abx.class`, `abx.dose_level` | Antibiotika-Selektion und Stress |
| `adherence` | skaliert effektive Antibiotika-Kill-Rate |
| `host.immune_strength` | skaliert Immunclearance und spätere `p_clearance` |
| `initial_state.resistant_fraction` | resistenter Anteil bei neuer Episode |
| `initial_state.dominant_genotype` | Seed-Klasse für resistente Startstämme |
| `initial_state.dominant_strain_name` | bevorzugter Name für den dominanten Seed |
| `initial_state.seed_genome` | bei Transmission geerbtes dominantes Genom |

Diese Felder werden transportiert, aber aktuell nicht direkt von der Mikro-Engine verwendet:

- `setting`
- `host.age_years`
- `host.history_flags`
- Hygiene, Isolation und Diagnostik wirken derzeit nur auf Makro-Ebene.

### Mikro -> Patient -> Makro

`population_to_response()` liefert:

| Response-Feld | Wirkung |
|---|---|
| `updated_state.resistant_fraction` | Anteil resistenter Bakterien in der aktuellen Within-Host-Population |
| `updated_state.dominant_genotype` | dominante Resistenzklasse `S`, `R1`, `R2`, `R3` |
| `updated_state.dominant_strain_name` | Name des größten Stamms |
| `updated_state.dominant_genome` | Genom des dominanten Stamms; wird bei Transmission vererbt |
| `derived_effects.relative_transmissibility` | skaliert Makro-Transmission über `transmission_multiplier_for_macro()` |
| `derived_effects.p_clearance` | tägliche Wahrscheinlichkeit für `CARRIER -> SUSCEPTIBLE` |
| `derived_effects.severity_modifier` | skaliert Makro-Aufenthaltsverlängerung und Mortalität |
| `population_stats.total_population`, `n_strains` | Logging und Plausibilitätsprüfung |

Wichtig: `compute_lethality()` existiert in `genome.py`, aber `population_to_response()` gibt aktuell keinen `lethality_modifier` zurück und `Patient` speichert kein solches Feld. Die Mortalitätskopplung läuft im aktuellen Code über `severity_modifier`.

---

## Genommodell

Das Genom ist ein Trait-Vektor, keine Sequenz. Ein "Gen" bedeutet hier: kontinuierliche Ausprägung einer Eigenschaft.

| Index | Gen | Hauptwirkung |
|---:|---|---|
| 0 | `GROWTH_BASE` | Basiswachstum |
| 1 | `METABOLIC_OPTIMIZATION` | reduziert Fitnesskosten von Resistenz |
| 2 | `EFFLUX_PUMPS` | Antibiotika-Efflux |
| 3 | `TARGET_MODIFICATION` | Schutz durch verändertes Antibiotika-Ziel |
| 4 | `PERMEABILITY_REDUCTION` | geringeres Eindringen von Antibiotika |
| 5 | `VIRULENCE` | erhöht Severity und Lethality |
| 6 | `STEALTH` | Immunevasion, senkt Clearance |
| 7 | `ADHESION` | erhöht Transmission und Severity |
| 8 | `MUTATION_RATE_MODIFIER` | stammspezifische Mutationsrate |
| 9 | `HGT_COMPETENCE` | stammspezifische HGT-Wahrscheinlichkeit |
| 10 | `DNA_REPAIR` | Reparaturkapazität |
| 11 | `DORMANCY_PROPENSITY` | Dormanz, Persistenz, Wachstumskosten |
| 12 | `STRESS_RESPONSE` | Stressantwort |
| 13 | `DAMAGE_TOLERANCE` | Toleranz gegen Schadenslast |

Die Genotyp-Klasse ist abgeleitet aus dem Mittelwert von `EFFLUX_PUMPS`, `TARGET_MODIFICATION` und `PERMEABILITY_REDUCTION`:

| Resistenzscore | Klasse |
|---:|---|
| `< 0.2` | `S` |
| `< 0.4` | `R1` |
| `< 0.7` | `R2` |
| `>= 0.7` | `R3` |

`resistant_fraction` ist nicht der dominante Genotyp, sondern der Populationsanteil mit Resistenzscore `>= 0.3`.

---

## Zentrale Formeln

### Fitness

```text
fitness = (growth_base - net_resistance_costs)
          * abx_survival
          * immune_survival
```

Resistenzkosten:

```text
raw_costs = efflux * 0.15 + target_mod * 0.12 + permeability * 0.08
net_costs = raw_costs * (1 - 0.8 * metabolic_optimization)
```

Antibiotika-Überleben:

```text
effective_kill = base_kill_rate * dose_multiplier * adherence
protection = weighted_resistance / 1.5, capped at 0.95
abx_survival = 1 - effective_kill * (1 - protection)
```

Immun-Überleben:

```text
base_clearance = 0.15 * immune_strength
evasion = stealth * 0.7
immune_survival = 1 - base_clearance * (1 - evasion)
```

### Selektion und Populationsdynamik

```text
relative_fitness = fitness / mean_fitness
selection_factor = relative_fitness ** selection_strength
growth = growth_rate_per_step * selection_factor * fitness
```

Dormanz reduziert Wachstum, kann aber Replikationsdruck, Schaden und Turnover senken. Schaden entsteht aus Basisverschleiß, Replikationsdruck und Umweltstress; Reparatur hängt von Reparatur-, Dormanz- und Stressgenen ab.

```text
net_growth = growth - death
population_next = population * (1 + net_growth)
```

Danach wird demografische Stochastik angewandt. Wenn die Gesamtpopulation über `carrying_capacity` liegt, werden alle Stämme proportional herunterskaliert.

### Mutation

```text
stress_multiplier = 1 + abx_stress * (stress_mutation_boost - 1)
effective_rate = base_mutation_rate * stress_multiplier
strain_rate = effective_rate * (0.5 + mutation_rate_modifier)
n_mutations ~ Poisson(strain_rate * NUM_GENES)
```

Bei Mutation entsteht ein neuer Stamm. Pro mutiertem Gen wird ein normalverteilter Wert mit Standardabweichung `mutation_std` addiert und auf `[0, 1]` begrenzt.

### Horizontaler Gentransfer

HGT läuft alle drei Mikro-Schritte.

```text
hgt_prob = base_hgt_rate * (0.5 + HGT_COMPETENCE)
```

Transferierbar sind:

- `EFFLUX_PUMPS`
- `TARGET_MODIFICATION`
- `PERMEABILITY_REDUCTION`
- `METABOLIC_OPTIMIZATION`
- `DORMANCY_PROPENSITY`
- `STRESS_RESPONSE`
- `DAMAGE_TOLERANCE`

Pro HGT-Ereignis wird ein Donor populationsgewichtet gewählt. Für jedes transferierbare Gen entscheidet `hgt_gene_transfer_prob`, ob es gemischt wird. Der Donorwert ersetzt den Empfängerwert nicht hart; er wird mit einem Blend-Faktor zwischen `0.3` und `0.7` eingemischt.

### Clearance

```text
if total_population <= 0: p_clearance = 1.0
elif total_population < clearance_threshold: p_clearance = 0.95
elif total_population < min_population: base_prob = 0.3
else: base_prob = 0.02 * (1 - total_population / carrying_capacity)

immune_mult = immune_strength
stealth_effect = 1 - avg_stealth * 0.5
p_clearance = clip(base_prob * immune_mult * stealth_effect, 0.001, 0.95)
```

---

## Parameterübersicht

### Ausführung und technische Struktur

| Parameter | Typ | Wert | Wirkung | Einschätzung |
|---|---|---:|---|---|
| `micro.workers` | technisch | `null` | nutzt CPU-Anzahl; beeinflusst nur Performance | nicht biologisch kalibrieren |
| `steps_per_day` | technisch / kalibrieren | `12` | Anzahl Mikro-Schritte pro Makro-Tag; mehr Schritte bedeuten mehr Selektions-, Mutations- und HGT-Gelegenheiten | primär numerische Zeitskala, aber mit Dynamik gekoppelt |
| `max_strains` | technisch | `40` | maximale Anzahl aktiver Stämme pro Episode | Stabilitäts-/Speicherparameter; Sensitivität prüfen |
| `strain_prune_threshold` | technisch / kalibrieren | `200` | entfernt Stämme unterhalb absoluter Populationsgröße | beeinflusst Erhalt seltener Resistenzvarianten; Sensitivität nötig |
| `founder_pool_size` | technisch / Modellannahme | `32` | Größe der globalen Founder-Bibliothek | Diversität der Startpopulation; nicht direkt messbar |
| `founder_pool_seed` | technisch | `1` | deterministische Founder-Erzeugung | Reproduzierbarkeit |
| `founder_pool_gene_noise_std` | Modellannahme | `0.02` | Rauschen um archetypische Founder-Genome | Trait-Modell-Skala, nicht direkt biologisch |
| `gene_presence_threshold` | technisch | `0.2` | Schwelle für Gen-Präsenz in Rohlogs | Logging-Schwelle, nicht Engine-Dynamik |

### Startpopulation neuer Episoden

Diese Werte stehen nicht in YAML, sondern in `StrainPopulation.create_initial()`.

| Parameter | Typ | Wert | Wirkung | Einschätzung |
|---|---|---:|---|---|
| `initial_population` | geschätzt / kalibrieren | `1e6` | Startgröße einer neuen Within-Host-Population | S. aureus-Lasten sind aus CFU-Studien ableitbar, aber Mapping von CFU/Swab auf Modellpopulation muss kalibriert werden |
| `n_susceptible_strains` | Modellannahme | `3` | Anzahl sensibler Seed-Stämme | modelliert Startdiversität |
| `n_resistant_strains` | Modellannahme | `2` | Anzahl resistenter Seed-Stämme bei `resistant_fraction > 0` | modelliert Startdiversität |
| `resistant_fraction` | geschätzt | aus Patient / Template | Anteil resistenter Startpopulation | aus Surveillance/Screening ableitbar, z. B. MRSA-Anteil |
| `dominant_genotype` | geschätzt / Modellannahme | aus Patient / Template | Resistenzklasse des dominanten Seeds | epidemiologisch grob ableitbar; Klassen `R1-R3` sind Modellabstraktion |
| `seed_genome` | Modellannahme | optional | übernimmt dominantes Quellgenom bei Transmission | wichtig für Strain-Kontinuität, nicht externer Parameter |

### Antibiotika-Profile

Diese Werte sind Code-Konstanten in `genome.py`, nicht YAML-Parameter.

| ABX-Klasse | Efflux | Target Mod | Permeability | Base Kill | Einschätzung |
|---|---:|---:|---:|---:|---|
| `none` | 0.0 | 0.0 | 0.0 | 0.0 | technisch |
| `beta_lactam` | 0.3 | 0.8 | 0.4 | 0.75 | online qualitativ ableitbar, quantitativ kalibrieren |
| `fluoroquinolone` | 0.6 | 0.7 | 0.3 | 0.80 | online qualitativ ableitbar, quantitativ kalibrieren |
| `aminoglycoside` | 0.4 | 0.5 | 0.6 | 0.70 | online qualitativ ableitbar, quantitativ kalibrieren |
| `macrolide` | 0.7 | 0.4 | 0.3 | 0.65 | online qualitativ ableitbar, quantitativ kalibrieren |
| `tetracycline` | 0.8 | 0.3 | 0.2 | 0.60 | online qualitativ ableitbar, quantitativ kalibrieren |
| `glycopeptide` | 0.2 | 0.9 | 0.5 | 0.85 | online qualitativ ableitbar, quantitativ kalibrieren |

Dosis-Multiplikatoren:

| Parameter | Typ | Wert | Wirkung |
|---|---|---:|---|
| `DOSE_MULTIPLIERS.low` | Modellannahme | `0.6` | schwächerer Kill |
| `DOSE_MULTIPLIERS.std` | Modellannahme | `1.0` | Referenz |
| `DOSE_MULTIPLIERS.high` | Modellannahme | `1.4` | stärkerer Kill |

Antibiotikamechanismen sind gut identifizierbar, aber die Zahlen in diesem Trait-Modell sind keine pharmakokinetischen Parameter. Sie müssen gegen beobachtete Resistenz- und Clearance-Dynamik kalibriert werden.

### Wachstum, Population und Clearance

| Parameter | Typ | Wert | Wirkung | Einschätzung |
|---|---|---:|---|---|
| `carrying_capacity` | geschätzt / kalibrieren | `5e8` | Obergrenze der Within-Host-Gesamtpopulation | CFU-Studien geben Größenordnungen; Modellskala kalibrieren |
| `min_population` | kalibrieren | `100` | unterhalb davon springt `base_prob` für Clearance auf `0.3` | effektive Schwelle, nicht direkt messbar |
| `clearance_threshold` | kalibrieren | `1000` | unterhalb davon ist `p_clearance = 0.95` | effektive Schwelle, mit Carrier-Dauer kalibrieren |
| `growth_rate_per_step` | geschätzt / kalibrieren | `0.18` | Wachstum pro Schritt vor Fitness- und Dormanzkorrektur | S. aureus-Generationszeiten geben Obergrenzen; per-step-Wert kalibrieren |
| `death_rate_per_step` | kalibrieren | `0.06` | basale Todesrate pro Schritt | effektiver Modellparameter |
| `selection_strength` | kalibrieren | `2.5` | Exponent auf relative Fitness | steuert, wie schnell Gewinner-Stämme sweepen |

### Mutation und HGT

| Parameter | Typ | Wert | Wirkung | Einschätzung |
|---|---|---:|---|---|
| `base_mutation_rate` | geschätzt / kalibrieren | `0.012` | Basismutation pro Gen und Schritt im Trait-Modell | echte S. aureus-Mutationsraten sind ableitbar; Mapping auf Trait-Schritte muss kalibriert werden |
| `mutation_std` | kalibrieren | `0.025` | Größe einer Trait-Veränderung | reine Modellskala |
| `stress_mutation_boost` | geschätzt / kalibrieren | `40.0` | Mutationsboost unter Antibiotikastress | SOS-/Stressmutagenese ist belegt; Faktor im Modell kalibrieren |
| `base_hgt_rate` | geschätzt / kalibrieren | `0.03` | Basiswahrscheinlichkeit eines HGT-Ereignisses pro eligiblem Stamm und Schritt | HGT bei S. aureus ist belegt; effektive Within-Host-Rate kalibrieren |
| `hgt_gene_transfer_prob` | kalibrieren | `0.25` | Wahrscheinlichkeit pro transferierbarem Gen innerhalb eines HGT-Ereignisses | Trait-Mischparameter, nicht direkt messbar |

### Schaden, Lebenszyklus und Dormanz

| Parameter | Typ | Wert | Wirkung | Einschätzung |
|---|---|---:|---|---|
| `base_damage_per_step` | kalibrieren | `0.004` | Hintergrundschaden pro Schritt | latent, nicht direkt beobachtbar |
| `replication_damage_factor` | kalibrieren | `0.03` | zusätzlicher Schaden durch Replikationsdruck | effektiver Modellparameter |
| `stress_damage_factor` | kalibrieren | `0.06` | zusätzlicher Schaden durch Umweltstress | effektiver Modellparameter |
| `repair_rate_per_step` | kalibrieren | `0.08` | Abbau von Schadenslast | effektiver Modellparameter |
| `age_mortality_scale` | kalibrieren | `0.001` | Einfluss des Linienalters auf Turnover | latent |
| `damage_mortality_scale` | kalibrieren | `0.025` | Einfluss der Schadenslast auf Turnover | latent |
| `lifecycle_half_life_steps` | kalibrieren | `200` | Zeitskala für Alterungsdruck | latent |
| `max_damage_load` | technisch / kalibrieren | `5.0` | Sättigung der Schadenslast | numerische Skala |
| `dormancy_growth_penalty` | geschätzt / kalibrieren | `0.55` | Wachstumsreduktion bei Dormanz | Persister-/Dormanzliteratur liefert Richtung, nicht diesen Koeffizienten |
| `synergy_repair_dormancy_bonus` | kalibrieren | `0.25` | Bonus für Reparatur × Dormanz | Modellinteraktion |
| `synergy_stress_tolerance_bonus` | kalibrieren | `0.20` | Bonus für Stressantwort × Toleranz | Modellinteraktion |

### Demografische Stochastik

| Parameter | Typ | Wert | Wirkung | Einschätzung |
|---|---|---:|---|---|
| `stochastic_threshold` | technisch / kalibrieren | `10000` | darunter Poisson-Sampling, darüber Normalrauschen | beeinflusst Aussterben seltener Stämme |
| `stochastic_noise_scale` | technisch / kalibrieren | `0.08` | Rauschstärke oberhalb der Schwelle | Sensitivität prüfen |

---

## Was online identifizierbar ist

Diese Größen können aus Online-Quellen sinnvoll abgeleitet werden:

| Bereich | Parameter / Zielgröße | Warum online identifizierbar | Mögliche Quellen |
|---|---|---|---|
| MRSA-/S.-aureus-Prävalenz | initiale Carrier-Fraktion, `community_carrier_fraction`, Zielprävalenz | Surveillance- und PPS-Daten messen Kolonisation/Infektion direkt | Swissnoso, ECDC EARS-Net, nationale AMR-Berichte |
| MRSA-Anteil / Resistenzanteil | `resistant_fraction`, `replacement_resistant_fraction` | Bei bestätigten MRSA-Carriern liegt MSSA-Anteil bei >90% der Patienten unter Nachweisgrenze (5%); nahezu klonale MRSA-Dominanz → 0.90 | Dall'Antonia et al. 2005, J Hosp Infect |
| Carriage-Dauer / Dekolonisierung | Zielwert für `p_clearance`, `immune_strength`, Clearance-Schwellen | Studien messen Dauer von MRSA-/S.-aureus-Trägerschaft | z. B. Clin. Infect. Dis. 32(10):1393, S.-aureus-Nasal-Carriage-Studien |
| Bakterielle Last | Größenordnung für `initial_population`, `carrying_capacity` | CFU/Swab-Studien berichten S.-aureus-Lasten von sehr niedrig bis Millionen und teils höher | Nasal-load-Studien, z. B. MRSA-Last unter Antibiotika |
| Generationszeit / Replikation | grobe Plausibilität für `steps_per_day`, `growth_rate_per_step` | genomische Studien schätzen S.-aureus-Teilungsraten im menschlichen Nasenraum | BMC Genomics 2019 zu S. aureus im Nasenraum |
| Basale Mutationsrate | Plausibilitätsanker für `base_mutation_rate` | Mutationsraten pro Nukleotid/Generation sind experimentell messbar | BMC Genomics 2019; Mutationsakkumulationsstudien |
| Antibiotika-Stress und SOS | Richtung/Größenordnung für `stress_mutation_boost` | S. aureus-SOS-Antwort und stressinduzierte Mutation sind belegt | Cirz et al. 2007; Ubeda et al. 2006; Reviews zu Stressmutagenese |
| HGT und mobile genetische Elemente | Plausibilität für `base_hgt_rate`, transferierbare Gene | SCCmec, Plasmide, Phagen und Transposons sind gut beschrieben | Reviews zu S. aureus-HGT und Plasmidmobilisierung |
| Dormanz / Persister | Richtung für `dormancy_growth_penalty`, Stress-/Toleranzblock | Persister und Dormanz erklären Antibiotikatoleranz | Reviews zu S.-aureus-Persistern und Dormanz |
| Antibiotikaklassen und Resistenzmechanismen | ABX-Profilstruktur | Wirkmechanismen sind gut bekannt | Mikrobiologie-/Pharmakologiequellen, Guidelines |

Wichtig: Diese Quellen identifizieren meist reale biologische Größen. Die Mikrosimulation verwendet aber ein komprimiertes Trait-Modell. Deshalb liefern Online-Quellen oft Zielgrößen und plausible Bereiche, aber nicht direkt den exakten YAML-Wert.

---

## Was wahrscheinlich kalibriert werden muss

Diese Parameter sollten nicht direkt aus Literatur übernommen werden:

| Parametergruppe | Parameter | Grund |
|---|---|---|
| Trait-Mutation | `base_mutation_rate`, `mutation_std`, `stress_mutation_boost` | echte Mutationen pro Base/Generation müssen auf 14 kontinuierliche Trait-Slots und 12 Schritte/Tag übersetzt werden |
| HGT | `base_hgt_rate`, `hgt_gene_transfer_prob` | reale HGT-Raten sind kontext-, Stamm-, Plasmid- und Nischen-abhängig; das Modell nutzt eine effektive Ereignisrate |
| Selektion | `selection_strength`, `growth_rate_per_step`, `death_rate_per_step` | bestimmt Sweep-Geschwindigkeit und Populationsstabilität im Modell |
| Population/Clearance | `carrying_capacity`, `min_population`, `clearance_threshold` | CFU-Werte sind messbar, aber die Modellpopulation ist eine normalisierte Binnenpopulation mit eigener Clearance-Logik |
| Schaden/Alterung | `base_damage_per_step`, `replication_damage_factor`, `stress_damage_factor`, `repair_rate_per_step`, `age_mortality_scale`, `damage_mortality_scale`, `lifecycle_half_life_steps`, `max_damage_load` | latente Variablen ohne direkte klinische Messgröße |
| Dormanz-Synergien | `dormancy_growth_penalty`, `synergy_repair_dormancy_bonus`, `synergy_stress_tolerance_bonus` | Persister-Dormanz ist qualitativ belegt, aber die Modellinteraktionen sind frei gewählt |
| Stochastik/Pruning | `stochastic_threshold`, `stochastic_noise_scale`, `strain_prune_threshold`, `max_strains` | beeinflusst seltene Varianten und numerische Stabilität; muss per Sensitivitätsanalyse abgesichert werden |
| ABX-Profilzahlen | `ABXProfile.*`, `DOSE_MULTIPLIERS.*` | qualitative Mechanismen sind bekannt, aber `base_kill_rate` und Schutzgewichte sind keine PK/PD-Parameter |

---

## Empfohlene Kalibrierungsziele

Die Mikrokalibrierung sollte nicht jeden Parameter einzeln gegen eine einzelne Zielgröße optimieren. Sinnvoller ist eine stufenweise Kalibrierung mit mehreren Outputs:

1. **Population stabilisieren**
   - Ziel: `mean_total_population` bleibt ohne ABX in plausibler Größenordnung.
   - Parameter: `carrying_capacity`, `growth_rate_per_step`, `death_rate_per_step`, `selection_strength`.

2. **Clearance / Carriage-Dauer treffen**
   - Ziel: mittlere `p_clearance` und simulierte Carrier-Dauer passen zu Literatur.
   - Parameter: `min_population`, `clearance_threshold`, `immune_strength`, `STEALTH`-Effekt, `carrying_capacity`.

3. **Resistenzentwicklung unter ABX**
   - Ziel: `mean_resistant_fraction`, `p50/p90_resistant_fraction`, Zeit bis `R2/R3`-Dominanz.
   - Parameter: `base_mutation_rate`, `mutation_std`, `stress_mutation_boost`, `selection_strength`, ABX-Profile.

4. **HGT/Diversität**
   - Ziel: `mean_n_strains`, `genotype_entropy`, Rohlogs `micro_strain_daily` und `micro_episode_gene_daily`.
   - Parameter: `base_hgt_rate`, `hgt_gene_transfer_prob`, `max_strains`, `strain_prune_threshold`.

5. **Persistenz unter Stress**
   - Ziel: resistente Minderheiten verschwinden nicht unrealistisch schnell, aber Population explodiert nicht.
   - Parameter: Dormanz-, Reparatur-, Schaden- und Stochastikblock.

---

## Quellenanker

- S. aureus-Mutationsrate und Replikation im menschlichen Nasenraum: <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5604-6>
- S. aureus-Nasal-Carriage und Kolonisationsmechanismen: <https://pmc.ncbi.nlm.nih.gov/articles/PMC6186810/>
- Nasale MRSA-Last unter Antibiotika: <https://pubmed.ncbi.nlm.nih.gov/18632184/>
- S. aureus-SOS-Antwort auf Ciprofloxacin: <https://pmc.ncbi.nlm.nih.gov/articles/PMC1797410/>
- Beta-Lactame, SOS-Antwort und Transfer von Virulenzfaktoren in S. aureus: <https://journals.asm.org/doi/abs/10.1128/jb.188.7.2726-2729.2006>
- Plasmid-Konjugation und Mobilisierung in Staphylococcus: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4993578/>
- S. aureus-Persister und Dormanz unter oxidativem Stress: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8865412/>
- Überblick zu Persister-Zellen und Dormanz: <https://journals.asm.org/doi/10.1128/aem.02636-13>
