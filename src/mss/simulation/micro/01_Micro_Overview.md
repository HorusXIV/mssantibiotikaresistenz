# Micro-Ebene der Simulation

## Ziel der Micro-Ebene

Die Micro-Ebene modelliert, was **innerhalb eines einzelnen Patienten** mit einer bakteriellen Population passiert. Im Gegensatz zur Makro-Ebene, die Patientenbewegungen, Transmission zwischen Patienten und Krankenhauslogik abbildet, simuliert die Micro-Ebene:

- Wachstum und Absterben von Bakterien
- Konkurrenz zwischen mehreren Staemmen
- Mutation
- Horizontalen Gentransfer (HGT)
- Selektion durch Antibiotika
- Selektion durch das Immunsystem
- Ableitung von makro-relevanten Groessen wie:
  - `p_clearance`
  - `relative_transmissibility`
  - `severity_modifier`
  - `lethality_modifier`

Die drei zentralen Dateien sind:

- `genome.py`: Genomdarstellung, Fitness-Helfer, Resistenzlogik
- `engine.py`: eigentliche Within-Host-Simulation
- `simulator.py`: Episodenverwaltung und Batch-Verarbeitung fuer viele Patienten

## Grundidee des Modells

Die Micro-Simulation arbeitet nicht mit einem einzelnen Bakterium, sondern mit einer **Population aus Staemmen**. Jeder Stamm hat:

- ein Genom als Vektor aus `NUM_GENES = 14` kontinuierlichen Genwerten
- eine Populationsgroesse
- ein Linienalter (`lineage_ages`)
- eine angesammelte Schadenslast (`damage_loads`)
- einen Namen (`strain_names`)

Ein Tag wird in `steps_per_day = 12` diskrete Schritte zerlegt. In jedem Schritt passiert in dieser Reihenfolge:

1. **Selection Step**: Wachstum und Tod haengen von Fitness, Immunstatus, Antibiotika und Lebenszyklus-Kosten ab.
2. **Mutation**: Neue Varianten koennen durch verrauschte Genveraenderungen entstehen.
3. **HGT**: In jedem dritten Schritt koennen Gene zwischen Staemmen uebertragen werden.
4. **Consolidation**: Sehr kleine Staemme werden entfernt, zu viele Staemme werden abgeschnitten.

Am Ende des Tages wird aus der finalen Population eine Micro-Response gebaut, die an den `Patient` zurueckgeht.

## Architektur und Datenfluss

### 1. Makro-Layer baut den Tageskontext

Der Makro-Simulator erzeugt pro Tag einen `PatientDailyContext`. Relevant dafuer sind vor allem:

- Krankenhaus und Station
- Isolationsstatus
- Antibiotika-Regime

Das passiert in `mss/simulation/macro/simulator.py` in `_build_context()`.

### 2. `patient.py` uebersetzt Makro-Zustand in einen Micro-Request

Die Klasse `Patient` ist die Bruecke zwischen Makro und Micro.

Wichtige Methoden:

- `update_context(ctx)`: uebernimmt Makro-Kontext und leitet z. B. `adherence` ab
- `make_micro_request(...)`: baut die Eingabe fuer die Micro-Simulation
- `apply_micro_response(resp)`: uebernimmt die Antwort der Micro-Simulation

Wichtig: Ein Request wird nur gebaut, wenn der Patient `HealthState.CARRIER` ist und eine `episode_id` hat.

### 3. `MicroSimulator` verarbeitet den Request

`mss/simulation/micro/simulator.py` kuemmert sich um:

- Wiederherstellung des bisherigen Episodenzustands
- Initialisierung einer neuen Population bei neuer Episode
- Aufruf von `simulate_day(...)`
- Rueckgabe einer Response
- Persistenz des EpisodeState zwischen Tagen

### EpisodeState und Persistenz

Der `EpisodeState` ist entscheidend, weil die Micro-Ebene **nicht jeden Tag bei null beginnt**.

Pro Episode werden gespeichert:

- `episode_id`
- `patient_id`
- `day`
- `population.genomes`
- `population.populations`
- `population.lineage_ages`
- `population.damage_loads`
- `population.strain_names`

Das bedeutet:

- Mutationen vom Vortag bleiben erhalten
- HGT-Ereignisse bleiben erhalten
- dominante Staemme koennen sich ueber mehrere Tage etablieren
- Alterungs- und Schadenseffekte akkumulieren ueber mehrere Tage

Ohne diesen Zustand waere die Simulation nur eine taegliche Einzelsimulation ohne echte Within-Host-Evolution.

### 4. Rueckwirkung auf die Makro-Ebene

Die Micro-Response beeinflusst die Makro-Ebene ueber `Patient.apply_micro_response()`:

- `resistant_fraction`
- `dominant_genotype`
- `dominant_strain_name`
- `relative_transmissibility`
- `p_clearance`
- `severity_modifier`
- `lethality_modifier`

Diese Werte wirken danach in der Makro-Ebene weiter:

- `p_clearance` steuert `C -> S` in `Patient.should_clear_today()`
- `relative_transmissibility` skaliert Transmission in `transmission_multiplier_for_macro()`
- `lethality_modifier` skaliert das Todesrisiko in `daily_death_risk_multiplier()`
- `dominant_genotype` und `resistant_fraction` werden bei Uebertragung auf neue Patienten vererbt

## Verbindung zu `patient.py`

### Relevante Patient-Felder mit direktem Einfluss auf die Micro-Ebene

### `state`

- Nur Traeger (`CARRIER`) werden an die Micro-Ebene geschickt.
- `SUSCEPTIBLE` erzeugt keinen Micro-Request.

### `episode_id`

- Identifiziert die laufende Besiedelungs-/Infektionsepisode.
- Dient in `MicroSimulator` als Schluessel fuer persistenten Zustand ueber mehrere Tage.

### `compliance`

- Wird in `update_context()` in `adherence` umgewandelt.
- ICU-Patienten erhalten dort aktuell einen Bonus von `+0.1`, gedeckelt auf `1.0`.
- Hoehere `adherence` verstaerkt die effektive Antibiotikawirkung.

### `immune_strength`

- Geht direkt in die Immun-Selektion ein.
- Hoehere Werte erhoehen die Immunclearance und spaeter auch `p_clearance`.

### `immune_status`

- Relevant sind aktuell `"normal"` und `"suppressed"`.
- `"suppressed"` reduziert die Immun-Clearance stark.

### `resistant_fraction`

- Wird bei neuen Episoden als Anteil resistenter Startpopulation interpretiert.
- Bestimmt, wie viel der Anfangspopulation aus resistenten Seed-Staemmen besteht.

### `dominant_genotype`

- Gibt die resistente Startklasse vor: `"S"`, `"R1"`, `"R2"` oder `"R3"`.
- Steuert, welches Seed-Genom bei resistenten Anfangsstaemmen erzeugt wird.

### `dominant_strain_name`

- Wenn vorhanden, wird dieser Name bevorzugt fuer den dominanten Seed-Stamm weiterverwendet.

### Patient-Felder ohne direkten Micro-Effekt in der aktuellen Implementierung

Diese Felder werden zwar in den Request geschrieben oder in `Patient` gepflegt, beeinflussen die Engine derzeit aber **nicht direkt**:

- `age_years`
- `vulnerability`
- `history_flags`
- `sociability`
- `is_isolated`
- `hospital_id`
- `treatment_phase`
- `severity_modifier`
- `lethality_modifier`

Davon wirken manche nur auf der Makro-Ebene:

- `sociability` beeinflusst Transmission im Makro-Layer
- `vulnerability` beeinflusst Makro-Suszeptibilitaet und Todesrisiko
- `is_isolated` beeinflusst Makro-Transmission und Entlassungslogik

## Verbindung zum Makro-Layer

### Direkte Einfluesse von Makro auf Micro

### Station (`Department`)

Der Request enthaelt `setting = self.department.value`, aber:

- das Feld wird aktuell in der Micro-Engine **nicht verwendet**
- indirekt wirkt die Station trotzdem:
  - auf ICU wird `adherence` in `Patient.update_context()` erhoeht
  - ICU hat im Makro eine hoehere Wahrscheinlichkeit fuer Antibiotika

### Antibiotika-Politik

Makro erzeugt ueber `_build_context()` taeglich ein `AntibioticRegimen`:

- `regimen.on`
- `regimen.abx_class`
- `regimen.dose_level`

Diese Werte gehen direkt in die Micro-Selektion ein und bestimmen:

- welche `ABXProfile` genutzt werden
- wie hoch die Kill-Rate ist
- welche Resistenzmechanismen Schutz geben

### Isolation, Hygiene und Diagnostik

Im `PatientDailyContext` existieren:

- `hygiene_level`
- `isolation_effectiveness`
- `diagnostic_speed`
- `is_isolated`

Aber in der aktuellen Implementierung gilt:

- `hygiene_level` wird **nicht** an Micro uebergeben
- `isolation_effectiveness` wird **nicht** an Micro uebergeben
- `diagnostic_speed` wird **nicht** an Micro weitergereicht
- `is_isolated` wird zwar im Patient gespeichert, beeinflusst die Micro-Engine aber nicht

Diese Groessen sind momentan also primaer Makro-Parameter.

### Direkte Rueckkopplung von Micro auf Makro

### `relative_transmissibility`

- wird aus dem dominanten Stamm berechnet
- skaliert in Makro den Beitrag eines Traegers zur Transmission

### `p_clearance`

- wird aus Gesamtpopulation, Immunsystem und Stealth berechnet
- wird am naechsten Makro-Tag fuer die Entscheidung `C -> S` verwendet

### `severity_modifier` und `lethality_modifier`

- werden aus dem dominanten Stamm abgeleitet
- koennen in Makro Schweregrad und Mortalitaet skalieren

### `dominant_genotype` und `resistant_fraction`

- werden bei Transmission auf neue Patienten uebernommen
- Makro kann dabei zusaetzlich kleine Mutation/Drift auf Resistenz und Genotyp anwenden

## Der Request an die Micro-Ebene

Ein Request aus `Patient.make_micro_request()` hat diese Struktur:

```python
{
    "schema_version": "1.0",
    "run_id": ...,
    "episode_id": ...,
    "patient_id": ...,
    "t_day": ...,
    "dt_days": 1,
    "setting": ...,
    "abx": {
        "on": ...,
        "class": ...,
        "dose_level": ...,
    },
    "adherence": ...,
    "host": {
        "age_years": ...,
        "immune_strength": ...,
        "immune_status": ...,
        "vulnerability": ...,
        "history_flags": ...,
    },
    "initial_state": {
        "resistant_fraction": ...,
        "dominant_genotype": ...,
        "dominant_strain_name": ...,
    },
    "seed": ...,
}
```

## Initialisierung einer neuen Episode

Wenn zu einer `episode_id` noch kein gespeicherter `EpisodeState` existiert, erzeugt `StrainPopulation.create_initial(...)` eine Startpopulation.

Wichtige Default-Werte:

- `initial_population = 1e6`
- `n_susceptible_strains = 3`
- `n_resistant_strains = 2`

Die Startpopulation wird dann so aufgeteilt:

- Anteil `1 - resistant_fraction` auf empfindliche Staemme
- Anteil `resistant_fraction` auf resistente Staemme

Die Genomwerte werden anschliessend noch leicht verrauscht, damit nicht alle Seed-Staemme exakt identisch sind.

Wichtig fuer die Erklaerung:

- `dominant_genotype` bestimmt, **welche Art resistenter Seed-Staemme** erzeugt wird
- `resistant_fraction` bestimmt, **wie viel der Anfangspopulation** resistent ist
- beides zusammen bestimmt also Startlage und Selektionsreserve

## Welche Request-Felder werden aktuell wirklich verwendet?

### Aktiv genutzt

- `episode_id`
- `patient_id`
- `t_day`
- `seed`
- `abx.on`
- `abx.class`
- `abx.dose_level`
- `adherence`
- `host.immune_strength`
- `host.immune_status`
- `initial_state.resistant_fraction`
- `initial_state.dominant_genotype`
- `initial_state.dominant_strain_name`

### Aktuell nur transportiert, aber nicht in der Engine genutzt

- `schema_version`
- `run_id`
- `dt_days`
- `setting`
- `host.age_years`
- `host.vulnerability`
- `host.history_flags`

Das ist wichtig fuer muendliche Erklaerungen: Im Datenmodell sind diese Felder schon vorgesehen, die aktuelle Engine verarbeitet aber nur einen Teil davon.

## Genom und "Allele"

### Wichtige Praezisierung

Im Code gibt es **keine expliziten diskreten Allele im klassischen mendelschen Sinn**. Es gibt also keine Objekte wie `"Allel A"` oder `"Allel a"`, keine Diploidie und keine festen Basenfolgen.

Stattdessen wird jedes "Gen" als **kontinuierlicher Zahlenwert zwischen `0.0` und `1.0`** modelliert.

Du kannst dir das so merken:

- `0.0` = Eigenschaft praktisch nicht ausgepraegt
- `1.0` = Eigenschaft sehr stark ausgepraegt

Wenn du nach "Allelen" gefragt wirst, ist die praezise Antwort fuer dieses Modell:

> Ein "Allel" ist hier keine diskrete Kategorie, sondern die aktuelle Auspraegung eines Gen-Slots als kontinuierlicher Wert im Genomvektor.

Mutation veraendert diesen Wert durch Gauss-Rauschen, HGT mischt Werte zwischen Donor und Empfaenger.

## Gene des Genoms

Das Genom hat `14` Slots. Die Indizes sind in `GeneIndex` definiert.

| Index | Gen | Bedeutung | Hauptwirkung |
|---|---|---|---|
| 0 | `GROWTH_BASE` | Basiswachstum | Hoeherer Grundwert steigert Fitness |
| 1 | `METABOLIC_OPTIMIZATION` | Kompensation metabolischer Kosten | Senkt Resistenzkosten |
| 2 | `EFFLUX_PUMPS` | Auspumpen von Antibiotika | Schuetzt vor mehreren ABX-Klassen |
| 3 | `TARGET_MODIFICATION` | Veraenderung des Angriffsziels | Starker Schutz bei zielgerichteten ABX |
| 4 | `PERMEABILITY_REDUCTION` | Geringere Zellpermeabilitaet | Senkt ABX-Eindringen |
| 5 | `VIRULENCE` | Virulenz | Erhoeht Severity und Lethality, leicht auch Transmission |
| 6 | `STEALTH` | Immunevasion | Reduziert Immun-Clearance und senkt `p_clearance` |
| 7 | `ADHESION` | Anhaftung/Haftfaehigkeit | Erhoeht Transmission und Severity |
| 8 | `MUTATION_RATE_MODIFIER` | Evolvierbarkeit | Erhoeht stammspezifische Mutationsrate |
| 9 | `HGT_COMPETENCE` | Aufnahmefaehigkeit fuer HGT | Erhoeht HGT-Wahrscheinlichkeit |
| 10 | `DNA_REPAIR` | Reparaturkapazitaet | Senkt Schadenslast und Alterungsdruck |
| 11 | `DORMANCY_PROPENSITY` | Neigung zur Dormanz | Weniger Wachstum, aber mehr Ueberleben unter Stress |
| 12 | `STRESS_RESPONSE` | Antwort auf Umweltstress | Unterstuetzt Ueberleben, Reparatur und Toleranz |
| 13 | `DAMAGE_TOLERANCE` | Toleranz gegen Schaeden | Senkt turnover-/schadensbedingten Tod |

## Wildtyp und resistente Startgenome

### `create_wild_type_genome()`

Der Wildtyp startet:

- mit gutem Basiswachstum
- mit niedrigen Resistenzwerten
- mit moderater Virulenz/Adhaesion/Stealth

Das ist ein eher empfindlicher, aber konkurrenzfaehiger Ausgangsstamm.

### `create_resistant_genome(resistance_level)`

Ein resistentes Genom:

- erhoeht `EFFLUX_PUMPS`
- erhoeht `TARGET_MODIFICATION`
- erhoeht `PERMEABILITY_REDUCTION`
- erhoeht Kompensations- und Stressgene
- senkt leicht `GROWTH_BASE`

Die Idee dahinter:

- Resistenz schuetzt vor Antibiotika
- Resistenz kostet aber typischerweise Wachstum
- Kompensationsgene federn einen Teil dieser Kosten ab

### Startgenome fuer `R1`, `R2`, `R3`

Bei der Initialisierung wird fuer resistente Seed-Staemme nicht nur ein generisches Resistenzlevel genutzt, sondern `_create_seed_genome_for_genotype()`:

- `S`: Wildtyp bzw. niedrige Resistenz
- `R1`: leichte Resistenz
- `R2`: mittlere Resistenz mit deutlich hoeheren Stress-/Persistenzwerten
- `R3`: starke Resistenz, groessere Kosten, aber starke Schutzmechanismen

## Genotyp-Klassifikation

`classify_genotype(genome)` bildet aus einem Genom eine diskrete Klasse:

- `S`
- `R1`
- `R2`
- `R3`

Die Klasse basiert auf dem Mittelwert aus:

- `EFFLUX_PUMPS`
- `TARGET_MODIFICATION`
- `PERMEABILITY_REDUCTION`

Schwellen:

- `< 0.2` -> `S`
- `< 0.4` -> `R1`
- `< 0.7` -> `R2`
- `>= 0.7` -> `R3`

Wichtig: Der Genotyp ist also **kein eigener Gen-Slot**, sondern eine **abgeleitete Klassifikation** aus drei Resistenzgenen.

## Resistenzanteil der Population

`compute_resistant_fraction(...)` berechnet den Anteil der Gesamtpopulation, deren Resistenzscore `>= 0.3` ist.

Das bedeutet:

- Ein Patient kann `dominant_genotype = "S"` haben
- gleichzeitig aber trotzdem einen kleinen resistenten Subanteil in `resistant_fraction`

Das ist biologisch sinnvoll, weil Dominanz und Teilpopulationen nicht dasselbe sind.

## Antibiotika-Logik

In `ABX_PROFILES` ist festgelegt, wie gut verschiedene Resistenzmechanismen gegen einzelne Antibiotikaklassen wirken.

| Klasse | Efflux | Target Mod | Permeability | Base Kill |
|---|---|---|---|---|
| `none` | 0.0 | 0.0 | 0.0 | 0.0 |
| `beta_lactam` | 0.3 | 0.8 | 0.4 | 0.75 |
| `fluoroquinolone` | 0.6 | 0.7 | 0.3 | 0.80 |
| `aminoglycoside` | 0.4 | 0.5 | 0.6 | 0.70 |
| `macrolide` | 0.7 | 0.4 | 0.3 | 0.65 |
| `tetracycline` | 0.8 | 0.3 | 0.2 | 0.60 |
| `glycopeptide` | 0.2 | 0.9 | 0.5 | 0.85 |

Zusaetzlich gibt es Dosis-Multiplikatoren:

- `low` -> `0.6`
- `std` -> `1.0`
- `high` -> `1.4`

Die effektive Kill-Rate ist:

```text
effective_kill = base_kill_rate * dose_multiplier * adherence
```

Der Resistenzschutz wird aus den drei Resistenzgenen aufaddiert, profilabhaengig gewichtet, dann normalisiert und auf `0.95` begrenzt.

Die resultierende Ueberlebenswahrscheinlichkeit unter ABX ist:

```text
survival = 1 - effective_kill * (1 - protection)
```

Interpretation:

- hohe Adherence -> Antibiotika wirken staerker
- hohe Dosis -> Antibiotika wirken staerker
- passende Resistenzmechanismen -> Schutz steigt

## Fitnesslogik

Die zentrale Fitnessfunktion ist:

```text
fitness = (growth_base - net_costs) * abx_survival * immune_survival
```

mit Begrenzung auf sinnvolle Werte.

### Resistenzkosten

`compute_resistance_costs(...)` verwendet `ResistanceCosts`:

- `efflux_pumps = 0.15`
- `target_modification = 0.12`
- `permeability_reduction = 0.08`

Die Rohkosten sind die gewichtete Summe der Resistenzgene.

`METABOLIC_OPTIMIZATION` reduziert diese Kosten um bis zu `80%`.

Interpretation:

- Resistenz ist nicht kostenlos
- Kompensationsgene machen resistente Staemme wieder konkurrenzfaehiger

### Immunsystem

`compute_immune_survival(...)` berechnet:

- Basisclearance = `0.15 * immune_strength`
- bei `immune_status == "suppressed"` nur `30%` davon
- `STEALTH` reduziert die Erkennung um bis zu `70%`

Hoeheres `STEALTH` bedeutet also:

- mehr Ueberleben unter Immunangriff
- spaeter auch geringere Clearance-Wahrscheinlichkeit des Patienten

## Ablauf eines Selection Steps

`selection_step(...)` ist der biologisch wichtigste Teil.

### 1. Fitness berechnen

Fuer jeden Stamm wird die Fitness unter aktuellen Umweltbedingungen berechnet.

### 2. Relative Selektion

Nicht die absolute Fitness allein ist entscheidend, sondern auch die Fitness relativ zum Populationsmittel:

```text
relative_fitness = fitness / mean_fitness
selection_factor = relative_fitness ** selection_strength
```

Hoeheres `selection_strength` bedeutet:

- Gewinner wachsen schneller
- Verlierer verlieren schneller Population

### 3. Wachstum

Grundwachstum:

```text
growth = growth_rate_per_step * selection_factor * fitness
```

Dann greifen Lebenszyklusgene ein:

- `DNA_REPAIR`
- `DORMANCY_PROPENSITY`
- `STRESS_RESPONSE`
- `DAMAGE_TOLERANCE`
- `METABOLIC_OPTIMIZATION`
- `STEALTH`

Es werden daraus drei Kapazitaeten abgeleitet:

- `repair_capacity`
- `dormancy_capacity`
- `tolerance_capacity`

### 4. Dormanz

`active_dormancy` steigt bei hohem Umweltstress.

Effekt:

- Wachstum sinkt
- Replikationsdruck sinkt
- Ueberleben unter Stress kann steigen

Das ist ein klassischer Trade-off:

- langsamere Expansion
- dafuer robustere Persistenz

### 5. Schadensaufbau und Reparatur

Schadenszuwachs kommt aus:

- Grundverschleiss
- Replikationsdruck
- Umweltstress

Schadensabbau kommt aus:

- `DNA_REPAIR`
- Dormanz-/Repair-Synergie
- Stress-/Toleranz-Synergie

Hohe `damage_loads` erhoehen spaeter den Tod.

### 6. Linienalter und Turnover

`lineage_ages` steigt mit der Zeit und mit Replikationsdruck.

Aus Alter und Schaden wird `turnover_pressure` berechnet. Das ist ein Zusatzterm auf die Todesrate.

### 7. Tod

Die Todesrate kombiniert:

- `death_rate_per_step`
- inverse Fitness
- Alterungs-/Schadensdruck

### 8. Demographische Stochastik

Kleine Populationen werden nicht nur deterministisch gerechnet:

- unter `stochastic_threshold` wird Poisson-Sampling verwendet
- darueber Normalrauschen mit Standardabweichung `sqrt(expected) * stochastic_noise_scale`

Dadurch koennen kleine Staemme zufaellig verschwinden.

### 9. Carrying Capacity

Wenn die Gesamtpopulation groesser als `carrying_capacity` wird, werden alle Populationen proportional nach unten skaliert.

## Mutation

`mutate_population(...)` erzeugt neue Staemme.

### Mutationsrate

Die effektive Mutationsrate pro Stamm ist:

```text
effective_rate = base_mutation_rate * stress_multiplier
strain_rate = effective_rate * (0.5 + MUTATION_RATE_MODIFIER)
```

mit:

```text
stress_multiplier = 1 + abx_stress * (stress_mutation_boost - 1)
```

Interpretation:

- Antibiotika-Stress erhoeht Mutationsraten
- Staemme mit hohem `MUTATION_RATE_MODIFIER` evolvieren schneller

### Art der Mutation

- Anzahl Mutationen pro Schritt: Poisson
- betroffene Gene: zufaellige Auswahl ohne Zuruecklegen
- Aenderung pro Gen: Gauss-Verteilung mit `mutation_std`
- neue Genwerte werden auf `[0.0, 1.0]` begrenzt

Wenn eine Mutation stattfindet und der Stamm gross genug ist:

- ein kleiner Teil der Population wird abgespalten
- daraus entsteht ein neuer Stamm
- Alter und Schadenslast werden teilweise "verjuengt"

## Horizontaler Gentransfer (HGT)

`horizontal_gene_transfer(...)` laesst Staemme Gene von anderen Staemmen uebernehmen.

### HGT-Wahrscheinlichkeit

Sie haengt ab von:

- `base_hgt_rate`
- `HGT_COMPETENCE`

Form:

```text
hgt_prob = base_hgt_rate * (0.5 + HGT_COMPETENCE)
```

### Transferierbare Gene

Nicht alle Gene werden per HGT uebertragen. Transferierbar sind:

- `EFFLUX_PUMPS`
- `TARGET_MODIFICATION`
- `PERMEABILITY_REDUCTION`
- `METABOLIC_OPTIMIZATION`
- `DORMANCY_PROPENSITY`
- `STRESS_RESPONSE`
- `DAMAGE_TOLERANCE`

Das ist biologisch sinnvoll, weil vor allem Resistenz- und Persistenzmechanismen mobil sind.

### Transfermechanik

- Donor wird populationsgewichtet gezogen
- pro transferierbarem Gen entscheidet `hgt_gene_transfer_prob`, ob es uebertragen wird
- der Empfaengerwert wird nicht hart ersetzt, sondern mit dem Donorwert gemischt

Auch hier entsteht ein neuer rekombinanter Teilstamm.

## Clearance-Wahrscheinlichkeit des Patienten

`compute_clearance_probability(...)` erzeugt die taegliche Wahrscheinlichkeit, dass der Patient im Makro-Layer wieder von `CARRIER` zu `SUSCEPTIBLE` wechselt.

Einflussfaktoren:

- Gesamtpopulation
- `clearance_threshold`
- `min_population`
- `carrying_capacity`
- `immune_strength`
- `immune_status`
- durchschnittliches `STEALTH` der Population

Logik:

- sehr kleine Population -> hohe Clearance
- grosse Population -> niedrige Clearance
- starkes Immunsystem -> hoehere Clearance
- hohe Stealth-Werte -> niedrigere Clearance

## Response der Micro-Ebene

`population_to_response(...)` erzeugt drei Bloecke:

### `updated_state`

- `resistant_fraction`
- `dominant_genotype`
- `dominant_strain_name`

### `derived_effects`

- `relative_transmissibility`
- `lethality_modifier`
- `severity_modifier`
- `p_clearance`

### `population_stats`

- `total_population`
- `n_strains`
- `dominant_strain_name`

Die `population_stats` werden aktuell in `Patient.apply_micro_response()` nicht genutzt, koennen aber fuer Debugging oder spaetere Erweiterungen hilfreich sein.

## Alle Parameter der Micro-Ebene und ihr Einfluss

Gemeint ist hier die `SimulationConfig` aus `engine.py`.

| Parameter | Default | Einfluss |
|---|---:|---|
| `steps_per_day` | `12` | Mehr Schritte pro Tag = feinere Dynamik, mehr Gelegenheiten fuer Selektion, Mutation und HGT |
| `max_strains` | `50` | Begrenzt, wie viele unterschiedliche Staemme gehalten werden; zu klein kann Diversitaet abschneiden |
| `carrying_capacity` | `1e9` | Obergrenze der Gesamtpopulation; kleiner = staerkerer logistischer Druck |
| `min_population` | `1e3` | Unterhalb davon wird Clearance deutlich wahrscheinlicher |
| `clearance_threshold` | `1e2` | Unterhalb davon ist Clearance fast sicher |
| `base_mutation_rate` | `0.01` | Grundwahrscheinlichkeit fuer Mutationen pro Gen und Schritt |
| `mutation_std` | `0.05` | Staerke einzelner Mutationsschritte; groesser = mutationale Spruenge groesser |
| `stress_mutation_boost` | `3.0` | Verstaerkt Mutation unter Antibiotikastress |
| `base_hgt_rate` | `0.02` | Grundchance fuer HGT pro Stamm und Schritt |
| `hgt_gene_transfer_prob` | `0.3` | Wahrscheinlichkeit, dass ein einzelnes transferierbares Gen wirklich uebernommen wird |
| `selection_strength` | `2.0` | Verstaerkt Unterschiede zwischen fitten und unfitten Staemmen |
| `growth_rate_per_step` | `0.3` | Oberes Wachstumslimit pro Schritt vor den Korrekturen |
| `death_rate_per_step` | `0.1` | Grundsterberate pro Schritt |
| `strain_prune_threshold` | `1.0` | Staemme unterhalb dieser Groesse werden entfernt |
| `base_damage_per_step` | `0.03` | Grundschaeden unabhaengig von Replikation und Stress |
| `replication_damage_factor` | `0.25` | Zusatzschaeden durch starke Replikation |
| `stress_damage_factor` | `0.20` | Zusatzschaeden durch unguenstige Umwelt |
| `repair_rate_per_step` | `0.08` | Wie stark Reparaturmechanismen Schaden wieder abbauen |
| `age_mortality_scale` | `0.015` | Einfluss des Linienalters auf Sterblichkeit |
| `damage_mortality_scale` | `0.20` | Einfluss der Schadenslast auf Sterblichkeit |
| `lifecycle_half_life_steps` | `36.0` | Skala, ab wann Linienalter relevant wird |
| `max_damage_load` | `5.0` | Obergrenze/Saettigung fuer Schadenslast |
| `dormancy_growth_penalty` | `0.25` | Wie stark Dormanz Wachstum reduziert |
| `synergy_repair_dormancy_bonus` | `0.25` | Zusatznutzen, wenn Reparatur und Dormanz gemeinsam hoch sind |
| `synergy_stress_tolerance_bonus` | `0.20` | Zusatznutzen, wenn Stressantwort und Toleranz gemeinsam hoch sind |
| `stochastic_threshold` | `1e5` | Unterhalb davon wird exakte Poisson-Stochastik genutzt |
| `stochastic_noise_scale` | `1.0` | Staerke des Rauschens bei grossen Populationen |

## Welche Parameter sind fuer das Verhalten besonders wichtig?

Wenn du die Dynamik schnell erklaeren musst, sind diese Parameter die Schluesselhebel:

### Resistenzentwicklung

- `base_mutation_rate`
- `mutation_std`
- `stress_mutation_boost`
- `base_hgt_rate`
- `hgt_gene_transfer_prob`

### Selektionsdruck

- `selection_strength`
- `growth_rate_per_step`
- `death_rate_per_step`
- Antibiotikaklasse
- Dosis
- Adherence
- `immune_strength`
- `immune_status`

### Persistenz und Clearance

- `min_population`
- `clearance_threshold`
- `carrying_capacity`
- `repair_rate_per_step`
- `damage_mortality_scale`
- `dormancy_growth_penalty`

## Was die Micro-Ebene explizit nicht macht

Die aktuelle Implementierung modelliert nicht:

- einzelne Bakterienzellen
- explizite DNA-Sequenzen
- diskrete klassische Allelobjekte
- Pharmakokinetik ueber den Tag
- unterschiedliche Gewebe-/Organraeume
- direkte Wirkung von Hygiene oder Isolation innerhalb des Wirts
- Nutzung von `age_years`, `history_flags` oder `vulnerability` in der Engine

## Ein typischer Tag eines Traeger-Patienten

1. Makro erkennt den Patienten als `CARRIER`.
2. Makro baut den `PatientDailyContext`.
3. `Patient.update_context()` setzt Kontext, Regime und `adherence`.
4. `Patient.make_micro_request()` baut den Tagesrequest.
5. `MicroSimulator` laedt den bisherigen EpisodeState oder erzeugt eine neue Population.
6. `simulate_day()` fuehrt 12 Schritte aus.
7. `population_to_response()` berechnet dominanten Stamm und Makro-Effekte.
8. `Patient.apply_micro_response()` schreibt die Werte in den Patient zurueck.
9. Makro nutzt diese Werte fuer Clearance, Transmission und weitere Entscheidungen.

## Merksaetze fuer die Erklaerung

- Die Makro-Ebene sagt: **In welchem Patienten und unter welchen Rahmenbedingungen?**
- Die Micro-Ebene sagt: **Welcher Stamm setzt sich innerhalb dieses Patienten durch und mit welchen Folgen?**
- Das Genom ist hier kein Sequenzmodell, sondern ein **Trait-Vektor**.
- "Allele" sind in diesem Modell am besten als **kontinuierliche Auspraegungen von Traits** zu verstehen.
- `dominant_genotype` ist eine **abgeleitete Resistenzklasse**, kein eigenes Gen.
- `resistant_fraction` beschreibt die gesamte Population, nicht nur den dominanten Stamm.

## Kurzfazit

Die Micro-Ebene ist ein zustandsbehaftetes Within-Host-Modell fuer bakterielle Evolution. Sie nimmt pro Traegerpatient die aktuellen Selektionsbedingungen aus Makro und `patient.py` entgegen, entwickelt ueber 12 Schritte die Bakterienpopulation weiter und liefert danach genau die Groessen zurueck, die die Makro-Ebene fuer Transmission, Clearance und klinische Schwere braucht.
