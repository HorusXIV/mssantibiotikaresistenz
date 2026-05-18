# Makro-Ebene: Dokumentation

## Übersicht

Die Makro-Ebene simuliert die **Ausbreitung von Antibiotikaresistenz auf Spitalebene**. Sie modelliert:
- Aufnahme, Entlassung und Verlegung von Patienten zwischen Spitälern
- Stochastische Übertragung von Keimen zwischen Patienten (S → C)
- Spontane Erholung von Trägern (C → S)
- Antibiotikaverordnung und Isolationsmassnahmen
- Den täglichen Datenaustausch mit der Mikro-Ebene

Der Zustand eines Patienten ist auf dieser Ebene binär: **SUSCEPTIBLE (S)** oder **CARRIER (C)**.

---

## Dateien

| Datei | Zweck |
|---|---|
| `config.py` | `SimulationConfig` — alle Parameter der Makro-Ebene |
| `simulator.py` | `MacroSimulator` — Tagesablauf, Formeln, Logik |
| `grid.py` | `HospitalDepartmentGrid`, `HospitalNetworkGrid` — räumliche Modellierung |
| `agents.py` | `PatientAgent` — Mesa-Wrapper für Patienten im Grid |

---

## Tagesablauf (`MacroSimulator.step`)

Jeden simulierten Tag werden die folgenden Schritte **in dieser Reihenfolge** ausgeführt:

```
1. Pro Spital (im Loop):
   a. Spontane Erholung (C → S)
   b. Kontext-Update + Isolationserkennung
   c. Mikro-Requests sammeln
2. Entlassung (inkl. Mortalität)       pro Spital
3. Aufnahme neuer Patienten            global (Poisson)
4. Mikro-Batch anwenden                alle Carrier global
5. Übertragung (S → C)                 pro Spital
6. Verlegungen zwischen Spitälern      global
```

Die Reihenfolge ist bewusst gewählt: Clearance und Kontext-Update laufen gemeinsam im Spital-Loop, bevor Entlassung und Aufnahme verarbeitet werden. Wer heute gesund wird (1a), ist bei der Entlassung (2) bereits als Susceptible sichtbar.

---

## Parameter

### Infektionskontrolle

#### `base_hygiene`
- **Typ:** float, [0, 1]
- **Einfluss:** Reduziert die Übertragungswahrscheinlichkeit für alle Patienten im Spital. Wirkt als `hygiene_factor = 1 - base_hygiene` direkt in der Übertragungsformel.
- **Formel:** `hazard = daily_contact_attempts × hygiene_factor × ...`
- **Beispiel:** `0.65` → 35% der Kontakte sind infektiös wirksam. `0.9` → nur 10%.
- **Verbindung:** Wird täglich via `PatientDailyContext.hygiene_level` an jeden Patienten übergeben. `_do_transmission` liest den Wert aus dem Kontext des empfänglichen Patienten.

#### `base_isolation_effectiveness`
- **Typ:** float, [0, 1]
- **Einfluss:** Wie stark Isolation die Übertragung reduziert. Wirkt auf **beide Seiten**: isolierte Carrier senden weniger, isolierte Susceptible empfangen weniger.
- **Formel:** `contrib_carrier *= (1 - isolation_effectiveness)` wenn Carrier isoliert; `hazard *= (1 - isolation_effectiveness)` wenn Susceptible isoliert.
- **Beispiel:** `0.65` → Isolation reduziert Übertragung um 65%.
- **Verbindung:** Wird via `PatientDailyContext.isolation_effectiveness` pro Patient gespeichert; `_do_transmission` liest den Wert aus dem Kontext des jeweiligen Patienten.

#### `base_diagnostic_speed`
- **Typ:** float, > 0
- **Einfluss:** Multiplikator auf die tägliche Erkennungswahrscheinlichkeit von Trägern. Beeinflusst, wie schnell Carrier erkannt und isoliert werden.
- **Formel:** `effective_detection_prob = min(1.0, carrier_isolation_probability × base_diagnostic_speed)`
- **Beispiel:** `1.0` (Standard) → Erkennungsrate unverändert. `0.5` → halbierte Rate. `2.0` → doppelte Rate (max. 1.0).
- **Verbindung:** Wird in `_build_context` berechnet und bestimmt ob `is_isolated = True` gesetzt wird, was dann via `PatientDailyContext` in `patient.is_isolated` landet.

---

### Übertragung

#### `base_transmission_rate`
- **Typ:** float
- **Einfluss:** Basis-Übertragungswahrscheinlichkeit pro Kontakt und pro Carrier. Multipliziert mit dem patientenspezifischen `transmission_multiplier_for_macro()` (Soziabilität × Transmissibilität des Stammes).
- **Formel:** `contrib = base_transmission_rate × patient.transmission_multiplier_for_macro() × proximity_weight`
- **Verbindung zu Patient:** `Patient.relative_transmissibility` (vom Mikro-Layer gesetzt) und `Patient.sociability` skalieren diesen Wert.

#### `daily_contact_attempts`
- **Typ:** float
- **Einfluss:** Anzahl effektiver täglicher Kontakte, die ein Patient mit anderen hat. Multipliziert die gesamte Übertragungsgefährdung.
- **Formel:** `hazard = daily_contact_attempts × hygiene_factor × (weighted_force / n_total) × susceptibility`
- **Hinweis:** `base_transmission_rate × daily_contact_attempts` ergibt den kombinierten Übertragungskoeffizienten. R0-Schätzung: `R0 ≈ daily_contact_attempts × (1 - hygiene) × base_transmission_rate × susceptibles / p_clearance`

**Die vollständige Übertragungsformel:**
```
weighted_force = Σ [ base_transmission_rate × transmissibility × proximity_weight × (1 - iso_eff)^(carrier_isolated) ]

hazard = daily_contact_attempts × (1 - hygiene) × (weighted_force / n_total) × susceptibility × (1 - iso_eff)^(susceptible_isolated)

p_colonize = 1 - exp(-max(0, hazard))
```

---

### Antibiotikaverordnung

#### `icu_abx_probability`
- **Typ:** float, [0, 1]
- **Einfluss:** Tägliche Wahrscheinlichkeit, dass ein ICU-Patient Antibiotika erhält.
- **Verbindung zum Mikro-Layer:** Das zugewiesene `AntibioticRegimen` (Klasse, Dosierung) wird via `PatientDailyContext` → `patient.regimen` gespeichert und fliesst in den Mikro-Request als `"abx": {"on": True, "class": ..., "dose_level": ...}`. Die Mikro-Ebene berechnet daraus den Selektionsdruck auf die Bakterienpopulation.

#### `ward_abx_probability`
- **Typ:** float, [0, 1]
- **Einfluss:** Tägliche Wahrscheinlichkeit für ABX auf der Normalstation. Typischerweise deutlich tiefer als ICU.
- **Verbindung:** Identisch wie `icu_abx_probability`, aber für Ward-Patienten.

> Antibiotika-Klasse und Dosierung werden **zufällig** aus `["beta_lactam", "fluoroquinolone", "glycopeptide", "macrolide", "aminoglycoside"]` bzw. `["low", "std", "high"]` gewählt. Die Verteilung ist gleichmässig.

---

### Erkennung und Isolation

#### `carrier_isolation_probability`
- **Typ:** float, [0, 1]
- **Einfluss:** Basiswahrscheinlichkeit, dass ein nicht-isolierter Carrier an einem Tag erkannt und isoliert wird. Wird durch `base_diagnostic_speed` skaliert.
- **Formel:** `effective_detection = min(1.0, carrier_isolation_probability × base_diagnostic_speed)`
- **Konsequenz bei Erkennung:** `patient.is_isolated = True` → Entlassung wird um `carrier_extension_days` hinausgezögert; Übertragungskoeffizient sinkt um `isolation_effectiveness`.
- **Verbindung zu Patient:** `Patient.is_isolated` bleibt `True` für alle Folgetage (wird nie zurückgesetzt, ausser der Patient wird entlassen).

---

### Übertragungsvererbung

Wenn ein Susceptible durch einen Carrier kolonisiert wird, erbt der neue Carrier den Stamm des Quell-Carriers — mit einer kleinen Zufallskomponente.

#### `transmission_mutation_probability`
- **Typ:** float, [0, 1]
- **Einfluss:** Wahrscheinlichkeit, dass der übertragene Stamm bei der Kolonisierung eine kleine Mutation erfährt.
- **Formel:** `if rng.random() < transmission_mutation_probability: delta = gauss(0, std)`

#### `transmission_resistance_mutation_std`
- **Typ:** float
- **Einfluss:** Standardabweichung des Gausschen Mutationsdeltas auf den `resistant_fraction`-Wert. Grösserer Wert → grössere Drifts bei jeder Übertragung.
- **Verbindung:** Beide Parameter steuern `_inherit_transmitted_state()` in `simulator.py`. Der resultierende `resistant_fraction` und `dominant_genotype` wird auf den neuen Carrier gesetzt und beim nächsten Tag an den Mikro-Layer übergeben.

---

### Verweildauer (Length of Stay)

Die Verweildauer wird bei Aufnahme einmalig aus einer **Log-Normalverteilung** gezogen.

#### `los_mean_ward`
- **Typ:** float (Tage)
- **Einfluss:** Mittlere Verweildauer auf der Normalstation. Bestimmt wann Patienten ihren geplanten Entlassungstag erreichen.

#### `los_mean_icu`
- **Typ:** float (Tage)
- **Einfluss:** Mittlere Verweildauer auf der ICU. Typischerweise länger als Ward.

#### `los_sigma`
- **Typ:** float
- **Einfluss:** Form-Parameter der Log-Normalverteilung. Grösser → breitere Streuung der Verweildauer.
- **Formel:** `mu = log(mean) - 0.5 × sigma²; LOS = max(1, round(lognormal(mu, sigma)))`

---

### Entlassung

Patienten werden nicht exakt am geplanten Entlassungstag entlassen, sondern mit einer **logistischen Wahrscheinlichkeit** nach Überschreitung dieses Tages.

#### `carrier_extension_days`
- **Typ:** float (Tage)
- **Einfluss:** Wird an zwei Stellen verwendet:
  1. **Bei Ersterkennung** (Transition nicht-isoliert → isoliert): `planned_discharge_day = erkennungstag + carrier_extension_days` (ohne Skalierung).
  2. **Rollierender Discharge-Loop**: Solange ein Carrier noch aktiv ist und seinen Entlassungstag überschritten hat, wird das Datum täglich neu gesetzt: `planned_discharge_day = heute + carrier_extension_days × severity_modifier`. Schwerere Fälle (höherer `severity_modifier`) bleiben entsprechend länger.

#### `base_mortality_rate`
- **Typ:** float, [0, 1]
- **Einfluss:** Tägliche Sterbewahrscheinlichkeit für alle Patienten. Wird für Carrier mit `severity_modifier` multipliziert.
- **Formel:** `mortality_rate = base_mortality_rate × severity_modifier` (nur für Carrier; für Susceptible = `base_mortality_rate`)
- **Basis:** Abgeleitet aus Schweizer Spitalmortalität ~2.5% pro Aufenthalt ÷ 5.5-Tage-LoS (BFS Medizinische Statistik 2022; Huang & Platt 2003). Standard: `0.0045` → ~0.5% Tageswahrscheinlichkeit.
- **Konsequenz:** Ein Patient wird via `discharge()` aus dem Spital entfernt, bevor die Entlassungslogistik greift. Carrier-Verstorbene haben keinen Entlassungs-Plot-Eintrag.

#### `discharge_logistic_k`
- **Typ:** float
- **Einfluss:** Steilheit der logistischen Entlassungskurve. Höherer Wert → abrupterer Übergang.
- **Formel:** `p_discharge = 1 / (1 + exp(-k × (days_over - t_half)))`

#### `discharge_logistic_t_half`
- **Typ:** float (Tage)
- **Einfluss:** Anzahl Tage nach dem Zieldatum, an denen die Entlassungswahrscheinlichkeit 50% beträgt. Bestimmt die typische Verzögerung.
- **Beispiel:** `3.0` → Im Durchschnitt 3 Tage nach dem Zieldatum werden Susceptible entlassen.

---

### Aufnahmen

#### `daily_admission_rate`
- **Typ:** float (Patienten/Tag)
- **Einfluss:** Poisson-Erwartungswert für neue Patienten pro Tag über alle Spitäler. `0` deaktiviert Aufnahmen.
- **Formel:** `n_new ~ Poisson(daily_admission_rate)`. Neue Patienten werden dem Spital mit der grössten freien Kapazität zugewiesen (gewichtet).

#### `community_carrier_fraction`
- **Typ:** float, [0, 1]
- **Einfluss:** Anteil der neu aufgenommenen Patienten, die bereits Carrier sind (Einschleppung aus der Gemeinschaft).
- **Verbindung:** Bestimmt den Anfangszustand neuer Patienten — `HealthState.CARRIER` mit vordefinierten Resistenzwerten.

#### `replacement_resistant_fraction`
- **Typ:** float, [0, 1]
- **Einfluss:** Resistenzgrad der aus der Gemeinschaft eingeschleppten Carrier. Initialisiert `patient.resistant_fraction` bei Aufnahme.
- **Verbindung zum Mikro-Layer:** Dieser Wert fliesst beim ersten Mikro-Request als `"initial_state": {"resistant_fraction": ...}` ein.

#### `replacement_dominant_genotype`
- **Typ:** str
- **Einfluss:** Initialer Genotyp für eingeschleppte Carrier (typisch `"S"` = sensitiver Wildtyp).
- **Verbindung:** Wird auf `patient.dominant_genotype` gesetzt und beim ersten Mikro-Request übergeben.

#### `max_occupancy_per_hospital`
- **Typ:** int
- **Einfluss:** Harte Kapazitätsgrenze pro Spital. Wenn alle Spitäler voll sind, werden keine neuen Patienten mehr aufgenommen.

---

### Verlegungen

#### `daily_transfer_rate`
- **Typ:** float, [0, 1]
- **Einfluss:** Tägliche Wahrscheinlichkeit pro Patient, in ein anderes Spital verlegt zu werden. `0` deaktiviert Verlegungen.
- **Zielwahl:** Gewichtet nach `(freie_kapazität / distanz)` — nahe und weniger belegte Spitäler werden bevorzugt.
- **Hinweis:** Der Patient behält seine Abteilung (Ward/ICU), seinen Zustand und seinen Mikro-Episode-Status.

---

### Grid-Parameter

Das Spital wird als zweidimensionales Grid modelliert. Patienten erhalten eine Position im Grid, welche die räumliche Nähe für die Übertragung bestimmt.

#### `dept_grid_cols`
- **Typ:** int
- **Einfluss:** Anzahl Spalten im Abteilungs-Grid. Bestimmt die Gesamtzahl möglicher Positionen: `cols × rows`.

#### `dept_grid_rows`
- **Typ:** int
- **Einfluss:** Gesamtzahl Reihen im Grid (Ward + ICU zusammen).

#### `dept_grid_icu_rows`
- **Typ:** int
- **Einfluss:** Anzahl Reihen, die als ICU gelten. Reihen 0 bis `icu_rows - 1` sind ICU, der Rest Ward.
- **Formel ICU-Anteil:** `icu_fraction = icu_rows / rows`. Bestimmt die Aufnahmewahrscheinlichkeit für ICU vs. Ward.

#### `network_grid_cols`
- **Typ:** int
- **Einfluss:** Anzahl Spalten im Spitalnetzwerk-Grid. Bestimmt die geometrischen Distanzen zwischen Spitälern für die Verlegungsgewichtung.

#### `proximity_decay_alpha`
- **Typ:** float
- **Einfluss:** Räumlicher Zerfallskoeffizient für die Übertragungsgewichtung. Bestimmt, wie stark die Übertragungswahrscheinlichkeit mit der Distanz abnimmt.
- **Formel:** `proximity_weight = exp(-alpha × chebyshev_distance)`
- **Beispiel:** `alpha = 0.5` → Nachbar (Distanz 1): Gewicht `≈ 0.61`. Distanz 2: `≈ 0.37`. `alpha = 2.0` → Nachbar: `≈ 0.14` (viel stärkerer Abfall).

---

## Verbindung zu `patient.py`

Die `Patient`-Klasse ist das zentrale Interface zwischen Makro- und Mikro-Ebene.

### Makro → Patient (täglich via `PatientDailyContext`)

```
MacroSimulator._build_context()
    └── PatientDailyContext(
            hospital_id, department,
            hygiene_level,           ← aus base_hygiene
            isolation_effectiveness, ← aus base_isolation_effectiveness
            diagnostic_speed,        ← aus base_diagnostic_speed
            is_isolated,             ← stochastisch (carrier_isolation_probability × diagnostic_speed)
            regimen                  ← stochastisch (icu/ward_abx_probability)
        )
        └── patient.update_context(ctx)
                ├── patient.is_isolated     = ctx.is_isolated
                ├── patient.regimen         = ctx.regimen
                ├── patient.department      = ctx.department
                └── patient.adherence       ← abgeleitet aus patient.compliance + department
```

### Patient → Makro (Multiplikatoren)

Der Patient gibt der Makro-Übertragungsformel zwei Grössen zurück:

| Methode | Bestandteile | Einfluss |
|---|---|---|
| `transmission_multiplier_for_macro()` | `sociability × relative_transmissibility` | Skaliert `base_transmission_rate` nach oben/unten |
| `susceptibility_multiplier_for_macro()` | `vulnerability × (1 / immune_strength) × prior_infection_flag` | Skaliert den Hazard nach oben/unten |

`relative_transmissibility` wird täglich vom Mikro-Layer aktualisiert — resistentere, besser angepasste Stämme können eine höhere Transmissibilität haben.

---

## Verbindung zum Mikro-Layer

### Makro → Mikro (täglicher Request)

Für jeden aktiven Carrier baut `patient.make_micro_request()` einen Dictionary-Request:

```python
{
    "episode_id": "...",
    "t_day": 42,
    "setting": "icu",           # department
    "abx": {
        "on": True,
        "class": "beta_lactam",
        "dose_level": "std"
    },
    "adherence": 0.85,
    "host": {
        "age_years": 65,
        "immune_strength": 0.75,
        "immune_status": "normal",
        "vulnerability": 1.0,
        "history_flags": ["prior_abx"]
    },
    "initial_state": {
        "resistant_fraction": 0.45,
        "dominant_genotype": "R2",
        "dominant_strain_name": ""
    }
}
```

Die Mikro-Ebene simuliert daraus die bakterielle Evolution für diesen Tag (12 Schritte).

### Mikro → Makro (tägliche Response)

Die Antwort des Mikro-Layers wird via `patient.apply_micro_response()` auf den Patienten geschrieben:

```python
{
    "updated_state": {
        "resistant_fraction": 0.72,     # → patient.resistant_fraction
        "dominant_genotype": "R3",      # → patient.dominant_genotype
        "dominant_strain_name": "..."
    },
    "derived_effects": {
        "relative_transmissibility": 1.3,  # → patient.relative_transmissibility
        "p_clearance": 0.005,              # → patient.p_clearance
        "severity_modifier": 1.4           # → patient.severity_modifier
    }
}
```

**Wichtig:** Nur die Makro-Ebene entscheidet über den S/C-Zustandswechsel. Der Mikro-Layer liefert lediglich `p_clearance` als Wahrscheinlichkeit zurück — die eigentliche C→S Entscheidung trifft `MacroSimulator` via `patient.should_clear_today()`.

---

## Zustandsübergänge auf Makro-Ebene

```
                  carrier_isolation_probability × diagnostic_speed
                  ┌─────────────────────────────────────────────────┐
                  │                                                   ↓
    S ──────────────────────────────────────────────────────────── C (isolated)
         p_colonize                                                    │
         (Hazard-Formel)                                               │ carrier_extension_days
                                                                       │ verlängert Aufenthalt
    S ←──────────────────────────────────────────────────────────── C
              p_clearance (vom Mikro-Layer, täglich aktualisiert)
```

| Übergang | Auslöser | Verantwortlich |
|---|---|---|
| S → C | `p_colonize` aus Hazard-Formel | `MacroSimulator._do_transmission()` |
| C → S | `patient.p_clearance` (vom Mikro) | `MacroSimulator.step()` → `patient.should_clear_today()` |
| C → C (isolated) | `effective_detection_prob` | `MacroSimulator._build_context()` |

---

## Zusammenfassung: Parameterübersicht

| Parameter | Gruppe | Direkte Wirkung |
|---|---|---|
| `base_hygiene` | Infektionskontrolle | Reduziert Hazard für alle Susceptible |
| `base_isolation_effectiveness` | Infektionskontrolle | Reduziert Beitrag isolierter Carrier + Hazard isolierter Susceptible |
| `base_diagnostic_speed` | Infektionskontrolle | Skaliert Erkennungsrate von Carriern |
| `base_transmission_rate` | Übertragung | Basis-Beta pro Kontakt und Carrier |
| `daily_contact_attempts` | Übertragung | Skaliert den gesamten Hazard |
| `icu_abx_probability` | ABX-Politik | Wahrscheinlichkeit für ABX in ICU → Mikro-Selektionsdruck |
| `ward_abx_probability` | ABX-Politik | Wahrscheinlichkeit für ABX auf Ward → Mikro-Selektionsdruck |
| `carrier_isolation_probability` | Isolation | Basisrate für tägliche Carrier-Erkennung |
| `transmission_mutation_probability` | Übertragungsvererbung | Chance auf Stamm-Drift bei Kolonisierung |
| `transmission_resistance_mutation_std` | Übertragungsvererbung | Grösse des Resistenz-Drifts bei Kolonisierung |
| `los_mean_ward` | Verweildauer | Mittlere Ward-Verweildauer (Log-Normal) |
| `los_mean_icu` | Verweildauer | Mittlere ICU-Verweildauer (Log-Normal) |
| `los_sigma` | Verweildauer | Streuung der Verweildauerverteilung |
| `carrier_extension_days` | Entlassung | Verlängerung des Aufenthalts bei Carrier-Erkennung (× `severity_modifier` im Discharge-Loop) |
| `base_mortality_rate` | Entlassung | Tägliche Basissterblichkeit; für Carrier × `severity_modifier` |
| `discharge_logistic_k` | Entlassung | Steilheit der Entlassungskurve |
| `discharge_logistic_t_half` | Entlassung | Verzögerung nach Zieldatum bis 50% Entlassungswahrscheinlichkeit |
| `daily_admission_rate` | Aufnahmen | Neue Patienten/Tag (Poisson) |
| `community_carrier_fraction` | Aufnahmen | Anteil eingeschleppter Carrier bei Aufnahme |
| `replacement_resistant_fraction` | Aufnahmen | Resistenzgrad eingeschleppter Carrier |
| `replacement_dominant_genotype` | Aufnahmen | Genotyp eingeschleppter Carrier |
| `max_occupancy_per_hospital` | Kapazität | Maximale Patientenzahl pro Spital |
| `daily_transfer_rate` | Verlegungen | Tägliche Verlegungswahrscheinlichkeit pro Patient |
| `dept_grid_cols` | Grid | Spalten im Abteilungs-Grid |
| `dept_grid_rows` | Grid | Reihen im Abteilungs-Grid (gesamt) |
| `dept_grid_icu_rows` | Grid | Davon ICU-Reihen (bestimmt ICU-Anteil) |
| `network_grid_cols` | Grid | Spalten im Spitalnetzwerk-Grid (Distanzen) |
| `proximity_decay_alpha` | Grid | Räumlicher Zerfallskoeffizient der Übertragungsgewichtung |
