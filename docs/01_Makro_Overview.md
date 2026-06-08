# Makro-Ebene (Spital-Netzwerk)

Die Makro-Ebene simuliert die **Ausbreitung von MRSA auf Spitalebene**: Aufnahme,
Entlassung und Verlegung von Patienten, stochastische Übertragung zwischen Patienten
(S → C), spontane Erholung von Trägern (C → S), Antibiotikaverordnung, Isolation und den
täglichen Austausch mit der Mikro-Ebene. Der Patientenzustand ist hier binär:
**SUSCEPTIBLE (S)** oder **CARRIER (C)**.

Code: [`config.py`](../src/mss/simulation/macro/config.py) (alle Parameter),
[`simulator.py`](../src/mss/simulation/macro/simulator.py) (Tagesablauf, Formeln),
[`grid.py`](../src/mss/simulation/macro/grid.py) (Räume und Netzwerk),
[`agents.py`](../src/mss/simulation/macro/agents.py) (Mesa-Wrapper).

> Werte, Typen und Quellen aller Parameter stehen in der Referenz
> [`config/01_Makro_Parameterübersicht.md`](../config/01_Makro_Parameterübersicht.md) (inklusive
> Kalibrierungen und der Tabelle zum gezielten Ein-/Ausschalten von Mechanismen).

## Tagesablauf (`MacroSimulator.step`)

Jeden Tag laufen die folgenden Schritte **in dieser Reihenfolge**:

```
1. Pro Spital:
   a. Spontane Erholung (C → S)
   b. Kontext-Update + Isolationserkennung
   c. Mikro-Requests sammeln
2. Entlassung (inkl. Mortalität)       pro Spital
3. Aufnahme neuer Patienten            global (Poisson)
4. Mikro-Batch anwenden                alle Carrier global
5. Übertragung (S → C)                 pro Spital
6. Verlegungen zwischen Spitälern      global
```

Die Reihenfolge ist bewusst gewählt: Wer heute gesund wird (1a), ist bei der Entlassung
(2) bereits als Susceptible sichtbar. Der detaillierte Phasen-Ablauf inklusive Mikro-
Kopplung steht in [`docs/system_overview/Flowchart_v1.mmd`](system_overview/Flowchart_v1.mmd).

## Mechanismen im Überblick

Jeder Mechanismus ist über Config-Parameter steuerbar (Werte/Typen/Quellen in
`config/01_Makro_Parameterübersicht.md`):

- **Infektionskontrolle**: `base_hygiene` senkt den Hazard für alle;
  `base_isolation_effectiveness` reduziert Senden und Empfangen isolierter Patienten;
  `base_diagnostic_speed` skaliert die Erkennungsrate.
- **Übertragung**: `base_transmission_rate` (β₀, kalibriert) pro Kontakt und Carrier,
  `daily_contact_attempts` skaliert den Gesamt-Hazard, `proximity_decay_alpha` gewichtet
  räumliche Nähe.
- **Antibiotika**: `icu_abx_probability` / `ward_abx_probability` bestimmen die tägliche
  Verordnungschance; Klasse und Dosis werden gleichverteilt gezogen und gehen als
  Selektionsdruck in die Mikro-Ebene.
- **Erkennung/Isolation**: `carrier_isolation_probability` × `base_diagnostic_speed`.
- **Verweildauer/Entlassung**: `los_mean_ward` / `los_mean_icu` / `los_sigma`
  (Log-Normal); `carrier_extension_days` verlängert den Aufenthalt erkannter Carrier; die
  logistische Entlassung steuern `discharge_logistic_k` / `discharge_logistic_t_half`.
- **Mortalität**: `base_mortality_rate` (für Carrier × `severity_modifier`), täglich
  ausgewertet, von Aufnahmen entkoppelt.
- **Patientenfluss**: `daily_admission_rate` (Poisson), `community_carrier_fraction` und
  `replacement_*` für eingeschleppte Träger, `max_occupancy_per_hospital`,
  `daily_transfer_rate` für Verlegungen.

### Übertragungsformel

```
weighted_force = Σ [ base_transmission_rate × transmissibility × proximity_weight × (1 - iso_eff)^(carrier_isolated) ]

hazard = daily_contact_attempts × (1 - hygiene) × (weighted_force / n_total)
         × susceptibility × (1 - iso_eff)^(susceptible_isolated)

p_colonize = 1 - exp(-max(0, hazard))
```

`proximity_weight = exp(-proximity_decay_alpha × chebyshev_distance)` gewichtet die Nähe
im Abteilungsgitter (eine Zelle ist ein Mehrbettzimmer). Bei Übertragung erbt der neue
Carrier den Stamm der Quelle exakt; die Resistenz-Drift übernimmt ab dem Folgetag die
Mikro-Ebene.

## Verbindung zu `patient.py`

Die `Patient`-Klasse ist die Schnittstelle zwischen Makro und Mikro.

**Makro → Patient** (täglich via `PatientDailyContext` in `_build_context()`):
Krankenhaus/Station, `hygiene_level`, `isolation_effectiveness`, `diagnostic_speed`,
`is_isolated` (stochastisch aus `carrier_isolation_probability × diagnostic_speed`) und
`regimen` (stochastisch aus den ABX-Wahrscheinlichkeiten). `patient.update_context()`
übernimmt das und leitet `adherence` aus `compliance` ab (ICU +0.1).

**Patient → Makro** (Multiplikatoren in die Übertragungsformel):

| Methode | Bestandteile | Wirkung |
|---|---|---|
| `transmission_multiplier_for_macro()` | `sociability × relative_transmissibility` | skaliert `base_transmission_rate` |
| `susceptibility_multiplier_for_macro()` | `1 / immune_strength` | skaliert den Hazard |

`relative_transmissibility` wird täglich von der Mikro-Ebene aktualisiert.

## Verbindung zur Mikro-Ebene

Für jeden aktiven Carrier baut `patient.make_micro_request()` einen Request (ABX-Regime,
`adherence`, Immunstärke, aktueller Resistenzzustand). Die Antwort schreibt
`patient.apply_micro_response()` zurück: `resistant_fraction`, `dominant_genotype`,
`relative_transmissibility`, `p_clearance` und `severity_modifier`. Den vollständigen
Within-Host-Ablauf beschreibt [`docs/02_Mikro_Overview.md`](02_Mikro_Overview.md).

**Wichtig:** Über den S/C-Zustandswechsel entscheidet allein die Makro-Ebene. Die Mikro
liefert nur `p_clearance`; die C→S-Entscheidung trifft `patient.should_clear_today()`.

## Zustandsübergänge


| Übergang | Auslöser | Verantwortlich |
|---|---|---|
| S → C | `p_colonize` aus der Hazard-Formel | `MacroSimulator._do_transmission()` |
| C → S | `patient.p_clearance` (vom Mikro) | `MacroSimulator.step()` → `should_clear_today()` |
| C → C (isoliert) | `effective_detection_prob` | `MacroSimulator._build_context()` |
