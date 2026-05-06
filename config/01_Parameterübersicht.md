# Parameterübersicht

Jeder Parameter ist einem der folgenden Typen zugeordnet:

- **geschätzt** – aus Literatur oder plausiblen Annahmen abgeleitet
- **kalibriert** – durch Simulation so bestimmt, dass eine beobachtete Zielgrösse reproduziert wird
- **nicht identifizierbar** – Sweep durchgeführt, Modell reagiert im plausiblen Bereich nicht signifikant; Wert physikalisch oder neutral motiviert
- **offen** – noch nicht bestimmt (Phase 3)

---

## Kalibrierungen

### Kalibrierung 1 — `base_transmission_rate` (β₀)

**Script:** `src/mss/cli/run_single_ward_calibration.py`
**Methode:** Analytische Rückrechnung aus Literatur-Zielwert, validiert mit stochastischer Einzelgitter-Simulation (Mikro deaktiviert)

Modellgleichung:

$$
\lambda = \beta_0 \cdot c \cdot (1 - H) \cdot \frac{I}{N}
$$

**Zielwert:** MRSA-Akquisitionsrate auf Station = 4.6–5.4 pro 1000 Patienten-Tage, d.h. λ ≈ 0.0046–0.0054 pro Tag — [PMC 5384532](https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/)

Gegebene Werte: c = 14, H = 0.75, I/N ≈ 0.017

$$
\beta_0 = \frac{\lambda}{c \cdot (1 - H) \cdot (I/N)} \approx 0.065 \text{ bis } 0.075
$$

**Ergebnis: β₀ = 0.07**

---

### Kalibrierung 2 — Phase 1: Makro-Transmissionsparameter

**Script:** `mss-sweep --sweep config/calibration/phase1_*.yml`
**Methode:** Stochastischer Parameter-Sweep (n\_seeds = 5, run\_days = 90, `use_micro = false`)
`use_micro = false` isoliert die Makro-Dynamik von unkalibriertem Mikro-Verhalten: Patienten behalten ihre Template-Standardwerte (`p_clearance`, `relative_transmissibility`, `severity_modifier`) unverändert. Methodisch analog zur analytischen β₀-Kalibrierung — stochastisch statt geschlossen.

#### `base_isolation_effectiveness`

| | |
|---|---|
| Sweep-Bereich | 0.40 – 1.00 (zwei Runden: grob + fein) |
| Zielgrösse | `acquisition_rate_per_1000 = 1.5` |
| Ziel-Begründung | Schweizer Niedrig-MRSA-Setting; plausibler Korridor 0.5–2.0/1000 Pat.-Tage (Swissnoso Annual Report 2022, ECDC EARS-Net 2023) |
| Sweep-Befund | Monotoner Abfall mit steigender Isolation. Im Bereich 0.95–1.0 erreichen alle Werte die Zielgrösse innerhalb des Simulationsrauschens (±0.3). |
| Mechanisch bester Wert | 0.9857 → acquisition\_rate = 1.5302 (Abstand 0.030) |
| **Gewählter Wert** | **0.98** — runder Wert im identifizierten Bereich; Unterschiede zwischen 0.97–1.0 sind Rauschen, nicht Signal |

#### `proximity_decay_alpha`

| | |
|---|---|
| Sweep-Bereich | 0.1 – 3.0 (zwei Runden) |
| Zielgrösse | `mean_prevalence = 0.02` |
| Ziel-Begründung | 1.5–4 % MRSA-Prävalenz in CH-Spitälern (Swissnoso 2022, Rohner et al. 2017) |
| Sweep-Befund | **Zielwert nicht erreichbar.** Community-Aufnahmen (`community_carrier_fraction × daily_admission_rate` ≈ 2.1/Tag) erzeugen einen strukturellen Prävalenz-Floor von ~11 %, unabhängig vom räumlichen Abklingparameter. Kurve flach ab α ≈ 1.5 (mean\_prevalence 0.113–0.124, reines Rauschen). |
| Schlussfolgerung | Parameter ist über diese Zielgrösse nicht identifizierbar. Wahl erfolgt physikalisch: `w = exp(−α·d)` → bei α = 1.5: Nachbar (d=1): 0.22, zwei Zellen (d=2): 0.05. Patienten interagieren primär mit unmittelbaren Nachbarn — realistisch für eine Spitalstation. |
| **Gewählter Wert** | **1.5** — Untergrenze des nicht identifizierbaren Bereichs; physikalisch plausibel |

---

### Kalibrierung 3 — Phase 2: Populations-Template-Parameter

**Script:** `mss-sweep --sweep config/calibration/phase2_*.yml`
**Methode:** Identisch zu Kalibrierung 2 (`use_micro = false`, n\_seeds = 5, run\_days = 90)
**Zielgrösse für alle Phase-2-Sweeps:** `acquisition_rate_per_1000 = 1.5` (Swissnoso 2022 / ECDC EARS-Net 2023)

**Übergreifender Befund:** Alle getesteten Template-Parameter zeigen im plausiblen Bereich nur schwachen, verrauschten Einfluss auf `acquisition_rate_per_1000`. Der Community-Aufnahmen-Floor dominiert die Prävalenz; die Akquisitionsrate variiert im Bereich 1.5–1.8 ohne klares Optimum. Neutralwerte (1.0) sind in allen Fällen vertretbar.

| Parameter | Sweep-Bereich | Mechanisch bester Wert | Befund | Gewählter Wert |
|---|---|---|---|---|
| `susceptible_template.vulnerability` | 0.4 – 2.5 | 2.2 → 1.506 | Nicht-monotone Kurve (Rauschen); kein Signal identifizierbar | **1.0** (Neutralwert) |
| `susceptible_template.immune_strength` | 0.3 – 2.0 | 1.03 → 1.573 | Schwacher Effekt; 1.0 liegt im Bereich des Minimums | **1.0** (gerundet von 1.03) |
| `susceptible_template.sociability` | 0.4 – 2.0 | 1.8 → 1.563 | Flache Kurve über gesamten Bereich (alle Werte ±0.15 um Ziel) | **1.0** (Neutralwert) |
| `carrier_template.sociability` | 0.4 – 2.0 | 1.4 → 1.534 | Insensitiv im Bereich 1.0–1.6 (Floor bei <1.0: identische Werte) | **1.0** (Neutralwert) |
| `carrier_template.immune_strength` | — | — | Kein Makro-Effekt (nur Mikro); Phase-3-Parameter | **0.75** (Schätzwert) |
| `carrier_template.vulnerability` | — | — | **Entfernt** — kein Makro- und kein Mikro-Effekt für Carrier | — |

---

## Populations-Parameter

`initial_department` ist bei Simulationsstart für alle Patienten auf `ward` gesetzt.

### Populationsgrösse

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `hospitals` | geschätzt | 6 | Anzahl Spitäler im simulierten Netzwerk | — |
| `carrier_count` | geschätzt | 20 (≈ 2 % von N) | Anzahl initial kolonisierter Patienten bei Simulationsstart | [PMC 5384532](https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/) |
| `susceptible_count` | geschätzt | 980 | Anzahl empfänglicher Patienten bei Simulationsstart | abgeleitet |

### Agenten-Initialwerte: Susceptible

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `compliance` | geschätzt | 0.80 | Durchschnittliche Händehygiene-Compliance im Akutspital | [Swissnoso Report 2023](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/8_Swissnoso_Publikationen/Swissnoso_Annual_Report_HAI_CH_2023_FULL_REPORT_EN_DOUBLE.pdf) |
| `vulnerability` | nicht identifizierbar | 1.0 | Suszeptibilitäts-Multiplikator (1.0 = Referenz). Sweep nicht-monoton, kein Signal; Neutralwert gewählt. | Kalibrierung 3 |
| `immune_strength` | nicht identifizierbar | 1.0 | Immunkraft-Multiplikator (1.0 = Referenz). Sweep-Optimum bei 1.03, gerundet. | Kalibrierung 3 |
| `immune_status` | geschätzt | normal | Immunstatus (normal / compromised) | Modellannahme |
| `sociability` | nicht identifizierbar | 1.0 | Kontaktintensität relativ zum Durchschnitt. Sweep zeigt flache Kurve; Neutralwert gewählt. | Kalibrierung 3 |

### Agenten-Initialwerte: Carrier

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `compliance` | geschätzt | 0.70 | Leicht geringere Compliance als bei Susceptible | Modellannahme |
| `immune_strength` | geschätzt | 0.75 | Leicht reduzierte Immunabwehr für persistente Carrier (kein Makro-Effekt; wirkt im Mikro auf Immun-Clearance). Phase-3-Kalibrierung ausstehend. | Template-Default; Kalibrierung 3 |
| `immune_status` | geschätzt | normal | Immunstatus (normal / compromised) | Modellannahme |
| `sociability` | nicht identifizierbar | 1.0 | Kontaktintensität. Sweep nicht identifizierbar im Bereich 1.0–1.6; Neutralwert gewählt. | Kalibrierung 3 |
| `resistant_fraction` | geschätzt | 0.073 | Anteil MRSA an S. aureus Isolaten (CH-Durchschnitt) | [ANRESIS Report 2023](https://www.anresis.ch/wp-content/uploads/2023/04/CAESAR_report_2023.pdf) |
| `dominant_genotype` | geschätzt | R2 | Dominanter Stamm bei Simulationsstart | Schweizer/Europäischer MRSA-Standardstamm, als Modellannahme |
| `relative_transmissibility` | geschätzt | 1.0 | Initialer Übertragbarkeits-Multiplikator (Referenzwert). β₀ wurde für HA-MRSA-Carrier kalibriert → Multiplikator = 1.0 ist definitional. Nach Tag 1 durch Mikro überschrieben. | [Grundmann et al. 2006, Lancet](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(06)68849-2/abstract) |
| `p_clearance` | geschätzt | 0.0039 | Tägliche Wahrscheinlichkeit der spontanen Dekolonisierung ohne Behandlung (≈ 1/255 Tage mittlere Tragezeit). | [Clin. Infect. Dis. 32(10)](https://academic.oup.com/cid/article/32/10/1393/465089) |
| `severity_modifier` | geschätzt | 1.2 | Skaliert `carrier_extension_days` (längerer Aufenthalt) und `base_mortality_rate` (erhöhtes Sterblichkeitsrisiko für Carrier). Wird täglich durch Mikro aktualisiert. | [Cosgrove et al. 2003, CID](https://academic.oup.com/cid/article/36/1/53/298827); [Huang & Platt 2003, Ann Intern Med](https://www.acpjournals.org/doi/10.7326/0003-4819-139-5_Part_1-200309020-00008) |

---

## Makro-Parameter

### Netzwerk & Gitter

Jedes Spital wird als Gitter aus `dept_grid_cols × dept_grid_rows` Zellen abgebildet. Die untersten `dept_grid_icu_rows` Reihen sind ICU. `network_grid_cols` bestimmt die räumliche Anordnung der Spitäler im Netzwerk.

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `dept_grid_cols` | geschätzt | 5 | Zimmer-Anordnung (5 × 4 = 20 Betten/Zimmer pro Ward) | — |
| `dept_grid_rows` | geschätzt | 4 | Gesamtreihen im Gitter | — |
| `dept_grid_icu_rows` | geschätzt | 1 | Ca. 25 % ICU-Betten-Anteil | — |
| `network_grid_cols` | geschätzt | 3 | 3 × 2-Cluster für 6 Spitäler | — |
| `max_occupancy_per_hospital` | geschätzt | 200 | Hard Cap für Aufnahmen und Transfers pro Spital | — |
| `proximity_decay_alpha` | nicht identifizierbar | 1.5 | Abklingkoeffizient α: w = exp(−α·d). Sweep zeigte strukturellen Prävalenz-Floor; Parameter physikalisch motiviert (Nachbar: 0.22, d=2: 0.05). | Kalibrierung 2 |

### Transmission

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `daily_contact_attempts` | geschätzt | 14 | Relevante Kontaktepisoden pro Patient und Tag | [PMC 9152759](https://pmc.ncbi.nlm.nih.gov/articles/PMC9152759/) |
| `base_hygiene` | geschätzt | 0.75 | Reduktion der Transmissionswahrscheinlichkeit durch Hygienemassnahmen | — |
| `base_transmission_rate` | kalibriert | 0.07 | Übertragungswahrscheinlichkeit pro relevantem Kontakt | Kalibrierung 1 |
| `carrier_isolation_probability` | geschätzt | 0.65 | Wahrscheinlichkeit, dass ein erkannter Carrier täglich isoliert wird | [HIS/IPS MRSA Guidelines](https://www.his.org.uk/media/djom3t4f/joint-healthcare-infection-society-his-and-infection-prevention-society-ips-guidelines-for-the-prevention-and-control-of-meticillin-resistant-staphylococcus-aureus-mrsa-in-healthcare-facilities.pdf) |
| `base_isolation_effectiveness` | kalibriert | 0.98 | Reduktion der Transmissionswahrscheinlichkeit durch Isolationsmassnahmen. Sweep-Optimum 0.9857; auf 0.98 gerundet (Unterschiede im Bereich 0.95–1.0 sind Simulationsrauschen). | Kalibrierung 2 |
| `base_diagnostic_speed` | geschätzt | 1.0 | Geschwindigkeit der MRSA-Diagnose (Skalierungsfaktor; 1.0 = Standard) | Vorerst belassen |

### Aufenthalt & Entlassung

Die Entlassung folgt einer logistischen Kurve nach dem geplanten Entlassungstag.

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `los_mean_ward` | geschätzt | 5.0 Tage | Mittlere Aufenthaltsdauer auf der Normalstation | [BFS Spitalstatistik 2023](https://www.bfs.admin.ch/bfs/de/home.assetdetail.34027817.html) |
| `los_mean_icu` | geschätzt | 2.8 Tage | Mittlere Aufenthaltsdauer auf der Intensivstation | [PubMed 10890670](https://pubmed.ncbi.nlm.nih.gov/10890670/) |
| `los_sigma` | geschätzt | 0.6 | Log-Normal-Parameter für LOS-Streuung | [FOPH Report](https://www.bag.admin.ch/dam/en/sd-web/3aJgqPNhp9yQ/sentinel-corona-bericht-maerz-2023.pdf) |
| `carrier_extension_days` | geschätzt | 12 Tage | Verlängerung bei erkannter MRSA-Kolonisierung durch Kontaktisolation; wird mit `severity_modifier` skaliert | [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S187603412030438X) |
| `discharge_logistic_k` | geschätzt | 1.0 | Steilheit der logistischen Entlassungskurve | Logistisches Approximationsmodell |
| `discharge_logistic_t_half` | geschätzt | 3.0 Tage | Halbwert der logistischen Entlassungskurve nach geplantem Entlassungstag | Logistisches Approximationsmodell |

### Patientenfluss

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `daily_admission_rate` | geschätzt | 126 | Tägliche Neuaufnahmen ins gesamte Netzwerk (zielt auf ~105 Patienten/Spital bei LOS = 5) | Abgeleitet aus Bettenbelegung und LOS |
| `daily_transfer_rate` | geschätzt | 0.004 | Tägliche Verlegungswahrscheinlichkeit zwischen Spitälern | [PMC 3727950](https://pmc.ncbi.nlm.nih.gov/articles/PMC3727950/) |
| `community_carrier_fraction` | geschätzt | 0.017 | Anteil MRSA-positiver Patienten bei Spitaleintritt (Community-acquired) | [PubMed 27658666](https://pubmed.ncbi.nlm.nih.gov/27658666/) |
| `replacement_resistant_fraction` | geschätzt | 0.2 | Anteil resistenter Stämme bei neu eintretenden Community-Carriern | [ECDC EARS-Net](https://www.ecdc.europa.eu/en/about-us/networks/disease-networks-and-laboratory-networks/ears-net-data) |
| `replacement_dominant_genotype` | geschätzt | S | Dominanter Genotyp bei neu eintretenden Community-Carriern | — |
| `base_mortality_rate` | geschätzt | 0.0045 | Tägliche Sterblichkeit für alle Patienten. Für Carrier skaliert mit `severity_modifier` (→ 0.0054/Tag). Herleitung: CH In-Hospital-Sterblichkeit ~2.5 % ÷ mittlere LOS 5.5 Tage = 0.0045/Tag. | [BFS Medizinische Statistik 2022](https://www.bfs.admin.ch/bfs/de/home.assetdetail.34027817.html); [Huang & Platt 2003, Ann Intern Med](https://www.acpjournals.org/doi/10.7326/0003-4819-139-5_Part_1-200309020-00008) |

### Antibiotika

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `icu_abx_probability` | geschätzt | 0.62 | Tägliche Wahrscheinlichkeit, dass ein ICU-Patient Antibiotika erhält | [SwissNoso PPS 2017](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/2_Punktpraevalenzstudie/Report_Point_Prevalence_Survey_2017_of_HAI_and_antimicrobial_use_in_Swiss_acute_care_hospitals.pdf) |
| `ward_abx_probability` | geschätzt | 0.326 | Tägliche Wahrscheinlichkeit, dass ein Ward-Patient Antibiotika erhält | [SwissNoso Annual Report 2023](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/8_Swissnoso_Publikationen/Swissnoso_Annual_Report_HAI_CH_2023_FULL_REPORT_EN_SINGLE.pdf) |

### Mutation & Resistenz (Makro-Ebene)

Diese Parameter steuern, ob und wie stark bei einer Transmission eine Resistenzveränderung auf Makro-Ebene ausgelöst wird (unabhängig von der Mikro-Simulation). Phase-3-Kalibrierung ausstehend.

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `transmission_mutation_probability` | offen | 0.25 | Wahrscheinlichkeit einer Resistenzmutation bei Transmission | Phase 3 |
| `transmission_resistance_mutation_std` | offen | 0.08 | Standardabweichung der `resistant_fraction`-Änderung bei Transmissionsmutation | Phase 3 |

---

## Mikro-Parameter

Die Mikro-Simulation modelliert die bakterielle Evolution innerhalb einzelner Patienten. Ein Makro-Tag wird in `steps_per_day` feinere Zeitschritte unterteilt. Bei Phase-1/2-Kalibrierungen ist der Mikro-Simulator deaktiviert (`use_micro = false`); Patienten behalten ihre Template-Standardwerte. Alle Mikro-Parameter sind Phase-3-Kalibrierung zugeordnet.

### Simulationsstruktur

| Parameter | Typ | Wert | Beschreibung |
|---|---|---|---|
| `steps_per_day` | geschätzt | 12 | Mikro-Zeitschritte pro Simulationstag |
| `founder_pool_size` | geschätzt | 32 | Anzahl vordefinierter Gründerstämme |
| `founder_pool_gene_noise_std` | geschätzt | 0.02 | Rauschen bei Geninitialisierung der Gründerstämme |
| `gene_presence_threshold` | geschätzt | 0.2 | Minimale Genexpressionsstärke für "Gen vorhanden" |
| `max_strains` | geschätzt | 40 | Maximale gleichzeitig aktive Stämme pro Patient |

### Bakterielle Populationsdynamik

| Parameter | Typ | Wert | Beschreibung |
|---|---|---|---|
| `carrying_capacity` | offen | 5 × 10⁸ | Maximale Bakterienpopulation pro Patient |
| `min_population` | offen | 100 | Minimale Populationsgrösse; darunter gilt Stamm als ausgestorben |
| `clearance_threshold` | offen | 1 000 | Schwellenwert für immunologische Klärung |
| `growth_rate_per_step` | offen | 0.18 | Wachstumsrate pro Mikro-Zeitschritt |
| `death_rate_per_step` | offen | 0.06 | Sterberate pro Mikro-Zeitschritt |
| `strain_prune_threshold` | offen | 200 | Pruning-Untergrenze; Stämme darunter werden entfernt |

### Mutation & Horizontaler Gentransfer

| Parameter | Typ | Wert | Beschreibung |
|---|---|---|---|
| `base_mutation_rate` | offen | 0.012 | Basis-Punktmutationsrate pro Mikro-Zeitschritt |
| `mutation_std` | offen | 0.025 | Standardabweichung der Mutationsrate |
| `stress_mutation_boost` | offen | 40 | Multiplikator auf Mutationsrate unter Antibiotikastress (SOS-Antwort; Literaturbereich: 10–100×) |
| `base_hgt_rate` | offen | 0.03 | Basisrate für horizontalen Gentransfer pro Zeitschritt |
| `hgt_gene_transfer_prob` | offen | 0.25 | Wahrscheinlichkeit der Genübertragung pro HGT-Ereignis |
| `selection_strength` | offen | 2.5 | Selektionsvorteil resistenter Stämme unter ABX-Druck |

### Schaden & Sterblichkeit

| Parameter | Typ | Wert | Beschreibung |
|---|---|---|---|
| `base_damage_per_step` | offen | 0.004 | Basale Schadensrate pro Zeitschritt |
| `replication_damage_factor` | offen | 0.03 | Zusätzlicher Schaden durch Replikationsstress |
| `stress_damage_factor` | offen | 0.06 | Zusätzlicher Schaden durch Antibiotikastress |
| `repair_rate_per_step` | offen | 0.08 | Reparaturrate pro Zeitschritt |
| `max_damage_load` | offen | 5.0 | Maximale Schadensakkumulation |
| `age_mortality_scale` | offen | 0.001 | Skalierungsfaktor altersabhängiger Mortalität |
| `damage_mortality_scale` | offen | 0.025 | Skalierungsfaktor schadensbedingter Mortalität |
| `lifecycle_half_life_steps` | offen | 200 | Halbwertszeit des Bakterienlebenszyklus in Mikro-Zeitschritten |

### Dormanz & Synergie

| Parameter | Typ | Wert | Beschreibung |
|---|---|---|---|
| `dormancy_growth_penalty` | offen | 0.55 | Wachstumsreduktion während Dormanz |
| `synergy_repair_dormancy_bonus` | offen | 0.25 | Reparaturbonus durch Dormanz-Synergie |
| `synergy_stress_tolerance_bonus` | offen | 0.20 | Stresstoleranzbonus durch Synergie-Effekte |

### Stochastik

| Parameter | Typ | Wert | Beschreibung |
|---|---|---|---|
| `stochastic_threshold` | offen | 10 000 | Populationsgrösse; darunter stochastische statt deterministische Simulation |
| `stochastic_noise_scale` | offen | 0.08 | Skalierung des stochastischen Rauschens bei kleinen Populationen |
