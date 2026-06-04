# Parameterübersicht

Jeder Parameter ist einem der folgenden Typen zugeordnet:

- **geschätzt** – aus Literatur oder plausiblen Annahmen abgeleitet
- **Kontextualisierungsparameter** – legt die strukturelle Skala und Topologie der Simulation fest (Gitter, Netzwerk, Kapazität). Modellierungstechnische Setup-Entscheidung, keine aus Daten geschätzte epidemiologische Grösse
- **kalibriert** – durch Simulation so bestimmt, dass eine beobachtete Zielgrösse reproduziert wird
- **Referenzwert** – neutraler Multiplikator (Standard 1.0). Beim Standardwert ohne eigenständigen Effekt, da mit einem kalibrierten Parameter (z.B. β₀) konfundiert; dient als Stellschraube für Sensitivitäts-Experimente
- **nicht identifizierbar** – Sweep durchgeführt, Modell reagiert im plausiblen Bereich nicht signifikant; Wert physikalisch oder neutral motiviert
- **offen** – noch nicht bestimmt (Mikro-Parameter)

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

### Kalibrierung 2 — `proximity_decay_alpha`

**Script:** `mss-sweep --sweep config/calibration/cal2_proximity_decay.yml`
**Methode:** Stochastischer Parameter-Sweep (`use_micro = false`, n\_seeds = 5, run\_days = 90).

**Zielgrösse:** `same_cell_transmission_fraction` — Anteil der In-Hospital-Übertragungen, deren Quelle in derselben Gitterzelle (Zimmernachbar) liegt. Diese Grösse wird **allein von α** bestimmt und ist **nicht** mit β₀/Isolation konfundiert (die ändern das Niveau der Übertragung, nicht die Distanzverteilung) — daher die eigene, isolierte Kalibrierungsachse. Prävalenz/Akquisition eignen sich **nicht** als Ziel (vom Community-Import-Floor dominiert).

| | |
|---|---|
| Sweep-Bereich | 0.8 – 1.2 (50 Seeds) |
| Zielgrösse | `same_cell_transmission_fraction = 0.25` |
| Ziel-Begründung | Zimmernachbar-attribuierbarer Anteil nosokomialer MRSA-Akquisition ~20–30 % (VA-LTCF ~36 % als Obergrenze; ICU-WGS ~20 % patient-to-patient; hierarchisches Raummodell: schwache räumliche Struktur als Untergrenze) |
| Sweep-Befund | Monotone Kurve (50 Seeds): 0.182 (α=0.8) → 0.292 (α=1.2). Zielwert 0.25 wird bei α ≈ 1.05 erreicht (α=1.03 → 0.240, α=1.09 → 0.257). |
| Einschränkung | Das Gitter ist eine Abstraktion (Zelle ≈ Zimmer/Bettzone); der glatte Distanz-Abfall ist eine Vereinfachung — die Literatur stützt v.a. den Zimmer-Effekt, nicht einen Gradienten über die ganze Station. |
| **Gewählter Wert** | **1.05** → same_cell_transmission_fraction ≈ 0.246 (Zimmernachbar-Anteil ~25 %) |

---

### Kalibrierung 3 — `base_isolation_effectiveness`

**Script:** `mss-sweep --sweep config/calibration/cal3_isolation_effectiveness.yml`
**Methode:** Stochastischer Parameter-Sweep (`use_micro = false`, n\_seeds = 30, run\_days = 90) mit Counterfactual-Baseline.

**Zielgrösse:** `acquisition_reduction` — relative Senkung der kumulativen In-Hospital-Übertragungen bei wirksamer gegenüber unwirksamer Isolation: 1 − n\_an / n\_aus. Im Verhältnis kürzen sich **β₀, Kontakte, Hygiene und Raumstruktur heraus** → eigene, β₀-unabhängige Achse (das absolute Niveau ist mit β₀ konfundiert und scheidet aus). Der Baseline-Lauf verwendet dieselbe Isolationspolitik, aber Wirksamkeit 0 (`base_isolation_effectiveness = 0`) — so kürzt sich die LOS-Verlängerung durch Erkennung heraus und das Verhältnis misst **nur** die Wirksamkeit. Pro Seed gepaart.

| | |
|---|---|
| Sweep-Bereich | 0.30 – 0.90 (30 Seeds) |
| Zielgrösse | `acquisition_reduction = 0.60` |
| Ziel-Begründung | Kontaktisolation senkt die MRSA-Akquisition populationsweit um ~20–60 % (Bandbreite über Studien); 0.60 = oberes Band. Über Übertragungsketten ist die populationsweite Reduktion grösser als die Per-Patient-Wirksamkeit — daher fällt der kalibrierte Per-Patient-Wert beim oberen Reduktionsband in den konsistenten Literaturbereich. |
| Sweep-Befund | Monotone Kurve (30 Seeds): 0.44 (eff=0.30) → 0.90 (eff=0.90). Zielwert 0.60 wird bei eff ≈ 0.47 erreicht. |
| Einschränkung | Konditional auf `carrier_isolation_probability = 0.65` (literaturbelegte Isolationsabdeckung) — Wirksamkeit und Abdeckung sind konfundiert, letztere bleibt fixiert. |
| **Gewählter Wert** | **0.47** → acquisition_reduction ≈ 0.60; liegt am unteren Rand der Per-Patient-Literatur (0.5–0.9) |

---

## Populations-Parameter

`initial_department` ist bei Simulationsstart für alle Patienten auf `ward` gesetzt.

### Populationsgrösse

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `hospitals` | Kontextualisierungsparameter | 6 | Anzahl Spitäler im simulierten Netzwerk | — |
| `carrier_count` | geschätzt | 20 (≈ 2 % von N) | Anzahl initial kolonisierter Patienten bei Simulationsstart | [PMC 5384532](https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/) |
| `susceptible_count` | geschätzt | 980 | Anzahl empfänglicher Patienten bei Simulationsstart | abgeleitet |

### Agenten-Initialwerte: Susceptible

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `compliance` | geschätzt | 0.80 | Allgemeine Verhaltens-Compliance des Patienten gegenüber medizinischen Massnahmen [0–1]. Wird täglich in `adherence` umgewandelt (1:1-Mapping), welches die effektive Antibiotikawirkung im Mikro-Layer skaliert (`effective_kill = base_kill_rate × dose × adherence`). Kein direkter Effekt auf Transmission. Wert orientiert sich an beobachteter Protokoll-Compliance im Akutspital. | [Swissnoso Report 2023](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/8_Swissnoso_Publikationen/Swissnoso_Annual_Report_HAI_CH_2023_FULL_REPORT_EN_DOUBLE.pdf) |
| `immune_strength` | Referenzwert | 1.0 | Immunkompetenz-Multiplikator. Steuert Makro-Suszeptibilität (`1/immune_strength`) und Mikro-Clearance, ist aber mit β₀ konfundiert → beim Standard 1.0 ohne eigenständigen Effekt. Sensitivitätsanalyse bestätigte Nicht-Identifizierbarkeit. | Modellannahme |
| `sociability` | Referenzwert | 1.0 | Kontaktintensitäts-Multiplikator (skaliert Transmission, sobald die Person Carrier wird). Mit β₀ konfundiert → beim Standard 1.0 ohne eigenständigen Effekt. | Modellannahme |

### Agenten-Initialwerte: Carrier

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `compliance` | geschätzt | 0.70 | Allgemeine Verhaltens-Compliance [0–1]; leicht tiefer als bei Susceptible (persistente Carrier zeigen erfahrungsgemäss schlechtere Therapietreue). Wirkt ausschliesslich via `adherence` auf die Antibiotikawirkung im Mikro-Layer. | Modellannahme |
| `immune_strength` | Referenzwert | 0.75 | Immunkompetenz-Multiplikator; leicht reduziert für persistente Carrier (wirkt im Mikro auf Immun-Clearance). Mit der Mikro-Basis-Clearancerate konfundiert → nicht eigenständig identifizierbar, daher als Referenzwert gesetzt. | Modellannahme |
| `sociability` | Referenzwert | 1.0 | Kontaktintensitäts-Multiplikator (skaliert Carrier-Transmission). Mit β₀ konfundiert → beim Standard 1.0 ohne eigenständigen Effekt. | Modellannahme |
| `resistant_fraction` | geschätzt | 0.90 | Anteil MRSA an Gesamt-S.aureus innerhalb eines bestätigten MRSA-Carriers. Bei >90% der MRSA-positiven Patienten liegt der MSSA-Anteil unter der Nachweisgrenze von 5% (30-Kolonie-Methode). Modelliert nahezu klonale MRSA-Dominanz mit kleiner sensitiver Restpopulation. | [Dall'Antonia et al. 2005, J Hosp Infect](https://www.sciencedirect.com/science/article/pii/S0195670105000691?via%3Dihub) |
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
| `dept_grid_cols` | Kontextualisierungsparameter | 5 | Zimmer-Anordnung (5 × 4 = 20 Betten/Zimmer pro Ward) | — |
| `dept_grid_rows` | Kontextualisierungsparameter | 4 | Gesamtreihen im Gitter | — |
| `dept_grid_icu_rows` | Kontextualisierungsparameter | 1 | Ca. 25 % ICU-Betten-Anteil | — |
| `network_grid_cols` | Kontextualisierungsparameter | 3 | 3 × 2-Cluster für 6 Spitäler | — |
| `max_occupancy_per_hospital` | Kontextualisierungsparameter | 200 | Hard Cap für Aufnahmen und Transfers pro Spital | — |
| `proximity_decay_alpha` | kalibriert | 1.05 | Abklingkoeffizient α: w = exp(−α·d). Kalibriert gegen den Zimmernachbar-Übertragungsanteil (`same_cell_transmission_fraction` ≈ 0.25). | Kalibrierung 2 |

### Transmission

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `daily_contact_attempts` | geschätzt | 14 | Relevante Kontaktepisoden pro Patient und Tag | [PMC 9152759](https://pmc.ncbi.nlm.nih.gov/articles/PMC9152759/) |
| `base_hygiene` | geschätzt | 0.75 | Reduktion der Transmissionswahrscheinlichkeit durch Hygienemassnahmen | — |
| `base_transmission_rate` | kalibriert | 0.07 | Übertragungswahrscheinlichkeit pro relevantem Kontakt | Kalibrierung 1 |
| `carrier_isolation_probability` | geschätzt | 0.65 | Wahrscheinlichkeit, dass ein erkannter Carrier täglich isoliert wird | [HIS/IPS MRSA Guidelines](https://www.his.org.uk/media/djom3t4f/joint-healthcare-infection-society-his-and-infection-prevention-society-ips-guidelines-for-the-prevention-and-control-of-meticillin-resistant-staphylococcus-aureus-mrsa-in-healthcare-facilities.pdf) |
| `base_isolation_effectiveness` | kalibriert | 0.47 | Per-Patient-Transmissionsreduktion für *isolierte* Carrier. Kalibriert über die **relative** Akquisitionsreduktion durch Isolation (Kalibrierung 3) — das absolute Niveau wäre mit β₀/proximity konfundiert. Ergebnis liegt am unteren Rand der Per-Patient-Literatur (0.5–0.9). | Kalibrierung 3; [Marshall et al. 2013, PLoS One](https://pmc.ncbi.nlm.nih.gov/articles/PMC4049803/); [ICU-Modellstudie, J Hosp Infect 2023](https://www.sciencedirect.com/science/article/abs/pii/S0195670123003742) |
| `base_diagnostic_speed` | Referenzwert | 1.0 | Skalierungsfaktor auf die Carrier-Erkennungsrate. Beim Standard 1.0 ohne eigenständigen Effekt (mit `carrier_isolation_probability` konfundiert); Stellschraube für Szenarien: 0.5 = Kultur (Nasenabstrich, langsam), 2.0 = PCR-Schnelltest. | Modellannahme |

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
| `replacement_resistant_fraction` | geschätzt | 0.90 | Anteil MRSA an Gesamt-S.aureus bei neu eintretenden Community-Carriern. Gleiche Begründung wie `resistant_fraction`: bestätigte Carrier haben nahezu klonale MRSA-Dominanz. Dieser Wert wirkt bei jeder neuen Community-Aufnahme (~2.1/Tag) als Initialisierung für Tag 1 der Mikro-Episode. | [Dall'Antonia et al. 2005, J Hosp Infect](https://www.sciencedirect.com/science/article/pii/S0195670105000691?via%3Dihub) |
| `replacement_dominant_genotype` | geschätzt | R2 | Dominanter Genotyp bei neu eintretenden Community-Carriern. CH-MRSA ist dominant der CC22/ST22-MRSA-IV-Klon (Barnim-Stamm); entspricht Resistenzklasse R2. Konsistent mit `carrier_template.dominant_genotype`. | [Warnke et al. 2014, PMC4213029](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4213029/) |
| `base_mortality_rate` | geschätzt | 0.0045 | Tägliche Sterblichkeit für alle Patienten. Für Carrier skaliert mit `severity_modifier` (→ 0.0054/Tag). Herleitung: CH In-Hospital-Sterblichkeit ~2.5 % ÷ mittlere LOS 5.5 Tage = 0.0045/Tag. | [BFS Medizinische Statistik 2022](https://www.bfs.admin.ch/bfs/de/home.assetdetail.34027817.html); [Huang & Platt 2003, Ann Intern Med](https://www.acpjournals.org/doi/10.7326/0003-4819-139-5_Part_1-200309020-00008) |

### Antibiotika

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `icu_abx_probability` | geschätzt | 0.62 | Tägliche Wahrscheinlichkeit, dass ein ICU-Patient Antibiotika erhält | [SwissNoso PPS 2017](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/2_Punktpraevalenzstudie/Report_Point_Prevalence_Survey_2017_of_HAI_and_antimicrobial_use_in_Swiss_acute_care_hospitals.pdf) |
| `ward_abx_probability` | geschätzt | 0.326 | Tägliche Wahrscheinlichkeit, dass ein Ward-Patient Antibiotika erhält | [SwissNoso Annual Report 2023](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/8_Swissnoso_Publikationen/Swissnoso_Annual_Report_HAI_CH_2023_FULL_REPORT_EN_SINGLE.pdf) |

---

## Mikro-Parameter

Die Mikro-Simulation modelliert die bakterielle Evolution innerhalb einzelner Patienten. Ein Makro-Tag wird in `steps_per_day` feinere Zeitschritte unterteilt. Bei den Makro-Kalibrierungen ist der Mikro-Simulator deaktiviert (`use_micro = false`); Patienten behalten ihre Template-Standardwerte. Die Mikro-Parameter sind noch nicht kalibriert (`offen`).

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
