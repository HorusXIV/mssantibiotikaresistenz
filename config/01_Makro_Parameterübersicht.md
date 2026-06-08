# Makro-Parameterübersicht

Die Parameter der Mikro-Ebene (Within-Host-Evolution) stehen in der eigenen Referenz `config/02_Mikro_Parameterübersicht.md`.

Jeder Parameter ist einem der folgenden Typen zugeordnet:

- **geschätzt** – aus Literatur oder plausiblen Annahmen abgeleitet
- **Kontextualisierungsparameter** – legt die strukturelle Skala und Topologie der Simulation fest (Gitter, Netzwerk, Kapazität). Modellierungstechnische Setup-Entscheidung, keine aus Daten geschätzte epidemiologische Grösse
- **kalibriert** – durch Simulation so bestimmt, dass eine beobachtete Zielgrösse reproduziert wird
- **Referenzwert** – neutraler Multiplikator (Standard 1.0). Beim Standardwert ohne eigenständigen Effekt, da mit einem kalibrierten Parameter (z.B. β₀) konfundiert; dient als Stellschraube für Sensitivitäts-Experimente
- **nicht identifizierbar** – Sweep durchgeführt, Modell reagiert im plausiblen Bereich nicht signifikant; Wert physikalisch oder neutral motiviert

---

## Mechanismen-Steuerung (Szenarien)

Jeder Modellmechanismus lässt sich rein über die Config gezielt abschalten, ohne
Code-Änderung. Das erlaubt Ablations- und Szenario-Läufe (z.B. "ohne Isolation",
"ohne Antibiotikadruck", "nur Makro"). Die Abschaltung erfolgt entweder über den
Master-Schalter `run.use_micro` oder über das Nullsetzen der jeweiligen Rate; die
Engine fängt den Nullfall sauber ab (kein Sonderzweig nötig).

| Mechanismus | Ebene | Schalter (Config) | Verhalten bei "aus" | Code |
|---|---|---|---|---|
| Mikro gesamt | – | `run.use_micro: false` | Carrier behalten Template-Defaults (`p_clearance`, `relative_transmissibility`, `severity_modifier`); keine Within-Host-Evolution | `run_coupled_simulation.py` (`macro.step(micro_simulator=None)`) |
| Mutation | Mikro | `micro.base_mutation_rate: 0` | keine neuen Mutanten (`Poisson(0)`) | `engine.py` `mutate_population` |
| Horizontaler Gentransfer | Mikro | `micro.base_hgt_rate: 0` | kein Gentransfer (`hgt_prob = 0`) | `engine.py` `horizontal_gene_transfer` |
| Selektion | Mikro | `micro.selection_strength: 0` | neutrale Selektion (`selection_factor = 1`, keine differenzielle Verstärkung) | `engine.py` `selection_step` |
| Transmission | Makro | `macro.base_transmission_rate: 0` | keine S→C-Übertragung im Spital | `simulator.py` `_do_transmission` |
| Erkennung/Isolation | Makro | `macro.carrier_isolation_probability: 0` | Carrier werden nie erkannt/isoliert | `simulator.py` `_build_context` |
| Isolationswirkung | Makro | `macro.base_isolation_effectiveness: 0` | Isolation findet statt, senkt aber die Transmission nicht | `simulator.py` `_do_transmission` |
| Antibiotika | Makro | `macro.icu_abx_probability: 0` + `macro.ward_abx_probability: 0` | kein ABX-Regime, kein Selektionsdruck im Mikro | `simulator.py` `_build_context` |
| Mortalität | Makro | `macro.base_mortality_rate: 0` | keine tägliche Sterblichkeit (Guard `> 0`) | `simulator.py` Discharge-Loop |
| Aufnahmen | Makro | `macro.daily_admission_rate: 0` | keine Neuaufnahmen **und** keine logistische Entlassung (geschlossene Kohorte) | `simulator.py` Admissions (Guard `<= 0`) |
| Community-Import | Makro | `macro.community_carrier_fraction: 0` | Neuaufnahmen sind alle Susceptible | `simulator.py` Admissions |
| Verlegungen | Makro | `macro.daily_transfer_rate: 0` | keine Inter-Spital-Verlegung (Guard `<= 0`) | `simulator.py` Transfer |

`seed` ist **kein** inhaltlicher Parameter, sondern eine rein technische Kontrolle für
Reproduzierbarkeit. Ein gekoppelter Lauf (`mss-run`) nutzt einen festen Seed; die
Robustheit gegenüber dem Zufall prüfen die Kalibrierungs- und Sweep-Werkzeuge über
Ensembles vieler Seeds (`mss-calibrate --n-runs`, `mss-sweep` mit `n_seeds`) mit
Perzentil-/Konfidenzbändern (siehe `docs/03_Modellverhalten_und_Methodik.md`).

---

## Kalibrierungen

> **Begriffliche Einordnung:** Die drei Kalibrierungen passen je einen Parameter an einen externen Literatur-**Zielwert** an. Das ist eine *Kalibrierung gegen Literaturanker*, keine *Validierung* gegen unabhängige, nicht zum Fitten genutzte Daten. cal1 prüft β₀ zusätzlich per Simulation gegen denselben λ-Zielwert (Konsistenz-Check, keine Out-of-Sample-Validierung).

### Kalibrierung 1: `base_transmission_rate` (β₀)

**Script:** `src/mss/cli/run_single_ward_calibration.py`
**Methode:** Analytische Rückrechnung aus Literatur-Zielwert, validiert mit stochastischer Einzelgitter-Simulation (Mikro deaktiviert)

Modellgleichung:

$$
\lambda = \beta_0 \cdot c \cdot (1 - H) \cdot \frac{I}{N}
$$

**Zielwert:** MRSA-Akquisitionsrate auf Station = 4.6–5.4 pro 1000 Patienten-Tage, d.h. λ ≈ 0.0046–0.0054 pro Tag, [PMC 5384532](https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/)

Gegebene Werte: c = 14, H = 0.75, I/N = 0.02 (Startanteil der Single-Ward-Kalibrierungspopulation: carrier_count 2 / 100, **nicht** der Community-Import-Wert 0.017)

$$
\beta_0 = \frac{\lambda}{c \cdot (1 - H) \cdot (I/N)} \approx 0.066 \text{ bis } 0.077
$$

**Ergebnis: β₀ = 0.07**

---

### Kalibrierung 2: `proximity_decay_alpha`

**Script:** `mss-sweep --sweep config/calibration/cal2_proximity_decay.yml`
**Methode:** Stochastischer Parameter-Sweep (`use_micro = false`, n\_seeds = 100, run\_days = 90).

**Zielgrösse:** `same_cell_transmission_fraction`, Anteil der In-Hospital-Übertragungen, deren Quelle in derselben Gitterzelle (Zimmernachbar) liegt. Diese Grösse ist **β₀-unabhängig** (β₀ ist ein globaler Skalar und kürzt sich in der normierten Quellenauswahl heraus), daher die eigene Kalibrierungsachse. Die Unabhängigkeit gilt **erster Ordnung**: Isolation und stammspezifische Transmissibilität sind nicht über alle Carrier uniform und die Carrier-Geographie ist endogen (Clustering), wirken also nur zweiter Ordnung. Prävalenz/Akquisition eignen sich **nicht** als Ziel (vom Community-Import-Floor dominiert).

| | |
|---|---|
| Sweep-Bereich | 0.8 – 1.2 (100 Seeds) |
| Zielgrösse | `same_cell_transmission_fraction = 0.25` |
| Ziel-Begründung | Modellgrösse: Anteil der nosokomialen MRSA-Übertragungen mit Quelle im selben Zimmer. Literatur-Korridor ~20–36 %: [Price et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3922217/) (ICU-WGS) 19 % patient-to-patient als untere Orientierung; [Stone et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22561709/) (VA-LTCF, PFGE) 36 % als Obergrenze; [Kong et al. 2013](https://doi.org/10.1186/1471-2334-13-449) schwache räumliche Struktur als Untergrenze. Zielwert 0.25 liegt im unteren Teil dieses Korridors, nahe dem ICU-WGS-Wert. |
| Sweep-Befund | Monotone Kurve (100 Seeds, explizite Gitterpunkte 0.80…1.20 in 0.05-Schritten): 0.190 (α=0.80) → 0.232 (α=1.00) → 0.248 (α=1.05) → 0.294 (α=1.20). Zielwert 0.25 wird beim Gitterpunkt α = 1.05 am besten getroffen. |
| Einschränkung | Eine Gitterzelle ist ein **Mehrbettzimmer**; der glatte Distanz-Abfall über die Station ist eine Modellvereinfachung, die v.a. den Zimmer-Effekt abbildet. |
| **Gewählter Wert** | **1.05** → same_cell_transmission_fraction = 0.248 ± 0.040 (±1σ über 100 Seeds); das ±1σ-Band [0.21, 0.29] liegt im Literatur-Korridor 20–36 % (Zimmernachbar-Anteil ~25 %) |

---

### Kalibrierung 3: `base_isolation_effectiveness`

**Script:** `mss-sweep --sweep config/calibration/cal3_isolation_effectiveness.yml`
**Methode:** Stochastischer Parameter-Sweep (`use_micro = false`, n\_seeds = 30, run\_days = 90) mit Counterfactual-Baseline.

**Zielgrösse:** `acquisition_reduction`, relative Senkung der kumulativen In-Hospital-Übertragungen bei wirksamer gegenüber unwirksamer Isolation: 1 − n\_an / n\_aus. Im Verhältnis kürzen sich **β₀, Kontakte, Hygiene und Raumstruktur heraus** → eigene, β₀-unabhängige Achse (das absolute Niveau ist mit β₀ konfundiert und scheidet aus). Der Baseline-Lauf verwendet dieselbe Isolationspolitik, aber Wirksamkeit 0 (`base_isolation_effectiveness = 0`), so kürzt sich die LOS-Verlängerung durch Erkennung heraus und das Verhältnis misst **nur** die Wirksamkeit. Pro Seed gepaart.

| | |
|---|---|
| Sweep-Bereich | 0.30 – 0.90 (30 Seeds) |
| Zielgrösse | `acquisition_reduction = 0.60` |
| Ziel-Begründung | Kontaktisolation senkt die MRSA-Akquisition populationsweit um ~20–60 % (Bandbreite über Studien); 0.60 = oberes Band. Über Übertragungsketten ist die populationsweite Reduktion grösser als die Per-Patient-Wirksamkeit, daher fällt der kalibrierte Per-Patient-Wert beim oberen Reduktionsband in den konsistenten Literaturbereich. |
| Sweep-Befund | Monotone Kurve über den Sweep-Bereich; der Zielwert 0.60 wird bei eff = 0.45 erreicht (acquisition_reduction = 0.599). |
| Einschränkung | Konditional auf `carrier_isolation_probability = 0.65` (literaturbelegte Isolationsabdeckung), Wirksamkeit und Abdeckung sind konfundiert, letztere bleibt fixiert. |
| **Gewählter Wert** | **0.45** → acquisition_reduction = 0.599 ± 0.059 (±1σ über 30 Seeds, Band [0.54, 0.66] umschliesst das Ziel 0.60). Der Per-Patient-Wert liegt knapp unter dem berichteten Bereich (0.5–0.9) und ist über die populationsweite Verstärkung der Übertragungsketten konsistent. |

---

## Populations-Parameter

`initial_department` ist bei Simulationsstart für alle Patienten auf `ward` gesetzt.

### Populationsgrösse

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `hospitals` | Kontextualisierungsparameter | 6 | Anzahl Spitäler im simulierten Netzwerk | - |
| `carrier_count` | geschätzt | 20 (≈ 2 % von N) | Anzahl initial kolonisierter Patienten bei Simulationsstart | [PMC 5384532](https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/) |
| `susceptible_count` | geschätzt | 980 | Anzahl empfänglicher Patienten bei Simulationsstart | abgeleitet |

### Agenten-Initialwerte: Susceptible

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `compliance` | geschätzt | 0.80 | Allgemeine Verhaltens-Compliance des Patienten gegenüber medizinischen Massnahmen [0–1]. Wird täglich in `adherence` umgewandelt (1:1-Mapping), welches die effektive Antibiotikawirkung im Mikro-Layer skaliert (`effective_kill = base_kill_rate × dose × adherence`). Kein direkter Effekt auf Transmission. Wert orientiert sich an beobachteter Protokoll-Compliance im Akutspital. | [Swissnoso Report 2023](https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/8_Swissnoso_Publikationen/Swissnoso_Annual_Report_HAI_CH_2023_FULL_REPORT_EN_DOUBLE.pdf) |
| `immune_strength` | Referenzwert | 1.0 | Immunkompetenz-Multiplikator. Steuert Makro-Suszeptibilität (`1/immune_strength`) und Mikro-Clearance, ist aber mit β₀ konfundiert → beim Standard 1.0 ohne eigenständigen Effekt (Multiplikator 1.0 = Identität). Eigenständige Identifizierbarkeit ist analytisch/definitional ausgeschlossen; dient als **Referenz-Stellschraube** für spätere Sensitivitäts-Experimente. | Modellannahme |
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
| `dept_grid_cols` | Kontextualisierungsparameter | 5 | Zimmer-Anordnung (5 × 4 = 20 Betten/Zimmer pro Ward) | - |
| `dept_grid_rows` | Kontextualisierungsparameter | 4 | Gesamtreihen im Gitter | - |
| `dept_grid_icu_rows` | Kontextualisierungsparameter | 1 | Ca. 25 % ICU-Betten-Anteil | - |
| `network_grid_cols` | Kontextualisierungsparameter | 3 | 3 × 2-Cluster für 6 Spitäler | - |
| `max_occupancy_per_hospital` | Kontextualisierungsparameter | 200 | Hard Cap für Aufnahmen und Transfers pro Spital | - |
| `proximity_decay_alpha` | kalibriert | 1.05 | Abklingkoeffizient α: w = exp(−α·d). Kalibriert gegen den Zimmernachbar-Übertragungsanteil (`same_cell_transmission_fraction` ≈ 0.25). | Kalibrierung 2; [Price et al. 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3922217/); [Stone et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22561709/); [Kong et al. 2013](https://doi.org/10.1186/1471-2334-13-449) |

### Transmission

| Parameter | Typ | Wert | Beschreibung | Quelle |
|---|---|---|---|---|
| `daily_contact_attempts` | geschätzt | 14 | Relevante Kontaktepisoden pro Patient und Tag | [PMC 9152759](https://pmc.ncbi.nlm.nih.gov/articles/PMC9152759/) |
| `base_hygiene` | geschätzt | 0.75 | Reduktion der Transmissionswahrscheinlichkeit durch Hygienemassnahmen | - |
| `base_transmission_rate` | kalibriert | 0.07 | Übertragungswahrscheinlichkeit pro relevantem Kontakt | Kalibrierung 1 |
| `carrier_isolation_probability` | geschätzt | 0.65 | Wahrscheinlichkeit, dass ein erkannter Carrier täglich isoliert wird | [HIS/IPS MRSA Guidelines](https://www.his.org.uk/media/djom3t4f/joint-healthcare-infection-society-his-and-infection-prevention-society-ips-guidelines-for-the-prevention-and-control-of-meticillin-resistant-staphylococcus-aureus-mrsa-in-healthcare-facilities.pdf) |
| `base_isolation_effectiveness` | kalibriert | 0.45 | Per-Patient-Transmissionsreduktion für *isolierte* Carrier. Kalibriert über die **relative** Akquisitionsreduktion durch Isolation (Kalibrierung 3, bester Gitterpunkt), das absolute Niveau wäre mit β₀/proximity konfundiert. Ergebnis liegt knapp unter dem berichteten Per-Patient-Bereich (0.5–0.9). | Kalibrierung 3; [Marshall et al. 2013, PLoS One](https://pmc.ncbi.nlm.nih.gov/articles/PMC4049803/); [ICU-Modellstudie, J Hosp Infect 2023](https://www.sciencedirect.com/science/article/abs/pii/S0195670123003742) |
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
