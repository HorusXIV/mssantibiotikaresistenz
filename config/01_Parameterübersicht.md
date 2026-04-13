# Parameterübersicht Makro-Kalibrierung

Diese Datei beschreibt die wichtigsten Parameter der Simulation und ordnet sie ein in:

- geschätzt: aus Literatur oder plausiblen Annahmen abgeleitet
- kalibriert: so angepasst, dass beobachtete Zielgrössen reproduziert werden

---

## 1. Population

### carrier_count
- Typ: geschätzt
- Beschreibung: Anzahl initial kolonisierter Patienten

- Grundlage:
  MRSA-Carriage bei Aufnahme etwa 1.7 %

- Quelle:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/

- Interpretation:
  Startzustand der Simulation

---

### susceptible_count
- Typ: geschätzt
- Beschreibung: Anzahl empfänglicher Patienten

- Definition:

$$
S = N - I
$$

- Hinweis:
  Abgeleiteter Parameter

---

## 2. Kontakte und Hygiene

### daily_contact_attempts ($c$)
- Typ: geschätzt
- Beschreibung: relevante Kontakte pro Person und Tag

- Literatur:
  - etwa 23 Kontaktepisoden pro Patient und Tag
  - etwa 14 verschiedene Kontakte pro Tag

- Gewählter Wert:

$$
c = 12
$$

- Begründung:
  Konservative und realistische Annahme

- Quelle:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC9152759/

---

### base_hygiene ($H$)
- Typ: geschätzt
- Beschreibung: Reduktion der Transmission durch Hygiene

- Wertebereich:

$$
0 \leq H \leq 1
$$

- Interpretation:
  - $H = 0$ keine Hygiene
  - $H = 1$ perfekte Vermeidung

- Gewählter Wert:

$$
H = 0.75
$$

- Begründung:
  Gute Standardhygiene im Spital

---

## 3. Transmission

### base_transmission_rate ($\beta_0$)
- Typ: kalibriert
- Beschreibung: Übertragungswahrscheinlichkeit pro Kontakt

---

## 4. Kalibrierungslogik

Modellgleichung:

$$
\lambda = \beta_0 \cdot c \cdot (1 - H) \cdot \frac{I}{N}
$$

Dabei gilt:
- $\lambda$ ist die Akquisitionsrate
- $I/N$ ist der Carrier-Anteil

---

### Zielwert aus Literatur

MRSA-Akquisition auf Station:

$$
4.6 \text{ bis } 5.4 \text{ pro 1000 patient-days}
$$

$$
\lambda \approx 0.0046 \text{ bis } 0.0054 \text{ pro Tag}
$$

- Quelle:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5384532/

---

### Kalibrierung von $\beta_0$

Gegeben:
- $c = 12$
- $H = 0.75$
- $I/N \approx 0.017$

Berechnung:

$$
\beta_0 = \frac{\lambda}{c \cdot (1-H) \cdot (I/N)}
$$

Ergebnis:

$$
\beta_0 \approx 0.09 \text{ bis } 0.10
$$

- Gewählter Bereich:

$$
\beta_0 = 0.08 \text{ bis } 0.10
$$

---

## 5. Aufenthalt und Abgang

### los_mean_ward
- Typ: geschätzt
- Beschreibung: Mittlere Aufenthaltsdauer auf der Normalstation

- Gewählter Wert:

$$
\text{los\_mean\_ward} = 5 \text{ Tage}
$$

- Quelle:
  https://www.bfs.admin.ch/bfs/de/home.assetdetail.34027817.html

---

### los_mean_icu
- Typ: geschätzt
- Beschreibung: Mittlere Aufenthaltsdauer auf der Intensivstation

- Gewählter Wert:

$$
\text{los\_mean\_icu} = 2.8 \text{ Tage}
$$

- Quelle:
  https://pubmed.ncbi.nlm.nih.gov/10890670/

---

### carrier_extension_days
- Typ: geschätzt
- Beschreibung: Verlängerung des Aufenthalts bei erkannter MRSA-Kolonisierung (Kontaktisolation)

- Gewählter Wert:

$$
\text{carrier\_extension\_days} = 12 \text{ Tage}
$$

- Begründung:
  Proxy-Studie zur Aufenthaltsverlängerung durch Isolation

- Quelle:
  https://www.sciencedirect.com/science/article/pii/S187603412030438X

---

### daily_transfer_rate
- Typ: geschätzt
- Beschreibung: Tägliche Wahrscheinlichkeit einer Verlegung in ein anderes Spital

- Berechnung:

$$
\text{daily\_transfer\_rate} = \frac{0.02}{\text{los\_mean\_ward}} = \frac{0.02}{5} = 0.004
$$

- Gewählter Wert:

$$
\text{daily\_transfer\_rate} = 0.004
$$

- Begründung:
  Proxy-Studie, skaliert auf mittlere Aufenthaltsdauer

- Quelle:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3727950/

---

## 6. Antibiotika

### icu_abx_probability
- Typ: geschätzt
- Beschreibung: Tägliche Wahrscheinlichkeit, dass ein ICU-Patient Antibiotika erhält

- Gewählter Wert:

$$
\text{icu\_abx\_probability} = 0.62
$$

- Quelle:
  SwissNoso Point Prevalence Survey 2017 — https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/2_Punktpraevalenzstudie/Report_Point_Prevalence_Survey_2017_of_HAI_and_antimicrobial_use_in_Swiss_acute_care_hospitals.pdf

---

### ward_abx_probability
- Typ: geschätzt
- Beschreibung: Tägliche Wahrscheinlichkeit, dass ein Ward-Patient Antibiotika erhält

- Gewählter Wert:

$$
\text{ward\_abx\_probability} = 0.326
$$

- Quelle:
  SwissNoso Annual Report HAI CH 2023 — https://www.swissnoso.ch/fileadmin/swissnoso/Dokumente/5_Forschung_und_Entwicklung/8_Swissnoso_Publikationen/Swissnoso_Annual_Report_HAI_CH_2023_FULL_REPORT_EN_SINGLE.pdf

---

## 7. Kolonisierung und Clearance

### community_carrier_fraction
- Typ: geschätzt
- Beschreibung: Anteil MRSA-positiver Patienten bei Spitaleintritt (community-acquired)

- Gewählter Wert:

$$
\text{community\_carrier\_fraction} = 0.017
$$

- Quelle:
  https://pubmed.ncbi.nlm.nih.gov/27658666/

---

### p_clearance
- Typ: geschätzt
- Beschreibung: Tägliche Wahrscheinlichkeit der spontanen Dekolonisierung ohne Behandlung

- Berechnung:

$$
\text{p\_clearance} = \frac{1}{255} \approx 0.0039
$$

- Begründung:
  Mittlere Carrierdauer 8.5 Monate ≈ 255 Tage

- Quelle:
  https://academic.oup.com/cid/article/32/10/1393/465089

---

### replacement_resistant_fraction
- Typ: geschätzt
- Beschreibung: Anteil resistenter Stämme bei neu eintretenden Community-Carriern

- Gewählter Wert:

$$
\text{replacement\_resistant\_fraction} = 0.2
$$

- Quellen:
  https://www.scirp.org/reference/referencespapers?referenceid=3296034
  https://www.ecdc.europa.eu/en/about-us/networks/disease-networks-and-laboratory-networks/ears-net-data

---

## 8. Zusammenfassung

| Parameter | Typ | Wert | Rolle |
|---|---|---|---|
| carrier_count | geschätzt | 20 | Startzustand |
| susceptible_count | geschätzt | 980 | abgeleitet |
| daily_contact_attempts | geschätzt | 12 | Kontaktstruktur |
| base_hygiene | geschätzt | 0.75 | Hygieneeffekt |
| base_transmission_rate | kalibriert | 0.09 | Transmission |
| los_mean_ward | geschätzt | 5 | Aufenthaltsdauer Ward |
| los_mean_icu | geschätzt | 2.8 | Aufenthaltsdauer ICU |
| carrier_extension_days | geschätzt | 12 | Isolationsverlängerung |
| daily_transfer_rate | geschätzt | 0.004 | Verlegungsrate |
| icu_abx_probability | geschätzt | 0.62 | ABX-Rate ICU |
| ward_abx_probability | geschätzt | 0.326 | ABX-Rate Ward |
| community_carrier_fraction | geschätzt | 0.017 | Aufnahmeprävalenz |
| p_clearance | geschätzt | 0.0039 | Dekolonisierungsrate |
| replacement_resistant_fraction | geschätzt | 0.2 | Resistenzanteil Community |
| daily_admission_rate | offen | — | noch nicht bestimmt |
| base_carrier_mortality_rate | offen | — | noch nicht bestimmt |

---

## 9. Wichtiger Hinweis

Die Ausgangsstudie liefert:
- Aufnahmeprävalenz etwa 1.7 %
- Akquisitionsrate etwa 5 pro 1000 patient-days

Sie liefert keine feste Verteilung auf Station

Daher gilt:
- Carrier-Anteil ist geschätzt
- $\beta_0$ ist kalibriert
