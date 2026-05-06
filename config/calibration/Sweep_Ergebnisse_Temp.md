# phase1_isolation_effectiveness:

  min: 0.40
  max: 0.95
  n_steps: 8

  Parameter-Sweep: macro.base_isolation_effectiveness
  Werte: [0.400, 0.479, 0.557, 0.636, 0.714, 0.793, 0.871, 0.950]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Modellannahme für Schweizer Niedrig-MRSA-Setting (plausibler Bereich 0.5–2.0/1000 Pat.-Tage; abgeleitet aus Swissnoso/EARS-Net Kontextdaten))
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  macro.base_isolation_effectiveness=0.4000  →  acquisition_rate_per_1000=3.5064  (Abstand=2.0064)
  macro.base_isolation_effectiveness=0.4786  →  acquisition_rate_per_1000=2.9638  (Abstand=1.4638)
  macro.base_isolation_effectiveness=0.5571  →  acquisition_rate_per_1000=2.8618  (Abstand=1.3618)
  macro.base_isolation_effectiveness=0.6357  →  acquisition_rate_per_1000=3.2196  (Abstand=1.7196)
  macro.base_isolation_effectiveness=0.7143  →  acquisition_rate_per_1000=2.4068  (Abstand=0.9068)
  macro.base_isolation_effectiveness=0.7929  →  acquisition_rate_per_1000=2.3589  (Abstand=0.8589)
  macro.base_isolation_effectiveness=0.8714  →  acquisition_rate_per_1000=2.1339  (Abstand=0.6339)
  macro.base_isolation_effectiveness=0.9500  →  acquisition_rate_per_1000=1.7944  (Abstand=0.2944)

BESTER WERT: macro.base_isolation_effectiveness = 0.9500  →  acquisition_rate_per_1000 = 1.7944  (Abstand zum Ziel: 0.2944)


Eingeschränkter:
  min: 0.95
  max: 1
  n_steps: 8

Parameter-Sweep: macro.base_isolation_effectiveness
  Werte: [0.950, 0.957, 0.964, 0.971, 0.979, 0.986, 0.993, 1.000]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Modellannahme für Schweizer Niedrig-MRSA-Setting (plausibler Bereich 0.5–2.0/1000 Pat.-Tage; abgeleitet aus Swissnoso/EARS-Net Kontextdaten))
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  macro.base_isolation_effectiveness=0.9500  →  acquisition_rate_per_1000=1.7944  (Abstand=0.2944)
  macro.base_isolation_effectiveness=0.9571  →  acquisition_rate_per_1000=1.8150  (Abstand=0.3150)
  macro.base_isolation_effectiveness=0.9643  →  acquisition_rate_per_1000=1.6728  (Abstand=0.1728)
  macro.base_isolation_effectiveness=0.9714  →  acquisition_rate_per_1000=1.4470  (Abstand=0.0530)
  macro.base_isolation_effectiveness=0.9786  →  acquisition_rate_per_1000=1.4143  (Abstand=0.0857)
  macro.base_isolation_effectiveness=0.9857  →  acquisition_rate_per_1000=1.5302  (Abstand=0.0302)
  macro.base_isolation_effectiveness=0.9929  →  acquisition_rate_per_1000=1.5741  (Abstand=0.0741)
  macro.base_isolation_effectiveness=1.0000  →  acquisition_rate_per_1000=1.5502  (Abstand=0.0502)

BESTER WERT: macro.base_isolation_effectiveness = 0.9857  →  acquisition_rate_per_1000 = 1.5302  (Abstand zum Ziel: 0.0302)

---

# phase1_proximity_decay

Parameter-Sweep: macro.proximity_decay_alpha
  Werte: [0.100, 0.443, 0.786, 1.129, 1.471, 1.814, 2.157, 2.500]
  Zielgroesse: mean_prevalence = 0.02  (Swissnoso 2022 / Rohner et al. 2017 (1.5–4 % MRSA-Praevalenz CH))
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  macro.proximity_decay_alpha=0.1000  →  mean_prevalence=0.1462  (Abstand=0.1262)
  macro.proximity_decay_alpha=0.4429  →  mean_prevalence=0.1370  (Abstand=0.1170)
  macro.proximity_decay_alpha=0.7857  →  mean_prevalence=0.1375  (Abstand=0.1175)
  macro.proximity_decay_alpha=1.1286  →  mean_prevalence=0.1205  (Abstand=0.1005)
  macro.proximity_decay_alpha=1.4714  →  mean_prevalence=0.1255  (Abstand=0.1055)
  macro.proximity_decay_alpha=1.8143  →  mean_prevalence=0.1135  (Abstand=0.0935)
  macro.proximity_decay_alpha=2.1571  →  mean_prevalence=0.1168  (Abstand=0.0968)
  macro.proximity_decay_alpha=2.5000  →  mean_prevalence=0.1155  (Abstand=0.0955)

BESTER WERT: macro.proximity_decay_alpha = 1.8143  →  mean_prevalence = 0.1135  (Abstand zum Ziel: 0.0935)


---

Parameter-Sweep: population.carrier_template.sociability
  Werte: [0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Swissnoso 2022 / ECDC EARS-Net 2023)
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  population.carrier_template.sociability=0.4000  →  acquisition_rate_per_1000=1.6348  (Abstand=0.1348)
  population.carrier_template.sociability=0.6000  →  acquisition_rate_per_1000=1.6348  (Abstand=0.1348)
  population.carrier_template.sociability=0.8000  →  acquisition_rate_per_1000=1.6348  (Abstand=0.1348)
  population.carrier_template.sociability=1.0000  →  acquisition_rate_per_1000=1.5722  (Abstand=0.0722)
  population.carrier_template.sociability=1.2000  →  acquisition_rate_per_1000=1.5923  (Abstand=0.0923)
  population.carrier_template.sociability=1.4000  →  acquisition_rate_per_1000=1.5343  (Abstand=0.0343)
  population.carrier_template.sociability=1.6000  →  acquisition_rate_per_1000=1.5343  (Abstand=0.0343)
  population.carrier_template.sociability=1.8000  →  acquisition_rate_per_1000=1.8709  (Abstand=0.3709)
  population.carrier_template.sociability=2.0000  →  acquisition_rate_per_1000=2.1946  (Abstand=0.6946)

BESTER WERT: population.carrier_template.sociability = 1.4000  →  acquisition_rate_per_1000 = 1.5343  (Abstand zum Ziel: 0.0343)


# phase2_carrier_sociaility.yml
Parameter-Sweep: population.carrier_template.sociability
  Werte: [0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Swissnoso 2022 / ECDC EARS-Net 2023)
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  population.carrier_template.sociability=0.4000  →  acquisition_rate_per_1000=1.6348  (Abstand=0.1348)
  population.carrier_template.sociability=0.6000  →  acquisition_rate_per_1000=1.6348  (Abstand=0.1348)
  population.carrier_template.sociability=0.8000  →  acquisition_rate_per_1000=1.6348  (Abstand=0.1348)
  population.carrier_template.sociability=1.0000  →  acquisition_rate_per_1000=1.5722  (Abstand=0.0722)
  population.carrier_template.sociability=1.2000  →  acquisition_rate_per_1000=1.5923  (Abstand=0.0923)
  population.carrier_template.sociability=1.4000  →  acquisition_rate_per_1000=1.5343  (Abstand=0.0343)
  population.carrier_template.sociability=1.6000  →  acquisition_rate_per_1000=1.5343  (Abstand=0.0343)
  population.carrier_template.sociability=1.8000  →  acquisition_rate_per_1000=1.8709  (Abstand=0.3709)
  population.carrier_template.sociability=2.0000  →  acquisition_rate_per_1000=2.1946  (Abstand=0.6946)

BESTER WERT: population.carrier_template.sociability = 1.4000  →  acquisition_rate_per_1000 = 1.5343  (Abstand zum Ziel: 0.0343)

# phase2_susceptible_immune_strength.yml

Parameter-Sweep: population.susceptible_template.immune_strength
  Werte: [0.300, 0.543, 0.786, 1.029, 1.271, 1.514, 1.757, 2.000]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Swissnoso 2022 / ECDC EARS-Net 2023)
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  population.susceptible_template.immune_strength=0.3000  →  acquisition_rate_per_1000=2.0855  (Abstand=0.5855)
  population.susceptible_template.immune_strength=0.5429  →  acquisition_rate_per_1000=1.6070  (Abstand=0.1070)
  population.susceptible_template.immune_strength=0.7857  →  acquisition_rate_per_1000=1.6320  (Abstand=0.1320)
  population.susceptible_template.immune_strength=1.0286  →  acquisition_rate_per_1000=1.5733  (Abstand=0.0733)
  population.susceptible_template.immune_strength=1.2714  →  acquisition_rate_per_1000=1.7268  (Abstand=0.2268)
  population.susceptible_template.immune_strength=1.5143  →  acquisition_rate_per_1000=1.5985  (Abstand=0.0985)
  population.susceptible_template.immune_strength=1.7571  →  acquisition_rate_per_1000=1.7303  (Abstand=0.2303)
  population.susceptible_template.immune_strength=2.0000  →  acquisition_rate_per_1000=1.7303  (Abstand=0.2303)

BESTER WERT: population.susceptible_template.immune_strength = 1.0286  →  acquisition_rate_per_1000 = 1.5733  (Abstand zum Ziel: 0.0733)


# phase2_susceptible_sociability.yml

Parameter-Sweep: population.susceptible_template.sociability
  Werte: [0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Swissnoso 2022 / ECDC EARS-Net 2023)
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  population.susceptible_template.sociability=0.4000  →  acquisition_rate_per_1000=1.5681  (Abstand=0.0681)
  population.susceptible_template.sociability=0.6000  →  acquisition_rate_per_1000=1.6474  (Abstand=0.1474)
  population.susceptible_template.sociability=0.8000  →  acquisition_rate_per_1000=1.7152  (Abstand=0.2152)
  population.susceptible_template.sociability=1.0000  →  acquisition_rate_per_1000=1.5733  (Abstand=0.0733)
  population.susceptible_template.sociability=1.2000  →  acquisition_rate_per_1000=1.5729  (Abstand=0.0729)
  population.susceptible_template.sociability=1.4000  →  acquisition_rate_per_1000=1.5652  (Abstand=0.0652)
  population.susceptible_template.sociability=1.6000  →  acquisition_rate_per_1000=1.6690  (Abstand=0.1690)
  population.susceptible_template.sociability=1.8000  →  acquisition_rate_per_1000=1.5630  (Abstand=0.0630)
  population.susceptible_template.sociability=2.0000  →  acquisition_rate_per_1000=1.7416  (Abstand=0.2416)

BESTER WERT: population.susceptible_template.sociability = 1.8000  →  acquisition_rate_per_1000 = 1.5630  (Abstand zum Ziel: 0.0630)

# phase2_susceptible_vulnerability.yml

Parameter-Sweep: population.susceptible_template.vulnerability
  Werte: [0.400, 0.700, 1.000, 1.300, 1.600, 1.900, 2.200, 2.500]
  Zielgroesse: acquisition_rate_per_1000 = 1.5  (Swissnoso 2022 / ECDC EARS-Net 2023)
  Seeds pro Rasterpunkt: 5
  Simulationsdauer: 90 Tage
  Mikro-Simulator: deaktiviert (use_micro: false)
  Config: /Users/damianszedalik/Documents/01_Coding/02_FHNW/07_Angewandte_Data_Science/mss/mssantibiotikaresistenz/config/simulation_realistic.yml

  population.susceptible_template.vulnerability=0.4000  →  acquisition_rate_per_1000=1.6186  (Abstand=0.1186)
  population.susceptible_template.vulnerability=0.7000  →  acquisition_rate_per_1000=1.5985  (Abstand=0.0985)
  population.susceptible_template.vulnerability=1.0000  →  acquisition_rate_per_1000=1.5699  (Abstand=0.0699)
  population.susceptible_template.vulnerability=1.3000  →  acquisition_rate_per_1000=1.5401  (Abstand=0.0401)
  population.susceptible_template.vulnerability=1.6000  →  acquisition_rate_per_1000=1.6545  (Abstand=0.1545)
  population.susceptible_template.vulnerability=1.9000  →  acquisition_rate_per_1000=1.6616  (Abstand=0.1616)
  population.susceptible_template.vulnerability=2.2000  →  acquisition_rate_per_1000=1.5058  (Abstand=0.0058)
  population.susceptible_template.vulnerability=2.5000  →  acquisition_rate_per_1000=1.6722  (Abstand=0.1722)

BESTER WERT: population.susceptible_template.vulnerability = 2.2000  →  acquisition_rate_per_1000 = 1.5058  (Abstand zum Ziel: 0.0058)
