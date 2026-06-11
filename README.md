# MSS: Antibiotikaresistenz-Simulation

MSS ist ein mehrskaliges Simulations-Framework zur Modellierung der Antibiotikaresistenz-Dynamik in Spitalumgebungen. Der Code ist um ein einziges `src/mss`-Paket herum organisiert, damit Simulationslogik, Domänenobjekte und CLI-Einstiegspunkte von Konfiguration, Dokumentation, Tests und generierten Artefakten getrennt sind.

> **Hier starten:** [`docs/00_Overview.md`](docs/00_Overview.md) ist der geführte Einstieg. Es verbindet Makro, Mikro, Patient, Config und Outputs ("ein Patient, ein Tag") und listet die Diagramme auf. Dieses README ist die Repository-Karte und die Ausführungsreferenz.

## Repository-Prinzipien

- Aller gepflegte Python-Code liegt unter `src/mss/`.
- Ordner auf oberster Ebene sind für Konfiguration, Dokumentation, Tests, Container und generierte Artefakte reserviert.
- Makro, Mikro, Domäne und CLI haben klar getrennte Verantwortlichkeiten.

## Repository-Struktur

```text
MSS/
├── Organizational/
│   ├── Mini_Challenge.md
│   └── Modulbeschreibung.md
├── config/
│   ├── calibration/
│   │   ├── cal1_simulation_single_ward.yml
│   │   ├── cal2_proximity_decay.yml
│   │   └── cal3_isolation_effectiveness.yml
│   ├── 01_Makro_Parameterübersicht.md
│   ├── 02_Mikro_Parameterübersicht.md
│   ├── simulation_abx.yml
│   ├── simulation_realistic.yml
│   └── template.yml
├── containers/
│   └── mss_image.def
├── docs/
│   ├── 00_Overview.md
│   ├── 01_Makro_Overview.md
│   ├── 02_Mikro_Overview.md
│   ├── 03_Modellverhalten_und_Methodik.md
│   ├── organizational/
│   └── system_overview/
│       ├── Flowchart_v0.mmd
│       ├── Flowchart_v1.mmd
│       ├── MindMap.mmd
│       ├── Pruned_MindMap.mmd
│       ├── amr_system_map.gexf
│       ├── amr_system_map_edges.csv
│       ├── amr_system_map_nodes.csv
│       └── build_gephi_graphs.py
├── logs/
│   └── *.out
├── notebooks/
│   └── output_explorer.ipynb
├── outputs/
│   └── <timestamp>_<name>/
│       ├── data/
│       │   └── *.parquet
│       └── plots/
│           └── *.png
├── src/
│   └── mss/
│       ├── __init__.py
│       ├── cli/
│       │   ├── __init__.py
│       │   ├── run_coupled_simulation.py
│       │   ├── run_parameter_sweep.py
│       │   ├── run_single_ward_calibration.py
│       │   └── visualize_results.py
│       ├── domain/
│       │   ├── __init__.py
│       │   └── patient.py
│       └── simulation/
│           ├── __init__.py
│           ├── macro/
│           │   ├── __init__.py
│           │   ├── agents.py
│           │   ├── config.py
│           │   ├── grid.py
│           │   └── simulator.py
│           └── micro/
│               ├── __init__.py
│               ├── engine.py
│               ├── genome.py
│               └── simulator.py
├── tests/
│   ├── __init__.py
│   ├── conftest.py
│   ├── test_run_coupled_simulation.py
│   └── integration_tests/
│       ├── __init__.py
│       ├── test_discharge_turnover.py
│       ├── test_grid.py
│       ├── test_macro_patient_integration.py
│       └── test_micro_patient_integration.py
├── .gitignore
├── .gitlab-ci.yml
├── .pre-commit-config.yaml
├── pyproject.toml
├── README.md
├── slurm_runner.sh
└── uv.lock
```

## Ordner-Übersicht

### `src/mss/`

Das einzige Source-Root für den gepflegten Anwendungscode.

### `src/mss/domain/`

Gemeinsame Domänenobjekte, die von Makro- und Mikro-Ebene genutzt werden.

- `patient.py`: kanonisches Patientenmodell, Enums, Behandlungszustand und der Austauschvertrag zwischen Makro und Mikro.

### `src/mss/simulation/macro/`

Logik der Spital-Netzwerk-Simulation.

- `config.py`: Konfigurations-Dataclass der Makro-Ebene.
- `agents.py`: Mesa-Agenten-Wrapper für die Spital-Gitter.
- `grid.py`: Abteilungsgitter und grobes Spital-Netzwerkgitter.
- `simulator.py`: Orchestrierung der Makro-Simulation, Aufnahmen, Verlegungen, Übertragung, Entlassung und Mikro-Kopplung.

### `src/mss/simulation/micro/`

Logik der bakteriellen Within-Host-Evolution.

- `genome.py`: Genomdarstellung, Resistenz-Traits und Fitness-Helfer.
- `engine.py`: Mikro-Konfiguration und die Stamm-Populations-Evolutions-Engine.
- `simulator.py`: Batch-Verarbeitungs-Schnittstelle und Verwaltung des Episoden-Lebenszyklus.
- `time_calibration.py`: Umrechnung schrittbezogener Mikro-Parameter auf neue reale Zeitschrittdefinitionen und kontrollierte Diagnose-Ensembles.

### `src/mss/cli/`

Ausführbare Einstiegspunkte, die die Anwendung aus tieferliegenden Modulen zusammensetzen.

- `run_coupled_simulation.py`: lädt die YAML-Konfiguration, führt die gekoppelte Makro/Mikro-Simulation aus und schreibt Parquet-Outputs. Stellt zudem `run_realistic_once()` bereit, das vom Sweep-Kalibrierungswerkzeug genutzt wird.
- `run_single_ward_calibration.py`: analytische β₀-Kalibrierung für eine Einzelstation. Leitet `base_transmission_rate` aus einer geschlossenen Formel ab und validiert sie per Simulation. `--n-runs > 1` führt ein stochastisches Ensemble über Seeds aus und aggregiert die Ergebnisse.
- `run_parameter_sweep.py`: strukturierte Parameter-Sweep-Kalibrierung. Variiert einen YAML-Parameter über ein definiertes Gitter, führt die Simulation für jeden Wert aus und plottet den Effekt auf eine Zielgrösse.
- `run_micro_time_calibration.py`: validiert die stündliche Mikro-Zeitskala als 12 einstündige Nachtfenster-Schritte und erzeugt Diagnose-Outputs bei alternativen Auflösungen.
- `visualize_results.py`: liest generierte Parquet-Outputs und schreibt Diagnose-Plots.

### `config/`

Laufzeit-Konfigurationsdateien. Diese sollen umgebungs- oder szenariospezifisch sein, nicht codespezifisch.

- `simulation_realistic.yml`: Haupt-Szenario der realistischen Simulation mit kalibrierten Parameterwerten.
- `simulation_abx.yml`: alternatives Szenario, abgestimmt auf antibiotika-fokussierte Läufe.
- `template.yml`: vollständig dokumentierte Referenzdatei, die jede unterstützte YAML-Variable mit Erklärungen auflistet. Zum Erstellen neuer Szenarien kopieren und anpassen.
- `01_Makro_Parameterübersicht.md`: Parameter-Referenztabelle für die Makro-Ebene, mit Typen (geschätzt / Kontextualisierungsparameter / kalibriert / Referenzwert / nicht identifizierbar), Quellen und Kalibrierergebnissen.
- `02_Mikro_Parameterübersicht.md`: Parameter-Referenz für die Mikro-Ebene (Within-Host-Evolution), inklusive Genmodell und Formeln.
- `calibration/`: eine Konfigurationsdatei pro Kalibrierung, in Ausführungsreihenfolge nummeriert.
  - `cal1_simulation_single_ward.yml`: Einzelstation-β₀-Kalibrierung (analytisch).
  - `cal2_proximity_decay.yml`: räumlicher Distanzabfall gegenüber dem Zimmernachbar-Übertragungsanteil (gleiche Zelle).
  - `cal3_isolation_effectiveness.yml`: Isolationswirksamkeit gegenüber der relativen Akquisitionsreduktion (Counterfactual-Baseline).
  - `cal4_micro_hourly_time_scale.yml`: Rezept für die Mikro-Zeitskalen-Ausrichtung auf 12 einstündige Schritte im Nachtfenster.

### `Organizational/`

Planungs- und Bewertungsdokumente auf Modulebene.

- `Mini_Challenge.md`: Aufgabenbeschreibung der Mini-Challenge-Komponente.
- `Modulbeschreibung.md`: Modulanforderungen, Lernziele und Bewertungskriterien.

### `docs/`

Projektdokumentation, Diagramme und analytische Materialien für Menschen. In dieser Reihenfolge lesen:

- `docs/00_Overview.md`: geführter Einstieg. Wie Makro, Mikro, Patient, Config und Outputs zusammenpassen ("ein Patient, ein Tag"), plus der Diagramm-Index. Hier starten.
- `docs/01_Makro_Overview.md`: Makro-Ebene (Spital-Netzwerk) mit Tagesablauf, Übertragung und Patient-Kopplung.
- `docs/02_Mikro_Overview.md`: Mikro-Ebene (Within-Host) mit Stamm-Populationen, Kopplung und Evolutionsmechanik.
- `docs/03_Modellverhalten_und_Methodik.md`: beobachtetes Modellverhalten (effektive gegenüber konfigurierten Raten, Stabilität) und die Methodik im Umgang mit nicht-plausiblen Ergebnissen.
- `docs/04_Micro_Time_Calibration.md`: Quellenanker, Formeln und Validierungsworkflow für stündliche Mikro-Schritte.
- `docs/system_overview/`: Mermaid-Diagramme, Systemkarten und das Hilfsskript zum Erzeugen der Graphen.
- `docs/organizational/`: reserviert für Planungs- oder Prozessdokumentation.

### `notebooks/`

Explorative Jupyter-Notebooks zur Inspektion generierter Outputs. Nicht Teil des importierbaren Pakets.

- `output_explorer.ipynb`: lädt den neuesten `outputs/<timestamp>`-Lauf und rendert Heatmaps mit Tages-Slider sowie Diagnoseansichten.

### `tests/`

Automatisierte Verifikation für die `src`-Struktur.

- `conftest.py`: stellt sicher, dass `src/` während der Testläufe im Importpfad liegt.
- `test_run_coupled_simulation.py`: Laden der Konfiguration und Aufbau der Startpopulation.
- `integration_tests/`: modulübergreifende Verhaltenstests für Makro, Mikro und Gitter-Interaktionen.

### `outputs/`

Generierte Simulationsartefakte. Das sind keine Quelldateien. Jeder Lauf erzeugt ein Unterverzeichnis mit Zeitstempel:

- `outputs/<timestamp>_<name>/data/`: Simulationsergebnis-Tabellen als Parquet-Dateien.
- `outputs/<timestamp>_<name>/plots/`: gerenderte Diagnose-Plots.

### `logs/`

Ausführungslogs, einschliesslich Slurm-Job-Ausgabe.

### `containers/`

Container-Definitionen und zugehörige Laufzeit-Assets.

## Namens- und Platzierungskonventionen

Diese Regeln für alle künftigen Ergänzungen verwenden:

- Allen Python-Anwendungscode unter `src/mss/` ablegen.
- Schichtübergreifende Geschäftsentitäten in `src/mss/domain/` ablegen.
- Simulations-Engines unter `src/mss/simulation/<layer>/` ablegen.
- Ausführbare Einstiegspunkte in `src/mss/cli/` ablegen.
- Szenario-YAML in `config/` ablegen, nie neben Code-Modulen.
- Diagramme, Architekturnotizen und generierte Karten in `docs/` ablegen.
- Generierte Daten nur in `outputs/` oder `logs/` ablegen, nie unter `src/` oder `tests/`.

Dateinamenskonventionen:

- `config.py` für reine Konfigurationsmodule verwenden.
- `simulator.py` für Orchestrierungsobjekte verwenden, die eine Ebene koordinieren.
- `engine.py` für Rechenkerne oder tieferliegende Simulationsmechanik verwenden.
- Singular-Namen für Domänenentitäten wie `patient.py` verwenden.
- Beschreibende Testdateinamen bevorzugen, die das geprüfte Verhalten widerspiegeln.

Importkonventionen:

- Aus `mss...` importieren, nicht aus relativen Top-Level-Ordnern.
- CLI-Module schlank halten und wiederverwendbare Logik nach `domain/` oder `simulation/` verschieben.
- Zirkuläre Abhängigkeiten vermeiden, indem `domain/` unabhängig von `cli/` bleibt.

## Projekt ausführen

Abhängigkeiten installieren:

```bash
uv sync
```

Gekoppelte Makro/Mikro-Simulation ausführen:

```bash
uv run mss-run --config config/simulation_realistic.yml
```

Plots aus vorhandenem Parquet-Output erzeugen (optional):

```bash
uv run mss-visualize --output-dir outputs/<timestamp>_<name>
```

Einzelstation-β₀-Kalibrierung ausführen (Kalibrierung 1):

```bash
uv run mss-calibrate --config config/calibration/cal1_simulation_single_ward.yml
# stochastisches Ensemble über viele Seeds:
uv run mss-calibrate --config config/calibration/cal1_simulation_single_ward.yml --n-runs 1000
```

Parameter-Sweep-Kalibrierung ausführen (Kalibrierungen 2 und 3):

```bash
uv run mss-sweep --sweep config/calibration/cal2_proximity_decay.yml
uv run mss-sweep --sweep config/calibration/cal3_isolation_effectiveness.yml
```

Mikro-Zeitskalen-Kalibrierung auf stündliche Schritte ausführen (Kalibrierung 4):

```bash
uv run mss-micro-time-calibrate --config config/simulation_realistic.yml --target-steps-per-day 12 --active-window-hours 12
```

Tests ausführen:

```bash
uv run pytest
```

### Reproduzierbarkeit und Seeds

Der Wert `run.seed` ist eine rein technische Kontrolle für die Reproduzierbarkeit, kein inhaltlicher Parameter zum Tunen. Ein gekoppelter Lauf (`mss-run`) ist bei gegebenem Seed deterministisch und nutzt einen einzelnen Seed. Die Robustheit gegenüber dem Zufall prüfen die Kalibrierungs- und Sweep-Werkzeuge, die Ensembles über viele Seeds ausführen (`mss-calibrate --n-runs N`, `mss-sweep` mit `n_seeds`) und Perzentil-/Konfidenzbänder berichten. Siehe `config/01_Makro_Parameterübersicht.md` und `docs/03_Modellverhalten_und_Methodik.md`.

### Grosse Outputs (git-lfs)

Das Verzeichnis `outputs/` wird über git-lfs getrackt und ist gross (mehrere GB an Lauf-Artefakten). Die Blobs sind optional herunterladbar. Um das Repository ohne sie zu klonen:

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone <repo-url>
```

Einzelne Läufe bei Bedarf nachladen:

```bash
git lfs pull --include="outputs/<timestamp>_<name>/**"
```

## Richtlinien zum Hinzufügen neuer Komponenten

Neue Funktionalität nach Verantwortung einordnen, nicht allein nach Feature-Namen.

### Eine neue Domänenentität hinzufügen

- In `src/mss/domain/` ablegen.
- Frei von CLI- und Dateisystem-Belangen halten.
- Aus `src/mss/domain/__init__.py` exportieren, falls Teil der Paket-API.

### Ein neues Makro-Verhalten hinzufügen

- Konfigurationsfelder kommen in `src/mss/simulation/macro/config.py`.
- Räumliche oder Topologie-Logik kommt in `src/mss/simulation/macro/grid.py`.
- Änderungen an der täglichen Orchestrierung kommen in `src/mss/simulation/macro/simulator.py`.

### Ein neues Mikro-Verhalten hinzufügen

- Genom- oder Trait-Berechnungen kommen in `src/mss/simulation/micro/genome.py`.
- Populations-Evolutionslogik kommt in `src/mss/simulation/micro/engine.py`.
- Batch-Ausführung, Persistenz oder Parallelität kommen in `src/mss/simulation/micro/simulator.py`.

### Einen neuen CLI-Befehl hinzufügen

- Ein neues Modul in `src/mss/cli/` erstellen.
- Einen Eintrag unter `[project.scripts]` in `pyproject.toml` hinzufügen.
- Das Argument-Parsing lokal im CLI-Modul halten und wiederverwendbare Logik aus tieferen Schichten importieren.

### Eine neue Szenario-Konfiguration hinzufügen

- Eine neue YAML-Datei unter `config/` anlegen.
- Nach dem Zweck des Szenarios benennen, nicht nach einem temporären Experiment.
- Aus der Dokumentation oder den Job-Runnern referenzieren, sobald sie ein unterstützter Workflow wird.

### Einen neuen Test hinzufügen

- Unit-artige Tests gehören in den oberen Bereich von `tests/`.
- Schicht- oder ablaufübergreifende Tests gehören in `tests/integration_tests/`.
- Testnamen am geprüften Modul oder Verhalten ausrichten.

## Skalierungshinweise

Diese Struktur ist auf kontrolliertes Wachstum ausgelegt:

- Neue Simulationsebenen können neben `macro/` und `micro/` unter `src/mss/simulation/` ergänzt werden.
- Weitere Domänenkonzepte können unter `src/mss/domain/` wachsen, ohne die Orchestrierung zu verunreinigen.
- CLI-Befehle können unabhängig wachsen, ohne dass Domänenmodule von Prozess-Belangen abhängig werden.
- Szenario-Wachstum bleibt auf `config/` beschränkt, statt Logik über Skripte zu duplizieren.

Wächst ein Unterpaket über etwa 5 bis 7 Dateien hinaus, nur dann ein fokussiertes Unterverzeichnis einführen, wenn es eine klarere Grenze schafft, zum Beispiel `src/mss/simulation/macro/policies/`.

## Wartungsempfehlungen

Um die Struktur langfristig gesund zu halten:

- Neue Top-Level-Code-Ordner ablehnen, ausser sie sind klar keine Source-Belange.
- Doppelte "Runner"-Logik in mehreren Dateien vermeiden.
- Mehrdeutige Module früh umbenennen, wenn sich ihre Verantwortung erweitert.
- Generierte Artefakte aus `src/`, `tests/` und `docs/` heraushalten.
- Dieses README aktualisieren, sobald ein neuer Top-Level-Ordner, ein neues Paket oder ein neuer öffentlicher Einstiegspunkt hinzukommt.
- Importe im Code-Review prüfen. `mss...` soll das Standard-Import-Root bleiben.

## Onboarding-Empfehlungen

Für neue Entwickelnde ist der schnellste Weg:

1. [`docs/00_Overview.md`](docs/00_Overview.md) für die geführte Tour und den dort empfohlenen Lesepfad lesen (Überblick -> Makro -> Mikro -> Modellverhalten).
2. Dieses README für die Repository-Karte und die Ausführungsbefehle überfliegen.
3. `src/mss/domain/patient.py` lesen, um den gemeinsamen Makro/Mikro-Vertrag zu verstehen, danach die beiden `simulator.py`-Dateien für die Orchestrierungsgrenzen.
4. `uv run pytest` ausführen, danach `uv run mss-run --config config/simulation_realistic.yml` und das erzeugte `outputs/`-Verzeichnis inspizieren.

Bei Onboarding-Reviews gilt eine Regel: Gehört eine Datei nicht unter `src/mss`, darf sie nur dann im Repository-Root liegen, wenn sie Konfiguration, Dokumentation, Automatisierung oder generierter Output ist.
