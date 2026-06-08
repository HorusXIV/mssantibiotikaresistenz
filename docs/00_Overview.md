# System-Überblick: Einstiegspunkte & Zusammenhänge

Dieser Überblick zeigt die Einstiegspunkte und wie die Teile zusammenhängen; die Details
stehen in der jeweils verlinkten Datei.

**Worum geht es?** MRSA (methicillin-resistenter *Staphylococcus aureus*) ist ein
Spitalkeim, der gegen viele Antibiotika resistent ist: Er besiedelt Patienten, verbreitet
sich über Kontakte und erschwert Behandlungen. Das Modell ist **agentenbasiert**, weil
räumliche Nähe (Zimmer, Station) und die individuelle Vorgeschichte jedes Patienten
zählen. Ein "Schritt" der Mikro-Ebene ist eine numerische Diskretisierung des Tages
(`steps_per_day`), kein festes biologisches Intervall.

Die Simulation ist **mehrskalig** und über die `Patient`-Klasse gekoppelt:

- **Makro** (`src/mss/simulation/macro/`): Spital-Netzwerk, Übertragung zwischen
  Patienten, Aufnahme/Entlassung/Verlegung, Isolation, Antibiotikaverordnung. Zeitschritt: ein Tag.
- **Mikro** (`src/mss/simulation/micro/`): Within-Host-Evolution der Bakterienpopulation
  eines Carriers (Selektion, Mutation, HGT). Zeitschritt: `steps_per_day` Schritte pro Tag.
- **Patient** (`src/mss/domain/patient.py`): das Bindeglied. Makro schreibt täglich den
  Kontext auf den Patienten, der Patient erzeugt daraus einen Mikro-Request, und die
  Mikro-Antwort fliesst über den Patienten zurück in die Makro-Dynamik.

---

## Wo fange ich an? (Einstiegspunkte)

| Ich will verstehen … | Lies | Diagramm |
|---|---|---|
| Repository-Struktur, wie man es startet | `README.md` | – |
| Den Tagesablauf und den Makro↔Mikro-Datenfluss | `docs/01_Makro_Overview.md` | `docs/system_overview/Flowchart_v1.mmd` |
| Was die Within-Host-Evolution macht | `docs/02_Mikro_Overview.md` | `docs/system_overview/Pruned_MindMap.mmd` |
| Was jeder Parameter bedeutet, Typ und Quelle | `config/01_Makro_Parameterübersicht.md`, `config/02_Mikro_Parameterübersicht.md` | – |
| Alle YAML-Variablen als kopierbare Referenz | `config/template.yml` | – |
| Wie man Mechanismen ein-/ausschaltet (Szenarien) | `config/01_Makro_Parameterübersicht.md` (Abschnitt "Mechanismen-Steuerung") | – |
| Wie kalibriert wird | `config/01_Makro_Parameterübersicht.md` (Kalibrierungen) | Ensemble-/Sweep-Plots in `outputs/` |
| Modellverhalten, effektive Raten, Umgang mit Unplausiblem | `docs/03_Modellverhalten_und_Methodik.md` | – |
| Den fachlichen Problemraum (Kontext der Mini-Challenge) | – | `docs/system_overview/MindMap.mmd` |
| Die kausalen Wirkungszusammenhänge (Gephi) | `docs/system_overview/build_gephi_graphs.py` | `docs/system_overview/amr_system_map.gexf` |

Empfohlener Lesepfad für Neueinsteiger: `README.md` -> dieses Dokument ->
`docs/system_overview/Flowchart_v1.mmd` -> `src/mss/domain/patient.py` ->
`docs/01_Makro_Overview.md` und `docs/02_Mikro_Overview.md`.

---

## Ein Patient, ein Tag

So fliessen die Daten durch einen Trägerpatienten, von der Aufnahme bis zur Entlassung.
Die genaue tägliche Phasen-Reihenfolge der Engine steht in
[`docs/01_Makro_Overview.md`](01_Makro_Overview.md#tagesablauf-macrosimulatorstep) und im
Flowchart.

1. **Aufnahme**: Der Patient kommt ins Spital (Susceptible oder als Community-Carrier),
   bekommt eine Gitterposition (Zimmer) und eine geplante Verweildauer (LOS).
2. **Clearance** (nur Carrier): Mit Wahrscheinlichkeit `p_clearance` (aus dem Mikro)
   wird er spontan wieder Susceptible.
3. **Kontext & Erkennung**: Makro setzt Hygiene, mögliche Isolation (Erkennung mit
   `carrier_isolation_probability × base_diagnostic_speed`) und ein Antibiotika-Regime;
   das alles landet via `PatientDailyContext` auf dem Patienten.
4. **Mikro-Episode** (nur Carrier, nur wenn `use_micro: true`): `make_micro_request()`
   liefert ABX, Adherence und Immunstärke an die Mikro-Engine. Diese simuliert
   `steps_per_day` Schritte Selektion -> Mutation -> HGT -> Konsolidierung.
5. **Antwort zurück**: `apply_micro_response()` schreibt `resistant_fraction`,
   `dominant_genotype`, `relative_transmissibility`, `p_clearance` und
   `severity_modifier` auf den Patienten. Über `severity_modifier` koppelt der Mikro an
   die Makro-Mortalität und Aufenthaltsverlängerung.
6. **Übertragung**: Ist der Patient Carrier, trägt er zur Kolonisierungsgefährdung
   seiner Nachbarn bei (distanzgewichtet, reduziert durch Hygiene/Isolation). Ist er
   Susceptible, kann er mit `p_colonize` selbst Carrier werden und erbt den Stamm der
   Quelle.
7. **Verlegung/Entlassung/Mortalität**: Der Patient kann verlegt, entlassen werden oder
   versterben; danach beginnt der nächste Tag mit dem aktualisierten Zustand.

Übersicht in `docs/system_overview/Flowchart_v1.mmd` zu sehen.

---

## Diagramm-Legende (was zeigt was)

Alle Diagramme liegen in `docs/system_overview/`.

| Datei | Zeigt | Wann nützlich |
|---|---|---|
| `Flowchart_v1.mmd` | Detaillierter Tagesablauf mit allen Phasen und dem Makro↔Mikro-Datenfluss | Um die Ausführungsreihenfolge und Kopplung zu verstehen |
| `Flowchart_v0.mmd` | Frühe, kompakte Fassung des Tagesablaufs (Vorläufer von v1, weniger Phasen) | Für einen schnellen ersten Eindruck; `Flowchart_v1.mmd` ist die aktuelle, vollständige Version |
| `Pruned_MindMap.mmd` | Was tatsächlich im Modell steckt (Spital, Patient, 14 Gene, 12 Schritte) | Um den Modellumfang zu überblicken |
| `MindMap.mmd` | Der fachliche Problemraum der Antibiotikaresistenz (auch nicht-modellierte Faktoren) | Um den Kontext und die Abgrenzung der Mini-Challenge einzuordnen |
| `amr_system_map.gexf` / `*_nodes.csv` / `*_edges.csv` | Kausaler Wirkungsgraph mit gerichteten, signierten Kanten | Für die Analyse der Wirkungszusammenhänge in Gephi |
| `build_gephi_graphs.py` | Baut die Gephi-Dateien aus den CSVs (benötigt `networkx`) | Um den Wirkungsgraph neu zu erzeugen |
