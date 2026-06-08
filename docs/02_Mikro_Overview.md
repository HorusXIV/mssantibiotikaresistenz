# Mikro-Ebene (Within-Host-Evolution)

Die Mikro-Ebene modelliert, was **innerhalb eines einzelnen Trägerpatienten** mit der
bakteriellen Population passiert: Wachstum, Tod, Konkurrenz zwischen Stämmen, Mutation,
horizontaler Gentransfer (HGT) und Selektion durch Antibiotika und Immunsystem. Daraus
leitet sie die makro-relevanten Grössen `p_clearance`, `relative_transmissibility` und
`severity_modifier` ab.

Code: [`genome.py`](../src/mss/simulation/micro/genome.py) (Genom, Fitness, Resistenz),
[`engine.py`](../src/mss/simulation/micro/engine.py) (Within-Host-Dynamik),
[`simulator.py`](../src/mss/simulation/micro/simulator.py) (Episodenverwaltung, Batch).

> Alle Parameter (mit Typ, Wert, Einheit), die vollständige 14-Gen-Tabelle, die
> ABX-Profile und die ausgeschriebenen Formeln stehen in der Referenz
> [`config/02_Mikro_Parameterübersicht.md`](../config/02_Mikro_Parameterübersicht.md).

## Grundidee

Modelliert wird nicht ein einzelnes Bakterium, sondern eine **Population aus Stämmen**.
Jeder Stamm hat ein Genom als Vektor aus `NUM_GENES = 14` kontinuierlichen Genwerten,
eine Populationsgrösse, ein Linienalter und eine angesammelte Schadenslast.

Ein Makro-Tag wird in `steps_per_day` Schritte zerlegt (Standard 12; ein Schritt ist
eine numerische Diskretisierung, kein festes biologisches Intervall, und sichert
numerische Stabilität bei schnellen Populationsdynamiken). Pro Schritt läuft in dieser
Reihenfolge:

1. **Selektion**: Wachstum und Tod aus Fitness, Immundruck, Antibiotika und
   Lebenszyklus-Kosten.
2. **Mutation**: neue Varianten durch verrauschte Genveränderungen.
3. **HGT** (jeder dritte Schritt): Gene wandern zwischen Stämmen.
4. **Konsolidierung**: sehr kleine Stämme werden entfernt, zu viele abgeschnitten.

## Architektur und Datenfluss

Die Mikro-Ebene beginnt **nicht jeden Tag bei null**: Der Episodenzustand bleibt über
Tage erhalten, sodass sich dominante Stämme etablieren und Alterung/Schaden akkumulieren.

```
Makro (_build_context) ── PatientDailyContext ──▶ Patient.update_context
Patient.make_micro_request ── request ──▶ MicroSimulator
   EpisodeState (Vortag) ──▶ simulate_day (steps_per_day Schritte)
MicroSimulator ── response ──▶ Patient.apply_micro_response ──▶ Makro
```

Die vier Stationen:

1. **Makro baut den Tageskontext** (`_build_context()`): Station, Isolation,
   Antibiotika-Regime.
2. **`patient.py` übersetzt** den Makro-Zustand in einen Request. Ein Request entsteht
   nur für `CARRIER` mit `episode_id`. `update_context()` leitet u.a. `adherence` ab.
3. **`MicroSimulator` rechnet**: lädt den `EpisodeState` (oder initialisiert eine neue
   Population), ruft `simulate_day()`, gibt die Response zurück und persistiert den
   Zustand bis zum nächsten Tag.
4. **Rückwirkung auf Makro** über `apply_micro_response()`.

Der `EpisodeState` speichert pro Episode `episode_id`, `patient_id`, `day` und die
Population (`genomes`, `populations`, `lineage_ages`, `damage_loads`, `strain_names`).

### Was von Makro/Patient einfliesst

| Eingang | Wirkung in Mikro |
|---|---|
| `state = CARRIER` + `episode_id` | Voraussetzung; nur dann läuft Mikro |
| `abx.on/class/dose_level` | wählt das ABX-Profil und den Selektionsdruck |
| `adherence` (aus `compliance`, ICU +0.1) | skaliert die effektive Kill-Rate |
| `immune_strength` | skaliert Immun-Clearance und spätere `p_clearance` |
| `initial_state.resistant_fraction` | Anteil resistenter Startpopulation (neue Episode) |
| `initial_state.dominant_genotype` | Klasse der resistenten Seed-Stämme (S/R1/R2/R3) |

Hygiene, Isolation und Diagnostik wirken ausschliesslich auf der Makro-Ebene.

### Was an Makro zurückgeht

| Ausgang | Wirkung in Makro |
|---|---|
| `p_clearance` | Entscheidung `CARRIER → SUSCEPTIBLE` (`should_clear_today()`) |
| `relative_transmissibility` | skaliert den Transmissionsbeitrag des Trägers |
| `severity_modifier` | skaliert Aufenthaltsverlängerung und Mortalität |
| `dominant_genotype`, `resistant_fraction` | werden bei Übertragung an neue Träger vererbt |

Bei Übertragung wird der Quellstamm exakt vererbt; Resistenz-Drift entsteht erst
within-host (Makro wendet keine Mutation an).

## Genommodell (konzeptionell)

Es gibt **keine diskreten Allele** im mendelschen Sinn. Jedes "Gen" ist ein
**kontinuierlicher Wert in `[0, 1]`** (0 = kaum ausgeprägt, 1 = stark ausgeprägt).
Mutation verändert diesen Wert durch Gauss-Rauschen, HGT mischt Werte zwischen Stämmen.

Die 14 Gene decken Wachstum, Resistenzmechanismen (Efflux, Zielmodifikation,
Permeabilität), Virulenz/Adhäsion/Stealth, Evolvierbarkeit (Mutationsrate, HGT-Kompetenz)
und Lebenszyklus (Reparatur, Dormanz, Stress, Toleranz) ab. Die vollständige Tabelle mit
Indizes und Wirkung steht in der [Mikro-Parameterreferenz](../config/02_Mikro_Parameterübersicht.md).

Zwei Begriffe, die oft verwechselt werden:

- **`dominant_genotype` (S/R1/R2/R3)** ist kein eigenes Gen, sondern eine **abgeleitete
  Klasse** aus dem Mittel der drei Resistenzgene (Efflux, Zielmodifikation,
  Permeabilität) über feste Schwellen.
- **`resistant_fraction`** ist der **Populationsanteil** mit Resistenzscore `>= 0.3`. Ein
  Patient kann `dominant_genotype = S` haben und trotzdem einen kleinen resistenten
  Subanteil tragen.

Bei einer neuen Episode erzeugt `StrainPopulation.create_initial()` aus
`resistant_fraction` (wie viel resistent) und `dominant_genotype` (welche Resistenzklasse)
eine leicht verrauschte Startpopulation (Default 3 sensible + 2 resistente Seed-Stämme).

## Mechanik pro Schritt (konzeptionell)

Die genauen Formeln stehen in der [Mikro-Parameterreferenz](../config/02_Mikro_Parameterübersicht.md);
hier die Wirkprinzipien.

**Selektion.** Fitness kombiniert Basiswachstum, Resistenzkosten (durch
`METABOLIC_OPTIMIZATION` gedämpft), Antibiotika-Überleben (Kill-Rate × Dosis ×
`adherence`, gemindert durch passende Resistenzgene) und Immun-Überleben (gemindert durch
`STEALTH`). Entscheidend ist die Fitness **relativ zum Populationsmittel**, potenziert mit
`selection_strength`: fitte Stämme wachsen schneller, schwache schrumpfen schneller.
Dormanz, Schadensaufbau/-reparatur und Linienalter wirken als Lebenszyklus-Kosten auf
Wachstum und Tod. Sehr kleine Stämme unterliegen demografischer Stochastik (Poisson) und
können zufällig aussterben; oberhalb der `carrying_capacity` wird proportional
herunterskaliert.

**Mutation.** Die effektive Rate steigt unter Antibiotikastress
(`stress_mutation_boost`) und mit dem Gen `MUTATION_RATE_MODIFIER`. Pro Ereignis wird ein
kleiner Teil des Stamms zu einem neuen, leicht veränderten Stamm abgespalten.

**HGT.** Mit einer von `HGT_COMPETENCE` abhängigen Wahrscheinlichkeit übernimmt ein Stamm
Gene von einem populationsgewichtet gezogenen Donor. Übertragbar sind vor allem
Resistenz- und Persistenzgene; der Donorwert wird eingemischt, nicht hart ersetzt.

**Clearance.** `compute_clearance_probability()` liefert die tägliche
`CARRIER → SUSCEPTIBLE`-Wahrscheinlichkeit: kleine Population und starkes Immunsystem
erhöhen sie, hohe `STEALTH`-Werte senken sie.

## Ein typischer Tag eines Trägers

1. Makro erkennt den Patienten als `CARRIER` und baut den `PatientDailyContext`.
2. `Patient.update_context()` setzt Kontext, Regime und `adherence`.
3. `Patient.make_micro_request()` baut den Tagesrequest.
4. `MicroSimulator` lädt den `EpisodeState` oder erzeugt eine neue Population.
5. `simulate_day()` führt `steps_per_day` Schritte aus (Selektion, Mutation, HGT,
   Konsolidierung).
6. `population_to_response()` bestimmt dominanten Stamm und Makro-Effekte.
7. `Patient.apply_micro_response()` schreibt die Werte zurück; Makro nutzt sie für
   Clearance, Transmission und Mortalität.

## Merksätze

- Makro fragt: **In welchem Patienten und unter welchen Rahmenbedingungen?** Mikro
  antwortet: **Welcher Stamm setzt sich durch und mit welchen Folgen?**
- Das Genom ist ein **Trait-Vektor**, kein Sequenzmodell.
- `dominant_genotype` ist eine **abgeleitete Resistenzklasse**, kein eigenes Gen.
- `resistant_fraction` beschreibt die **gesamte Population**, nicht nur den dominanten
  Stamm.
- Der Episodenzustand **persistiert über Tage**; das macht echte Within-Host-Evolution
  erst möglich.
