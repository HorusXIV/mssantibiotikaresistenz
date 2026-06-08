# Modellverhalten und Methodik

Zwei Aspekte ergeben sich erst am fertigen Modell, nicht durch Vorab-Festlegung:

1. **Effektive vs. konfigurierte Raten**: Welche Grössen sind Eingaben, und welche
   ergeben sich erst emergent aus dem Zusammenspiel der Mechanismen?
2. **Methodik**: Wie gehen wir mit nicht-plausiblen Ergebnissen um und wie nutzen wir
   das Modell, um Einsichten in Mechanismen zu gewinnen, statt jeden Parameter einzeln
   zu rechtfertigen?

Die genannten Zahlen stammen aus echten Läufen (Single-Ward-Ensemble über die
cal1-Konfiguration sowie der realistische 365-Tage-Lauf in
`outputs/20260607_173453`), nicht aus Annahmen.

---

## Effektive vs. konfigurierte Raten

Ein wiederkehrendes Missverständnis ist, dass die in der YAML gesetzten Raten gleich
den beobachteten Raten sind. Tatsächlich sind viele beobachtete Grössen **emergent**:
Sie entstehen aus balancierenden Schleifen mehrerer Mechanismen und haben teils eine
andere Zeitkonstante als die Eingaben.

### Beispiel 1: Akquisitionsrate (konfiguriert ≈ effektiv, geschlossenes System)

Im Single-Ward-Setup (cal1, `use_micro = false`, geschlossene Kohorte) wurde `base_transmission_rate`
(β₀) analytisch so gewählt, dass die anfängliche Akquisitionsrate den Literaturanker
λ ≈ 0.0046–0.0054 pro Patient-Tag (4.6–5.4 pro 1000 Patient-Tage) trifft; gewählt:
β₀ = 0.07 -> λ(0) ≈ 0.0049 pro Tag.

Das Ensemble bestätigt diese effektive Rate: Über 40 Seeds entspricht die mittlere Zahl
der Erstinfektionen an Tag 1 dem theoretischen Erwartungswert λ(0)·S₀. Hier stimmen
konfigurierte und effektive Rate also **per Konstruktion** überein, was als
Konsistenz-Check dient.

### Beispiel 2: Prävalenz (emergent, nicht konfiguriert)

Im realistischen Netzwerk-Lauf (365 Tage, Mikro aktiv) ist die MRSA-Prävalenz **keine
Eingabe**. Gemessen (zweite Jahreshälfte): im Mittel **4.7 %** (sd 0.94 Prozentpunkte,
Spannweite 2.6–6.7 %), entsprechend ~44 Carriern bei ~930 Patienten. Dieser Wert ergibt
sich aus dem Gleichgewicht von

- Zustrom: Community-Import (1.7 % von ~126 Aufnahmen/Tag ≈ 2.1 Carrier/Tag) plus
  In-Hospital-Übertragung,
- Abstrom: Entlassung, Mortalität und spontane Dekolonisierung.

Die konfigurierte tägliche Hazard-Rate eines Susceptible (β₀ · c · (1−H) · I/N) liefert
bei der beobachteten Prävalenz nur ~0.011 pro Tag *vor* Isolation und stammspezifischer
Transmissibilität; die **effektive** Übertragung liegt darunter, weil 65 % der Carrier
erkannt und ihre Übertragung um 45 % gesenkt wird. Prävalenz ist damit ein effektives
Resultat, kein Stellhebel.

### Beispiel 3: Clearance-Zeitkonstante vs. Verweildauer (balancierende Schleife)

Die konfigurierte spontane Dekolonisierung `p_clearance` ≈ 0.0039 pro Tag entspricht
einer mittleren Tragezeit von ~255 Tagen. Die mittlere Verweildauer beträgt aber nur
~5 Tage (Ward) bzw. ~2.8 Tage (ICU). Folge: Innerhalb eines typischen Aufenthalts wird
ein Carrier mit weit grösserer Wahrscheinlichkeit **entlassen** als spontan geklärt. Die
*effektive* In-Hospital-Verweildauer als Carrier wird also durch LOS und
`carrier_extension_days` bestimmt, nicht durch die 255-Tage-Clearance. Die konfigurierte
Clearance-Rate wirkt erst über viele wiederholte Aufenthalte spürbar.

Genau das ist eine balancierende Schleife mit **anderer Zeitkonstante**: Die schnelle
Schleife (Aufnahme/Entlassung, Tage) dominiert die langsame Schleife (spontane
Clearance, Monate). Die beobachtete Prävalenz-Schwankung (sd ~0.9 Prozentpunkte) ist
daher überwiegend stochastisches Rauschen der schnellen Schleife und teils realistisch,
nicht ein Zeichen numerischer Instabilität.

---

## Methodik: Umgang mit nicht-plausiblen Ergebnissen und Einsichtsgewinn

Nicht jeder Parameter lässt sich einzeln gegen eine Literaturzahl rechtfertigen; viele
sind effektive Modellparameter (siehe Typen in `config/01_Makro_Parameterübersicht.md` und
`config/02_Mikro_Parameterübersicht.md`). Statt Einzelrechtfertigung nutzen
wir das Modell systematisch, um Mechanismen zu verstehen. Das Vorgehen:

1. **Ablation**: Mechanismen gezielt abschalten (siehe Tabelle "Mechanismen-Steuerung"
   in `config/01_Makro_Parameterübersicht.md`). `use_micro = false` isoliert die reine
   Makro-Dynamik; einzelne Raten auf 0 isolieren Mutation, HGT, Isolation, Transfer usw.
   So lässt sich zuordnen, welcher Mechanismus ein auffälliges Ergebnis verursacht.
2. **Ensemble statt Einzel-Seed**: Wo der Zufall geprüft wird, geschieht das über
   Ensembles vieler Seeds mit Perzentil-/Konfidenzbändern (Kalibrierung über `--n-runs`,
   Sweep über `n_seeds`), nicht über einen Einzellauf. Ein gekoppelter `mss-run` nutzt
   einen festen Seed; der Seed ist eine rein technische Kontrolle.
3. **Sweeps zur Konfundierungs-Diagnose**: Ein strukturierter Parameter-Sweep
   (`mss-sweep`) zeigt, ob eine Zielgrösse monoton und identifizierbar auf einen
   Parameter reagiert oder ob mehrere Parameter konfundiert sind.
4. **Effektive vs. konfigurierte Raten gegenprüfen**: Beobachtete Grössen werden mit den
   Eingaben verglichen (siehe oben), um emergente Effekte von Eingaben zu trennen.

### Fallbeispiel A: Schwankende Prävalenz ist kein Bug

Die Prävalenz schwankt im realistischen Lauf um ~0.9 Prozentpunkte. Vorgehen zur
Einordnung: (a) `use_micro = false` zeigt, dass die Schwankung schon ohne Mikro auftritt,
also aus der Makro-Aufnahme/Entlassung/Übertragung stammt; (b) der Vergleich der
Zeitkonstanten (Tage vs. Monate, siehe Beispiel 3) erklärt sie als balancierende
Schleife; (c) die Schwankung bleibt über den Lauf in einem schmalen Band und läuft nicht
weg (kein Aufschaukeln). Schlussfolgerung: stochastisches Gleichgewicht, teils
realistisch, kein numerischer Defekt.

### Fallbeispiel B: Konfundierung bei der Isolationswirkung

Das absolute Niveau der Isolationswirksamkeit ist mit β₀, Kontakten, Hygiene und
Raumstruktur konfundiert; eine direkte Kalibrierung des Absolutwerts wäre nicht
identifizierbar. Vorgehen: In Kalibrierung 3 wird stattdessen die **relative**
Akquisitionsreduktion (mit vs. ohne wirksame Isolation, pro Seed gepaart) als Zielachse
verwendet, in der sich die konfundierten Faktoren herauskürzen. So wird ein
identifizierbarer, effektiver Wert bestimmt, ohne einen nicht-identifizierbaren
Absolutwert zu erfinden. Nicht-identifizierbare Multiplikatoren werden konsequent als
Referenzwert (Standard 1.0) deklariert statt scheinpräzise gefittet.
