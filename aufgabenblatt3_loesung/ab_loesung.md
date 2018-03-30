## Aufgabe 1

* Lohnt es sich das Programm zu parallelisieren und wie lange braucht es dann?
* Skalierbarkeit abschätzen
* Kosten abschätzen

## Aufgabe 2

![Berechnung2](https://github.com/KingMus/hpc-aufgaben/blob/master/aufgabenblatt3_loesung/scr_bilder/AB3_Aufgabe2.png)

##### Benötigte Formeln:

> Speedup = T_para / T_seq <br> <br>
> Effizienz = Speedup / Anzahl Recheneinheiten <br>

##### Aussage zu Diagrammen:

> **Testreihe 1** zeigt einen sinkenden Effizienzanstieg umso mehr Recheneinheiten verwendet werden.
Dies kann an einem dadurch entstehenden Overhead erklärt werden. <br>
> **Testreihe 2** zeigt einen einen Effizienzrückgang nachdem eine gewisse Menge an Recheneinheiten überschritten wurde. <br>
> **Testreihe 3** zeigt ebenfalls einen sinkenden Effiziensanstieg, allerdings bewegen sich die Messungen näher am idealen Speedup. <br>

## Aufgabe 3

![Berechnung3](https://github.com/KingMus/hpc-aufgaben/blob/master/aufgabenblatt3_loesung/scr_bilder/AB3_Aufgabe3.png)

##### Benötigte Formeln:

> t_seq = T_seq / T_gesamt <br>
> t_para = T_para / T_gesamt <br> <br>
> t_para + t_seq = 1 <br> <br>
> Speedup = 1 / (t_seq + (1 / Anzahl Recheneinheiten) * t_para) <br>
> oder <br>
> Speedup = T_gesamt / (T_seq + (1 / Anzahl Recheneinheiten) * T_para) <br>

##### Aussage zu Diagrammen:

> **Testreihe 1** zeigt einen kleinen Anstieg, bleibt aber danach auf einem Level.
Ein weiteres Erhöhen der Recheneinheiten wird (durch das Diagramm abschätzbar) wahrscheinlich nichts mehr bewirken. <br>
> **Testreihe 2** zeigt dasselbe wie das Diagramm zu Testreihe 1. <br>

## Aufgabe 4 - Gustafsons Gesetz

![Berechnung4](https://github.com/KingMus/hpc-aufgaben/blob/master/aufgabenblatt3_loesung/scr_bilder/AB3_Aufgabe4.png)

##### Benötigte Formeln:

> t_seq = T_seq / T_gesamt <br>
> t_para = T_para / T_gesamt <br> <br>
> t_para + t_seq = 1 <br> <br>
> Speedup = t_seq + Anzahl Recheneinheiten * t_para <br>

##### Aussage zu Diagrammen:

> **Testreihe 1** zeigt einen kleinen Anstieg, bleibt aber danach auf einem Level.
Ein weiteres Erhöhen der Recheneinheiten wird (durch das Diagramm abschätzbar) wahrscheinlich nichts mehr bewirken. <br>
> **Testreihe 2** zeigt dasselbe wie das Diagramm zu Testreihe 1. <br>
