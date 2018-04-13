# Übung 1. zur Vorlesung Paralleles Rechnen

### Aufgabe: 1 Hello World (3 P)

a) done <br>
b) Mit #pragma parallel gelöst <br>
c) Die Reihenfolge der Threads wechselt <br>
d) Mit #pragma parallel sections gelöst <br>

### Aufgabe: 2 Errors (2 P)

a) Die Threadanzahl lässt sich nicht über *nthreads* ändern, einfügen von *num_threads* bei *#prama omp parallel private* löst das Problem. Außerdem verursacht die unterste Barriere möglicherweise einen Deadlock.<br>
b) Beide Threads üben einen *lock* auf jeweils *a* und *b* aus. Dadurch können beide im jeweils nächsten Schritt nicht mehr auf die andere Variable zugreifen. <br>

### Aufgabe: 3 Kreiszahl π (5 P)

a) done <br>
b) mit *#pragma omp parallel for* über der for-Schleife gelöst <br>
c) *reduction ( +:globalCount)* eingefügt <br>
d) mit if-Anweisung gelöst (möglicherweise overengineered, nachfragen) <br>
e) mit eigener Variable *numthreads* und der Anweisung *#pragma num_threads(numthreads)* kann eine feste Anzahl an Threads vorgegeben werden. Unterbunden kann es werden, indem eine Konsoleneingabe gefordert wird die deaktiviert werden kann bzw. einen Standardvalue setzt

# Übung 3. zur Vorlesung Paralleles Rechnen

Siehe [Übungsblatt3-Lösung](https://github.com/KingMus/hpc-aufgaben/blob/master/aufgabenblatt3_loesung/ab3_loesung.md)

