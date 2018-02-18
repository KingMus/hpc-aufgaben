# Übung 1. zur Vorlesung Paralleles Rechnen

### Aufgabe: 1 Hello World (3 P)

a) done <br>
b) Mit #pragma parallel gelöst <br>
c) Die Reihenfolge der Threads wechselt <br>
d) Mit #pragma parallel sections gelöst <br>

### Aufgabe: 2 Errors (2 P)

a) Die Threadanzahl lässt sich nicht über *nthreads* ändern, einfügen von *num_threads* bei *#prama omp parallel private* löst das Problem. <br>
b) Beide Threads üben einen *lock* auf jeweils *a* und *b* aus. Dadurch können beide im jeweils nächsten Schritt nicht mehr auf die andere Variable zugreifen. <br>

### Aufgabe: 3 Kreiszahl π (5 P)

a) done <br>
b) mit *#pragma omp parallel for* über der for-Schleife gelöst <br>
c) *reduction ( +:globalCount)* eingefügt <br>
d) mit if-Anweisung gelöst (möglicherweise overengineered, nachfragen) <br>
e) mit eigener Variable *numthreads* und der Anweisung *#pragma num_threads(numthreads)* kann eine feste Anzahl an Threads vorgegeben werden. Unterbunden kann es werden, indem eine Konsoleneingabe gefordert wird die deaktiviert werden kann bzw. einen Standardvalue setzt

# Übung 2. zur Vorlesung Paralleles Rechnen

### Aufgabe: 1 Game of life (3 P)
### Aufgabe: 2 Domain Decomposition (6 P)
### Aufgabe: 3 Optional (5 P)

# Übung 3. zur Vorlesung Paralleles Rechnen

### Aufgabe: 1 Kennzahlen (1.5 P)
### Aufgabe: 2 Speedup (3.5 P)
### Aufgabe: 3 Amdahls Gesetz (2.5 P)
### Aufgabe: 4 Gustafsons Gesetz (2.5 P)

# Übung 4. zur Vorlesung Paralleles Rechnen

### Aufgabe: 1 Game of life (20 (+20) P)
