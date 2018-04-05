# hpc-aufgaben
Aufgaben der Vorlesung High Performance Computing 2018

<hr>

#### Kompilieren

Build in Eclipse. <br>
Run-Befehl aus Folien im Terminal im Ordnerverzeichnis.

#### Branches

In den Branches finden sich verschiedene Lösungsansätze, die mehr oder weniger weit ausgebaut sind.
Der Name gibt in etwa vor, bei welcher Aufgabe der Branch Relevanz hat.

###### gol_parallelfor

AB2 umgesetzt mit omp_parallel_for. Entspricht nicht unbedingt dem gewünschten Ergebnis,
da keine eigenen Dateien für jeden Thread geschrieben werden, funktioniert aber an sich sehr gut.

###### gol_sectionsplitting

AB2 umgesetzt mit omp_section. Je nach Auslegung des Arbeitsblattes eine mögliche Lösung.
Keine Möglichkeit der flexiblen Threadanzahl (immer 4), dafür aber mit Ausgabe.
Branch geschlossen, findet sich im Master-Gol.

###### gol_mpi_comments

Die für AB4 umgesetzte Aufgabe ist mit Comments versehen, die den Code erklären. Zumindest versuchen sie es.
Wenn auch schlecht. Egal. Im Code sind außerdem regelmäßig Code-Stücke, die Infos ausgeben an denen der Ablauf der Anwendung
nachvollzogen werden kann so gut es geht. Im Master befindet sich der indentische Code in "schlank".
