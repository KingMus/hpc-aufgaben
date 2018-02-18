#Aufgabe 1

a) done
b) done
c) Die Reihenfolge der Threads wechselt
d) done

#Aufgabe 2

a) - Threadanzahl lässt sich nicht über nthreads ändern, einfügen von num_threads bei #prama omp parallel private löst das Problem
b) beide Threads üben einen lock auf jeweils a und b aus. dadurch können beide im jeweils nächsten schritt nicht mehr auf die jeweils andere variable zugreifen.
