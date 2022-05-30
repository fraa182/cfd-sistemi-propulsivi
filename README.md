# CFD Sistemi Propulsivi (Larocca)
Repository del progetto di gruppo 

# Come funziona il tutto
- Ogni volta che facciamo esercitazione o che si modifica il codice, chi modifica il codice deve caricarlo nella repository facendo "add file" e poi "upload files", oppure direttamente modificando i file caricati copiando e incollando il nuovo codice modificato. In questo modo evitiamo che si incasini tutto e possiamo gestire gran parte delle cose da qua.
- Se si vogliono scaricare i file in locale per runnarli basta fare "code" e poi "download zip". Una volta scaricati basta aprirli nel progetto con codeblocks o altri IDE.

**Come è organizzata la repository**
- Nella cartella "images" ci vanno le immagini che generiamo con VisIt e che vanno eventualmente nella relazione (disegni e schemi che non c'entrano proprio coi grafici che vuole Larocca basta che li mettiamo nell'overleaf)
- Nella cartella main invece ci vanno semplicemente i vari codici sorgente di subroutine, moduli, funzioni e main.

# Come compilare
- Scaricare la cartella e aprirla nel terminale
- Scrivere "gfortran *.f90 *.FOR" e premere invio, poi scrivere "a.exe" se su Windows, oppure "./a.out" se su UNIX (Mac o Linux)
- La prima volta darà errore su "variabili.mod", ricompilare una seconda volta e va
