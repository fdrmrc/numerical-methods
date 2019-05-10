# numerical-methods
Raccolta esami ed esercizi per il corso di Metodi Numerici per le Equazioni Differenziali (aggiornata A.A. 2017/2018) 

Tutti i file sono scritti in Octave, perciò normalmente se lanciati in MatLab potrebbero dare degli errori di sintassi. Questi sono:

1. Gli incrementi del tipo a++, invece che con a=a+1

2. La chiusura di un ciclo for con endfor, invece che con end

3. La chiusura di un while con endwhile, invece che con end

4. BVP: ogni funzione da azzerare è scritta con una sintassi che in MatLab (versione R2018b e precedenti) probabilmente da errore. Infatti le condizioni ai bordi sono imposte direttamente nella definizione della funzione (vedi riga F=(u)). L'alternativa è quella di cambiare la prima e l'ultima riga delle matrici/termini noti in modo da forzare le condizioni ai bordi.
