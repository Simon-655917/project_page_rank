CC=gcc																# imposta gcc come compilatore
CFLAGS=-std=c11 -Wall -g -O3 -pthread								# flag per il compilatore
LDLIBS=-lm -lrt -pthread											# linka le librerie

# elenco degli eseguibili da creare
EXECS=pagerank														# nome dell'eseguibili da creare

all: $(EXECS)														# alla chiamata di make costruisce gli eseguibili da creare nella variabile EXECS

# regola per la creazione degli eseguibili utilizzando xerrori.o
$(EXECS): pagerank.o xerrori.o										# crea l'eseguibile pagerank collegando pagerank.o e xerrori.o

    #   - $(CC): compilatore utilizzato per compilare il file, 
    #   - $(LDFLAGS): contiene percorsi per speifiche librerie, 
    #   - -o $@: specifica il nome del file di output $@ prende il target ovvero pagerank e creerà un eseguibile di nome pagerank, 
    #   - $^: prende sia pagerank.o sia xerrori.o e garantisce che tutti i prerequisiti siano passati al linker, 
    #   - $(LDLIBS): contiene le librerie con cui il programma deve essere linkato
	
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)	
# regola per la creazione di file oggetto che dipendono da xerrori.h
%.o: %.c xerrori.h
	$(CC) $(CFLAGS) -c $<

# esempio di target che non corrisponde a una compilazione
# ma esegue la cancellazione dei file oggetto e degli eseguibili
clean:
	rm -f *.o $(EXECS) pagerank.o

# crea file zip della lezione
zip:
	zip threads.zip *.c *.h *.py makefile