#define _GNU_SOURCE   // avverte che usiamo le estensioni GNU 
#include <stdio.h>    // permette di usare scanf printf etc ...
#include <math.h>     // per il valore assoluto in virgola mobile
#include <stdlib.h>   // conversioni stringa/numero exit() etc ...
#include <stdbool.h>  // gestisce tipo bool (per variabili booleane)
#include <assert.h>   // permette di usare la funzione assert
#include <string.h>   // funzioni di confronto/copia/etc di stringhe
#include <signal.h>
#include <unistd.h>
#include <errno.h>    // necessaria per usare errno
#include "xerrori.h"

#define QUI __LINE__,__FILE__
#define Buf_size 10

typedef struct{
    int *nodi;
    int num_de;
}dead;

typedef struct{
    int nodo;
    double pr;
}page_r;

typedef struct {
    int from;
    int to;
} Arco;

typedef struct {
    int *nodi_entranti;
    int  num_in;
} inmap;

typedef struct {
    int N;        // Numero dei nodi del grafo
    int *out;     // Array con il numero di archi uscenti da ogni nodo
    inmap *in;    // Array con gli insiemi di nodi entranti in ogni nodo
    dead de;       // Array con i nodi senza archi uscenti
    int N_archi;     // numero di archi validi
} grafo;
 
// per i thread produttori
typedef struct {
  Arco *buffer;                     // puntatore a un buffer condiviso
  int *ppindex;                     // puntatore per tenere traccia della posizione corrente nel buffer da parte del produttore
  pthread_mutex_t *pmutex_buf;      // puntatore a un mutex per sincronizzare l'accesso al buffer
  sem_t *sem_free_slots;            // puntatore a un semaforo per tracciare il numero di slot liberi nel buffer
  sem_t *sem_data_items;            // puntatore a un semaforo che indica il numero di elementi disponibili nel buffer 
  FILE *nomefile;                   // nome del file da cui i produttori leggono i dati da inserire nel buffer.
} dati_produttori;

// per i thread consumatori 
typedef struct {
  Arco *buffer;                     // puntatore a un buffer condiviso
  int *pcindex;                     // puntatore per tenere traccia della posizione corrente nel buffer da parte del consumatore
  pthread_mutex_t *pmutex_buf;      // puntatore a un mutex per sincronizzare l'accesso al buffer
  pthread_mutex_t *pmutex_grafo;    // puntatore ad un mutex per accedere al grafo
  pthread_mutex_t *pmutex_archi;    // puntatore ad un mutex per aggiornare il numero di archi validi
  sem_t *sem_free_slots;            // puntatore a un semaforo per tracciare il numero di slot liberi nel buffer
  sem_t *sem_data_items;            // puntatore a un semaforo che indica il numero di elementi disponibili nel buffer 
  grafo g;                          // grafo di output
} dati_consumatori;

typedef struct {
    double *x, *x_1;                // array al tempo t e t+1 
    double t1, *t2;                 // t2 e' un puntatore in quanto andro' a modificare il valore dopo la dichiarazione
    int *current_idx;               // Indice condiviso per il prossimo elemento da elaborare
    double *e;                      // ad e ogni thread gli andra' ad aggiungere l'errore di x nella posizione current_idx
    grafo *g;                       // grafo g che serve per il calcolo di t3
    double d;                       // serve per il calcolo di t3
    pthread_mutex_t *pmutex_e;      // puntatore ad un mutex per accedere alla var e
    pthread_mutex_t *pmutex_idx;    // puntatore ad un mutex per accedere alla var idx
    sem_t *sem_work_to_do;          // puntatore a un semaforo per tracciare il numero di slot liberi in x
    sem_t *sem_work_done;           // puntatore a un semaforo per tracciare il numero di elementi calcolati
    bool *continue_working;         // flag di controllo al ThreadData
    int *numiter;                   // numero d'iterazioni
} ThreadData;

typedef struct {
    int *exit;                      // puntatore a un intero che serve come flag di uscita. Se il valore puntato è 0, il thread continua a eseguire; se è 1, il thread termina.
    int *iter;                      // Un puntatore a un intero che tiene traccia del numero di iterazioni effettuate.
    double *x;                      // Un puntatore a un array di double che salva il contenuto dell'array aggiornato applicando il pagerank fino a quando viene chiamato il segnale
    int N;                          // grandezza vettore x
} ThreadSignal;

// funzione che libera la memoria allocata per il grafo
void free_grafo(grafo *g) {
    if (g->out) free(g->out);
    if (g->in) {
        for (int i = 0; i < g->N; i++) {
            free(g->in[i].nodi_entranti);
            
        }
        free(g->in);
    }
    if (g->de.nodi) free(g->de.nodi);
    
}

//funzione che inizializza il grafo
void inizializza_grafo(grafo *g, int num_nodi, int *num_archi) {
    g->N = num_nodi;
    g->out = calloc(num_nodi, sizeof(int)); //creo un array di dimensione n nodi
    if (g->out == NULL) {
    perror("Allocazione fallita per out");
    exit(EXIT_FAILURE);
    }            

    g->in = malloc(num_nodi * sizeof(inmap));
    if (g->in == NULL) {
    perror("Allocazione fallita per in");
    free(g->out);  // Liberare la memoria già allocata prima dell'errore
    exit(EXIT_FAILURE);
    }   
    for (int i = 0; i < num_nodi; i++) {
        g->in[i].nodi_entranti = NULL;  // Nessun arco entrante inizialmente
        g->in[i].num_in = 0;            // numero di nodi entranti (utile per far meno realloc possibili)
    }

    g->de.nodi = NULL;
    g->de.num_de = 0;
    g->N_archi = *num_archi;            // inserisco inizialmente il numero totale di archi, ogni qual volta vi è un arco non valido lo decremento
}

// funzione che fa un controllo se l'arco è valido
bool arco_valido(grafo *g, int from, int to) {

    if (from == to) return false;                           // return false se from = to
    inmap *dest = &g->in[to];
    for (int i = 0; i < dest->num_in; i++) {
        if (dest->nodi_entranti[i] == from) return false;   // return false se trova un arco duplicato
    }
    return true;                                            // return 1 se l'arco e' valido
}


// funzione utilizzata ogni qual volta si deve aggiungere un arco
void aggiungi_arco(grafo *g, int from, int to) {
    // Incrementa il numero di archi uscenti da 'from'
    g->out[from]++;

    // Aggiungi 'from' alla lista di nodi entranti per 'to'
    inmap *dest = &g->in[to];

    // alloco 10 locazioni di tipo intero così da ridurre i realloc per poi andare a fare una realloc finale per ridurre gli eventuali spazi in esubero
    if (dest->num_in % 10 == 0) {
        int *temp = realloc(dest->nodi_entranti, (dest->num_in + 10) * sizeof(int));
        if (temp == NULL) {
            perror("Memoria insufficiente per aggiungere l'arco");
            free(temp);
            return;  // Uscita precoce per evitare la perdita di puntatore
        }
        dest->nodi_entranti = temp;
    }

    dest->nodi_entranti[dest->num_in] = from;
    dest->num_in++;
    
}


FILE *grafo_start(const char *filename, grafo *g){
    FILE *file = fopen(filename, "r");                      // apro il file e lo inserisco nel puntatore "file"
    if (!file) {
        perror("Errore nell'apertura del file");            // esci se vi e' stato un errore nell'apertura del file
        exit(EXIT_FAILURE);
    }

    char line[1024];                                        // prendo l'intera riga

    int r, c, num_archi;                                    // nodo a, nodo b, numero di archi

     // Legge riga per riga fino a trovare una riga valida che non inizia con '%'
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '%') continue;  // Salta le righe di commento

        // Tentativo di leggere i numeri dalla riga
        if (sscanf(line, "%d %d %d", &r, &c, &num_archi) == 3) {
            break;  // Se ha letto correttamente i tre numeri, esce dal ciclo
        } else {
            fprintf(stderr, "Errore nella lettura del numero di nodi e archi\n");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    // controllo se il numero di nodi in entrata è diverso da quelli in uscita
    if (r != c) {
        fprintf(stderr, "Il numero di righe e colonne non corrisponde (r != c)\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    inizializza_grafo(g, r, &num_archi);            // chiamo la funzione di inizializzazione del grafo

    return file;                                    // ritorno il file per leggere gli archi 
}

// funzione che riempie l'array de.nodi con tutti quei nodi che non hanno archi uscenti
void DE(grafo *g) {
    for (int i = 0; i < g->N; i++) {
        if (g->out[i] == 0) {
            if (g->de.num_de % 10 == 0) {                                               // allocazione a spazi da 10 per diminuire il numero di realloc
                int *temp = realloc(g->de.nodi, (g->de.num_de + 10) * sizeof(int));
                if (temp == NULL) {
                    perror("Memoria insufficiente per aggiungere dead-end node");
                    return;  // Uscita precoce per evitare la perdita di puntatore
                }
                g->de.nodi = temp;
            } 
            g->de.nodi[g->de.num_de] = i;
            g->de.num_de++;
        }
    }

    int *temp = realloc(g->de.nodi, g->de.num_de * sizeof(int));                        // realloc per eliminare gli spazi vuoti allocati
    if (temp == NULL) {
        perror("Memoria insufficiente per ridimensionare dead-end nodes");
        return;  // Uscita precoce per evitare la perdita di puntatore
    }
    g->de.nodi = temp;
}

// funzione di stampa del grafo
void stampa_grafo(grafo *g){
     printf("Number of nodes: %d\n",g->N);
    printf("Number of dead-end nodes: %d\n",g->de.num_de);
    printf("Number of valid arcs: %d\n",g->N_archi);

}

// funzione eseguita dai thread producer
void *pbody(void *arg)
{  
  dati_produttori *a = (dati_produttori *)arg; // a = arg -> tipo produttori

  
  // apre il file e termina se non riesce
  FILE *f = a->nomefile;                                        // apro il file
  Arco n;                                                       // var utilizzata per prendere il valore dal file e poi essere inserita nel buffer
  do {
    int e = fscanf(f,"%d %d",&n.from, &n.to);
    if(e!=2) break;
    xsem_wait(a->sem_free_slots,QUI);                           // Decrementa il semaforo degli slot liberi
    xpthread_mutex_lock(a->pmutex_buf,QUI);                     // Blocca il mutex prima di accedere al buffer
    a->buffer[*(a->ppindex) % Buf_size].from = n.from -1;       // Inserisce il numero letto nel buffer alla posizione corrente dell'indice
    a->buffer[*(a->ppindex) % Buf_size].to = n.to -1; 
    *(a->ppindex) +=1;                                          // Incrementa l'indice del produttore
    xpthread_mutex_unlock(a->pmutex_buf,QUI);                   // Sblocca il mutex dopo aver accesso al buffer
    xsem_post(a->sem_data_items,QUI);                           // Incrementa il semaforo degli archi disponibili nel buffer
  } while(true);
  //puts("produttore sta per finire");
  fclose(f);
  return (void*)0;                            // Termina il thread
}

void *cbody(void *arg) {
    dati_consumatori *a = (dati_consumatori *)arg;

    Arco n;
    do {
        xsem_wait(a->sem_data_items,__LINE__,__FILE__);
        xpthread_mutex_lock(a->pmutex_buf,QUI);
        n = a->buffer[*(a->pcindex) % Buf_size];
        *(a->pcindex) +=1;
        xpthread_mutex_unlock(a->pmutex_buf,QUI);
        xsem_post(a->sem_free_slots,__LINE__,__FILE__);

        if(n.from == -1 && n.to == -1) break;

        xpthread_mutex_lock(a->pmutex_grafo,QUI); // Blocca il mutex del grafo prima di verificare e aggiungere l'arco
        if(arco_valido(&a->g, n.from, n.to)) {
            aggiungi_arco(&a->g, n.from, n.to);
        } else {
            xpthread_mutex_lock(a->pmutex_archi,QUI);
            a->g.N_archi -= 1;
            xpthread_mutex_unlock(a->pmutex_archi,QUI);
        }
        xpthread_mutex_unlock(a->pmutex_grafo,QUI);
    } while(true);
    return (void*)0;
}


double teleporting(double d, double N){                      
    return (1-d)/N;
}

double dead_end(grafo *g, double d, double *page_arr){      //eventuale ottimizzazione: al posto di chiedere tutto g chiedo solo il vettore e la dimensione
    double x = d/((double)(g->N));
    double somma = 0;
    for(int i = 0; i < g->de.num_de; i++){
        somma += page_arr[g->de.nodi[i]];
    }

    return x*somma;
}

double vet_y(double d, grafo *g, int j, double *page_arr){
    double somma = 0;
    for(int i = 0; i < g->in[j].num_in; i++){
        somma += ((page_arr[g->in[j].nodi_entranti[i]])/(g->out[g->in[j].nodi_entranti[i]]));
    }
    return (d*somma);
}

int comparePageRank(const void *a, const void *b) {
    const page_r *pa = (const page_r *)a;
    const page_r *pb = (const page_r *)b;

    if (pa->pr < pb->pr)
        return 1;  // Restituisce -1 se il primo elemento deve venire prima del secondo
    else if (pa->pr > pb->pr)
        return -1;   // Restituisce 1 se il primo elemento deve venire dopo il secondo
    return 0;       // Restituisce 0 se sono uguali
}

// funzione utilizzata da qsort per l'ordinamento dei nodi del pagerank
void stampa_page_rank(double *x, int n, int numiter, int k){            
    printf("Converged after %d iterations\n",numiter);
    if(k>n) k = n;

    double somma = 0;
    page_r *arr = malloc(n*sizeof(page_r));
    for(int i = 0; i < n; i++){
        arr[i].nodo = (i);
        arr[i].pr = x[i];
        somma += x[i];
    }
    printf("Sum of ranks: %f   (should be 1)\n",somma);
    qsort(arr, n, sizeof(page_r), comparePageRank);

    printf("Top %d nodes:\n",k);
    for(int i = 0; i < k; i++){
        printf("%d %f\n",arr[i].nodo,arr[i].pr);
    }
    free(arr);
}

void stampa_signal(double *x, int n, int numiter, int k){
    fprintf(stderr,"iterazione: %d\n",numiter);             // stampa il numero d'iterazione corrente

    double somma = 0;
    page_r *arr = malloc(n*sizeof(page_r));                 // creo un array che possa contenere tutti i nodi di tipo page_r
    for(int i = 0; i < n; i++){
        arr[i].nodo = i;                                    // inserisco il nodo e il relativo pagerank
        arr[i].pr = x[i];
        somma += x[i];                                      // somma di tutti i pagerank
    }
    qsort(arr, n, sizeof(page_r), comparePageRank);         // funzione di ordinamento chiamata su arr

    // stampo i migliori k nodi
    fprintf(stderr,"Top %d node:\n",k);
    for(int i = 0; i < k; i++){
        fprintf(stderr,"%d %f\n",arr[i].nodo,arr[i].pr);
    }
    free(arr);
}


void *fun_segnale(void *arg){
    ThreadSignal *d = (ThreadSignal*)arg;               // converto il puntatore arg in una struttura ThreadSignal
    sigset_t mask;                                      // variabile di tipo sigset_t utilizzata per rappresentare un insieme di segnali
    sigemptyset(&mask);                                 // inizializzo maschera segnali rimuovendo tutti i segnali da essa
    sigaddset(&mask, SIGUSR1);                          // aggiungo SIGUSR1 a maschera segnali, quindi unico segnale che attiverà la funzione
    
    int s;                                              // variabile che salva il tipo di segnale ricevuto
    while(*d->exit == 0){                               // continua finchè d->exit == 0
        int e = sigwait(&mask,&s);                      //  La funzione sigwait attende che uno dei segnali specificati nella maschera mask venga inviato al thread. Il segnale ricevuto viene memorizzato in s
        if(e != 0) termina("Errore: sigwait\n");
        if(s == SIGUSR1 &&*d->exit == 0){               // nel caso in cui il segnale è SIGUSR1 e d->exit == 0 allora chiama la funzione stampa_signal
            stampa_signal(d->x,d->N,*(d->iter),1);
        }
    }
    return (void*)0;                                    // sostituisce pthread_exit
}

void *aux_thread(void *arg){
    ThreadData *a = (ThreadData *)arg;
    int local_idx;

    while (true){

        sem_wait(a->sem_work_to_do);                    // Aspetta che ci sia lavoro da fare

        if (!(*(a->continue_working))) break;           // Uscire dal ciclo se non si deve più continuare a lavorare
    
        pthread_mutex_lock(a->pmutex_idx);              //mutex dell'indice
        if (*(a->current_idx) >= a->g->N) {             // se l'indice supera la grandezza dell'array
            pthread_mutex_unlock(a->pmutex_idx);        // sblocco l'indice
            sem_post(a->sem_work_done);                 // Segnala che questo thread ha finito
            continue;                                   // passo alla prossima iterazione e si mette in wait
        }

        local_idx = *(a->current_idx);                  // salvo in una var locale prima di incrementare
        (*(a->current_idx))++;                          // incremento current_idx
        pthread_mutex_unlock(a->pmutex_idx);            // sblocco current_idx

        double y = 0;
            y = vet_y(a->d, a->g, local_idx, a->x);     // calcolo t3

            a->x_1[local_idx] = a->t1 + *(a->t2) + y;               // in x_1 metto il valore calcolato attraverso il pagerank
            xpthread_mutex_lock(a->pmutex_e, QUI);                  // lock del mutex dell'errore
            *(a->e) += fabs(a->x_1[local_idx] - a->x[local_idx]);   // errore locale = differenza in valore assoluto tra x_1 e x
            xpthread_mutex_unlock(a->pmutex_e, QUI);                // unlock del mutex dell'errore
    }

    return (void*)0;                                                // sostituisce pthread_exit
}

double *pagerank(grafo *g, double d, double eps, int maxiter, int taux, int *numiter){
   
    double e = 0;                                                   // variabile che salverà l'errore ad ogni iterazione
    bool continue_working = true;                                   // booleana che impostata a true farà continuare i thread ausiliari, al contrario a false li fermerà

    // dichiarazione 2 vettori e primo elemento della formula
    double *x, *x_1, t1;          

    // alloco g->N spazi di dimesione double a x , x_1 , y 
    x = malloc((g->N)*sizeof(double));                              // array x
    if (x == NULL) {
    perror("Allocazione fallita per in");
    free(x);
    exit(EXIT_FAILURE);
    }

    x_1 = malloc((g->N)*sizeof(double));                            // array x_1
    if (x_1 == NULL) {
    perror("Allocazione fallita per in");
    free(x_1);
    free(x);
    exit(EXIT_FAILURE);
    }

    for(int i = 0; i < g->N; i++){                                  // riempio x al tempo 1
        x[i] = 1.0/g->N;
    }   

    t1 = teleporting(d,((double)(g->N)));                           // calcolo il primo elemento

    // Inizializzazione delle variabili di sincronizzazione
    pthread_mutex_t mue = PTHREAD_MUTEX_INITIALIZER;                // inizializzazione mutex errore
    pthread_mutex_t muidx = PTHREAD_MUTEX_INITIALIZER;              // inizializzazione mutex indice array
    sem_t sem_work_to_do, sem_work_done;                            // inizializzo i semafori impostandoli entrambi a 0
    xsem_init(&sem_work_to_do, 0, 0,__LINE__,__FILE__);
    xsem_init(&sem_work_done, 0, 0,__LINE__,__FILE__);

    ThreadData ac[taux];                                            // creo taux thread, taux -> passato dalla linea di comando
    pthread_t cons[taux];                                           // id dei thread ausiliari

    double t2_val = 0;                                              // variabile per il calcolo del secondo elemento della formula da elaborare ad ogni inizio del tempo t
    int current_idx = 0;                                            // indice dell'array x da passare ai thread
    int t_exit = 0;                                                 // t_exit impostato a 0 per far si che il thread segnale rimanga attivo
    
    ThreadSignal signal_thread_data = {                             // creazione thread segnale
    .N = g->N,
    .x = x,
    .exit = &t_exit,
    .iter = numiter
    };
    
    pthread_t signal_thread;                                        // id thread segnale
     
    sigset_t mask;                              // dichiarazione di mask di tipo sigset_t usata per rappresentare un insieme di segnali
    sigfillset(&mask);                          // inizializza la maschera dei segnali mask aggiungendo tutti i segnali disponibili.
    sigdelset(&mask, SIGQUIT);                  // Rimuove il segnale SIGQUIT dalla maschera.
    pthread_sigmask(SIG_BLOCK, &mask, NULL);    // Applica la maschera dei segnali al thread corrente. SIG_BLOCK indica che i segnali specificati nella mask devono essere bloccati. Quindi, tutti i segnali tranne SIGQUIT sono bloccati.

    // Crea il thread che gestirà SIGUSR1
    if (pthread_create(&signal_thread, NULL, fun_segnale, (void *)&signal_thread_data)) {
    perror("Failed to create the signal handler thread");
    exit(EXIT_FAILURE);
    }

    // Creazione dei thread
    for (int i = 0; i < taux; i++) {
        ac[i].x = x;
        ac[i].x_1 = x_1;
        ac[i].t1 = t1;
        ac[i].t2 = &t2_val;
        ac[i].d = d;
        ac[i].e = &e;
        ac[i].current_idx = &current_idx;
        ac[i].pmutex_e = &mue;
        ac[i].pmutex_idx = &muidx;
        ac[i].sem_work_to_do = &sem_work_to_do;
        ac[i].sem_work_done = &sem_work_done;
        ac[i].g = g;
        ac[i].continue_working = &continue_working;
        ac[i].numiter = numiter;
            
        if (pthread_create(&cons[i], NULL, aux_thread, (void *)&ac[i]) != 0) {
            perror("Failed to create thread");
            exit(EXIT_FAILURE);
        }
    }

    // Ciclo principale che controlla la convergenza e coordina i thread ed esce appena continue working e' false
    do{        
        
        current_idx = 0;
        e=0;
        t2_val = dead_end(g, d, x);               // calcolo il secondo elemento

        // Sblocca i thread per iniziare a lavorare
        for (int i = 0; i < ((g->N)+taux); i++) {       // incremento il semaforo sem_work_to_do g->n + taux affinchè tutti i thread riescano a rimettersi in uno stato di wait e a segnalare che non vi sono più elementi da calcolare
            sem_post(&sem_work_to_do);
        }

        // Attendi che tutti i thread completino il loro lavoro
        for (int i = 0; i < taux; i++) {
            sem_wait(&sem_work_done);
        }

        // Aggiornamento del vettore x con i nuovi valori calcolati
        memcpy(x, x_1, g->N * sizeof(double));

        if (e < eps || *numiter >= maxiter) {               // il valore di e viene modificato dai thread
        
            continue_working = false;                       // se e < eps o si superano le max iterazioni allo continue_working = false e si esce dal while 
            for (int i = 0; i < taux; i++) {
                sem_post(&sem_work_to_do);                  // si risbloccano i thread per farli terminare avendo continue_working = false
            }
            
        }
        (*numiter)++;                                       // incremento il numero d'iterazioni
    }while(continue_working);
    

     // Attesa della terminazione di tutti i thread
    for (int i = 0; i < taux; i++) {
        pthread_join(cons[i], NULL);
    }

    // kill del segnale
    t_exit = 1;
    pthread_kill(signal_thread, SIGUSR1);
    pthread_join(signal_thread,NULL);

    // Distruggi i mutex e le condizioni
    pthread_mutex_destroy(&mue);
    pthread_mutex_destroy(&muidx);

    // Distruggi i semafori
    sem_destroy(&sem_work_to_do);
    sem_destroy(&sem_work_done);

    // Libera le risorse allocate
    free(x_1);

    return x; // x viene restituito e gestito esternamente
}

int main(int argc, char *argv[]){
    int opt;
    int k = 3, m = 100, t = 3;
    double d = 0.9, e = 1.0e-7;
    

    while ((opt = getopt(argc, argv, "k:m:d:e:t:")) != -1) {                                // getopt cercherà negli argomenti passati al main le lettere -k, -m, ecc. Quando gli argomenti specificati finisco restituisce -1
        switch (opt) {
            case 'k':
                k = atoi(optarg);                                                           // optarg è una variabile gloabe al quale viene inserito il valore letto da getopt
                break;
            case 'm':
                m = atoi(optarg);
                break;
            case 'd':
                d = atof(optarg);
                break;
            case 'e':
                e = atof(optarg);
                break;
            case 't':
                t = atoi(optarg);
                break;
            default: 
                fprintf(stderr, "Usage: %s [-k K] [-m M] [-d D] [-e E] [-t T] infile\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // optind è l'indice utilizzato da getopt per scorrere gli argomenti opzionali di argv
    // se è >= di argc è un errore in quanto il programma si aspetta anche argomenti come il nome del file e il nome del programma
    if (optind >= argc) {                                                                      
        fprintf(stderr, "Expected infile after options\n");
        exit(EXIT_FAILURE);
    }

    int tc = t;                                                 // numero di thread
    assert(tc>0);                                               // thread devono essere positivi

    grafo g;                                                    // creo la variabile grafo
    FILE *f = grafo_start(argv[optind],&g);                     // inizializzo grafo, a grafo_start gli passo argv[optind] e la locazione di memoria di g, optind -> indice primo elemento non opzionale quindi il file di testo
   
    inmap *dest;
    for(int i = 0; i < g.N; i++){
        dest = &g.in[i];
        int *temp = realloc(dest->nodi_entranti, (dest->num_in) * sizeof(int));
                if (temp == NULL) {
                    perror("realloc fallita");
                    break;;  // Uscita precoce per evitare la perdita di puntatore
                }
                dest->nodi_entranti = temp;
    }

    // buffer produttori-consumatori
    Arco buffer[Buf_size];                                      // creo un beffer di tipo Arco di dimensione Buf_size inizializzata all'inizio del codice
    int pindex=0, cindex=0;                                     // creo gli indice per il produttore e il consumatore che scorrono il buffer
    pthread_mutex_t mucbuf = PTHREAD_MUTEX_INITIALIZER;         // creo il mutex per il buffer in modo da non creare concorrenza tra il produttore e tra i consumatori
    pthread_mutex_t mucgrafo = PTHREAD_MUTEX_INITIALIZER;       // creo il mutex per il grafo in modo da non creare concorrenza tra i consumatori
    pthread_mutex_t muagrafo = PTHREAD_MUTEX_INITIALIZER;       // creo il mutex per modificare il valore degli archi validi 
    sem_t sem_free_slots, sem_data_items;                       // creo i semafori per gestire il produttore e i consumatori
    xsem_init(&sem_free_slots,0,Buf_size,__LINE__,__FILE__);    // inizializzo il semaforo degli elementi liberi del buffer alla grandezza del buffer
    xsem_init(&sem_data_items,0,0,__LINE__,__FILE__);           // inizializzo a 0 il semaforo degli elementi disponibili nel buffer

    // dati per i thread
    dati_produttori ap;                                         // creo ap di tipo dati_produttori
    dati_consumatori ac[tc];                                    // creo un array di tipo dati_consumatori di garndezza tc ovvero t
    pthread_t prod;                                             // id thread produttore
    pthread_t cons[tc];                                         // id thread consumatori

    // associo le variabili ai campi del produttore
    ap.buffer = buffer;
    ap.ppindex = &pindex;
    ap.pmutex_buf = &mucbuf;
    ap.sem_free_slots = &sem_free_slots;
    ap.sem_data_items = &sem_data_items;
    ap.nomefile = f;

    // creo il produttore
    pthread_create(&prod, NULL, pbody, (void *)&ap);            // pthread_create(indice thread, attributi (NULL: default), funzione per ogni thread, dato utilizzato dal thread)

    // associo le variabili ai campi dei consumatori e per ogni iterazione creo un consumatore
    for (int i = 0; i < tc; i++) {
    ac[i].buffer = buffer;
    ac[i].pcindex = &cindex;
    ac[i].pmutex_buf = &mucbuf;
    ac[i].pmutex_grafo = &mucgrafo;
    ac[i].pmutex_archi = &muagrafo;
    ac[i].sem_free_slots = &sem_free_slots;
    ac[i].sem_data_items = &sem_data_items;
    ac[i].g = g;
    
    pthread_create(&cons[i], NULL, cbody, (void *)&ac[i]);
    }
    
    // attendo il produttore
    pthread_join(prod, NULL);                               // attendo che il thread prod termini, null -> non mi interessa il valore di ritorno
    for (int i = 0; i < tc; i++) {
    xsem_wait(&sem_free_slots, QUI);                        // attende che vi siano slot liberi
    xpthread_mutex_lock(&mucbuf, QUI);                      // lock del mutex buffer
    buffer[pindex % Buf_size] = (Arco){-1, -1};             // inserisce -1 , -1 in modo da comunicare ai consumatori che si possono arrestare 
    pindex++;                                               // incremento dell'indice
    xpthread_mutex_unlock(&mucbuf, QUI);                    // unlock del mutex buffer
    xsem_post(&sem_data_items, QUI);                        // incremento sem_data_items in modo da far leggere -1 , -1 ai consumatori
    }

    // comunico ai consumatori che possono terminare
    // attendo i consumatori 
    for (int i = 0; i < tc; i++) {
    pthread_join(cons[i], NULL);
    }

    // deallocazione, saluti, etc....
    pthread_mutex_destroy(&mucbuf);
    pthread_mutex_destroy(&mucgrafo);
    sem_destroy(&sem_free_slots);
    sem_destroy(&sem_data_items);

    DE(&g);                                                         // aggiungo i dead-end al grafo

    stampa_grafo(&g);                                               // Stampa il grafo

    int numiter = 0;                                                // numiter salva il numero di iterazioni 
    
    double *arr_page_rank = pagerank(&g, d, e , m, t, &numiter);    // chiamo la funzione pagerank

    stampa_page_rank(arr_page_rank,g.N,numiter,k);                  // stampo i primi k elementi dell'array arr_page_rank

    //pulizia memoria
    free_grafo(&g);
    free(arr_page_rank);

    return 0;
}
