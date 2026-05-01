#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <float.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

static inline double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}

typedef struct {
    long long osm_id;   /* identifiant OSM original (grand entier) */
    double    lat;      /* latitude  en degrés décimaux            */
    double    lon;      /* longitude en degrés décimaux            */
} Node;

/* CSRGraph : le graphe stocké en format Compressed Sparse Row */
typedef struct {
    int     n_nodes;
    int     n_edges;
    int    *row_ptr;   /* tableau des débuts d'adjacence, taille n_nodes+1 */
    int    *adj;       /* tableau des destinations, taille n_edges          */
    double *weights;   /* tableau des poids en mètres, taille n_edges       */
} CSRGraph;

/*
 * HeapNode : un élément de notre file de priorité (tas min).
 * On stocke le nœud et sa distance courante.
 * La priorité = f_score (dist pour Dijkstra, dist+heuristique pour A*).
 */
typedef struct {
    int    node;     /* indice du nœud 0..N-1 */
    double f_score;  /* clé de priorité (plus petit = plus prioritaire) */
} HeapNode;

/*
 * MinHeap : file de priorité min implémentée comme un tas binaire.
 * Un tas binaire est un tableau où parent(i) = (i-1)/2,
 * enfant_gauche(i) = 2*i+1, enfant_droit(i) = 2*i+2.
 * La règle du tas : parent.f_score <= enfants.f_score (c'est un MIN-tas).
 */
typedef struct {
    HeapNode *data;   /* tableau des éléments */
    int       size;   /* nombre d'éléments actuels */
    int       cap;    /* capacité allouée */
} MinHeap;

/*
 * ShortestPathResult : résultat renvoyé par Dijkstra et A*.
 * Contient la distance, le chemin, et des stats de debug.
 */
typedef struct {
    double  dist;           /* distance totale en mètres (-1 si inaccessible) */
    int    *path;           /* tableau des nœuds du chemin (source..target)   */
    int     path_len;       /* nombre de nœuds dans le chemin                 */
    int     nodes_explored; /* nombre de nœuds extraits du tas (pour benchmark) */
    int     relaxations;    /* nombre de fois qu'on a mis à jour une distance  */
    double  time_ms;        /* temps de calcul en millisecondes                */
} ShortestPathResult;


/* Table de hashage*/

#define HASH_EMPTY  (-1LL)

typedef struct {
    long long key;   /* osm_id, HASH_EMPTY si case vide */
    int       val;   /* indice local 0..N-1             */
} HashEntry;

typedef struct {
    HashEntry *table;
    int        cap;    /* taille du tableau, toujours une puissance de 2 */
    int        size;
} HashMap;

static HashMap *map_create(int capacity) {
    int cap = 1;
    while (cap < capacity * 2) cap <<= 1;  /* puissance de 2 >= 2*capacity */
    HashMap *m = malloc(sizeof(HashMap));
    m->table   = malloc(cap * sizeof(HashEntry));
    m->cap     = cap;
    m->size    = 0;
    for (int i = 0; i < cap; i++) m->table[i].key = HASH_EMPTY;
    return m;
}

/* Calcule l'indice de départ dans le tableau pour une clé donnée. */
static inline int map_slot(HashMap *m, long long key) {
    unsigned long long h = (unsigned long long)key * 2654435761ULL;
    return (int)(h & (unsigned)(m->cap - 1));
}

/* Insère (key, val) dans la table. Utilise le sondage linéaire en cas de collision. */
static void map_insert(HashMap *m, long long key, int val) {
    int i = map_slot(m, key);
    while (m->table[i].key != HASH_EMPTY && m->table[i].key != key)
        i = (i + 1) & (m->cap - 1);
    m->table[i].key = key;
    m->table[i].val = val;
    m->size++;
}

/* Recherche key dans la table. Retourne val ou -1 si absent. */
static int map_get(HashMap *m, long long key) {
    int i = map_slot(m, key);
    while (m->table[i].key != HASH_EMPTY) {
        if (m->table[i].key == key) return m->table[i].val;
        i = (i + 1) & (m->cap - 1);
    }
    return -1;
}

static void map_free(HashMap *m) { free(m->table); free(m); }

static Node *read_nodes(const char *path, int *out_n, HashMap **out_map) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    int count = 0;
    char line[256];
    if (!fgets(line, sizeof(line), f)) { fprintf(stderr, "fichier vide\n"); exit(1); }
    while (fgets(line, sizeof(line), f)) count++;
    rewind(f);
    if (!fgets(line, sizeof(line), f)) { fprintf(stderr, "fichier vide\n"); exit(1); }

    Node    *nodes = malloc(count * sizeof(Node));
    HashMap *map   = map_create(count);
    int idx = 0;

    while (fgets(line, sizeof(line), f)) {
        long long id; double lat, lon;
        if (sscanf(line, "%lld,%lf,%lf", &id, &lat, &lon) != 3) continue;
        nodes[idx].osm_id = id;
        nodes[idx].lat    = lat;
        nodes[idx].lon    = lon;
        map_insert(map, id, idx);
        idx++;
    }
    fclose(f);
    *out_n = idx; *out_map = map;
    return nodes;
}

typedef struct { int u, v; double w; } RawEdge;

static RawEdge *read_edges(const char *path, HashMap *map, int *out_m) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    int count = 0;
    char line[256];
    if (!fgets(line, sizeof(line), f)) { fprintf(stderr, "fichier vide\n"); exit(1); }
    while (fgets(line, sizeof(line), f)) count++;
    rewind(f);
    if (!fgets(line, sizeof(line), f)) { fprintf(stderr, "fichier vide\n"); exit(1); }

    RawEdge *edges = malloc(count * sizeof(RawEdge));
    int idx = 0, skipped = 0;

    while (fgets(line, sizeof(line), f)) {
        long long u_osm, v_osm; double w;
        if (sscanf(line, "%lld,%lld,%lf", &u_osm, &v_osm, &w) != 3) continue;
        int u = map_get(map, u_osm);
        int v = map_get(map, v_osm);
        if (u < 0 || v < 0) { skipped++; continue; }
        edges[idx].u = u; edges[idx].v = v; edges[idx].w = w;
        idx++;
    }
    fclose(f);
    if (skipped)
        fprintf(stderr, "  [warn] %d arêtes ignorées\n", skipped);
    *out_m = idx;
    return edges;
}

static CSRGraph *build_csr(int n, RawEdge *edges, int m) {
    CSRGraph *g = malloc(sizeof(CSRGraph));
    g->n_nodes = n; g->n_edges = m;
    g->row_ptr = calloc(n + 1, sizeof(int));
    g->adj     = malloc(m * sizeof(int));
    g->weights = malloc(m * sizeof(double));

    /* Passe 1 */
    for (int i = 0; i < m; i++) g->row_ptr[edges[i].u + 1]++;
    /* Passe 2 */
    for (int i = 1; i <= n; i++) g->row_ptr[i] += g->row_ptr[i - 1];
    /* Passe 3 */
    int *cursor = malloc(n * sizeof(int));
    memcpy(cursor, g->row_ptr, n * sizeof(int));
    for (int i = 0; i < m; i++) {
        int u = edges[i].u, pos = cursor[u]++;
        g->adj[pos] = edges[i].v; g->weights[pos] = edges[i].w;
    }
    free(cursor);
    return g;
}

static void csr_free(CSRGraph *g) {
    free(g->row_ptr); free(g->adj); free(g->weights); free(g);
}


/* MIN-HEAP */

static MinHeap *heap_create(int initial_cap) {
    MinHeap *h = malloc(sizeof(MinHeap));
    h->data    = malloc(initial_cap * sizeof(HeapNode));
    h->size    = 0;
    h->cap     = initial_cap;
    return h;
}

static void heap_free(MinHeap *h) { free(h->data); free(h); }

/* heap_push : ajoute (node, f_score) dans le tas.
 * On place le nouvel élément à la fin du tableau (indice size),
 * puis on le "fait remonter" (sift-up) tant qu'il est plus petit que son parent.
 */
static void heap_push(MinHeap *h, int node, double f_score) {
    /* Agrandir si nécessaire */
    if (h->size == h->cap) {
        h->cap *= 2;
        h->data = realloc(h->data, h->cap * sizeof(HeapNode));
    }
    /* Placer à la fin */
    int i = h->size++;
    h->data[i].node    = node;
    h->data[i].f_score = f_score;

    /* Sift-up : remonter tant que plus petit que le parent */
    while (i > 0) {
        int parent = (i - 1) / 2;
        if (h->data[parent].f_score <= h->data[i].f_score) break;
        /* Échange avec le parent */
        HeapNode tmp   = h->data[parent];
        h->data[parent] = h->data[i];
        h->data[i]      = tmp;
        i = parent;
    }
}

/* heap_pop : retire et retourne l'élément de priorité minimale (la racine).
 * On remplace la racine par le dernier élément,
 * puis on "fait descendre" (sift-down) pour restaurer la propriété du tas.
 */
static HeapNode heap_pop(MinHeap *h) {
    HeapNode result = h->data[0];           /* sauvegarde du minimum */
    h->data[0] = h->data[--h->size];        /* dernier élément → racine */

    /* Sift-down : descendre tant qu'un enfant est plus petit */
    int i = 0;
    for (;;) {
        int left  = 2 * i + 1;
        int right = 2 * i + 2;
        int smallest = i;

        if (left  < h->size && h->data[left].f_score  < h->data[smallest].f_score)
            smallest = left;
        if (right < h->size && h->data[right].f_score < h->data[smallest].f_score)
            smallest = right;

        if (smallest == i) break;  /* déjà à la bonne place */

        HeapNode tmp      = h->data[i];
        h->data[i]        = h->data[smallest];
        h->data[smallest] = tmp;
        i = smallest;
    }
    return result;
}

static double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371000.0;  /* rayon moyen de la Terre en mètres */
    double phi1 = lat1 * M_PI / 180.0;
    double phi2 = lat2 * M_PI / 180.0;
    double dphi = (lat2 - lat1) * M_PI / 180.0;
    double dlam = (lon2 - lon1) * M_PI / 180.0;
    double a = sin(dphi/2)*sin(dphi/2)
             + cos(phi1)*cos(phi2)*sin(dlam/2)*sin(dlam/2);
    return 2.0 * R * asin(sqrt(a));
}

static ShortestPathResult dijkstra(const CSRGraph *g, int source, int target) {
    ShortestPathResult res = {-1.0, NULL, 0, 0, 0, 0.0};
    double t0 = now_ms();

    int n = g->n_nodes;

    /* dist[v] = meilleure distance connue de source à v */
    double *dist = malloc(n * sizeof(double));
    /* prev[v] = prédécesseur de v sur le meilleur chemin (-1 si inconnu) */
    int    *prev = malloc(n * sizeof(int));

    /* Initialisation : tout à +infini */
    for (int i = 0; i < n; i++) { dist[i] = DBL_MAX; prev[i] = -1; }
    dist[source] = 0.0;

    /* Créer le tas avec une capacité initiale raisonnable */
    MinHeap *heap = heap_create(1024);
    heap_push(heap, source, 0.0);

    while (heap->size > 0) {
        HeapNode cur = heap_pop(heap);
        int    u = cur.node;
        double f = cur.f_score;
        res.nodes_explored++;

        /* Entrée périmée : une meilleure distance a déjà été trouvée pour u */
        if (f > dist[u]) continue;

        /* Cible atteinte : on peut s'arrêter */
        if (u == target) break;

        /* Relaxation des voisins de u */
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            int    v        = g->adj[k];
            double new_dist = dist[u] + g->weights[k];
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                prev[v] = u;
                heap_push(heap, v, new_dist);
                res.relaxations++;
            }
        }
    }

    heap_free(heap);

    /* Reconstruction du chemin si target est accessible */
    if (dist[target] < DBL_MAX) {
        res.dist = dist[target];
        /* Compter la longueur du chemin en remontant */
        int len = 0;
        for (int v = target; v != -1; v = prev[v]) len++;
        res.path     = malloc(len * sizeof(int));
        res.path_len = len;
        /* Remplir le chemin à l'envers puis inverser */
        int idx = len - 1;
        for (int v = target; v != -1; v = prev[v]) res.path[idx--] = v;
    }

    free(dist);
    free(prev);

    double t1 = now_ms();
    res.time_ms = t1 - t0;
    return res;
}


/* A* */
static ShortestPathResult astar(const CSRGraph *g,
                                 const Node *nodes,
                                 int source, int target) {
    ShortestPathResult res = {-1.0, NULL, 0, 0, 0, 0.0};
    double t0 = now_ms();

    int n = g->n_nodes;

    /* g_score[v] = dist réelle connue de source à v (équivalent de dist[] dans Dijkstra) */
    double *g_score = malloc(n * sizeof(double));
    int    *prev    = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) { g_score[i] = DBL_MAX; prev[i] = -1; }
    g_score[source] = 0.0;

    /* Coordonnées de la cible (pré-calculées une fois pour toutes) */
    double tlat = nodes[target].lat;
    double tlon = nodes[target].lon;

    MinHeap *heap = heap_create(1024);
    /* Pour la source : f = g(source) + h(source, target) = 0 + haversine(...) */
    double h_source = haversine(nodes[source].lat, nodes[source].lon, tlat, tlon);
    heap_push(heap, source, h_source);

    while (heap->size > 0) {
        HeapNode cur = heap_pop(heap);
        int    u = cur.node;
        /* cur.f_score est f = g+h, mais on compare avec g_score[u] + h(u) */
        /* Entrée périmée : vérification sur g_score, pas f_score */
        double h_u = haversine(nodes[u].lat, nodes[u].lon, tlat, tlon);
        if (cur.f_score > g_score[u] + h_u) continue;

        res.nodes_explored++;

        if (u == target) break;  /* cible atteinte */

        /* Relaxation : identique à Dijkstra */
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            int    v        = g->adj[k];
            double new_g    = g_score[u] + g->weights[k];
            if (new_g < g_score[v]) {
                g_score[v] = new_g;
                prev[v]    = u;
                /* Clé du tas = f = g(v) + h(v, target)     ← seule différence avec Dijkstra */
                double h_v = haversine(nodes[v].lat, nodes[v].lon, tlat, tlon);
                heap_push(heap, v, new_g + h_v);
                res.relaxations++;
            }
        }
    }

    heap_free(heap);

    if (g_score[target] < DBL_MAX) {
        res.dist = g_score[target];
        int len = 0;
        for (int v = target; v != -1; v = prev[v]) len++;
        res.path     = malloc(len * sizeof(int));
        res.path_len = len;
        int idx = len - 1;
        for (int v = target; v != -1; v = prev[v]) res.path[idx--] = v;
    }

    free(g_score);
    free(prev);

    double t1 = now_ms();
    res.time_ms = t1 - t0;
    return res;
}

/*  DIJKSTRA "ALL-NODES" + GRAPHE INVERSE */
static void dijkstra_all(const CSRGraph *g, int source, double *dist) {
    int n = g->n_nodes;
    for (int i = 0; i < n; i++) dist[i] = DBL_MAX;
    dist[source] = 0.0;

    MinHeap *heap = heap_create(1024);
    heap_push(heap, source, 0.0);

    while (heap->size > 0) {
        HeapNode cur = heap_pop(heap);
        int    u = cur.node;
        double f = cur.f_score;
        if (f > dist[u]) continue;  /* entrée périmée */

        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            int    v        = g->adj[k];
            double new_dist = dist[u] + g->weights[k];
            if (new_dist < dist[v]) {
                dist[v] = new_dist;
                heap_push(heap, v, new_dist);
            }
        }
    }
    heap_free(heap);
}


/* Construction du graphe inverse 
* À partir d'un graphe CSR direct, on construit son graphe inverse.
*/

static CSRGraph *build_reverse_csr(const CSRGraph *g) {
    int n = g->n_nodes;
    int m = g->n_edges;

    RawEdge *rev = malloc(m * sizeof(RawEdge));
    int idx = 0;
    for (int u = 0; u < n; u++) {
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            rev[idx].u = g->adj[k];
            rev[idx].v = u;
            rev[idx].w = g->weights[k];
            idx++;
        }
    }
    CSRGraph *gr = build_csr(n, rev, m);
    free(rev);
    return gr;
}


/* ALT : A* AVEC LANDMARKS (PRÉTRAITEMENTS) */

typedef struct {
    int       num_landmarks;
    int      *landmarks;     /* tableau des indices des landmarks (taille K) */
    double  **dist_from;     /* dist_from[i][v] = d(landmark_i, v) — forward  */
    double  **dist_to;       /* dist_to[i][v]   = d(v, landmark_i) — backward */
    double    prep_time_ms;  /* temps de prétraitement total                  */
    size_t    memory_bytes;  /* mémoire utilisée (estimation)                 */
} ALTData;

static void select_landmarks_farthest(const CSRGraph *g, int K,
                                      int *landmarks, unsigned seed) {
    int n = g->n_nodes;

    int seed_node = 0, max_deg = -1;
    for (int v = 0; v < n; v++) {
        int d = g->row_ptr[v + 1] - g->row_ptr[v];
        if (d > max_deg) { max_deg = d; seed_node = v; }
    }

    char *in_main = calloc(n, sizeof(char));
    int  *queue   = malloc(n * sizeof(int));
    int qhead = 0, qtail = 0;
    in_main[seed_node] = 1;
    queue[qtail++] = seed_node;
    while (qhead < qtail) {
        int u = queue[qhead++];
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            int v = g->adj[k];
            if (!in_main[v]) { in_main[v] = 1; queue[qtail++] = v; }
        }
    }
    int main_size = qtail;
    printf("    [info] plus grande composante atteignable : %d / %d nœuds (%.1f%%)\n",
           main_size, n, 100.0 * main_size / n);
    free(queue);

    /* On collecte les nœuds de la grande composante dans un tableau,
     * pour pouvoir piocher facilement dedans. */
    int *main_nodes = malloc(main_size * sizeof(int));
    int idx = 0;
    for (int v = 0; v < n; v++) if (in_main[v]) main_nodes[idx++] = v;

    /* ── Étape 2 : sélection farthest dans la grande composante ── */
    double *min_dist = malloc(n * sizeof(double));
    double *tmp_dist = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) min_dist[i] = DBL_MAX;

    srand(seed);

    /* Premier landmark : un nœud aléatoire DANS la grande composante */
    int first = main_nodes[rand() % main_size];
    landmarks[0] = first;

    dijkstra_all(g, first, tmp_dist);
    for (int v = 0; v < n; v++)
        if (tmp_dist[v] < min_dist[v]) min_dist[v] = tmp_dist[v];

    /* Landmarks suivants : argmax de min_dist, restreint à la grande composante */
    for (int i = 1; i < K; i++) {
        int    best_v   = -1;
        double best_val = -1.0;
        for (int j = 0; j < main_size; j++) {
            int v = main_nodes[j];
            if (min_dist[v] == DBL_MAX) continue;
            if (min_dist[v] > best_val) {
                best_val = min_dist[v];
                best_v   = v;
            }
        }
        if (best_v < 0) {
            fprintf(stderr,
                "  [warn] sélection farthest : impossible de placer le "
                "landmark %d\n", i);
            best_v = main_nodes[rand() % main_size];
        }
        landmarks[i] = best_v;

        dijkstra_all(g, best_v, tmp_dist);
        for (int v = 0; v < n; v++)
            if (tmp_dist[v] < min_dist[v]) min_dist[v] = tmp_dist[v];
    }

    free(min_dist); free(tmp_dist);
    free(main_nodes); free(in_main);
}

static ALTData *alt_preprocess(const CSRGraph *g, const CSRGraph *g_rev,
                                int K, unsigned seed) {
    double t0 = now_ms();

    ALTData *alt = malloc(sizeof(ALTData));
    alt->num_landmarks = K;
    alt->landmarks     = malloc(K * sizeof(int));
    alt->dist_from     = malloc(K * sizeof(double *));
    alt->dist_to       = malloc(K * sizeof(double *));

    int n = g->n_nodes;

    /* 1) Sélection */
    printf("    Sélection des %d landmarks (farthest) ...\n", K);
    select_landmarks_farthest(g, K, alt->landmarks, seed);

    /* 2) Distances FORWARD : un Dijkstra depuis chaque landmark sur G */
    printf("    Calcul des distances forward (G) ...\n");
    for (int i = 0; i < K; i++) {
        alt->dist_from[i] = malloc(n * sizeof(double));
        dijkstra_all(g, alt->landmarks[i], alt->dist_from[i]);
    }

    /* 3) Distances BACKWARD : un Dijkstra depuis chaque landmark sur G' */
    /* dist_to[i][v] = d(v, ℓ_i) = distance v→ℓ_i sur G
     * Or distance v→ℓ_i sur G = distance ℓ_i→v sur G' (graphe inverse)
     * Donc : un Dijkstra depuis ℓ_i sur G' nous donne directement dist_to. */
    printf("    Calcul des distances backward (G') ...\n");
    for (int i = 0; i < K; i++) {
        alt->dist_to[i] = malloc(n * sizeof(double));
        dijkstra_all(g_rev, alt->landmarks[i], alt->dist_to[i]);
    }


    double t1 = now_ms();
    alt->prep_time_ms = t1 - t0;
    alt->memory_bytes = (size_t)2 * K * n * sizeof(double)
                      + K * sizeof(int);
    return alt;
}

static void alt_free(ALTData *alt) {
    for (int i = 0; i < alt->num_landmarks; i++) {
        free(alt->dist_from[i]);
        free(alt->dist_to[i]);
    }
    free(alt->dist_from);
    free(alt->dist_to);
    free(alt->landmarks);
    free(alt);
}

static inline double alt_h(const ALTData *alt, int v,
                            const double *C_fwd, const double *C_bwd) {
    double best = 0.0;
    int K = alt->num_landmarks;
    for (int i = 0; i < K; i++) {
        double dlv_fwd = alt->dist_from[i][v];
        double dlv_bwd = alt->dist_to[i][v];

        /* Forward : d(ℓ, t) − d(ℓ, v) */
        if (dlv_fwd != DBL_MAX && C_fwd[i] != DBL_MAX) {
            double f = C_fwd[i] - dlv_fwd;
            if (f > best) best = f;
        }
        /* Backward : d(v, ℓ) − d(t, ℓ) */
        if (dlv_bwd != DBL_MAX && C_bwd[i] != DBL_MAX) {
            double b = dlv_bwd - C_bwd[i];
            if (b > best) best = b;
        }
    }
    return best;
}


/*  ALT : LA REQUÊTE */
static ShortestPathResult alt_query(const CSRGraph *g,
                                     const ALTData *alt,
                                     int source, int target) {
    ShortestPathResult res = {-1.0, NULL, 0, 0, 0, 0.0};
    double t0 = now_ms();

    int n = g->n_nodes;
    int K = alt->num_landmarks;

    /* Constantes de l'heuristique pour cette cible */
    double *C_fwd = malloc(K * sizeof(double));
    double *C_bwd = malloc(K * sizeof(double));
    for (int i = 0; i < K; i++) {
        C_fwd[i] = alt->dist_from[i][target];
        C_bwd[i] = alt->dist_to[i][target];
    }

    double *g_score = malloc(n * sizeof(double));
    int    *prev    = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) { g_score[i] = DBL_MAX; prev[i] = -1; }
    g_score[source] = 0.0;

    MinHeap *heap = heap_create(1024);
    double h_source = alt_h(alt, source, C_fwd, C_bwd);
    heap_push(heap, source, h_source);

    while (heap->size > 0) {
        HeapNode cur = heap_pop(heap);
        int    u = cur.node;
        double h_u = alt_h(alt, u, C_fwd, C_bwd);
        if (cur.f_score > g_score[u] + h_u) continue;

        res.nodes_explored++;
        if (u == target) break;

        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            int    v     = g->adj[k];
            double new_g = g_score[u] + g->weights[k];
            if (new_g < g_score[v]) {
                g_score[v] = new_g;
                prev[v]    = u;
                double h_v = alt_h(alt, v, C_fwd, C_bwd);
                heap_push(heap, v, new_g + h_v);
                res.relaxations++;
            }
        }
    }
    heap_free(heap);

    if (g_score[target] < DBL_MAX) {
        res.dist = g_score[target];
        int len = 0;
        for (int v = target; v != -1; v = prev[v]) len++;
        res.path     = malloc(len * sizeof(int));
        res.path_len = len;
        int idx = len - 1;
        for (int v = target; v != -1; v = prev[v]) res.path[idx--] = v;
    }

    free(g_score); free(prev);
    free(C_fwd);   free(C_bwd);

    double t1 = now_ms();
    res.time_ms = t1 - t0;
    return res;
}


/* CONTRACTION HIERARCHIES (CH) */

typedef struct {
    int    to;
    double w;
    int    via;   /* -1 = arête originale ; sinon indice du nœud contracté milieu */
} CHEdge;

typedef struct {
    int        n_nodes;

    /* Adjacence dynamique : un tableau par nœud, qui grandit. */
    CHEdge   **out;
    int       *out_len;
    int       *out_cap;
    CHEdge   **in;
    int       *in_len;
    int       *in_cap;

    int       *level;       /* level[v] = ordre de contraction (0..N-1) ; -1 pas encore */
    int       *contracted;  /* 1 si v déjà contracté */
    int        next_level;  /* compteur d'ordre */
} CHGraph;


/* Init structures dynamiques avec une capacité initiale par nœud
 * égale au degré original (in et out), pour limiter les realloc. */
static CHGraph *ch_init(const CSRGraph *g, const CSRGraph *g_rev) {
    int n = g->n_nodes;
    CHGraph *ch = malloc(sizeof(CHGraph));
    ch->n_nodes    = n;
    ch->out        = malloc(n * sizeof(CHEdge *));
    ch->out_len    = calloc(n, sizeof(int));
    ch->out_cap    = malloc(n * sizeof(int));
    ch->in         = malloc(n * sizeof(CHEdge *));
    ch->in_len     = calloc(n, sizeof(int));
    ch->in_cap     = malloc(n * sizeof(int));
    ch->level      = malloc(n * sizeof(int));
    ch->contracted = calloc(n, sizeof(int));
    ch->next_level = 0;
    for (int v = 0; v < n; v++) ch->level[v] = -1;

    /* Copie des arêtes originales depuis les CSR */
    for (int u = 0; u < n; u++) {
        int d_out = g->row_ptr[u + 1] - g->row_ptr[u];
        int cap   = d_out > 4 ? d_out : 4;
        ch->out_cap[u] = cap;
        ch->out[u]     = malloc(cap * sizeof(CHEdge));
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            CHEdge e = { .to = g->adj[k], .w = g->weights[k], .via = -1 };
            ch->out[u][ch->out_len[u]++] = e;
        }

        int d_in = g_rev->row_ptr[u + 1] - g_rev->row_ptr[u];
        cap      = d_in > 4 ? d_in : 4;
        ch->in_cap[u] = cap;
        ch->in[u]     = malloc(cap * sizeof(CHEdge));
        /* Dans g_rev, les voisins sortants de u sont les "prédécesseurs" de u dans G */
        for (int k = g_rev->row_ptr[u]; k < g_rev->row_ptr[u + 1]; k++) {
            CHEdge e = { .to = g_rev->adj[k], .w = g_rev->weights[k], .via = -1 };
            ch->in[u][ch->in_len[u]++] = e;
        }
    }
    return ch;
}

static void ch_free_dynamic(CHGraph *ch) {
    int n = ch->n_nodes;
    for (int v = 0; v < n; v++) { free(ch->out[v]); free(ch->in[v]); }
    free(ch->out);     free(ch->out_len); free(ch->out_cap);
    free(ch->in);      free(ch->in_len);  free(ch->in_cap);
    /* level et contracted restent valides après ; libérés à part */
}


/* Ajouter une arête (u → v, weight w, via mid) dans les deux listes. */
static inline void ch_push_out(CHGraph *ch, int u, int to, double w, int via) {
    if (ch->out_len[u] == ch->out_cap[u]) {
        ch->out_cap[u] = ch->out_cap[u] * 2 + 1;
        ch->out[u] = realloc(ch->out[u], ch->out_cap[u] * sizeof(CHEdge));
    }
    ch->out[u][ch->out_len[u]++] = (CHEdge){ .to = to, .w = w, .via = via };
}
static inline void ch_push_in(CHGraph *ch, int v, int from, double w, int via) {
    if (ch->in_len[v] == ch->in_cap[v]) {
        ch->in_cap[v] = ch->in_cap[v] * 2 + 1;
        ch->in[v] = realloc(ch->in[v], ch->in_cap[v] * sizeof(CHEdge));
    }
    ch->in[v][ch->in_len[v]++] = (CHEdge){ .to = from, .w = w, .via = via };
}

/*
 * Si une arête u→to existe déjà avec un poids ≥ new_w, on la met à jour.
 * Sinon on en ajoute une nouvelle. (Évite l'accumulation de raccourcis
 * redondants entre les mêmes nœuds.) Renvoie 1 si on a fait une modif,
 * 0 si l'arête existante était meilleure (donc rien fait).
 */
static int ch_add_or_update(CHGraph *ch, int u, int to, double w, int via) {
    /* Recherche d'une arête existante u→to */
    for (int k = 0; k < ch->out_len[u]; k++) {
        if (ch->out[u][k].to == to) {
            if (w < ch->out[u][k].w) {
                ch->out[u][k].w   = w;
                ch->out[u][k].via = via;
                /* Trouver l'arête correspondante dans in[to] et la mettre à jour */
                for (int j = 0; j < ch->in_len[to]; j++) {
                    if (ch->in[to][j].to == u) {
                        ch->in[to][j].w   = w;
                        ch->in[to][j].via = via;
                        break;
                    }
                }
                return 1;
            }
            return 0;  /* arête existante déjà aussi bonne ou meilleure */
        }
    }
    /* Nouvelle arête */
    ch_push_out(ch, u, to, w, via);
    ch_push_in(ch, to, u, w, via);
    return 1;
}

typedef struct {
    int       n;
    double   *value;
    int      *epoch;
    int       cur_epoch;
    MinHeap  *heap;
} WitnessBuf;

static WitnessBuf *witness_create(int n) {
    WitnessBuf *wb = malloc(sizeof(WitnessBuf));
    wb->n         = n;
    wb->value     = malloc(n * sizeof(double));
    wb->epoch     = calloc(n, sizeof(int));
    wb->cur_epoch = 0;
    wb->heap      = heap_create(64);
    return wb;
}
static void witness_free(WitnessBuf *wb) {
    free(wb->value); free(wb->epoch);
    heap_free(wb->heap);
    free(wb);
}
static inline void witness_reset(WitnessBuf *wb) {
    wb->cur_epoch++;
    wb->heap->size = 0;
}
static inline double witness_get(const WitnessBuf *wb, int v) {
    return (wb->epoch[v] == wb->cur_epoch) ? wb->value[v] : DBL_MAX;
}
static inline void witness_set(WitnessBuf *wb, int v, double d) {
    wb->epoch[v] = wb->cur_epoch;
    wb->value[v] = d;
}


/* CH : RECHERCHE DE TÉMOINS, CONTRACTION, ORDRE */

/* Limites de la recherche de témoins (à régler selon le graphe) */
#define WITNESS_HOP_LIMIT     5      /* nombre max d'arêtes du témoin */
#define WITNESS_NODE_LIMIT  500      /* nombre max de pop avant abandon */

static void witness_search(const CHGraph *ch, WitnessBuf *wb,
                            int u, int forbidden,
                            double max_dist, int hop_limit) {
    witness_reset(wb);
    witness_set(wb, u, 0.0);
    /* On stocke (node, f) dans le tas. Pour le hop, on utilise un tableau
     * séparé : on n'a pas la place dans HeapNode. Approche simple : on
     * recalcule à la volée — pas vraiment besoin, on tracke avec un tableau
     * éphémère. Trop lourd. Variante : on borne juste sur node_limit. */
    heap_push(wb->heap, u, 0.0);

    int popped = 0;
    while (wb->heap->size > 0 && popped < WITNESS_NODE_LIMIT) {
        HeapNode cur = heap_pop(wb->heap);
        int x = cur.node;
        if (cur.f_score > witness_get(wb, x)) continue;
        if (cur.f_score > max_dist) break;  /* tas trié → tout le reste est pire */

        popped++;
        (void)hop_limit;  /* hop_limit non géré strictement ici, on borne par node_limit */

        for (int k = 0; k < ch->out_len[x]; k++) {
            int    y = ch->out[x][k].to;
            if (ch->contracted[y]) continue;     /* skip nœuds déjà contractés */
            if (y == forbidden)    continue;     /* interdit : c'est v */
            double new_d = cur.f_score + ch->out[x][k].w;
            if (new_d > max_dist)  continue;     /* dépasse → pas de témoin via ce chemin */
            if (new_d < witness_get(wb, y)) {
                witness_set(wb, y, new_d);
                heap_push(wb->heap, y, new_d);
            }
        }
    }
}


/*
 * Calcule (ou simule) les raccourcis qui devraient être ajoutés en
 * contractant v, et les renvoie dans le tableau (alloué par l'appelant)
 * shortcuts[] de capacité cap. Le nombre de raccourcis produit est
 * renvoyé via *out_count. Si le tableau est NULL, on se contente de
 * compter (mode "simulation" pour calculer ED).
 */
typedef struct { int u, w; double weight; } CHShortcut;

static int ch_simulate_or_apply(CHGraph *ch, WitnessBuf *wb, int v,
                                 CHShortcut *shortcuts, int cap,
                                 int apply) {
    int n_shortcuts = 0;

    for (int i = 0; i < ch->in_len[v]; i++) {
        int    u  = ch->in[v][i].to;
        double wu = ch->in[v][i].w;
        if (ch->contracted[u]) continue;

        /* max_dCH pour cette source */
        double max_dCH = 0.0;
        int has_target = 0;
        for (int j = 0; j < ch->out_len[v]; j++) {
            int    w  = ch->out[v][j].to;
            if (ch->contracted[w]) continue;
            if (w == u) continue;
            double ww = ch->out[v][j].w;
            double dCH = wu + ww;
            if (!has_target || dCH > max_dCH) max_dCH = dCH;
            has_target = 1;
        }
        if (!has_target) continue;

        /* Recherche de témoins depuis u */
        witness_search(ch, wb, u, v, max_dCH, WITNESS_HOP_LIMIT);

        /* Décision pour chaque cible w */
        for (int j = 0; j < ch->out_len[v]; j++) {
            int    w  = ch->out[v][j].to;
            if (ch->contracted[w]) continue;
            if (w == u) continue;
            double ww = ch->out[v][j].w;
            double dCH = wu + ww;
            double dwit = witness_get(wb, w);

            /* Tolérance numérique : si égal à dCH, témoin existe */
            if (dwit <= dCH + 1e-9) continue;

            /* Pas de témoin → raccourci nécessaire */
            if (apply) {
                ch_add_or_update(ch, u, w, dCH, v);
            }
            if (shortcuts && n_shortcuts < cap) {
                shortcuts[n_shortcuts].u = u;
                shortcuts[n_shortcuts].w = w;
                shortcuts[n_shortcuts].weight = dCH;
            }
            n_shortcuts++;
        }
    }
    return n_shortcuts;
}


/* Compter les voisins NON contractés de v, dans les deux sens. */
static int ch_count_active_neighbors(const CHGraph *ch, int v) {
    int c = 0;
    for (int k = 0; k < ch->in_len[v]; k++)
        if (!ch->contracted[ch->in[v][k].to]) c++;
    for (int k = 0; k < ch->out_len[v]; k++)
        if (!ch->contracted[ch->out[v][k].to]) c++;
    return c;
}


/* Edge difference : nbre de raccourcis simulés − degré actif. */
static int ch_edge_difference(CHGraph *ch, WitnessBuf *wb, int v) {
    int sc = ch_simulate_or_apply(ch, wb, v, NULL, 0, /*apply=*/0);
    int dg = ch_count_active_neighbors(ch, v);
    return sc - dg;
}


/* Contracter effectivement v : ajouter les raccourcis et marquer contracté. */
static void ch_contract_node(CHGraph *ch, WitnessBuf *wb, int v) {
    /* On applique les raccourcis. La capacité n'est pas vraiment limitée
     * ici puisqu'on passe NULL — on ne stocke pas dans un buffer externe. */
    ch_simulate_or_apply(ch, wb, v, NULL, 0, /*apply=*/1);
    ch->contracted[v] = 1;
    ch->level[v]      = ch->next_level++;
}


/* CH : BOUCLE DE CONTRACTION + PLIAGE EN CSR */

typedef struct {
    int     n_nodes;
    int     n_edges;
    int    *row_ptr;
    int    *adj;
    double *weights;
    int    *via;       /* taille n_edges */
} CSRGraphCH;

static void csr_ch_free(CSRGraphCH *g) {
    free(g->row_ptr); free(g->adj); free(g->weights); free(g->via); free(g);
}

/*
 * CHIndex : tout ce dont la requête CH a besoin.
 * Construit en une fois après prétraitement, peut servir des milliers de
 * requêtes ensuite.
 */
typedef struct {
    int          n_nodes;
    int         *level;        /* level[v] = ordre de contraction */
    CSRGraphCH  *up;           /* graphe augmenté direct (toutes les arêtes) */
    CSRGraphCH  *up_rev;       /* graphe augmenté inversé */
    int          num_shortcuts;
    double       prep_time_ms;
    size_t       memory_bytes;
} CHIndex;


/*
 * Boucle de contraction principale.
 * Stratégie : tas min sur edge difference, mises à jour paresseuses.
 *   - Initial : ED(v) pour tous les v actifs → tas.
 *   - Boucle : pop v, recalcule ED(v), compare au top du tas.
 *     Si v n'est plus le minimum, on le repousse. Sinon on contracte.
 *   - Après contraction : repousser les voisins actifs avec leur nouvelle ED.
 */
static void ch_contract_all(CHGraph *ch, WitnessBuf *wb, int verbose) {
    int n = ch->n_nodes;
    MinHeap *pq = heap_create(n);

    /* Initialisation : on calcule ED pour tous les nœuds */
    if (verbose) printf("    Initialisation des priorités (%d nœuds) ...\n", n);
    for (int v = 0; v < n; v++) {
        int ed = ch_edge_difference(ch, wb, v);
        heap_push(pq, v, (double)ed);
    }

    int contracted_count = 0;
    int progress_step    = n / 20 > 1 ? n / 20 : 1;

    while (pq->size > 0) {
        HeapNode top = heap_pop(pq);
        int v = top.node;
        if (ch->contracted[v]) continue;

        /* Recalcul de la priorité */
        int ed_now = ch_edge_difference(ch, wb, v);

        /* Comparer au minimum courant du tas (sans le pop) */
        if (pq->size > 0 && (double)ed_now > pq->data[0].f_score) {
            heap_push(pq, v, (double)ed_now);
            continue;
        }

        /* Contracter v */
        ch_contract_node(ch, wb, v);
        contracted_count++;

        if (verbose && contracted_count % progress_step == 0) {
            printf("    [%d/%d] contractés (%.0f%%)\n",
                   contracted_count, n,
                   100.0 * contracted_count / n);
        }

        /* Pousser les voisins actifs avec une nouvelle ED. Les anciennes
         * entrées resteront dans le tas mais seront détectées au pop. */
        /* Voisins out actifs */
        for (int k = 0; k < ch->out_len[v]; k++) {
            int y = ch->out[v][k].to;
            if (!ch->contracted[y]) {
                int ed_y = ch_edge_difference(ch, wb, y);
                heap_push(pq, y, (double)ed_y);
            }
        }
        /* Voisins in actifs */
        for (int k = 0; k < ch->in_len[v]; k++) {
            int y = ch->in[v][k].to;
            if (!ch->contracted[y]) {
                int ed_y = ch_edge_difference(ch, wb, y);
                heap_push(pq, y, (double)ed_y);
            }
        }
    }
    heap_free(pq);

    if (verbose) printf("    Contraction terminée : %d nœuds.\n", contracted_count);
}


/*
 * Convertit le graphe dynamique CHGraph (post-contraction) en CSR.
 * On copie toutes les arêtes (originales + raccourcis) dans un CSR
 * direct. Si reverse=1, on inverse le sens (pour le graphe inverse).
 *
 * Les arêtes incidentes à des nœuds isolés (sans voisin) sont gardées
 * comme un "trou" dans row_ptr (différence row_ptr[v+1]-row_ptr[v] = 0).
 */
static CSRGraphCH *ch_to_csr(const CHGraph *ch, int reverse) {
    int n = ch->n_nodes;

    /* Compter le total d'arêtes (= somme des out_len, c'est exact car on
     * a maintenu out et in en symétrie) */
    int m = 0;
    for (int u = 0; u < n; u++) m += ch->out_len[u];

    CSRGraphCH *g = malloc(sizeof(CSRGraphCH));
    g->n_nodes = n;
    g->n_edges = m;
    g->row_ptr = calloc(n + 1, sizeof(int));
    g->adj     = malloc(m * sizeof(int));
    g->weights = malloc(m * sizeof(double));
    g->via     = malloc(m * sizeof(int));

    if (!reverse) {
        /* CSR direct : facile, on suit ch->out tel quel */
        int idx = 0;
        for (int u = 0; u < n; u++) {
            g->row_ptr[u] = idx;
            for (int k = 0; k < ch->out_len[u]; k++) {
                g->adj[idx]     = ch->out[u][k].to;
                g->weights[idx] = ch->out[u][k].w;
                g->via[idx]     = ch->out[u][k].via;
                idx++;
            }
        }
        g->row_ptr[n] = idx;
    } else {
        /* CSR inverse : pour chaque (u → v) on insère (v → u).
         * On utilise un build en 3 passes comme build_csr. */
        for (int u = 0; u < n; u++)
            for (int k = 0; k < ch->out_len[u]; k++)
                g->row_ptr[ ch->out[u][k].to + 1 ]++;
        for (int i = 1; i <= n; i++) g->row_ptr[i] += g->row_ptr[i - 1];

        int *cursor = malloc(n * sizeof(int));
        memcpy(cursor, g->row_ptr, n * sizeof(int));
        for (int u = 0; u < n; u++) {
            for (int k = 0; k < ch->out_len[u]; k++) {
                int v   = ch->out[u][k].to;
                int pos = cursor[v]++;
                g->adj[pos]     = u;                  /* l'inverse part de v vers u */
                g->weights[pos] = ch->out[u][k].w;
                g->via[pos]     = ch->out[u][k].via;
            }
        }
        free(cursor);
    }
    return g;
}


/*
 * Prétraitement CH complet :
 *   - copie le graphe original dans CHGraph (dynamique)
 *   - boucle de contraction
 *   - convertit en CSR (direct + inverse)
 *   - retourne un CHIndex prêt à servir les requêtes
 */
static CHIndex *ch_preprocess(const CSRGraph *g, const CSRGraph *g_rev,
                                int verbose) {
    double t0 = now_ms();

    int n = g->n_nodes;
    if (verbose) printf("    Initialisation graphe dynamique CH ...\n");
    CHGraph *ch = ch_init(g, g_rev);

    int initial_edges = 0;
    for (int u = 0; u < n; u++) initial_edges += ch->out_len[u];

    WitnessBuf *wb = witness_create(n);

    if (verbose) printf("    Contraction (peut être long) ...\n");
    ch_contract_all(ch, wb, verbose);

    /* Compter les raccourcis ajoutés */
    int total_edges = 0;
    for (int u = 0; u < n; u++) total_edges += ch->out_len[u];
    int n_shortcuts = total_edges - initial_edges;

    /* Pliage en CSR */
    if (verbose) printf("    Pliage en CSR (direct + inverse) ...\n");
    CHIndex *idx = malloc(sizeof(CHIndex));
    idx->n_nodes        = n;
    idx->level          = malloc(n * sizeof(int));
    memcpy(idx->level, ch->level, n * sizeof(int));
    idx->up             = ch_to_csr(ch, 0);
    idx->up_rev         = ch_to_csr(ch, 1);
    idx->num_shortcuts  = n_shortcuts;

    witness_free(wb);
    ch_free_dynamic(ch);
    /* level et contracted dans ch sont libérés ici, mais level a été copié */
    free(ch->level); free(ch->contracted); free(ch);

    double t1 = now_ms();
    idx->prep_time_ms = t1 - t0;

    /* Estimation mémoire : 2 CSR (chacun adj+weights+via+row_ptr) + level */
    size_t bytes = 0;
    bytes += (size_t)idx->up->n_edges       * (sizeof(int) + sizeof(double) + sizeof(int));
    bytes += (size_t)idx->up_rev->n_edges   * (sizeof(int) + sizeof(double) + sizeof(int));
    bytes += (size_t)(n + 1) * 2 * sizeof(int);  /* row_ptr × 2 */
    bytes += (size_t)n * sizeof(int);            /* level */
    idx->memory_bytes = bytes;

    return idx;
}

static void ch_index_free(CHIndex *idx) {
    free(idx->level);
    csr_ch_free(idx->up);
    csr_ch_free(idx->up_rev);
    free(idx);
}


/* CH : REQUÊTE BIDIRECTIONNELLE ET UNPACKING */


/* Cherche dans le CSR direct l'arête (u → v) de poids minimum.
 * Retourne via via_out (peut valoir -1 si arête originale).
 * Renvoie la position k de l'arête, -1 si pas trouvée. */
static int find_edge(const CSRGraphCH *g, int u, int v, int *via_out) {
    int    best_k = -1;
    double best_w = DBL_MAX;
    for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
        if (g->adj[k] == v && g->weights[k] < best_w) {
            best_w = g->weights[k];
            best_k = k;
        }
    }
    if (best_k >= 0) *via_out = g->via[best_k];
    return best_k;
}


/* Tampon dynamique pour stocker un chemin en cours de reconstruction. */
typedef struct { int *data; int len; int cap; } IntVec;
static void iv_init(IntVec *v) { v->data = NULL; v->len = 0; v->cap = 0; }
static void iv_push(IntVec *v, int x) {
    if (v->len == v->cap) {
        v->cap = v->cap * 2 + 16;
        v->data = realloc(v->data, v->cap * sizeof(int));
    }
    v->data[v->len++] = x;
}


/* Unpack récursif d'une arête (u → v) du graphe augmenté.
 * On écrit la séquence de nœuds dans path en n'incluant PAS u au début
 * (on suppose que u a déjà été émis par l'appelant) ; on inclut v à la fin. */
static void unpack_edge(const CSRGraphCH *g, int u, int v, IntVec *path) {
    int via;
    int k = find_edge(g, u, v, &via);
    if (k < 0) {
        /* Sécurité : arête introuvable (ne devrait pas arriver) */
        iv_push(path, v);
        return;
    }
    if (via < 0) {
        /* Arête originale : on émet directement v */
        iv_push(path, v);
    } else {
        /* Raccourci : (u → via) puis (via → v) */
        unpack_edge(g, u,   via, path);
        unpack_edge(g, via, v,   path);
    }
}


/*
 * Requête CH bidirectionnelle.
 * Retourne distance + chemin reconstruit (avec unpacking complet).
 *
 * On utilise deux MinHeap, deux tableaux dist[], deux tableaux prev[]
 * et deux tableaux prev_via[] (pour pouvoir retrouver via lors du
 * unpacking, sans scan inutile).
 */
static ShortestPathResult ch_query(const CHIndex *idx, int source, int target) {
    ShortestPathResult res = {-1.0, NULL, 0, 0, 0, 0.0};
    double t0 = now_ms();

    int n = idx->n_nodes;
    const int *level = idx->level;

    /* Cas trivial */
    if (source == target) {
        res.dist     = 0.0;
        res.path     = malloc(sizeof(int));
        res.path[0]  = source;
        res.path_len = 1;
        double t1 = now_ms();
        res.time_ms = t1 - t0;
        return res;
    }

    double *df = malloc(n * sizeof(double));
    double *db = malloc(n * sizeof(double));
    int    *pf = malloc(n * sizeof(int));
    int    *pb = malloc(n * sizeof(int));
    int    *vf = malloc(n * sizeof(int));
    int    *vb = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        df[i] = DBL_MAX; db[i] = DBL_MAX;
        pf[i] = -1;      pb[i] = -1;
        vf[i] = -1;      vb[i] = -1;
    }
    df[source] = 0.0;
    db[target] = 0.0;

    MinHeap *hf = heap_create(1024);
    MinHeap *hb = heap_create(1024);
    heap_push(hf, source, 0.0);
    heap_push(hb, target, 0.0);

    double best_meet = DBL_MAX;
    int    meet_node = -1;

    while (hf->size > 0 || hb->size > 0) {
        /* Critère d'arrêt : la plus petite clé restante ≥ best_meet
         * ⇒ aucune amélioration future possible. */
        double top_f = (hf->size > 0) ? hf->data[0].f_score : DBL_MAX;
        double top_b = (hb->size > 0) ? hb->data[0].f_score : DBL_MAX;
        double min_top = (top_f < top_b) ? top_f : top_b;
        if (min_top >= best_meet) break;

        /* On avance celui dont le top est le plus petit */
        if (top_f <= top_b) {
            HeapNode cur = heap_pop(hf);
            int u = cur.node;
            if (cur.f_score > df[u]) continue;
            res.nodes_explored++;

            /* Mise à jour du best_meet via dist_bwd[u] */
            if (db[u] < DBL_MAX) {
                double total = df[u] + db[u];
                if (total < best_meet) { best_meet = total; meet_node = u; }
            }

            /* Relaxation : suivre les arêtes "montantes" dans idx->up */
            const CSRGraphCH *G = idx->up;
            for (int k = G->row_ptr[u]; k < G->row_ptr[u + 1]; k++) {
                int v = G->adj[k];
                if (level[v] <= level[u]) continue;       /* pas montant : ignore */
                double new_d = df[u] + G->weights[k];
                if (new_d < df[v]) {
                    df[v] = new_d;
                    pf[v] = u;
                    vf[v] = G->via[k];
                    heap_push(hf, v, new_d);
                    res.relaxations++;
                }
            }
        } else {
            HeapNode cur = heap_pop(hb);
            int u = cur.node;
            if (cur.f_score > db[u]) continue;
            res.nodes_explored++;

            if (df[u] < DBL_MAX) {
                double total = df[u] + db[u];
                if (total < best_meet) { best_meet = total; meet_node = u; }
            }

            const CSRGraphCH *G = idx->up_rev;
            for (int k = G->row_ptr[u]; k < G->row_ptr[u + 1]; k++) {
                int v = G->adj[k];
                if (level[v] <= level[u]) continue;
                double new_d = db[u] + G->weights[k];
                if (new_d < db[v]) {
                    db[v] = new_d;
                    pb[v] = u;
                    vb[v] = G->via[k];
                    heap_push(hb, v, new_d);
                    res.relaxations++;
                }
            }
        }
    }
    heap_free(hf); heap_free(hb);

    if (meet_node < 0 || best_meet == DBL_MAX) {
        free(df); free(db); free(pf); free(pb); free(vf); free(vb);
        double t1 = now_ms();
        res.time_ms = t1 - t0;
        return res;
    }

    res.dist = best_meet;

    /* Reconstruction du chemin avec unpacking */

    /* Étape 1 : récupérer la séquence d'arêtes (a→b avec via) du chemin. */
    /*   forward part : meet → source en remontant pf, à inverser
     *   backward part : meet → target en suivant pb */

    /* On collecte d'abord (en sens inverse) les paires (predecesseur, via)
     * du segment forward, puis on l'inverse. */
    typedef struct { int from, to, via; } EdgeRec;
    int  cap_e = 64;
    EdgeRec *edges = malloc(cap_e * sizeof(EdgeRec));
    int n_e = 0;

    /* Forward : remonter pf du meet à la source */
    {
        /* On construit en sens inverse, puis on inverse */
        int  cap_tmp = 64; int n_tmp = 0;
        EdgeRec *tmp = malloc(cap_tmp * sizeof(EdgeRec));
        int v = meet_node;
        while (pf[v] != -1) {
            int u = pf[v];
            if (n_tmp == cap_tmp) { cap_tmp *= 2; tmp = realloc(tmp, cap_tmp * sizeof(EdgeRec)); }
            tmp[n_tmp].from = u; tmp[n_tmp].to = v; tmp[n_tmp].via = vf[v];
            n_tmp++;
            v = u;
        }
        /* Insérer dans edges en ordre inversé (du source vers meet) */
        for (int i = n_tmp - 1; i >= 0; i--) {
            if (n_e == cap_e) { cap_e *= 2; edges = realloc(edges, cap_e * sizeof(EdgeRec)); }
            edges[n_e++] = tmp[i];
        }
        free(tmp);
    }

    /* Backward : suivre pb du meet vers la target */
    {
        int v = meet_node;
        while (pb[v] != -1) {
            int u = pb[v];
            if (n_e == cap_e) { cap_e *= 2; edges = realloc(edges, cap_e * sizeof(EdgeRec)); }
            edges[n_e].from = v; edges[n_e].to = u; edges[n_e].via = vb[v];
            n_e++;
            v = u;
        }
    }

    /* Étape 2 : unpacker chaque arête. On émet d'abord source, puis pour
     * chaque arête on délègue à unpack_edge qui ajoute les nœuds intermédiaires
     * et l'extrémité. */
    IntVec path; iv_init(&path);
    iv_push(&path, source);
    for (int i = 0; i < n_e; i++) {
        if (edges[i].via < 0) {
            iv_push(&path, edges[i].to);
        } else {
            unpack_edge(idx->up, edges[i].from, edges[i].to, &path);
        }
    }
    free(edges);

    res.path     = path.data;
    res.path_len = path.len;

    free(df); free(db); free(pf); free(pb); free(vf); free(vb);

    double t1 = now_ms();
    res.time_ms = t1 - t0;
    return res;
}

typedef struct {
    const char *nodes_path;
    const char *edges_path;
    const char *csv_out;
    int         num_landmarks;
    int         num_queries;
    unsigned    seed;
    int         skip_ch;     /* permet de désactiver CH si trop lent */
} Config;

static void print_usage(const char *prog) {
    fprintf(stderr,
        "Usage : %s [options] [nodes.csv] [edges.csv]\n"
        "Options :\n"
        "  --landmarks N    nombre de landmarks ALT (défaut 16)\n"
        "  --queries   N    nombre de requêtes pour le benchmark (défaut 100)\n"
        "  --csv FILE       écrire les résultats détaillés dans FILE\n"
        "  --seed S         graine du générateur (défaut 42)\n"
        "  --no-ch          désactiver CH (utile en dev sur gros graphe)\n"
        "  --help           afficher cette aide\n",
        prog);
}

static Config parse_args(int argc, char *argv[]) {
    Config c = {
        .nodes_path     = "nodes.csv",
        .edges_path     = "edges.csv",
        .csv_out        = NULL,
        .num_landmarks  = 16,
        .num_queries    = 100,
        .seed           = 42,
        .skip_ch        = 0,
    };
    int positional = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--help") || !strcmp(a, "-h")) {
            print_usage(argv[0]); exit(0);
        } else if (!strcmp(a, "--landmarks") && i + 1 < argc) {
            c.num_landmarks = atoi(argv[++i]);
        } else if (!strcmp(a, "--queries") && i + 1 < argc) {
            c.num_queries = atoi(argv[++i]);
        } else if (!strcmp(a, "--csv") && i + 1 < argc) {
            c.csv_out = argv[++i];
        } else if (!strcmp(a, "--seed") && i + 1 < argc) {
            c.seed = (unsigned)atoi(argv[++i]);
        } else if (!strcmp(a, "--no-ch")) {
            c.skip_ch = 1;
        } else if (a[0] == '-') {
            fprintf(stderr, "Option inconnue : %s\n", a);
            print_usage(argv[0]); exit(1);
        } else {
            if (positional == 0)      c.nodes_path = a;
            else if (positional == 1) c.edges_path = a;
            else { print_usage(argv[0]); exit(1); }
            positional++;
        }
    }
    return c;
}


/* Validation d'un chemin reconstruit */

static int validate_path(const CSRGraph *g,
                          const int *path, int len,
                          double expected_dist,
                          const char *who) {
    if (len <= 0) {
        fprintf(stderr, "  [ERR] %s : chemin vide\n", who); return 0;
    }
    if (len == 1) {
        if (fabs(expected_dist) > 1e-3) {
            fprintf(stderr, "  [ERR] %s : chemin de 1 nœud mais dist=%.3f\n",
                    who, expected_dist);
            return 0;
        }
        return 1;
    }
    double sum = 0.0;
    for (int i = 0; i + 1 < len; i++) {
        int u = path[i], v = path[i + 1];
        if (u < 0 || u >= g->n_nodes || v < 0 || v >= g->n_nodes) {
            fprintf(stderr, "  [ERR] %s : nœud hors-borne à i=%d\n", who, i);
            return 0;
        }
        /* Cherche la meilleure arête u→v dans G (on prend la plus courte en cas de doublons) */
        double best = DBL_MAX;
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++)
            if (g->adj[k] == v && g->weights[k] < best) best = g->weights[k];
        if (best == DBL_MAX) {
            fprintf(stderr, "  [ERR] %s : arête %d→%d inexistante dans G\n",
                    who, u, v);
            return 0;
        }
        sum += best;
    }
    if (fabs(sum - expected_dist) > 1e-3) {
        fprintf(stderr,
            "  [ERR] %s : somme des poids=%.3f mais dist annoncée=%.3f\n",
            who, sum, expected_dist);
        return 0;
    }
    return 1;
}


/*Statistiques (moyenne, p50, p95, débit)*/

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da < db) ? -1 : (da > db);
}

typedef struct {
    double mean_ms;
    double p50_ms;
    double p95_ms;
    double throughput;   /* requêtes/s */
    double total_ms;
    int    count;
} TimeStats;

/*
 * Calcule les stats à partir d'un tableau de temps (en ms).
 */
static TimeStats compute_stats(double *times_ms, int n) {
    TimeStats s = {0};
    if (n <= 0) return s;
    s.count = n;
    double sum = 0;
    for (int i = 0; i < n; i++) sum += times_ms[i];
    s.mean_ms  = sum / n;
    s.total_ms = sum;
    qsort(times_ms, n, sizeof(double), cmp_double);
    int i50 = (int)ceil(0.50 * n) - 1; if (i50 < 0) i50 = 0;
    int i95 = (int)ceil(0.95 * n) - 1; if (i95 < 0) i95 = 0;
    s.p50_ms = times_ms[i50];
    s.p95_ms = times_ms[i95];
    s.throughput = (sum > 0) ? (1000.0 * n / sum) : 0.0;
    return s;
}


/* Classification des requêtes par distance */
/*
 * On classe les requêtes en 3 catégories selon la distance Dijkstra :
 *   SHORT  : distance < 20 km   (intra-ville/agglo)
 *   MEDIUM : 20 - 80 km          (inter-villes courtes)
 *   LONG   : > 80 km             (régional)
 *
 * Ces seuils sont calibrés pour un graphe régional comme la Champagne-Ardenne.
 */
typedef enum { CAT_SHORT = 0, CAT_MEDIUM = 1, CAT_LONG = 2, CAT_N = 3 } Category;
static const char *CAT_NAMES[CAT_N] = {"SHORT", "MEDIUM", "LONG"};

static Category classify_distance(double dist_m) {
    if (dist_m < 20000.0) return CAT_SHORT;
    if (dist_m < 80000.0) return CAT_MEDIUM;
    return CAT_LONG;
}


/* Main */

int main(int argc, char *argv[]) {
    Config cfg = parse_args(argc, argv);

    printf("=== Route Planning : Dijkstra | A* | ALT | CH ===\n");
    printf("    nodes=%s edges=%s\n", cfg.nodes_path, cfg.edges_path);
    printf("    landmarks=%d queries=%d seed=%u%s\n",
           cfg.num_landmarks, cfg.num_queries, cfg.seed,
           cfg.skip_ch ? " (CH désactivé)" : "");
    if (cfg.csv_out) printf("    csv=%s\n", cfg.csv_out);
    printf("\n");

    /* Chargement et CSR */
    int      n_nodes;
    HashMap *map;
    printf("[1/4] Lecture de %s ...\n", cfg.nodes_path);
    Node *nodes = read_nodes(cfg.nodes_path, &n_nodes, &map);
    printf("      %d noeuds chargés\n", n_nodes);

    int      n_edges;
    printf("[2/4] Lecture de %s ...\n", cfg.edges_path);
    RawEdge *raw_edges = read_edges(cfg.edges_path, map, &n_edges);
    printf("      %d arêtes chargées\n", n_edges);

    printf("[3/4] Construction CSR direct + inverse ...\n");
    CSRGraph *g     = build_csr(n_nodes, raw_edges, n_edges);
    CSRGraph *g_rev = build_reverse_csr(g);
    printf("      G  : %d noeuds, %d arêtes\n", g->n_nodes, g->n_edges);
    printf("      G' : %d noeuds, %d arêtes\n", g_rev->n_nodes, g_rev->n_edges);

    free(raw_edges);
    map_free(map);

    /* Mémoire CSR de base : row_ptr (×2 G+G') + adj + weights (×2) */
    size_t mem_csr = (size_t)2 * (n_nodes + 1) * sizeof(int)
                   + (size_t)2 * n_edges * (sizeof(int) + sizeof(double))
                   + (size_t)n_nodes * sizeof(Node);
    printf("      Mémoire graphes : %.1f Mo\n", mem_csr / (1024.0 * 1024.0));


    /* Prétraitement ALT et CH */
    printf("\n[4/4] Prétraitements\n");
    printf("  > ALT (%d landmarks)\n", cfg.num_landmarks);
    ALTData *alt = alt_preprocess(g, g_rev, cfg.num_landmarks, cfg.seed);
    printf("    Prétraitement ALT : %.1f ms | mémoire %.1f Mo\n",
           alt->prep_time_ms, alt->memory_bytes / (1024.0 * 1024.0));

    CHIndex *ch = NULL;
    if (!cfg.skip_ch) {
        printf("  > CH (contraction complète, edge-difference)\n");
        ch = ch_preprocess(g, g_rev, /*verbose=*/1);
        printf("    Prétraitement CH  : %.2f s | mémoire %.1f Mo\n",
               ch->prep_time_ms / 1000.0,
               ch->memory_bytes / (1024.0 * 1024.0));
        printf("    Raccourcis ajoutés : %d (graphe augmenté = %d arêtes, "
               "facteur ×%.2f)\n",
               ch->num_shortcuts, ch->up->n_edges,
               (double)ch->up->n_edges / g->n_edges);
    }


    /* Génération du lot de requêtes */
    int N = cfg.num_queries;
    srand(cfg.seed);
    int *src = malloc(N * sizeof(int));
    int *tgt = malloc(N * sizeof(int));
    for (int q = 0; q < N; q++) {
        src[q] = rand() % n_nodes;
        tgt[q] = rand() % n_nodes;
    }


    /* Stockage des temps et résultats par requête */
    #define ALGO_N 4
    const char *ALGO_NAMES[ALGO_N] = {"Dijkstra", "A*", "ALT", "CH"};
    double *times[ALGO_N];
    int    *explor[ALGO_N];
    int    *relax[ALGO_N];
    double *dists[ALGO_N];
    for (int a = 0; a < ALGO_N; a++) {
        times[a]  = malloc(N * sizeof(double));
        explor[a] = malloc(N * sizeof(int));
        relax[a]  = malloc(N * sizeof(int));
        dists[a]  = malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) { times[a][i] = 0.0; dists[a][i] = -1.0; }
    }
    Category *cats = malloc(N * sizeof(Category));
    int       valid_q = 0;          /* nb de requêtes où src est connecté à tgt */
    int       cat_count[CAT_N] = {0, 0, 0};


    /*Benchmark*/
    printf("\n=== Benchmark sur %d requêtes ===\n", N);
    int errors = 0;
    int progress = N / 20 > 1 ? N / 20 : 1;

    for (int q = 0; q < N; q++) {
        int s = src[q], t = tgt[q];

        ShortestPathResult r_dij = dijkstra(g, s, t);
        ShortestPathResult r_ast = astar(g, nodes, s, t);
        ShortestPathResult r_alt = alt_query(g, alt, s, t);
        ShortestPathResult r_ch  = (ch ? ch_query(ch, s, t)
                                       : (ShortestPathResult){-1, NULL, 0, 0, 0, 0});

        ShortestPathResult *all[ALGO_N] = {&r_dij, &r_ast, &r_alt, &r_ch};
        for (int a = 0; a < ALGO_N; a++) {
            times[a][q]  = all[a]->time_ms;
            explor[a][q] = all[a]->nodes_explored;
            relax[a][q]  = all[a]->relaxations;
            dists[a][q]  = all[a]->dist;
        }

        /* Vérification : Dijkstra est la référence */
        if (r_dij.dist < 0) {
            cats[q] = CAT_SHORT;
            goto cleanup_q;
        }
        valid_q++;
        cats[q] = classify_distance(r_dij.dist);
        cat_count[cats[q]]++;

        /* Vérification distance et chemin pour chaque algo non-Dijkstra */
        for (int a = 1; a < ALGO_N; a++) {
            if (a == 3 && !ch) continue;
            if (all[a]->dist < 0) {
                fprintf(stderr,
                    "  [ERR] q=%d %s dit inaccessible alors que Dijkstra a trouvé\n",
                    q, ALGO_NAMES[a]);
                errors++; continue;
            }
            if (fabs(all[a]->dist - r_dij.dist) > 1e-3) {
                fprintf(stderr,
                    "  [ERR] q=%d %s dist=%.3f vs Dijkstra=%.3f\n",
                    q, ALGO_NAMES[a], all[a]->dist, r_dij.dist);
                errors++;
            }
            /* Validation chemin */
            if (!validate_path(g, all[a]->path, all[a]->path_len,
                                all[a]->dist, ALGO_NAMES[a])) {
                errors++;
            }
        }
        if (!validate_path(g, r_dij.path, r_dij.path_len,
                             r_dij.dist, "Dijkstra")) {
            errors++;
        }

cleanup_q:
        free(r_dij.path); free(r_ast.path); free(r_alt.path);
        if (ch) free(r_ch.path);

        if ((q + 1) % progress == 0)
            printf("    [%d/%d] requêtes traitées (%.0f%%)\n",
                   q + 1, N, 100.0 * (q + 1) / N);
    }

    if (errors == 0)
        printf("\n  [OK] %d requêtes vérifiées (distances + chemins)\n", valid_q);
    else
        printf("\n  [%d ERREURS] détectées (voir stderr)\n", errors);

    printf("\n=== Récapitulatif global (%d requêtes valides sur %d) ===\n",
           valid_q, N);
    printf("  Algo     |  moy (ms) |  p50 (ms) |  p95 (ms) |   req/s  | nœuds expl. moy | xDijk\n");
    printf("  ---------+-----------+-----------+-----------+----------+-----------------+-------\n");

    /* On calcule les stats sur les requêtes VALIDES uniquement.
     * Pour ça on filtre dans des tableaux temporaires. */
    double *tmp_times = malloc(N * sizeof(double));
    long    sum_explor[ALGO_N] = {0};
    TimeStats stats[ALGO_N];

    for (int a = 0; a < ALGO_N; a++) {
        if (a == 3 && !ch) continue;
        int k = 0;
        for (int q = 0; q < N; q++) {
            if (dists[0][q] < 0) continue;
            tmp_times[k++] = times[a][q];
            sum_explor[a] += explor[a][q];
        }
        stats[a] = compute_stats(tmp_times, k);
    }
    free(tmp_times);

    for (int a = 0; a < ALGO_N; a++) {
        if (a == 3 && !ch) continue;
        double speedup = (stats[a].mean_ms > 0)
                        ? stats[0].mean_ms / stats[a].mean_ms : 0.0;
        printf("  %-8s | %9.3f | %9.3f | %9.3f | %8.0f | %15ld | x%.1f\n",
               ALGO_NAMES[a],
               stats[a].mean_ms, stats[a].p50_ms, stats[a].p95_ms,
               stats[a].throughput,
               valid_q > 0 ? sum_explor[a] / valid_q : 0,
               speedup);
    }


    /* STATS PAR CATÉGORIE (SHORT / MEDIUM / LONG) */
    printf("\n=== Stats par catégorie de distance ===\n");
    printf("  Seuils : SHORT < 20km , MEDIUM 20-80km , LONG > 80km\n");
    printf("  Effectifs : SHORT=%d  MEDIUM=%d  LONG=%d\n",
           cat_count[CAT_SHORT], cat_count[CAT_MEDIUM], cat_count[CAT_LONG]);

    double *tmp = malloc(N * sizeof(double));
    for (int c = 0; c < CAT_N; c++) {
        if (cat_count[c] == 0) {
            printf("\n  [%s] aucune requête dans cette catégorie\n", CAT_NAMES[c]);
            continue;
        }
        printf("\n  [%s] %d requêtes\n", CAT_NAMES[c], cat_count[c]);
        printf("    Algo     |  moy (ms) |  p50 (ms) |  p95 (ms) |   req/s  | nœuds expl.\n");
        printf("    ---------+-----------+-----------+-----------+----------+------------\n");
        for (int a = 0; a < ALGO_N; a++) {
            if (a == 3 && !ch) continue;
            int    k = 0;
            long   sumE = 0;
            for (int q = 0; q < N; q++) {
                if (dists[0][q] < 0)            continue;
                if ((int)cats[q] != c)          continue;
                tmp[k++] = times[a][q];
                sumE += explor[a][q];
            }
            TimeStats st = compute_stats(tmp, k);
            printf("    %-8s | %9.3f | %9.3f | %9.3f | %8.0f | %11ld\n",
                   ALGO_NAMES[a],
                   st.mean_ms, st.p50_ms, st.p95_ms, st.throughput,
                   k > 0 ? sumE / k : 0);
        }
    }
    free(tmp);


    /*COÛTS DE PRÉTRAITEMENT*/
    printf("\n=== Prétraitement et mémoire ===\n");
    printf("  Graphe CSR (G + G' + nodes) : %.1f Mo\n", mem_csr / (1024.0 * 1024.0));
    printf("  ALT : %.1f ms , %.1f Mo (%d landmarks)\n",
           alt->prep_time_ms, alt->memory_bytes / (1024.0 * 1024.0),
           alt->num_landmarks);
    if (ch) {
        printf("  CH  : %.2f s , %.1f Mo (%d raccourcis ; ratio ×%.2f)\n",
               ch->prep_time_ms / 1000.0,
               ch->memory_bytes / (1024.0 * 1024.0),
               ch->num_shortcuts,
               (double)ch->up->n_edges / g->n_edges);
    }


    /* EXPORT CSV */
    if (cfg.csv_out) {
        FILE *fcsv = fopen(cfg.csv_out, "w");
        if (!fcsv) {
            perror(cfg.csv_out);
        } else {
            fprintf(fcsv,
                "query_id,source,target,category,algo,"
                "dist_m,nodes_explored,relaxations,time_ms\n");
            for (int q = 0; q < N; q++) {
                if (dists[0][q] < 0) continue;
                const char *catstr = CAT_NAMES[cats[q]];
                for (int a = 0; a < ALGO_N; a++) {
                    if (a == 3 && !ch) continue;
                    fprintf(fcsv, "%d,%d,%d,%s,%s,%.3f,%d,%d,%.6f\n",
                            q, src[q], tgt[q], catstr, ALGO_NAMES[a],
                            dists[a][q], explor[a][q], relax[a][q],
                            times[a][q]);
                }
            }
            fclose(fcsv);
            printf("\n  [CSV] résultats détaillés écrits dans %s\n", cfg.csv_out);
        }
    }


    /* ── Nettoyage ── */
    for (int a = 0; a < ALGO_N; a++) {
        free(times[a]); free(explor[a]); free(relax[a]); free(dists[a]);
    }
    free(cats); free(src); free(tgt);
    alt_free(alt);
    if (ch) ch_index_free(ch);
    csr_free(g_rev);
    csr_free(g);
    free(nodes);

    printf("\n=== Terminé ===\n");
    return (errors > 0) ? 1 : 0;
}