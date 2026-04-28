/*
 * route_planning.c
 * Étapes 2, 3 et 4 — CSR + Dijkstra + A*
 *
 * Compilation : gcc -O2 -o route_planning route_planning.c -lm
 * Usage       : ./route_planning nodes.csv edges.csv
 *
 * Ce fichier contient tout en un seul endroit :
 *   - Chargement des CSV et construction du graphe CSR  (étape 2)
 *   - Algorithme de Dijkstra avec min-heap              (étape 3)
 *   - Algorithme A* avec heuristique Haversine          (étape 4)
 *   - Benchmark automatique sur 10 paires aléatoires
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>    /* pour mesurer le temps d'exécution */
#include <float.h>   /* pour DBL_MAX (infini initial des distances) */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/* ─── Chrono haute précision (nanoseconde) ─────────────────────────────
 * clock() est trop grossier (10 ms sous Linux). On utilise CLOCK_MONOTONIC
 * pour mesurer les requêtes courtes (< 1 ms) sans les voir comme "0".
 */
static inline double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 1 — STRUCTURES DE DONNÉES
 * ═══════════════════════════════════════════════════════════════════════════ */

/*
 * Node : informations géographiques d'un nœud après renumérotation.
 * On stocke la latitude et longitude pour pouvoir calculer
 * l'heuristique Haversine dans A*.
 */
typedef struct {
    long long osm_id;   /* identifiant OSM original (grand entier) */
    double    lat;      /* latitude  en degrés décimaux            */
    double    lon;      /* longitude en degrés décimaux            */
} Node;

/*
 * CSRGraph : le graphe stocké en format Compressed Sparse Row.
 *
 * Rappel de la structure (voir README) :
 *   Pour lister les voisins du nœud i, on fait :
 *       for (k = row_ptr[i]; k < row_ptr[i+1]; k++)
 *           voisin = adj[k], poids = weights[k]
 *
 *   row_ptr a une taille de n_nodes + 1.
 *   adj et weights ont une taille de n_edges.
 */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 2 — TABLE DE HACHAGE (pour renuméroter les IDs OSM)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Les IDs OSM sont de grands entiers (ex: 3703570445).
 * On ne peut pas les utiliser directement comme indices de tableau.
 * On crée donc une table de hachage qui fait la correspondance :
 *   osm_id (long long) → indice local (int, 0..N-1)
 */

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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 3 — LECTURE DES CSV ET CONSTRUCTION CSR
 * ═══════════════════════════════════════════════════════════════════════════ */

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

/*
 * build_csr : construit la structure CSR en 3 passes.
 *
 * Passe 1 — On compte combien d'arêtes partent de chaque nœud.
 *   row_ptr[u+1]++ pour chaque arête (u→v).
 *   Après la passe 1, row_ptr[i+1] contient le degré sortant de i.
 *
 * Passe 2 — Somme préfixe (prefix sum).
 *   row_ptr[i] = row_ptr[0] + row_ptr[1] + ... + row_ptr[i-1]
 *   Maintenant row_ptr[i] est l'indice de début des arêtes de i dans adj[].
 *
 * Passe 3 — On remplit adj[] et weights[].
 *   On utilise un tableau "cursor" (copie de row_ptr) comme compteur
 *   qui avance au fur et à mesure qu'on place les arêtes.
 */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 4 — MIN-HEAP (FILE DE PRIORITÉ)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Un tas binaire min est un arbre binaire complet stocké dans un tableau.
 * La propriété fondamentale : chaque nœud a une clé <= celles de ses enfants.
 * Donc le minimum est TOUJOURS à l'indice 0.
 *
 * Deux opérations principales :
 *   heap_push : ajoute un élément → le place à la fin, puis "remonte"
 *   heap_pop  : retire le minimum → remplace la racine par le dernier
 *               élément, puis "descend" pour restaurer la propriété du tas
 *
 * Exemple de tas avec les valeurs [1, 3, 2, 7, 5, 4] :
 *         1          ← indice 0 (minimum, toujours en haut)
 *        / \
 *       3   2        ← indices 1, 2
 *      / \ / \
 *     7  5 4         ← indices 3, 4, 5
 *
 * Relations dans le tableau : parent(i) = (i-1)/2
 *                              enfant_gauche(i) = 2*i + 1
 *                              enfant_droit(i)  = 2*i + 2
 */

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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 5 — HEURISTIQUE HAVERSINE (pour A*)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Calcule la distance à vol d'oiseau entre deux points GPS.
 * C'est la distance réelle en ligne droite sur la Terre (sphère).
 *
 * Cette distance est TOUJOURS <= distance réelle sur route.
 * C'est ce qu'on appelle une heuristique "admissible" pour A* :
 * elle ne surestime jamais le coût réel, donc A* reste correct.
 *
 * Si on divisait par la vitesse max du réseau, on obtiendrait
 * une heuristique en secondes (pour un critère temps plutôt que distance).
 * On reste ici en mètres pour la cohérence avec les poids du graphe.
 */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 6 — DIJKSTRA
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Dijkstra trouve le plus court chemin depuis 'source' vers 'target'.
 *
 * IDÉE GÉNÉRALE :
 *   On maintient pour chaque nœud v une "distance courante" dist[v].
 *   Au départ : dist[source] = 0, dist[tous les autres] = +infini.
 *   On utilise un tas min pour toujours traiter en premier
 *   le nœud non-traité avec la plus petite distance connue.
 *
 * INVARIANT DE L'ALGORITHME :
 *   Quand on extrait un nœud u du tas, dist[u] est définitivement optimal.
 *   (C'est pourquoi Dijkstra ne fonctionne qu'avec des poids positifs !)
 *
 * DÉROULEMENT ÉTAPE PAR ÉTAPE :
 *   1. Initialiser dist[source] = 0, tout le reste = DBL_MAX
 *   2. Pousser (source, 0) dans le tas
 *   3. Tant que le tas n'est pas vide :
 *      a. Extraire u = nœud avec la plus petite distance dans le tas
 *      b. Si dist extrait > dist[u], ce nœud est "périmé" → ignorer (*)
 *      c. Si u == target → on a trouvé le plus court chemin, arrêter
 *      d. Pour chaque voisin v de u :
 *         - Calculer new_dist = dist[u] + poids(u,v)
 *         - Si new_dist < dist[v] : mettre à jour dist[v], pousser (v, new_dist)
 *
 * (*) Note sur les "entrées périmées" :
 *   On utilise un tas "paresseux" : au lieu de modifier une entrée existante
 *   dans le tas (opération coûteuse), on ajoute simplement une nouvelle entrée
 *   avec la distance mise à jour. Les anciennes entrées restent dans le tas
 *   mais seront ignorées grâce au test "if (f > dist[u]) continue".
 *
 * RECONSTRUCTION DU CHEMIN :
 *   On maintient un tableau prev[v] = nœud depuis lequel on a atteint v
 *   avec la meilleure distance. À la fin, on remonte de target à source.
 */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 7 — A* (A-ÉTOILE)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * A* est une amélioration de Dijkstra qui utilise une heuristique h(v)
 * pour "guider" la recherche vers la cible.
 *
 * DIFFÉRENCE AVEC DIJKSTRA :
 *   Dijkstra trie les nœuds par dist[v]           (coût depuis source)
 *   A*       trie les nœuds par dist[v] + h(v,t)  (coût estimé total)
 *                                         ↑ heuristique = estimation du reste
 *
 * INTUITION :
 *   Dijkstra explore dans un cercle autour de la source (tous les nœuds
 *   à distance d, puis tous les nœuds à distance d+ε, etc.).
 *   A* explore dans une ellipse orientée vers la cible : il "sait"
 *   approximativement où aller et évite d'explorer dans la mauvaise direction.
 *
 * POURQUOI H(V) = HAVERSINE EST CORRECT ?
 *   Pour que A* trouve bien le plus court chemin, h doit être "admissible" :
 *       h(v, t) <= vraie distance(v, t)    pour tout v
 *   La distance à vol d'oiseau est toujours <= distance réelle sur route.
 *   Donc Haversine est admissible. ✓
 *
 * LE CODE EST QUASI-IDENTIQUE À DIJKSTRA, avec une seule différence :
 *   La clé du tas f = dist[v] + h(v, target) au lieu de dist[v].
 *
 * NOTATION CLASSIQUE EN INFORMATIQUE :
 *   g(v) = dist[v]         = coût réel depuis la source
 *   h(v) = haversine(v, t) = estimation du coût restant
 *   f(v) = g(v) + h(v)     = estimation du coût total du chemin passant par v
 */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 8 — (libre — l'ancien print_result et verify_same_dist ont été
 * remplacés par les utilitaires du benchmark, voir section 16)
 * ═══════════════════════════════════════════════════════════════════════════ */


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 9 — DIJKSTRA "ALL-NODES" + GRAPHE INVERSE
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Pour ALT et CH, on a besoin de deux choses qu'on n'avait pas :
 *
 *   1) Un Dijkstra qui calcule les distances depuis une source vers TOUS
 *      les nœuds atteignables (pas d'arrêt anticipé sur une cible).
 *      → c'est ce qu'on utilise pour pré-calculer les distances aux
 *        landmarks (ALT) et pour la recherche de témoins (CH).
 *
 *   2) Le graphe inverse G' : si (u → v) est dans G, alors (v → u) est
 *      dans G'. Indispensable pour :
 *        - calculer d(v, ℓ) (distance VERS un landmark) en faisant un
 *          Dijkstra depuis ℓ sur G',
 *        - faire la recherche en arrière dans le bidirectionnel CH.
 */

/* ─── 9.1 Dijkstra qui remplit le tableau dist[] depuis source ────────── */
/*
 * Variante "pas d'arrêt" du Dijkstra de la section 6.
 * À la sortie : dist[v] = plus court chemin source → v (ou DBL_MAX si v
 * inaccessible). Aucun chemin n'est reconstruit ici.
 *
 * Complexité : O((N + M) log N).
 */
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


/* ─── 9.2 Construction du graphe inverse ──────────────────────────────── */
/*
 * À partir d'un graphe CSR direct, on construit son graphe inverse.
 * L'astuce : on part d'une liste de RawEdge inversées (v, u, w),
 * puis on appelle build_csr() qui sait déjà construire un CSR.
 */
static CSRGraph *build_reverse_csr(const CSRGraph *g) {
    int n = g->n_nodes;
    int m = g->n_edges;

    RawEdge *rev = malloc(m * sizeof(RawEdge));
    int idx = 0;
    for (int u = 0; u < n; u++) {
        for (int k = g->row_ptr[u]; k < g->row_ptr[u + 1]; k++) {
            rev[idx].u = g->adj[k];   /* l'arc inverse part de v */
            rev[idx].v = u;
            rev[idx].w = g->weights[k];
            idx++;
        }
    }
    CSRGraph *gr = build_csr(n, rev, m);
    free(rev);
    return gr;
}


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 10 — ALT : A* AVEC LANDMARKS (PRÉTRAITEMENT)
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * IDÉE
 *   A* avec Haversine est limité : la distance à vol d'oiseau sous-estime
 *   très fortement le vrai coût routier (autoroutes obliques, montagnes,
 *   contournements...). Du coup l'heuristique reste faible et A* explore
 *   encore beaucoup de nœuds inutiles.
 *
 *   ALT utilise un PRÉ-CALCUL de distances exactes vers/depuis quelques
 *   nœuds repères (landmarks), et combine ça avec l'inégalité triangulaire :
 *
 *       d(v, t)  ≥  d(ℓ, t)  − d(ℓ, v)        (forward, depuis landmark)
 *       d(v, t)  ≥  d(v, ℓ)  − d(t, ℓ)        (backward, vers landmark)
 *
 *   On prend le max sur tous les landmarks ℓ et sur les deux variantes :
 *
 *       h_ALT(v, t) = max_ℓ max( d(ℓ, t)−d(ℓ, v),
 *                                 d(v, ℓ)−d(t, ℓ),
 *                                 0 )
 *
 *   Cette heuristique reste ADMISSIBLE (le max d'admissibles l'est).
 *   Et comme elle est plus serrée que la haversine, A* explore moins.
 *
 *   Note : le sujet ne donne que la version "forward". On ajoute la
 *   "backward" car le graphe routier est dirigé (sens uniques) : la
 *   version bidirectionnelle est nettement plus informative.
 *
 * COÛTS
 *   - Prétraitement : 2K Dijkstra "all-nodes" (K = nombre de landmarks),
 *     soit O(K · (N + M) log N).
 *   - Mémoire       : 2K · N · 8 octets.
 *   - Requête       : un A* normal mais avec h_ALT au lieu de Haversine.
 *
 * SÉLECTION DES LANDMARKS — STRATÉGIE "FARTHEST"
 *   Choisir K landmarks de façon à les disperser au maximum :
 *     ℓ_0 : un nœud aléatoire
 *     ℓ_i : le nœud le plus éloigné des landmarks déjà choisis
 *           (au sens : max sur les v du min sur les ℓ_j déjà pris de d(ℓ_j, v))
 *
 *   Bien meilleur qu'un tirage aléatoire, presque aussi bon que
 *   les méthodes plus sophistiquées (avoid, planar...) selon la littérature.
 */

typedef struct {
    int       num_landmarks;
    int      *landmarks;     /* tableau des indices des landmarks (taille K) */
    double  **dist_from;     /* dist_from[i][v] = d(landmark_i, v) — forward  */
    double  **dist_to;       /* dist_to[i][v]   = d(v, landmark_i) — backward */
    double    prep_time_ms;  /* temps de prétraitement total                  */
    size_t    memory_bytes;  /* mémoire utilisée (estimation)                 */
} ALTData;


/*
 * Sélection "farthest" des landmarks, RESTREINTE à la plus grande
 * composante connexe accessible.
 *
 * Pourquoi cette restriction ?
 *   Un graphe routier OSM contient souvent de petites composantes
 *   isolées (parkings privés, enclaves mal connectées, erreurs de
 *   tagging). Si on tire le premier landmark dans une de ces micro-
 *   composantes (ex. 10 nœuds), il ne peut atteindre que ces 10 nœuds
 *   et l'heuristique ALT vaut 0 partout ailleurs → ALT ≡ Dijkstra.
 *
 * Algorithme :
 *   1. BFS depuis un nœud "central" (on prend celui de plus haut degré
 *      sortant comme heuristique simple — il est probablement dans la
 *      grande composante). On marque tous les nœuds atteints.
 *   2. On ne pioche les landmarks que parmi les nœuds atteints.
 */
static void select_landmarks_farthest(const CSRGraph *g, int K,
                                      int *landmarks, unsigned seed) {
    int n = g->n_nodes;

    /* ── Étape 1 : trouver la plus grande composante (par BFS) ── */
    /* Heuristique : on part du nœud de plus haut degré sortant. Sur un
     * graphe routier réel, ce nœud est forcément dans la grande
     * composante (carrefour important). */
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


/*
 * Prétraitement complet ALT :
 *   - sélectionne K landmarks (farthest)
 *   - calcule dist_from[i][v] = d(ℓ_i, v) (Dijkstra sur le graphe direct)
 *   - calcule dist_to[i][v]   = d(v, ℓ_i) (Dijkstra sur le graphe inverse)
 */
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


/* ─── Heuristique ALT, à plat et inlinable ───────────────────────────────
 * Pour une requête (s, t) figée, on précalcule un tableau de constantes
 *   C_fwd[i] = dist_from[i][t]
 *   C_bwd[i] = dist_to[i][t]
 * Du coup, l'évaluation de h(v) ne fait qu'un parcours linéaire des K
 * landmarks, sans rappeler aucune fonction lourde. C'est important :
 * h(v) est appelé une fois par relaxation, donc des millions de fois.
 */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 11 — ALT : LA REQUÊTE
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Code identique à A* (section 7), avec deux différences :
 *   - On précalcule C_fwd[i] = d(ℓ_i, t) et C_bwd[i] = d(t, ℓ_i)
 *     (constantes pour toute la requête).
 *   - On remplace haversine() par alt_h(v, C_fwd, C_bwd).
 *
 * Le reste — file de priorité, vérification d'entrée périmée, relaxation,
 * reconstruction de chemin — est strictement le même.
 */
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
        if (cur.f_score > g_score[u] + h_u) continue;  /* périmé */

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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 12 — CONTRACTION HIERARCHIES (CH) : STRUCTURES DE DONNÉES
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * VUE D'ENSEMBLE
 *   CH est une approche très différente de Dijkstra / A* / ALT.
 *   On accepte de payer cher en prétraitement pour obtenir des requêtes
 *   ultra-rapides. Idéal quand on doit répondre à BEAUCOUP de requêtes.
 *
 *   Le prétraitement CH se déroule ainsi :
 *     1. On choisit un ordre de contraction sur les nœuds (du moins
 *        important au plus important).
 *     2. On contracte les nœuds un par un, dans cet ordre. "Contracter"
 *        un nœud v signifie : pour chaque paire (u, w) telle que u→v→w,
 *        on regarde si ce passage par v est un plus court chemin unique.
 *        Si oui, on ajoute un RACCOURCI u→w (weight = w(u,v)+w(v,w)).
 *     3. Le graphe final = arêtes originales + raccourcis.
 *
 *   À la requête, on fait un Dijkstra BIDIRECTIONNEL "qui ne descend
 *   jamais" : on ne suit que les arêtes qui montent dans la hiérarchie
 *   (level[u] < level[v]).
 *
 * STRUCTURES DYNAMIQUES
 *   Pendant le prétraitement, on ajoute des arêtes au fur et à mesure :
 *   le CSR (statique) ne convient plus. On utilise donc des listes
 *   d'adjacence dynamiques (un tableau par nœud, qui grandit avec realloc).
 *   À la fin de la phase de contraction, on REPLIE tout en CSR pour la
 *   phase requête (qui exige de la rapidité maximale).
 *
 *   Chaque arête CH stocke :
 *     - to     : nœud cible
 *     - w      : poids
 *     - via    : -1 si arête originale, sinon indice du nœud "milieu"
 *                contracté qui a généré ce raccourci. C'est ce qui permet
 *                la reconstruction du vrai chemin (unpacking).
 *
 * REMARQUE : on stocke aussi les arêtes ENTRANTES (in[v]). C'est nécessaire
 * pour, lors de la contraction de v, énumérer les paires (u, v) → (v, w)
 * sans devoir parcourir tout le graphe. C'est une duplication mémoire,
 * mais limitée à la phase de prétraitement (libérée ensuite).
 */

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


/* ───────────────────────────────────────────────────────────────────────
 * Buffer global réutilisable pour la recherche de témoins.
 * On utilise l'astuce des "epochs" pour éviter de réinitialiser un
 * tableau de taille N à chaque mini-Dijkstra.
 *
 *   dist_epoch[v] == cur_epoch  →  dist_value[v] est valide
 *   sinon                        →  dist_value[v] est considéré comme +∞
 *
 * On incrémente cur_epoch entre deux recherches : tous les anciens
 * marquages deviennent automatiquement "périmés" sans avoir besoin de
 * boucle de reset.
 * ─────────────────────────────────────────────────────────────────────── */
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 13 — CH : RECHERCHE DE TÉMOINS, CONTRACTION, ORDRE
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * 13.1 RECHERCHE DE TÉMOINS
 *   Pour décider si u→v→w doit donner un raccourci :
 *     - On calcule dCH = w(u,v) + w(v,w).
 *     - On cherche s'il existe un autre chemin de u à w, n'utilisant ni v
 *       ni les nœuds déjà contractés, de longueur ≤ dCH.
 *     - Si oui (TÉMOIN trouvé) : pas de raccourci.
 *     - Sinon : on doit ajouter le raccourci u→w avec poids dCH.
 *
 *   Contrainte d'efficacité : on LIMITE la recherche
 *     - en distance (on ne dépasse pas dCH),
 *     - en nombre de "hops" (on ne va pas trop loin en nombre d'arêtes).
 *   Si la recherche se termine sans avoir prouvé l'existence d'un témoin,
 *   on est conservateur et on ajoute le raccourci. Cela peut produire
 *   QUELQUES raccourcis superflus, mais la correction n'est pas affectée.
 *
 * 13.2 OPTIMISATION — DIJKSTRA "ONE-TO-MANY" PAR SOURCE
 *   Quand on contracte v, pour un même prédécesseur u (in-neighbor),
 *   il peut exister plusieurs successeurs w. Plutôt que de relancer une
 *   recherche par cible, on fait UNE SEULE recherche depuis u, et on lit
 *   ensuite dist[w] pour chaque w. Gain net.
 *
 * 13.3 EDGE DIFFERENCE
 *   Pour décider QUEL nœud contracter en premier, on calcule, pour chaque
 *   candidat v, le nombre de raccourcis qui SERAIENT ajoutés s'il était
 *   contracté maintenant. La priorité (clé du tas) est :
 *
 *      ED(v) = (#raccourcis simulés) − (in_degree(v) + out_degree(v))
 *
 *   Plus c'est petit (voire négatif), plus on préfère contracter v tôt.
 *   On utilise un tas avec mises à jour PARESSEUSES : quand on extrait v,
 *   on recalcule sa priorité ; si ce n'est plus le minimum, on le repousse
 *   et on continue.
 */

/* Limites de la recherche de témoins (à régler selon le graphe) */
#define WITNESS_HOP_LIMIT     5      /* nombre max d'arêtes du témoin */
#define WITNESS_NODE_LIMIT  500      /* nombre max de pop avant abandon */


/*
 * Dijkstra borné depuis u, sur le graphe non-contracté, en évitant le
 * nœud forbidden (= v en cours de contraction).
 *   - max_dist     : on ne pousse jamais une distance > max_dist.
 *   - hop_limit    : on s'arrête de relâcher à partir de hop_limit hops.
 *   - node_limit   : on extrait au plus node_limit nœuds.
 *
 * À la sortie, witness_get(wb, w) renvoie soit la distance trouvée pour w,
 * soit DBL_MAX si on ne l'a pas atteint dans la limite.
 *
 * Note : on ne suit que les arêtes vers des nœuds NON contractés. Sinon,
 * un témoin pourrait passer par un raccourci... ce qui est en fait une
 * bonne chose ! Donc on autorise les nœuds non-contractés y compris ceux
 * qu'on contractera plus tard, mais pas ceux DÉJÀ contractés (hors hierarchie).
 */
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
 *
 * Algorithme :
 *   pour chaque u ∈ in_neighbors(v), non contracté :
 *     pour chaque w ∈ out_neighbors(v), non contracté, w ≠ u :
 *       dCH = w(u,v) + w(v,w)
 *     max_dCH = max sur tous les w
 *     witness_search(u, forbidden=v, max_dist=max_dCH)
 *     pour chaque w :
 *       si witness_get(w) > w(u,v) + w(v,w) : raccourci nécessaire
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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 14 — CH : BOUCLE DE CONTRACTION + PLIAGE EN CSR
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * 14.1 BOUCLE PRINCIPALE
 *   On gère l'ordre de contraction avec une file de priorité (tas min)
 *   sur l'edge difference. Mises à jour PARESSEUSES :
 *     - Quand on extrait v du tas, on recalcule ED(v).
 *     - Si ED(v) > priorité minimale courante du tas → v n'est plus le
 *       meilleur candidat, on le repousse avec sa nouvelle priorité.
 *     - Sinon → on le contracte.
 *     - Après contraction, on repousse les voisins actifs avec une
 *       nouvelle ED (les anciennes entrées resteront dans le tas, mais
 *       seront détectées comme périmées à leur tour).
 *
 * 14.2 PLIAGE EN CSR
 *   Une fois la contraction terminée, on recopie tout le graphe
 *   augmenté (arêtes originales + raccourcis) dans une structure CSR
 *   "lourde" qui inclut un champ via[] supplémentaire.
 *   On construit aussi son inverse pour la recherche en arrière.
 *
 *   Lors de la requête, on filtre dynamiquement les arêtes selon la
 *   condition CH : on ne suit que les arêtes (u→v) avec level[u]<level[v].
 *   C'est un test "if" par arête en plus, mais ça évite de gérer deux
 *   CSR séparés (up et down).
 */

/*
 * CSRGraphCH : variante du CSR qui ajoute un champ via par arête.
 * via[k] = -1 si l'arête est originale, sinon indice du nœud milieu
 * contracté qui a généré ce raccourci. Indispensable pour l'unpacking.
 */
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
        if (ch->contracted[v]) continue;       /* entrée périmée */

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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 15 — CH : REQUÊTE BIDIRECTIONNELLE ET UNPACKING
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * 15.1 RECHERCHE BIDIRECTIONNELLE
 *   On lance DEUX Dijkstras en parallèle :
 *     - "forward"  depuis source, sur idx->up,     sens normal,
 *     - "backward" depuis target, sur idx->up_rev, sens normal aussi.
 *   DANS LES DEUX, on ne suit que les arêtes qui montent dans la
 *   hiérarchie : level[destination] > level[origine].
 *
 *   À chaque itération, on pop le nœud de plus faible f_score parmi les
 *   deux tas. On entretient une variable best_meet = min sur les nœuds
 *   visités v de (dist_fwd[v] + dist_bwd[v]).
 *
 *   ARRÊT : on s'arrête dès que la plus petite clé restante (parmi les
 *   deux tas) est ≥ best_meet. À ce moment, plus aucune amélioration
 *   n'est possible.
 *
 * 15.2 UNPACKING
 *   Une fois la distance et le nœud de jonction m connus, on reconstruit
 *   le chemin réel :
 *     - forward part  : m → ... → s (en remontant prev_fwd), inversée.
 *     - backward part : m → ... → t (en suivant prev_bwd directement).
 *   On obtient une suite d'arêtes (a → b) du graphe AUGMENTÉ, qui peuvent
 *   être :
 *     - originales (via = -1) → on les émet telles quelles,
 *     - des raccourcis (via = mid) → on remplace récursivement par
 *       (a → mid) puis (mid → b).
 *
 *   Les sous-arêtes sont retrouvées par scan dans idx->up.
 */


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
    int    *vf = malloc(n * sizeof(int));   /* via de l'arête prev_fwd[v] → v */
    int    *vb = malloc(n * sizeof(int));   /* via de l'arête v → prev_bwd[v] */
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

    /* ─── Reconstruction du chemin avec unpacking ─── */

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


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 16 — UTILITAIRES POUR LE BENCHMARK
 * ═══════════════════════════════════════════════════════════════════════════
 *
 * Ce qu'on met ici (en plus du main) :
 *   - parsing de la ligne de commande (--landmarks, --queries, --csv, --seed)
 *   - validation d'un chemin reconstruit (vérifie que la séquence est un
 *     vrai chemin du graphe et que la somme des poids matche la distance)
 *   - calcul de statistiques sur des temps : moyenne, p50, p95
 *   - classification des requêtes en SHORT / MEDIUM / LONG selon la distance
 *     Dijkstra (référence) — réponse à la question "le gain est-il stable ?"
 *   - export CSV des résultats détaillés (1 ligne par requête × algo)
 */


/* ─── 16.1 Configuration de la ligne de commande ───────────────────────── */

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


/* ─── 16.2 Validation d'un chemin reconstruit ──────────────────────────── */
/*
 * Vérifie que (path[0], path[1], ..., path[len-1]) est un vrai chemin
 * dans le graphe G, et que la somme des poids sur ce chemin est
 * (à tolérance près) égale à expected_dist.
 *
 * Ce check est essentiel pour CH : si un raccourci est mal "unpacké",
 * la distance peut être correcte mais le chemin retourné peut être
 * faux (sauter des nœuds, ou utiliser une arête inexistante du graphe
 * original). On ne s'en rendrait pas compte en regardant juste la dist.
 *
 * Renvoie 1 si OK, 0 sinon (et écrit la cause sur stderr).
 */
static int validate_path(const CSRGraph *g,
                          const int *path, int len,
                          double expected_dist,
                          const char *who) {
    if (len <= 0) {
        fprintf(stderr, "  [ERR] %s : chemin vide\n", who); return 0;
    }
    if (len == 1) {
        /* source == target : ok ssi expected_dist ~= 0 */
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
        /* Cherche la meilleure arête u→v dans G (on prend la plus courte
         * en cas de doublons, peu probable ici) */
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


/* ─── 16.3 Statistiques (moyenne, p50, p95, débit) ─────────────────────── */

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
 * Modifie le tableau (qsort en place) — c'est OK car on l'a déjà
 * exploité par ailleurs. Si non, faire une copie au préalable.
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
    /* Percentile par "nearest-rank" : indice = ceil(p/100 × n) − 1 */
    int i50 = (int)ceil(0.50 * n) - 1; if (i50 < 0) i50 = 0;
    int i95 = (int)ceil(0.95 * n) - 1; if (i95 < 0) i95 = 0;
    s.p50_ms = times_ms[i50];
    s.p95_ms = times_ms[i95];
    s.throughput = (sum > 0) ? (1000.0 * n / sum) : 0.0;
    return s;
}


/* ─── 16.4 Classification des requêtes par distance ────────────────────── */
/*
 * On classe les requêtes en 3 catégories selon la distance Dijkstra :
 *   SHORT  : distance < 5 km
 *   MEDIUM : entre 5 km et 30 km
 *   LONG   : > 30 km
 *
 * Ces seuils sont raisonnables pour un graphe routier régional.
 * À adapter pour une ville (mettre ex. 1 km / 5 km) ou un pays
 * (ex. 50 km / 200 km).
 */
typedef enum { CAT_SHORT = 0, CAT_MEDIUM = 1, CAT_LONG = 2, CAT_N = 3 } Category;
static const char *CAT_NAMES[CAT_N] = {"SHORT", "MEDIUM", "LONG"};

static Category classify_distance(double dist_m) {
    if (dist_m < 5000.0)  return CAT_SHORT;
    if (dist_m < 30000.0) return CAT_MEDIUM;
    return CAT_LONG;
}


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 17 — PROGRAMME PRINCIPAL
 * ═══════════════════════════════════════════════════════════════════════════ */

int main(int argc, char *argv[]) {
    Config cfg = parse_args(argc, argv);

    printf("=== Route Planning : Dijkstra | A* | ALT | CH ===\n");
    printf("    nodes=%s edges=%s\n", cfg.nodes_path, cfg.edges_path);
    printf("    landmarks=%d queries=%d seed=%u%s\n",
           cfg.num_landmarks, cfg.num_queries, cfg.seed,
           cfg.skip_ch ? " (CH désactivé)" : "");
    if (cfg.csv_out) printf("    csv=%s\n", cfg.csv_out);
    printf("\n");

    /* ── [1/4] Chargement et CSR ── */
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


    /* ── [4/4] Prétraitement ALT et CH ── */
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


    /* ── Génération du lot de requêtes ── */
    int N = cfg.num_queries;
    srand(cfg.seed);
    int *src = malloc(N * sizeof(int));
    int *tgt = malloc(N * sizeof(int));
    for (int q = 0; q < N; q++) {
        src[q] = rand() % n_nodes;
        tgt[q] = rand() % n_nodes;
    }


    /* ── Stockage des temps et résultats par requête ── */
    /* Indexation : [algo][query]. ALGO = 0 Dijkstra, 1 A*, 2 ALT, 3 CH. */
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


    /* ── Benchmark ── */
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
            cats[q] = CAT_SHORT;     /* arbitraire : non comptabilisé */
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
        /* Validation chemin Dijkstra aussi (sanity check) */
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


    /* ─────────────────────────────────────────────────────────────────
     * RÉCAPITULATIF GLOBAL : moyenne / p50 / p95 / débit / nœuds
     * ───────────────────────────────────────────────────────────────── */
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
            if (dists[0][q] < 0) continue;     /* Dijkstra inaccessible : skip */
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


    /* ─────────────────────────────────────────────────────────────────
     * STATS PAR CATÉGORIE (SHORT / MEDIUM / LONG)
     * ───────────────────────────────────────────────────────────────── */
    printf("\n=== Stats par catégorie de distance ===\n");
    printf("  Seuils : SHORT < 5km , MEDIUM 5-30km , LONG > 30km\n");
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


    /* ─────────────────────────────────────────────────────────────────
     * COÛTS DE PRÉTRAITEMENT
     * ───────────────────────────────────────────────────────────────── */
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


    /* ─────────────────────────────────────────────────────────────────
     * EXPORT CSV (1 ligne par requête × algo)
     * ───────────────────────────────────────────────────────────────── */
    if (cfg.csv_out) {
        FILE *fcsv = fopen(cfg.csv_out, "w");
        if (!fcsv) {
            perror(cfg.csv_out);
        } else {
            fprintf(fcsv,
                "query_id,source,target,category,algo,"
                "dist_m,nodes_explored,relaxations,time_ms\n");
            for (int q = 0; q < N; q++) {
                if (dists[0][q] < 0) continue;     /* on n'écrit que les valides */
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