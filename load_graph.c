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
    fgets(line, sizeof(line), f);           /* skip header */
    while (fgets(line, sizeof(line), f)) count++;
    rewind(f);
    fgets(line, sizeof(line), f);           /* skip header again */

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
    fgets(line, sizeof(line), f);
    while (fgets(line, sizeof(line), f)) count++;
    rewind(f);
    fgets(line, sizeof(line), f);

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
    clock_t t0 = clock();
    

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

    clock_t t1 = clock();
    res.time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
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
    clock_t t0 = clock();

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

    clock_t t1 = clock();
    res.time_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    return res;
}


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 8 — AFFICHAGE ET VÉRIFICATION
 * ═══════════════════════════════════════════════════════════════════════════ */

static void print_result(const char *algo_name,
                          const ShortestPathResult *r,
                          int source, int target) {
    printf("  %-10s | ", algo_name);
    if (r->dist < 0) {
        printf("INACCESSIBLE\n");
        return;
    }
    printf("dist = %8.1f m | "
           "nœuds explorés = %7d | "
           "relaxations = %7d | "
           "temps = %.3f ms\n",
           r->dist, r->nodes_explored, r->relaxations, r->time_ms);
}

/*
 * Vérifie que Dijkstra et A* donnent la même distance.
 * Si ce n'est pas le cas, il y a un bug dans A*.
 */
static void verify_same_dist(const ShortestPathResult *dijk,
                               const ShortestPathResult *astar_res) {
    if (dijk->dist < 0 && astar_res->dist < 0) return;  /* tous deux inaccessibles */
    if (dijk->dist < 0 || astar_res->dist < 0) {
        printf("  [ERREUR] Un algo dit inaccessible, l'autre non !\n");
        return;
    }
    double diff = fabs(dijk->dist - astar_res->dist);
    if (diff > 1e-3) {  /* tolérance 1 mm */
        printf("  [ERREUR] Distances différentes : Dijkstra=%.3f, A*=%.3f (diff=%.6f)\n",
               dijk->dist, astar_res->dist, diff);
    } else {
        printf("  [OK] Distances identiques (diff < 1mm)\n");
    }
}


/* ═══════════════════════════════════════════════════════════════════════════
 * SECTION 9 — PROGRAMME PRINCIPAL
 * ═══════════════════════════════════════════════════════════════════════════ */

int main(int argc, char *argv[]) {
    const char *nodes_path = (argc > 1) ? argv[1] : "nodes.csv";
    const char *edges_path = (argc > 2) ? argv[2] : "edges.csv";

    printf("=== Route Planning — Étapes 2, 3, 4 ===\n\n");

    /* ── Chargement des données ── */
    int      n_nodes;
    HashMap *map;
    printf("[1/3] Lecture de %s ...\n", nodes_path);
    Node *nodes = read_nodes(nodes_path, &n_nodes, &map);
    printf("      %d noeuds chargés\n", n_nodes);

    int      n_edges;
    printf("[2/3] Lecture de %s ...\n", edges_path);
    RawEdge *raw_edges = read_edges(edges_path, map, &n_edges);
    printf("      %d arêtes chargées\n", n_edges);

    printf("[3/3] Construction CSR ...\n");
    CSRGraph *g = build_csr(n_nodes, raw_edges, n_edges);
    printf("      CSR prêt : %d noeuds, %d arêtes\n\n", g->n_nodes, g->n_edges);

    free(raw_edges);
    map_free(map);

    /* ── Vérification basique de la structure CSR ── */
    if (g->row_ptr[g->n_nodes] != g->n_edges) {
        fprintf(stderr, "ERREUR CSR : row_ptr[N]=%d != n_edges=%d\n",
                g->row_ptr[g->n_nodes], g->n_edges);
        return 1;
    }

    /* ── Benchmark sur N_QUERIES paires aléatoires ── */
    /*
     * On génère des paires (source, target) avec une graine fixe (srand(42))
     * pour avoir des résultats reproductibles. C'est important pour le rapport :
     * les chiffres doivent être les mêmes d'une exécution à l'autre.
     */
    #define N_QUERIES 10
    srand(42);  /* graine fixe → résultats reproductibles */

    printf("=== Benchmark sur %d requêtes aléatoires ===\n", N_QUERIES);
    printf("  Algo       | Distance         | Nœuds explorés    | Relaxations    | Temps\n");
    printf("  -----------+------------------+--------------------+----------------+----------\n");

    double total_time_dijk  = 0.0;
    double total_time_astar = 0.0;
    int    total_explored_dijk  = 0;
    int    total_explored_astar = 0;

    for (int q = 0; q < N_QUERIES; q++) {
        int source = rand() % n_nodes;
        int target = rand() % n_nodes;
        printf("\nRequête %2d : nœud %d → nœud %d\n", q+1, source, target);

        /* Dijkstra */
        ShortestPathResult r_dijk = dijkstra(g, source, target);
        print_result("Dijkstra", &r_dijk, source, target);

        /* A* */
        ShortestPathResult r_astar = astar(g, nodes, source, target);
        print_result("A*", &r_astar, source, target);

        /* Vérification : les deux doivent trouver la même distance */
        verify_same_dist(&r_dijk, &r_astar);

        /* Gain A* vs Dijkstra */
        if (r_dijk.nodes_explored > 0) {
            double gain = 100.0 * (1.0 - (double)r_astar.nodes_explored
                                          / r_dijk.nodes_explored);
            printf("  Gain A*   : %.1f%% de nœuds explorés en moins\n", gain);
        }

        total_time_dijk      += r_dijk.time_ms;
        total_time_astar     += r_astar.time_ms;
        total_explored_dijk  += r_dijk.nodes_explored;
        total_explored_astar += r_astar.nodes_explored;

        /* Libérer les chemins alloués */
        free(r_dijk.path);
        free(r_astar.path);
    }

    /* ── Récapitulatif ── */
    printf("\n=== Récapitulatif sur %d requêtes ===\n", N_QUERIES);
    printf("  Dijkstra : temps moyen = %.2f ms | nœuds/requête = %d\n",
           total_time_dijk / N_QUERIES, total_explored_dijk / N_QUERIES);
    printf("  A*       : temps moyen = %.2f ms | nœuds/requête = %d\n",
           total_time_astar / N_QUERIES, total_explored_astar / N_QUERIES);
    double speedup = total_time_dijk / (total_time_astar > 0 ? total_time_astar : 1);
    printf("  Accélération A* vs Dijkstra : x%.2f\n", speedup);

    /* ── Nettoyage ── */
    free(nodes);
    csr_free(g);

    printf("\n=== Terminé ===\n");
    return 0;
}