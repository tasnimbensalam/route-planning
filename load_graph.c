/*
 * load_graph.c
 * Étape 2 — Lecture des CSV + représentation CSR
 *
 * Lecture de nodes.csv et edges.csv produits par l'étape Python,
 * renumérotation contiguë des noeuds (0..N-1),
 * puis construction de la représentation CSR du graphe orienté pondéré.
 *
 * Compilation : gcc -O2 -o load_graph load_graph.c -lm
 * Usage       : ./load_graph nodes.csv edges.csv
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/* ─────────────────────────────────────────────────────────────────────────
 * Structures
 * ───────────────────────────────────────────────────────────────────────── */

/* Informations géographiques d'un noeud (après renumérotation) */
typedef struct {
    long long osm_id;   /* identifiant OSM original            */
    double    lat;
    double    lon;
} Node;

/*
 * Représentation CSR (Compressed Sparse Row) du graphe orienté.
 *
 * Pour le noeud i, ses voisins se trouvent dans :
 *   adj[ row_ptr[i] .. row_ptr[i+1] - 1 ]
 * avec les poids correspondants dans :
 *   weights[ row_ptr[i] .. row_ptr[i+1] - 1 ]
 */
typedef struct {
    int     n_nodes;    /* nombre de noeuds                    */
    int     n_edges;    /* nombre d'arêtes                     */
    int    *row_ptr;    /* taille n_nodes + 1                  */
    int    *adj;        /* taille n_edges  — noeud destination */
    double *weights;    /* taille n_edges  — poids (mètres)    */
} CSRGraph;

/* ─────────────────────────────────────────────────────────────────────────
 * Table de hachage minimaliste  osm_id -> indice local
 * (adressage ouvert, facteur de charge ≤ 0.7)
 * ───────────────────────────────────────────────────────────────────────── */
#define HASH_EMPTY  (-1LL)

typedef struct {
    long long key;   /* osm_id, HASH_EMPTY si case libre */
    int       val;   /* indice local 0..N-1              */
} HashEntry;

typedef struct {
    HashEntry *table;
    int        cap;    /* puissance de 2                  */
    int        size;
} HashMap;

static HashMap *map_create(int capacity) {
    /* arrondir à la puissance de 2 supérieure */
    int cap = 1;
    while (cap < capacity * 2) cap <<= 1;

    HashMap *m = malloc(sizeof(HashMap));
    m->table   = malloc(cap * sizeof(HashEntry));
    m->cap     = cap;
    m->size    = 0;
    for (int i = 0; i < cap; i++) m->table[i].key = HASH_EMPTY;
    return m;
}

static inline int map_slot(HashMap *m, long long key) {
    /* FNV-1a simplifié pour entiers 64 bits */
    unsigned long long h = (unsigned long long)key * 2654435761ULL;
    return (int)(h & (unsigned)(m->cap - 1));
}

static void map_insert(HashMap *m, long long key, int val) {
    int i = map_slot(m, key);
    while (m->table[i].key != HASH_EMPTY && m->table[i].key != key)
        i = (i + 1) & (m->cap - 1);
    m->table[i].key = key;
    m->table[i].val = val;
    m->size++;
}

/* Retourne l'indice local, ou -1 si absent */
static int map_get(HashMap *m, long long key) {
    int i = map_slot(m, key);
    while (m->table[i].key != HASH_EMPTY) {
        if (m->table[i].key == key) return m->table[i].val;
        i = (i + 1) & (m->cap - 1);
    }
    return -1;
}

static void map_free(HashMap *m) {
    free(m->table);
    free(m);
}

/* ─────────────────────────────────────────────────────────────────────────
 * Lecture de nodes.csv
 *   Format attendu : id,lat,lon  (première ligne = en-tête)
 * ───────────────────────────────────────────────────────────────────────── */
static Node *read_nodes(const char *path, int *out_n, HashMap **out_map) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    /* Compter les lignes (hors en-tête) pour allouer */
    int count = 0;
    char line[256];
    fgets(line, sizeof(line), f); /* skip header */
    while (fgets(line, sizeof(line), f)) count++;
    rewind(f);
    fgets(line, sizeof(line), f); /* skip header again */

    Node    *nodes = malloc(count * sizeof(Node));
    HashMap *map   = map_create(count);
    int idx = 0;

    while (fgets(line, sizeof(line), f)) {
        long long id;
        double lat, lon;
        if (sscanf(line, "%lld,%lf,%lf", &id, &lat, &lon) != 3) continue;
        nodes[idx].osm_id = id;
        nodes[idx].lat    = lat;
        nodes[idx].lon    = lon;
        map_insert(map, id, idx);
        idx++;
    }
    fclose(f);

    *out_n   = idx;
    *out_map = map;
    return nodes;
}

/* ─────────────────────────────────────────────────────────────────────────
 * Lecture de edges.csv  (deux passes pour construire le CSR)
 *   Format attendu : u,v,distance_m  (première ligne = en-tête)
 * ───────────────────────────────────────────────────────────────────────── */

/* Arête brute (avant CSR) */
typedef struct { int u, v; double w; } RawEdge;

static RawEdge *read_edges(const char *path, HashMap *map,
                            int *out_m, int n_nodes) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    /* Compter les lignes valides */
    int count = 0;
    char line[256];
    fgets(line, sizeof(line), f); /* skip header */
    while (fgets(line, sizeof(line), f)) count++;
    rewind(f);
    fgets(line, sizeof(line), f);

    RawEdge *edges = malloc(count * sizeof(RawEdge));
    int idx = 0, skipped = 0;

    while (fgets(line, sizeof(line), f)) {
        long long u_osm, v_osm;
        double w;
        if (sscanf(line, "%lld,%lld,%lf", &u_osm, &v_osm, &w) != 3) continue;

        int u = map_get(map, u_osm);
        int v = map_get(map, v_osm);
        if (u < 0 || v < 0) { skipped++; continue; }

        edges[idx].u = u;
        edges[idx].v = v;
        edges[idx].w = w;
        idx++;
    }
    fclose(f);

    if (skipped)
        fprintf(stderr, "  [warn] %d arêtes ignorées (noeuds introuvables)\n",
                skipped);

    *out_m = idx;
    return edges;
}

/* ─────────────────────────────────────────────────────────────────────────
 * Construction de la représentation CSR
 *
 * Algorithme en 3 passes :
 *   1. Compter le degré sortant de chaque noeud → row_ptr (décalé de 1)
 *   2. Préfixe cumulé → row_ptr définitif
 *   3. Remplir adj[] et weights[] en utilisant une copie temporaire de row_ptr
 * ───────────────────────────────────────────────────────────────────────── */
static CSRGraph *build_csr(int n, RawEdge *edges, int m) {
    CSRGraph *g = malloc(sizeof(CSRGraph));
    g->n_nodes  = n;
    g->n_edges  = m;
    g->row_ptr  = calloc(n + 1, sizeof(int));
    g->adj      = malloc(m * sizeof(int));
    g->weights  = malloc(m * sizeof(double));

    /* Passe 1 : degrés sortants */
    for (int i = 0; i < m; i++)
        g->row_ptr[edges[i].u + 1]++;

    /* Passe 2 : somme préfixe */
    for (int i = 1; i <= n; i++)
        g->row_ptr[i] += g->row_ptr[i - 1];

    /* Passe 3 : remplissage (curseur temporaire) */
    int *cursor = malloc(n * sizeof(int));
    memcpy(cursor, g->row_ptr, n * sizeof(int));

    for (int i = 0; i < m; i++) {
        int u   = edges[i].u;
        int pos = cursor[u]++;
        g->adj[pos]     = edges[i].v;
        g->weights[pos] = edges[i].w;
    }
    free(cursor);
    return g;
}

/* ─────────────────────────────────────────────────────────────────────────
 * Libération mémoire
 * ───────────────────────────────────────────────────────────────────────── */
static void csr_free(CSRGraph *g) {
    free(g->row_ptr);
    free(g->adj);
    free(g->weights);
    free(g);
}

/* ─────────────────────────────────────────────────────────────────────────
 * Vérifications rapides
 * ───────────────────────────────────────────────────────────────────────── */
static void validate(const CSRGraph *g, const Node *nodes) {
    /* 1. Cohérence row_ptr */
    if (g->row_ptr[g->n_nodes] != g->n_edges) {
        fprintf(stderr, "ERREUR : row_ptr[N] = %d ≠ n_edges = %d\n",
                g->row_ptr[g->n_nodes], g->n_edges);
        exit(1);
    }

    /* 2. Poids positifs, indices valides */
    int neg = 0, oob = 0;
    for (int i = 0; i < g->n_edges; i++) {
        if (g->weights[i] < 0)          neg++;
        if (g->adj[i] < 0 || g->adj[i] >= g->n_nodes) oob++;
    }
    if (neg) fprintf(stderr, "  [warn] %d arêtes à poids négatif\n", neg);
    if (oob) fprintf(stderr, "  [warn] %d arêtes hors bornes\n",     oob);

    /* 3. Exemple : voisins du noeud 0 */
    printf("\n── Exemple : voisins du noeud 0 (osm_id=%lld, lat=%.5f, lon=%.5f) ──\n",
           nodes[0].osm_id, nodes[0].lat, nodes[0].lon);
    int deg0 = g->row_ptr[1] - g->row_ptr[0];
    printf("  degré sortant = %d\n", deg0);
    int show = deg0 < 5 ? deg0 : 5;
    for (int k = g->row_ptr[0]; k < g->row_ptr[0] + show; k++)
        printf("  -> noeud %d  (%.1f m)\n", g->adj[k], g->weights[k]);
    if (deg0 > show) printf("  ... (%d autres)\n", deg0 - show);
}

/* ─────────────────────────────────────────────────────────────────────────
 * Programme principal
 * ───────────────────────────────────────────────────────────────────────── */
int main(int argc, char *argv[]) {
    const char *nodes_path = (argc > 1) ? argv[1] : "nodes.csv";
    const char *edges_path = (argc > 2) ? argv[2] : "edges.csv";

    printf("=== Étape 2 : Chargement du graphe + construction CSR ===\n\n");

    /* ── 1. Noeuds ── */
    int      n_nodes;
    HashMap *map;
    printf("[1/3] Lecture de %s …\n", nodes_path);
    Node *nodes = read_nodes(nodes_path, &n_nodes, &map);
    printf("      %d noeuds chargés (renumérotés 0..%d)\n", n_nodes, n_nodes - 1);

    /* ── 2. Arêtes ── */
    int      n_edges;
    printf("[2/3] Lecture de %s …\n", edges_path);
    RawEdge *edges = read_edges(edges_path, map, &n_edges, n_nodes);
    printf("      %d arêtes chargées\n", n_edges);

    /* ── 3. CSR ── */
    printf("[3/3] Construction de la représentation CSR …\n");
    CSRGraph *g = build_csr(n_nodes, edges, n_edges);
    printf("      CSR construit : %d noeuds, %d arêtes\n",
           g->n_nodes, g->n_edges);

    /* Empreinte mémoire approx. */
    long mem_bytes = (long)(n_nodes + 1) * sizeof(int)      /* row_ptr  */
                   + (long)n_edges       * sizeof(int)       /* adj      */
                   + (long)n_edges       * sizeof(double)    /* weights  */
                   + (long)n_nodes       * sizeof(Node);     /* coords   */
    printf("      Mémoire CSR ≈ %.1f Mo\n", mem_bytes / 1e6);

    /* ── Validation ── */
    validate(g, nodes);

    /* ── Nettoyage ── */
    free(edges);
    free(nodes);
    map_free(map);
    csr_free(g);

    printf("\n=== Terminé ===\n");
    return 0;
}
