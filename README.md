

# Route Planning — Moteur de plus court chemin sur graphe routier réel

> Projet de Master 1 — Structure de données avancées  
> Université Sorbonne Paris Nord · Année 2025–2026  
> Encadrant : M. Olivier Bodini  
> Étudiants : Hassan Charaf · Tesnime Bensalem · Melissa Hachemi

---

## Table des matières

- [Présentation](#présentation)
- [Algorithmes implémentés](#algorithmes-implémentés)
- [Structure du projet](#structure-du-projet)
- [Prérequis](#prérequis)
- [Données d'entrée](#données-dentrée)
- [Compilation](#compilation)
- [Utilisation](#utilisation)
- [Résultats de référence](#résultats-de-référence)
- [Détails techniques](#détails-techniques)

---

## Présentation

Ce projet implémente un **moteur de route planning** capable de répondre efficacement à des requêtes de plus court chemin `(s, t)` sur un graphe routier réel extrait d'OpenStreetMap.

Le graphe utilisé couvre la région **Champagne-Ardenne** et compte :
- **905 283 nœuds** (intersections)
- **1 753 308 arêtes** (segments de route orientés)
- **67,8 Mo** d'empreinte mémoire totale (structure CSR + graphe inverse + coordonnées)

Quatre algorithmes sont étudiés et comparés sur 1 000 requêtes aléatoires reproductibles (`--seed 42`).

---

## Algorithmes implémentés

| Algorithme | Prétraitement | Mémoire auxiliaire | Speedup moyen | Nœuds explorés (moy.) |
|---|---|---|---|---|
| **Dijkstra** (référence) | — | — | ×1,0 | 457 459 |
| **A\*** (heuristique haversine) | — | — | ×1,6 | 119 395 |
| **ALT** (A\* + Landmarks, K=16) | ~5,8 s | 221 Mo | ×2,9 | 39 259 |
| **CH** (Contraction Hierarchies) | ~65 s | 119 Mo | ×14,3 | 1 344 |

### Dijkstra
Algorithme de référence avec arrêt anticipé sur la cible. Utilise un tas binaire min avec mises à jour paresseuses (lazy deletion). Complexité : O((|V| + |E|) log |V|).

### A\*
Extension de Dijkstra avec heuristique **haversine** (distance à vol d'oiseau vers la cible). L'heuristique est admissible par construction — la distance réelle est toujours supérieure à la ligne droite sur le globe. Aucun prétraitement requis.

### ALT — A\* avec Landmarks
Heuristique plus serrée basée sur l'**inégalité triangulaire** et un ensemble de points repères (landmarks) présélectionnés. Le prétraitement calcule les distances exactes depuis/vers chaque landmark vers tous les nœuds du graphe. La sélection utilise la stratégie **farthest-first** restreinte à la composante connexe principale (99,1 % des nœuds sur ce graphe).

Le nombre de landmarks `K` est configurable via `--landmarks`. L'étude montre que **K=8 offre le meilleur compromis** sur ce graphe : au-delà, le coût d'évaluation de l'heuristique (linéaire en K) et les cache misses sur les tables de distance (> 200 Mo à K=16) dégradent les performances malgré un nombre de nœuds explorés plus faible.

### CH — Contraction Hierarchies
Approche de prétraitement intensif : les nœuds sont contractés un par un selon leur **edge difference**, en ajoutant des raccourcis conservateurs pour préserver les plus courts chemins. La requête est ensuite un Dijkstra bidirectionnel ne suivant que les arêtes **montantes** dans la hiérarchie. CH est quasiment insensible à la longueur du trajet (×1,7 de SHORT à LONG, contre ×20 pour Dijkstra).

---

## Structure du projet

```
.
├── load_graph.c        # Source unique — tous les algorithmes, benchmark et export CSV
├── extract.py          # Script Python d'extraction OSM → nodes.csv + edges.csv
├── nodes.csv           # Nœuds du graphe (id, latitude, longitude)
├── edges.csv           # Arêtes du graphe (u, v, distance_m)
└── README.md
```

Tout le code C est contenu dans un **fichier source unique** structuré en sections numérotées pour faciliter la lecture et les extensions. Les paramètres globaux (seuils de catégories, nombre de landmarks par défaut, limite de la recherche de témoins CH) sont centralisés en tête de fichier.

---

## Prérequis

**Compilateur**
- GCC ≥ 9 (ou Clang ≥ 10)
- Flag recommandé : `-O2`

**Extraction des données (optionnel)**
- Python 3.8+
- Bibliothèque `osmium` pour lire les fichiers `.osm.pbf`

**Aucune dépendance externe** n'est requise pour la compilation du moteur C (bibliothèque standard uniquement : `math.h`, `time.h`).

---

## Données d'entrée

Le moteur attend deux fichiers CSV dans le répertoire courant (ou passés en argument) :

**`nodes.csv`** — une ligne par nœud, sans espace :
```
id,latitude,longitude
123456789,49.2441,4.0337
...
```

**`edges.csv`** — une ligne par arête orientée, poids en mètres :
```
u,v,distance
123456789,987654321,142.7
...
```

### Génération depuis OpenStreetMap

```bash
# Télécharger une région sur https://download.geofabrik.de/
python extract.py champagne-ardenne-latest.osm.pbf
# Produit nodes.csv et edges.csv
```

Le script filtre les voies carrossables (`motorway`, `trunk`, `primary`, `secondary`, `tertiary`, `residential`, `unclassified` et leurs liens), respecte les sens uniques et gère les ronds-points.

---

## Compilation

```bash
# Compilation standard (recommandée)
gcc -O2 -Wall -Wextra -o load_graph load_graph.c -lm
---

## Utilisation

```
./load_graph [OPTIONS]
```

### Options

| Option | Défaut | Description |
|---|---|---|
| `--nodes <fichier>` | `nodes.csv` | Fichier des nœuds |
| `--edges <fichier>` | `edges.csv` | Fichier des arêtes |
| `--queries <N>` | `1000` | Nombre de requêtes du benchmark |
| `--seed <S>` | `42` | Graine aléatoire (reproductibilité) |
| `--landmarks <K>` | `16` | Nombre de landmarks ALT |
| `--no-ch` | — | Désactive CH (prétraitement long) |
| `--csv <fichier>` | — | Export détaillé des résultats par requête |

### Exemples

```bash
# Benchmark complet avec paramètres par défaut
./load_graph

# Benchmark rapide sans CH (pas de prétraitement de 65 s)
./load_graph --no-ch --queries 500

# Étude de l'effet du nombre de landmarks
./load_graph --no-ch --landmarks 8
./load_graph --no-ch --landmarks 16
./load_graph --no-ch --landmarks 32

# Export CSV pour analyse externe (Python, R…)
./load_graph --csv resultats.csv
```
### commande pour lance 
compiler gcc -O2 -o load_graph load_graph.c -lm        
lancer  ./load_graph nodes.csv edges.csv  

### Sortie console

```
=== Benchmark sur 1000 requêtes ===
    [50/1000] requêtes traitées (5%)
    ...
  [OK] 984 requêtes vérifiées (distances + chemins)

=== Récapitulatif global (984 requêtes valides sur 1000) ===
  Algo     |  moy (ms) |  p50 (ms) |  p95 (ms) |   req/s  | nœuds expl. moy | xDijk
  ---------+-----------+-----------+-----------+----------+-----------------+-------
  Dijkstra |    55.950 |    55.120 |   104.560 |       18 |          457459 | x1.0
  A*       |    34.200 |    24.670 |    97.450 |       29 |          119395 | x1.6
  ALT      |    19.620 |    14.670 |    57.500 |       51 |           39259 | x2.9
  CH       |     3.915 |     3.810 |     5.420 |      255 |            1344 | x14.3

=== Stats par catégorie de distance ===
  Seuils : SHORT < 20km , MEDIUM 20-80km , LONG > 80km
  ...

=== Prétraitement et mémoire ===
  Graphe CSR (G + G' + nodes) : 67.8 Mo
  ALT : 5775.0 ms , 221.0 Mo (16 landmarks)
  CH  : 65.00 s , 119.2 Mo (1812084 raccourcis ; ratio ×2.03)
```

---

## Résultats de référence

Résultats obtenus sur **Champagne-Ardenne** (905 283 nœuds, 1 753 308 arêtes), 1 000 requêtes uniformes, graine 42, machine AMD / GCC -O2.

### Performances globales (984 requêtes valides)

| Algorithme | Moy. (ms) | p50 (ms) | p95 (ms) | Débit (req/s) |
|---|---|---|---|---|
| Dijkstra | 55,95 | 55,12 | 104,56 | 18 |
| A\* | 34,20 | 24,67 | 97,45 | 29 |
| ALT (K=16) | 19,62 | 14,67 | 57,50 | 51 |
| **CH** | **3,92** | **3,81** | **5,42** | **255** |

### Performances par catégorie de distance

| Algorithme | SHORT < 20 km (24 req.) | MEDIUM 20–80 km (252 req.) | LONG > 80 km (708 req.) |
|---|---|---|---|
| Dijkstra | 3,58 ms | 21,96 ms | 69,83 ms |
| A\* | 2,45 ms | 7,68 ms | 44,72 ms |
| ALT (K=16) | 2,44 ms | 7,04 ms | 24,68 ms |
| CH | 2,53 ms | 3,20 ms | 4,22 ms |

> **Observation clé :** CH est quasi-insensible à la distance (facteur ×1,7 de SHORT à LONG). Sur les requêtes courtes, A\* et ALT-8 sont compétitifs — le surcoût du bidirectionnel ne se rentabilise pas en dessous de ~20 km.

### Effet du nombre de landmarks K (ALT)

| K | Prétr. | Mémoire | SHORT | MEDIUM | LONG |
|---|---|---|---|---|---|
| 4 | 1,4 s | 55,3 Mo | 1,57 ms | 4,46 ms | 16,62 ms |
| **8** | **2,4 s** | **110,5 Mo** | **1,80 ms** | **4,63 ms** | **15,65 ms** |
| 16 | 5,8 s | 221,0 Mo | 2,44 ms | 7,04 ms | 24,68 ms |
| 32 | 10,7 s | 442,0 Mo | 3,70 ms | 12,06 ms | 43,06 ms |

Au-delà de K=8, le temps d'évaluation de l'heuristique (O(K) par relaxation) et les cache misses sur les tables de distance inversent le gain algorithmique.

---

## Détails techniques

### Représentation mémoire — CSR (Compressed Sparse Row)

Le graphe est stocké dans trois tableaux contigus :
- `row_ptr[v]` — indice de début des voisins de `v` dans `adj`
- `adj[k]` — nœud cible de la k-ième arête
- `weights[k]` — poids (distance haversine en mètres) de la k-ième arête

Un **CSR inverse** (graphe transposé) est construit par une procédure en trois passes pour les algorithmes nécessitant un parcours en arrière (ALT, CH).

### Tas binaire avec lazy deletion

Plutôt qu'un tas indexé avec `decrease-key`, le projet utilise des **mises à jour paresseuses** : les entrées obsolètes sont ignorées à l'extraction (test `f_score > dist[u]`). Plus simple à implémenter et très performant en pratique.

### Validation

Chaque résultat est vérifié sur deux critères :
1. **Distance** : correspondance à 1 mm près avec Dijkstra (référence)
2. **Chemin** : la séquence de nœuds constitue un chemin réel du graphe, avec une somme de poids égale à la distance annoncée (vérification critique pour CH avec dépaquetage récursif des raccourcis)

Aucune erreur détectée sur les 984 requêtes valides. Validé avec `-fsanitize=address,undefined` (zéro fuite mémoire, zéro comportement indéfini).

### Paramètre CH — Recherche de témoins

La recherche de témoins (witness search) lors de la contraction est bornée à **500 extractions** (`WITNESS_NODE_LIMIT`). Cette borne conservatrice peut générer quelques raccourcis superflus mais garantit la correction et maintient un temps de prétraitement raisonnable.

---

*Rapport complet disponible dans `rapport_SDA.pdf`.*
