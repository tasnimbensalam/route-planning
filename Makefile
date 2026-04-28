# ─── Makefile pour le projet Route Planning ──────────────────────────
#
# Cibles principales :
#   make           → compile en mode release (O2)
#   make fast      → compile avec O3 + march=native (pour les bench finaux)
#   make debug     → compile avec ASan + UBSan (pour traquer les bugs)
#   make run       → compile puis lance avec les valeurs par défaut
#   make bench     → compile puis lance avec 1000 requêtes + export CSV
#   make extract   → relance l'extraction OSM (si .pbf présent)
#   make clean     → supprime le binaire et les CSV générés

CC      = gcc
CFLAGS  = -O2 -Wall -Wextra
LDFLAGS = -lm
SRC     = load_graph.c
BIN     = load_graph

# ─── Cible par défaut : compilation release ────────────────────────────
all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

# ─── Compilation pour les mesures finales du rapport ──────────────────
# -O3 active toutes les optimisations agressives
# -march=native utilise le jeu d'instructions natif du CPU (SSE/AVX)
# -DNDEBUG désactive les éventuels asserts
fast: $(SRC)
	$(CC) -O3 -march=native -DNDEBUG -o $(BIN) $< $(LDFLAGS)

# ─── Compilation debug avec sanitizers ─────────────────────────────────
# Détecte les fuites mémoire, les accès hors-borne, les comportements UB
# Ralentit l'exécution mais indispensable quand on touche au code
debug: $(SRC)
	$(CC) -O1 -g -fsanitize=address,undefined -Wall -Wextra \
	      -o $(BIN)_dbg $< $(LDFLAGS)

# ─── Lancements pratiques ──────────────────────────────────────────────
run: $(BIN)
	./$(BIN) nodes.csv edges.csv

bench: $(BIN)
	./$(BIN) --queries 1000 --csv results.csv nodes.csv edges.csv

# Bench quand CH est trop lent (en cours de dev sur gros graphe)
bench-fast: $(BIN)
	./$(BIN) --no-ch --queries 1000 --csv results.csv nodes.csv edges.csv

# ─── Extraction OSM (nécessite osmium en Python) ──────────────────────
extract:
	python3 extract.py

# ─── Génération des figures pour le rapport ───────────────────────────
figures: results.csv
	python3 plot_results.py results.csv

# ─── Nettoyage ─────────────────────────────────────────────────────────
clean:
	rm -f $(BIN) $(BIN)_dbg results.csv figures/*.png

distclean: clean
	rm -f nodes.csv edges.csv

.PHONY: all fast debug run bench bench-fast extract figures clean distclean
