import osmium
import csv
import math

# ─── Formule Haversine ───────────────────────────────────────────────
def haversine(lat1, lon1, lat2, lon2):
    R = 6_371_000  # rayon de la Terre en mètres
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    # CORRECTION : ** 2 (puissance) et non * 2 (multiplication)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))

# ─── Handler OSM ────────────────────────────────────────────────────
class RoadHandler(osmium.SimpleHandler):
    def __init__(self):
        super().__init__()
        self.nodes = {}   # id_osm -> (lat, lon)
        self.edges = []   # liste de (u, v)

    def node(self, n):
        self.nodes[n.id] = (n.location.lat, n.location.lon)

    def way(self, w):
        highway = w.tags.get("highway")
        if highway not in (
            "motorway", "trunk", "primary", "secondary",
            "tertiary", "residential", "unclassified",
            "motorway_link", "trunk_link", "primary_link"
        ):
            return

        nodes = list(w.nodes)
        is_oneway = w.tags.get("oneway", "no") in ("yes", "1", "true")
        is_roundabout = w.tags.get("junction") == "roundabout"
        # Les voies de type motorway sont toujours à sens unique
        is_oneway = is_oneway or is_roundabout or highway in ("motorway", "motorway_link")

        for i in range(len(nodes) - 1):
            u, v = nodes[i].ref, nodes[i + 1].ref
            self.edges.append((u, v))           # sens direct
            if not is_oneway:
                self.edges.append((v, u))        # sens inverse (double sens)

# ─── Extraction ──────────────────────────────────────────────────────
handler = RoadHandler()
handler.apply_file("champagne-ardenne-260412.osm.pbf", locations=True)
print(f"Noeuds chargés      : {len(handler.nodes)}")
print(f"Arêtes extraites    : {len(handler.edges)}")

# ─── Calcul des poids (distance en mètres) ──────────────────────────
weighted_edges = []
skipped = 0
for u, v in handler.edges:
    if u in handler.nodes and v in handler.nodes:
        lat1, lon1 = handler.nodes[u]
        lat2, lon2 = handler.nodes[v]
        dist = haversine(lat1, lon1, lat2, lon2)
        weighted_edges.append((u, v, dist))
    else:
        skipped += 1

print(f"Arêtes pondérées    : {len(weighted_edges)}")
print(f"Arêtes ignorées     : {skipped}  (noeuds manquants)")

# ─── Nettoyage : garder uniquement les noeuds utilisés ──────────────
used_node_ids = set()
for u, v, _ in weighted_edges:
    used_node_ids.add(u)
    used_node_ids.add(v)

print(f"Noeuds utilisés     : {len(used_node_ids)}")

# ─── Sauvegarde CSV ──────────────────────────────────────────────────
with open("nodes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "lat", "lon"])
    for nid in used_node_ids:
        lat, lon = handler.nodes[nid]
        w.writerow([nid, lat, lon])

with open("edges.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["u", "v", "distance_m"])
    w.writerows(weighted_edges)

print("Sauvegarde terminée : nodes.csv + edges.csv")