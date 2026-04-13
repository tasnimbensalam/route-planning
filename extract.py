import osmium
import csv

class RoadHandler(osmium.SimpleHandler):
    def __init__(self):
        super().__init__()
        self.nodes = {}   # id -> (lat, lon)
        self.edges = []   # (u, v, weight_m)

    def node(self, n):
        self.nodes[n.id] = (n.location.lat, n.location.lon)

    def way(self, w):
        # Filtrer uniquement les routes
        highway = w.tags.get("highway")
        if highway not in (
            "motorway", "trunk", "primary", "secondary",
            "tertiary", "residential", "unclassified",
            "motorway_link", "trunk_link", "primary_link"
        ):
            return

        # Créer les arêtes entre noeuds consécutifs
        nodes = list(w.nodes)
        for i in range(len(nodes) - 1):
            u, v = nodes[i].ref, nodes[i+1].ref
            self.edges.append((u, v))

handler = RoadHandler()
handler.apply_file("champagne-ardenne-260412.osm.pbf", locations=True)

print(f"Noeuds chargés : {len(handler.nodes)}")
print(f"Arêtes extraites : {len(handler.edges)}")
import math

def haversine(lat1, lon1, lat2, lon2):
    R = 6_371_000  # rayon terre en mètres
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
    return 2 * R * math.asin(math.sqrt(a))

weighted_edges = []
for u, v in handler.edges:
    if u in handler.nodes and v in handler.nodes:
        lat1, lon1 = handler.nodes[u]
        lat2, lon2 = handler.nodes[v]
        dist = haversine(lat1, lon1, lat2, lon2)
        weighted_edges.append((u, v, dist))



import csv

# Sauvegarder les noeuds
with open("nodes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "lat", "lon"])
    for nid, (lat, lon) in handler.nodes.items():
        w.writerow([nid, lat, lon])

# Sauvegarder les arêtes
with open("edges.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["u", "v", "distance_m"])
    w.writerows(weighted_edges)

print(f"Sauvegardé : {len(weighted_edges)} arêtes pondérées")