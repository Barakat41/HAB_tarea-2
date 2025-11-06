#!/usr/bin/env python3
import argparse
import csv
import sys
from collections import defaultdict
import networkx as nx
import numpy as np
from scipy.stats import hypergeom

# -----------------------
# Utilidades
# -----------------------

def read_network(path, score_threshold=400):
    """
    Lee una red desde fichero:
    - 2 columnas (nodeA, nodeB) o 3 columnas (nodeA, nodeB, score)
    - Si hay una columna numérica en la 3a posición, se filtra por score_threshold.
    Devuelve un networkx.Graph

    path: ruta al fichero
    score_threshold: si hay scores, filtrar por este umbral
    """
    G = nx.Graph()
    with open(path, 'r', encoding='utf-8') as f:
        sniffer = csv.Sniffer()
        sample = f.read(2048)
        f.seek(0)
        try:
            dialect = sniffer.sniff(sample)
            delim = dialect.delimiter
        except Exception:
            # fallback
            delim = '\t' if '\t' in sample else ',' if ',' in sample else None
            if delim is None:
                delim = None  # will use split()

        reader = csv.reader(f, delimiter=delim) if delim else (line.split() for line in f)
        for row in reader:
            if not row:
                continue
            # limpiar espacios
            row = [c.strip() for c in row if c is not None]
            if len(row) < 2:
                continue
            # ignorar líneas que parezcan encabezado
            if row[0].lower() in ('protein1', 'node1', 'gene', 'protein') or row[1].lower() in ('protein2', 'node2', 'gene', 'protein'):
                continue
            node1, node2 = row[0], row[1]
            # si hay una 3a columna numérica, interpretarla como score
            if len(row) >= 3:
                try:
                    score = float(row[2])
                    if score >= score_threshold:
                        G.add_edge(node1, node2)
                except ValueError:
                    # 3a columna no es numérica -> asumimos edgelist simple
                    G.add_edge(node1, node2)
            else:
                G.add_edge(node1, node2)
    return G

def read_seeds(path):
    """
    Lee un fichero de semillas: una por línea o separadas por comas.
    Si no se recibe path, devuelve las semillas por defecto (ENO1, PGK1, HK2)

    path: ruta al fichero de semillas
    """
    default = {'ENO1', 'PGK1', 'HK2'}
    if not path:
        return default
    seeds = set()
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.replace(',', '\n').split()]
            for p in parts:
                if p:
                    seeds.add(p)
    return seeds if seeds else default

# -----------------------
# Función DIAMOnD
# -----------------------

def diamond(G, seeds, X, alpha=1, verbose=True):
    """
    Versión simple de DIAMOnD:
    G: networkx graph
    seeds: iterable de nodos iniciales (sus nombres deben existir en G)
    X: número de nodos a añadir
    alpha: integer weight for seeds (simple: s = alpha * |seeds|)
    Devuelve la lista de nodos añadidos en orden.
    """
    # validar semillas y avisar si faltan
    seeds_present = set(s for s in seeds if s in G)
    missing = set(seeds) - seeds_present
    if missing:
        print("WARNING: some seed nodes not found in network:", ", ".join(sorted(missing)))
    if not seeds_present:
        print("ERROR: none of the seed nodes are present in the network. Exiting.")
        return []

    # estructura de datos
    cluster = set(seeds_present)  # nodos del cluster (semillas + añadidos)
    added = []
    nodes = set(G.nodes())
    N = len(nodes)  # tamaño del universo para hipergeométrica

    if verbose:
        print(f"Network nodes: {N}, seeds in network: {len(seeds_present)}")

    # precalcular vecinos y grados
    neighbors = {n: set(G.neighbors(n)) for n in nodes}
    degrees = {n: G.degree(n) for n in nodes}

    # bucle principal
    while len(added) < X:
        best_node = None
        best_p = 1.0
        # tamaño efectivo de la muestra (s). Usamos alpha * |cluster|, cap a N
        s = int(alpha * len(cluster))
        if s < 1:
            s = len(cluster)
        if s > N:
            s = N

        # candidatos: todos los nodos que no están en el cluster
        candidates = nodes - cluster
        if not candidates:
            if verbose:
                print("No more candidate nodes available.")
            break

        for node in candidates:
            k = degrees.get(node, 0)
            if k == 0:
                continue
            # kb = numero vecinos en el cluster
            kb = sum(1 for nb in neighbors.get(node, ()) if nb in cluster)
            # usamos la cola superior de la hipergeométrica: P(X >= kb)
            # hipergeom.sf(kb-1, N, s, k) es exacto y estable
            try:
                pval = hypergeom.sf(kb - 1, N, s, k)
            except Exception:
                # en caso de problemas numéricos, saltar
                pval = 1.0
            # escoger el mínimo p
            if pval < best_p:
                best_p = pval
                best_node = node
        if best_node is None:
            if verbose:
                print("No node selected in this iteration (maybe isolated nodes).")
            break
        # añadir elegido
        added.append(best_node)
        cluster.add(best_node)
        # actualizar estructuras auxiliares: (vecinos y grados no cambian)
        if verbose:
            print(f"Added node: {best_node} (p={best_p:.2e}) — total added: {len(added)}/{X}")
    return added

# -----------------------
# CLI
# -----------------------
def parse_args():
    """
    Parseo de argumentos de ejecución
    """
    parser = argparse.ArgumentParser(description="Simple DIAMOnD script")
    parser.add_argument("network", help="Network file (edgelist or STRING-like TSV)")
    parser.add_argument("-s", "--seeds", help="Seeds file (one gene per line). If omitted, uses ENO1,PGK1,HK2", default=None)
    parser.add_argument("-n", "--num", help="Number of DIAMOnD nodes to add (X)", type=int, default=200)
    parser.add_argument("-a", "--alpha", help="Seed weight alpha (>=1)", type=int, default=1)
    parser.add_argument("-o", "--out", help="Output file for added nodes", default=None)
    parser.add_argument("--score-threshold", help="Score threshold (default 400)", type=float, default=400.0)
    return parser.parse_args()

def main():
    args = parse_args()
    print("INFO: reading network...")
    G = read_network(args.network, score_threshold=args.score_threshold)
    print(f"INFO: network loaded. Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    seeds = read_seeds(args.seeds)
    print("INFO: using seeds:", ", ".join(sorted(seeds)))

    print("INFO: running DIAMOnD...")
    added = diamond(G, seeds, X=args.num, alpha=max(1, args.alpha), verbose=True)

    out_file = args.out if args.out else f"output.txt"
    with open(out_file, 'w', encoding='utf-8') as fo:
        for node in added:
            fo.write(node + "\n")
    print("DONE: results saved to", out_file)

if __name__ == "__main__":
    main()
