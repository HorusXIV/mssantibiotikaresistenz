import pandas as pd
import networkx as nx
from pathlib import Path

BASE_DIR = Path(__file__).parent
out_dir = BASE_DIR


def build_system_map(nodes_df, edges_df) -> nx.DiGraph:
    """Build the AMR system map from the CSV files in System_Overview."""
    G = nx.DiGraph(name="AMR_System_Map")

    for row in nodes_df.to_dict("records"):
        G.add_node(
            row["Id"],
            label=row.get("Label", row["Id"]),
            layer=row.get("layer", ""),
            category=row.get("category", ""),
        )

    for row in edges_df.to_dict("records"):
        G.add_edge(
            row["Source"],
            row["Target"],
            sign=row.get("sign", "+"),
            weight=float(row.get("weight", 1.0)),
            relation=row.get("relation", ""),
        )

    return G


def export_for_gephi(G, stem):
    gexf_path = out_dir / f"{stem}.gexf"
    nx.write_gexf(G, gexf_path)

    nodes = []
    for nid, attrs in G.nodes(data=True):
        row = {"Id": nid, "Label": attrs.get("label", nid)}
        for key, value in attrs.items():
            if key != "label":
                row[key] = value
        nodes.append(row)
    nodes_df = pd.DataFrame(nodes)

    edges = []
    for source, target, attrs in G.edges(data=True):
        row = {
            "Source": source,
            "Target": target,
            "Type": "Directed" if G.is_directed() else "Undirected",
        }
        for key, value in attrs.items():
            row[key] = value
        edges.append(row)
    edges_df = pd.DataFrame(edges)

    nodes_path = out_dir / f"{stem}_nodes.csv"
    edges_path = out_dir / f"{stem}_edges.csv"
    nodes_df.to_csv(nodes_path, index=False, encoding="utf-8")
    edges_df.to_csv(edges_path, index=False, encoding="utf-8")

    print("Wrote:", gexf_path, nodes_path, edges_path)


if __name__ == "__main__":
    nodes_df = pd.read_csv(BASE_DIR / "amr_system_map_nodes.csv")
    edges_df = pd.read_csv(BASE_DIR / "amr_system_map_edges.csv")
    system_map = build_system_map(nodes_df, edges_df)
    export_for_gephi(system_map, "amr_system_map")
