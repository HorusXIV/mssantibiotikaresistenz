import pandas as pd
import networkx as nx
import numpy as np
from pathlib import Path

out_dir = Path(".")  # oder z.B. Path("./out")

def build_system_map():
    """
    Konzeptuelle Systemkarte (gerichteter Einfluss-/Fluss-Graph) für AMR:
    - Layer: macro / micro / external
    - Kopplung: Makro -> Mikro (Selektionsdruck, unvollständige Therapie),
                Mikro -> Makro (Resistenzanteil, Therapieversagen, Transmission)
    """
    G = nx.DiGraph(name="AMR_System_Map")

    nodes = [
        # External context / boundary
        {"id":"Community", "label":"Community", "layer":"external", "category":"context"},
        {"id":"Regulator", "label":"Regulator/Guidelines", "layer":"external", "category":"context"},

        # Macro layer (hospital & patients)
        {"id":"Hospitals", "label":"Hospitals (Macro)", "layer":"macro", "category":"entity"},
        {"id":"Patients", "label":"Patients", "layer":"macro", "category":"entity"},
        {"id":"Admissions", "label":"Admissions", "layer":"macro", "category":"process"},
        {"id":"IntraContact", "label":"Intra-hospital contacts", "layer":"macro", "category":"process"},
        {"id":"Screening", "label":"Screening & testing", "layer":"macro", "category":"intervention"},
        {"id":"Isolation", "label":"Isolation / quarantine", "layer":"macro", "category":"intervention"},
        {"id":"Hygiene", "label":"Hygiene compliance", "layer":"macro", "category":"intervention"},
        {"id":"AbxUse", "label":"Antibiotic prescribing", "layer":"macro", "category":"intervention"},
        {"id":"IncompleteTx", "label":"Incomplete treatment", "layer":"macro", "category":"behavior"},
        {"id":"Transfers", "label":"Patient transfers", "layer":"macro", "category":"process"},
        {"id":"InterHospSpread", "label":"Inter-hospital spread", "layer":"macro", "category":"outcome"},
        {"id":"ColonPressure", "label":"Colonization pressure", "layer":"macro", "category":"state"},
        {"id":"TransProb", "label":"Transmission probability", "layer":"macro", "category":"state"},
        {"id":"DetectProb", "label":"Detection probability", "layer":"macro", "category":"state"},
        {"id":"TreatFailure", "label":"Treatment failure risk", "layer":"macro", "category":"outcome"},
        {"id":"InfectDuration", "label":"Infection/colonization duration", "layer":"macro", "category":"state"},

        # Micro layer (within-host evolution)
        {"id":"WithinHost", "label":"Within-host evolution", "layer":"micro", "category":"process"},
        {"id":"SelPressure", "label":"Selection pressure (antibiotics)", "layer":"micro", "category":"mechanism"},
        {"id":"MutationHGT", "label":"Mutation / HGT", "layer":"micro", "category":"mechanism"},
        {"id":"ResFrac", "label":"Resistant fraction", "layer":"micro", "category":"state"},
        {"id":"FitnessCost", "label":"Fitness cost (no antibiotics)", "layer":"micro", "category":"mechanism"},
        {"id":"TransFitness", "label":"Transmission fitness", "layer":"micro", "category":"state"},
    ]

    # edges: (src, tgt, sign, weight, relation)
    edges = [
        # Boundary / context influences
        ("Community", "Admissions", "+", 1.0, "Community inflow to hospital"),
        ("Regulator", "Screening", "+", 0.6, "Guidelines affect screening policy"),
        ("Regulator", "AbxUse", "-", 0.6, "Stewardship can reduce unnecessary use"),

        # Macro structure
        ("Hospitals", "Admissions", "+", 1.0, "Hospitals receive admissions"),
        ("Patients", "Admissions", "+", 1.0, "Patients are admitted"),
        ("Admissions", "IntraContact", "+", 0.8, "More admissions → more contacts"),
        ("IntraContact", "TransProb", "+", 0.8, "Contacts drive transmission probability"),
        ("Hygiene", "TransProb", "-", 0.9, "Hygiene reduces transmission"),
        ("Isolation", "TransProb", "-", 0.9, "Isolation reduces transmission"),
        ("Screening", "DetectProb", "+", 0.9, "Testing increases detection"),
        ("DetectProb", "Isolation", "+", 0.8, "Detection enables isolation"),

        # Inter-hospital spread
        ("Transfers", "InterHospSpread", "+", 1.0, "Transfers drive between-hospital spread"),
        ("InterHospSpread", "ColonPressure", "+", 0.7, "Imported cases increase pressure"),
        ("ColonPressure", "TransProb", "+", 0.6, "Higher pressure increases spread"),

        # Linking macro to micro (inputs)
        ("AbxUse", "SelPressure", "+", 1.0, "Antibiotics create selection"),
        ("IncompleteTx", "WithinHost", "+", 0.9, "Incomplete therapy selects resistance"),
        ("SelPressure", "WithinHost", "+", 1.0, "Selection acts within host"),
        ("MutationHGT", "WithinHost", "+", 0.7, "New variants arise"),

        # Micro dynamics (outputs)
        ("WithinHost", "ResFrac", "+", 1.0, "Evolution changes resistant fraction"),
        ("FitnessCost", "ResFrac", "-", 0.4, "Cost reduces resistance without antibiotics"),
        ("ResFrac", "TreatFailure", "+", 0.9, "Resistance increases failure risk"),
        ("ResFrac", "TransFitness", "+", 0.6, "If low cost, fitness remains high"),
        ("TransFitness", "TransProb", "+", 0.4, "Fitter strains transmit better"),

        # Macro outcomes feedback
        ("TreatFailure", "InfectDuration", "+", 0.7, "Failure prolongs duration"),
        ("InfectDuration", "ColonPressure", "+", 0.8, "Longer duration increases pressure"),
        ("TransProb", "ColonPressure", "+", 0.8, "Transmission increases pressure (feedback)"),

        # Optional: outbreak feedback to policy tightening
        ("ColonPressure", "Screening", "+", 0.4, "High pressure triggers more screening"),
        ("ColonPressure", "Isolation", "+", 0.3, "High pressure triggers more isolation"),
    ]

    for n in nodes:
        G.add_node(n["id"], label=n["label"], layer=n["layer"], category=n["category"])

    for src, tgt, sign, w, rel in edges:
        G.add_edge(src, tgt, sign=sign, weight=float(w), relation=rel)

    return G

def build_transfer_network(n_hospitals=10, density=0.22, seed=7):
    """
    Reines Spital-zu-Spital Transfernetz (gerichtet, gewichtet).
    """
    rng = np.random.default_rng(seed)
    G = nx.DiGraph(name="Hospital_Transfer_Network")

    hospital_ids = [f"H{i:02d}" for i in range(1, n_hospitals+1)]
    for hid in hospital_ids:
        G.add_node(hid, label=f"Hospital {hid}", layer="macro", category="hospital")

    for i in range(n_hospitals):
        for j in range(n_hospitals):
            if i == j:
                continue
            if rng.random() < density:
                w = float(rng.integers(1, 30))  # toy: transfers per time unit
                G.add_edge(hospital_ids[i], hospital_ids[j], weight=w, relation="patient_transfer")
    return G

def export_for_gephi(G, stem):
    # GEXF for Gephi
    gexf_path = out_dir / f"{stem}.gexf"
    nx.write_gexf(G, gexf_path)

    # Optional CSV export
    nodes = []
    for nid, attrs in G.nodes(data=True):
        row = {"Id": nid, "Label": attrs.get("label", nid)}
        for k, v in attrs.items():
            if k != "label":
                row[k] = v
        nodes.append(row)
    nodes_df = pd.DataFrame(nodes)

    edges = []
    for u, v, attrs in G.edges(data=True):
        row = {"Source": u, "Target": v, "Type": "Directed" if G.is_directed() else "Undirected"}
        for k, val in attrs.items():
            row[k] = val
        edges.append(row)
    edges_df = pd.DataFrame(edges)

    nodes_path = out_dir / f"{stem}_nodes.csv"
    edges_path = out_dir / f"{stem}_edges.csv"
    nodes_df.to_csv(nodes_path, index=False, encoding="utf-8")
    edges_df.to_csv(edges_path, index=False, encoding="utf-8")

    print("Wrote:", gexf_path, nodes_path, edges_path)

if __name__ == "__main__":
    system_map = build_system_map()
    transfer_net = build_transfer_network()

    export_for_gephi(system_map, "amr_system_map")
    export_for_gephi(transfer_net, "amr_transfer_network")