"""Join browser residue picking to the same WT alignment and colors used in exports."""
import csv
from collections import defaultdict

import maxitia_mapping as core


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def geometry(path, chain):
    # Use the polymer label ID as the picking key, even if author IDs have offsets
    # or insertion codes. The downloaded model retains all original identifiers.
    atoms = core.parse_structure_cif(path, chain)["atoms"]
    fields = "group_PDB id type_symbol label_atom_id label_comp_id label_asym_id label_seq_id auth_seq_id Cartn_x Cartn_y Cartn_z"
    lines = ["data_selected_chain", "loop_", *("_atom_site." + f for f in fields.split())]
    seen = set()
    for a in atoms:
        key = (a["label_seq_id"], a["atom_name"])
        if key in seen:
            continue
        seen.add(key)
        lines.append(" ".join(map(str, ("ATOM", len(seen), a["element"], a["atom_name"],
                                        a["comp_id"], "A", a["label_seq_id"], a["label_seq_id"], a["x"], a["y"], a["z"]))))
    return "\n".join(lines + ["#", ""])


def single_rows(folder):
    rows = []
    for r in read_rows(folder / "structure_mapped_residue_effects.csv"):
        if not int(r["n_genotype_pass"]) or r["alignment_status"] != "exact":
            continue
        kind = core.substitution_display(r)["visualization_class"]
        rows.append(dict(position=int(r["sequence_position_1based"]), variant=r["mutation_labels"].split(";")[0],
                         multi=r["effect_attribution"] == "multi_substitution_confounded",
                         status=kind, ratios=[core.to_float(r[k + "_fold_vs_WT_pass"]) if kind == "QUANTITATIVE"
                                             else None for k in ("KD","kon","koff")],
                         support=f'{r["n_active_quantitative_pass"]} PASS active · {r["n_inactive"]} inactive',
                         note=r["phenotype_consensus"]))
    return rows


def evidence_rows(wt, rows, batch=False):
    effects = []
    for r in rows:
        if r["variant_category"] not in ("single_substitution", "multi_substitution"):
            continue
        support = (f'{r["n_pass_active"]} PASS active · {r["n_inactive"]} inactive' if batch else
                   f'{r["batches_with_fold"]}/{r["batches_observed"]} observed batches with a fold')
        for i,(a,b) in enumerate(zip(wt,r["aa_sequence"]),1):
            if a != b:
                effects.append(dict(position=i, variant=r["variant"], multi=r["variant_category"]=="multi_substitution",
                                    status=r["structure_status"], support=support,
                                    note=r["phenotype"] if batch else r["interpretation"],
                                    ratios=[r[k + ("_fold" if batch else "_fold_center")] for k in ("KD","kon","koff")]))
    return effects


def scope(label, table, effects):
    grouped = defaultdict(list)
    for item in effects:
        grouped[item["position"]].append(item)
    residues = []
    for r in table:
        p = int(r["sequence_position_1based"])
        if r["alignment_status"] != "exact" or not r["structure_label_seq_id"]:
            continue
        residues.append(dict(pos=p, aa=r["sequence_aa"], model=int(r["structure_label_seq_id"]),
                             author=f'{r["structure_auth_chain_id"]} {r["structure_residue_number"]}{r["structure_insertion_code"]}',
                             color=r["visualization_color_hex"], status=r["visualization_class"],
                             selected=r.get("selected_substitution",""), energy=core.to_float(r.get("display_ddG_kcal_mol")),
                             effects=sorted(grouped[p], key=lambda e: (e["multi"], e["variant"])) ))
    return dict(label=label, residues=residues)


def demo_structure(home):
    # Illustrative four-residue peptide in the same frame as the synthetic wells.
    names = ["ALA","CYS","ASP","GLU"]
    fields = "group_PDB id label_entity_id label_asym_id auth_asym_id label_seq_id auth_seq_id label_atom_id type_symbol label_comp_id Cartn_x Cartn_y Cartn_z B_iso_or_equiv"
    lines = ["data_demo", "loop_", "_entity_poly_seq.entity_id","_entity_poly_seq.num","_entity_poly_seq.mon_id"]
    lines += [f"1 {i} {aa}" for i, aa in enumerate(names,1)]
    lines += ["#", "loop_", *("_atom_site." + f for f in fields.split())]
    serial = 0
    for i,aa in enumerate(names,1):
        sign = (-1)**i
        for atom, x,y,z in [("N",-1.16,.20,0),("CA",0,.95,.18),("C",1.22,.12,.15),
                            ("O",1.3,-1.08,.35),("CB",.08,1.6,1.55)]:
            serial += 1
            lines.append(f"ATOM {serial} 1 A A {i} {i} {atom} {atom[0]} {aa} {3.72*i+x:.3f} {sign*y:.3f} {z:.3f} 90")
    path = home / "demo_structure.cif"
    path.write_text("\n".join(lines+["#",""]))
    return path


def write_html(path, viewer):
    """Portable result with the exact browser map, data, renderer and notices."""
    import json
    from html import escape
    from pathlib import Path
    assets = Path(__file__).resolve().parent
    def asset(name):
        return (assets / name).read_text(encoding="utf-8")
    payload = json.dumps(viewer, ensure_ascii=True, allow_nan=False, separators=(",", ":")).replace("<", "\\u003c")
    notices = escape(asset("vendor/3Dmol-LICENSE.txt") + "\n" + asset("vendor/3Dmol-min.js.LICENSE.txt"))
    path.write_text(f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MAXTIA · Structure explorer</title><style>{asset("style.css")}</style></head><body>
<main><h1>MAXTIA structure map</h1><section id="modelPanel" class="panel"></section>
<details><summary>Software licenses</summary><pre>{notices}</pre></details></main>
<script>{asset("vendor/3Dmol-min.js")}</script><script>{asset("molecular-view.js")}</script><script>{asset("explorer.js")}</script>
<script type="application/json" id="study-data">{payload}</script>
<script>structureExplorer.show(JSON.parse(document.getElementById("study-data").textContent));</script>
</body></html>""", encoding="utf-8")
