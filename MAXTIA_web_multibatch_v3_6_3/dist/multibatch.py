"""Analyze independent LSA/NGS runs against a shared protein reference."""
import csv
import shutil
import json
import math
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import bridge
import structure_explorer as explorer
import maxitia_mapping as core

VERSION = "3.6.3"
METRICS = {"KD": "KD_M", "kon": "kon_M_inv_s", "koff": "koff_s_inv"}


def csv_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def demo_batches(home, include_cluster=False):
    batches = []
    for i, start in enumerate([1, 5]):
        folder = home / f"sample-{i+1}"
        folder.mkdir()
        base = .001 * (i + 1)
        wells = [(1, 1, "C", "WT", "ACDE", "None", base),
                 (2, 1, "C", "D3A", "ACAE", "None", base * (2 if i else .5)),
                 (3, 1, "C", "E4A", "ACDA", "Inactive" if i else "None", None if i else base),
                 (4, 1, "C", "A1G" if i else "C2G", "GCDE" if i else "AGDE", "None", base*.5),
                 (1, 2, "C", "C2A", "AADE", "None", base*.6)]
        if include_cluster:
            wells.extend([(2, 2, "C", "A1V", "VCDE", "None", base*3),
                          (3, 2, "C", "D3N", "ACNE", "None", base*2.5)])
        files = bridge.synthetic_inputs(folder, wells, start)
        batches.append(dict(id=f"Batch {i+1}", condition=["Buffer A", "Buffer B"][i],
                            block=f"P{start:02d}-P{start+3:02d}", paths=files))
    return batches


def prepare(batches, supplied_wt, shared_chain):
    references = {bridge.wt_input(supplied_wt)} if supplied_wt else set()
    ids, lsa_hashes, partners, chains = set(), set(), set(), set()
    for batch in batches:
        label = batch["id"].strip()
        if not label or label.casefold() in ids:
            raise ValueError("Enter a distinct, nonempty name for every batch.")
        ids.add(label.casefold())
        batch["id"] = label
        files = batch["paths"]
        digest = core.sha256(Path(files["lsa"]))
        if digest in lsa_hashes:
            raise ValueError(f"{label} reuses an identical LSA workbook. Select one independent export per batch.")
        lsa_hashes.add(digest)
        policy = json.loads(Path(files["policy"]).read_text(encoding="utf-8-sig")) if files.get("policy") else {}
        if not isinstance(policy, dict):
            raise ValueError(f"{label} requires a QC policy JSON object.")
        meta = policy.get("experiment_metadata", {})
        if not isinstance(meta, dict):
            raise ValueError(f"{label} has invalid experiment metadata in the QC policy.")
        batch["metadata"] = dict(meta)
        if "temperature_kelvin" in policy:
            batch["metadata"]["temperature_kelvin"] = policy["temperature_kelvin"]
        chains.add(str(shared_chain or policy.get("structure_chain_id", "A")))
        if meta.get("analyte_identity"):
            partners.add(str(meta["analyte_identity"]).strip().casefold())
        config = core.deep_merge(core.DEFAULT_CONFIG, policy)
        if "M" in config["accepted_ngs_call_codes"]:
            raise ValueError(f"{label}: the M call is ineligible for sequence attribution.")
        if policy.get("wt_sequence"):
            references.add(bridge.wt_input(policy["wt_sequence"]))
        if batch.get("block"):
            config.update(ngs_plate_block=batch["block"], ngs_target_plate_start=1)
        try:
            ngs, _ = core.parse_ngs(Path(files["ngs"]), config)
        except ValueError as exc:
            raise ValueError(f"{label}: {exc}") from exc
        # Only unambiguous, trusted, WT-labeled wells can establish study identity.
        wells = Counter(r["well"] for r in ngs)
        candidates = {r["aa_sequence"] for r in ngs if r["call_code"] == "C"
                      and wells[r["well"]] == 1 and r["mutation_labels"] == ["WT"]
                      and core.AA_RE.fullmatch(r["aa_sequence"])}
        references.update(candidates)
    if len(partners) > 1:
        raise ValueError("The policies identify different binding partners. Use one partner per study.")
    if len(references) != 1:
        raise ValueError("The batches have conflicting WT sequences. Check the protein, barcode blocks, and WT calls."
                         if references else "Enter a shared WT sequence. No trusted WT sequence could be established.")
    if len(chains) != 1:
        raise ValueError("The batch policies select different structure chains. Select the shared chain explicitly.")
    return references.pop(), chains.pop()


def positive(value):
    number = core.to_float(value)
    return number if number is not None and math.isfinite(number) and number > 0 else None


def log_median(values):
    return 2 ** statistics.median(math.log2(x) for x in values) if values else None


def batch_roi_columns(rows):
    return {
        key + "_by_batch": json.dumps({row["batch_id"]: row[key] for row in rows}, ensure_ascii=False)
        for key in ("observed_rois", "median_rois")
    }


def describe(rows, available, across_conditions=False):
    qualified = [r for r in rows if r["KD_fold"] is not None]
    kd = [r["KD_fold"] for r in qualified]
    n_inactive = sum(r["n_inactive"] for r in rows)
    active_seen = any(r["active_evidence"] for r in rows)
    mixed = active_seen and n_inactive > 0
    opposite = any(x < .8 for x in kd) and any(x > 1.25 for x in kd)
    caution = []
    if mixed:
        caution.append("Active/Inactive discordance")
    if opposite:
        caution.append("Opposing KD effects")
    if not qualified:
        caution.append("No eligible WT-normalized estimate")
    elif len(qualified) == 1:
        caution.append("Single quantitative batch")
    if any(not r["usable_wt"] for r in rows):
        caution.append("Batch without eligible WT")
    if any(r["variant_qc_status"] != "PASS" for r in rows):
        caution.append("Includes QC-limited observations")
    output = dict(batches_observed=len(rows), batches_total=available, batches_with_fold=len(qualified),
                  conditions_with_fold=len({r["condition_id"] for r in qualified}),
                  n_inactive_wells=n_inactive, n_pass_active_wells=sum(r["n_pass_active"] for r in rows),
                  active_inactive_discordance=mixed, opposing_KD_directions=opposite,
                  interpretation="; ".join(caution) or "Repeated quantitative evidence")
    for metric in METRICS:
        by_condition = defaultdict(list)
        for r in qualified:
            if r[f"{metric}_fold"] is not None:
                by_condition[r["condition_id"]].append(r[f"{metric}_fold"])
        vals = [x for ys in by_condition.values() for x in ys]
        centers = [log_median(xs) for xs in by_condition.values()] if across_conditions else vals
        output[f"{metric}_fold_center"] = None if mixed or opposite else log_median(centers)
        output[f"{metric}_fold_min"] = min(vals, default=None)
        output[f"{metric}_fold_max"] = max(vals, default=None)
    by_condition_ddg = defaultdict(list)
    for r in qualified:
        if r.get("ddG_kcal_mol") is not None:
            by_condition_ddg[r["condition_id"]].append(r["ddG_kcal_mol"])
    energies = [v for vals in by_condition_ddg.values() for v in vals]
    if across_conditions:
        energies = [statistics.median(vals) for vals in by_condition_ddg.values()]
    ddg = statistics.median(energies) if energies and not (mixed or opposite) else None
    output["ddG_kcal_mol"] = ddg
    statuses = {r.get("structure_status") for r in rows}
    if ddg is not None:
        kind = "QUANTITATIVE"
    elif mixed or opposite or "REVIEW" in statuses:
        kind = "REVIEW"
    elif statuses == {"INACTIVE_ONLY"}:
        kind = "INACTIVE_ONLY"
    else:
        kind = "NO_PASS_DATA"
    output["structure_status"] = kind
    return output


def analyze_study(request, home):
    batches = demo_batches(home, include_cluster=True) if request.get("demo") else request["batches"]
    if len(batches) < 2:
        raise ValueError("Add at least two independent batches.")
    wt, chain = prepare(batches, request.get("wt", ""), request.get("chain"))
    partner = (request.get("partner") or "").strip()
    if request.get("demo"):
        partner = "Example analyte"
        request = {**request, "shared": {"cif": explorer.demo_structure(home)}}
    if not partner:
        raise ValueError("Enter the common binding partner for this study.")
    output = home / "study"
    output.mkdir()
    long_rows, audits, well_rows, structures = [], [], [], []
    batch_rankings = []
    view_jobs = []
    condition_ids = {}
    counts_by_condition = Counter()
    for index, batch in enumerate(batches, 1):
        label = batch["id"]
        condition = batch.get("condition", "").strip()
        # Unspecified conditions cannot silently count as matched replicates.
        signature = (condition.casefold() or f"\0{index}", json.dumps(batch["metadata"], sort_keys=True))
        cid = condition_ids.setdefault(signature, f"C{len(condition_ids)+1}")
        counts_by_condition[cid] += 1
        batch_dir = output / f"batch_{index:02d}"
        batch_dir.mkdir()
        paths = {**request.get("shared", {}), **batch["paths"]}
        try:
            bridge.analyze(dict(paths=paths, demo=False, block=batch.get("block"), wt=wt,
                                chain=chain, study_batch=label), batch_dir, pending_views=view_jobs)
        except (ValueError, KeyError, OSError, zipfile.BadZipFile) as exc:
            raise ValueError(f"{label}: {exc}") from exc
        folder = batch_dir / "results"
        q = json.loads((folder / "qc_summary.json").read_text())
        variants = csv_rows(folder / "variant_summary.csv")
        records = csv_rows(folder / "integrated_records.csv")
        sets = sorted({r["analyte_set_name"].strip() for r in records if r["analyte_set_name"].strip()})
        if len(sets) > 1:
            raise ValueError(f"{label} contains several analyte sets. Supply a single-set export for this batch.")
        wt_control = q["wt_internal_control"]
        wt_var = next((v for v in variants if v["aa_sequence"] == wt), {})
        usable = wt_var.get("variant_qc_status") == "PASS" and wt_control["n_pass_WT"] > 0
        batch_rankings.append(core.ranked_variant_tables(variants, heading_level=4) if usable else
                              "Rankings are unavailable because this batch has no eligible WT baseline.")
        base = dict(batch_id=label, condition_id=cid, condition=condition or "Unspecified")
        baseline = {m:positive(wt_control["parameter_summary"][col]["median"]) for m,col in METRICS.items()}
        map_info = q["ngs_plate_mapping"]
        audit = dict(**base, partner=partner, folder=batch_dir.name, analyte_set_names=sets,
                     conditions_from_policy=batch["metadata"], source_files=q["source_files"],
                     ngs_plate_mapping=map_info, join=q["key_join"], qc=q["qc"]["analysis_eligibility_counts"],
                     WT_status=wt_var.get("variant_qc_status", "Unavailable"), usable_wt=bool(usable),
                     n_pass_WT=wt_control["n_pass_WT"], wt_medians=baseline, parser_warnings=q["parser_warnings"])
        audits.append(audit)
        for r in records:
            well_rows.append({**base, **r})
        for v in variants:
            if int(v["n_genotype_pass"]) == 0 or v["aa_sequence"] == wt:
                continue
            row = {**base, "aa_sequence":v["aa_sequence"], "variant":v["mutation_labels"],
                   "variant_category":v["variant_category"], "variant_qc_status":v["variant_qc_status"],
                   "observed_rois":v["observed_rois"], "median_rois":v["median_rois"],
                   "median_Rmax_RU":core.to_float(v["median_Rmax_RU_pass"]),
                   "n_pass_active":int(v["n_active_quantitative_pass"]), "n_inactive":int(v["n_inactive"]),
                   "n_genotype_pass":int(v["n_genotype_pass"]), "active_evidence":(
                        int(v["n_active_quantified"]) + int(v["n_active_qc_limited"]) > 0),
                   "usable_wt":bool(usable), "phenotype":v["phenotype_consensus"]}
            for metric, col in METRICS.items():
                val = positive(v[f"median_{col}_pass"])
                row[f"median_{col}"] = val
                row[f"WT_median_{col}"] = baseline[metric]
                row[f"{metric}_fold"] = (val / baseline[metric] if val and baseline[metric]
                                         and usable and v["variant_qc_status"] == "PASS" else None)
            row["ddG_kcal_mol"] = core.to_float(v["delta_delta_G_kcal_mol_pass"]) if row["KD_fold"] is not None else None
            display = core.substitution_display(v)
            row["structure_status"] = display["visualization_class"]
            if row["structure_status"] == "QUANTITATIVE" and row["ddG_kcal_mol"] is None:
                row["structure_status"] = "NO_PASS_DATA"
            long_rows.append(row)
        mapping = folder / "sequence_structure_mapping.csv"
        if mapping.exists() and not structures:
            structures = csv_rows(mapping)
        # The enclosing study archive carries the original files directly.
        (batch_dir / "results.zip").unlink()
    groups = defaultdict(list)
    for r in long_rows:
        groups[r["aa_sequence"]].append(r)
    overview, by_condition = [], []
    for seq, rows in sorted(groups.items(), key=lambda item: item[1][0]["variant"]):
        labels = ";".join(f"{aa}{i}{bb}" for i,(aa,bb) in enumerate(zip(wt,seq),1) if aa != bb)
        ident = dict(variant=labels or rows[0]["variant"], aa_sequence=seq, variant_category=rows[0]["variant_category"])
        overview.append({**ident, **batch_roi_columns(rows), **describe(rows, len(batches), across_conditions=True)})
        strata = defaultdict(list)
        for r in rows:
            strata[r["condition_id"]].append(r)
        for cid, subset in strata.items():
            by_condition.append({**ident, "condition_id":cid, "condition":subset[0]["condition"],
                                 **batch_roi_columns(subset), **describe(subset, counts_by_condition[cid])})
    core.write_csv(output/"variant_by_batch.csv", long_rows)
    core.write_csv(output/"integrated_wells.csv", well_rows)
    core.write_csv(output/"variant_by_condition.csv", by_condition)
    core.write_csv(output/"variant_overview.csv", overview)
    ranked = study_rankings(overview)
    core.write_csv(output / "study_ranked_variants.csv", ranked)
    structure_views, scopes, residue_tables = [], [], {}
    clusters = None
    viewer = None
    if structures:
        coordinates = {int(r["sequence_position_1based"]): r for r in structures}
        condition_labels = {a["condition_id"]: a["condition"] for a in audits}
        groups_to_map = {"study": overview}
        for cid in counts_by_condition:
            groups_to_map[cid] = [r for r in by_condition if r["condition_id"] == cid]
        mapped_conditions, positions_by_scope = [], {}
        for key, rows in groups_to_map.items():
            direct, confounded = [], []
            for r in rows:
                for pos, (aa, mutant) in enumerate(zip(wt, r["aa_sequence"]), 1):
                    point = coordinates.get(pos)
                    if aa == mutant or point is None or point["alignment_status"] != "exact":
                        continue
                    if r["variant_category"] == "single_substitution":
                        direct.append({**r, **point})
                    elif r["variant_category"] == "multi_substitution":
                        confounded.append({**point, "wt_aa": aa, "n_genotype_pass": 1,
                                           "effect_attribution": "multi_substitution_confounded"})
            if key == "study":
                core.write_csv(output / "study_structure_mapping.csv", direct)
                members, clusters = study_clusters(direct, structures)
                core.write_csv(output / "study_cluster_members.csv", members)
                core.write_csv(output / "study_clusters.csv", clusters)
            else:
                mapped_conditions.extend(direct)
            positions = structure_positions(direct)
            core.include_confounded_positions(positions, confounded)
            positions_by_scope[key] = positions
        core.write_csv(output / "condition_structure_mapping.csv", mapped_conditions)
        if view_jobs:
            model_path = Path(request["shared"]["cif"])
            shutil.copy2(model_path, output / "structure_model.cif")
            all_positions = [r for job in view_jobs for r in job["positions"]]
            all_positions.extend(r for group in positions_by_scope.values() for r in group)
            scale = core.structure_scale(all_positions, [r for job in view_jobs for r in job["variants"]])
            (output/"structure_color_scale.json").write_text(json.dumps(scale, indent=2))
            for key, positions in positions_by_scope.items():
                core.apply_structure_scale(positions, scale)
                table = core.all_residue_colors(structures, positions)
                residue_tables[key] = table
                core.write_csv(output/f"residue_colors_{key}.csv", table)
                svg = core.render_dependency_free_structure_svg(structures, positions, core.DEFAULT_CONFIG, scale)
                (output/f"structure_{key}.svg").write_text(svg)
                (output/f"structure_{key}.cxc").write_text(core.render_chimerax_overview_script(positions, core.DEFAULT_CONFIG))
                if key == "study":
                    (output/"structure_study.pml").write_text(core.render_pymol_overview_script(positions, core.DEFAULT_CONFIG))
                label = "All batches" if key == "study" else key+" / "+condition_labels[key]
                scopes.append(explorer.scope(label, table, explorer.evidence_rows(wt, groups_to_map[key])))
                if key != "study":
                    structure_views.append(dict(label=label, svg=svg))
            by_folder = {a["folder"]: a["batch_id"] for a in audits}
            for job in view_jobs:
                folder = job["out_dir"]
                q = json.loads((folder/"qc_summary.json").read_text())
                q["structure"]["visualization"].update(core.write_structure_views(job, scale))
                (folder/"qc_summary.json").write_text(json.dumps(q, ensure_ascii=False, indent=2))
                core.write_run_manifest(folder, q, job["config"])
                label = by_folder[folder.parent.name]
                batch_effects = explorer.evidence_rows(wt, [r for r in long_rows if r["batch_id"] == label], batch=True)
                scopes.append(explorer.scope("Batch / "+label, core.all_residue_colors(job["mapping"], job["positions"]), batch_effects))
            viewer = dict(geometry=explorer.geometry(model_path, chain), scopes=scopes)
            explorer.write_html(output / "structure_study.html", viewer)
    manifest = dict(version=VERSION, single_batch_engine=core.MODULE_VERSION, binding_partner=partner,
                    wt_sequence=wt, batches=audits)
    (output/"study_manifest.json").write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
    report = render_report(partner, wt, audits, overview, by_condition, residue_tables, batch_rankings, ranked, clusters)
    (output/"analysis_report.md").write_text(report)
    with zipfile.ZipFile(home/"results.zip","w",zipfile.ZIP_DEFLATED) as archive:
        for f in sorted(output.rglob("*")):
            if f.is_file():
                archive.write(f, f.relative_to(output))
    return json.dumps(dict(batches=len(batches), wells=len(well_rows), variants=len(overview),
                           repeated=sum(r["batches_with_fold"]>=2 for r in overview),
                           overview=overview, by_condition=by_condition, batch_rows=long_rows, structure_views=structure_views, viewer=viewer, clusters=clusters, report=report))


def structure_positions(entries):
    variants = []
    for row in entries:
        fields = core.structure_fields(row["structure_status"], row["ddG_kcal_mol"])
        variants.append({**row, **fields, "wt_aa": row["sequence_aa"], "mutation_label": row["variant"],
                         "variant_qc_status": "PASS" if row["structure_status"] == "QUANTITATIVE" else "REVIEW",
                         "phenotype_consensus": "mixed_active_inactive" if row["active_inactive_discordance"] else "",
                         "KD_fold_vs_WT_pass": row["KD_fold_center"]})
    return core.build_structure_visualization_positions(variants, core.DEFAULT_CONFIG)


def study_rankings(overview):
    eligible = [r for r in overview if r["variant_category"] == "single_substitution"
                and r["KD_fold_center"] is not None]
    result = []
    for label, rows in (("Improvement", sorted((r for r in eligible if r["KD_fold_center"] < .8),
                                               key=lambda r: (r["KD_fold_center"], r["variant"]))),
                        ("Loss", sorted((r for r in eligible if r["KD_fold_center"] > 1.25),
                                        key=lambda r: (-r["KD_fold_center"], r["variant"])))):
        result.extend(dict(direction=label, rank=i, **row) for i, row in enumerate(rows[:5], 1))
    return result


def study_clusters(entries, mapping):
    grouped = defaultdict(list)
    for r in entries:
        grouped[int(r["sequence_position_1based"])].append(r)
    summary = []
    for pos, variants in sorted(grouped.items()):
        qualified = [r for r in variants if r["KD_fold_center"] is not None]
        worst = max(qualified, key=lambda r: r["KD_fold_center"], default=None)
        inactive_only = all(r["structure_status"] == "INACTIVE_ONLY" for r in variants)
        first = variants[0]
        summary.append(dict(sequence_position_1based=pos, wt_aa=first["sequence_aa"],
                            worst_KD_fold_vs_WT=worst["KD_fold_center"] if worst else None,
                            worst_KD_mutation=worst["variant"] if worst else None,
                            variants_with_pass_measurement=len(qualified),
                            variants_with_any_inactive=len(variants) if inactive_only else 0,
                            plddt_ca=first.get("plddt_ca"), residue_sasa_A2=first.get("residue_sasa_A2"),
                            sidechain_sasa_A2=first.get("sidechain_sasa_A2")))
    return core.build_candidate_spatial_clusters(summary, mapping, core.DEFAULT_CONFIG)


def render_report(partner, wt, audits, overview, by_condition, residue_tables, batch_rankings, ranked, clusters):
    lines = ["# Multi-batch binding analysis", "", f"Binding partner: {partner}. Shared WT length: {len(wt)} residues.",
             "", "## Batch coverage", "", "| Batch | Condition | Barcode block | Matched wells | PASS WT wells | WT KD (nM) | WT QC |",
             "|---|---|---|---:|---:|---:|---|"]
    def cell(v):
        return str(v).replace("|", "\\|").replace("\n", " ")
    for a in audits:
        kd = a["wt_medians"]["KD"]
        lines.append("| " + " | ".join(map(cell, [a["batch_id"],a["condition_id"]+" / "+a["condition"],
                     a["ngs_plate_mapping"]["selected_source_plate_block"] or "Original numbering",
                     a["join"]["joined_lsa_rows"],a["n_pass_WT"],f"{kd*1e9:.4g}" if kd else "—",a["WT_status"]]))+" |")
    lines += ["", "## Integrated evidence", "", f"- Distinct sequence-verified variants excluding WT: {len(overview)}.",
              f"- Observed in one batch only: {sum(r['batches_observed']==1 for r in overview)}.",
              f"- Quantified in two or more batches: {sum(r['batches_with_fold']>=2 for r in overview)}.",
              f"- Variants with opposing KD effects: {sum(r['opposing_KD_directions'] for r in overview)}.",
              f"- Variants with active and inactive evidence: {sum(r['active_inactive_discordance'] for r in overview)}.",
              "", "## Interpretation", "",
              "A fold is the variant median divided by the WT median in the same batch. Both must pass QC. "
              "Within each condition, the center is the median log2 fold transformed back to a ratio. "
              "The study center uses the same calculation across condition centers, giving each represented condition equal weight.",
              "", "Ranges span eligible batch effects. Lower KD indicates stronger binding. KD folds below 0.80 and above 1.25 "
              "mark opposing directions when both occur. Opposing effects or active/inactive discordance suppress the center. "
              "PASS active and Inactive counts, QC-limited records, and single-batch variants remain traceable.",
              "", "Missing observations carry no effect estimate. Different conditions can modify the relative effect of a variant. "
              "Coverage varies by variant, so study centers describe the sampled conditions only. "
              "These descriptive results do not test statistical significance or separate batch effects from condition effects.",
              "", "See variant_overview.csv, variant_by_condition.csv, variant_by_batch.csv and integrated_wells.csv. "
              "Each numbered batch folder contains its full single-batch reports and QC tables."]
    def fold(val):
        return "—" if val is None else f"{val:.3g}"
    lines += ["", "## Integrated improvements and losses", "",
              "Top five eligible single substitutions per direction, ordered by study KD fold (<0.80 or >1.25).",
              "", "| Direction | Variant | KD fold | kon fold | koff fold | ΔΔG (kcal/mol) | Eligible / observed batches | Conditions with fold | Observed ROIs by batch |",
              "|---|---|---:|---:|---:|---:|---|---:|---|"]
    for r in ranked:
        rois = ", ".join(f"{name} ({ids})" for name,ids in json.loads(r["observed_rois_by_batch"]).items())
        values = [r["direction"], r["variant"], *(fold(r[m+"_fold_center"]) for m in METRICS),
                  fold(r["ddG_kcal_mol"]), f"{r['batches_with_fold']} / {r['batches_observed']}", r["conditions_with_fold"], rois]
        lines.append("| "+" | ".join(map(cell, values))+" |")
    if not ranked:
        lines += ["", "No eligible single substitutions passed the study ranking thresholds."]
    lines += ["", "## Integrated spatial clusters", ""]
    if clusters is None:
        lines.append("Supply a shared structure to identify spatial clusters.")
    else:
        lines.append("Candidates have a study KD fold ≥2 for at least one single substitution, or exclusively inactive-only single substitutions at that position. Exact WT matches with Cα coordinates within 8 Å are linked into connected components. Inactive-only evidence is qualitative. Singletons are listed separately from multi-residue regions.")
        if clusters:
            lines += ["", "| Group | Type | WT residues | KD-loss positions | Inactive-only positions | Maximum study KD fold |",
                      "|---|---|---|---:|---:|---:|"]
            for c in clusters:
                values = [c["cluster_id"], "Region" if c["is_multi_position_cluster"] else "Singleton",
                          ", ".join(c["member_labels"]), c["qc_supported_strong_loss_positions"],
                          c["inactive_without_pass_kd_review_positions"], fold(c["maximum_observed_KD_fold_vs_WT"])]
                lines.append("| "+" | ".join(map(cell, values))+" |")
        else:
            lines += ["", "No mapped positions met the study cluster criteria."]
        lines += ["", "See study_cluster_members.csv and study_clusters.csv. Spatial proximity alone cannot establish a binding interface."]
    lines += ["", "## Strongest improvements and losses by batch"]
    for audit, rankings in zip(audits, batch_rankings):
        lines += ["", f"### {cell(audit['batch_id'])} / {cell(audit['condition'])}", "", rankings]
    lines += ["", "## Residues by structure classification", ""]
    if residue_tables:
        labels = {a["condition_id"]: a["condition"] for a in audits}
        for cid, table in residue_tables.items():
            heading = "All batches" if cid == "study" else f"{cid} / {cell(labels[cid])}"
            lines += [f"### {heading}", "", core.report_residue_classes(table), ""]
    else:
        lines.append(core.report_residue_classes(None))
    lines += ["", "## Variant effects across sampled conditions", "",
              "| Variant | Batches observed / total | Eligible batches | Conditions with fold | KD fold | KD range | Interpretation |",
              "|---|---|---:|---:|---:|---|---|"]
    for row in overview:
        values = [row["variant"], f"{row['batches_observed']} / {row['batches_total']}", row["batches_with_fold"],
                  row["conditions_with_fold"], fold(row["KD_fold_center"]),
                  fold(row["KD_fold_min"])+" – "+fold(row["KD_fold_max"]), row["interpretation"]]
        lines.append("| "+" | ".join(map(cell, values))+" |")
    lines += ["", "## Variant effects by condition", "",
              "| Variant | Condition | Observed / eligible batches | KD fold | KD range | kon fold | koff fold | Interpretation |",
              "|---|---|---|---:|---|---:|---:|---|"]
    for row in by_condition:
        values = [row["variant"], row["condition_id"]+" / "+row["condition"],
                  f"{row['batches_observed']} / {row['batches_with_fold']}", fold(row["KD_fold_center"]),
                  fold(row["KD_fold_min"])+" – "+fold(row["KD_fold_max"]), fold(row["kon_fold_center"]),
                  fold(row["koff_fold_center"]), row["interpretation"]]
        lines.append("| "+" | ".join(map(cell, values))+" |")
    if not by_condition:
        lines.append("No sequence-verified non-WT variants.")
    return "\n".join(lines)+"\n"
