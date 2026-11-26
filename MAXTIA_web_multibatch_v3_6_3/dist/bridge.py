"""Browser adapter for the MAXTIA structure mapping engine."""
import csv
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import maxitia_mapping as core
import structure_explorer as explorer


def synthetic_inputs(home, wells=None, source_start=1):
    """Four artificial wells in the same two formats accepted by the CLI."""
    lsa = home / "demo.xlsx"
    ngs = home / "demo.csv"
    headers = list(core.LSA_FIELDS)
    sheet = ET.Element("worksheet", xmlns=core.MAIN_NS)
    data = ET.SubElement(sheet, "sheetData")
    table = [headers]
    if wells is None:
        wells = [
            (1, 1, "C", "WT", "ACDE", "None", 0.001),
            (1, 2, "C", "D3A", "ACAE", "None", 0.0025),
            (1, 3, "C", "E4A", "ACDA", "Inactive", None),
            (1, 4, "M", "D3A", "ACAE", "None", 0.001),
        ]
    for index, (plate, well, call, label, seq, flag, off) in enumerate(wells, 1):
        record = {"ROI": index, "Ligand Name": f"{plate}-{well}", "Flagged Data": flag,
                  "Analyte Set Name": "Demo", "ka (1/Ms)": 10000, "kd (1/s)": off,
                  "KD (M)": off / 10000 if off else None, "Rmax": 100, "Res sd": 1}
        for key in headers:
            if key.endswith(" Error") and record.get(key[:-6]):
                record[key] = 0.01 * record[key[:-6]]
        table.append([record.get(key) for key in headers])
    for index, values in enumerate(table, 1):
        row = ET.SubElement(data, "row", r=str(index))
        for column, value in enumerate(values):
            if value is None:
                continue
            cell = ET.SubElement(row, "c", r=f"{chr(65+column)}{index}")
            if isinstance(value, str):
                cell.set("t", "inlineStr")
                ET.SubElement(ET.SubElement(cell, "is"), "t").text = value
            else:
                ET.SubElement(cell, "v").text = str(value)
    with zipfile.ZipFile(lsa, "w") as book:
        book.writestr("xl/workbook.xml", f'<workbook xmlns="{core.MAIN_NS}" xmlns:r="{core.DOC_REL_NS}"><sheets><sheet name="Kinetic Snapshot" sheetId="1" r:id="r1"/></sheets></workbook>')
        book.writestr("xl/_rels/workbook.xml.rels", f'<Relationships xmlns="{core.PKG_REL_NS}"><Relationship Id="r1" Target="worksheets/sheet1.xml"/></Relationships>')
        book.writestr("xl/worksheets/sheet1.xml", ET.tostring(sheet))
    with ngs.open("w", newline="") as file:
        writer = csv.writer(file)
        for plate, well, call, label, seq, flag, off in wells:
            writer.writerow([f"P{plate+source_start-1:02d}-{well:02d}", call, 120,
                             0.01 if call == "C" else 0.4, 0.02, "", repr([label]), seq])
    return {"lsa": lsa, "ngs": ngs}


def wt_input(text):
    if not text or not text.strip():
        return ""
    lines = text.strip().splitlines()
    if lines[0].startswith(">"):
        lines = lines[1:]
    clean = "".join("".join(lines).split()).upper()
    if not core.AA_RE.fullmatch(clean):
        raise ValueError("Enter one wild-type sequence using standard one-letter amino acid codes.")
    return clean


def analyze(request, home, pending_views=None):
    inputs = {key: Path(name) for key, name in request["paths"].items()}
    if request["demo"]:
        inputs = synthetic_inputs(home)
        inputs["cif"] = explorer.demo_structure(home)
    policy = json.loads(inputs["policy"].read_text(encoding="utf-8-sig")) if "policy" in inputs else {}
    if not isinstance(policy, dict) or {"inputs", "config_overrides", "out_dir", "lsa_xlsx", "ngs_csv"} & policy.keys():
        raise ValueError("Select a QC policy JSON file. Run configurations containing file paths are unsupported.")
    config = core.deep_merge(core.DEFAULT_CONFIG, policy)
    if request.get("study_batch"):
        config["experiment_metadata"] = policy.get("experiment_metadata", {}).copy()
        config["analysis_name"] = request["study_batch"]
    temperature = policy.get("experiment_metadata", {}).get("temperature_C")
    if temperature is not None:
        kelvin = float(temperature) + 273.15
        if kelvin <= 0 or ("temperature_kelvin" in policy and abs(float(policy["temperature_kelvin"]) - kelvin) > 0.1):
            raise ValueError("Check the Celsius and Kelvin temperatures in the QC policy.")
        config["temperature_kelvin"] = kelvin
    if not (core.to_float(config["temperature_kelvin"]) or 0) > 0:
        raise ValueError("Enter a finite positive absolute temperature in kelvin.")
    if request.get("block"):
        config["ngs_plate_block"] = request["block"]
        config["ngs_target_plate_start"] = 1
    config["require_explicit_ngs_plate_block_for_superset"] = True
    config["wt_sequence"] = wt_input(request.get("wt") or config.get("wt_sequence", ""))
    chain = request.get("chain") or config["structure_chain_id"]
    if not isinstance(chain, str) or not re.fullmatch(r"[A-Za-z0-9_]+", chain):
        raise ValueError("Enter a chain ID using letters, digits, or underscores.")
    config["structure_chain_id"] = chain
    config["structure_visualization"]["render_static_pngs_if_pymol_available"] = False
    config["data_flagging_settings_confirmed_for_experiment"] = policy.get("data_flagging_settings_confirmed_for_experiment") is True
    config["inactive_response_threshold_confirmed_for_experiment"] = policy.get("inactive_response_threshold_confirmed_for_experiment") is True
    if "cif" not in inputs and ("full" in inputs or "summary" in inputs):
        raise ValueError("Select the structure CIF file from the same model as the AlphaFold 3 JSON files.")
    for key, required in [("full", {"pae", "contact_probs", "token_chain_ids", "token_res_ids"}), ("summary", {"ptm"})]:
        if key in inputs:
            obj = json.loads(inputs[key].read_text())
            if not isinstance(obj, dict) or not required.issubset(obj):
                raise ValueError(f"Check that the selected file is the AlphaFold 3 {key} JSON for this model.")
            if key == "full" and chain not in obj["token_chain_ids"]:
                raise ValueError("The selected chain ID was not found in the full data JSON.")
    results = home / "results"
    report = core.run(inputs["lsa"], inputs["ngs"], results, config, inputs.get("cif"),
                      inputs.get("full"), inputs.get("summary"), chain, pending_views=pending_views)
    if not report["input_shape"]["lsa_kinetic_rows"] or not report["key_join"]["joined_lsa_rows"]:
        raise ValueError("No LSA rows matched the NGS records. Check the input files and primer barcode block.")
    if "cif" in inputs and not report["wt_sequence"]["length"]:
        raise ValueError("The NGS records do not establish a wild-type sequence. Enter the WT sequence for structure mapping.")
    (results / "effective_policy.json").write_text(json.dumps(config, indent=2))
    # The CLI's legacy report describes its original apparatus settings in the
    # past tense. Keep file-derived observations distinct from policy assumptions.
    def settings_note(line):
        markers = ("The confirmed Carterra", "The run used a", "The run enabled",
                   "`Data <50% Rmax`, `kd Below Limit`", "The 10 RU Data Filter")
        return not line.lstrip("- ").startswith(markers)
    qc_path = results / "qc_report.md"
    qc_path.write_text("\n".join(filter(settings_note, qc_path.read_text().splitlines())) +
                       "\n\nConsult effective_policy.json for configured filter, highlighting, and fit settings. "
                       "Instrument flags are preserved as exported.\n")
    report["provisional_assumptions"] = list(filter(settings_note, report["provisional_assumptions"]))
    (results / "qc_summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2))
    structure = report.get("structure", {})
    viewer = None
    if pending_views is None and (results / "residue_colors.csv").exists():
        table = explorer.read_rows(results / "residue_colors.csv")
        viewer = dict(geometry=explorer.geometry(inputs["cif"], chain),
                      scopes=[explorer.scope("Experiment", table, explorer.single_rows(results))])
        explorer.write_html(results / "structure_explorer.html", viewer)
    core.write_run_manifest(results, report, config)
    with zipfile.ZipFile(home / "results.zip", "w", zipfile.ZIP_DEFLATED) as archive:
        for file in sorted(results.iterdir()):
            archive.write(file, file.name)
    return json.dumps({"viewer": viewer, "joined": report["key_join"]["joined_lsa_rows"],
                       "rows": report["input_shape"]["lsa_kinetic_rows"],
                       "quantitative": report["qc"]["analysis_eligibility_counts"].get("QUANTITATIVE", 0),
                       "inactive": report["qc"]["analysis_eligibility_counts"].get("CENSORED_INACTIVE", 0),
                       "mapped": structure.get("mapped_positions"),
                       "block": report["ngs_plate_mapping"]["selected_source_plate_block"],
                       "warnings": report.get("parser_warnings", []),
                       "report": (results / "analysis_report.md").read_text()})
