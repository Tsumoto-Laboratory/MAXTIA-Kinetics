#!/usr/bin/env python3
"""MAXTIA/Carterra LSA kinetic results + headerless NGS library output integrator.

The parser intentionally reads XLSX cell XML directly. This avoids a known issue in
the supplied Carterra workbook where a malformed style color prevents openpyxl from
opening an otherwise readable workbook.

Python 3.11+; standard library only.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path, PurePosixPath
from typing import Any, Iterable
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"m": MAIN_NS, "r": DOC_REL_NS, "p": PKG_REL_NS}
MODULE_VERSION = "3.6.3"
GAS_CONSTANT_KCAL = 1.987204258e-3
STRUCTURE_COLORS = {
    "NO_DATA": "#808080",
    "REVIEW": "#800080",
    "INACTIVE_ONLY": "#FFFF00",
    "NO_PASS_DATA": "#FF00FF",
}
STRUCTURE_CLASS_LABELS = {
    "NO_DATA": "No data",
    "REVIEW": "Review",
    "INACTIVE_ONLY": "Inactive only",
    "NO_PASS_DATA": "No PASS data",
}

DEFAULT_CONFIG: dict[str, Any] = {
    "qc_policy_name": "carterra_1.9.2_ngs_defined_evidence_preserving_v3_4",
    "analysis_name": "MAXTIA Kinetic-to-Structure Analysis",
    "carterra_kinetics_software_version": "1.9.2.4463",
    "data_flagging_settings_confirmed_for_experiment": True,
    "lsa_sheet": "Kinetic Snapshot",
    "kinetic_ligand_pattern": r"^(?P<plate>\d+)-(?P<well>\d+)$",
    "ngs_well_plate_prefix": "P",
    "ngs_plate_width": 2,
    "ngs_well_width": 2,
    "ngs_primer_barcode_type": "",
    "ngs_plate_block": None,
    "ngs_target_plate_start": 1,
    "require_explicit_ngs_plate_block_for_superset": True,
    "ngs_has_header": False,
    "ngs_columns": [
        "well",
        "call_code",
        "total_valid_read_count",
        "inter_cluster_mixing_score",
        "intra_cluster_mixing_score",
        "dna_sequence",
        "mutation_labels_raw",
        "aa_sequence",
    ],
    "ngs_cluster_identity_threshold": 0.98,
    "ngs_inter_cluster_mixing_threshold": 0.10,
    "ngs_intra_cluster_mixing_threshold": 0.20,
    "ngs_score_threshold_comparison": "strict_less_than",
    "ngs_score_boundary_handling": "unknown_keep_authoritative_call_and_review",
    "ngs_read_count_definition": "total_valid_reads",
    "ngs_call_code_is_authoritative": True,
    "ngs_pcr_induced_mutation_possible": True,
    "accepted_ngs_call_codes": ["C"],
    "lsa_pass_flags": ["None"],
    "lsa_review_flags": [],
    "lsa_complex_flags": ["Complex"],
    "lsa_inactive_flags": ["Inactive"],
    "lsa_fail_flags": [],
    "relative_error_review_threshold": 0.5,
    "relative_error_fail_threshold": None,
    "require_fit_errors_for_active": True,
    "zero_fit_error_is_fail": True,
    "require_residual_sd_for_active": True,
    "residual_sd_over_rmax_fail_threshold": 0.10,
    "rmax_review_threshold_RU": None,
    "inactive_default_min_binding_response_RU": 20.0,
    "inactive_response_threshold_confirmed_for_experiment": True,
    "carterra_data_filter": {
        "exclude_low_or_no_binding_enabled": True,
        "exclude_low_or_no_binding_threshold_RU": 10.0,
        "exclude_poor_1_to_1_fits_enabled": True,
        "poor_fit_residual_sd_percent_maxbound": 10.0,
        "fitting_assist_concentration_selective_enabled": False,
        "fitting_assist_minimum_n_concentrations": 3,
        "fitting_assist_residual_sd_percent_maxbound": 10.0,
    },
    "carterra_highlighting": {
        "response_below_RU_enabled": True,
        "response_below_RU": 20.0,
        "residual_sd_percent_rmax_enabled": True,
        "residual_sd_percent_rmax": 10.0,
        "excluded_group_highlighting_enabled": True,
        "manual_excluded_group_used": False,
        "data_below_50_percent_rmax_enabled": True,
        "kd_below_limit_enabled": True,
        "km_below_limit_enabled": True,
    },
    "fit_parameter_bounds": {
        "kon_M_inv_s": {"min": 100.0, "max": None},
        "koff_s_inv": {"min": 1e-5, "max": None},
        "Rmax_RU": {"min": 0.0, "max": 20000.0},
    },
    "fit_boundary_relative_tolerance": 1e-4,
    "fit_boundary_absolute_tolerance": 1e-12,
    "kd_identity_relative_tolerance": 1e-6,
    "per_fit_kd_error_propagation_check_enabled": False,
    "kd_error_propagation_relative_tolerance": 0.10,
    "require_positive_rmax": True,
    "minimum_active_pass_replicates_for_quantitation": 1,
    "mixed_active_inactive_policy": "use_pass_active_review_exclude_inactive_from_loss_evidence",
    "replicate_consistency_review_fold": None,
    "replicate_consistency_fail_fold": None,
    "ngs_min_read_count_review": None,
    "colony_pcr_failure_inference_enabled": False,
    "require_explicit_substitution_labels": True,
    "temperature_kelvin": 298.15,
    "experiment_metadata": {
        "manual_reference": "Carterra Kinetics Data Analysis Software User Manual, May 2024",
        "fit_model": "1:1 Langmuir",
        "fit_scope": "concentration series global fit",
        "rmax_mode": "Global Rmax",
        "temperature_C": 25.0,
        "association_time_s": 300.0,
        "dissociation_time_s": 300.0,
        "analyte_concentrations_M": [6.25e-8, 1.25e-7, 2.5e-7, 5e-7, 1e-6, 2e-6],
        "blank_subtraction": True,
        "double_reference": True,
        "ligand_molecular_weight_kDa": 15.0,
        "analyte_molecular_weight_kDa": 26.0,
        "binding_stoichiometry": "1:1",
        "immobilization_level_reliable": False,
        "expression_or_spot_validity_control": False,
    },
    "visualization_colors": {
        "PASS": "#2E7D32",
        "REVIEW": "#F9A825",
        "FAIL": "#C62828",
        "INACTIVE": "#757575",
        "MIXED_REVIEW": "#7B1FA2",
    },
    "structure_visualization": {
        "enabled": True,
        "render_static_pngs_if_pymol_available": False,
        "static_png_width": 1400,
        "static_png_height": 1000,
        "strong_improvement_max_KD_fold": 0.5,
        "improvement_max_KD_fold": 0.8,
        "wt_like_max_KD_fold": 1.25,
        "moderate_loss_max_KD_fold": 2.0,
        "candidate_cluster_ca_cutoff_A": 8.0,
        "colors": STRUCTURE_COLORS.copy(),
    },
    "wt_sequence": "",
    "structure_chain_id": "A",
    "sasa_probe_radius_A": 1.4,
    "sasa_sphere_points": 240,
}

LSA_FIELDS = {
    "ROI": "roi",
    "Bay": "bay",
    "Position": "position",
    "Ligand Name": "ligand_name",
    "Analyte Set Name": "analyte_set_name",
    "Flagged Data": "lsa_flag",
    "ka (1/Ms)": "kon_M_inv_s",
    "ka (1/Ms) Error": "kon_fit_error_M_inv_s",
    "kd (1/s)": "koff_s_inv",
    "kd (1/s) Error": "koff_fit_error_s_inv",
    "KD (M)": "KD_M",
    "KD (M) Error": "KD_fit_error_M",
    "Rmax": "Rmax_RU",
    "Rmax Error": "Rmax_fit_error_RU",
    "Res sd": "residual_sd_RU",
}

NUMERIC_LSA_FIELDS = {
    "kon_M_inv_s",
    "kon_fit_error_M_inv_s",
    "koff_s_inv",
    "koff_fit_error_s_inv",
    "KD_M",
    "KD_fit_error_M",
    "Rmax_RU",
    "Rmax_fit_error_RU",
    "residual_sd_RU",
}

SUBSTITUTION_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")
AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}

VDW_RADII_A = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}


def deep_merge(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def col_number(cell_reference: str) -> int:
    match = re.match(r"([A-Z]+)", cell_reference)
    if not match:
        raise ValueError(f"Invalid cell reference: {cell_reference!r}")
    result = 0
    for character in match.group(1):
        result = result * 26 + ord(character) - 64
    return result


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    return [
        "".join(node.text or "" for node in item.iter(f"{{{MAIN_NS}}}t"))
        for item in root.findall("m:si", NS)
    ]


def _workbook_sheets(archive: zipfile.ZipFile) -> dict[str, str]:
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    relation_targets = {
        relation.attrib["Id"]: relation.attrib["Target"]
        for relation in relations.findall("p:Relationship", NS)
    }
    result: dict[str, str] = {}
    for sheet in workbook.findall("m:sheets/m:sheet", NS):
        relation_id = sheet.attrib[f"{{{DOC_REL_NS}}}id"]
        target = relation_targets[relation_id].lstrip("/")
        path = str(PurePosixPath("xl") / target)
        result[sheet.attrib["name"]] = path
    return result


def _cell_value(cell: ET.Element, strings: list[str]) -> Any:
    cell_type = cell.attrib.get("t")
    value_node = cell.find("m:v", NS)
    inline_node = cell.find("m:is", NS)
    if value_node is None and inline_node is None:
        return None
    if cell_type == "s":
        return strings[int(value_node.text or "0")]
    if cell_type == "inlineStr":
        return "".join(
            node.text or "" for node in inline_node.iter(f"{{{MAIN_NS}}}t")
        )
    if cell_type == "b":
        return value_node is not None and value_node.text == "1"
    if cell_type in {"str", "e"}:
        return value_node.text if value_node is not None else ""
    raw = value_node.text if value_node is not None else ""
    try:
        number = float(raw)
        return int(number) if number.is_integer() else number
    except (TypeError, ValueError):
        return raw


def read_xlsx_sheet_rows(path: Path, sheet_name: str) -> list[tuple[int, dict[int, Any]]]:
    """Read an XLSX worksheet without parsing the style table."""
    with zipfile.ZipFile(path) as archive:
        strings = _shared_strings(archive)
        sheets = _workbook_sheets(archive)
        if sheet_name not in sheets:
            raise ValueError(
                f"Worksheet {sheet_name!r} not found. Available: {', '.join(sheets)}"
            )
        root = ET.fromstring(archive.read(sheets[sheet_name]))
        rows: list[tuple[int, dict[int, Any]]] = []
        for row in root.findall("m:sheetData/m:row", NS):
            row_number = int(row.attrib.get("r", len(rows) + 1))
            cells = {
                col_number(cell.attrib["r"]): _cell_value(cell, strings)
                for cell in row.findall("m:c", NS)
            }
            rows.append((row_number, cells))
        return rows


def normalize_header(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip())


def read_cif_loop(path: Path, category_prefix: str) -> list[dict[str, str]]:
    """Read one simple mmCIF loop category using only the standard library."""
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        if lines[index].strip() != "loop_":
            index += 1
            continue
        index += 1
        headers: list[str] = []
        while index < len(lines) and lines[index].strip().startswith("_"):
            headers.append(lines[index].strip())
            index += 1
        if not headers or not all(header.startswith(category_prefix) for header in headers):
            while index < len(lines):
                stripped = lines[index].strip()
                if stripped == "#":
                    index += 1
                    break
                if stripped == "loop_" or stripped.startswith("_"):
                    break
                index += 1
            continue

        records: list[dict[str, str]] = []
        pending: list[str] = []
        while index < len(lines):
            stripped = lines[index].strip()
            if stripped == "#":
                index += 1
                break
            if stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break
            if stripped:
                if stripped.startswith(";"):
                    raise ValueError(
                        f"Multiline semicolon values are not supported in {category_prefix}"
                    )
                pending.extend(shlex.split(stripped, posix=True))
                while len(pending) >= len(headers):
                    row = pending[: len(headers)]
                    del pending[: len(headers)]
                    records.append(dict(zip(headers, row)))
            index += 1
        if pending:
            raise ValueError(f"Incomplete mmCIF row in category {category_prefix}")
        return records
    raise ValueError(f"mmCIF category {category_prefix!r} not found in {path}")


def cif_int(value: str | None) -> int | None:
    if value in (None, "", ".", "?"):
        return None
    return int(str(value))


def cif_float(value: str | None) -> float | None:
    if value in (None, "", ".", "?"):
        return None
    return float(str(value))


def parse_structure_cif(path: Path, chain_id: str) -> dict[str, Any]:
    sequence_rows = read_cif_loop(path, "_entity_poly_seq.")
    atom_rows = read_cif_loop(path, "_atom_site.")
    chain_atoms = [
        row
        for row in atom_rows
        if row.get("_atom_site.group_PDB") == "ATOM"
        and (
            row.get("_atom_site.auth_asym_id") == chain_id
            or row.get("_atom_site.label_asym_id") == chain_id
        )
        and row.get("_atom_site.pdbx_PDB_model_num", "1") == "1"
    ]
    if not chain_atoms:
        available = sorted(
            {
                row.get("_atom_site.auth_asym_id") or row.get("_atom_site.label_asym_id")
                for row in atom_rows
                if row.get("_atom_site.group_PDB") == "ATOM"
            }
        )
        raise ValueError(f"Chain {chain_id!r} not found. Available chains: {available}")

    entity_ids = {row.get("_atom_site.label_entity_id") for row in chain_atoms}
    if len(entity_ids) != 1:
        raise ValueError(f"Chain {chain_id!r} maps to multiple entities: {entity_ids}")
    entity_id = next(iter(entity_ids))
    polymer_rows = sorted(
        [row for row in sequence_rows if row.get("_entity_poly_seq.entity_id") == entity_id],
        key=lambda row: int(row["_entity_poly_seq.num"]),
    )
    polymer_sequence = "".join(
        THREE_TO_ONE.get(row["_entity_poly_seq.mon_id"].upper(), "X")
        for row in polymer_rows
    )

    atoms: list[dict[str, Any]] = []
    residues: dict[int, dict[str, Any]] = {}
    for row in chain_atoms:
        label_seq_id = cif_int(row.get("_atom_site.label_seq_id"))
        if label_seq_id is None:
            continue
        atom = {
            "element": row.get("_atom_site.type_symbol", "").upper(),
            "atom_name": row.get("_atom_site.label_atom_id", ""),
            "comp_id": row.get("_atom_site.label_comp_id", "").upper(),
            "label_seq_id": label_seq_id,
            "auth_seq_id": cif_int(row.get("_atom_site.auth_seq_id")),
            "label_chain_id": row.get("_atom_site.label_asym_id", ""),
            "auth_chain_id": row.get("_atom_site.auth_asym_id", ""),
            "x": cif_float(row.get("_atom_site.Cartn_x")),
            "y": cif_float(row.get("_atom_site.Cartn_y")),
            "z": cif_float(row.get("_atom_site.Cartn_z")),
            "plddt": cif_float(row.get("_atom_site.B_iso_or_equiv")),
        }
        if None in (atom["x"], atom["y"], atom["z"]):
            continue
        atoms.append(atom)
        residue = residues.setdefault(
            label_seq_id,
            {
                "label_seq_id": label_seq_id,
                "auth_seq_id": atom["auth_seq_id"],
                "auth_chain_id": atom["auth_chain_id"],
                "insertion_code": row.get("_atom_site.pdbx_PDB_ins_code", ""),
                "comp_id": atom["comp_id"],
                "one_letter": THREE_TO_ONE.get(atom["comp_id"], "X"),
                "atoms": [],
            },
        )
        residue["atoms"].append(atom)
    observed_sequence = "".join(residues[index]["one_letter"] for index in sorted(residues))
    return {
        "chain_id": chain_id,
        "entity_id": entity_id,
        "polymer_sequence": polymer_sequence,
        "observed_sequence": observed_sequence,
        "polymer_rows": polymer_rows,
        "residues": residues,
        "atoms": atoms,
    }


def global_alignment(reference: str, target: str) -> dict[str, Any]:
    """Needleman–Wunsch mapping from experimental sequence positions to structure positions."""
    match_score, mismatch_score, gap_score = 2, -1, -2
    rows, cols = len(reference) + 1, len(target) + 1
    score = [[0] * cols for _ in range(rows)]
    trace = [[""] * cols for _ in range(rows)]
    for i in range(1, rows):
        score[i][0] = i * gap_score
        trace[i][0] = "U"
    for j in range(1, cols):
        score[0][j] = j * gap_score
        trace[0][j] = "L"
    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = score[i - 1][j - 1] + (
                match_score if reference[i - 1] == target[j - 1] else mismatch_score
            )
            up = score[i - 1][j] + gap_score
            left = score[i][j - 1] + gap_score
            best = max(diagonal, up, left)
            score[i][j] = best
            trace[i][j] = "D" if best == diagonal else ("U" if best == up else "L")
    aligned_reference: list[str] = []
    aligned_target: list[str] = []
    i, j = len(reference), len(target)
    while i or j:
        direction = trace[i][j]
        if direction == "D":
            aligned_reference.append(reference[i - 1])
            aligned_target.append(target[j - 1])
            i -= 1
            j -= 1
        elif direction == "U":
            aligned_reference.append(reference[i - 1])
            aligned_target.append("-")
            i -= 1
        else:
            aligned_reference.append("-")
            aligned_target.append(target[j - 1])
            j -= 1
    aligned_reference.reverse()
    aligned_target.reverse()

    mapping: dict[int, int | None] = {}
    reference_position = 0
    target_position = 0
    matches = mismatches = gaps_in_reference = gaps_in_target = 0
    for ref_aa, target_aa in zip(aligned_reference, aligned_target):
        if ref_aa != "-":
            reference_position += 1
        if target_aa != "-":
            target_position += 1
        if ref_aa == "-":
            gaps_in_reference += 1
        elif target_aa == "-":
            gaps_in_target += 1
            mapping[reference_position] = None
        else:
            mapping[reference_position] = target_position
            if ref_aa == target_aa:
                matches += 1
            else:
                mismatches += 1
    denominator = max(1, matches + mismatches)
    return {
        "mapping": mapping,
        "aligned_reference": "".join(aligned_reference),
        "aligned_target": "".join(aligned_target),
        "matches": matches,
        "mismatches": mismatches,
        "gaps_in_reference": gaps_in_reference,
        "gaps_in_target": gaps_in_target,
        "identity_over_aligned_residues": matches / denominator,
    }


def fibonacci_sphere(number_of_points: int) -> list[tuple[float, float, float]]:
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    points: list[tuple[float, float, float]] = []
    for index in range(number_of_points):
        y = 1.0 - 2.0 * (index + 0.5) / number_of_points
        radius = math.sqrt(max(0.0, 1.0 - y * y))
        angle = golden_angle * index
        points.append((math.cos(angle) * radius, y, math.sin(angle) * radius))
    return points


def squared_distance(a: dict[str, Any], b: dict[str, Any]) -> float:
    return (
        (a["x"] - b["x"]) ** 2
        + (a["y"] - b["y"]) ** 2
        + (a["z"] - b["z"]) ** 2
    )


def compute_sasa(
    atoms: list[dict[str, Any]], probe_radius: float, sphere_points: int
) -> dict[int, dict[str, float]]:
    unit_points = fibonacci_sphere(sphere_points)
    expanded_radii = [
        VDW_RADII_A.get(atom["element"], 1.70) + probe_radius for atom in atoms
    ]
    neighbors: list[list[int]] = [[] for _ in atoms]
    for left in range(len(atoms)):
        for right in range(left + 1, len(atoms)):
            cutoff = expanded_radii[left] + expanded_radii[right]
            if squared_distance(atoms[left], atoms[right]) < cutoff * cutoff:
                neighbors[left].append(right)
                neighbors[right].append(left)
    by_residue: dict[int, dict[str, float]] = defaultdict(
        lambda: {"residue_sasa_A2": 0.0, "sidechain_sasa_A2": 0.0}
    )
    for index, atom in enumerate(atoms):
        radius = expanded_radii[index]
        accessible = 0
        for unit_x, unit_y, unit_z in unit_points:
            point_x = atom["x"] + radius * unit_x
            point_y = atom["y"] + radius * unit_y
            point_z = atom["z"] + radius * unit_z
            occluded = False
            for neighbor_index in neighbors[index]:
                neighbor = atoms[neighbor_index]
                neighbor_radius = expanded_radii[neighbor_index]
                if (
                    (point_x - neighbor["x"]) ** 2
                    + (point_y - neighbor["y"]) ** 2
                    + (point_z - neighbor["z"]) ** 2
                    < neighbor_radius * neighbor_radius
                ):
                    occluded = True
                    break
            if not occluded:
                accessible += 1
        area = 4.0 * math.pi * radius * radius * accessible / sphere_points
        residue_id = atom["label_seq_id"]
        by_residue[residue_id]["residue_sasa_A2"] += area
        if atom["atom_name"] not in BACKBONE_ATOMS:
            by_residue[residue_id]["sidechain_sasa_A2"] += area
    return dict(by_residue)


def compute_geometric_contacts(residues: dict[int, dict[str, Any]]) -> dict[int, dict[str, int]]:
    output = {
        residue_id: {"ca_contact_count_8A_nonlocal": 0, "heavy_contact_count_5A_nonlocal": 0}
        for residue_id in residues
    }
    residue_ids = sorted(residues)
    for left_index, left_id in enumerate(residue_ids):
        left = residues[left_id]
        left_ca = next((atom for atom in left["atoms"] if atom["atom_name"] == "CA"), None)
        for right_id in residue_ids[left_index + 1 :]:
            if abs(left_id - right_id) <= 2:
                continue
            right = residues[right_id]
            right_ca = next((atom for atom in right["atoms"] if atom["atom_name"] == "CA"), None)
            if left_ca and right_ca and squared_distance(left_ca, right_ca) <= 64.0:
                output[left_id]["ca_contact_count_8A_nonlocal"] += 1
                output[right_id]["ca_contact_count_8A_nonlocal"] += 1
            has_heavy_contact = any(
                squared_distance(left_atom, right_atom) <= 25.0
                for left_atom in left["atoms"]
                for right_atom in right["atoms"]
            )
            if has_heavy_contact:
                output[left_id]["heavy_contact_count_5A_nonlocal"] += 1
                output[right_id]["heavy_contact_count_5A_nonlocal"] += 1
    return output


def load_af3_token_features(
    full_data_path: Path | None, chain_id: str
) -> tuple[dict[int, dict[str, float]], dict[str, Any]]:
    if full_data_path is None:
        return {}, {}
    data = json.loads(full_data_path.read_text(encoding="utf-8"))
    chain_ids = data.get("token_chain_ids", [])
    residue_ids = data.get("token_res_ids", [])
    pae = data.get("pae", [])
    contact_probs = data.get("contact_probs", [])
    if not (len(chain_ids) == len(residue_ids) == len(pae)):
        raise ValueError("AF3 token arrays have inconsistent lengths")
    chain_indices = [index for index, value in enumerate(chain_ids) if value == chain_id]
    features: dict[int, dict[str, float]] = {}
    for index in chain_indices:
        residue_id = int(residue_ids[index])
        same_chain = chain_indices
        local = [other for other in same_chain if abs(int(residue_ids[other]) - residue_id) <= 5]
        nonlocal_indices = [
            other for other in same_chain if abs(int(residue_ids[other]) - residue_id) > 2
        ]
        pae_values = [float(pae[index][other]) for other in same_chain]
        local_pae_values = [float(pae[index][other]) for other in local]
        contact_values = [float(contact_probs[index][other]) for other in nonlocal_indices]
        features[residue_id] = {
            "pae_mean_to_selected_chain": statistics.mean(pae_values),
            "pae_max_to_selected_chain": max(pae_values),
            "pae_local_mean_window_5": statistics.mean(local_pae_values),
            "predicted_contact_mean_nonlocal": (
                statistics.mean(contact_values) if contact_values else 0.0
            ),
            "predicted_contact_count_ge_0_8_nonlocal": sum(
                value >= 0.8 for value in contact_values
            ),
        }
    return features, {
        "token_count_all_chains": len(chain_ids),
        "token_count_selected_chain": len(chain_indices),
        "selected_chain_token_residue_ids": [int(residue_ids[index]) for index in chain_indices],
    }


def finite_stats(values: Iterable[float | None]) -> dict[str, float | int | None]:
    clean = [float(value) for value in values if value is not None and math.isfinite(value)]
    if not clean:
        return {"count": 0, "min": None, "median": None, "mean": None, "max": None}
    return {
        "count": len(clean),
        "min": min(clean),
        "median": statistics.median(clean),
        "mean": statistics.mean(clean),
        "max": max(clean),
    }


def build_structure_mapping(
    wt_sequence: str,
    structure_cif: Path,
    full_data_json: Path | None,
    summary_confidences_json: Path | None,
    chain_id: str,
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    structure = parse_structure_cif(structure_cif, chain_id)
    alignment = global_alignment(wt_sequence, structure["polymer_sequence"])
    sasa = compute_sasa(
        structure["atoms"],
        float(config["sasa_probe_radius_A"]),
        int(config["sasa_sphere_points"]),
    )
    contacts = compute_geometric_contacts(structure["residues"])
    token_features, token_summary = load_af3_token_features(full_data_json, chain_id)
    mapping_rows: list[dict[str, Any]] = []
    for sequence_position in range(1, len(wt_sequence) + 1):
        structure_position = alignment["mapping"].get(sequence_position)
        residue = structure["residues"].get(structure_position) if structure_position else None
        ca_atom = (
            next((atom for atom in residue["atoms"] if atom["atom_name"] == "CA"), None)
            if residue
            else None
        )
        atom_plddts = [atom["plddt"] for atom in residue["atoms"] if atom["plddt"] is not None] if residue else []
        structure_aa = residue["one_letter"] if residue else ""
        status = (
            "gap_in_structure"
            if residue is None
            else ("exact" if wt_sequence[sequence_position - 1] == structure_aa else "mismatch")
        )
        mapping_rows.append(
            {
                "sequence_position_1based": sequence_position,
                "sequence_aa": wt_sequence[sequence_position - 1],
                "structure_chain_id": chain_id,
                "structure_auth_chain_id": (residue.get("auth_chain_id") or chain_id) if residue else chain_id,
                "structure_insertion_code": (residue["insertion_code"] if residue and residue["insertion_code"] not in (".","?") else ""),
                "structure_label_seq_id": residue["label_seq_id"] if residue else "",
                "structure_auth_seq_id": residue["auth_seq_id"] if residue else "",
                "structure_residue_number": residue["auth_seq_id"] if residue else "",
                "structure_comp_id": residue["comp_id"] if residue else "",
                "structure_aa": structure_aa,
                "alignment_status": status,
                "ca_x_A": ca_atom["x"] if ca_atom else None,
                "ca_y_A": ca_atom["y"] if ca_atom else None,
                "ca_z_A": ca_atom["z"] if ca_atom else None,
                "plddt_ca": ca_atom["plddt"] if ca_atom else None,
                "plddt_mean_all_atoms": statistics.mean(atom_plddts) if atom_plddts else None,
                **sasa.get(structure_position, {}),
                **contacts.get(structure_position, {}),
                **token_features.get(structure_position, {}),
            }
        )
    confidence_summary = (
        json.loads(summary_confidences_json.read_text(encoding="utf-8"))
        if summary_confidences_json
        else {}
    )
    plddt_values = [row["plddt_ca"] for row in mapping_rows]
    structure_summary = {
        "source_files": {
            "structure_cif": {"name": structure_cif.name, "sha256": sha256(structure_cif)},
            "full_data_json": (
                {"name": full_data_json.name, "sha256": sha256(full_data_json)}
                if full_data_json
                else None
            ),
            "summary_confidences_json": (
                {"name": summary_confidences_json.name, "sha256": sha256(summary_confidences_json)}
                if summary_confidences_json
                else None
            ),
        },
        "chain_id": chain_id,
        "entity_id": structure["entity_id"],
        "experimental_sequence_length": len(wt_sequence),
        "structure_polymer_sequence_length": len(structure["polymer_sequence"]),
        "structure_observed_residue_count": len(structure["residues"]),
        "structure_atom_count": len(structure["atoms"]),
        "structure_polymer_sequence": structure["polymer_sequence"],
        "sequence_alignment": {
            key: alignment[key]
            for key in [
                "matches",
                "mismatches",
                "gaps_in_reference",
                "gaps_in_target",
                "identity_over_aligned_residues",
                "aligned_reference",
                "aligned_target",
            ]
        },
        "mapped_positions": sum(row["alignment_status"] != "gap_in_structure" for row in mapping_rows),
        "exact_positions": sum(row["alignment_status"] == "exact" for row in mapping_rows),
        "plddt_ca": finite_stats(plddt_values),
        "plddt_ca_counts": {
            "lt_70": sum(value is not None and value < 70 for value in plddt_values),
            "70_to_lt_90": sum(value is not None and 70 <= value < 90 for value in plddt_values),
            "ge_90": sum(value is not None and value >= 90 for value in plddt_values),
        },
        "sasa_method": {
            "algorithm": "Shrake-Rupley Fibonacci sphere",
            "probe_radius_A": float(config["sasa_probe_radius_A"]),
            "sphere_points_per_atom": int(config["sasa_sphere_points"]),
        },
        "af3_summary_confidences": confidence_summary,
        "af3_full_data": token_summary,
    }
    return mapping_rows, structure_summary


def merge_residue_effects_with_structure(
    effects: list[dict[str, Any]], mapping_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    by_position = {row["sequence_position_1based"]: row for row in mapping_rows}
    merged: list[dict[str, Any]] = []
    for effect in effects:
        structure = by_position.get(effect["sequence_position_1based"], {})
        merged.append({**effect, **structure})
    return merged


def summarize_mapped_positions(mapped_effects: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for effect in mapped_effects:
        if effect["effect_attribution"] == "direct_single_substitution":
            groups[int(effect["sequence_position_1based"])].append(effect)
    output: list[dict[str, Any]] = []
    for position, effects in sorted(groups.items()):
        measured = [
            effect for effect in effects if effect.get("KD_fold_vs_WT_pass") is not None
        ]
        best = (
            min(measured, key=lambda effect: float(effect["KD_fold_vs_WT_pass"]))
            if measured
            else None
        )
        worst = (
            max(measured, key=lambda effect: float(effect["KD_fold_vs_WT_pass"]))
            if measured
            else None
        )
        first = effects[0]
        active_replicates = sum(
            int(effect.get("n_active_quantitative_pass") or 0) for effect in effects
        )
        inactive_replicates = sum(int(effect.get("n_inactive") or 0) for effect in effects)
        inactive_loss_evidence_replicates = sum(
            int(effect.get("n_inactive_loss_of_binding_evidence") or 0)
            for effect in effects
        )
        inactive_excluded_by_active_override = sum(
            int(effect.get("n_inactive_excluded_by_pass_active_override") or 0)
            for effect in effects
        )
        informative_replicates = active_replicates + inactive_replicates
        loss_evidence_denominator = (
            active_replicates + inactive_loss_evidence_replicates
        )
        output.append(
            {
                "sequence_position_1based": position,
                "wt_aa": first["wt_aa"],
                "single_substitution_variants": len(effects),
                "variants_with_pass_measurement": len(measured),
                "variants_with_any_inactive": sum(
                    int(effect.get("n_inactive") or 0) > 0 for effect in effects
                ),
                "consistently_inactive_variants": sum(
                    effect.get("phenotype_consensus") == "consistently_inactive"
                    for effect in effects
                ),
                "mixed_active_inactive_variants": sum(
                    effect.get("phenotype_consensus") == "mixed_active_inactive"
                    for effect in effects
                ),
                "active_replicates": active_replicates,
                "inactive_replicates_raw": inactive_replicates,
                "inactive_replicates_loss_of_binding_evidence": (
                    inactive_loss_evidence_replicates
                ),
                "inactive_replicates_excluded_by_pass_active_override": (
                    inactive_excluded_by_active_override
                ),
                "position_inactive_fraction_raw": safe_ratio(
                    inactive_replicates, informative_replicates
                ),
                "position_inactive_fraction_for_loss_evidence": safe_ratio(
                    inactive_loss_evidence_replicates,
                    loss_evidence_denominator,
                ),
                "plddt_ca": first.get("plddt_ca"),
                "residue_sasa_A2": first.get("residue_sasa_A2"),
                "sidechain_sasa_A2": first.get("sidechain_sasa_A2"),
                "best_KD_mutation": best["mutation_label"] if best else None,
                "best_KD_fold_vs_WT": best["KD_fold_vs_WT_pass"] if best else None,
                "worst_KD_mutation": worst["mutation_label"] if worst else None,
                "worst_KD_fold_vs_WT": worst["KD_fold_vs_WT_pass"] if worst else None,
            }
        )
    return output


def build_candidate_spatial_clusters(
    position_summary: list[dict[str, Any]],
    mapping_rows: list[dict[str, Any]],
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Group structure-mapped loss candidates by connected C-alpha proximity.

    A position enters the graph when it has a QC-supported KD fold at or above
    the configured strong-loss boundary, or when it has Inactive evidence and
    no PASS KD measurement. The latter remains explicitly review-only.
    """

    visualization = config.get("structure_visualization", {})
    cutoff = float(visualization.get("candidate_cluster_ca_cutoff_A", 8.0))
    loss_threshold = float(visualization.get("moderate_loss_max_KD_fold", 2.0))
    mapping_by_position = {
        int(row["sequence_position_1based"]): row for row in mapping_rows
    }
    candidates: list[dict[str, Any]] = []
    for position_record in position_summary:
        position = int(position_record["sequence_position_1based"])
        structure = mapping_by_position.get(position)
        if structure is None:
            continue
        coordinates = [
            to_float(structure.get("ca_x_A")),
            to_float(structure.get("ca_y_A")),
            to_float(structure.get("ca_z_A")),
        ]
        if any(value is None for value in coordinates):
            continue
        worst_fold = to_float(position_record.get("worst_KD_fold_vs_WT"))
        pass_count = int(position_record.get("variants_with_pass_measurement") or 0)
        inactive_count = int(position_record.get("variants_with_any_inactive") or 0)
        if worst_fold is not None and worst_fold >= loss_threshold:
            evidence_class = "QC_SUPPORTED_STRONG_KD_LOSS"
        elif pass_count == 0 and inactive_count > 0:
            evidence_class = "INACTIVE_WITHOUT_PASS_KD_REVIEW"
        else:
            continue
        candidates.append(
            {
                **position_record,
                "evidence_class": evidence_class,
                "ca_x_A": coordinates[0],
                "ca_y_A": coordinates[1],
                "ca_z_A": coordinates[2],
            }
        )

    adjacency: dict[int, set[int]] = {
        index: set() for index in range(len(candidates))
    }
    for left in range(len(candidates)):
        for right in range(left + 1, len(candidates)):
            distance = math.dist(
                [
                    float(candidates[left]["ca_x_A"]),
                    float(candidates[left]["ca_y_A"]),
                    float(candidates[left]["ca_z_A"]),
                ],
                [
                    float(candidates[right]["ca_x_A"]),
                    float(candidates[right]["ca_y_A"]),
                    float(candidates[right]["ca_z_A"]),
                ],
            )
            if distance <= cutoff:
                adjacency[left].add(right)
                adjacency[right].add(left)

    components: list[list[int]] = []
    unvisited = set(adjacency)
    while unvisited:
        seed = min(unvisited)
        stack = [seed]
        component: list[int] = []
        unvisited.remove(seed)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency[current]):
                if neighbor in unvisited:
                    unvisited.remove(neighbor)
                    stack.append(neighbor)
        components.append(
            sorted(
                component,
                key=lambda index: int(candidates[index]["sequence_position_1based"]),
            )
        )
    components.sort(
        key=lambda component: int(candidates[component[0]]["sequence_position_1based"])
    )

    member_rows: list[dict[str, Any]] = []
    cluster_summaries: list[dict[str, Any]] = []
    for cluster_index, component in enumerate(components, start=1):
        cluster_id = f"C{cluster_index}"
        members = [candidates[index] for index in component]
        labels = [
            f"{member['wt_aa']}{int(member['sequence_position_1based'])}"
            for member in members
        ]
        folds = [
            float(member["worst_KD_fold_vs_WT"])
            for member in members
            if to_float(member.get("worst_KD_fold_vs_WT")) is not None
        ]
        sasas = [
            float(member["residue_sasa_A2"])
            for member in members
            if to_float(member.get("residue_sasa_A2")) is not None
        ]
        plddts = [
            float(member["plddt_ca"])
            for member in members
            if to_float(member.get("plddt_ca")) is not None
        ]
        cluster_summaries.append(
            {
                "cluster_id": cluster_id,
                "member_positions": [
                    int(member["sequence_position_1based"]) for member in members
                ],
                "member_labels": labels,
                "member_count": len(members),
                "is_multi_position_cluster": len(members) >= 2,
                "qc_supported_strong_loss_positions": sum(
                    member["evidence_class"] == "QC_SUPPORTED_STRONG_KD_LOSS"
                    for member in members
                ),
                "inactive_without_pass_kd_review_positions": sum(
                    member["evidence_class"]
                    == "INACTIVE_WITHOUT_PASS_KD_REVIEW"
                    for member in members
                ),
                "ca_connectivity_cutoff_A": cutoff,
                "maximum_observed_KD_fold_vs_WT": max(folds) if folds else None,
                "median_residue_sasa_A2": statistics.median(sasas) if sasas else None,
                "mean_plddt_ca": statistics.fmean(plddts) if plddts else None,
            }
        )
        for member in members:
            neighbor_distances = []
            for other in members:
                if other is member:
                    continue
                neighbor_distances.append(
                    math.dist(
                        [
                            float(member["ca_x_A"]),
                            float(member["ca_y_A"]),
                            float(member["ca_z_A"]),
                        ],
                        [
                            float(other["ca_x_A"]),
                            float(other["ca_y_A"]),
                            float(other["ca_z_A"]),
                        ],
                    )
                )
            member_rows.append(
                {
                    "cluster_id": cluster_id,
                    "is_multi_position_cluster": len(members) >= 2,
                    "sequence_position_1based": int(
                        member["sequence_position_1based"]
                    ),
                    "wt_aa": member["wt_aa"],
                    "evidence_class": member["evidence_class"],
                    "worst_KD_mutation": member.get("worst_KD_mutation"),
                    "worst_KD_fold_vs_WT": member.get("worst_KD_fold_vs_WT"),
                    "variants_with_any_inactive": member.get(
                        "variants_with_any_inactive"
                    ),
                    "plddt_ca": member.get("plddt_ca"),
                    "residue_sasa_A2": member.get("residue_sasa_A2"),
                    "sidechain_sasa_A2": member.get("sidechain_sasa_A2"),
                    "nearest_cluster_member_ca_distance_A": (
                        min(neighbor_distances) if neighbor_distances else None
                    ),
                    "ca_x_A": member["ca_x_A"],
                    "ca_y_A": member["ca_y_A"],
                    "ca_z_A": member["ca_z_A"],
                }
            )
    return member_rows, cluster_summaries


def structure_fields(kind: str, ddg: float | None = None) -> dict[str, Any]:
    raw = to_float(ddg) if kind == "QUANTITATIVE" else None
    if raw is None and kind == "QUANTITATIVE":
        kind = "NO_PASS_DATA"
    return dict(visualization_class=kind, display_ddG_kcal_mol=raw,
                normalized_ddG_for_color=None, visualization_color_hex=STRUCTURE_COLORS.get(kind))


def structure_scale(positions, variants=()):
    """Endpoints of the residue maps, with a fallback for variant-only signs."""
    selected = [r["display_ddG_kcal_mol"] for r in positions if r["display_ddG_kcal_mol"] is not None]
    other = [r["display_ddG_kcal_mol"] for r in variants if r["display_ddG_kcal_mol"] is not None]
    return {"negative_endpoint_kcal_mol": min((v for v in selected if v < 0),
                                               default=min((v for v in other if v < 0), default=None)),
            "positive_endpoint_kcal_mol": max((v for v in selected if v > 0),
                                               default=max((v for v in other if v > 0), default=None))}


def apply_structure_scale(rows, scale):
    for row in rows:
        raw = row["display_ddG_kcal_mol"]
        if raw is None:
            continue
        if raw < 0:
            value = max(-1.0, raw / abs(scale["negative_endpoint_kcal_mol"]))
        elif raw > 0:
            value = min(1.0, raw / scale["positive_endpoint_kcal_mol"])
        else:
            value = 0.0
        faint = int(255 * (1 - abs(value)) + 0.5)
        rgb = (faint, faint, 255) if value < 0 else (255, faint, faint)
        row["normalized_ddG_for_color"] = value
        row["visualization_color_hex"] = "#{:02X}{:02X}{:02X}".format(*rgb)


def substitution_display(effect: dict[str, Any]) -> dict[str, Any]:
    phenotype = effect.get("phenotype_consensus")
    status = effect.get("variant_qc_status")
    energy = to_float(effect.get("delta_delta_G_kcal_mol_pass"))
    if phenotype == "consistently_inactive":
        return structure_fields("INACTIVE_ONLY")
    if (phenotype == "mixed_active_inactive" or effect.get("replicate_consistency_qc_status") == "REVIEW"
            or (status == "REVIEW" and (energy is not None or int(effect.get("n_review_fits") or 0)))):
        return structure_fields("REVIEW")
    if status == "PASS":
        return structure_fields("QUANTITATIVE", energy)
    return structure_fields("NO_PASS_DATA")


def build_structure_visualization_records(
    mapped_effects: list[dict[str, Any]], config: dict[str, Any]
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for effect in mapped_effects:
        if effect.get("effect_attribution") != "direct_single_substitution" or not int(effect.get("n_genotype_pass") or 0):
            continue
        kd_fold = to_float(effect.get("KD_fold_vs_WT_pass"))
        output.append(
            {
                "mutation_label": effect.get("mutation_label", ""),
                "observed_rois": effect.get("observed_rois", ""),
                "median_rois": effect.get("median_rois", ""),
                "sequence_position_1based": effect.get("sequence_position_1based"),
                "wt_aa": effect.get("wt_aa", ""),
                "mutant_aa": effect.get("mutant_aa", ""),
                "structure_chain_id": effect.get("structure_chain_id", ""),
                "structure_auth_chain_id": effect.get("structure_auth_chain_id", ""),
                "structure_insertion_code": effect.get("structure_insertion_code", ""),
                "structure_residue_number": effect.get(
                    "structure_residue_number", ""
                ),
                "alignment_status": effect.get("alignment_status", ""),
                "KD_fold_vs_WT_pass": kd_fold,
                "log2_KD_fold_vs_WT_pass": (
                    math.log2(kd_fold) if kd_fold is not None and kd_fold > 0 else None
                ),
                "kon_fold_vs_WT_pass": effect.get("kon_fold_vs_WT_pass"),
                "koff_fold_vs_WT_pass": effect.get("koff_fold_vs_WT_pass"),
                "median_KD_M_pass": effect.get("median_KD_M_pass"),
                "n_active_quantitative_pass": effect.get(
                    "n_active_quantitative_pass", 0
                ),
                "n_inactive": effect.get("n_inactive", 0),
                "phenotype_consensus": effect.get("phenotype_consensus", ""),
                "variant_qc_status": effect.get("variant_qc_status", ""),
                **substitution_display(effect),
                "plddt_ca": effect.get("plddt_ca"),
                "residue_sasa_A2": effect.get("residue_sasa_A2"),
            }
        )
    return sorted(
        output,
        key=lambda row: (
            int(row["sequence_position_1based"] or 0),
            str(row["mutation_label"]),
        ),
    )


def build_structure_visualization_positions(
    variant_records: list[dict[str, Any]], config: dict[str, Any]
) -> list[dict[str, Any]]:
    groups: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for record in variant_records:
        groups[int(record["sequence_position_1based"])].append(record)
    output: list[dict[str, Any]] = []
    for position, records in sorted(groups.items()):
        usable = [r for r in records if r["visualization_class"] == "QUANTITATIVE"]
        selected = min(usable, key=lambda r: (-abs(r["display_ddG_kcal_mol"]), r["mutation_label"]), default=None)
        statuses = {r["visualization_class"] for r in records}
        if selected is not None:
            fields = structure_fields("QUANTITATIVE", selected["display_ddG_kcal_mol"])
        elif "REVIEW" in statuses:
            fields = structure_fields("REVIEW")
        elif statuses == {"INACTIVE_ONLY"}:
            fields = structure_fields("INACTIVE_ONLY")
        else:
            fields = structure_fields("NO_PASS_DATA")
        measured = [
            record
            for record in records
            if to_float(record.get("KD_fold_vs_WT_pass")) is not None
        ]
        best = (
            min(measured, key=lambda row: float(row["KD_fold_vs_WT_pass"]))
            if measured
            else None
        )
        worst = (
            max(measured, key=lambda row: float(row["KD_fold_vs_WT_pass"]))
            if measured
            else None
        )
        first = records[0]
        output.append(
            {
                "sequence_position_1based": position,
                "wt_aa": first.get("wt_aa", ""),
                "structure_chain_id": first.get("structure_chain_id", ""),
                "structure_auth_chain_id": first.get("structure_auth_chain_id", ""),
                "structure_insertion_code": first.get("structure_insertion_code", ""),
                "structure_residue_number": first.get(
                    "structure_residue_number", ""
                ),
                "single_substitution_variants": len(records),
                "pass_variants": len(usable),
                "review_variants": sum(
                    record.get("variant_qc_status") != "PASS" for record in records
                ),
                "position_visualization_class": fields["visualization_class"],
                "selected_substitution": selected["mutation_label"] if selected is not None else "",
                **fields,
                "best_KD_mutation": best.get("mutation_label") if best else None,
                "best_KD_fold_vs_WT": (
                    best.get("KD_fold_vs_WT_pass") if best else None
                ),
                "worst_KD_mutation": worst.get("mutation_label") if worst else None,
                "worst_KD_fold_vs_WT": (
                    worst.get("KD_fold_vs_WT_pass") if worst else None
                ),
                "has_mixed_active_inactive_variant": any(
                    record.get("phenotype_consensus") == "mixed_active_inactive"
                    for record in records
                ),
                "plddt_ca": first.get("plddt_ca"),
                "residue_sasa_A2": first.get("residue_sasa_A2"),
            }
        )
    return output


def include_confounded_positions(positions, mapped_effects):
    present = {int(r["sequence_position_1based"]) for r in positions}
    for r in mapped_effects:
        pos = int(r["sequence_position_1based"])
        if pos in present or r.get("alignment_status") != "exact":
            continue
        if r.get("effect_attribution") == "multi_substitution_confounded" and int(r.get("n_genotype_pass") or 0):
            positions.append(dict(sequence_position_1based=pos, wt_aa=r["wt_aa"], structure_chain_id=r["structure_chain_id"],
                                  structure_residue_number=r["structure_residue_number"], selected_substitution="",
                                  structure_auth_chain_id=r.get("structure_auth_chain_id",""),
                                  structure_insertion_code=r.get("structure_insertion_code",""),
                                  position_visualization_class="NO_PASS_DATA", attribution="multi_substitution_only",
                                  **structure_fields("NO_PASS_DATA")))
            present.add(pos)
    positions.sort(key=lambda r: int(r["sequence_position_1based"]))


def all_residue_colors(mapping, tested):
    lookup = {int(row["sequence_position_1based"]): row for row in tested}
    rows = []
    for mapped in mapping:
        pos = int(mapped["sequence_position_1based"])
        default = dict(position_visualization_class="NO_DATA", selected_substitution="",
                       **structure_fields("NO_DATA"))
        result = lookup.get(pos, default) if mapped.get("alignment_status") == "exact" else default
        rows.append({**mapped, **result})
    return rows


def report_residue_classes(rows: list[dict[str, Any]] | None) -> str:
    if not rows:
        return "Residue classifications are unavailable. Supply a structure and enable structure visualization to list residues."
    groups = {key: [] for key in STRUCTURE_COLORS}
    for row in sorted(rows, key=lambda r: int(r["sequence_position_1based"])):
        kind = row["position_visualization_class"]
        if kind in groups:
            groups[kind].append(f"{row['sequence_aa']}{row['sequence_position_1based']}")
    chain = next((r["structure_chain_id"] for r in rows if r.get("structure_chain_id")), "")
    chain_label = f" (chain {chain})" if chain else ""
    lines = [f"Residues use the input WT sequence numbered from 1{chain_label}. Each position has one structure color.",
             "", "| Structure class | Count | WT residues |", "|---|---:|---|"]
    for kind, name in STRUCTURE_CLASS_LABELS.items():
        labels = groups[kind]
        lines.append(f"| {name} | {len(labels)} | {', '.join(labels) or 'None'} |")
    if any(r.get("alignment_status") != "exact" for r in rows):
        lines.extend(["", "No data also includes WT positions without an exact structure sequence match."])
    return "\n".join(lines)


def render_dependency_free_structure_svg(
    mapping_rows: list[dict[str, Any]],
    position_records: list[dict[str, Any]],
    config: dict[str, Any],
    scale: dict[str, Any],
) -> str:
    coordinate_rows = [
        row
        for row in mapping_rows
        if all(to_float(row.get(axis)) is not None for axis in ("ca_x_A", "ca_y_A", "ca_z_A"))
    ]
    if not coordinate_rows:
        raise ValueError("No C-alpha coordinates are available for SVG rendering")

    tested_by_position = {
        int(record["sequence_position_1based"]): record
        for record in position_records
    }
    label_candidates = sorted(
        position_records,
        key=lambda row: (
            bool(row.get("has_mixed_active_inactive_variant")),
            row.get("position_visualization_class") == "INACTIVE_ONLY",
            to_float(row.get("worst_KD_fold_vs_WT")) or -1.0,
        ),
        reverse=True,
    )[:12]
    labeled_positions = {
        int(record["sequence_position_1based"]) for record in label_candidates
    }

    width = 1260
    height = 525
    panel_width = 390
    panel_height = 360
    panel_top = 70
    panel_margin = 24
    projections = [
        ("Front (X/Y)", "ca_x_A", "ca_y_A"),
        ("Top (X/Z)", "ca_x_A", "ca_z_A"),
        ("Side (Y/Z)", "ca_y_A", "ca_z_A"),
    ]

    def project(
        row: dict[str, Any],
        x_axis: str,
        y_axis: str,
        panel_left: float,
    ) -> tuple[float, float]:
        x_values = [float(item[x_axis]) for item in coordinate_rows]
        y_values = [float(item[y_axis]) for item in coordinate_rows]
        x_min, x_max = min(x_values), max(x_values)
        y_min, y_max = min(y_values), max(y_values)
        x_span = max(x_max - x_min, 1.0)
        y_span = max(y_max - y_min, 1.0)
        scale = min(
            (panel_width - 2 * panel_margin) / x_span,
            (panel_height - 2 * panel_margin) / y_span,
        )
        x_padding = (panel_width - x_span * scale) / 2
        y_padding = (panel_height - y_span * scale) / 2
        x = panel_left + x_padding + (float(row[x_axis]) - x_min) * scale
        y = panel_top + panel_height - y_padding - (float(row[y_axis]) - y_min) * scale
        return x, y

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<defs><linearGradient id="ddg"><stop offset="0%" stop-color="#0000FF"/><stop offset="50%" stop-color="#FFFFFF"/><stop offset="100%" stop-color="#FF0000"/></linearGradient></defs>',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<text x="30" y="32" font-family="Arial, sans-serif" font-size="20" font-weight="bold" fill="#1F2937">Binding ΔΔG on the predicted structure</text>',
        '<text x="30" y="52" font-family="Arial, sans-serif" font-size="12" fill="#4B5563">The PASS substitution with the largest absolute ΔΔG sets each residue color.</text>',
    ]
    sorted_rows = sorted(coordinate_rows, key=lambda row: int(row["sequence_position_1based"]))
    for panel_index, (title, x_axis, y_axis) in enumerate(projections):
        panel_left = 20 + panel_index * 410
        lines.append(
            f'<rect x="{panel_left}" y="{panel_top}" width="{panel_width}" height="{panel_height}" rx="8" fill="#FAFAFA" stroke="#D1D5DB"/>'
        )
        lines.append(
            f'<text x="{panel_left + 12}" y="{panel_top + 22}" font-family="Arial, sans-serif" font-size="13" font-weight="bold" fill="#374151">{escape(title)}</text>'
        )
        projected = [project(row, x_axis, y_axis, panel_left) for row in sorted_rows]
        point_text = " ".join(f"{x:.2f},{y:.2f}" for x, y in projected)
        lines.append(
            f'<polyline points="{point_text}" fill="none" stroke="#808080" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>'
        )
        for row, (x, y) in zip(sorted_rows, projected):
            position = int(row["sequence_position_1based"])
            record = tested_by_position.get(position)
            radius = 5.5 if record else 1.5
            color = record["visualization_color_hex"] if record else STRUCTURE_COLORS["NO_DATA"]
            aa = str(row.get("sequence_aa", ""))
            title = f"{aa}{position} · " + (record["position_visualization_class"] if record else "NO_DATA")
            if record and record.get("display_ddG_kcal_mol") is not None:
                title += f" · {record['selected_substitution']} · {record['display_ddG_kcal_mol']:.4g} kcal/mol"
            stroke = "#27323A" if record else color
            lines.append(f'<circle data-position="{position}" cx="{x:.2f}" cy="{y:.2f}" r="{radius}" '
                         f'fill="{color}" stroke="{stroke}" stroke-width="0.8"><title>{escape(title)}</title></circle>')
            if panel_index == 0 and position in labeled_positions:
                label = f"{record.get('wt_aa', '')}{position}"
                lines.append(
                    f'<text x="{x + 7:.2f}" y="{y - 7:.2f}" font-family="Arial, sans-serif" font-size="10" font-weight="bold" fill="#111827">{escape(label)}</text>'
                )

    def endpoint(key):
        energy = scale[key + "_endpoint_kcal_mol"]
        return f"{energy:+.3g} kcal/mol" if energy is not None else "no observations"
    lines.extend([
        '<text x="30" y="461" font-family="Arial, sans-serif" font-size="14">Normalized binding ΔΔG (color only)</text>',
        '<rect x="30" y="471" width="330" height="12" fill="url(#ddg)" stroke="#67737C"/>',
        '<g font-family="Arial, sans-serif" font-size="12" fill="#263640">'
        '<text x="30" y="501">−1</text><text x="195" y="501" text-anchor="middle">0 (WT)</text>'
        '<text x="360" y="501" text-anchor="end">+1</text></g>',
    ])
    lines.append(f'<text x="30" y="519" font-family="Arial, sans-serif" font-size="12">Blue: {endpoint("negative")} · Red: {endpoint("positive")}</text>')
    for i, (kind, label) in enumerate(STRUCTURE_CLASS_LABELS.items()):
        x = 410 + i * 207
        lines.append(f'<circle cx="{x}" cy="476" r="6" fill="{STRUCTURE_COLORS[kind]}" stroke="#475569"/>')
        lines.append(f'<text x="{x+13}" y="481" font-family="Arial, sans-serif" font-size="14">{label}</text>')
    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def safe_viewer_identifier(value: Any, fallback: str) -> str:
    text = str(value or "")
    return text if re.fullmatch(r"[A-Za-z0-9_+.-]+", text) else fallback


def pymol_color_name(classification: str) -> str:
    return "ksm_" + re.sub(r"[^a-z0-9]+", "_", classification.lower()).strip("_")


def hex_to_rgb_fraction(color: str) -> tuple[float, float, float]:
    value = color.lstrip("#")
    if not re.fullmatch(r"[0-9A-Fa-f]{6}", value):
        value = "BDBDBD"
    return tuple(int(value[index : index + 2], 16) / 255 for index in (0, 2, 4))


def pymol_residue_selection(record: dict[str, Any], object_name: str) -> str:
    chain = safe_viewer_identifier(record.get("structure_auth_chain_id") or record.get("structure_chain_id"), "")
    residue = safe_viewer_identifier(
        record.get("structure_residue_number"),
        str(record.get("sequence_position_1based", "")),
    )
    residue += safe_viewer_identifier(record.get("structure_insertion_code"), "")
    chain_clause = f" and chain {chain}" if chain else ""
    return f"{object_name}{chain_clause} and resi {residue}"


def render_pymol_overview_script(
    position_records: list[dict[str, Any]], config: dict[str, Any]
) -> str:
    colors = {r["visualization_color_hex"] for r in position_records} | {STRUCTURE_COLORS["NO_DATA"]}
    lines = [
        "# Generated by MAXTIA kinetic-to-structure mapping.",
        "# Run from the result directory: pymol visualize_overview.pml",
        "reinitialize",
        "load structure_model.cif, kinetic_structure",
        "hide everything, kinetic_structure",
        "show cartoon, kinetic_structure",
        "color 0x808080, kinetic_structure",
        "set cartoon_transparency, 0, kinetic_structure",
        "set stick_radius, 0.18",
    ]
    for color in sorted(colors):
        red, green, blue = hex_to_rgb_fraction(str(color))
        lines.append(
            f"set_color {pymol_color_name(color)}, "
            f"[{red:.4f}, {green:.4f}, {blue:.4f}]"
        )
    for record in position_records:
        selection_name = f"position_{int(record['sequence_position_1based'])}"
        expression = pymol_residue_selection(record, "kinetic_structure")
        color = record["visualization_color_hex"]
        lines.extend(
            [
                f"select {selection_name}, {expression}",
                f"color {pymol_color_name(color)}, {selection_name}",
                f"show sticks, {selection_name}",
            ]
        )
    lines.extend(
        [
            "select tested_positions, "
            + (" or ".join(
                f"position_{int(record['sequence_position_1based'])}"
                for record in position_records
            ) or "none"),
            "orient kinetic_structure",
            "zoom kinetic_structure, 4",
            "scene overview, store",
            "deselect",
        ]
    )
    return "\n".join(lines) + "\n"


def render_pymol_variant_scenes_script(
    variant_records: list[dict[str, Any]], config: dict[str, Any]
) -> str:
    colors = {r["visualization_color_hex"] for r in variant_records} | {STRUCTURE_COLORS["NO_DATA"]}
    lines = [
        "# Generated by MAXTIA kinetic-to-structure mapping.",
        "# Use the PyMOL Scene menu or PageDown/PageUp to inspect variants.",
        "reinitialize",
        "load structure_model.cif, kinetic_structure",
        "hide everything, kinetic_structure",
        "show cartoon, kinetic_structure",
        "set cartoon_transparency, 0, kinetic_structure",
        "set stick_radius, 0.20",
    ]
    for color in sorted(colors):
        red, green, blue = hex_to_rgb_fraction(str(color))
        lines.append(
            f"set_color {pymol_color_name(color)}, "
            f"[{red:.4f}, {green:.4f}, {blue:.4f}]"
        )
    scene_definitions: list[dict[str, str]] = []
    for index, record in enumerate(variant_records, start=1):
        mutation = safe_viewer_identifier(record.get("mutation_label"), f"variant_{index}")
        selection_name = "variant_" + re.sub(r"[^A-Za-z0-9_]+", "_", mutation)
        selection = pymol_residue_selection(record, "kinetic_structure")
        classification = str(record["visualization_class"])
        kd_fold = to_float(record.get("KD_fold_vs_WT_pass"))
        fold_label = f"{kd_fold:.3g}x KD" if kd_fold is not None else classification
        label_text = f"{mutation} {fold_label}".replace('"', "'")
        lines.extend(
            [
                f"select {selection_name}, {selection}",
                f"label {selection_name} and name CA, \"{label_text}\"",
            ]
        )
        scene_definitions.append(
            {
                "mutation": mutation,
                "selection_name": selection_name,
                "color_name": pymol_color_name(record["visualization_color_hex"]),
            }
        )
    lines.extend(["hide sticks, kinetic_structure", "hide labels, kinetic_structure"])
    for scene in scene_definitions:
        selection_name = scene["selection_name"]
        lines.extend(
            [
                "hide sticks, kinetic_structure",
                "hide labels, kinetic_structure",
                "color 0x808080, kinetic_structure",
                f"color {scene['color_name']}, {selection_name}",
                f"show sticks, {selection_name}",
                f"show labels, {selection_name} and name CA",
                f"orient {selection_name}",
                f"zoom {selection_name}, 10",
                f"scene {scene['mutation']}, store",
            ]
        )
    if scene_definitions:
        first_scene = safe_viewer_identifier(
            variant_records[0].get("mutation_label"), "variant_1"
        )
        lines.append(f"scene {first_scene}, recall")
    lines.append("deselect")
    return "\n".join(lines) + "\n"


def render_chimerax_overview_script(
    position_records: list[dict[str, Any]], config: dict[str, Any]
) -> str:
    lines = [
        "# Generated by MAXTIA kinetic-to-structure mapping.",
        "# Run in a fresh ChimeraX session, from the result folder.",
        "open structure_model.cif id #1",
        "hide #1 atoms",
        "cartoon #1",
        f"color #1 {STRUCTURE_COLORS['NO_DATA']} target acs",
    ]
    for record in position_records:
        chain = safe_viewer_identifier(record.get("structure_auth_chain_id") or record.get("structure_chain_id"), "")
        residue = safe_viewer_identifier(
            record.get("structure_residue_number"),
            str(record.get("sequence_position_1based", "")),
        )
        residue += safe_viewer_identifier(record.get("structure_insertion_code"), "")
        target = f"#1/{chain}:{residue}" if chain else f"#1:{residue}"
        lines.append(f"color {target} {record['visualization_color_hex']} target acs")
        lines.append(f"show {target} atoms")
    lines.extend(["view #1", "lighting soft"])
    return "\n".join(lines) + "\n"


def render_visualization_readme(config: dict[str, Any]) -> str:
    return """# Structure visualization outputs

The visualization files are generated deterministically from QC-filtered single-substitution records. Multi-substitution variants are retained in the analytical CSV files but are not assigned as single-residue effects.

## Files

- `structure_model.cif`: local copy of the analyzed structure
- `structure_visualization_variants.csv`: one row per single substitution, with KD/kon/koff folds, QC status, residue mapping, and display color
- `structure_visualization_positions.csv`: position-level overview that preserves substitution-dependent directionality
- `structure_candidate_cluster_members.csv`: structure-connected strong-loss and review-only Inactive candidates
- `visualize_overview.pml`: PyMOL position overview
- `visualize_variant_scenes.pml`: one PyMOL scene per single substitution
- `visualize_overview.cxc`: ChimeraX position overview
- `structure_overview.svg`: dependency-free Cα projections that open in a web browser without PyMOL

## PyMOL

Open a terminal in this result directory and run:

```bash
pymol visualize_overview.pml
pymol visualize_variant_scenes.pml
```

In the variant-scene view, use the Scene menu or PageDown/PageUp to move between substitutions.

## ChimeraX

Open a terminal in this result directory and run:

```bash
chimerax visualize_overview.cxc
```

## Binding ΔΔG colors

The binding energy change is RT ln(KD variant / KD WT), in kcal/mol. The most extreme absolute PASS energy among single substitutions sets the residue color. Equal magnitudes use ascending substitution label. The signed, unclipped energy and selected substitution are recorded in residue_colors.csv. Only the display is normalized. Negative values are divided by the magnitude of the most negative selected residue energy, and positive values by the largest positive selected energy. Blue is −1, white is zero, and red is +1 on this dimensionless scale. The endpoints in structure_color_scale.json are in kcal/mol. All batch and condition maps in a study use the same endpoints. If a sign occurs only in individual substitution scenes, its variant extreme supplies that endpoint. Individual scenes saturate at the residue endpoints.

Gray indicates no attributable experimental observation. Purple indicates review, yellow indicates exclusively inactive evidence, and magenta indicates observations without a usable PASS estimate. QC colors apply when no eligible single substitution supplies an energy. A multiple substitution alone cannot supply an independent residue energy.
"""


def render_static_images_with_pymol(
    out_dir: Path,
    variant_records: list[dict[str, Any]],
    config: dict[str, Any],
) -> dict[str, Any]:
    visualization = config.get("structure_visualization", {})
    if not visualization.get("render_static_pngs_if_pymol_available", False):
        return {"status": "DISABLED", "files": []}
    pymol = shutil.which("pymol")
    if pymol is None:
        return {"status": "SKIPPED_PYMOL_NOT_FOUND", "files": []}
    width = int(visualization.get("static_png_width", 1400))
    height = int(visualization.get("static_png_height", 1000))
    rendered: list[str] = []
    overview_name = "structure_effect_overview.png"
    overview_command = [
        pymol,
        "-cq",
        "visualize_overview.pml",
        "-d",
        "bg_color white",
        "-d",
        "set ray_opaque_background, off",
        "-d",
        f"png {overview_name}, width={width}, height={height}, ray=1",
        "-d",
        "quit",
    ]
    try:
        subprocess.run(
            overview_command,
            cwd=out_dir,
            check=True,
            capture_output=True,
            text=True,
            timeout=180,
        )
        rendered.append(overview_name)
        measured_pass = [
            record
            for record in variant_records
            if record.get("variant_qc_status") == "PASS"
            and to_float(record.get("KD_fold_vs_WT_pass")) is not None
        ]
        scene_records: list[dict[str, Any]] = []
        if measured_pass:
            scene_records.append(
                min(
                    measured_pass,
                    key=lambda row: float(row["KD_fold_vs_WT_pass"]),
                )
            )
            scene_records.append(
                max(
                    measured_pass,
                    key=lambda row: float(row["KD_fold_vs_WT_pass"]),
                )
            )
        review_records = [
            record
            for record in variant_records
            if record.get("visualization_class") == "REVIEW"
        ]
        if review_records:
            scene_records.append(review_records[0])
        unique_scenes: list[dict[str, Any]] = []
        seen: set[str] = set()
        for record in scene_records:
            mutation = str(record.get("mutation_label", ""))
            if mutation and mutation not in seen:
                unique_scenes.append(record)
                seen.add(mutation)
        if unique_scenes:
            variant_command = [pymol, "-cq", "visualize_variant_scenes.pml"]
            for record in unique_scenes:
                mutation = safe_viewer_identifier(record["mutation_label"], "variant")
                selection_name = "variant_" + re.sub(
                    r"[^A-Za-z0-9_]+", "_", mutation
                )
                classification = str(record["visualization_class"])
                kd_fold = to_float(record.get("KD_fold_vs_WT_pass"))
                fold_label = (
                    f"{kd_fold:.3g}x KD" if kd_fold is not None else classification
                )
                label_text = f"{mutation} {fold_label}".replace('"', "'")
                filename = f"variant_{mutation}.png"
                variant_command.extend(
                    [
                        "-d",
                        f"scene {mutation}, recall",
                        "-d",
                        "hide labels, kinetic_structure",
                        "-d",
                        f"label {selection_name} and name CA, \"{label_text}\"",
                        "-d",
                        f"show labels, {selection_name} and name CA",
                        "-d",
                        "bg_color white",
                        "-d",
                        f"png {filename}, width={width}, height={height}, ray=1",
                    ]
                )
                rendered.append(filename)
            variant_command.extend(["-d", "quit"])
            subprocess.run(
                variant_command,
                cwd=out_dir,
                check=True,
                capture_output=True,
                text=True,
                timeout=180,
            )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
        return {
            "status": "FAILED",
            "files": [name for name in rendered if (out_dir / name).exists()],
            "error_type": type(error).__name__,
        }
    return {"status": "RENDERED", "files": rendered, "pymol": Path(pymol).name}


def to_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def parse_lsa(path: Path, config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    rows = read_xlsx_sheet_rows(path, config["lsa_sheet"])
    required = {"ROI", "Ligand Name", "ka (1/Ms)", "kd (1/s)", "KD (M)"}
    header_row_number: int | None = None
    headers_by_col: dict[int, str] = {}
    for row_number, cells in rows:
        candidate = {column: normalize_header(value) for column, value in cells.items()}
        if required.issubset(set(candidate.values())):
            header_row_number = row_number
            headers_by_col = candidate
            break
    if header_row_number is None:
        raise ValueError(f"Could not locate the kinetic header in {config['lsa_sheet']!r}")

    warnings: list[str] = []
    records: list[dict[str, Any]] = []
    pattern = re.compile(config["kinetic_ligand_pattern"])
    for row_number, cells in rows:
        if row_number <= header_row_number:
            continue
        raw = {
            header: cells.get(column)
            for column, header in headers_by_col.items()
            if header
        }
        if not any(value not in (None, "") for value in raw.values()):
            continue
        record = {output: raw.get(source) for source, output in LSA_FIELDS.items()}
        record["source_lsa_row"] = row_number
        for field in NUMERIC_LSA_FIELDS:
            record[field] = to_float(record.get(field))
        # Backward-compatible aliases. Carterra calls these columns "Error";
        # the manual does not identify them as standard errors.
        record["kon_se_M_inv_s"] = record.get("kon_fit_error_M_inv_s")
        record["koff_se_s_inv"] = record.get("koff_fit_error_s_inv")
        record["KD_se_M"] = record.get("KD_fit_error_M")
        record["Rmax_se_RU"] = record.get("Rmax_fit_error_RU")
        record["roi"] = "" if record["roi"] is None else str(record["roi"]).strip()
        record["ligand_name"] = str(record.get("ligand_name") or "").strip()
        record["lsa_flag"] = str(record.get("lsa_flag") or "").strip()
        match = pattern.fullmatch(record["ligand_name"])
        if match:
            plate = int(match.group("plate"))
            well = int(match.group("well"))
            record["normalized_well"] = (
                f"{config['ngs_well_plate_prefix']}"
                f"{plate:0{int(config['ngs_plate_width'])}d}-"
                f"{well:0{int(config['ngs_well_width'])}d}"
            )
        else:
            record["normalized_well"] = ""
            warnings.append(
                f"LSA row {row_number}: ligand key {record['ligand_name']!r} did not match the configured pattern"
            )
        records.append(record)
    return records, warnings


def parse_mutation_labels(raw: str) -> tuple[list[str], str | None]:
    if raw.strip() == "":
        return [], None
    try:
        value = ast.literal_eval(raw)
    except (SyntaxError, ValueError) as exc:
        return [], f"Could not parse mutation list: {exc}"
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        return [], "Mutation field is not a list of strings"
    return [item.strip().upper() for item in value], None


def parse_ngs_plate_block(value: Any) -> tuple[int, int] | None:
    if value in (None, ""):
        return None
    match = re.fullmatch(
        r"\s*P?(\d+)\s*-\s*P?(\d+)\s*", str(value), re.IGNORECASE
    )
    if not match:
        raise ValueError(
            "ngs_plate_block must use a range such as 'P01-P04', 'P05-P08', or 'P09-P12'"
        )
    start, end = (int(match.group(1)), int(match.group(2)))
    if start < 1 or end < start:
        raise ValueError("ngs_plate_block must be an increasing positive plate range")
    return start, end


def parse_ngs_well(value: str, config: dict[str, Any]) -> tuple[int, int] | None:
    prefix = re.escape(str(config["ngs_well_plate_prefix"]))
    match = re.fullmatch(rf"{prefix}(\d+)-(\d+)", value, re.IGNORECASE)
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


def parse_ngs(path: Path, config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str]]:
    column_names = list(config["ngs_columns"])
    plate_block = parse_ngs_plate_block(config.get("ngs_plate_block"))
    target_plate_start = int(config.get("ngs_target_plate_start", 1))
    if target_plate_start < 1:
        raise ValueError("ngs_target_plate_start must be a positive integer")
    records: list[dict[str, Any]] = []
    warnings: list[str] = []
    selected_source_plates: set[int] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        for source_row, values in enumerate(reader, start=1):
            if source_row == 1 and config.get("ngs_has_header"):
                continue
            if len(values) != len(column_names):
                warnings.append(
                    f"NGS row {source_row}: expected {len(column_names)} columns, found {len(values)}"
                )
                values = (values + [""] * len(column_names))[: len(column_names)]
            record = dict(zip(column_names, values))
            record["source_ngs_row"] = source_row
            source_well = str(record.get("well") or "").strip().upper()
            parsed_well = parse_ngs_well(source_well, config)
            if parsed_well is None:
                warnings.append(
                    f"NGS row {source_row}: well key {source_well!r} did not match the configured format"
                )
                if plate_block is not None:
                    continue
                mapped_well = source_well
                source_plate_number = None
                mapped_plate_number = None
            else:
                source_plate_number, well_number = parsed_well
                if plate_block is not None:
                    source_start, source_end = plate_block
                    if not source_start <= source_plate_number <= source_end:
                        continue
                    mapped_plate_number = (
                        target_plate_start + source_plate_number - source_start
                    )
                else:
                    mapped_plate_number = source_plate_number
                selected_source_plates.add(source_plate_number)
                mapped_well = (
                    f"{config['ngs_well_plate_prefix']}"
                    f"{mapped_plate_number:0{int(config['ngs_plate_width'])}d}-"
                    f"{well_number:0{int(config['ngs_well_width'])}d}"
                )
            record["source_well"] = source_well
            record["well"] = mapped_well
            record["source_plate_number"] = source_plate_number
            record["mapped_plate_number"] = mapped_plate_number
            record["call_code"] = str(record.get("call_code") or "").strip().upper()
            total_valid_read_count = to_float(
                record.get("total_valid_read_count", record.get("read_count"))
            )
            record["total_valid_read_count"] = total_valid_read_count
            # Deprecated alias retained so older downstream readers do not break.
            record["read_count"] = total_valid_read_count
            inter_cluster_score = to_float(
                record.get("inter_cluster_mixing_score", record.get("metric_1"))
            )
            intra_cluster_score = to_float(
                record.get("intra_cluster_mixing_score", record.get("metric_2"))
            )
            record["inter_cluster_mixing_score"] = inter_cluster_score
            record["intra_cluster_mixing_score"] = intra_cluster_score
            # Deprecated aliases retained so older downstream readers do not break.
            record["metric_1"] = inter_cluster_score
            record["metric_2"] = intra_cluster_score
            record["dna_sequence"] = re.sub(
                r"\s+", "", str(record.get("dna_sequence") or "").upper()
            )
            record["aa_sequence"] = re.sub(
                r"\s+", "", str(record.get("aa_sequence") or "").upper()
            )
            labels, error = parse_mutation_labels(
                str(record.get("mutation_labels_raw") or "")
            )
            record["mutation_labels"] = labels
            record["mutation_parse_error"] = error or ""
            records.append(record)
    if plate_block is not None:
        expected_plates = set(range(plate_block[0], plate_block[1] + 1))
        if not records:
            raise ValueError(
                f"No NGS rows were found in configured plate block P{plate_block[0]:02d}-P{plate_block[1]:02d}"
            )
        if selected_source_plates != expected_plates:
            missing = sorted(expected_plates - selected_source_plates)
            raise ValueError(
                "Configured NGS plate block is incomplete; missing source plates: "
                + ", ".join(f"P{plate:02d}" for plate in missing)
            )
    return records, warnings


def validate_ngs_plate_scope(
    lsa_records: list[dict[str, Any]],
    ngs_records: list[dict[str, Any]],
    config: dict[str, Any],
) -> None:
    if parse_ngs_plate_block(config.get("ngs_plate_block")) is not None:
        return
    if not config.get("require_explicit_ngs_plate_block_for_superset", True):
        return
    lsa_plates = {
        parsed[0]
        for record in lsa_records
        if (parsed := parse_ngs_well(record.get("normalized_well", ""), config))
    }
    ngs_source_plates = {
        int(record["source_plate_number"])
        for record in ngs_records
        if record.get("source_plate_number") is not None
    }
    if lsa_plates and len(ngs_source_plates) > len(lsa_plates):
        raise ValueError(
            "The NGS file contains more source plates than the LSA experiment. "
            "Set ngs_plate_block (for example P01-P04, P05-P08, or P09-P12) "
            "according to the Primer Barcode type before joining."
        )


def infer_wt_sequence(
    ngs_records: list[dict[str, Any]], configured_sequence: str
) -> tuple[str, dict[str, Any]]:
    configured = re.sub(r"\s+", "", configured_sequence.upper())
    candidates = [
        record["aa_sequence"]
        for record in ngs_records
        if record["mutation_labels"] == ["WT"] and record["aa_sequence"]
    ]
    counts = Counter(candidates)
    if configured:
        sequence = configured
        source = "config"
    elif counts:
        top = counts.most_common()
        if len(top) > 1 and top[0][1] == top[1][1]:
            raise ValueError("WT sequence inference is ambiguous; provide wt_sequence in config")
        sequence = top[0][0]
        source = "inferred_from_ngs_WT_labels"
    else:
        sequence = ""
        source = "unavailable"
    return sequence, {
        "source": source,
        "sequence": sequence,
        "length": len(sequence),
        "wt_labeled_rows_with_sequence": len(candidates),
        "distinct_wt_labeled_sequences": len(counts),
        "selected_sequence_support": counts.get(sequence, 0),
    }


def add_issue(
    fail_issues: list[str], review_issues: list[str], code: str, severity: str
) -> None:
    target = fail_issues if severity == "fail" else review_issues
    if code not in target:
        target.append(code)


def validate_sequence_and_mutations(
    wt_sequence: str,
    aa_sequence: str,
    labels: list[str],
    require_explicit_labels: bool,
) -> tuple[bool, str, list[str]]:
    issues: list[str] = []
    if not aa_sequence:
        return False, "missing_sequence", ["AA_SEQUENCE_MISSING"]
    if not AA_RE.fullmatch(aa_sequence):
        issues.append("AA_SEQUENCE_INVALID_CHARACTERS")
    if not wt_sequence:
        issues.append("WT_SEQUENCE_UNAVAILABLE")
        return False, "unknown", issues
    if len(aa_sequence) != len(wt_sequence):
        issues.append("SEQUENCE_LENGTH_MISMATCH")
    if labels == ["WT"]:
        category = "WT"
        if aa_sequence != wt_sequence:
            issues.append("WT_LABEL_SEQUENCE_MISMATCH")
    elif not labels:
        category = "unlabeled_or_indel"
        if require_explicit_labels:
            issues.append("MUTATION_LABEL_MISSING_OR_INDEL")
        if aa_sequence == wt_sequence:
            issues.append("WT_SEQUENCE_WITHOUT_WT_LABEL")
    else:
        category = "single_substitution" if len(labels) == 1 else "multi_substitution"
        expected = list(wt_sequence)
        seen_positions: set[int] = set()
        for label in labels:
            match = SUBSTITUTION_RE.fullmatch(label)
            if not match:
                issues.append("MUTATION_LABEL_NOT_SUBSTITUTION")
                continue
            source, position_text, target = match.groups()
            position = int(position_text)
            if position in seen_positions:
                issues.append("DUPLICATE_MUTATION_POSITION")
                continue
            seen_positions.add(position)
            index = position - 1
            if index < 0 or index >= len(wt_sequence):
                issues.append("MUTATION_POSITION_OUT_OF_RANGE")
                continue
            if wt_sequence[index] != source:
                issues.append("MUTATION_SOURCE_AA_MISMATCH")
            expected[index] = target
            if index >= len(aa_sequence) or aa_sequence[index] != target:
                issues.append("MUTATION_TARGET_AA_MISMATCH")
        if len(aa_sequence) == len(wt_sequence) and "".join(expected) != aa_sequence:
            issues.append("UNLABELED_SEQUENCE_DIFFERENCE")
    return not issues, category, issues


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or denominator == 0:
        return None
    result = numerator / denominator
    return result if math.isfinite(result) else None


def value_at_bound(
    value: float | None,
    bound: float | None,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> bool:
    if value is None or bound is None:
        return False
    return math.isclose(
        float(value),
        float(bound),
        rel_tol=relative_tolerance,
        abs_tol=absolute_tolerance,
    )


def max_min_fold(values: Iterable[float | None]) -> float | None:
    clean = [
        float(value)
        for value in values
        if value is not None and math.isfinite(value) and value > 0
    ]
    return max(clean) / min(clean) if len(clean) >= 2 else None


def join_and_qc(
    lsa_records: list[dict[str, Any]],
    ngs_records: list[dict[str, Any]],
    wt_sequence: str,
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    ngs_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in ngs_records:
        ngs_groups[record["well"]].append(record)

    accepted_codes = set(config["accepted_ngs_call_codes"])
    pass_flags = set(config["lsa_pass_flags"])
    review_flags = set(config["lsa_review_flags"])
    complex_flags = set(config.get("lsa_complex_flags", ["Complex"]))
    inactive_flags = set(config["lsa_inactive_flags"])
    fail_flags = set(config["lsa_fail_flags"])
    review_threshold = to_float(config.get("relative_error_review_threshold"))
    fail_threshold = to_float(config.get("relative_error_fail_threshold"))
    if (
        review_threshold is not None
        and fail_threshold is not None
        and fail_threshold < review_threshold
    ):
        raise ValueError("relative_error_fail_threshold must be >= review threshold")
    kd_tolerance = float(config["kd_identity_relative_tolerance"])
    boundary_relative_tolerance = float(config["fit_boundary_relative_tolerance"])
    boundary_absolute_tolerance = float(config["fit_boundary_absolute_tolerance"])
    parameter_bounds = dict(config.get("fit_parameter_bounds", {}))
    integrated: list[dict[str, Any]] = []

    for kinetic in lsa_records:
        fail_issues: list[str] = []
        review_issues: list[str] = []
        kinetic_fail_issues: list[str] = []
        kinetic_review_issues: list[str] = []

        def add_kinetic_issue(code: str, severity: str) -> None:
            add_issue(kinetic_fail_issues, kinetic_review_issues, code, severity)
            add_issue(fail_issues, review_issues, code, severity)

        candidates = ngs_groups.get(kinetic["normalized_well"], [])
        if not candidates:
            ngs: dict[str, Any] = {}
            add_issue(fail_issues, review_issues, "NGS_KEY_NOT_FOUND", "fail")
        else:
            ngs = candidates[0]
            if len(candidates) > 1:
                add_issue(fail_issues, review_issues, "NGS_DUPLICATE_KEY", "fail")

        genotype_qc_status = "PASS"
        if not candidates or len(candidates) > 1:
            genotype_qc_status = "FAIL"
        if ngs.get("call_code") not in accepted_codes:
            code = ngs.get("call_code") or "MISSING"
            add_issue(fail_issues, review_issues, f"NGS_CALL_CODE_{code}", "fail")
            genotype_qc_status = "FAIL"
        if ngs.get("mutation_parse_error"):
            add_issue(fail_issues, review_issues, "MUTATION_LIST_PARSE_ERROR", "fail")
            genotype_qc_status = "FAIL"

        inter_cluster_score = to_float(ngs.get("inter_cluster_mixing_score"))
        intra_cluster_score = to_float(ngs.get("intra_cluster_mixing_score"))
        ngs_read_count = to_float(
            ngs.get("total_valid_read_count", ngs.get("read_count"))
        )
        inter_cluster_threshold = float(
            config.get("ngs_inter_cluster_mixing_threshold", 0.10)
        )
        intra_cluster_threshold = float(
            config.get("ngs_intra_cluster_mixing_threshold", 0.20)
        )
        if ngs_read_count is not None and ngs_read_count <= 0:
            ngs_clone_score_status = "NOT_INTERPRETABLE_ZERO_READ_COUNT"
            if ngs.get("call_code") == "M":
                ngs_score_call_consistency = "CONSISTENT_M_CALL_ZERO_READ_COUNT"
            else:
                ngs_score_call_consistency = "FAIL_C_CALL_ZERO_READ_COUNT"
                add_issue(
                    fail_issues,
                    review_issues,
                    "NGS_C_CALL_WITH_ZERO_READ_COUNT",
                    "fail",
                )
                genotype_qc_status = "FAIL"
        elif inter_cluster_score is None or intra_cluster_score is None:
            ngs_clone_score_status = "NOT_EVALUATED_MISSING_SCORE"
            ngs_score_call_consistency = "NOT_EVALUATED"
        elif not (
            0 <= inter_cluster_score <= 1 and 0 <= intra_cluster_score <= 1
        ):
            ngs_clone_score_status = "FAIL_SCORE_OUTSIDE_0_TO_1"
            ngs_score_call_consistency = "FAIL_INVALID_SCORE_RANGE"
            add_issue(
                fail_issues,
                review_issues,
                "NGS_CLUSTER_MIXING_SCORE_OUTSIDE_0_TO_1",
                "fail",
            )
            genotype_qc_status = "FAIL"
        else:
            scores_strictly_pass = (
                inter_cluster_score < inter_cluster_threshold
                and intra_cluster_score < intra_cluster_threshold
            )
            scores_at_boundary = (
                math.isclose(inter_cluster_score, inter_cluster_threshold)
                or math.isclose(intra_cluster_score, intra_cluster_threshold)
            )
            if scores_strictly_pass:
                ngs_clone_score_status = "PASS_PROVISIONAL"
            elif scores_at_boundary:
                ngs_clone_score_status = "BORDERLINE_AT_REPORTED_PRECISION"
            else:
                ngs_clone_score_status = "FAIL_PROVISIONAL"
            if ngs.get("call_code") == "C" and not scores_strictly_pass:
                ngs_score_call_consistency = "REVIEW_C_CALL_NOT_STRICTLY_BELOW_THRESHOLDS"
                add_issue(
                    fail_issues,
                    review_issues,
                    "NGS_C_CALL_NOT_STRICTLY_BELOW_PROVISIONAL_SCORE_THRESHOLDS",
                    "review",
                )
            elif ngs.get("call_code") == "M" and scores_strictly_pass:
                ngs_score_call_consistency = (
                    "M_CALL_DESPITE_PASSING_SCORES_OTHER_PIPELINE_CRITERIA_POSSIBLE"
                )
            else:
                ngs_score_call_consistency = "CONSISTENT_WITH_PROVISIONAL_THRESHOLDS"

        sequence_valid, category, sequence_issues = validate_sequence_and_mutations(
            wt_sequence,
            ngs.get("aa_sequence", ""),
            ngs.get("mutation_labels", []),
            bool(config["require_explicit_substitution_labels"]),
        )
        for issue in sequence_issues:
            add_issue(fail_issues, review_issues, issue, "fail")
        sequence_qc_status = "PASS" if sequence_valid else "FAIL"

        lsa_flag = kinetic.get("lsa_flag", "")
        is_inactive = lsa_flag in inactive_flags
        is_complex = lsa_flag in complex_flags
        assay_qc_status = "INTERPRETABLE"
        if lsa_flag in fail_flags:
            add_kinetic_issue(f"LSA_FLAG_{lsa_flag.upper()}", "fail")
            assay_qc_status = "FAIL"
        elif is_complex:
            add_kinetic_issue("LSA_FLAG_COMPLEX", "fail")
            assay_qc_status = "FAIL_QUANTITATIVE_COMPLEX"
        elif lsa_flag in review_flags:
            add_kinetic_issue(f"LSA_FLAG_{lsa_flag.upper()}", "review")
            assay_qc_status = "REVIEW"
        elif is_inactive:
            assay_qc_status = "INTERPRETABLE_INACTIVE"
        elif lsa_flag not in pass_flags:
            add_kinetic_issue(
                f"LSA_FLAG_UNRECOGNIZED_{lsa_flag or 'BLANK'}", "review"
            )
            assay_qc_status = "REVIEW"

        required_metrics = ["kon_M_inv_s", "koff_s_inv", "KD_M"]
        required_values_present = all(kinetic.get(field) is not None for field in required_metrics)
        required_values_positive = required_values_present and all(
            kinetic[field] > 0 for field in required_metrics
        )
        calculated_kd = safe_ratio(kinetic.get("koff_s_inv"), kinetic.get("kon_M_inv_s"))
        kd_relative_difference = None
        kd_error_propagated_M = None
        kd_error_propagation_relative_difference = None
        relative_errors: dict[str, float | None] = {
            "kon_relative_fit_error": None,
            "koff_relative_fit_error": None,
            "KD_relative_fit_error": None,
            "Rmax_relative_fit_error": None,
        }
        bound_hits = {
            "kon_at_lower_fit_bound": False,
            "kon_at_upper_fit_bound": False,
            "koff_at_lower_fit_bound": False,
            "koff_at_upper_fit_bound": False,
            "Rmax_at_lower_fit_bound": False,
            "Rmax_at_upper_fit_bound": False,
        }
        residual_sd_over_rmax = None

        if is_inactive:
            binding_outcome = "INACTIVE_NO_DETECTABLE_BINDING"
            kinetic_quantification_status = "NOT_QUANTIFIABLE_INACTIVE"
            if any(kinetic.get(field) is not None for field in required_metrics):
                add_kinetic_issue("INACTIVE_WITH_PARTIAL_KINETIC_VALUES", "review")
        else:
            for field in required_metrics:
                value = kinetic.get(field)
                if value is None:
                    add_kinetic_issue(f"{field.upper()}_MISSING", "fail")
                elif value <= 0:
                    add_kinetic_issue(f"{field.upper()}_NONPOSITIVE", "fail")
            if config["require_positive_rmax"]:
                rmax = kinetic.get("Rmax_RU")
                if rmax is None:
                    add_kinetic_issue("RMAX_MISSING", "fail")
                elif rmax <= 0:
                    add_kinetic_issue("RMAX_NONPOSITIVE", "fail")

            if calculated_kd is not None and kinetic.get("KD_M") not in (None, 0):
                kd_relative_difference = abs(calculated_kd - kinetic["KD_M"]) / abs(
                    kinetic["KD_M"]
                )
                if kd_relative_difference > kd_tolerance:
                    add_kinetic_issue("KD_IDENTITY_MISMATCH", "fail")

            error_specs = [
                (
                    "kon_relative_fit_error",
                    "kon_M_inv_s",
                    "kon_fit_error_M_inv_s",
                ),
                (
                    "koff_relative_fit_error",
                    "koff_s_inv",
                    "koff_fit_error_s_inv",
                ),
                ("KD_relative_fit_error", "KD_M", "KD_fit_error_M"),
                ("Rmax_relative_fit_error", "Rmax_RU", "Rmax_fit_error_RU"),
            ]
            legacy_error_fields = {
                "kon_fit_error_M_inv_s": "kon_se_M_inv_s",
                "koff_fit_error_s_inv": "koff_se_s_inv",
                "KD_fit_error_M": "KD_se_M",
                "Rmax_fit_error_RU": "Rmax_se_RU",
            }
            fit_errors: dict[str, float | None] = {}
            for name, value_field, error_field in error_specs:
                error = kinetic.get(error_field)
                if error is None and legacy_error_fields[error_field] in kinetic:
                    error = kinetic.get(legacy_error_fields[error_field])
                fit_errors[error_field] = error
                value = kinetic.get(value_field)
                if error is None:
                    if config.get("require_fit_errors_for_active", True):
                        add_kinetic_issue(f"{error_field.upper()}_MISSING", "fail")
                    continue
                if error < 0:
                    add_kinetic_issue(f"{error_field.upper()}_NEGATIVE", "fail")
                    continue
                if error == 0 and config.get("zero_fit_error_is_fail", True):
                    add_kinetic_issue(
                        f"{error_field.upper()}_ZERO_BOUND_OR_UNAVAILABLE", "fail"
                    )
                relative = safe_ratio(error, abs(value) if value is not None else None)
                relative_errors[name] = relative
                if (
                    relative is not None
                    and fail_threshold is not None
                    and relative >= fail_threshold
                ):
                    add_kinetic_issue(
                        f"{name.upper()}_GE_{fail_threshold:g}", "fail"
                    )
                elif (
                    relative is not None
                    and review_threshold is not None
                    and relative >= review_threshold
                ):
                    add_kinetic_issue(
                        f"{name.upper()}_GE_{review_threshold:g}", "review"
                    )

            propagation_inputs = [
                kinetic.get("kon_M_inv_s"),
                kinetic.get("koff_s_inv"),
                kinetic.get("KD_M"),
                fit_errors.get("kon_fit_error_M_inv_s"),
                fit_errors.get("koff_fit_error_s_inv"),
                fit_errors.get("KD_fit_error_M"),
            ]
            kd_error_propagation_qc_status = (
                "NOT_EVALUATED_EXPORT_ERROR_SEMANTICS_UNKNOWN"
            )
            if config.get("per_fit_kd_error_propagation_check_enabled", False) and all(
                value is not None and value > 0 for value in propagation_inputs
            ):
                kon, koff, kd_value, kon_error, koff_error, kd_error = propagation_inputs
                kd_error_propagated_M = kd_value * math.sqrt(
                    (koff_error / koff) ** 2 + (kon_error / kon) ** 2
                )
                kd_error_propagation_relative_difference = abs(
                    kd_error_propagated_M - kd_error
                ) / kd_error
                propagation_tolerance = float(
                    config["kd_error_propagation_relative_tolerance"]
                )
                if kd_error_propagation_relative_difference > propagation_tolerance:
                    add_kinetic_issue("KD_ERROR_PROPAGATION_MISMATCH", "fail")
                    kd_error_propagation_qc_status = "FAIL"
                else:
                    kd_error_propagation_qc_status = "PASS"

            bound_specs = [
                ("kon", "kon_M_inv_s"),
                ("koff", "koff_s_inv"),
                ("Rmax", "Rmax_RU"),
            ]
            for label, field in bound_specs:
                bounds = parameter_bounds.get(field, {}) or {}
                for side in ("lower", "upper"):
                    bound_key = "min" if side == "lower" else "max"
                    bound = to_float(bounds.get(bound_key))
                    hit_key = f"{label}_at_{side}_fit_bound"
                    hit = value_at_bound(
                        kinetic.get(field),
                        bound,
                        boundary_relative_tolerance,
                        boundary_absolute_tolerance,
                    )
                    bound_hits[hit_key] = hit
                    if hit:
                        add_kinetic_issue(
                            f"{label.upper()}_AT_{side.upper()}_FIT_BOUND", "fail"
                        )

            rmax = kinetic.get("Rmax_RU")
            rmax_review_threshold = to_float(config.get("rmax_review_threshold_RU"))
            if (
                rmax is not None
                and rmax > 0
                and rmax_review_threshold is not None
                and rmax < rmax_review_threshold
            ):
                add_kinetic_issue(
                    f"RMAX_LT_{rmax_review_threshold:g}_RU_PROVISIONAL", "review"
                )

            residual_sd = kinetic.get("residual_sd_RU")
            if residual_sd is None:
                if config.get("require_residual_sd_for_active", True):
                    add_kinetic_issue("RESIDUAL_SD_MISSING", "fail")
            elif residual_sd < 0:
                add_kinetic_issue("RESIDUAL_SD_NEGATIVE", "fail")
            elif rmax is not None and rmax > 0:
                residual_sd_over_rmax = safe_ratio(residual_sd, rmax)
                residual_threshold = to_float(
                    config.get("residual_sd_over_rmax_fail_threshold")
                )
                if (
                    residual_threshold is not None
                    and residual_sd_over_rmax is not None
                    and residual_sd_over_rmax > residual_threshold
                ):
                    add_kinetic_issue(
                        f"RESIDUAL_SD_OVER_RMAX_GT_{residual_threshold:g}",
                        "fail",
                    )

            any_bound_hit = any(bound_hits.values())
            if is_complex:
                binding_outcome = (
                    "COMPLEX_FIT_QUANTIFIED"
                    if required_values_positive
                    else "COMPLEX_FIT_UNQUANTIFIED"
                )
                kinetic_quantification_status = (
                    "QUALITATIVE_ONLY_COMPLEX"
                    if required_values_positive
                    else "NOT_QUANTIFIABLE"
                )
            elif required_values_positive and any_bound_hit:
                binding_outcome = "ACTIVE_BOUND_CENSORED"
                kinetic_quantification_status = "NOT_QUANTIFIABLE_BOUND"
            elif required_values_positive:
                binding_outcome = "ACTIVE_QUANTIFIED"
                kinetic_quantification_status = "QUANTIFIED"
            else:
                binding_outcome = "ACTIVE_UNQUANTIFIED_INVALID"
                kinetic_quantification_status = "NOT_QUANTIFIABLE"

        kinetic_fit_qc_status = (
            "NOT_APPLICABLE_INACTIVE"
            if is_inactive
            else (
                "FAIL"
                if kinetic_fail_issues
                else ("REVIEW" if kinetic_review_issues else "PASS")
            )
        )

        ngs_read_qc_status = "NOT_EVALUATED_NO_THRESHOLD"
        low_read_interpretation = ""
        min_read_count = to_float(config.get("ngs_min_read_count_review"))
        if min_read_count is not None and ngs.get("read_count") is not None:
            if float(ngs["read_count"]) < min_read_count:
                ngs_read_qc_status = "LOW"
                if kinetic_fit_qc_status == "PASS":
                    low_read_interpretation = "KINETICS_SUPPORTS_USABLE_MEASUREMENT"
                elif config.get("colony_pcr_failure_inference_enabled", False):
                    low_read_interpretation = (
                        "PCR_OR_EXPRESSION_FAILURE_SUSPECTED_NOT_PROVEN"
                    )
                    add_issue(
                        fail_issues,
                        review_issues,
                        "LOW_NGS_READ_AND_POOR_OR_INACTIVE_LSA_PCR_OR_EXPRESSION_FAILURE_SUSPECTED",
                        "review",
                    )
            else:
                ngs_read_qc_status = "PASS"

        status = "FAIL" if fail_issues else ("REVIEW" if review_issues else "PASS")
        visualization_class = (
            "FAIL"
            if status == "FAIL"
            else (
                "INACTIVE"
                if binding_outcome == "INACTIVE_NO_DETECTABLE_BINDING"
                else status
            )
        )
        visualization_color = config.get("visualization_colors", {}).get(
            visualization_class, ""
        )
        genotype_and_sequence_valid = genotype_qc_status == "PASS" and sequence_valid
        phenotype_informative = genotype_and_sequence_valid and (
            binding_outcome in {
                "ACTIVE_QUANTIFIED",
                "INACTIVE_NO_DETECTABLE_BINDING",
            }
            or binding_outcome.startswith("COMPLEX_FIT")
        )
        if not genotype_and_sequence_valid:
            analysis_eligibility = "EXCLUDE_UNCERTAIN_GENOTYPE_OR_SEQUENCE"
        elif binding_outcome == "INACTIVE_NO_DETECTABLE_BINDING":
            analysis_eligibility = "CENSORED_INACTIVE"
        elif binding_outcome.startswith("COMPLEX_FIT"):
            analysis_eligibility = "QUALITATIVE_BINDING_ONLY_COMPLEX"
        elif kinetic_fit_qc_status == "FAIL":
            analysis_eligibility = "EXCLUDE_UNRELIABLE_KINETIC_FIT"
        elif binding_outcome == "ACTIVE_QUANTIFIED" and status == "PASS":
            analysis_eligibility = "QUANTITATIVE"
        elif status == "REVIEW":
            analysis_eligibility = "REVIEW"
        else:
            analysis_eligibility = "EXCLUDE_TECHNICAL"
        labels = ngs.get("mutation_labels", [])
        integrated.append(
            {
                **kinetic,
                "ngs_source_well": ngs.get("source_well", ngs.get("well", "")),
                "ngs_well": ngs.get("well", ""),
                "source_ngs_row": ngs.get("source_ngs_row", ""),
                "ngs_call_code": ngs.get("call_code", ""),
                "ngs_total_valid_read_count": ngs_read_count,
                "ngs_read_count": ngs_read_count,
                "ngs_inter_cluster_mixing_score": inter_cluster_score,
                "ngs_intra_cluster_mixing_score": intra_cluster_score,
                "ngs_clone_score_status": ngs_clone_score_status,
                "ngs_score_call_consistency": ngs_score_call_consistency,
                "ngs_metric_1": ngs.get("metric_1"),
                "ngs_metric_2": ngs.get("metric_2"),
                "ngs_read_qc_status": ngs_read_qc_status,
                "low_read_interpretation": low_read_interpretation,
                "dna_sequence": ngs.get("dna_sequence", ""),
                "mutation_labels": "|".join(labels),
                "mutation_labels_json": json.dumps(labels, ensure_ascii=False),
                "mutation_count": 0 if labels == ["WT"] else len(labels),
                "variant_category": category,
                "aa_sequence": ngs.get("aa_sequence", ""),
                "aa_length": len(ngs.get("aa_sequence", "")),
                "genotype_qc_status": genotype_qc_status,
                "sequence_qc_status": sequence_qc_status,
                "assay_qc_status": assay_qc_status,
                "kinetic_fit_qc_status": kinetic_fit_qc_status,
                "kinetic_fit_fail_reasons": ";".join(kinetic_fail_issues),
                "kinetic_fit_review_reasons": ";".join(kinetic_review_issues),
                "binding_outcome": binding_outcome,
                "kinetic_quantification_status": kinetic_quantification_status,
                "phenotype_informative": phenotype_informative,
                "analysis_eligibility": analysis_eligibility,
                "sequence_mapping_valid": sequence_valid,
                "KD_calculated_from_rates_M": calculated_kd,
                "KD_identity_relative_difference": kd_relative_difference,
                "KD_fit_error_propagated_from_rates_M": kd_error_propagated_M,
                "KD_error_propagation_relative_difference": (
                    kd_error_propagation_relative_difference
                ),
                "KD_error_propagation_qc_status": (
                    kd_error_propagation_qc_status
                    if not is_inactive
                    else "NOT_APPLICABLE_INACTIVE"
                ),
                **relative_errors,
                "kon_relative_se": relative_errors["kon_relative_fit_error"],
                "koff_relative_se": relative_errors["koff_relative_fit_error"],
                "KD_relative_se": relative_errors["KD_relative_fit_error"],
                "Rmax_relative_se": relative_errors["Rmax_relative_fit_error"],
                "residual_sd_over_Rmax": residual_sd_over_rmax,
                **bound_hits,
                "saturation_qc_status": "UNAVAILABLE_IN_KINETIC_SNAPSHOT_EXPORT",
                "inactive_response_threshold_RU_manual_default": config.get(
                    "inactive_default_min_binding_response_RU"
                ),
                "inactive_response_threshold_confirmed_for_experiment": bool(
                    config.get(
                        "inactive_response_threshold_confirmed_for_experiment",
                        False,
                    )
                ),
                "qc_status": status,
                "overall_qc_status": status,
                "visualization_class": visualization_class,
                "visualization_color": visualization_color,
                "qc_fail_reasons": ";".join(fail_issues),
                "qc_review_reasons": ";".join(review_issues),
                "technical_qc_status": status,
                "structure_mapping_ready": genotype_and_sequence_valid,
            }
        )
    return integrated


def median_or_none(values: Iterable[float | None]) -> float | None:
    clean = [float(value) for value in values if value is not None and math.isfinite(value)]
    return statistics.median(clean) if clean else None


def mean_or_none(values: Iterable[float | None]) -> float | None:
    clean = [float(value) for value in values if value is not None and math.isfinite(value)]
    return statistics.mean(clean) if clean else None


def sample_sd_or_none(values: Iterable[float | None]) -> float | None:
    clean = [float(value) for value in values if value is not None and math.isfinite(value)]
    return statistics.stdev(clean) if len(clean) >= 2 else None


def roi_list(records: list[dict[str, Any]]) -> str:
    """List original Carterra labels once and retain evidence of unlabeled rows."""
    labels = {"" if r.get("roi") is None else str(r["roi"]).strip() for r in records}
    missing = "" in labels
    labels.discard("")
    ordered = sorted(labels, key=lambda s: (0, int(s), s) if s.isdecimal() else (1, 0, s))
    if missing:
        ordered.append("unreported")
    return ", ".join(ordered)


def summarize_variants(
    records: list[dict[str, Any]], wt_sequence: str, config: dict[str, Any]
) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        if record["sequence_mapping_valid"] and record["aa_sequence"]:
            groups[record["aa_sequence"]].append(record)

    strict_wt = [
        record
        for record in records
        if record["aa_sequence"] == wt_sequence
        and record["genotype_qc_status"] == "PASS"
        and record["binding_outcome"] == "ACTIVE_QUANTIFIED"
        and record["kinetic_fit_qc_status"] == "PASS"
    ]
    wt_medians = {
        field: median_or_none(record[field] for record in strict_wt)
        for field in ["kon_M_inv_s", "koff_s_inv", "KD_M", "Rmax_RU", "residual_sd_RU"]
    }
    temperature = float(config["temperature_kelvin"])
    minimum_active_pass = int(
        config.get("minimum_active_pass_replicates_for_quantitation", 1)
    )
    summaries: list[dict[str, Any]] = []
    for aa_sequence, group in groups.items():
        genotype_valid = [
            record for record in group if record["genotype_qc_status"] == "PASS"
        ]
        active = [
            record
            for record in genotype_valid
            if record["binding_outcome"] == "ACTIVE_QUANTIFIED"
        ]
        inactive = [
            record
            for record in genotype_valid
            if record["binding_outcome"] == "INACTIVE_NO_DETECTABLE_BINDING"
        ]
        complex_fit = [
            record
            for record in genotype_valid
            if record["binding_outcome"].startswith("COMPLEX_FIT")
        ]
        active_qc_limited = [
            record
            for record in genotype_valid
            if record["binding_outcome"]
            in {
                "ACTIVE_BOUND_CENSORED",
                "ACTIVE_UNQUANTIFIED_INVALID",
            }
            or record["binding_outcome"].startswith("COMPLEX_FIT")
        ]
        strict = [
            record for record in active if record["kinetic_fit_qc_status"] == "PASS"
        ]
        strict_for_summary = strict if len(strict) >= minimum_active_pass else []
        informative_count = len(strict) + len(inactive)
        any_active_observation = bool(active or active_qc_limited)
        if informative_count == 0 and not any_active_observation:
            phenotype_consensus = "unresolved"
        elif inactive and not strict and not any_active_observation:
            phenotype_consensus = "consistently_inactive"
        elif strict and not inactive:
            phenotype_consensus = "consistently_active"
        elif inactive and (strict or any_active_observation):
            phenotype_consensus = "mixed_active_inactive"
        else:
            phenotype_consensus = "unresolved"

        if phenotype_consensus == "mixed_active_inactive":
            variant_qc_status = "REVIEW"
            replicate_interpretation = (
                "USE_PASS_ACTIVE_FOR_PRIMARY_QUANTITATION_KEEP_INACTIVE_AS_DISCORDANCE"
                if strict_for_summary
                else "MIXED_WITHOUT_PASS_ACTIVE_FIT"
            )
        elif phenotype_consensus == "consistently_inactive":
            variant_qc_status = "REVIEW"
            replicate_interpretation = (
                "CONSISTENT_NO_DETECTABLE_RESPONSE_EXPRESSION_OR_SPOT_VALIDITY_UNCONTROLLED"
            )
        elif phenotype_consensus == "consistently_active" and strict_for_summary:
            variant_qc_status = "PASS"
            replicate_interpretation = "QUANTIFIED_FROM_PASS_ACTIVE"
        else:
            variant_qc_status = "REVIEW"
            replicate_interpretation = "UNRESOLVED_OR_NO_PASS_ACTIVE_FIT"
        medians = {
            field: median_or_none(record[field] for record in strict_for_summary)
            for field in ["kon_M_inv_s", "koff_s_inv", "KD_M", "Rmax_RU", "residual_sd_RU"]
        }
        means = {
            field: mean_or_none(record[field] for record in strict_for_summary)
            for field in ["kon_M_inv_s", "koff_s_inv", "KD_M", "Rmax_RU", "residual_sd_RU"]
        }
        sample_sds = {
            field: sample_sd_or_none(record[field] for record in strict_for_summary)
            for field in ["kon_M_inv_s", "koff_s_inv", "KD_M", "Rmax_RU", "residual_sd_RU"]
        }
        kd_sd_propagated_from_rates = None
        if (
            all(
                value is not None and value > 0
                for value in [
                    means["KD_M"],
                    means["kon_M_inv_s"],
                    means["koff_s_inv"],
                ]
            )
            and all(
                value is not None and value >= 0
                for value in [
                    sample_sds["kon_M_inv_s"],
                    sample_sds["koff_s_inv"],
                ]
            )
        ):
            kd_sd_propagated_from_rates = means["KD_M"] * math.sqrt(
                (sample_sds["koff_s_inv"] / means["koff_s_inv"]) ** 2
                + (sample_sds["kon_M_inv_s"] / means["kon_M_inv_s"]) ** 2
            )

        pass_active_override = bool(strict_for_summary and inactive)
        inactive_loss_evidence_count = (
            len(inactive) if phenotype_consensus == "consistently_inactive" else 0
        )
        inactive_excluded_by_active_override = (
            len(inactive) if pass_active_override else 0
        )
        if pass_active_override:
            primary_phenotype = "ACTIVE_SUPPORTED_WITH_INACTIVE_DISCORDANCE"
            inactive_handling = (
                "RETAIN_RAW_BUT_EXCLUDE_FROM_LOSS_OF_BINDING_EVIDENCE"
            )
        elif phenotype_consensus == "consistently_inactive":
            primary_phenotype = "NO_DETECTABLE_RESPONSE_UNRESOLVED_CAUSE"
            inactive_handling = (
                "RETAIN_AS_QUALITATIVE_LOSS_EVIDENCE_EXPRESSION_OR_SURFACE_FAILURE_NOT_EXCLUDED"
            )
        elif phenotype_consensus == "consistently_active":
            primary_phenotype = "ACTIVE_SUPPORTED"
            inactive_handling = "NOT_APPLICABLE"
        else:
            primary_phenotype = "UNRESOLVED"
            inactive_handling = "NO_OVERRIDE_WITHOUT_PASS_ACTIVE"
        replicate_folds = {
            "KD_max_min_fold_pass": max_min_fold(
                record["KD_M"] for record in strict_for_summary
            ),
            "kon_max_min_fold_pass": max_min_fold(
                record["kon_M_inv_s"] for record in strict_for_summary
            ),
            "koff_max_min_fold_pass": max_min_fold(
                record["koff_s_inv"] for record in strict_for_summary
            ),
        }
        observed_replicate_folds = [
            value for value in replicate_folds.values() if value is not None
        ]
        consistency_review_fold = to_float(
            config.get("replicate_consistency_review_fold")
        )
        consistency_fail_fold = to_float(
            config.get("replicate_consistency_fail_fold")
        )
        if len(strict_for_summary) < 2:
            replicate_consistency_qc_status = "NOT_ASSESSABLE_N_LT_2"
        elif consistency_review_fold is None and consistency_fail_fold is None:
            replicate_consistency_qc_status = "REPORTED_NO_ACCEPTANCE_THRESHOLD"
        elif consistency_fail_fold is not None and any(
            value >= consistency_fail_fold for value in observed_replicate_folds
        ):
            replicate_consistency_qc_status = "FAIL"
            variant_qc_status = "FAIL"
        elif consistency_review_fold is not None and any(
            value >= consistency_review_fold for value in observed_replicate_folds
        ):
            replicate_consistency_qc_status = "REVIEW"
            if variant_qc_status == "PASS":
                variant_qc_status = "REVIEW"
        else:
            replicate_consistency_qc_status = "PASS"
        mutation_label_sets = sorted(
            {record["mutation_labels"] for record in group if record["mutation_labels"]}
        )
        category = Counter(record["variant_category"] for record in group).most_common(1)[0][0]
        kd_fold = safe_ratio(medians["KD_M"], wt_medians["KD_M"])
        ddg = (
            GAS_CONSTANT_KCAL * temperature * math.log(kd_fold)
            if kd_fold is not None and kd_fold > 0
            else None
        )
        variant_visualization_class = (
            "MIXED_REVIEW"
            if phenotype_consensus == "mixed_active_inactive"
            else variant_qc_status
        )
        summaries.append(
            {
                "variant_category": category,
                "mutation_labels": ";".join(mutation_label_sets),
                "aa_sequence": aa_sequence,
                "aa_length": len(aa_sequence),
                "observed_rois": roi_list(genotype_valid),
                "median_rois": roi_list(strict_for_summary),
                "n_total": len(group),
                "n_technical_pass": sum(
                    record["technical_qc_status"] == "PASS" for record in group
                ),
                "n_technical_review": sum(
                    record["technical_qc_status"] == "REVIEW" for record in group
                ),
                "n_technical_fail": sum(
                    record["technical_qc_status"] == "FAIL" for record in group
                ),
                "n_genotype_pass": len(genotype_valid),
                "n_review_fits": sum(r["kinetic_fit_qc_status"] == "REVIEW" for r in genotype_valid),
                "n_active_quantified": len(active),
                "n_active_quantitative_pass": len(strict),
                "n_active_qc_limited": len(active_qc_limited),
                "n_inactive": len(inactive),
                "n_complex": len(complex_fit),
                "n_informative_active_or_inactive": informative_count,
                "inactive_fraction": safe_ratio(len(inactive), informative_count),
                "phenotype_consensus": phenotype_consensus,
                "primary_phenotype": primary_phenotype,
                "inactive_handling": inactive_handling,
                "n_inactive_loss_of_binding_evidence": inactive_loss_evidence_count,
                "n_inactive_excluded_by_pass_active_override": (
                    inactive_excluded_by_active_override
                ),
                "variant_qc_status": variant_qc_status,
                "visualization_class": variant_visualization_class,
                "visualization_color": config.get("visualization_colors", {}).get(
                    variant_visualization_class, ""
                ),
                "replicate_interpretation": replicate_interpretation,
                "replicate_consistency_qc_status": replicate_consistency_qc_status,
                "minimum_active_pass_replicates_required": minimum_active_pass,
                **replicate_folds,
                "median_kon_M_inv_s_pass": medians["kon_M_inv_s"],
                "median_koff_s_inv_pass": medians["koff_s_inv"],
                "median_KD_M_pass": medians["KD_M"],
                "median_Rmax_RU_pass": medians["Rmax_RU"],
                "median_residual_sd_RU_pass": medians["residual_sd_RU"],
                "mean_kon_M_inv_s_pass": means["kon_M_inv_s"],
                "mean_koff_s_inv_pass": means["koff_s_inv"],
                "mean_KD_M_pass": means["KD_M"],
                "mean_Rmax_RU_pass": means["Rmax_RU"],
                "mean_residual_sd_RU_pass": means["residual_sd_RU"],
                "sample_sd_kon_M_inv_s_pass": sample_sds["kon_M_inv_s"],
                "sample_sd_koff_s_inv_pass": sample_sds["koff_s_inv"],
                "sample_sd_KD_M_observed_pass": sample_sds["KD_M"],
                "sample_sd_Rmax_RU_pass": sample_sds["Rmax_RU"],
                "sample_sd_residual_sd_RU_pass": sample_sds["residual_sd_RU"],
                "sample_sd_KD_M_propagated_from_rates_pass": (
                    kd_sd_propagated_from_rates
                ),
                "kon_fold_vs_WT_pass": safe_ratio(medians["kon_M_inv_s"], wt_medians["kon_M_inv_s"]),
                "koff_fold_vs_WT_pass": safe_ratio(medians["koff_s_inv"], wt_medians["koff_s_inv"]),
                "KD_fold_vs_WT_pass": kd_fold,
                "delta_delta_G_kcal_mol_pass": ddg,
            }
        )
    return sorted(
        summaries,
        key=lambda row: (
            0 if row["variant_category"] == "WT" else 1,
            row["mutation_labels"],
            row["aa_sequence"],
        ),
    )


def residue_effects(variant_summaries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for variant in variant_summaries:
        if variant["variant_category"] not in {"single_substitution", "multi_substitution"}:
            continue
        label_sets = [item for item in variant["mutation_labels"].split(";") if item]
        labels = label_sets[0].split("|") if label_sets else []
        for label in labels:
            match = SUBSTITUTION_RE.fullmatch(label)
            if not match:
                continue
            source, position, target = match.groups()
            output.append(
                {
                    "sequence_position_1based": int(position),
                    "wt_aa": source,
                    "mutant_aa": target,
                    "mutation_label": label,
                    "effect_attribution": (
                        "direct_single_substitution"
                        if variant["variant_category"] == "single_substitution"
                        else "multi_substitution_confounded"
                    ),
                    "structure_chain_id": "",
                    "structure_residue_number": "",
                    **variant,
                }
            )
    return sorted(output, key=lambda row: (row["sequence_position_1based"], row["mutation_label"]))


def select_qc_examples(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    targets = [
        "PASS",
        "LSA_FLAG_COMPLEX",
        "INACTIVE_NO_DETECTABLE_BINDING",
        "NGS_CALL_CODE_M",
        "AA_SEQUENCE_MISSING",
        "SEQUENCE_LENGTH_MISMATCH",
        "KON_RELATIVE_FIT_ERROR_GE_0.5",
        "RESIDUAL_SD_OVER_RMAX_GT_0.1",
        "KOFF_AT_LOWER_FIT_BOUND",
        "RMAX_NONPOSITIVE",
    ]
    examples: list[dict[str, Any]] = []
    for target in targets:
        desired_status = (
            "PASS"
            if target == "PASS"
            else (
                "REVIEW"
                if target
                in {
                    "KON_RELATIVE_FIT_ERROR_GE_0.5",
                }
                else ("PASS" if target == "INACTIVE_NO_DETECTABLE_BINDING" else "FAIL")
            )
        )
        ordered_records = sorted(
            records, key=lambda record: record["qc_status"] != desired_status
        )
        for record in ordered_records:
            issues = set(
                item
                for field in ["qc_fail_reasons", "qc_review_reasons"]
                for item in record[field].split(";")
                if item
            )
            if (
                (target == "PASS" and record["qc_status"] == "PASS")
                or target in issues
                or record["binding_outcome"] == target
            ):
                examples.append(
                    {
                        "example_issue": target,
                        "normalized_well": record["normalized_well"],
                        "ligand_name": record["ligand_name"],
                        "qc_status": record["qc_status"],
                        "qc_fail_reasons": record["qc_fail_reasons"],
                        "qc_review_reasons": record["qc_review_reasons"],
                        "lsa_flag": record["lsa_flag"],
                        "ngs_call_code": record["ngs_call_code"],
                        "genotype_qc_status": record["genotype_qc_status"],
                        "binding_outcome": record["binding_outcome"],
                        "analysis_eligibility": record["analysis_eligibility"],
                        "mutation_labels": record["mutation_labels"],
                        "aa_length": record["aa_length"],
                        "kon_M_inv_s": record["kon_M_inv_s"],
                        "koff_s_inv": record["koff_s_inv"],
                        "KD_M": record["KD_M"],
                        "Rmax_RU": record["Rmax_RU"],
                    }
                )
                break
    return examples


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not records:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for record in records:
        for key in record:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def counts(values: Iterable[Any]) -> dict[str, int]:
    return dict(sorted(Counter(str(value) for value in values).items()))


def summarize_internal_wt_control(
    records: list[dict[str, Any]], wt_sequence: str
) -> dict[str, Any]:
    accepted = [
        record
        for record in records
        if record.get("aa_sequence") == wt_sequence
        and record.get("genotype_qc_status") == "PASS"
        and record.get("sequence_qc_status") == "PASS"
        and record.get("binding_outcome") == "ACTIVE_QUANTIFIED"
        and record.get("kinetic_fit_qc_status") == "PASS"
    ]
    parameters = ["kon_M_inv_s", "koff_s_inv", "KD_M", "Rmax_RU"]
    parameter_summary = {
        field: {
            **finite_stats(record.get(field) for record in accepted),
            "max_min_fold": max_min_fold(record.get(field) for record in accepted),
        }
        for field in parameters
    }
    by_plate: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in accepted:
        match = re.match(r"^(P\d+)-", str(record.get("normalized_well", "")))
        by_plate[match.group(1) if match else "UNKNOWN"].append(record)
    per_plate = []
    for plate, group in sorted(by_plate.items()):
        per_plate.append(
            {
                "well_plate_proxy": plate,
                "n_pass_WT": len(group),
                **{
                    f"median_{field}": median_or_none(
                        record.get(field) for record in group
                    )
                    for field in parameters
                },
            }
        )
    return {
        "status": "BASELINE_ONLY_NO_PREDEFINED_ACCEPTANCE_RANGE",
        "n_pass_WT": len(accepted),
        "parameter_summary": parameter_summary,
        "per_well_plate_proxy": per_plate,
        "caution": (
            "Plate is inferred from the well key and is not a confirmed experimental batch. "
            "No historical WT acceptance range was supplied, so these values are descriptive only."
        ),
    }


def make_summary(
    lsa_path: Path,
    ngs_path: Path,
    lsa_records: list[dict[str, Any]],
    ngs_records: list[dict[str, Any]],
    integrated: list[dict[str, Any]],
    variants: list[dict[str, Any]],
    wt_summary: dict[str, Any],
    parser_warnings: list[str],
    config: dict[str, Any],
) -> dict[str, Any]:
    lsa_keys = {record["normalized_well"] for record in lsa_records if record["normalized_well"]}
    ngs_key_counts = Counter(record["well"] for record in ngs_records)
    joined = sum(record["ngs_well"] != "" for record in integrated)
    issue_counts: Counter[str] = Counter()
    kinetic_issue_counts: Counter[str] = Counter()
    for record in integrated:
        for field in ["qc_fail_reasons", "qc_review_reasons"]:
            issue_counts.update(item for item in record[field].split(";") if item)
        for field in ["kinetic_fit_fail_reasons", "kinetic_fit_review_reasons"]:
            kinetic_issue_counts.update(
                item for item in record[field].split(";") if item
            )
    kd_diffs = [
        record["KD_identity_relative_difference"]
        for record in integrated
        if record["KD_identity_relative_difference"] is not None
    ]
    plate_block = parse_ngs_plate_block(config.get("ngs_plate_block"))
    source_to_target_plate_map = {
        f"P{int(record['source_plate_number']):02d}": f"P{int(record['mapped_plate_number']):02d}"
        for record in ngs_records
        if record.get("source_plate_number") is not None
        and record.get("mapped_plate_number") is not None
    }
    return {
        "source_files": {
            "lsa_xlsx": {"name": lsa_path.name, "sha256": sha256(lsa_path)},
            "ngs_csv": {"name": ngs_path.name, "sha256": sha256(ngs_path)},
        },
        "input_shape": {
            "lsa_kinetic_rows": len(lsa_records),
            "ngs_rows": len(ngs_records),
            "ngs_columns": len(config["ngs_columns"]),
        },
        "ngs_plate_mapping": {
            "primer_barcode_type": config.get("ngs_primer_barcode_type") or None,
            "selected_source_plate_block": (
                f"P{plate_block[0]:02d}-P{plate_block[1]:02d}"
                if plate_block
                else None
            ),
            "target_plate_start": int(config.get("ngs_target_plate_start", 1)),
            "source_to_target_plate_map": dict(
                sorted(source_to_target_plate_map.items())
            ),
            "source_well_preserved_in_integrated_records": True,
        },
        "key_join": {
            "lsa_unique_normalized_keys": len(lsa_keys),
            "ngs_unique_keys": len(ngs_key_counts),
            "ngs_duplicate_key_count": sum(count > 1 for count in ngs_key_counts.values()),
            "joined_lsa_rows": joined,
            "unmatched_lsa_rows": len(integrated) - joined,
            "ngs_keys_outside_lsa_scope": len(set(ngs_key_counts) - lsa_keys),
        },
        "observed_codes": {
            "lsa_flags": counts(record["lsa_flag"] for record in lsa_records),
            "ngs_call_codes_all_rows": counts(record["call_code"] for record in ngs_records),
            "ngs_call_codes_joined": counts(record["ngs_call_code"] for record in integrated),
            "ngs_clone_score_status_joined": counts(
                record["ngs_clone_score_status"] for record in integrated
            ),
            "ngs_score_call_consistency_joined": counts(
                record["ngs_score_call_consistency"] for record in integrated
            ),
        },
        "qc": {
            "status_counts": counts(record["qc_status"] for record in integrated),
            "genotype_qc_counts": counts(
                record["genotype_qc_status"] for record in integrated
            ),
            "kinetic_fit_qc_counts": counts(
                record["kinetic_fit_qc_status"] for record in integrated
            ),
            "binding_outcome_counts": counts(
                record["binding_outcome"] for record in integrated
            ),
            "analysis_eligibility_counts": counts(
                record["analysis_eligibility"] for record in integrated
            ),
            "structure_mapping_ready_rows": sum(record["structure_mapping_ready"] for record in integrated),
            "sequence_mapping_valid_rows": sum(record["sequence_mapping_valid"] for record in integrated),
            "issue_counts": dict(sorted(issue_counts.items())),
            "kinetic_issue_counts": dict(sorted(kinetic_issue_counts.items())),
            "fit_boundary_hit_rows": sum(
                any(
                    record[field]
                    for field in [
                        "kon_at_lower_fit_bound",
                        "kon_at_upper_fit_bound",
                        "koff_at_lower_fit_bound",
                        "koff_at_upper_fit_bound",
                        "Rmax_at_lower_fit_bound",
                        "Rmax_at_upper_fit_bound",
                    ]
                )
                for record in integrated
            ),
            "saturation_qc_available_from_snapshot": False,
            "max_KD_identity_relative_difference": max(kd_diffs) if kd_diffs else None,
        },
        "variants": {
            "sequence_valid_unique_variants": len(variants),
            "categories": counts(variant["variant_category"] for variant in variants),
            "phenotype_consensus": counts(
                variant["phenotype_consensus"] for variant in variants
            ),
            "primary_phenotype": counts(
                variant["primary_phenotype"] for variant in variants
            ),
            "variant_qc_status": counts(
                variant["variant_qc_status"] for variant in variants
            ),
            "unique_variants_with_pass_measurement": sum(
                variant["n_active_quantitative_pass"] > 0 for variant in variants
            ),
            "unique_variants_with_inactive_observation": sum(
                variant["n_inactive"] > 0 for variant in variants
            ),
            "consistently_inactive_variants": sum(
                variant["phenotype_consensus"] == "consistently_inactive"
                for variant in variants
            ),
            "mixed_active_inactive_variants": sum(
                variant["phenotype_consensus"] == "mixed_active_inactive"
                for variant in variants
            ),
        },
        "wt_sequence": wt_summary,
        "wt_internal_control": summarize_internal_wt_control(
            integrated, wt_summary.get("sequence", "")
        ),
        "parser_warnings": parser_warnings,
        "active_config": config,
        "provisional_assumptions": [
            "The headerless NGS columns 4 and 5 are inter-cluster and intra-cluster mixing scores. Deprecated metric_1/metric_2 aliases are retained only for backward compatibility.",
            "NGS column 3 is total_valid_read_count. The provisional Clone rule is inter-cluster score <0.1 and intra-cluster score <0.2 after 98% identity clustering. Zero-total-valid-read rows are not interpreted as score passes. Boundary handling is unknown, so the supplied C/M call remains authoritative and boundary cases are REVIEW.",
            "NGS call code C is accepted and M is excluded because M fails the confirmed within-cluster homology reliability criterion.",
            "The Carterra export labels the uncertainty columns as Error; the manual does not define their exact statistical estimator, so Error/value is reported as relative fitting error rather than assumed to be a standard error.",
            "Relative fitting error >= 0.5 is REVIEW under the default policy; no automatic FAIL threshold is used until the Error statistic is defined and validated.",
            "Carterra Complex fits are excluded from quantitative summaries because the manual defines them as residual SD greater than 10% of fitted Rmax and instructs users to omit them from Stats means.",
            "C plus a valid sequence plus LSA Inactive is retained as a censored no-detectable-response phenotype, not a kinetic fit failure.",
            "Variant medians and fold changes use kinetic-fit PASS active records only; inactive observations are retained as dropout evidence and never converted to an artificial KD value.",
            "For a mixed active/inactive variant with at least one PASS active replicate, active kinetics are the primary quantitative evidence; the variant remains REVIEW and inactive records are retained but excluded from loss-of-binding hotspot evidence.",
            "All-inactive variants remain REVIEW because this experiment has no independent expression, immobilization, or spot-validity control.",
            "The run used a 10 RU low/no-binding Data Filter and a separate 20 RU response Highlighting threshold; neither is an Rmax threshold.",
            "The run enabled residual SD 10% of Rmax, Data <50% Rmax, kd-below-limit, km-below-limit, and Excluded Group highlighting; no sample was manually assigned to Excluded Group.",
            "The manual's KD propagation-of-errors formula is applied at replicate aggregation using rate means and sample SDs, not to validate per-fit Error columns.",
            "Kinetic Snapshot does not expose the software flag for observed binding below 50% of fitted Rmax; saturation adequacy requires raw sensorgrams or a richer export.",
        ],
    }


def render_qc_report(summary: dict[str, Any]) -> str:
    def fmt(value: Any, specification: str = ".4g") -> str:
        return format(float(value), specification) if value is not None else "n.a."

    join = summary["key_join"]
    qc = summary["qc"]
    wt = summary["wt_sequence"]
    wt_control = summary["wt_internal_control"]
    wt_kd = wt_control["parameter_summary"]["KD_M"]
    wt_plate_lines = "\n".join(
        f"- {row['well_plate_proxy']}: n={row['n_pass_WT']}, median KD={fmt(row['median_KD_M'])} M"
        for row in wt_control["per_well_plate_proxy"]
        if row["median_KD_M"] is not None
    )
    status_lines = "\n".join(
        f"- {key}: {value}" for key, value in qc["status_counts"].items()
    )
    kinetic_status_lines = "\n".join(
        f"- {key}: {value}" for key, value in qc["kinetic_fit_qc_counts"].items()
    )
    outcome_lines = "\n".join(
        f"- {key}: {value}" for key, value in qc["binding_outcome_counts"].items()
    )
    issue_lines = "\n".join(
        f"- {key}: {value}" for key, value in qc["issue_counts"].items()
    )
    kinetic_issue_lines = "\n".join(
        f"- {key}: {value}" for key, value in qc["kinetic_issue_counts"].items()
    )
    plate_mapping = summary.get("ngs_plate_mapping", {})
    plate_pairs = plate_mapping.get("source_to_target_plate_map", {})
    plate_pair_text = (
        ", ".join(f"{source}→{target}" for source, target in plate_pairs.items())
        or "No explicit remapping"
    )
    barcode_text = plate_mapping.get("primer_barcode_type") or "not supplied"
    selected_block_text = (
        plate_mapping.get("selected_source_plate_block") or "not selected"
    )
    return f"""# MAXTIA Kinetic–NGS Integrated QC Report

## Key results

- LSA kinetic rows: {summary['input_shape']['lsa_kinetic_rows']}
- NGS rows: {summary['input_shape']['ngs_rows']}
- LSA→NGS one-to-one joins: {join['joined_lsa_rows']} / {summary['input_shape']['lsa_kinetic_rows']}
- Unmatched LSA rows: {join['unmatched_lsa_rows']}
- Structure-mapping-ready rows (NGS C + valid sequence): {qc['structure_mapping_ready_rows']}
- Unique variants with a valid sequence: {summary['variants']['sequence_valid_unique_variants']}
- Unique variants with at least one PASS measurement: {summary['variants']['unique_variants_with_pass_measurement']}

## NGS plate-block mapping

- Primer Barcode type: {barcode_text}
- Selected NGS source block: {selected_block_text}
- Source-to-LSA plate mapping: {plate_pair_text}
- The original NGS well is preserved in `ngs_source_well`; the remapped join key is stored in `ngs_well`.

## QC status

{status_lines}

### Kinetic-fit QC (independent of genotype and sequence)

{kinetic_status_lines}

- Rows reaching a fitting boundary: {qc['fit_boundary_hit_rows']}
- Saturation adequacy (`Data < 50% Rmax`) available from this snapshot: {qc['saturation_qc_available_from_snapshot']}

## Binding phenotype

{outcome_lines}

- Unique variants with an Inactive observation: {summary['variants']['unique_variants_with_inactive_observation']}
- Variants for which every valid replicate is Inactive: {summary['variants']['consistently_inactive_variants']}
- Variants with mixed Active and Inactive observations: {summary['variants']['mixed_active_inactive_variants']}

## WT candidate

- Source: {wt['source']}
- Length: {wt['length']} aa
- NGS rows supporting the same sequence with a `WT` label: {wt['selected_sequence_support']}
- Distinct sequences carrying a WT label: {wt['distinct_wt_labeled_sequences']}
- Candidate sequence: `{wt['sequence']}`

### Internal WT control (technical PASS records only)

- Status: {wt_control['status']}
- PASS WT count: {wt_control['n_pass_WT']}
- Median KD: {fmt(wt_kd['median'])} M
- KD minimum–maximum: {fmt(wt_kd['min'])}–{fmt(wt_kd['max'])} M
- KD maximum/minimum fold: {fmt(wt_kd['max_min_fold'], '.3g')}×

{wt_plate_lines}

The plate grouping above is a proxy inferred from the well key, not actual batch metadata. Because no predefined WT acceptance range is available, it is used only for descriptive statistics.

## Row counts by issue code

{issue_lines}

### Kinetic-fit issue codes

{kinetic_issue_lines}

## Interpretation notes

- The NGS file has no header. Numeric columns 4 and 5 are `inter_cluster_mixing_score` and `intra_cluster_mixing_score`, respectively. The previous names `metric_1` and `metric_2` remain only as backward-compatible aliases.
- NGS column 3 is the total valid read count (`total_valid_read_count`). The provisional Clone rule is inter-cluster score <0.1 and intra-cluster score <0.2 after 98% identity clustering. A 0/0 score is not interpreted as PASS when the total valid read count is zero. Because boundary handling is not yet confirmed, the final C/M call is authoritative and boundary cases remain REVIEW.
- Only NGS code `C` is accepted. Code `M` means that the within-cluster homology criterion was not met, so the genotype is not attributed to a residue.
- Under the Carterra manual definition, `Complex` is an automatic poor-fit indication for residual SD greater than 10% of fitted Rmax. It is excluded as FAIL from quantitative kinetics, while the binding observation itself is retained as qualitative evidence.
- The confirmed Carterra Kinetics 1.9.2.4463 settings use a 10 RU low/no-binding Data Filter and 20 RU response Highlighting. The software `Inactive` flag is preserved, and these thresholds are not treated as interchangeable.
- An `Inactive` record with NGS `M` cannot be assigned to a genotype and is excluded from variant- and residue-level analysis.
- The manual does not define the exported `Error` statistic, so it is not assumed to be a standard error. The default policy marks `Error/value` ≥50% for REVIEW only and does not use it for automatic FAIL.
- The manual's propagation-of-errors equation is applied to replicate means and SDs. It is not used to recalculate the `KD Error` of an individual fit.
- Quantitative kinetics are assigned FAIL when kon, koff, or Rmax reaches a configured fitting boundary, or when Error is zero; raw values are preserved.
- The 10 RU Data Filter and 20 RU Highlighting thresholds are based on observed response, not fitted Rmax. The current export cannot independently reproduce response-based decisions.
- `Data <50% Rmax`, `kd Below Limit`, and `km Below Limit` highlighting were enabled, but their individual states are not exposed in the text columns of the `Kinetic Snapshot`. No manual Excluded Group was used.
- When Active and Inactive observations are mixed and at least one PASS Active record exists, Active records provide the primary quantitative evidence. The overall variant remains REVIEW and Inactive raw records are preserved, but they are excluded from loss-of-binding and hotspot evidence.
- A variant for which every valid replicate is Inactive is described only as `consistent no-detectable response`; it is not called a definitive binding knockout until expression, folding, and spot failure have been excluded.
- The `Kinetic Snapshot` alone cannot reconstruct the insufficient-saturation flag for observed binding below 50% of fitted Rmax. Raw sensorgrams or a richer export are required.
- Variants whose sequence length differs from the WT candidate, or that lack a substitution annotation, are excluded from the initial residue-substitution mapping. Indel-aware alignment should be implemented separately in a later stage.
"""


def render_structure_report(summary: dict[str, Any]) -> str:
    def fmt(value: Any, specification: str = ".3g") -> str:
        return format(float(value), specification) if value is not None else "n.a."

    structure = summary["structure"]
    alignment = structure["sequence_alignment"]
    plddt = structure["plddt_ca"]
    confidence = structure.get("af3_summary_confidences", {})
    analysis_name = summary.get("active_config", {}).get(
        "analysis_name", "Kinetic-to-Structure Analysis"
    )
    position_lines = []
    for item in structure.get("mapped_single_substitution_positions", []):
        best = (
            f"{item['best_KD_mutation']} ({float(item['best_KD_fold_vs_WT']):.3g}×)"
            if item.get("best_KD_mutation")
            else "n.a."
        )
        worst = (
            f"{item['worst_KD_mutation']} ({float(item['worst_KD_fold_vs_WT']):.3g}×)"
            if item.get("worst_KD_mutation")
            else "n.a."
        )
        position_lines.append(
            f"| {item['wt_aa']}{item['sequence_position_1based']} | "
            f"{item['single_substitution_variants']} | {item['variants_with_pass_measurement']} | "
            f"{item['variants_with_any_inactive']} | "
            f"{float(item['position_inactive_fraction_for_loss_evidence'] or 0):.3f} | "
            f"{fmt(item.get('plddt_ca'), '.2f')} | "
            f"{fmt(item.get('residue_sasa_A2'), '.1f')} | {best} | {worst} |"
        )
    position_table = "\n".join(position_lines) or (
        "| n.a. | 0 | 0 | 0 | 0.000 | n.a. | n.a. | n.a. | n.a. |"
    )
    exact_one_to_one = (
        alignment["mismatches"] == 0
        and alignment["gaps_in_reference"] == 0
        and alignment["gaps_in_target"] == 0
        and structure["exact_positions"] == structure["experimental_sequence_length"]
    )
    if exact_one_to_one:
        mapping_note = (
            f"Chain {structure['chain_id']} maps one-to-one to experimental positions "
            f"1–{structure['experimental_sequence_length']}. Both `structure_auth_seq_id` "
            "and `structure_label_seq_id` in `sequence_structure_mapping.csv` record "
            "the structure numbering."
        )
    else:
        mapping_note = (
            "The alignment contains mismatches or gaps. Use "
            "`sequence_structure_mapping.csv` rather than assuming a numbering offset."
        )
    positions = structure.get("mapped_single_substitution_positions", [])
    measured = [
        item
        for item in positions
        if item.get("best_KD_fold_vs_WT") is not None
        or item.get("worst_KD_fold_vs_WT") is not None
    ]
    highlight_lines: list[str] = []
    substitution_rich = sorted(
        (item for item in measured if item["single_substitution_variants"] >= 2),
        key=lambda item: item["single_substitution_variants"],
        reverse=True,
    )
    if substitution_rich:
        labels = ", ".join(
            f"{item['wt_aa']}{item['sequence_position_1based']} "
            f"({item['single_substitution_variants']} substitutions)"
            for item in substitution_rich[:5]
        )
        highlight_lines.append(f"- Substitution-rich positions: {labels}.")
    best_candidates = [
        item for item in measured if item.get("best_KD_fold_vs_WT") is not None
    ]
    if best_candidates:
        best = min(best_candidates, key=lambda item: float(item["best_KD_fold_vs_WT"]))
        highlight_lines.append(
            f"- Lowest observed PASS-summary KD fold: {best['best_KD_mutation']} "
            f"({float(best['best_KD_fold_vs_WT']):.3g}× WT) at "
            f"pLDDT {fmt(best.get('plddt_ca'), '.2f')} and residue SASA "
            f"{fmt(best.get('residue_sasa_A2'), '.1f')} Å²."
        )
    worst_candidates = [
        item for item in measured if item.get("worst_KD_fold_vs_WT") is not None
    ]
    if worst_candidates:
        worst = max(
            worst_candidates, key=lambda item: float(item["worst_KD_fold_vs_WT"])
        )
        highlight_lines.append(
            f"- Highest observed PASS-summary KD fold: {worst['worst_KD_mutation']} "
            f"({float(worst['worst_KD_fold_vs_WT']):.3g}× WT) at "
            f"pLDDT {fmt(worst.get('plddt_ca'), '.2f')} and residue SASA "
            f"{fmt(worst.get('residue_sasa_A2'), '.1f')} Å²."
        )
    mixed = [item for item in positions if item["mixed_active_inactive_variants"]]
    if mixed:
        labels = ", ".join(
            f"{item['wt_aa']}{item['sequence_position_1based']}" for item in mixed
        )
        highlight_lines.append(
            f"- Positions with mixed Active/Inactive variant evidence: {labels}. "
            "These remain REVIEW and are not definitive loss-of-binding evidence."
        )
    if not highlight_lines:
        highlight_lines.append("- No QC-supported single-substitution highlights are available.")
    highlights = "\n".join(highlight_lines)
    multi_position_clusters = [
        cluster
        for cluster in structure.get("candidate_spatial_clusters", [])
        if cluster.get("is_multi_position_cluster")
    ]
    if multi_position_clusters:
        cluster_lines = []
        for cluster in multi_position_clusters:
            labels = ", ".join(cluster["member_labels"])
            maximum_fold = cluster.get("maximum_observed_KD_fold_vs_WT")
            maximum_fold_text = (
                f"{float(maximum_fold):.3g}×"
                if maximum_fold is not None
                else "n.a."
            )
            cluster_lines.append(
                f"- {cluster['cluster_id']}: {labels}; "
                f"{cluster['qc_supported_strong_loss_positions']} QC-supported strong-loss "
                f"position(s), {cluster['inactive_without_pass_kd_review_positions']} "
                "Inactive-without-PASS-KD review position(s); "
                f"maximum observed KD fold {maximum_fold_text}. Members form a connected "
                f"Cα graph at ≤{float(cluster['ca_connectivity_cutoff_A']):.1f} Å."
            )
        cluster_text = "\n".join(cluster_lines)
    else:
        cluster_text = "- No multi-position candidate loss cluster met the configured criteria."
    return f"""# {analysis_name}: Structure-Mapping Validation

## Sequence and numbering

- Chain: {structure['chain_id']}
- Experimental WT sequence length: {structure['experimental_sequence_length']} aa
- Structure polymer sequence length: {structure['structure_polymer_sequence_length']} aa
- Residues observed in the structure: {structure['structure_observed_residue_count']}
- Alignment matches: {alignment['matches']} aa
- Mismatches: {alignment['mismatches']}
- Gaps in the experimental sequence: {alignment['gaps_in_reference']}
- Gaps in the structure sequence: {alignment['gaps_in_target']}
- Experimental residues mapped to the structure: {structure['mapped_positions']}

{mapping_note}

## AlphaFold confidence

- pTM: {confidence.get('ptm', 'n.a.')}
- ranking score: {confidence.get('ranking_score', 'n.a.')}
- clash: {confidence.get('has_clash', 'n.a.')}
- fraction disordered: {confidence.get('fraction_disordered', 'n.a.')}
- Mean Cα pLDDT: {fmt(plddt.get('mean'), '.2f')}
- Median Cα pLDDT: {fmt(plddt.get('median'), '.2f')}
- Cα pLDDT range: {fmt(plddt.get('min'), '.2f')}–{fmt(plddt.get('max'), '.2f')}
- Residue counts at pLDDT <70 / 70–<90 / ≥90: {structure['plddt_ca_counts']['lt_70']} / {structure['plddt_ca_counts']['70_to_lt_90']} / {structure['plddt_ca_counts']['ge_90']}

## Structural features

- Residue and side-chain SASA: Shrake–Rupley, probe {structure['sasa_method']['probe_radius_A']} Å, {structure['sasa_method']['sphere_points_per_atom']} surface points per atom
- Structural contacts: Cα 8 Å and heavy-atom 5 Å, excluding sequence neighbors within ±2 positions
- PAE: chain mean, local ±5-residue mean, and chain maximum
- AF3 predicted contact probability: mean over nonadjacent residues and count of contacts ≥0.8

`structure_mapped_residue_effects.csv` combines the kon, koff, and KD effect of each kinetic substitution with these structural features in one row. Multi-mutants are marked at each affected residue but remain `multi_substitution_confounded` rather than being assigned as single-residue effects.

## Single-substitution position summary

| WT position | Unique substitutions | Substitutions with PASS KD | Variants with Inactive evidence | Loss-evidence Inactive fraction | Cα pLDDT | Residue SASA (Å²) | Lowest KD fold | Highest KD fold |
|---|---:|---:|---:|---:|---:|---:|---|---|
{position_table}

## Automated evidence highlights

{highlights}

## Candidate spatial loss clusters

{cluster_text}

Inactive-without-PASS-KD members remain REVIEW because expression, folding, and spot-validity failure have not been excluded. Connected-component membership does not mean that every pair of residues is within the distance cutoff.

These highlights are deterministic summaries, not interface assignments. A single-entity structure can support local-environment and spatial-clustering analysis but cannot establish direct contact with a binding partner.
"""


def format_affinity(value_molar: Any) -> str:
    value = to_float(value_molar)
    if value is None:
        return "n.a."
    if value < 1e-6:
        return f"{value * 1e9:.3g} nM"
    if value < 1e-3:
        return f"{value * 1e6:.3g} µM"
    return f"{value:.3g} M"


def ranked_variant_tables(variants: list[dict[str, Any]], heading_level: int = 2) -> str:
    single_pass = [
        variant
        for variant in variants
        if variant.get("variant_category") == "single_substitution"
        and variant.get("variant_qc_status") == "PASS"
        and to_float(variant.get("KD_fold_vs_WT_pass")) is not None
    ]
    improvements = sorted(
        (
            variant
            for variant in single_pass
            if float(variant["KD_fold_vs_WT_pass"]) < 0.8
        ),
        key=lambda variant: float(variant["KD_fold_vs_WT_pass"]),
    )[:5]
    losses = sorted(
        (
            variant
            for variant in single_pass
            if float(variant["KD_fold_vs_WT_pass"]) > 1.25
        ),
        key=lambda variant: float(variant["KD_fold_vs_WT_pass"]),
        reverse=True,
    )[:5]

    def number(raw: Any, suffix: str = "") -> str:
        value = to_float(raw)
        return f"{value:.3g}{suffix}" if value is not None else "—"

    def variant_table(rows: list[dict[str, Any]]) -> str:
        if not rows:
            return "No QC-supported variants met this criterion."
        lines = [
            "| Variant | Observed ROIs | ROIs used for medians | PASS active n | Median KD | "
            "Median koff (s⁻¹) | Median kon (M⁻¹ s⁻¹) | Median Rmax (RU) | KD fold vs WT | kon fold | koff fold |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in rows:
            observed = (row.get("observed_rois") or "—").replace("|", "\\|")
            used = (row.get("median_rois") or "—").replace("|", "\\|")
            lines.append(
                f"| {row['mutation_labels']} | {observed} | {used} | {row['n_active_quantitative_pass']} | "
                f"{format_affinity(row.get('median_KD_M_pass'))} | "
                f"{number(row.get('median_koff_s_inv_pass'))} | {number(row.get('median_kon_M_inv_s_pass'))} | "
                f"{number(row.get('median_Rmax_RU_pass'))} | {number(row.get('KD_fold_vs_WT_pass'), '×')} | "
                f"{number(row.get('kon_fold_vs_WT_pass'), '×')} | {number(row.get('koff_fold_vs_WT_pass'), '×')} |"
            )
        return "\n".join(lines)

    heading = "#" * heading_level
    return "\n\n".join([
        f"{heading} Strongest QC-supported single-substitution improvements",
        "Observed ROIs identify sequence-confirmed occurrences in this batch. "
        "The median columns use the indicated active measurements that passed kinetic QC.",
        variant_table(improvements),
        f"{heading} Strongest QC-supported single-substitution losses",
        variant_table(losses),
    ])


def render_analysis_report(
    summary: dict[str, Any],
    variants: list[dict[str, Any]],
    position_summary: list[dict[str, Any]],
    config: dict[str, Any],
    residue_colors: list[dict[str, Any]] | None = None,
) -> str:
    analysis_name = config.get("analysis_name", "MAXTIA Kinetic-to-Structure Analysis")
    wt_summary = summary["wt_internal_control"]["parameter_summary"]["KD_M"]
    mixed_variants = [
        variant
        for variant in variants
        if variant.get("phenotype_consensus") == "mixed_active_inactive"
    ]
    mixed_lines = (
        "\n".join(
            f"- {variant['mutation_labels']}: {variant['n_active_quantitative_pass']} "
            f"PASS Active and {variant['n_inactive']} Inactive; overall "
            f"{variant['variant_qc_status']}."
            for variant in mixed_variants
        )
        or "- None."
    )
    structure = summary.get("structure")
    if structure:
        alignment = structure["sequence_alignment"]
        mean_plddt = structure["plddt_ca"].get("mean")
        multi_position_clusters = [
            cluster
            for cluster in structure.get("candidate_spatial_clusters", [])
            if cluster.get("is_multi_position_cluster")
        ]
        cluster_labels = (
            "; ".join(
                f"{cluster['cluster_id']} ({', '.join(cluster['member_labels'])})"
                for cluster in multi_position_clusters
            )
            or "none"
        )
        structure_lines = f"""- Chain: {structure['chain_id']}
- Sequence matches/mismatches/gaps: {alignment['matches']}/{alignment['mismatches']}/{alignment['gaps_in_reference'] + alignment['gaps_in_target']}
- Mapped experimental positions: {structure['mapped_positions']} / {structure['experimental_sequence_length']}
- Mean Cα pLDDT: {format(float(mean_plddt), '.2f') if mean_plddt is not None else 'n.a.'}
- Single-substitution positions summarized: {len(position_summary)}
- Candidate multi-position spatial loss clusters: {cluster_labels}"""
    else:
        structure_lines = "- No structure input was supplied."
    eligibility = summary["qc"]["analysis_eligibility_counts"]
    wt_fold = wt_summary.get("max_min_fold")
    wt_fold_text = f"{float(wt_fold):.3g}×" if wt_fold is not None else "n.a."
    return f"""# {analysis_name}: Deterministic Analysis Report

This report was generated by the command-line module without an LLM.

## Data integrity and QC

- LSA rows: {summary['input_shape']['lsa_kinetic_rows']}
- NGS rows: {summary['input_shape']['ngs_rows']}
- One-to-one LSA–NGS joins: {summary['key_join']['joined_lsa_rows']}
- Unmatched LSA rows: {summary['key_join']['unmatched_lsa_rows']}
- Quantitative rows: {eligibility.get('QUANTITATIVE', 0)}
- Review rows: {eligibility.get('REVIEW', 0)}
- Excluded uncertain genotype or sequence: {eligibility.get('EXCLUDE_UNCERTAIN_GENOTYPE_OR_SEQUENCE', 0)}
- Qualitative-only Complex rows: {eligibility.get('QUALITATIVE_BINDING_ONLY_COMPLEX', 0)}
- Censored Inactive rows: {eligibility.get('CENSORED_INACTIVE', 0)}

## WT baseline

- PASS WT records: {summary['wt_internal_control']['n_pass_WT']}
- Median WT KD: {format_affinity(wt_summary.get('median'))}
- WT KD maximum/minimum fold: {wt_fold_text}

{ranked_variant_tables(variants)}

## Mixed Active/Inactive evidence

{mixed_lines}

## Structure validation

{structure_lines}

## Residues by structure classification

{report_residue_classes(residue_colors)}

## Interpretation limits

- Variant rankings use kinetic-fit PASS Active records only.
- Inactive records are never converted to an artificial KD value.
- Mixed Active/Inactive variants remain REVIEW.
- The Kinetic Snapshot cannot independently reconstruct the `Data <50% Rmax` saturation flag.
- Without expression, immobilization, or spot-validity controls, Inactive observations are not definitive binding knockouts.
- A single-entity structure does not establish a binding interface.
"""


def write_run_manifest(
    out_dir: Path, summary: dict[str, Any], config: dict[str, Any]
) -> None:
    output_files: dict[str, Any] = {}
    for path in sorted(out_dir.iterdir()):
        if path.is_file() and path.name != "run_manifest.json":
            output_files[path.name] = {
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
    config_bytes = json.dumps(
        config, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    manifest = {
        "module_version": MODULE_VERSION,
        "qc_policy_name": config.get("qc_policy_name"),
        "analysis_name": config.get("analysis_name"),
        "runtime": {
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "executable_name": Path(sys.executable).name,
        },
        "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "source_files": summary.get("source_files", {}),
        "structure_source_files": summary.get("structure", {}).get(
            "source_files", {}
        ),
        "output_files": output_files,
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def write_structure_views(job, scale):
    folder, mapping, variants, positions, config = (job[k] for k in ("out_dir", "mapping", "variants", "positions", "config"))
    apply_structure_scale(positions, scale)
    apply_structure_scale(variants, scale)
    write_csv(folder / "structure_visualization_variants.csv", variants)
    write_csv(folder / "structure_visualization_positions.csv", positions)
    write_csv(folder / "residue_colors.csv", all_residue_colors(mapping, positions))
    (folder / "structure_color_scale.json").write_text(json.dumps(scale, indent=2))
    (folder / "visualize_overview.pml").write_text(render_pymol_overview_script(positions, config))
    (folder / "visualize_variant_scenes.pml").write_text(render_pymol_variant_scenes_script(variants, config))
    (folder / "visualize_overview.cxc").write_text(render_chimerax_overview_script(positions, config))
    (folder / "structure_overview.svg").write_text(render_dependency_free_structure_svg(mapping, positions, config, scale))
    (folder / "VISUALIZATION.md").write_text(render_visualization_readme(config))
    return dict(color_scale=scale, dependency_free_svg=dict(status="GENERATED", file="structure_overview.svg",
                representation="three orthographic C-alpha projections"),
                static_rendering=render_static_images_with_pymol(folder, variants, config))


def run(
    lsa_path: Path,
    ngs_path: Path,
    out_dir: Path,
    config: dict[str, Any],
    structure_cif: Path | None = None,
    full_data_json: Path | None = None,
    summary_confidences_json: Path | None = None,
    structure_chain_id: str | None = None,
    pending_views: list | None = None,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    lsa_records, lsa_warnings = parse_lsa(lsa_path, config)
    ngs_records, ngs_warnings = parse_ngs(ngs_path, config)
    validate_ngs_plate_scope(lsa_records, ngs_records, config)
    wt_sequence, wt_summary = infer_wt_sequence(ngs_records, config.get("wt_sequence", ""))
    integrated = join_and_qc(lsa_records, ngs_records, wt_sequence, config)
    variants = summarize_variants(integrated, wt_sequence, config)
    residues = residue_effects(variants)
    examples = select_qc_examples(integrated)
    summary = make_summary(
        lsa_path,
        ngs_path,
        lsa_records,
        ngs_records,
        integrated,
        variants,
        wt_summary,
        lsa_warnings + ngs_warnings,
        config,
    )
    mapping_rows: list[dict[str, Any]] = []
    mapped_effects: list[dict[str, Any]] = []
    position_summary: list[dict[str, Any]] = []
    candidate_cluster_members: list[dict[str, Any]] = []
    visualization_variants: list[dict[str, Any]] = []
    visualization_positions: list[dict[str, Any]] = []
    report_colors = None
    if structure_cif:
        mapping_rows, structure_summary = build_structure_mapping(
            wt_sequence,
            structure_cif,
            full_data_json,
            summary_confidences_json,
            structure_chain_id or str(config["structure_chain_id"]),
            config,
        )
        mapped_effects = merge_residue_effects_with_structure(residues, mapping_rows)
        position_summary = summarize_mapped_positions(mapped_effects)
        candidate_cluster_members, candidate_cluster_summary = (
            build_candidate_spatial_clusters(position_summary, mapping_rows, config)
        )
        if config.get("structure_visualization", {}).get("enabled", True):
            visualization_variants = build_structure_visualization_records(
                [e for e in mapped_effects if e.get("alignment_status") == "exact"], config
            )
            if pending_views is not None and not any(v["aa_sequence"] == wt_sequence and v["variant_qc_status"] == "PASS" for v in variants):
                for record in visualization_variants:
                    if record["visualization_class"] == "QUANTITATIVE":
                        record.update(structure_fields("NO_PASS_DATA"))
            visualization_positions = build_structure_visualization_positions(
                visualization_variants, config
            )
            include_confounded_positions(visualization_positions, mapped_effects)
            report_colors = all_residue_colors(mapping_rows, visualization_positions)
            structure_summary["visualization"] = {
                "variant_record_count": len(visualization_variants),
                "position_record_count": len(visualization_positions),
                "variant_class_counts": counts(
                    row["visualization_class"] for row in visualization_variants
                ),
                "position_class_counts": counts(
                    row["position_visualization_class"]
                    for row in visualization_positions
                ),
            }
        structure_summary["mapped_single_substitution_positions"] = position_summary
        structure_summary["candidate_spatial_clusters"] = candidate_cluster_summary
        summary["structure"] = structure_summary
    write_csv(out_dir / "integrated_records.csv", integrated)
    write_csv(out_dir / "variant_summary.csv", variants)
    write_csv(out_dir / "residue_effects.csv", residues)
    write_csv(out_dir / "qc_examples.csv", examples)
    if structure_cif:
        write_csv(out_dir / "sequence_structure_mapping.csv", mapping_rows)
        write_csv(out_dir / "structure_mapped_residue_effects.csv", mapped_effects)
        write_csv(out_dir / "position_phenotype_summary.csv", position_summary)
        write_csv(
            out_dir / "structure_candidate_cluster_members.csv",
            candidate_cluster_members,
        )
        if "visualization" in summary["structure"]:
            shutil.copy2(structure_cif, out_dir / "structure_model.cif")
            job = dict(out_dir=out_dir, mapping=mapping_rows, variants=visualization_variants,
                       positions=visualization_positions, config=config)
            if pending_views is None:
                scale = structure_scale(visualization_positions, visualization_variants)
                summary["structure"]["visualization"].update(write_structure_views(job, scale))
            else:
                pending_views.append(job)
    (out_dir / "qc_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (out_dir / "qc_report.md").write_text(render_qc_report(summary), encoding="utf-8")
    (out_dir / "analysis_report.md").write_text(
        render_analysis_report(summary, variants, position_summary, config, report_colors),
        encoding="utf-8",
    )
    if structure_cif:
        (out_dir / "structure_mapping_report.md").write_text(
            render_structure_report(summary), encoding="utf-8"
        )
    write_run_manifest(out_dir, summary, config)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Integrate Carterra LSA kinetics with MAXTIA NGS library calls"
    )
    parser.add_argument("--lsa-xlsx", required=True, type=Path)
    parser.add_argument("--ngs-csv", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--structure-cif", type=Path)
    parser.add_argument("--full-data-json", type=Path)
    parser.add_argument("--summary-confidences-json", type=Path)
    parser.add_argument("--structure-chain")
    parser.add_argument(
        "--ngs-plate-block",
        help="Primer-Barcode-dependent NGS source block, for example P05-P08",
    )
    parser.add_argument(
        "--ngs-target-plate-start",
        type=int,
        help="First LSA plate number to which the selected NGS block is mapped",
    )
    parser.add_argument("--ngs-primer-barcode-type")
    args = parser.parse_args()
    config = dict(DEFAULT_CONFIG)
    if args.config:
        config = deep_merge(
            config, json.loads(args.config.read_text(encoding="utf-8"))
        )
    if args.ngs_plate_block is not None:
        config["ngs_plate_block"] = args.ngs_plate_block
    if args.ngs_target_plate_start is not None:
        config["ngs_target_plate_start"] = args.ngs_target_plate_start
    if args.ngs_primer_barcode_type is not None:
        config["ngs_primer_barcode_type"] = args.ngs_primer_barcode_type
    summary = run(
        args.lsa_xlsx,
        args.ngs_csv,
        args.out_dir,
        config,
        structure_cif=args.structure_cif,
        full_data_json=args.full_data_json,
        summary_confidences_json=args.summary_confidences_json,
        structure_chain_id=args.structure_chain,
    )
    print(json.dumps(summary["key_join"], ensure_ascii=False))
    print(json.dumps(summary["qc"]["status_counts"], ensure_ascii=False))
    if "structure" in summary:
        print(json.dumps(summary["structure"]["sequence_alignment"], ensure_ascii=False))


if __name__ == "__main__":
    main()
