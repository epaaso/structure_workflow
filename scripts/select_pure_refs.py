#!/usr/bin/env python3

import argparse
import csv
import os
import sys
from typing import Dict, Iterable, List, Optional, Tuple


def detect_field(fieldnames: Iterable[str], wanted: List[str]) -> Optional[str]:
    normalized = [f.strip() for f in fieldnames if f is not None]
    for w in wanted:
        for f in normalized:
            if f.lower() == w.lower():
                return f
    return None


def read_panel(panel_path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    with open(panel_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in panel file: {panel_path}")
        sample_field = detect_field(reader.fieldnames, ["sample", "iid", "id"])
        pop_field = detect_field(reader.fieldnames, ["pop", "population"])
        super_field = detect_field(reader.fieldnames, ["super_pop", "superpop", "super_population"])
        if not sample_field or not pop_field or not super_field:
            raise ValueError(
                "Panel file is missing required columns (sample/pop/super_pop). "
                f"Found: {reader.fieldnames}"
            )
        pop_map: Dict[str, str] = {}
        super_map: Dict[str, str] = {}
        for row in reader:
            sample = (row.get(sample_field) or "").strip()
            pop = (row.get(pop_field) or "").strip()
            super_pop = (row.get(super_field) or "").strip()
            if not sample:
                continue
            pop_map[sample] = pop
            super_map[sample] = super_pop
    return pop_map, super_map


def read_ancestry(ancestry_path: str) -> Tuple[Dict[str, Dict[str, str]], str]:
    with open(ancestry_path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No header found in ancestry file: {ancestry_path}")
        iid_field = detect_field(reader.fieldnames, ["IID", "iid", "sample", "id"])
        if not iid_field:
            raise ValueError(
                "Ancestry file is missing an IID column. "
                f"Found: {reader.fieldnames}"
            )
        rows: Dict[str, Dict[str, str]] = {}
        for row in reader:
            iid = (row.get(iid_field) or "").strip()
            if not iid:
                continue
            rows[iid] = row
    return rows, iid_field


def select_group(
    label: str,
    ids: List[str],
    ancestry: Dict[str, Dict[str, str]],
    component_field: str,
    min_value: float,
    max_count: int,
) -> Tuple[List[Tuple[str, float]], int, int]:
    selected: List[Tuple[str, float]] = []
    missing = 0
    for iid in ids:
        row = ancestry.get(iid)
        if not row:
            missing += 1
            continue
        if component_field not in row:
            raise ValueError(
                f"Ancestry file missing component column '{component_field}' for {label}."
            )
        try:
            val = float(row[component_field])
        except (TypeError, ValueError):
            continue
        if val >= min_value:
            selected.append((iid, val))
    selected.sort(key=lambda x: x[1], reverse=True)
    if max_count and max_count > 0:
        selected = selected[:max_count]
    return selected, missing, len(ids)


def write_list(path: str, rows: List[Tuple[str, float]]) -> None:
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, "w", newline="") as f:
        for iid, _ in rows:
            f.write(f"{iid}\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Select high-ancestry reference samples from ADMIXTURE/STRUCTURE output. "
            "Outputs one IID per line for each reference group."
        )
    )
    parser.add_argument("--panel", required=True, help="1KG panel file (integrated_call_samples...panel)")
    parser.add_argument("--ancestry", required=True, help="all_ancestry.tsv from the pipeline")
    parser.add_argument("--mxl-pop", default="MXL", help="Population code for MXL refs")
    parser.add_argument("--eur-super", default="EUR", help="Super-population code for EUR refs")
    parser.add_argument("--afr-pop", default="YRI", help="Population code for AFR/YRI refs")
    parser.add_argument("--afr-super", default="AFR", help="Super-population code for AFR refs (used if --afr-pop is empty)")
    parser.add_argument("--mxl-comp", default="MXL", help="Ancestry component column for MXL")
    parser.add_argument("--eur-comp", default="EUR", help="Ancestry component column for EUR")
    parser.add_argument("--afr-comp", default="AFR", help="Ancestry component column for AFR/YRI")
    parser.add_argument("--min-mxl", type=float, default=0.85, help="Minimum MXL component to keep")
    parser.add_argument("--min-eur", type=float, default=0.90, help="Minimum EUR component to keep")
    parser.add_argument("--min-afr", type=float, default=0.90, help="Minimum AFR component to keep")
    parser.add_argument("--max-per-group", type=int, default=0, help="Optional cap per group (0 = no cap)")
    parser.add_argument("--out-mxl", required=True, help="Output list for MXL refs")
    parser.add_argument("--out-eur", required=True, help="Output list for EUR refs")
    parser.add_argument("--out-afr", required=True, help="Output list for AFR/YRI refs")
    args = parser.parse_args()

    pop_map, super_map = read_panel(args.panel)
    ancestry, _ = read_ancestry(args.ancestry)

    mxl_ids = [iid for iid, pop in pop_map.items() if pop == args.mxl_pop]
    eur_ids = [iid for iid, super_pop in super_map.items() if super_pop == args.eur_super]
    afr_ids: List[str]
    if args.afr_pop:
        afr_ids = [iid for iid, pop in pop_map.items() if pop == args.afr_pop]
    else:
        afr_ids = [iid for iid, super_pop in super_map.items() if super_pop == args.afr_super]

    mxl_selected, mxl_missing, mxl_total = select_group(
        "MXL", mxl_ids, ancestry, args.mxl_comp, args.min_mxl, args.max_per_group
    )
    eur_selected, eur_missing, eur_total = select_group(
        "EUR", eur_ids, ancestry, args.eur_comp, args.min_eur, args.max_per_group
    )
    afr_selected, afr_missing, afr_total = select_group(
        "AFR", afr_ids, ancestry, args.afr_comp, args.min_afr, args.max_per_group
    )

    summary = [
        ("MXL", mxl_total, mxl_missing, len(mxl_selected), args.min_mxl),
        ("EUR", eur_total, eur_missing, len(eur_selected), args.min_eur),
        ("AFR", afr_total, afr_missing, len(afr_selected), args.min_afr),
    ]
    for label, total, missing, kept, threshold in summary:
        print(
            f"{label}: total={total}, missing_in_ancestry={missing}, kept={kept}, min={threshold}",
            file=sys.stderr,
        )

    if not mxl_selected:
        raise SystemExit("No MXL samples passed the threshold; lower --min-mxl or set --max-per-group.")
    if not eur_selected:
        raise SystemExit("No EUR samples passed the threshold; lower --min-eur or set --max-per-group.")
    if not afr_selected:
        raise SystemExit("No AFR/YRI samples passed the threshold; lower --min-afr or set --max-per-group.")

    write_list(args.out_mxl, mxl_selected)
    write_list(args.out_eur, eur_selected)
    write_list(args.out_afr, afr_selected)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
