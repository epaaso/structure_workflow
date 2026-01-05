#!/usr/bin/env python3

import argparse
import csv
import glob
import gzip
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple


GT_SPLIT_RE = re.compile(r"[\/|]")


def normalize_chrom(chrom: str) -> str:
    c = chrom.strip()
    c_lower = c.lower()
    if c_lower.startswith("chr"):
        c = c[3:]
    return c


def detect_field(header: List[str], wanted: List[str]) -> Optional[str]:
    normalized_to_original = {h.strip(): h for h in header}
    for w in wanted:
        for key, original in normalized_to_original.items():
            if key.lower() == w.lower():
                return original
    return None


@dataclass
class Target:
    rsid: str
    chrom: str
    pos: int
    a1: str
    a2: str


@dataclass
class TargetStats:
    rsid: str
    chrom: str
    pos: int
    a1: str
    a2: str
    vcf_id: str = ""
    vcf_ref: str = ""
    vcf_alt: str = ""
    n_records_seen: int = 0
    n_samples_with_record: int = 0
    n_samples_called_gt: int = 0
    n_samples_missing_gt: int = 0
    n_samples_with_a1: int = 0
    n_samples_with_a2: int = 0
    notes: Set[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.notes is None:
            self.notes = set()


def read_aim_csv(path: str) -> List[Target]:
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in AIM file: {path}")

        rsid_field = detect_field(reader.fieldnames, ["SNP rsID", "rsID", "rsid", "SNP"])
        chr_field = detect_field(reader.fieldnames, ["chr", "chrom", "chromosome"])
        pos_field = detect_field(reader.fieldnames, ["position", "pos", "bp"])
        a1_field = detect_field(reader.fieldnames, ["A1", "allele1"])
        a2_field = detect_field(reader.fieldnames, ["A2", "allele2"])

        missing = [
            name
            for name, field in [
                ("rsid", rsid_field),
                ("chr", chr_field),
                ("position", pos_field),
                ("A1", a1_field),
                ("A2", a2_field),
            ]
            if field is None
        ]
        if missing:
            raise ValueError(
                "AIM CSV is missing required columns; need at least: SNP rsID/rsid, chr, position, A1, A2. "
                f"Missing: {', '.join(missing)}. Found: {reader.fieldnames}"
            )

        targets: List[Target] = []
        for row in reader:
            rsid = (row.get(rsid_field) or "").strip()
            chrom = normalize_chrom((row.get(chr_field) or "").strip())
            pos_raw = (row.get(pos_field) or "").strip()
            a1 = (row.get(a1_field) or "").strip().upper()
            a2 = (row.get(a2_field) or "").strip().upper()

            if not chrom or not pos_raw:
                continue
            try:
                pos = int(float(pos_raw))
            except ValueError:
                continue

            if rsid in ("", "."):
                rsid = f"{chrom}:{pos}"

            targets.append(Target(rsid=rsid, chrom=chrom, pos=pos, a1=a1, a2=a2))

    # de-duplicate by chrom/pos (keep first)
    seen: Set[Tuple[str, int]] = set()
    deduped: List[Target] = []
    for t in targets:
        key = (t.chrom, t.pos)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(t)
    return deduped


def iter_vcf_paths(vcf_globs: List[str]) -> List[str]:
    paths: List[str] = []
    for pattern in vcf_globs:
        matches = glob.glob(pattern)
        paths.extend(matches)
    # stable order
    paths = sorted(set(paths))
    return paths


def parse_gt_allele_indexes(gt: str) -> Optional[List[int]]:
    gt = gt.strip()
    if gt in (".", "./.", ".|."):
        return None
    parts = GT_SPLIT_RE.split(gt)
    idxs: List[int] = []
    for p in parts:
        if p == "." or p == "":
            return None
        try:
            idxs.append(int(p))
        except ValueError:
            return None
    return idxs


def scan_vcf(
    vcf_path: str,
    targets: Dict[Tuple[str, int], TargetStats],
    remaining_keys: Set[Tuple[str, int]],
) -> None:
    # Stream the bgzipped/gz VCF; no .tbi required.
    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as f:
        sample_names: List[str] = []
        gt_index: Optional[int] = None
        format_keys: List[str] = []

        for line in f:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_names = cols[9:]
                continue
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom_raw = parts[0]
            pos = int(parts[1])
            vid = parts[2] if len(parts) > 2 else "."
            ref = parts[3] if len(parts) > 3 else ""
            alt = parts[4] if len(parts) > 4 else ""
            fmt = parts[8] if len(parts) > 8 else ""
            samples = parts[9:] if len(parts) > 9 else []

            chrom = normalize_chrom(chrom_raw)
            key = (chrom, pos)
            if key not in targets:
                continue

            st = targets[key]
            st.n_records_seen += 1
            st.vcf_id = vid
            st.vcf_ref = ref
            st.vcf_alt = alt

            alts = alt.split(",") if alt else []
            allele_bases: List[str] = [ref] + alts
            if any(len(a) != 1 for a in allele_bases if a):
                st.notes.add("non_snv")

            if st.a1 and st.a1 not in allele_bases:
                st.notes.add("A1_not_in_REF_ALT")
            if st.a2 and st.a2 not in allele_bases:
                st.notes.add("A2_not_in_REF_ALT")

            format_keys = fmt.split(":") if fmt else []
            gt_index = None
            for i, k in enumerate(format_keys):
                if k == "GT":
                    gt_index = i
                    break
            if gt_index is None:
                st.notes.add("no_GT")

            # If no samples in this VCF, still consider record found.
            if not samples or not sample_names:
                remaining_keys.discard(key)
                continue

            # Some VCFs can be single-sample with missing #CHROM sample name(s);
            # fall back to counting sample columns.
            if not sample_names:
                sample_names = [f"sample{i+1}" for i in range(len(samples))]

            for sample_name, sample_field in zip(sample_names, samples):
                st.n_samples_with_record += 1
                if gt_index is None:
                    st.n_samples_missing_gt += 1
                    continue

                fields = sample_field.split(":")
                if gt_index >= len(fields):
                    st.n_samples_missing_gt += 1
                    continue

                gt = fields[gt_index]
                idxs = parse_gt_allele_indexes(gt)
                if idxs is None:
                    st.n_samples_missing_gt += 1
                    continue

                st.n_samples_called_gt += 1
                observed: Set[str] = set()
                for allele_idx in idxs:
                    if allele_idx < 0:
                        continue
                    if allele_idx >= len(allele_bases):
                        st.notes.add("GT_out_of_range")
                        continue
                    base = allele_bases[allele_idx].upper()
                    observed.add(base)

                if st.a1 and st.a1.upper() in observed:
                    st.n_samples_with_a1 += 1
                if st.a2 and st.a2.upper() in observed:
                    st.n_samples_with_a2 += 1

            remaining_keys.discard(key)


def write_report(path: str, stats: Iterable[TargetStats]) -> None:
    fieldnames = [
        "rsid",
        "chr",
        "position",
        "A1",
        "A2",
        "found_in_vcf",
        "vcf_id",
        "vcf_ref",
        "vcf_alt",
        "n_records_seen",
        "n_samples_with_record",
        "n_samples_called_gt",
        "n_samples_missing_gt",
        "A1_present",
        "A2_present",
        "n_samples_with_A1",
        "n_samples_with_A2",
        "notes",
    ]

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for st in stats:
            found = st.n_records_seen > 0
            if not found:
                st.notes.add("site_not_in_vcf")
            w.writerow(
                {
                    "rsid": st.rsid,
                    "chr": st.chrom,
                    "position": st.pos,
                    "A1": st.a1,
                    "A2": st.a2,
                    "found_in_vcf": "1" if st.n_records_seen > 0 else "0",
                    "vcf_id": st.vcf_id,
                    "vcf_ref": st.vcf_ref,
                    "vcf_alt": st.vcf_alt,
                    "n_records_seen": st.n_records_seen,
                    "n_samples_with_record": st.n_samples_with_record,
                    "n_samples_called_gt": st.n_samples_called_gt,
                    "n_samples_missing_gt": st.n_samples_missing_gt,
                    # If the site is absent from variant-only VCFs, we cannot distinguish
                    # homozygous reference from not-called; represent presence as unknown.
                    "A1_present": ("1" if st.n_samples_with_a1 > 0 else "0") if found else ".",
                    "A2_present": ("1" if st.n_samples_with_a2 > 0 else "0") if found else ".",
                    "n_samples_with_A1": st.n_samples_with_a1,
                    "n_samples_with_A2": st.n_samples_with_a2,
                    "notes": ";".join(sorted(st.notes)) if st.notes else "",
                }
            )


def main(argv: List[str]) -> int:
    p = argparse.ArgumentParser(
        description=(
            "Given an AIM SNP list (chr/pos/A1/A2), scan one or more bgzipped VCFs and report which AIM alleles "
            "are observed in sample genotypes. Works without .tbi indexes by streaming VCF.gz."
        )
    )
    p.add_argument("--aim", required=True, help="Path to AIM CSV (must include columns: chr, position, A1, A2, rsid)")
    p.add_argument(
        "--vcf",
        required=True,
        action="append",
        help="VCF path or glob. Repeatable. Example: --vcf '/datos/migccl/vcfs/*.vcf.gz'",
    )
    p.add_argument("--out", required=True, help="Output TSV path")

    args = p.parse_args(argv)

    vcf_paths = iter_vcf_paths(args.vcf)
    if not vcf_paths:
        print(f"No VCFs matched: {args.vcf}", file=sys.stderr)
        return 2

    targets_list = read_aim_csv(args.aim)
    if not targets_list:
        print(f"No targets parsed from AIM file: {args.aim}", file=sys.stderr)
        return 2

    targets: Dict[Tuple[str, int], TargetStats] = {}
    for t in targets_list:
        key = (t.chrom, t.pos)
        targets[key] = TargetStats(rsid=t.rsid, chrom=t.chrom, pos=t.pos, a1=t.a1, a2=t.a2)

    remaining_keys: Set[Tuple[str, int]] = set(targets.keys())

    for vcf_path in vcf_paths:
        scan_vcf(vcf_path, targets, remaining_keys)

    ordered_stats = sorted(targets.values(), key=lambda s: (int(s.chrom) if s.chrom.isdigit() else 10**9, s.chrom, s.pos))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    write_report(args.out, ordered_stats)

    n_total = len(ordered_stats)
    n_found = sum(1 for s in ordered_stats if s.n_records_seen > 0)
    n_a1_present = sum(1 for s in ordered_stats if s.n_records_seen > 0 and s.n_samples_with_a1 > 0)
    n_a2_present = sum(1 for s in ordered_stats if s.n_records_seen > 0 and s.n_samples_with_a2 > 0)
    print(f"Targets: {n_total}")
    print(f"Found in VCFs (by chr:pos): {n_found}")
    print(f"A1 present in any sample (among found sites): {n_a1_present}")
    print(f"A2 present in any sample (among found sites): {n_a2_present}")
    print(f"Wrote: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
