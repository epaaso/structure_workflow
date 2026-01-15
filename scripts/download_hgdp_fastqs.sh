#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  download_hgdp_fastqs.sh --samples <samples.txt> --outdir <dir> [--download] [--jobs N]

What it does:
  - For each HGDP sample ID (e.g. HGDP00455), queries ENA for read runs with sample_alias="HGDP00455".
  - Extracts the ENA-provided fastq_ftp field (one or more FASTQ URLs).
  - Writes a deduplicated URL list and an optional download manifest.
  - If --download is provided, downloads all FASTQs with resume support.

Notes:
  - This downloads raw FASTQ and can be very large.
  - Requires: curl, awk, sort, xargs, gzip.

Examples:
  # list URLs only
  ./download_hgdp_fastqs.sh --samples inputs/biaka_hgdp.samples.txt --outdir /datos/migccl/ancestry_refs/BiakaHGDP

  # download with 4 parallel transfers
  ./download_hgdp_fastqs.sh --samples inputs/biaka_hgdp.samples.txt --outdir /datos/migccl/ancestry_refs/BiakaHGDP --download --jobs 4
EOF
}

SAMPLES=""
OUTDIR=""
DO_DOWNLOAD=0
JOBS=4

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples)
      SAMPLES="$2"; shift 2;;
    --outdir)
      OUTDIR="$2"; shift 2;;
    --download)
      DO_DOWNLOAD=1; shift;;
    --jobs)
      JOBS="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "$SAMPLES" || -z "$OUTDIR" ]]; then
  usage
  exit 2
fi

if [[ ! -f "$SAMPLES" ]]; then
  echo "Samples file not found: $SAMPLES" >&2
  exit 2
fi

mkdir -p "$OUTDIR" "$OUTDIR/fastq" "$OUTDIR/meta"

URLS_TMP="$OUTDIR/meta/fastq_urls.tmp.txt"
URLS_OUT="$OUTDIR/meta/fastq_urls.txt"
MANIFEST_OUT="$OUTDIR/meta/fastq_manifest.tsv"
: > "$URLS_TMP"
: > "$MANIFEST_OUT"

ena_query() {
  local sample="$1"
  # ENA API: TSV with fields (run_accession, fastq_ftp). fastq_ftp is semi-colon separated.
  # limit=0 means no limit.
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=sample_alias%3D%22${sample}%22&fields=run_accession,fastq_ftp&format=tsv&limit=0"
}

while IFS= read -r sample; do
  sample="${sample//$'\r'/}"
  [[ -z "$sample" ]] && continue

  echo "[ENA] ${sample}" >&2

  # Output TSV has header. Extract run_accession + fastq_ftp.
  # Some samples can have multiple runs/lanes.
  ena_query "$sample" \
    | awk -F'\t' 'NR==1{next} $1!=""{print $1"\t"$2}' \
    | while IFS=$'\t' read -r run fastq_ftp; do
        [[ -z "$run" ]] && continue
        [[ -z "$fastq_ftp" ]] && continue

        # fastq_ftp: host/path1.gz;host/path2.gz (no scheme). Convert to ftp://...
        IFS=';' read -ra parts <<< "$fastq_ftp"
        for p in "${parts[@]}"; do
          [[ -z "$p" ]] && continue
          url="ftp://${p}"
          echo "$url" >> "$URLS_TMP"
          printf '%s\t%s\t%s\n' "$sample" "$run" "$url" >> "$MANIFEST_OUT"
        done
      done

done < "$SAMPLES"

# Deduplicate URLs
sort -u "$URLS_TMP" > "$URLS_OUT"
rm -f "$URLS_TMP"

echo "Wrote URL list: $URLS_OUT" >&2
echo "Wrote manifest: $MANIFEST_OUT" >&2
echo "FASTQ URLs: $(wc -l < "$URLS_OUT")" >&2

if [[ "$DO_DOWNLOAD" -eq 0 ]]; then
  echo "(dry-run) Not downloading. Re-run with --download to fetch FASTQs." >&2
  exit 0
fi

# Download all URLs into OUTDIR/fastq (keep basename)
# Use curl resume (-C -) and follow redirects (-L). Parallelize with xargs.
cat "$URLS_OUT" \
  | xargs -P "$JOBS" -I {} bash -lc 'u="$1"; outdir="$2"; f="$(basename "$u")"; dest="$outdir/fastq/$f"; if [[ -f "$dest" ]]; then if gzip -t "$dest" >/dev/null 2>&1; then echo "[SKIP] $f (ok)" >&2; exit 0; else echo "[RESUME] $f (corrupt or incomplete)" >&2; if ! curl -fL -C - -o "$dest" "$u"; then echo "[REDOWNLOAD] $f (resume failed)" >&2; rm -f "$dest"; fi; if [[ -f "$dest" ]] && gzip -t "$dest" >/dev/null 2>&1; then echo "[OK] $f" >&2; exit 0; else echo "[REDOWNLOAD] $f (still corrupt)" >&2; rm -f "$dest"; fi; fi; fi; echo "[DL] $f" >&2; curl -fL -o "$dest" "$u"' _ {} "$OUTDIR"

echo "Done. Files in: $OUTDIR/fastq" >&2
