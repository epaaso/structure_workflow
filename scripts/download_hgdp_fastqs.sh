#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  download_hgdp_fastqs.sh --samples <samples.txt> --outdir <dir> [--download] [--jobs N] [--quick-check] [--use-existing-urls] [--ftp]

What it does:
  - For each HGDP sample ID (e.g. HGDP00455), queries ENA for read runs with sample_alias="HGDP00455".
  - Extracts the ENA-provided fastq_ftp field (one or more FASTQ URLs).
  - Writes a deduplicated URL list and an optional download manifest.
  - If --download is provided, downloads all FASTQs with resume support.

Notes:
  - This downloads raw FASTQ and can be very large.
  - Requires: curl, awk, sort, xargs, gzip. (lftp if --ftp is used)
  - --quick-check only verifies the gzip header is present (very fast, least strict).
  - --use-existing-urls skips ENA queries and uses <outdir>/meta/fastq_urls.txt.
  - --ftp uses lftp (pget) instead of curl for downloading.

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
QUICK_CHECK=0
USE_EXISTING_URLS=0
USE_FTP=0

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
    --quick-check)
      QUICK_CHECK=1; shift;;
    --use-existing-urls)
      USE_EXISTING_URLS=1; shift;;
    --ftp)
      USE_FTP=1; shift;;
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

is_gzip_ok() {
  local f="$1"
  if [[ "$QUICK_CHECK" -eq 1 ]]; then
    head -c 2 "$f" | od -An -t x1 | tr -d ' \n' | grep -qi '^1f8b'
  else
    gzip -t "$f" >/dev/null 2>&1
  fi
}

ena_query() {
  local sample="$1"
  # ENA API: TSV with fields (run_accession, fastq_ftp). fastq_ftp is semi-colon separated.
  # limit=0 means no limit.
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=sample_alias%3D%22${sample}%22&fields=run_accession,fastq_ftp&format=tsv&limit=0"
}

if [[ "$USE_EXISTING_URLS" -eq 1 ]]; then
  if [[ ! -f "$URLS_OUT" ]]; then
    echo "URL list not found: $URLS_OUT (disable --use-existing-urls to regenerate)" >&2
    exit 2
  fi
else
  : > "$URLS_TMP"
  : > "$MANIFEST_OUT"
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
fi

echo "Wrote URL list: $URLS_OUT" >&2
echo "Wrote manifest: $MANIFEST_OUT" >&2
echo "FASTQ URLs: $(wc -l < "$URLS_OUT")" >&2

if [[ "$DO_DOWNLOAD" -eq 0 ]]; then
  echo "(dry-run) Not downloading. Re-run with --download to fetch FASTQs." >&2
  exit 0
fi

URLS_DL="$OUTDIR/meta/fastq_urls.to_download.txt"
: > "$URLS_DL"

while IFS= read -r url; do
  [[ -z "$url" ]] && continue
  f="$(basename "$url")"
  dest="$OUTDIR/fastq/$f"
  if [[ -f "$dest" ]] && is_gzip_ok "$dest"; then
    echo "[SKIP] $f (already downloaded)" >&2
    continue
  fi
  echo "$url" >> "$URLS_DL"
done < "$URLS_OUT"

if [[ ! -s "$URLS_DL" ]]; then
  echo "All FASTQ files already present and valid." >&2
  exit 0
fi

# Download all URLs into OUTDIR/fastq (keep basename)
# Use curl resume (-C -) and follow redirects (-L). Parallelize with xargs.
cat "$URLS_DL" \
  | xargs -P "$JOBS" -I {} bash -lc 'set -euo pipefail; u="$1"; outdir="$2"; f="$(basename "$u")"; dest="$outdir/fastq/$f"; quick="$3"; use_ftp="$4"; is_gzip_ok() { if [[ "$quick" -eq 1 ]]; then head -c 2 "$dest" | od -An -t x1 | tr -d " \\n" | grep -qi "^1f8b"; else gzip -t "$dest" >/dev/null 2>&1; fi; }; do_dl_resume() { if [[ "$use_ftp" -eq 1 ]]; then lftp -c "pget -n 4 -c \"$u\" -o \"$dest\""; else curl -fL -C - -o "$dest" "$u"; fi; }; do_dl_fresh() { rm -f "$dest"; if [[ "$use_ftp" -eq 1 ]]; then lftp -c "pget -n 4 -c \"$u\" -o \"$dest\""; else curl -fL -o "$dest" "$u"; fi; }; if [[ -f "$dest" ]]; then echo "[RESUME] $f" >&2; else echo "[DL] $f" >&2; fi; if ! do_dl_resume; then echo "[REDOWNLOAD] $f (resume failed)" >&2; do_dl_fresh; fi; if ! is_gzip_ok; then echo "[REDOWNLOAD] $f (corrupt after download)" >&2; do_dl_fresh; is_gzip_ok; fi' _ {} "$OUTDIR" "$QUICK_CHECK" "$USE_FTP"

echo "Done. Files in: $OUTDIR/fastq" >&2
