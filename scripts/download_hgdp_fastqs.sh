#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  download_hgdp_fastqs.sh --samples <samples.txt> --outdir <dir> [--download] [--jobs N] [--quick-check] [--no-check] [--force] [--use-existing-urls] [--ftp] [--progress] [--verbose]
  download_hgdp_fastqs.sh --fastqs <runs.txt> --outdir <dir> [--download] [--jobs N] [--quick-check] [--no-check] [--force] [--use-existing-urls] [--ftp] [--progress] [--verbose]

What it does:
  - For each HGDP sample ID (e.g. HGDP00455), queries ENA for read runs with sample_alias="HGDP00455".
  - For each run accession (e.g. ERR757835), queries ENA for fastq_ftp URLs.
  - Extracts the ENA-provided fastq_ftp field (one or more FASTQ URLs).
  - Writes a deduplicated URL list and an optional download manifest.
  - If --download is provided, downloads all FASTQs with resume support.

Notes:
  - This downloads raw FASTQ and can be very large.
  - Requires: curl, awk, sort, xargs, gzip. (lftp if --ftp is used)
  - --quick-check only verifies the gzip header is present (very fast, least strict).
  - --no-check disables FASTQ integrity checks (fastest, not recommended).
  - --force always download even if the file already exists.
  - --use-existing-urls skips ENA queries and uses <outdir>/meta/fastq_urls.txt.
  - --ftp uses lftp (pget) instead of curl for downloading.
  - --progress prints a simple per-file progress counter.
  - --verbose prints each URL before download.

Examples:
  # list URLs only
  ./download_hgdp_fastqs.sh --samples inputs/biaka_hgdp.samples.txt --outdir /datos/migccl/ancestry_refs/BiakaHGDP

  # download with 4 parallel transfers
  ./download_hgdp_fastqs.sh --samples inputs/biaka_hgdp.samples.txt --outdir /datos/migccl/ancestry_refs/BiakaHGDP --download --jobs 4
EOF
}

SAMPLES=""
FASTQS=""
OUTDIR=""
DO_DOWNLOAD=0
JOBS=4
QUICK_CHECK=0
NO_CHECK=0
FORCE=0
USE_EXISTING_URLS=0
USE_FTP=0
PROGRESS=0
VERBOSE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples)
      SAMPLES="$2"; shift 2;;
    --fastqs)
      FASTQS="$2"; shift 2;;
    --outdir)
      OUTDIR="$2"; shift 2;;
    --download)
      DO_DOWNLOAD=1; shift;;
    --jobs)
      JOBS="$2"; shift 2;;
    --quick-check)
      QUICK_CHECK=1; shift;;
    --no-check)
      NO_CHECK=1; shift;;
    --force)
      FORCE=1; shift;;
    --use-existing-urls)
      USE_EXISTING_URLS=1; shift;;
    --ftp)
      USE_FTP=1; shift;;
    --progress)
      PROGRESS=1; shift;;
    --verbose)
      VERBOSE=1; shift;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown arg: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -n "$SAMPLES" && -n "$FASTQS" ]]; then
  echo "Provide only one of --samples or --fastqs" >&2
  usage
  exit 2
fi

if [[ -z "$SAMPLES" && -z "$FASTQS" ]]; then
  usage
  exit 2
fi

if [[ -n "$SAMPLES" ]]; then
  if [[ ! -f "$SAMPLES" ]]; then
    echo "Samples file not found: $SAMPLES" >&2
    exit 2
  fi
fi

if [[ -n "$FASTQS" ]]; then
  if [[ ! -f "$FASTQS" ]]; then
    echo "FASTQ run list file not found: $FASTQS" >&2
    exit 2
  fi
fi

mkdir -p "$OUTDIR" "$OUTDIR/fastq" "$OUTDIR/meta"

URLS_TMP="$OUTDIR/meta/fastq_urls.tmp.txt"
URLS_OUT="$OUTDIR/meta/fastq_urls.txt"
MANIFEST_OUT="$OUTDIR/meta/fastq_manifest.tsv"
URLS_SRC="$URLS_OUT"
MANIFEST_SRC="$MANIFEST_OUT"

is_gzip_ok() {
  local f="$1"
  if [[ "$NO_CHECK" -eq 1 ]]; then
    return 0
  fi
  if [[ "$QUICK_CHECK" -eq 1 ]]; then
    head -c 2 "$f" | od -An -t x1 | tr -d ' \n' | grep -qi '^1f8b'
  else
    gzip -t "$f" >/dev/null 2>&1
  fi
}

ena_query_sample() {
  local sample="$1"
  # ENA API: TSV with fields (run_accession, sample_alias, fastq_ftp). fastq_ftp is semi-colon separated.
  # limit=0 means no limit.
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=sample_alias%3D%22${sample}%22&fields=run_accession,sample_alias,fastq_ftp&format=tsv&limit=0"
}

ena_query_run() {
  local run="$1"
  # ENA API: TSV with fields (run_accession, sample_alias, fastq_ftp). fastq_ftp is semi-colon separated.
  # limit=0 means no limit.
  curl -fsSL "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22${run}%22&fields=run_accession,sample_alias,fastq_ftp&format=tsv&limit=0"
}

filter_urls_by_runs() {
  local urls_in="$1"
  local manifest_in="$2"
  local runs_file="$3"
  local urls_out="$4"
  local manifest_out="$5"

  declare -A RUNS=()
  while IFS= read -r run; do
    run="${run//$'\r'/}"
    [[ -z "$run" ]] && continue
    RUNS["$run"]=1
  done < "$runs_file"

  : > "$urls_out"
  : > "$manifest_out"

  if [[ -s "$manifest_in" ]]; then
    while IFS=$'\t' read -r sample run url; do
      [[ -z "$run" || -z "$url" ]] && continue
      if [[ -n "${RUNS[$run]:-}" ]]; then
        echo "$url" >> "$urls_out"
        printf '%s\t%s\t%s\n' "$sample" "$run" "$url" >> "$manifest_out"
      fi
    done < "$manifest_in"
  else
    while IFS= read -r url; do
      [[ -z "$url" ]] && continue
      local base run
      base="$(basename "$url")"
      if [[ "$base" =~ ^([A-Za-z]+[0-9]+) ]]; then
        run="${BASH_REMATCH[1]}"
      else
        run=""
      fi
      if [[ -n "$run" && -n "${RUNS[$run]:-}" ]]; then
        echo "$url" >> "$urls_out"
        printf 'UNKNOWN\t%s\t%s\n' "$run" "$url" >> "$manifest_out"
      fi
    done < "$urls_in"
  fi

  if [[ -s "$urls_out" ]]; then
    sort -u -o "$urls_out" "$urls_out"
  fi
}

if [[ "$USE_EXISTING_URLS" -eq 1 ]]; then
  if [[ ! -f "$URLS_OUT" ]]; then
    echo "URL list not found: $URLS_OUT (disable --use-existing-urls to regenerate)" >&2
    exit 2
  fi
  if [[ -n "$FASTQS" ]]; then
    URLS_SRC="$OUTDIR/meta/fastq_urls.filtered.txt"
    MANIFEST_SRC="$OUTDIR/meta/fastq_manifest.filtered.tsv"
    filter_urls_by_runs "$URLS_OUT" "$MANIFEST_OUT" "$FASTQS" "$URLS_SRC" "$MANIFEST_SRC"
  fi
else
  : > "$URLS_TMP"
  : > "$MANIFEST_OUT"

  if [[ -n "$SAMPLES" ]]; then
    while IFS= read -r sample; do
      sample="${sample//$'\r'/}"
      [[ -z "$sample" ]] && continue

      echo "[ENA] sample ${sample}" >&2

      # Output TSV has header. Extract run_accession + sample_alias + fastq_ftp.
      # Some samples can have multiple runs/lanes.
      ena_query_sample "$sample" \
        | awk -F'\t' 'NR==1{next} $1!=""{print $1"\t"$2"\t"$3}' \
        | while IFS=$'\t' read -r run sample_alias fastq_ftp; do
            [[ -z "$run" ]] && continue
            [[ -z "$fastq_ftp" ]] && continue

            # fastq_ftp: host/path1.gz;host/path2.gz (no scheme). Convert to ftp://...
            IFS=';' read -ra parts <<< "$fastq_ftp"
            for p in "${parts[@]}"; do
              [[ -z "$p" ]] && continue
              url="ftp://${p}"
              echo "$url" >> "$URLS_TMP"
              printf '%s\t%s\t%s\n' "${sample_alias:-$sample}" "$run" "$url" >> "$MANIFEST_OUT"
            done
          done

    done < "$SAMPLES"
  else
    while IFS= read -r run; do
      run="${run//$'\r'/}"
      [[ -z "$run" ]] && continue

      echo "[ENA] run ${run}" >&2

      ena_query_run "$run" \
        | awk -F'\t' 'NR==1{next} $1!=""{print $1"\t"$2"\t"$3}' \
        | while IFS=$'\t' read -r run_accession sample_alias fastq_ftp; do
            [[ -z "$run_accession" ]] && continue
            [[ -z "$fastq_ftp" ]] && continue

            IFS=';' read -ra parts <<< "$fastq_ftp"
            for p in "${parts[@]}"; do
              [[ -z "$p" ]] && continue
              url="ftp://${p}"
              echo "$url" >> "$URLS_TMP"
              printf '%s\t%s\t%s\n' "${sample_alias:-UNKNOWN}" "$run_accession" "$url" >> "$MANIFEST_OUT"
            done
          done
    done < "$FASTQS"
  fi

# Deduplicate URLs
sort -u "$URLS_TMP" > "$URLS_OUT"
rm -f "$URLS_TMP"
fi

if [[ "$USE_EXISTING_URLS" -eq 1 ]]; then
  echo "Using URL list: $URLS_OUT" >&2
  echo "Using manifest: $MANIFEST_OUT" >&2
  if [[ "$URLS_SRC" != "$URLS_OUT" ]]; then
    echo "Filtered URL list: $URLS_SRC" >&2
    echo "Filtered manifest: $MANIFEST_SRC" >&2
  fi
else
  echo "Wrote URL list: $URLS_OUT" >&2
  echo "Wrote manifest: $MANIFEST_OUT" >&2
fi
echo "FASTQ URLs: $(wc -l < "$URLS_SRC" | tr -d ' ')" >&2

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
  if [[ "$FORCE" -eq 0 ]]; then
    if [[ -f "$dest" ]] && is_gzip_ok "$dest"; then
      echo "[SKIP] $f (already downloaded)" >&2
      continue
    fi
  fi
  echo "$url" >> "$URLS_DL"
done < "$URLS_SRC"

if [[ ! -s "$URLS_DL" ]]; then
  echo "All FASTQ files already present and valid." >&2
  exit 0
fi

# Download all URLs into OUTDIR/fastq (keep basename)
# Use curl resume (-C -) and follow redirects (-L). Parallelize with xargs.
TOTAL_DL="$(wc -l < "$URLS_DL" | tr -d ' ')"
nl -ba "$URLS_DL" \
  | xargs -P "$JOBS" -I {} bash -lc 'set -euo pipefail; line="$1"; outdir="$2"; quick="$3"; use_ftp="$4"; progress="$5"; verbose="$6"; total="$7"; force="$8"; IFS=$'"'"'\t'"'"' read -r idx u <<< "$line"; f="$(basename "$u")"; dest="$outdir/fastq/$f"; is_gzip_ok() { if [[ "$quick" -eq 1 ]]; then head -c 2 "$dest" | od -An -t x1 | tr -d " \\n" | grep -qi "^1f8b"; else gzip -t "$dest" >/dev/null 2>&1; fi; }; do_dl_resume() { if [[ "$use_ftp" -eq 1 ]]; then lftp -c "pget -n 4 -c \"$u\" -o \"$dest\""; else curl -fL -C - -o "$dest" "$u"; fi; }; do_dl_fresh() { rm -f "$dest"; if [[ "$use_ftp" -eq 1 ]]; then lftp -c "pget -n 4 -c \"$u\" -o \"$dest\""; else curl -fL -o "$dest" "$u"; fi; }; if [[ "$progress" -eq 1 ]]; then if [[ -f "$dest" ]]; then echo "[$idx/$total] RESUME $f" >&2; else echo "[$idx/$total] DL $f" >&2; fi; else if [[ -f "$dest" ]]; then echo "[RESUME] $f" >&2; else echo "[DL] $f" >&2; fi; fi; if [[ "$verbose" -eq 1 ]]; then echo "URL: $u" >&2; fi; if [[ "$force" -eq 1 ]]; then echo "[FORCE] $f" >&2; do_dl_fresh; else if ! do_dl_resume; then echo "[REDOWNLOAD] $f (resume failed)" >&2; do_dl_fresh; fi; if ! is_gzip_ok; then echo "[REDOWNLOAD] $f (corrupt after download)" >&2; do_dl_fresh; is_gzip_ok; fi; fi' _ {} "$OUTDIR" "$QUICK_CHECK" "$USE_FTP" "$PROGRESS" "$VERBOSE" "$TOTAL_DL" "$FORCE"

echo "Done. Files in: $OUTDIR/fastq" >&2
