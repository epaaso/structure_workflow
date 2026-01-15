# Ancestry Inference Workflow

This Nextflow pipeline performs ancestry inference for VCF samples using 1000 Genomes trio references (hg38). It automates the process of downloading reference data, merging with sample data, performing quality control, and running ancestry inference tools.

## Features

- **Reference Data**: Automatically downloads 1000 Genomes Project reference data (hg38).
- **Methods**: Supports both **ADMIXTURE** and **STRUCTURE** for ancestry inference.
- **QC & Pruning**: Includes steps for filtering (MAF, missingness) and Linkage Disequilibrium (LD) pruning.
- **Parallelization**: efficient processing of chromosomes and samples using Nextflow.

## Requirements

- [Nextflow](https://www.nextflow.io/)
- [Conda](https://docs.conda.io/en/latest/) (recommended for dependency management)
- `bwa`, `samtools`, and `bcftools` (only required if `--build_ref_vcfs` is enabled)

Dependencies are defined in `environment.yml`.

## Usage

Run the pipeline with Nextflow:

```bash
nextflow run main.nf --input "path/to/your/vcfs/*.vcf.gz" --outdir results
```

### Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Glob pattern for input VCF files (must be bgzipped and indexed) | `vcfs/*.vcf.gz` |
| `--outdir` | Directory for output results | `results` |
| `--method` | Inference method to use: `admixture` or `structure` | `admixture` |
| `--k` | Number of ancestral populations (K) | `3` |
| `--chroms` | Chromosomes to analyze (e.g., "1-22", "chr1,chr2") | `1-22` |
| `--exome_bed` | BED file of target regions used to subset the **reference** 1KG VCFs before merging (chr prefix is auto-normalized) | `exome_hg38.bed` |
| `--maf` | Minor Allele Frequency filter | `0.05` |
| `--geno` | Genotype missingness filter | `0.05` |
| `--mind` | Sample missingness filter (passed to PLINK `--mind`) | `0.8` |
| `--ld_window` | LD pruning window (passed to `--indep-pairwise`) | `50` |
| `--ld_step` | LD pruning step size (passed to `--indep-pairwise`) | `5` |
| `--ld_r2` | LD pruning r² threshold (passed to `--indep-pairwise`) | `0.2` |
| `--r2_sweep` | Enable LD r² sweep (PCA + ADMIXTURE checks) | `false` |
| `--r2_values` | Comma-separated r² values for sweep | `0.01,0.02,0.05,0.1,0.2,0.3` |
| `--r2_ld_window` | LD window used during r² sweep (passed to `--indep-pairwise`) | `500` |
| `--r2_ld_step` | LD step used during r² sweep (passed to `--indep-pairwise`) | `10` |
| `--r2_pca_components` | Number of PCA components computed during r² sweep | `10` |
| `--seed` | Random seed for ADMIXTURE runs | `42` |
| `--ref_limit` | Max samples per reference group from the 1KG panel (`0` or `all` = no limit) | `40` |
| `--ref_mxl` | Custom MXL reference list (one IID per line); requires `--ref_eur` and `--ref_afr` | unset |
| `--ref_eur` | Custom EUR reference list (one IID per line); requires `--ref_mxl` and `--ref_afr` | unset |
| `--ref_afr` | Custom AFR reference list (one IID per line); requires `--ref_mxl` and `--ref_eur` | unset |
| `--ref_mxl_pop` | 1KG population code used to select MXL references from the panel | `MXL` |
| `--ref_eur_super` | 1KG super-population used to select EUR references from the panel | `EUR` |
| `--ref_afr_super` | 1KG super-population used to select AFR references from the panel (unless `--ref_afr_pop` is set) | `AFR` |
| `--ref_afr_pop` | Specific AFR population to use (e.g., `YRI`); overrides `--ref_afr_super` | unset |
| `--build_ref_vcfs` | Build per-sample VCFs from paired-end FASTQs (runs only the FASTQ→VCF routine) | `false` |
| `--fastq_dir` | Base directory containing per-cohort `fastq/` folders | `/datos/migccl/ancestry_refs` |
| `--fastq_pattern` | FASTQ glob pattern (relative to each `fastq/` folder) | `*_{1,2}.fastq.gz` |
| `--fastq_manifest` | Optional manifest file (`sample_id<TAB>read1<TAB>read2`) | unset |
| `--ref_fasta` | Path to hg38 FASTA (used for FASTQ conversion); if unset, download from `--ref_fasta_url` | unset |
| `--ref_fasta_url` | hg38 FASTA URL (used when `--ref_fasta` is unset) | `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz` |
| `--bwa` | Path to `bwa` | `bwa` |
| `--samtools` | Path to `samtools` | `samtools` |
| `--bcftools` | Path to `bcftools` | `bcftools` |
| `--plink2` | Path to `plink2` binary (used for VCF→PLINK conversion + QC + default pruning) | `plink2` |
| `--plink` | Path to PLINK 1.9 binary (used for PCA, STRUCTURE conversion, and r² sweep pruning) | `/home/epaaso/bin/plink` |
| `--admixture` | Path to `admixture` binary | `admixture` |
| `--structure` | Path to `structure` binary | `bin/structure` |
| `--max_cpus` | Per-process CPU limit (also used by `nextflow.config` to set default `process.cpus`) | `18` |
| `--max_memory` | Per-process memory limit (also used by `nextflow.config` to set default `process.memory`) | `120 GB` |
| `--burnin` | STRUCTURE `BURNIN` iterations (only used when `--method structure`) | `1000` |
| `--numreps` | STRUCTURE `NUMREPS` iterations (only used when `--method structure`) | `1000` |
| `--panel_url` | URL to the 1KG panel file used to create reference sample lists | `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel` |
| `--vcf_base_url` | Base URL used to download 1KG per-chromosome hg38 VCFs | `https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes` |

## Output

The pipeline produces results in the specified output directory (default: `results/`):

- **ancestry/**: Contains the output from ADMIXTURE or STRUCTURE (Q files, P files).
- **plots/**: Visualization of ancestry proportions (if plotting is enabled/implemented).
- **merged_all/**: Merged VCF/PLINK datasets of samples and reference.
- **pca/**: PCA results (eigenvectors, eigenvalues).
- **r2_sweep/**: LD r² sweep outputs (pruned sets, PCA, ADMIXTURE CV logs).
- **reference/fastq_vcfs/**: Per-sample VCFs generated from FASTQs (when `--build_ref_vcfs` is enabled).

## Directory Structure

- `main.nf`: The main Nextflow pipeline script.
- `nextflow.config`: Configuration file for the pipeline.
- `environment.yml`: Conda environment definition.
- `bin/`: Contains executable scripts or binaries (e.g., `structure`).

## Utilities

### Optional: build reference VCFs from paired-end FASTQs

If you have reference FASTQs (paired-end, `_1` / `_2` suffixes), you can generate per-sample VCFs. This mode only runs the FASTQ→VCF routine:

```bash
nextflow run main.nf \
  --method admixture \
  --build_ref_vcfs true \
  --fastq_dir /datos/migccl/ancestry_refs \
  --ref_fasta_url https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz \
  --outdir results_fastq_refs
```

The FASTQ-derived VCFs land under `reference/fastq_vcfs/` inside your `--outdir`. Use those VCFs as input for the main pipeline in a separate run.
If your files end in `.fq.gz`, set `--fastq_pattern "*_{1,2}.fq.gz"`.
If `reference/genome/ref.fa` already exists in your `--outdir`, the workflow will reuse it and skip downloading/indexing.

You can also provide a manifest (`sample_id<TAB>read1<TAB>read2`) to control sample IDs:

```bash
nextflow run main.nf \
  --build_ref_vcfs true \
  --fastq_manifest /path/to/fastq_manifest.tsv \
  --ref_fasta /path/to/hg38.fa \
  --outdir results_fastq_refs
```

Then run the main workflow using the generated VCFs:

```bash
nextflow run main.nf \
  --input "results_fastq_refs/reference/fastq_vcfs/*.vcf.gz" \
  --outdir results_from_fastqs
```

### Check which AIM alleles are present in your samples

If you have a CSV of AIMs with at least these columns: `chr`, `position`, `A1`, `A2`, and `SNP rsID` (or `rsid`), you can scan your bgzipped VCFs and output a TSV summary:

```bash
python scripts/aim_alleles_from_vcfs.py \
	--aim path/to/AIM_latin.csv \
	--vcf "/datos/migccl/vcfs/*.vcf.gz" \
	--out results/aim_alleles_from_vcfs.tsv
```

The output includes, for each AIM SNP, whether `A1` and/or `A2` is observed in any sample genotype, plus basic counters and notes when the AIM alleles do not match the VCF `REF/ALT`.

### Select “pure” reference samples (from stage1 ancestry output)

After you run the pipeline once (stage1) and produce an `all_ancestry.tsv` (from either ADMIXTURE or STRUCTURE), you can select a stricter set of “pure” reference samples by thresholding ancestry proportions and writing new reference lists (one IID per line). These lists can then be fed back into the pipeline via `--ref_mxl/--ref_eur/--ref_afr` for a second run (stage2).

Script: `scripts/select_pure_refs.py`

Inputs:

- `--panel`: 1KG panel file (e.g. `integrated_call_samples_v3.20130502.ALL.panel`)
- `--ancestry`: pipeline `all_ancestry.tsv` (must have an IID column and component columns matching `--*-comp`)

Outputs:

- `--out-mxl`, `--out-eur`, `--out-afr`: text files with one IID per line

Parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--panel` | 1KG panel file used to map IID → population / super-population | required |
| `--ancestry` | `all_ancestry.tsv` from the pipeline | required |
| `--mxl-pop` | Panel population code used to define the MXL candidate set | `MXL` |
| `--eur-super` | Panel super-population code used to define the EUR candidate set | `EUR` |
| `--afr-pop` | Panel population code used to define the AFR candidate set (if non-empty) | `YRI` |
| `--afr-super` | Panel super-population code for AFR (used only when `--afr-pop` is empty) | `AFR` |
| `--mxl-comp` | Column name in `all_ancestry.tsv` for the MXL ancestry component | `MXL` |
| `--eur-comp` | Column name in `all_ancestry.tsv` for the EUR ancestry component | `EUR` |
| `--afr-comp` | Column name in `all_ancestry.tsv` for the AFR ancestry component | `AFR` |
| `--min-mxl` | Minimum MXL proportion to keep a sample | `0.85` |
| `--min-eur` | Minimum EUR proportion to keep a sample | `0.90` |
| `--min-afr` | Minimum AFR proportion to keep a sample | `0.90` |
| `--max-per-group` | Optional cap per group after sorting by ancestry proportion (`0` = no cap) | `0` |
| `--out-mxl` | Output IID list for the selected MXL refs | required |
| `--out-eur` | Output IID list for the selected EUR refs | required |
| `--out-afr` | Output IID list for the selected AFR/YRI refs | required |

Example (stage1 → pick pure refs → stage2):

```bash
python3 scripts/select_pure_refs.py \
	--panel results_pure_stage1/reference/meta/integrated_call_samples_v3.20130502.ALL.panel \
	--ancestry results_pure_stage1/plots/admixture/all_ancestry.tsv \
	--mxl-pop MXL --eur-super EUR --afr-pop YRI \
	--min-mxl 0.5 --min-eur 0.95 --min-afr 0.95 --max-per-group 40 \
	--out-mxl results_pure_stage1/reference/meta/mxl_pure.samples \
	--out-eur results_pure_stage1/reference/meta/eur_pure.samples \
	--out-afr results_pure_stage1/reference/meta/yri_pure.samples

nohup nextflow run main.nf -resume --method admixture \
	--input "results/samples/indexed/*.vcf.gz" \
	--outdir results_pure_stage2 \
	--ref_mxl results_pure_stage1/reference/meta/mxl_pure.samples \
	--ref_eur results_pure_stage1/reference/meta/eur_pure.samples \
	--ref_afr results_pure_stage1/reference/meta/yri_pure.samples \
	> results_pure
```

### Download HGDP FASTQs (ENA/FTP)

To download raw HGDP FASTQ files for a list of HGDP sample IDs (e.g., `HGDP00455`) using ENA metadata (the `fastq_ftp` field), use `scripts/download_hgdp_fastqs.sh`.

Inputs:

- A sample list file with one HGDP ID per line (example: `inputs/biaka_hgdp.samples.txt`)

This writes a deduplicated URL list and manifest under your output directory:

- `<outdir>/meta/fastq_urls.txt`
- `<outdir>/meta/fastq_manifest.tsv`

List URLs only (recommended first):

```bash
./scripts/download_hgdp_fastqs.sh \
	--samples inputs/biaka_hgdp.samples.txt \
	--outdir /datos/migccl/ancestry_refs/BiakaHGDP
```

Download (large; uses resume, parallelizable):

```bash
./scripts/download_hgdp_fastqs.sh \
	--samples inputs/biaka_hgdp.samples.txt \
	--outdir /datos/migccl/ancestry_refs/BiakaHGDP \
	--download --jobs 4
```
