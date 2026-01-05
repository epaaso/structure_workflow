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
| `--maf` | Minor Allele Frequency filter | `0.05` |
| `--geno` | Genotype missingness filter | `0.05` |
| `--r2_sweep` | Enable LD r² sweep (PCA + ADMIXTURE checks) | `false` |
| `--r2_values` | Comma-separated r² values for sweep | `0.01,0.02,0.05,0.1,0.2,0.3` |

## Output

The pipeline produces results in the specified output directory (default: `results/`):

- **ancestry/**: Contains the output from ADMIXTURE or STRUCTURE (Q files, P files).
- **plots/**: Visualization of ancestry proportions (if plotting is enabled/implemented).
- **merged_all/**: Merged VCF/PLINK datasets of samples and reference.
- **pca/**: PCA results (eigenvectors, eigenvalues).
- **r2_sweep/**: LD r² sweep outputs (pruned sets, PCA, ADMIXTURE CV logs).

## Directory Structure

- `main.nf`: The main Nextflow pipeline script.
- `nextflow.config`: Configuration file for the pipeline.
- `environment.yml`: Conda environment definition.
- `bin/`: Contains executable scripts or binaries (e.g., `structure`).

## Utilities

### Check which AIM alleles are present in your samples

If you have a CSV of AIMs with at least these columns: `chr`, `position`, `A1`, `A2`, and `SNP rsID` (or `rsid`), you can scan your bgzipped VCFs and output a TSV summary:

```bash
python scripts/aim_alleles_from_vcfs.py \
	--aim path/to/AIM_latin.csv \
	--vcf "/datos/migccl/vcfs/*.vcf.gz" \
	--out results/aim_alleles_from_vcfs.tsv
```

The output includes, for each AIM SNP, whether `A1` and/or `A2` is observed in any sample genotype, plus basic counters and notes when the AIM alleles do not match the VCF `REF/ALT`.
