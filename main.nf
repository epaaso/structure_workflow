nextflow.enable.dsl=2

/*
 * Pipeline: ancestry inference for VCF samples using 1000 Genomes trio references (hg38)
 * Inputs:  params.input -> glob for sample VCFs (bgzipped + indexed) e.g. 'vcfs/*.vcf.gz'
 * Outputs: per-sample admixture proportions and intermediate PLINK/VCF artifacts under params.outdir
 * Tools required on PATH: bcftools, admixture, plink2 (params.plink2 can override path)
 */

// Default parameters
params.input   = "vcfs/*.vcf.gz"
params.outdir  = "results"
params.chroms  = "1-22"      // comma list or ranges like 1-22 or chr1,chr2
params.k       = 3           // number of ancestries
params.maf     = 0.05
params.geno    = 0.05
params.mind    = 0.8
params.ld_window = 50
params.ld_step   = 5
params.ld_r2     = 0.2
params.seed      = 42
params.panel_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
params.vcf_base_url = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes"
params.plink2    = "plink2"
params.plink     = "/home/epaaso/bin/plink"
params.admixture = "admixture"
params.structure = "${workflow.projectDir}/bin/structure"
params.method    = "admixture" // "admixture" or "structure"
params.max_cpus  = 16          // per-process max; overridden by --max_cpus
params.max_memory = '120 GB'   // per-process max; overridden by --max_memory
params.burnin    = 1000
params.numreps   = 1000
params.exome_bed = "${workflow.projectDir}/exome_hg38.bed"

def parseChroms(value) {
    def tokens = value.toString().split(/[, ]+/).findAll { it }
    if (!tokens) {
        return (1..22).collect { "chr${it}" }
    }
    def expanded = tokens.collectMany { tok ->
        if (tok.contains("-")) {
            def (start, end) = tok.split("-")*.replaceFirst(/^chr/, "")*.toInteger()
            (start..end).collect { it }
        } else {
            [tok.replaceFirst(/^chr/, "")]
        }
    }
    return expanded.collect { "chr${it}" }
}

workflow {
    def chroms = parseChroms(params.chroms)
    def plink_ready
    def plink_pca

    // 1. Download Panel and Create Sample Lists
    DOWNLOAD_PANEL()
    CREATE_REF_LISTS(DOWNLOAD_PANEL.out)
    
    // Split outputs
    CREATE_REF_LISTS.out
        .map { mxl, eur, afr, all -> tuple(mxl, eur, afr) }
        .set { ref_lists_separate }
        
    CREATE_REF_LISTS.out
        .map { mxl, eur, afr, all -> all }
        .set { ref_list_all }

    // 2. Download 1KG VCFs per chromosome
    Channel.from(chroms).set { chrom_ch }
    DOWNLOAD_1KG_VCF(chrom_ch)

    // 3. Subset VCFs for all reference samples
    DOWNLOAD_1KG_VCF.out
        .combine(ref_list_all)
        .set { subset_requests }

    SUBSET_REF_VCF(subset_requests)
    // Output: tuple(chrom, vcf, tbi)

    FILTER_REF_CHROM(SUBSET_REF_VCF.out, file(params.exome_bed))

    // Sample VCFs
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { vcf ->
            def id = vcf.name.replaceAll(/\.vcf(\.gz|\.bgz)?$/, "")
            tuple(id, vcf)
        }
        .set { samples }

    INDEX_SAMPLE(samples)
    
    // Collect sample artifacts for joint reference merge
    INDEX_SAMPLE.out
        .map { id, vcf, tbi -> [ vcf, tbi ] }
        .collect()
        .map { items ->
            def files = (items && items[0] instanceof List) ? items.collectMany { it } : items
            [ files as List ]
        }
        .set { all_sample_files }

    // Merge all samples with ref per chromosome
    FILTER_REF_CHROM.out
        .combine(all_sample_files)
        .set { merge_all_requests }

    MERGE_ALL_SAMPLES_CHROM(merge_all_requests)
    
    MERGE_ALL_SAMPLES_CHROM.out
        .map { vcf, tbi -> vcf }
        .collect()
        .set { all_merged_chrom_vcfs }

    CONCAT_ALL_SAMPLES(all_merged_chrom_vcfs)
    .set { merged_all }
    
    def plink_all = PREPARE_PLINK(merged_all)
    plink_ready = plink_all
    plink_pca = plink_all

    if (params.method == "structure") {
        CONVERT_TO_STRUCTURE_ALL(plink_pca, ref_lists_separate)
        RUN_STRUCTURE_ALL(CONVERT_TO_STRUCTURE_ALL.out)
        SUMMARIZE_STRUCTURE_ALL(RUN_STRUCTURE_ALL.out, ref_lists_separate)
        SUMMARIZE_STRUCTURE_ALL.out
            .map { id, tsv -> tsv }
            .collect()
            .set { structure_ancestry_files }

        PLOT_ANCESTRY(structure_ancestry_files, ref_lists_separate)
    } else {
        RUN_ADMIXTURE(plink_ready)
        .set { admixture_outputs }

        SUMMARIZE_Q(admixture_outputs, ref_lists_separate)
        .map { id, tsv -> tsv }
        .collect()
        .set { all_tsvs }
        
        PLOT_ANCESTRY(all_tsvs, ref_lists_separate)
    }
    
    RUN_PCA_ALL(plink_pca)
    .set { pca_results }
    
    PLOT_PCA_ALL(pca_results, ref_lists_separate)
}

process MERGE_ALL_SAMPLES_CHROM {
    tag { chrom }
    
    input:
    tuple val(chrom), path(ref_vcf), path(ref_tbi), path(sample_files)

    output:
    tuple path("all_samples_ref_${chrom}.vcf.gz"), path("all_samples_ref_${chrom}.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    
    REF_IN="${ref_vcf}"
    if [ ! -f "\$REF_IN" ]; then
        echo "ERROR: Expected reference VCF '\$REF_IN' not found in task work dir" >&2
        ls -la >&2
        exit 1
    fi

    # Ensure Ref has 'chr' prefix to match samples (assuming samples have chr)
    if ! bcftools view -h "\$REF_IN" | grep -q "##contig=<ID=chr"; then
        for i in {1..22} X Y MT; do echo "\$i chr\$i" >> add_chr.txt; done
        bcftools annotate --threads ${task.cpus} --rename-chrs add_chr.txt "\$REF_IN" -Oz -o ref_renamed.vcf.gz
        bcftools index -t ref_renamed.vcf.gz
        REF_VCF=ref_renamed.vcf.gz
    else
        REF_VCF="\$REF_IN"
    fi
    
    # Determine which sites exist in *any* sample for this chromosome.
    # We intentionally drop variants present only in the reference panel to avoid batch effects.
    # The sample VCFs and their .tbi indices are staged via `sample_files`.
    ls *.vcf.gz \\
      | grep -v '^ref_' \\
      | sort -V > sample_vcfs.list
    
    if [ ! -s sample_vcfs.list ]; then
        echo "ERROR: No sample VCFs were staged for merge on ${chrom}" >&2
        ls -la >&2
        exit 1
    fi

    # bcftools query cannot query multiple VCFs with different sample sets at once,
    # so we query each sample VCF and take the union of positions.
    : > sample_sites.raw.tsv
    while read -r vcf; do
        bcftools query -f '%CHROM\\t%POS\\n' -r ${chrom} "\$vcf"
    done < sample_vcfs.list >> sample_sites.raw.tsv

    sort -u -k1,1 -k2,2n sample_sites.raw.tsv > sample_sites.tsv
    rm -f sample_sites.raw.tsv

    # Merge
    # We do NOT use --missing-to-ref to avoid batch effects.
    # Missing data will be handled by downstream tools (Plink/Structure).
    bcftools merge --threads ${task.cpus} \$REF_VCF -l sample_vcfs.list -r ${chrom} -Oz -o merged_all_${chrom}.vcf.gz
    bcftools index -t -f merged_all_${chrom}.vcf.gz

    # Drop reference-only sites (keep only sites observed in samples)
    bcftools view --threads ${task.cpus} -T sample_sites.tsv merged_all_${chrom}.vcf.gz -Oz -o all_samples_ref_${chrom}.vcf.gz
    bcftools index -t -f all_samples_ref_${chrom}.vcf.gz

    rm -f merged_all_${chrom}.vcf.gz merged_all_${chrom}.vcf.gz.csi merged_all_${chrom}.vcf.gz.tbi
    """
}

process CONCAT_ALL_SAMPLES {
    publishDir "${params.outdir}/merged_all", mode: 'copy'
    
    input:
    path chrom_vcfs

    output:
    tuple val("all_samples"), path("all_samples_ref.vcf.gz"), path("all_samples_ref.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    ls ${chrom_vcfs} | sort -V > file_list.txt
    
    bcftools concat --threads ${task.cpus} -f file_list.txt -Oz -o all_samples_ref.vcf.gz
    bcftools index -t -f all_samples_ref.vcf.gz
    """
}

process PREPARE_PLINK {
    tag { sample_id }
    publishDir "${params.outdir}/plink_all", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_vcf), path(merged_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}_pruned.bed"),
          path("${sample_id}_pruned.bim"),
          path("${sample_id}_pruned.fam")

    script:
    """
    set -euo pipefail
    plink=${params.plink2}
    mem_mb=${task.memory.toMega()}

    # Sort variants first as requested by plink2 if chromosomes are split
    \$plink --threads ${task.cpus} --memory \$mem_mb --vcf ${merged_vcf} --chr ${params.chroms} --max-alleles 2 --make-pgen --sort-vars --out ${sample_id}_sorted
    
    # Convert to BED and continue
    \$plink --threads ${task.cpus} --memory \$mem_mb --pfile ${sample_id}_sorted --make-bed --out ${sample_id}_base
    \$plink --threads ${task.cpus} --memory \$mem_mb --bfile ${sample_id}_base --snps-only just-acgt --max-alleles 2 --make-bed --out ${sample_id}_snps
    \$plink --threads ${task.cpus} --memory \$mem_mb --bfile ${sample_id}_snps --maf ${params.maf} --geno ${params.geno} --mind ${params.mind} --make-bed --out ${sample_id}_clean
    \$plink --threads ${task.cpus} --memory \$mem_mb --bfile ${sample_id}_clean --set-all-var-ids @:# --rm-dup exclude-all --make-bed --out ${sample_id}_dedup
    \$plink --threads ${task.cpus} --memory \$mem_mb --bfile ${sample_id}_dedup --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} --bad-ld --out ${sample_id}_prune
    \$plink --threads ${task.cpus} --memory \$mem_mb --bfile ${sample_id}_dedup --extract ${sample_id}_prune.prune.in --make-bed --out ${sample_id}_pruned
    """
}


process RUN_PCA_ALL {
    publishDir "${params.outdir}/pca_all", mode: 'copy'
    tag { sample_id }
    // No conda here to avoid path issues

    input:
    tuple val(sample_id), path(bed), path(bim), path(fam)

    output:
    tuple val(sample_id), path("${sample_id}.eigenvec"), path("${sample_id}.eigenval")

    script:
    """
    set -euo pipefail
    mem_mb=${task.memory.toMega()}
    # Use plink 1.9 for PCA
    ${params.plink} --threads ${task.cpus} --memory \$mem_mb --bfile ${bed.baseName} --pca header --out ${sample_id}
    """
}

process PLOT_PCA_ALL {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    tuple val(sample_id), path(eigenvec), path(eigenval)
    tuple path(mxl_list), path(eur_list), path(afr_list)

    output:
    path "all_samples_pca.png"

    script:
    """
    python3 - <<'PY'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

eigenvec = "${eigenvec}"
mxl_file = "${mxl_list}"
eur_file = "${eur_list}"
afr_file = "${afr_list}"

# Load reference samples
ref_groups = {}
for pop, fpath in [("MXL", mxl_file), ("EUR", eur_file), ("AFR", afr_file)]:
    with open(fpath, 'r') as f:
        ref_groups[pop] = set(line.strip() for line in f if line.strip())

# Read eigenvec
try:
    df = pd.read_csv(eigenvec, delim_whitespace=True)
except:
    df = pd.read_csv(eigenvec, sep='\\t')

# Rename columns if needed
col_map = {c: c.replace('#', '') for c in df.columns}
df.rename(columns=col_map, inplace=True)

# Ensure we have IID
if 'IID' not in df.columns:
    if len(df.columns) >= 3:
        df.rename(columns={df.columns[0]: 'IID'}, inplace=True)

def get_pop(iid):
    for pop, samples in ref_groups.items():
        if iid in samples:
            return pop
    # If not reference, it's a sample
    return "Sample"

df['Population'] = df['IID'].apply(get_pop)

# Plot
plt.figure(figsize=(12, 10))
sns.scatterplot(data=df, x='PC1', y='PC2', hue='Population', style='Population', s=100)

# Add labels for samples
for i, row in df.iterrows():
    if row['Population'] == "Sample":
        plt.text(row['PC1']+0.002, row['PC2']+0.002, row['IID'], fontsize=9)

plt.title('PCA: All Samples vs Reference Populations')
plt.grid(True, alpha=0.3)
plt.savefig("all_samples_pca.png")
PY
    """
}

process PLOT_ANCESTRY {
    publishDir "${params.outdir}/plots", mode: 'copy'
    conda 'environment.yml'

    input:
    path ancestry_files
    tuple path(mxl_list), path(eur_list), path(afr_list)

    output:
    path "ancestry_plot.png"
    path "all_ancestry.tsv"

    script:
    def ancestryList = ancestry_files instanceof List ? ancestry_files : [ ancestry_files ]
    def ancestryJson = groovy.json.JsonOutput.toJson(ancestryList.collect { it.toString() })

    """
    python3 - <<'PY'
import json
import pandas as pd
import matplotlib.pyplot as plt
import sys

mxl_file = "${mxl_list}"
eur_file = "${eur_list}"
afr_file = "${afr_list}"
input_files = json.loads('''${ancestryJson}''')

if not input_files:
    print("No ancestry files provided", file=sys.stderr)
    sys.exit(1)

# Load reference samples to exclude them
ref_samples = set()
for fpath in [mxl_file, eur_file, afr_file]:
    with open(fpath, 'r') as f:
        for line in f:
            if line.strip():
                ref_samples.add(line.strip())

# Combine all files
dfs = []
for f in input_files:
    df = pd.read_csv(f, sep='\\t')
    dfs.append(df)

if not dfs:
    print("No ancestry files found", file=sys.stderr)
    sys.exit(1)

full_df = pd.concat(dfs, ignore_index=True)
full_df.to_csv("all_ancestry.tsv", sep='\\t', index=False)

# Filter out reference samples
plot_df = full_df[~full_df['IID'].isin(ref_samples)].copy()

if plot_df.empty:
    print("No non-reference samples found for plotting", file=sys.stderr)
    # Fallback to full_df if empty (e.g. testing with only refs)
    plot_df = full_df

# Plotting
plot_df.set_index('IID', inplace=True)

# Sort by the first population column
pop_cols = plot_df.columns
plot_df.sort_values(by=list(pop_cols), ascending=False, inplace=True)

# Calculate figure width: 0.2 inches per sample (more space since fewer samples)
width = max(10, min(200, len(plot_df)*0.2))

ax = plot_df.plot(kind='bar', stacked=True, figsize=(width, 6), width=1.0)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Ancestry Admixture (Samples Only)")
plt.xlabel("Sample")
plt.ylabel("Proportion")
plt.tight_layout()
plt.savefig("ancestry_plot.png", dpi=300)
PY
    """
}


// Filter reference VCF to exome regions before concatenation
process FILTER_REF_CHROM {
    tag { chrom }

    input:
    tuple val(chrom), path(vcf), path(tbi)
    path exome_bed

    output:
    tuple val(chrom), path("ref_exome_${chrom}.vcf.gz"), path("ref_exome_${chrom}.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    
    # Check VCF headers for chr style
    if bcftools view -h ${vcf} | grep -q "##contig=<ID=chr"; then
        # VCF has chr prefix. Ensure BED has chr prefix.
        awk '{if(\$1 !~ /^chr/) print "chr"\$0; else print \$0}' ${exome_bed} > targets.bed
    else
        # VCF does not have chr prefix. Ensure BED doesn't.
        sed 's/^chr//' ${exome_bed} > targets.bed
    fi
    
    bcftools view --threads ${task.cpus} -T targets.bed ${vcf} -Oz -o ref_exome_${chrom}.vcf.gz
    bcftools index -t ref_exome_${chrom}.vcf.gz
    """
}

// Concatenate reference chromosomes into single multi-chrom VCF


// Index sample VCFs
process INDEX_SAMPLE {
    publishDir "${params.outdir}/samples/indexed", mode: 'copy'
    tag { sample_id }

    input:
    tuple val(sample_id), path(sample_vcf)

    output:
    tuple val(sample_id), path(sample_vcf), path("${sample_vcf.name}.tbi")

    script:
    """
    set -euo pipefail
    bcftools index -t -f ${sample_vcf}
    """
}

// Merge sample VCF with aggregated reference VCF


// Convert merged VCF to PLINK and prune variants (mirrors notebook QC)

// Run ADMIXTURE on pruned PLINK set
process RUN_ADMIXTURE {
    publishDir "${params.outdir}/admixture", mode: 'copy'
    tag { sample_id }

    input:
    tuple val(sample_id),
          path(bed),
          path(bim),
          path(fam)

    output:
    tuple val(sample_id),
          path("${sample_id}.Q"),
          path("${sample_id}.P"),
          path(fam)

    script:
    """
    set -euo pipefail
    export OMP_NUM_THREADS=${task.cpus}
    admixture_bin=${params.admixture}
    \$admixture_bin -j${task.cpus} --seed=${params.seed} --cv ${bed} ${params.k} > ${sample_id}_admixture.log
    mv ${bed.baseName}.${params.k}.Q ${sample_id}.Q
    mv ${bed.baseName}.${params.k}.P ${sample_id}.P
    """
}

// Extract sample ancestry percentages from ADMIXTURE output
process SUMMARIZE_Q {
    publishDir "${params.outdir}/ancestry", mode: 'copy'
    tag { sample_id }

    input:
    tuple val(sample_id), path(qfile), path(pfile), path(fam)
    tuple path(mxl_list), path(eur_list), path(afr_list)

    output:
    tuple val(sample_id), path("${sample_id}_ancestry.tsv")

    script:
    """
    set -euo pipefail
    python3 - <<'PY'
import csv
import sys

qfile = "${qfile}"
famfile = "${fam}"
outfile = "${sample_id}_ancestry.tsv"
mxl_file = "${mxl_list}"
eur_file = "${eur_list}"
afr_file = "${afr_list}"

def read_ids(path):
    with open(path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

ref_groups = {
    "MXL": read_ids(mxl_file),
    "EUR": read_ids(eur_file),
    "AFR": read_ids(afr_file),
}

# Read Q file
q_rows = []
with open(qfile, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        q_rows.append([float(x) for x in line.split()])

# Read FAM file
iids = []
with open(famfile, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            iids.append(parts[1])

if len(q_rows) != len(iids):
    sys.exit(f"Q rows ({len(q_rows)}) != FAM IIDs ({len(iids)}), check ADMIXTURE inputs.")

num_k = len(q_rows[0]) if q_rows else 0
if num_k == 0:
    sys.exit("Empty Q file")

if any(len(row) != num_k for row in q_rows):
    sys.exit("Inconsistent number of columns in Q file")

# Calculate mean Q values for each reference group to identify columns
col_means = {pop: [0.0] * num_k for pop in ref_groups}
col_counts = {pop: 0 for pop in ref_groups}

for idx, iid in enumerate(iids):
    for pop, samples in ref_groups.items():
        if iid in samples:
            for k in range(num_k):
                col_means[pop][k] += q_rows[idx][k]
            col_counts[pop] += 1

for pop in ref_groups:
    if col_counts[pop] > 0:
        col_means[pop] = [x / col_counts[pop] for x in col_means[pop]]
    else:
        print(f"Warning: No reference samples found for {pop} in FAM/Q", file=sys.stderr)

# Assign each population to the column where it has the highest average proportion
# We sort populations by their max column value to assign the most distinct ones first
pop_max_vals = []
for pop in ref_groups:
    if col_counts[pop] > 0:
        means = col_means[pop]
        max_val = max(means)
        best_k = means.index(max_val)
        pop_max_vals.append((max_val, pop, best_k))

# Sort by max value descending to assign clear clusters first
pop_max_vals.sort(key=lambda x: x[0], reverse=True)

col_to_pop = {}
for max_val, pop, k in pop_max_vals:
    if k in col_to_pop:
        print(f"Warning: Population {pop} maps to same column {k} as {col_to_pop[k]}", file=sys.stderr)
        continue
    col_to_pop[k] = pop

# Generate labels
labels = []
for k in range(num_k):
    if k in col_to_pop:
        labels.append(col_to_pop[k])
    else:
        labels.append(f"K{k+1}")

print(f"Inferred Labels: {labels}", file=sys.stderr)

# Write output
with open(outfile, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\\t')
    header = ["IID"] + labels
    writer.writerow(header)
    for iid, q in zip(iids, q_rows):
        writer.writerow([iid] + q)
PY
    """
}

process CONVERT_TO_STRUCTURE_ALL {
    tag { sample_id }
    conda 'environment.yml'

    input:
    tuple val(sample_id), path(bed), path(bim), path(fam)
    tuple path(mxl_list), path(eur_list), path(afr_list)

    output:
    tuple val(sample_id), path("${sample_id}.recode.strct_in"), path(bim), path(fam)

    script:
    """
    set -euo pipefail
    # Use plink 1.9 for structure conversion
    ${params.plink} --bfile ${bed.baseName} --recode structure --out ${sample_id}

    # Add PopData and PopFlag columns
    python3 - <<'PY'
import sys

infile = "${sample_id}.recode.strct_in"
outfile = "${sample_id}.recode.strct_in.tmp"
famfile = "${fam}"
mxl_file = "${mxl_list}"
eur_file = "${eur_list}"
afr_file = "${afr_list}"

# Load reference samples
ref_groups = {}
for pop, fpath in [("MXL", mxl_file), ("EUR", eur_file), ("AFR", afr_file)]:
    with open(fpath, 'r') as f:
        ref_groups[pop] = [line.strip() for line in f if line.strip()]

# Map to integers: MXL=1, EUR=2, AFR=3
# POPFLAG: 1 for references, 0 for others
pop_map = {}
flag_map = {}

for pop, samples in ref_groups.items():
    pid = 0
    if pop == "MXL": pid = 1
    elif pop == "EUR": pid = 2
    elif pop == "AFR": pid = 3
    for s in samples:
        pop_map[s] = pid
        flag_map[s] = 1

# Read FAM file to get IIDs
with open(famfile, 'r') as f:
    fam_lines = [l.strip().split() for l in f if l.strip()]
    # FAM: FID IID ...
    iids = [l[1] for l in fam_lines if len(l) >= 2]

with open(infile, 'r') as fin, open(outfile, 'w') as fout:
    lines = fin.readlines()
    
    # Write headers
    fout.write(lines[0])
    fout.write(lines[1])
    
    # Process data lines
    data_lines = lines[2:]
    
    for i, line in enumerate(data_lines):
        parts = line.strip().split()
        if len(parts) < 2: continue
        
        # Use IID from FAM if possible to be sure, or from file
        if i < len(iids):
            iid = iids[i]
        else:
            iid = parts[1]
        
        # Determine PopID and PopFlag
        pid = pop_map.get(iid, 0) # 0 for unknown/samples
        flag = flag_map.get(iid, 0) # 0 for samples
        
        genotypes = parts[2:]
        
        # Structure format: Label PopID PopFlag Genotypes
        new_line = f"{iid} {pid} {flag} " + " ".join(genotypes) + "\\n"
        fout.write(new_line)

PY
    mv ${sample_id}.recode.strct_in.tmp ${sample_id}.recode.strct_in
    """
}

process RUN_STRUCTURE_ALL {
    publishDir "${params.outdir}/structure", mode: 'copy'
    tag { sample_id }
    conda 'environment.yml'

    input:
    tuple val(sample_id), path(infile), path(bim), path(fam)

    output:
    tuple val(sample_id), path("${sample_id}_structure_out_f"), path(fam)

    script:
    """
    set -euo pipefail
    
    # Count individuals and loci
    num_inds=\$(wc -l < ${fam})
    num_loci=\$(wc -l < ${bim})
    
    # Create mainparams
    cp ${infile} input.txt
    
    cat <<EOF > mainparams
#define MAXPOPS    ${params.k}
#define BURNIN    ${params.burnin}
#define NUMREPS   ${params.numreps}
#define INFILE    input.txt
#define OUTFILE   output.txt
#define NUMINDS    \$num_inds
#define NUMLOCI    \$num_loci
#define PLOIDY     2
#define MISSING    -9
#define ONEROWPERIND 1
#define LABEL     1
#define POPDATA   1
#define POPFLAG   1
#define LOCDATA   0
#define PHENOTYPE 0
#define EXTRACOLS 0
#define MARKERNAMES 1
#define RECESSIVEALLELES 0
#define MAPDISTANCES 1
#define USEPOPINFO 1
#define PHASED   0
#define MARKOVPHASE 0
EOF

    touch extraparams
    
    ${params.structure}
    
    mv output.txt_f ${sample_id}_structure_out_f
    """
}

process SUMMARIZE_STRUCTURE_ALL {
    publishDir "${params.outdir}/ancestry", mode: 'copy'
    tag { sample_id }
    conda 'environment.yml'

    input:
    tuple val(sample_id), path(struct_out), path(fam)
    tuple path(mxl_list), path(eur_list), path(afr_list)

    output:
    tuple val(sample_id), path("${sample_id}_structure_ancestry.tsv")

    script:
    """
    set -euo pipefail
    python3 - <<'PY'
import csv
import sys

sample_id = "${sample_id}"
struct_file = "${struct_out}"
outfile = f"{sample_id}_structure_ancestry.tsv"
mxl_file = "${mxl_list}"
eur_file = "${eur_list}"
afr_file = "${afr_list}"

# Load reference lists
ref_map = {} # IID -> PopLabel
for pop, fpath in [("MXL", mxl_file), ("EUR", eur_file), ("AFR", afr_file)]:
    with open(fpath, 'r') as f:
        for line in f:
            if line.strip():
                ref_map[line.strip()] = pop

# Parse Structure Output
data = {} # IID -> [Q1, Q2, Q3]
pop_cluster_map = {} # PopID (int) -> [Prob_C1, Prob_C2, ...]

with open(struct_file, 'r') as f:
    lines = f.readlines()

# 1. Parse "Given Pop" table
in_table = False
for line in lines:
    if "Given    Inferred Clusters" in line:
        in_table = True
        continue
    if in_table:
        if "----------------" in line:
            in_table = False
            continue
        if ":" in line:
            parts = line.strip().split(':')
            try:
                pop_id = int(parts[0].strip())
                # parts[1] has probs and count
                vals = parts[1].strip().split()
                # The last one is count, remove it
                probs = [float(x) for x in vals[:-1]]
                pop_cluster_map[pop_id] = probs
            except:
                pass

if not pop_cluster_map:
    # Fallback if table not found (e.g. Unsupervised)
    # Use previous logic or fail?
    # For now, let's assume Supervised mode produces this table.
    # If not, we might need the old parser.
    pass

# Determine K
num_k = 3
if 1 in pop_cluster_map:
    num_k = len(pop_cluster_map[1])

# 2. Parse Individuals
in_data = False
for line in lines:
    if "Inferred ancestry of individuals:" in line:
        in_data = True
        continue
    if in_data:
        if "Estimated Allele Frequencies" in line:
            break
        if not line.strip(): continue
        if "Label" in line and "(%Miss)" in line: continue
        
        parts = line.strip().split()
        if not parts[0].isdigit(): continue
        
        if ":" not in parts: continue
        colon_idx = parts.index(":")
        
        iid = parts[1]
        
        try:
            pop_id = int(parts[colon_idx-1])
        except:
            continue
            
        if pop_id == 0:
            # Target Sample: Read Q values directly
            # After colon: Q1 Q2 Q3
            # Be careful with parsing, sometimes there are extra spaces
            q_part = parts[colon_idx+1:]
            # Take first K values
            q_vals = [float(x) for x in q_part[:num_k]]
            data[iid] = q_vals
        else:
            # Reference Sample: Use the table mapping
            if pop_id in pop_cluster_map:
                data[iid] = pop_cluster_map[pop_id]
            else:
                data[iid] = [0.0] * num_k

# Determine Labels
cluster_pops = {} # K -> [Pops]
for pid, probs in pop_cluster_map.items():
    if pid == 0: continue
    if not probs: continue
    max_p = max(probs)
    best_k = probs.index(max_p)
    if max_p > 0.4: # Threshold
        if best_k not in cluster_pops: cluster_pops[best_k] = []
        pname = {1:"MXL", 2:"EUR", 3:"AFR"}.get(pid, "?")
        cluster_pops[best_k].append(pname)

final_labels = []
for k in range(num_k):
    if k in cluster_pops:
        final_labels.append("-".join(cluster_pops[k]))
    else:
        final_labels.append(f"K{k+1}")

# Write Output
with open(outfile, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\\t')
    writer.writerow(["IID"] + final_labels)
    for iid, q in data.items():
        writer.writerow([iid] + q)
PY
    """
}

// Download the 1000 Genomes panel file
process DOWNLOAD_PANEL {
    publishDir "${params.outdir}/reference/meta", mode: 'copy'

    output:
    path "integrated_call_samples_v3.20130502.ALL.panel"

    script:
    """
    set -euo pipefail
    wget -O integrated_call_samples_v3.20130502.ALL.panel "${params.panel_url}"
    """
}

// Create sample lists for MXL, EUR, AFR
process CREATE_REF_LISTS {
    publishDir "${params.outdir}/reference/meta", mode: 'copy'

    input:
    path panel

    output:
    tuple path("mxl_40.samples"), path("eur_40.samples"), path("afr_40.samples"), path("all_ref.samples")

    script:
    """
    set -euo pipefail
    
    # MXL
    awk 'BEGIN{FS="\\t"; c=0} \$2=="MXL" {print \$1; c++; if(c==40) exit}' ${panel} > mxl_40.samples
    
    # EUR (Super-pop)
    awk 'BEGIN{FS="\\t"; c=0} \$3=="EUR" {print \$1; c++; if(c==40) exit}' ${panel} > eur_40.samples
    
    # AFR (Super-pop)
    awk 'BEGIN{FS="\\t"; c=0} \$3=="AFR" {print \$1; c++; if(c==40) exit}' ${panel} > afr_40.samples
    
    cat mxl_40.samples eur_40.samples afr_40.samples > all_ref.samples
    """
}

// Download 1KG VCF per chromosome
process DOWNLOAD_1KG_VCF {
    publishDir "${params.outdir}/reference/raw_vcfs", mode: 'copy'
    tag { chrom }

    input:
    val chrom

    output:
    tuple val(chrom), path("${chrom}.1kg.vcf.gz"), path("${chrom}.1kg.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    # Handle chr prefix if present in chrom val
    c=\$(echo ${chrom} | sed 's/^chr//')
    
    vcf_url="${params.vcf_base_url}/ALL.chr\${c}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    tbi_url="\${vcf_url}.tbi"
    
    wget -O chr\${c}.1kg.vcf.gz "\$vcf_url"
    wget -O chr\${c}.1kg.vcf.gz.tbi "\$tbi_url"
    """
}

// Subset VCF to specific population
process SUBSET_REF_VCF {
    publishDir "${params.outdir}/reference/subset_vcfs", mode: 'copy'
    tag { chrom }

    input:
    tuple val(chrom), path(vcf), path(tbi), path(ref_list)

    output:
    tuple val(chrom), path("ref_${chrom}.vcf.gz"), path("ref_${chrom}.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    c=\$(echo ${chrom} | sed 's/^chr//')
    
    bcftools view --threads ${task.cpus} -S ${ref_list} -Oz -o ref_${chrom}.vcf.gz ${vcf}
    bcftools index -t ref_${chrom}.vcf.gz
    """
}
