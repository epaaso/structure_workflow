import gzip
import sys

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    merged = []
    current_start, current_end = intervals[0]
    for i in range(1, len(intervals)):
        next_start, next_end = intervals[i]
        if next_start <= current_end:
            current_end = max(current_end, next_end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = next_start, next_end
    merged.append((current_start, current_end))
    return merged

intervals_by_chrom = {}

with gzip.open('ncbiRefSeqCurated.txt.gz', 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[2]
        
        # Filter for standard chromosomes if needed, or keep all
        if "_" in chrom: continue # Skip random contigs for simplicity
        
        exon_starts = parts[9].split(',')
        exon_ends = parts[10].split(',')
        
        for s, e in zip(exon_starts, exon_ends):
            if s and e:
                try:
                    start = int(s)
                    end = int(e)
                    if chrom not in intervals_by_chrom:
                        intervals_by_chrom[chrom] = []
                    intervals_by_chrom[chrom].append((start, end))
                except ValueError:
                    continue

with open('exome_hg38.bed', 'w') as out:
    # Sort chromosomes naturally if possible, or just alphabetically
    for chrom in sorted(intervals_by_chrom.keys()):
        merged = merge_intervals(intervals_by_chrom[chrom])
        for start, end in merged:
            out.write(f"{chrom}\t{start}\t{end}\n")
