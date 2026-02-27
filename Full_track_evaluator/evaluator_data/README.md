# evaluator_data/


This directory contains the script and intermediary files needed to create the sequence file: `enformer_borzoi_test_seqs.parquet'`

These are the only files required to run the pre-built container and can be downloaded from Zenodo (LINK):
- `enformer_borzoi_test_seqs.parquet` — genomic sequences for each test region
- `ENCFF972GVB.bigWig` — measured K562 DNase-seq signal (can be downloaded from here: https://www.encodeproject.org/experiments/ENCSR000EKS/)

The remaining files are included so you can follow along with how the input parquet was generated from scratch (see `sequence_design.py`).

---

## Directory Structure

```
evaluator_data/
├── enformer_borzoi_test_seqs.parquet         # ← Required to run evaluator
├── ENCFF972GVB.bigWig                         # ← Required to run evaluator
├── enformer_human_hg38_sequences.bed          # Enformer train/test/val BED file (hg38)
├── sequences_human_borzoi.bed                 # Borzoi train/test/val BED file (hg38)
├── enformer_human_hg38_sequences_test196kb.bed         # Enformer test regions extended to 196kb
├── enformer_human_hg38_sequences_test196kb_sorted.bed  # Sorted
├── enformer_human_hg38_sequences_test196kb_merged.bed  # Merged (85 regions)
├── sequences_human_borzoi_test.bed            # Borzoi test regions (fold3)
├── sequences_human_borzoi_test_sorted.bed     # Sorted
├── sequences_human_borzoi_test_merged.bed     # Merged (46 regions)
└── enformer_borzoi_test_merged.bed            # Intersection of Enformer and Borzoi test regions
```

---

## How the parquet was generated

The input parquet is produced by `sequence_design.py` in three steps.

### Step 1 — Filter BED files

`filter_enformer_bed()` keeps only `label == "test"` entries from the Enformer BED file and extends each region by 32,768 bp on each side, giving 196,608 bp Enformer training windows.

`filter_borzoi_bed()` keeps only `fold3` entries from the Borzoi BED file (fold 3 is the Borzoi test set).

```bash
python sequence_design.py  # run Steps 1 only
```

### Step 2 — Sort, merge, and intersect with bedtools

```bash
bedtools sort  -i enformer_human_hg38_sequences_test196kb.bed \
  | bedtools merge > enformer_human_hg38_sequences_test196kb_merged.bed
# 85 regions

bedtools sort  -i sequences_human_borzoi_test.bed \
  | bedtools merge > sequences_human_borzoi_test_merged.bed
# 46 regions

bedtools intersect \
  -a enformer_human_hg38_sequences_test196kb_merged.bed \
  -b sequences_human_borzoi_test_merged.bed \
  > enformer_borzoi_test_merged.bed
```

### Step 3 — Fetch sequences from hg38

Uncomment `pull_sequences_from_bed()` at the bottom of `sequence_design.py` and re-run. This fetches the DNA sequence for each region in `enformer_borzoi_test_merged.bed` from the hg38 FASTA and saves the result to `enformer_borzoi_test_seqs.parquet`.

The parquet contains columns: `chrom`, `start`, `end`, `sequence`, `region`.