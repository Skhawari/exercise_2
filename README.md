# Snakemake Workflow – Exercise 2 (Applied Sequence Analysis)

This repository contains a reproducible Snakemake workflow for automating common NGS preprocessing steps using Bowtie2 and Samtools. The workflow was developed as part of the "Applied Sequence Analysis" course at FU Berlin.

## 💡 What this workflow does

For each paired-end FASTQ sample, it performs:

1. **Indexing** of the reference genome with `bowtie2-build`
2. **Mapping** with `bowtie2`
3. **Conversion** from SAM to BAM (`samtools view`)
4. **Sorting** of BAM files (`samtools sort`)
5. **Indexing** of sorted BAM files (`samtools index`)
6. **Read count statistics** with `samtools idxstats`
7. **Aggregation** of per-sample read counts into a summary table

## 📁 Repository structure

```
.
├── config/             # Contains config.yaml with paths and parameters
├── resources/          # Contains reference genome (reference.fa)
├── results/            # Output data: sam, bam, sorted bam, idxstats
├── workflow/
│   ├── envs/           # Conda environment YAMLs (bowtie2.yaml, samtools.yaml)
│   ├── notebooks/      # Optional Jupyter notebooks
│   ├── report/         # Workflow documentation and reports
│   ├── rules/          # Modular Snakemake rule files
│   ├── scripts/        # Helper scripts (e.g., aggregate_idxstats.py)
│   ├── Snakefile       # Main Snakemake workflow
│   └── samples.csv     # Sample information (paths to FASTQ files)
├── .gitignore
├── LICENSE.md
└── README.md
```

## ⚙️ How to run

First, activate conda and run:

```bash
conda config --set channel_priority strict
snakemake --use-conda --cores 4
```

## ✅ Output

Final summary file (after complete run):

```
results/summary/idxstats_summary.tsv
```

## 🔬 Developed by

- Dr. Sandro Andreotti – Lecturer, Freie Universität Berlin
- Michael Riethmüller – Master's Student in Bioinformatics, Freie Universität Berlin
- Sajjad Khawari – Master's Student in Bioinformatics, Freie Universität Berlin