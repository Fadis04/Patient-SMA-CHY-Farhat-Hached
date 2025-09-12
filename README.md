SMN1/SMN2 TARGETED NANOPORE SEQUENCING ANALYSIS PIPELINE

CHU Farhat Hached, Sousse – Tunisia
Comprehensive Molecular Diagnosis of Spinal Muscular Atrophy (SMA) Using Long-Read Sequencing

📌 PROJECT OVERVIEW

This repository contains sequencing data and analysis pipelines for the molecular diagnosis of Spinal Muscular Atrophy (SMA) patients from CHU Farhat Hached, Sousse (Tunisia).

We employed Oxford Nanopore Technologies (ONT) targeted sequencing of the SMN1/SMN2 locus, combined with advanced bioinformatics, to:

Detect pathogenic variants

Identify copy number variations (CNVs)

Characterize SMN1–SMN2 hybridization patterns

🧬 KEY FINDINGS

8 SMA patients analyzed using targeted Nanopore sequencing

Pathogenic deletions identified in the SMN1 gene

Detection of the SMN1/SMN2 discriminatory variant (c.840C>T at chr5:70,951,946)

Use of a masked reference approach to improve SMN1-specific variant calling

Identification of pathogenic variants not detected by conventional methods

📂 REPOSITORY STRUCTURE
Patient-SMA-CHY-Farhat-Hached/
├── analysis_scripts/
│   ├── smn_analysis_pipeline.sh    # Main analysis pipeline
│   └── variant_annotation.py       # Variant annotation
├── processed_data/
│   ├── masked_reference/           # SMN2-masked reference
│   ├── variant_calls/              # Clair3 VCF outputs
│   └── coverage_analysis/          # Depth statistics
├── raw_data/
│   └── les_barcodes_par_nanopore/  # Barcoded patient data
├── results/
│   ├── pathogenic_variants.csv     # Curated pathogenic variants
│   ├── smn1_smn2_discrimination/   # PSV analysis
│   └── quality_metrics/            # QC reports
└── references/
    ├── SMN1_exons.bed              # SMN1 target regions
    ├── SMN2_exons.bed              # SMN2 target regions
    └── GRCh38_chr5.fa              # Reference genome

⚙️ WORKFLOW
Raw Data Preparation

Target region: SMN1 gene (chr5:70,910,000–70,960,000)

Platform: ONT PromethION

Samples: 8 patients (barcode08–barcode15)

Coverage: >100X per patient

Analysis Pipeline

Reference preparation with SMN2 masking (chr5:70,048,000–70,078,000)

Read alignment using minimap2

Variant calling with Clair3

Structural variant detection

SMN1/SMN2 discrimination analysis

Pathogenic variant validation

Masked Reference Strategy

The SMN2 region was replaced with N’s to force read alignment to SMN1, improving variant calling specificity and reducing mapping ambiguity.

🔬 GOALS

Provide a robust molecular diagnosis workflow for SMA using long-read sequencing

Detect SMN1-specific pathogenic variants (including deletions in homopolymer regions)

Compare Nanopore-based results with conventional diagnostic methods

Build an open bioinformatics resource for SMA research and clinical genomics

👨‍💻 AUTHOR

Fadi Slimi – Bioinformatician
📧 Email: fadi.slimi@insat.ucar.tn

🔗 LinkedIn: www.linkedin.com/in/fadi-slimi
