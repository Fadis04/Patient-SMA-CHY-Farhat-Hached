SMN1/SMN2 Targeted Nanopore Sequencing Analysis Pipeline

CHU Farhat Hached, Sousse – Tunisia
Comprehensive Molecular Diagnosis of Spinal Muscular Atrophy Using Long-Read Sequencing

📖 Project Overview

This repository provides data and analysis pipelines for the molecular diagnosis of Spinal Muscular Atrophy (SMA) in patients from CHU Farhat Hached, Sousse, Tunisia.
The project leverages Oxford Nanopore Technologies (ONT) targeted sequencing of the SMN1/SMN2 locus, combined with advanced bioinformatics workflows, to:

Detect pathogenic variants

Identify copy number variations (CNVs)

Characterize SMN1–SMN2 hybridization patterns

🧬 Key Findings

8 SMA patients analyzed using targeted ONT sequencing

Pathogenic deletions identified in the SMN1 gene

Detection of the SMN1/SMN2 discriminatory variant at chr5:70,951,946

Implementation of a masked reference strategy to improve variant calling accuracy

Identification of pathogenic variants missed by conventional approaches

🧪 Methods
🔹 Targeted Sequencing Design

Target region: SMN1 gene (chr5:70,910,000–70,960,000)

Platform: ONT PromethION

Samples: 8 patients (barcode08–barcode15)

Coverage: >100X per patient

🔹 Bioinformatics Pipeline
# Main analysis workflow
1. Reference preparation with SMN2 masking (chr5:70,048,000–70,078,000)  
2. Read alignment using minimap2  
3. Variant calling with Clair3  
4. Structural variant detection  
5. SMN1/SMN2 discrimination analysis  
6. Pathogenic variant validation  

🔹 Masked Reference Strategy

To enhance SMN1-specific detection, the SMN2 region was masked with N’s, forcing reads to align preferentially to SMN1 and reducing mapping ambiguity.

🔬 Key Results
1. SMN1/SMN2 Discriminatory Variant – c.840C>T

Position: chr5:70,951,946

Patient barcode08:

Total reads: 702

C (SMN1-specific): 697 (99%)

T (SMN2-specific): 5 (1%)

Deletions: 11

✅ Confirms the canonical SMN1 (C) vs SMN2 (T) marker with strong enrichment for SMN1 reads.

2. Pathogenic Deletion – c.584del (p.Pro195LeufsTer18)

Position: chr5:70,942,825

Variant: SMN1(NM_001297715.1):c.584del

Consequence: Frameshift → Pathogenic (VarSome)

Patient barcode08:

Total reads: 137

C reference: 136 (99%)

Deletion events: 91 (66%)

⚠️ Automated callers showed variable detection due to the homopolymer context (CCCCC), highlighting the need for manual review.

3. Additional Pathogenic Variants

(Analysis in progress – results will be added soon)

🧫 Comparative Analysis
Variant Detection Challenges

Homopolymer artifacts: c.584del inconsistently detected in CCCCC repeats

Structural variants: Larger deletions required IGV manual validation

Reference bias: Impacted variant calling in identical SMN1/SMN2 regions

Literature Context

“Long-read sequencing enables phased variant detection and improves SMA carrier screening” – Genome Medicine, 2025

“Comprehensive SMN1/SMN2 variant characterization requires multi-method approaches” – Human Mutation, 2023

📊 Repository Structure
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

🚀 Usage
Run Complete Analysis
bash scripts/smn_analysis_pipeline.sh

Process Individual Patients
python scripts/variant_calling.py --patient barcode08 --threads 16

Data Access
tar -xzf SMA_analysis/align.tar.gz

👥 Authors

Fadi Slimi – Bioinformatician

📧 Email: fadi.slimi@insat.ucar.tn

🔗 LinkedIn: www.linkedin.com/in/fadi-slimi
GitHub Issues: Project Issues

All patient data is anonymized.

Last updated: September 2024
Analysis pipeline version: 2.0
