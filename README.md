SMN1/SMN2 Targeted Nanopore Sequencing Analysis Pipeline
CHU Farhat Hached Sousse, Tunisia
Comprehensive Molecular Diagnosis of Spinal Muscular Atrophy Using Long-Read Sequencing

📖 Project Overview
This repository contains data and analysis pipelines for the molecular diagnosis of Spinal Muscular Atrophy (SMA) patients from CHU Farhat Hached Sousse, Tunisia. The study employs Oxford Nanopore Technologies (ONT) targeted sequencing of the SMN1/SMN2 locus combined with advanced bioinformatics analysis to detect pathogenic variants, copy number variations, and characterize SMN1-SMN2 hybridization patterns.

🧬 Key Findings
8 SMA patients analyzed using targeted Nanopore sequencing

Pathogenic deletions identified in the SMN1 gene

Discriminatory variant analysis at position chr5:70,951,946 (SMN1 vs SMN2)

Masked reference approach to improve variant calling accuracy

Multiple pathogenic variants undetected by conventional methods

🧪 Methods
Targeted Sequencing Design
Target Region: SMN1 gene (chr5:70,910,000-70,960,000)

Technology: Oxford Nanopore PromethION

Barcodes: 8 patients (barcode08-barcode15)

Coverage: >100X per patient

Bioinformatics Pipeline
bash
# Main analysis steps
1. Reference preparation with SMN2 masking (chr5:70,048,000-70,078,000)
2. Read alignment using minimap2
3. Variant calling with Clair3
4. Structural variant detection
5. SMN1/SMN2 discrimination analysis
6. Pathogenic variant validation
Masked Reference Strategy
To improve SMN1-specific variant detection, we created a masked reference genome where the SMN2 region (chr5:70,048,000-70,078,000) was replaced with N's, forcing reads to align to SMN1 and reducing mapping ambiguity.

🔬 Key Results
1. SMN1/SMN2 Discriminatory Variant (chr5:70,951,946)
Patient barcode08 results:

text
Position: chr5:70,951,946
Total reads: 702
C (SMN1-specific): 697 (99%) 
T (SMN2-specific): 5 (1%)
Deletions: 11
This position represents the canonical SMN1 (C) vs SMN2 (T) discriminatory variant (c.840C>T). The 99% C support indicates successful SMN1-targeted sequencing with minimal SMN2 cross-amplification.

2. Pathogenic Deletion Detection (chr5:70,942,825)
Pathogenic Variant:

Position: chr5:70,942,825

Variant: SMN1(NM_001297715.1):c.584del

Consequence: p.(Pro195LeufsTer18)

Pathogenicity: Pathogenic (VarSome)

Patient barcode08 quantification:

text
Total reads: 137
C reference: 136 (99%)
Deletion events: 91 (66%)
This frameshift deletion was consistently detected across multiple patients but showed variable detection efficiency in automated variant callers, likely due to the homopolymer context (CCCCC).

3. Additional Pathogenic Findings
??

🧫 Comparative Analysis
Variant Detection Challenges
Our analysis revealed several limitations in current variant detection pipelines:

Homopolymer Artifacts: Variants in repetitive regions (e.g., c.584del in CCCCC context) showed inconsistent detection

Structural Variants: Large deletions required manual IGV validation despite high coverage

SMN1/SMN2 Discrimination: Reference bias affected variant calling in identical regions

Comparison with Existing Literature
Our findings align with recent advancements in SMA molecular diagnostics:

"Long-read sequencing enables phased variant detection and improves SMA carrier screening" - Genome Medicine 2025

*"Comprehensive SMN1/SMN2 variant characterization requires multi-method approaches"* - Human Mutation 2023

📊 Repository Structure
text
Patient-SMA-CHY-Farhat-Hached/
├── analysis_scripts/
│   ├── smn_analysis_pipeline.sh    # Main analysis pipeline
│   └── variant_annotation.py       # Variant annotation
├── processed_data/
│   ├── masked_reference/           # SMN2-masked reference
│   ├── variant_calls/              # Clair3 VCF outputs
│   └── coverage_analysis/          Depth statistics
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
Pipeline Execution
bash
# Run complete analysis
bash scripts/smn_analysis_pipeline.sh

# Process individual patients
python scripts/variant_calling.py --patient barcode08 --threads 16
Data Access
Compressed alignment files require extraction:

bash
tar -xzf SMA_analysis/align.tar.gz
👥 Authors
Fadi Slimi - Bioinformatician 

📧 Contact
For questions regarding this dataset or analysis pipeline:

Email: fadi.slimi@insat.ucar.tn

GitHub Issues: Project Issues

All patient data is anonymized.

Last updated: September 2024
Analysis pipeline version: 2.0

