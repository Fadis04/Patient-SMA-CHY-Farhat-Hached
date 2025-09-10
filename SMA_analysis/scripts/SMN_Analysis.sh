#!/usr/bin/env bash
# ==============================================================================
# SMN1/SMN2 Targeted Nanopore Sequencing Analysis Pipeline
# Version: 2.0
# Author: Fadi Slimi
# Institution: CHU Farhat Hached Sousse, Tunisia
# 
# Description: Comprehensive pipeline for SMN1/SMN2 variant calling and analysis
# using Oxford Nanopore Technologies data with Clair3.
# ==============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# ==============================================================================
# CONFIGURATION SECTION - EDIT THESE VARIABLES
# ==============================================================================

# Sample and Directory Configuration
export PATIENT="barcode08"
export WD="/home/wassef/Patients_Moez_Gribaa/analysis_${PATIENT}"
export READS_DIR="/home/wassef/Patients_Moez_Gribaa/les_barcodes_par_nanopore"
export REF_DIR="/home/wassef/Patients_Moez_Gribaa/ref_gen"

# Reference Files
export REF_FASTA="${REF_DIR}/GRCh38_chr5.fa"
export SMN1_BED="/home/wassef/Patients_Moez_Gribaa/SMN1_exons.bed"
export SMN2_BED="/home/wassef/Patients_Moez_Gribaa/SMN2_exons.bed"

# SMN2 Region for Masking (GRCh38 coordinates)
export SMN2_CHR="5"
export SMN2_START=70048000
export SMN2_END=70078000

# Region of Interest (SMN1/SMN2 locus)
export ROI="5:70910000-70960000"
export EX7_START=70946065  # Exon 7 coordinates
export EX7_END=70946176

# Analysis Parameters
export PLOIDY=4            # Adjust according to total SMN copies
export THREADS=16
export CLAIR3_MODELS="/home/wassef/Patients_Moez_Gribaa/Clair3/models/ont"

# Derived Paths
export READS="${READS_DIR}/${PATIENT}/merged_${PATIENT}.fastq"
export REF_MASKED="${WD}/ref/GRCh38_chr5_masked_smn2.fa"
export REF_DICT="${WD}/ref/GRCh38_chr5_masked_smn2.dict"

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

print_header() {
    echo "================================================================"
    echo "$1"
    echo "================================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $2"
}

check_dependencies() {
    local tools=("minimap2" "samtools" "picard" "python3" "bcftools")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "ERROR: $tool not found in PATH"
            exit 1
        fi
    done
    echo "[✓] All dependencies found"
}

setup_directory_structure() {
    local subdirs=(
        "fastq" "ref" "align" "phasing" "clair3_out" 
        "sv" "outtables" "scripts" "logs" "qc"
    )
    
    for dir in "${subdirs[@]}"; do
        mkdir -p "${WD}/${dir}"
    done
    echo "[✓] Directory structure created in ${WD}"
}

prepare_masked_reference() {
    print_header "REFERENCE PREPARATION" "Creating masked reference for SMN2 region"
    
    # Index original reference
    if [[ ! -f "${REF_FASTA}.fai" ]]; then
        samtools faidx "${REF_FASTA}"
    fi

    # Create masked reference using Python
    python3 - <<PYTHON_SCRIPT
from Bio import SeqIO
from Bio.Seq import Seq

ref_path = "${REF_FASTA}"
output_path = "${REF_MASKED}"
mask_chr = "${SMN2_CHR}"
start_pos = ${SMN2_START} - 1  # Convert to 0-based
end_pos = ${SMN2_END}

print(f"Masking {mask_chr}:{start_pos+1}-{end_pos} in reference")

records = []
for record in SeqIO.parse(ref_path, "fasta"):
    if record.id == mask_chr:
        print(f"Processing chromosome {mask_chr}")
        seq_list = list(record.seq)
        # Mask the SMN2 region
        for i in range(start_pos, end_pos):
            if i < len(seq_list):
                seq_list[i] = 'N'
        masked_seq = Seq(''.join(seq_list))
        record.seq = masked_seq
    records.append(record)

SeqIO.write(records, output_path, "fasta")
print(f"Masked reference saved to: {output_path}")
PYTHON_SCRIPT

    # Index masked reference
    samtools faidx "${REF_MASKED}"
    
    # Create sequence dictionary if it doesn't exist
    if [[ ! -f "${REF_DICT}" ]]; then
        picard CreateSequenceDictionary \
            R="${REF_MASKED}" \
            O="${REF_DICT}" \
            > "${WD}/logs/picard_dict.log" 2>&1
    fi
}

run_alignment() {
    print_header "READ ALIGNMENT" "Mapping reads to masked reference"
    
    local out_bam="${WD}/align/${PATIENT}.GRCh38.sorted.bam"
    local roi_bam="${WD}/align/${PATIENT}.SMN_ROI.bam"
    local roi_rg_bam="${WD}/align/${PATIENT}.SMN_ROI.RG.bam"

    # Map reads with minimap2
    echo "[${PATIENT}] Mapping reads with minimap2..."
    minimap2 -ax map-ont --secondary=yes -t ${THREADS} "${REF_MASKED}" "${READS}" \
        | samtools sort -@8 -m 2G -o "${out_bam}" -
    
    samtools index "${out_bam}"

    # Extract Region of Interest
    echo "[${PATIENT}] Extracting ROI: ${ROI}"
    samtools view -b "${out_bam}" "${ROI}" -o "${roi_bam}"
    samtools index "${roi_bam}"

    # Add read groups
    picard AddOrReplaceReadGroups \
        I="${roi_bam}" \
        O="${roi_rg_bam}" \
        RGID=1 \
        RGLB="lib1" \
        RGPL="ONT" \
        RGPU="unit1" \
        RGSM="${PATIENT}" \
        > "${WD}/logs/picard_rg.log" 2>&1
    
    samtools index "${roi_rg_bam}"
    echo "[✓] Alignment completed: ${roi_rg_bam}"
}

run_variant_calling() {
    print_header "VARIANT CALLING" "Running Clair3 for variant detection"
    
    local clair3_out="${WD}/clair3_out/${PATIENT}"
    local roi_rg_bam="${WD}/align/${PATIENT}.SMN_ROI.RG.bam"

    export PATH="/home/wassef/Patients_Moez_Gribaa/Clair3:${PATH}"

    ./run_clair3.sh \
        --bam_fn="${roi_rg_bam}" \
        --ref_fn="${REF_MASKED}" \
        --threads=8 \
        --platform="ont" \
        --model_path="${CLAIR3_MODELS}" \
        --output="${clair3_out}" \
        --sample_name="${PATIENT}" \
        --gvcf \
        --print_ref_calls \
        --snp_min_af=0.1 \
        --indel_min_af=0.1 \
        --include_all_ctgs \
        --chr="5" \
        > "${WD}/logs/clair3.log" 2>&1

    echo "[✓] Variant calling completed: ${clair3_out}"
}

calculate_coverage() {
    print_header "COVERAGE ANALYSIS" "Calculating exon coverage metrics"
    
    local roi_rg_bam="${WD}/align/${PATIENT}.SMN_ROI.RG.bam"
    local depth_file="${WD}/outtables/${PATIENT}_exon7.depth.txt"
    local coverage_summary="${WD}/qc/${PATIENT}_coverage_summary.txt"

    # Calculate depth for exon 7
    echo "[${PATIENT}] Calculating exon 7 coverage..."
    samtools depth -r "5:${EX7_START}-${EX7_END}" "${roi_rg_bam}" > "${depth_file}"

    # Calculate coverage statistics
    awk '{
        sum += $3; 
        if ($3 > 0) covered++; 
        total++;
        if ($3 < 10) low_cov++;
    } END {
        avg = sum / total;
        coverage_pct = (covered / total) * 100;
        low_cov_pct = (low_cov / total) * 100;
        print "Sample: ${PATIENT}";
        print "Region: Exon 7 (5:${EX7_START}-${EX7_END})";
        print "Average coverage: " avg;
        print "Coverage percentage: " coverage_pct "%";
        print "Bases with <10x coverage: " low_cov " (" low_cov_pct "%)";
        print "Total bases: " total;
    }' "${depth_file}" > "${coverage_summary}"

    cat "${coverage_summary}"
    echo "[✓] Coverage analysis saved: ${coverage_summary}"
}

generate_qc_report() {
    print_header "QUALITY CONTROL" "Generating QC report"
    
    local qc_report="${WD}/qc/${PATIENT}_qc_report.txt"
    local roi_rg_bam="${WD}/align/${PATIENT}.SMN_ROI.RG.bam"

    {
        echo "SMN1/SMN2 Analysis QC Report"
        echo "============================"
        echo "Sample: ${PATIENT}"
        echo "Date: $(date)"
        echo ""
        
        echo "Read Statistics:"
        echo "----------------"
        samtools flagstat "${roi_rg_bam}" | head -5
        echo ""
        
        echo "Mapping Quality:"
        echo "----------------"
        samtools stats "${roi_rg_bam}" | grep "^SN" | grep -E "(raw total sequences|error rate|average length|insert size)"
        echo ""
        
        echo "Coverage Summary:"
        echo "-----------------"
        cat "${WD}/qc/${PATIENT}_coverage_summary.txt"
        
    } > "${qc_report}"

    echo "[✓] QC report generated: ${qc_report}"
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

main() {
    print_header "SMN ANALYSIS PIPELINE" "Starting analysis for ${PATIENT}"
    
    # Initialize
    check_dependencies
    setup_directory_structure
    
    # Analysis steps
    prepare_masked_reference
    run_alignment
    run_variant_calling
    calculate_coverage
    generate_qc_report
    
    print_header "PIPELINE COMPLETED" "Analysis finished successfully for ${PATIENT}"
    echo "Results available in: ${WD}"
    echo "VCF files: ${WD}/clair3_out/${PATIENT}"
    echo "QC reports: ${WD}/qc/"
}

# Execute main function and log output
main 2>&1 | tee "${WD}/logs/pipeline_${PATIENT}_$(date +%Y%m%d_%H%M%S).log"
