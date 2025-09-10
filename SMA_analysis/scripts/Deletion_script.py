import pysam
import csv
import glob
import os
from statistics import mean

# ================================
# Paths and Parameters
# ================================
bam_folder = "/home/wassef/Patients_Moez_Gribaa/SMA_analysis/align/"
bam_files = glob.glob(os.path.join(bam_folder, "*SMN_ROI.RG.bam"))

roi_chr = "5"
roi_start = 70910000
roi_end = 70960000

ref_fasta = "/home/wassef/Patients_Moez_Gribaa/SMA_analysis/ref/GRCh38_chr5_masked_smn2.fa"

AF_threshold = 0.2   # allele proportion threshold (0–1)
MinDepth = 10        # minimum total depth to consider

output_file = "/home/wassef/Patients_Moez_Gribaa/SMA_analysis/SMN_DEL_quantitative.csv"

# ================================
# Functions
# ================================
def compute_max_depth(bam, chr, start, end):
    """Compute maximum coverage in ROI for normalization."""
    depths = [p.n for p in bam.pileup(chr, start, end, truncate=True, stepper="all")]
    return max(depths) if depths else 1  # avoid division by zero

# ================================
# Main Script
# ================================
with open(output_file, "w", newline="") as csvfile:
    fieldnames = [
        "Patient", "Chromosome", "Position", "Ref_Base", "Variant_Type",
        "Total_Depth", "DEL_Count", "+_Strand", "-_Strand",
        "Allele_Frequency(%)", "AF_proportion", "Coverage_proportion", "Zygosity_proportion"
    ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for bam_path in bam_files:
        patient = os.path.basename(bam_path).split(".")[0]
        bam = pysam.AlignmentFile(bam_path, "rb")
        fasta = pysam.FastaFile(ref_fasta)

        max_depth = compute_max_depth(bam, roi_chr, roi_start, roi_end)

        del_dict = {}  # store deletion counts per position

        for read in bam.fetch(roi_chr, roi_start, roi_end):
            ref_pos = read.reference_start
            read_pos = 0
            cigartuples = read.cigartuples
            strand = "-" if read.is_reverse else "+"

            if cigartuples is None:
                continue

            for op, length in cigartuples:
                # 2 = deletion
                if op == 2:
                    for i in range(length):
                        pos = ref_pos + i + 1
                        if roi_start <= pos <= roi_end:
                            if pos not in del_dict:
                                del_dict[pos] = {"DEL_Count":0, "+":0, "-":0, "Total_Depth":0}
                            del_dict[pos]["DEL_Count"] += 1
                            del_dict[pos][strand] += 1

                # Count all reference-consuming positions for total depth
                if op in [0, 7, 8, 2, 3]:
                    for i in range(length):
                        pos = ref_pos + i + 1
                        if roi_start <= pos <= roi_end:
                            if pos not in del_dict:
                                del_dict[pos] = {"DEL_Count":0, "+":0, "-":0, "Total_Depth":0}
                            del_dict[pos]["Total_Depth"] += 1

                # Update positions
                if op in [0, 2, 3, 7, 8]:
                    ref_pos += length
                if op in [0, 1, 4, 7, 8]:
                    read_pos += length

        # Write quantitative deletion info
        for pos, counts in del_dict.items():
            if counts["DEL_Count"] > 0 and counts["Total_Depth"] >= MinDepth:
                AF_prop = counts["DEL_Count"] / counts["Total_Depth"]
                AF_percent = round(AF_prop * 100, 2)
                coverage_prop = counts["Total_Depth"] / max_depth
                Zygosity_prop = AF_prop

                # Only include deletions above AF_threshold
                if AF_prop >= AF_threshold:
                    ref_base = fasta.fetch(roi_chr, pos-1, pos)
                    writer.writerow({
                        "Patient": patient,
                        "Chromosome": roi_chr,
                        "Position": pos,
                        "Ref_Base": ref_base,
                        "Variant_Type": "DEL",
                        "Total_Depth": counts["Total_Depth"],
                        "DEL_Count": counts["DEL_Count"],
                        "+_Strand": counts["+"] ,
                        "-_Strand": counts["-"] ,
                        "Allele_Frequency(%)": AF_percent,
                        "AF_proportion": round(AF_prop, 4),
                        "Coverage_proportion": round(coverage_prop, 4),
                        "Zygosity_proportion": round(Zygosity_prop, 4)
                    })

print(f"Optimized quantitative deletions saved to: {output_file}")
