import pysam
import csv
import glob
import os

# Paths
bam_folder = "/home/wassef/Patients_Moez_Gribaa/SMA_analysis/align/"
bam_files = glob.glob(os.path.join(bam_folder, "*SMN_ROI.RG.bam"))

roi_chr = "5"
roi_start = 70910000
roi_end = 70960000

ref_fasta = "/home/wassef/Patients_Moez_Gribaa/SMA_analysis/ref/GRCh38_chr5_masked_smn2.fa"

# Allele frequency threshold
AF_threshold = 20.0  # percent

output_file = "/home/wassef/Patients_Moez_Gribaa/SMA_analysis/SMN_DEL_full_span.csv"

# Open CSV for writing
with open(output_file, "w", newline="") as csvfile:
    fieldnames = [
        "Patient", "Chromosome", "Position", "Ref_Base",
        "Variant_Type", "Total_Depth", "DEL_Count",
        "+_Strand", "-_Strand", "Allele_Frequency(%)"
    ]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for bam_path in bam_files:
        patient = os.path.basename(bam_path).split(".")[0]
        bam = pysam.AlignmentFile(bam_path, "rb")
        fasta = pysam.FastaFile(ref_fasta)

        # Create a dictionary to hold deletion counts per position
        del_dict = {}  # key = pos, value = dict with counts

        for read in bam.fetch(roi_chr, roi_start, roi_end):
            ref_pos = read.reference_start
            read_pos = 0
            cigartuples = read.cigartuples  # list of (operation,length)
            strand = "-" if read.is_reverse else "+"

            if cigartuples is None:
                continue

            for op, length in cigartuples:
                # CIGAR op 2 = deletion
                if op == 2:
                    for i in range(length):
                        pos = ref_pos + i + 1
                        if roi_start <= pos <= roi_end:
                            if pos not in del_dict:
                                del_dict[pos] = {"DEL_Count":0, "+":0, "-":0, "Total_Depth":0}
                            del_dict[pos]["DEL_Count"] += 1
                            del_dict[pos][strand] += 1
                # Count matches/mismatches for total depth
                if op in [0, 7, 8, 2, 3]:  # M, =, X, D, N
                    for i in range(length):
                        pos = ref_pos + i + 1
                        if roi_start <= pos <= roi_end:
                            if pos not in del_dict:
                                del_dict[pos] = {"DEL_Count":0, "+":0, "-":0, "Total_Depth":0}
                            del_dict[pos]["Total_Depth"] += 1
                # Update reference position
                if op in [0, 2, 3, 7, 8]:  # consumes reference
                    ref_pos += length
                # Update read position
                if op in [0, 1, 4, 7, 8]:  # consumes query/read
                    read_pos += length

        # Write positions with deletion above AF threshold
        for pos, counts in del_dict.items():
            if counts["DEL_Count"] > 0:
                AF = counts["DEL_Count"] / counts["Total_Depth"] * 100 if counts["Total_Depth"] > 0 else 0
                if AF >= AF_threshold:
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
                        "Allele_Frequency(%)": round(AF, 2)
                    })

print(f"Full-span high-value deletions saved to: {output_file}")
