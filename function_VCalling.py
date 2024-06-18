import os
import pysam

#######PERFORM VARIANT CALLING TO A SINGLE BAM#######
# This script takes as input a BAM file and a VCF file and performs
# variant calling for the BAM file. The resulting variants are 
# displayed in a text file.

# Specify and set the path to the working directory
directory_path = "/scratch/125-emmer/MartaOrtigas/pysam"
os.chdir(directory_path)

def variant_calling(bam_file, vcf_file, output_file):
    # Open BAM and VCF files
    bam = pysam.AlignmentFile(bam_file, "rb", index_filename=bam_file + ".csi")
    vcf = pysam.VariantFile(vcf_file)

    # Prepare output file
    with open(output_file, "w") as out:
        # Iterate over variants in VCF file
        for record in vcf.fetch():
            chrom = record.chrom
            pos = record.pos
            ref = record.ref
            alt = record.alts[0]

            # Initialize allele counts
            ref_count = 0
            alt_count = 0

            # Initialize sequence list
            sequences = []

            # Coverage counters before and after applying the filters to see the variation in the values
            coverage_prev = 0 # Counter of coverage before filters
            coverage_pos = 0 # Counter of coverage after filters

            # Fetch reads overlapping the variant position
            for read in bam.fetch(chrom, pos - 1, pos):
                coverage_prev += 1

                # Filtering criteria
                if read.is_duplicate == False:
                    if read.mapping_quality > 20:
                        base_qualities = read.query_qualities
                        # Check if the position is within the valid range of base_qualities
                        if 0 <= pos - read.reference_start - 1 < len(base_qualities):
                            if base_qualities[pos - read.reference_start - 1] > 20:
                                if pos + 2 < read.reference_end:
                                    if pos - 2 > read.reference_start:
                                        coverage_pos += 1

                                        # Check if read supports reference allele
                                        if read.query_sequence[pos - read.reference_start - 1] == ref:
                                            ref_count += 1
                                            seq = read.query_sequence
                                            sequences.append(seq + " " + "(R)")  # Append "R" for reference allele

                                        # Check if read supports alternate allele
                                        elif read.query_sequence[pos - read.reference_start - 1] == alt:
                                            alt_count += 1
                                            seq = read.query_sequence
                                            sequences.append(seq + " " + "(A)")  # Append "A" for alternate allele
                        else:
                            # If position is out of range, skip processing this read
                            continue

            # Write annotation to output file
            if ref_count > 0 or alt_count > 0:
                # Join sequences into a single string separated by whitespace
                sequence_str = " ".join(sequences) if len(sequences) > 0 else ""
                out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{ref_count},{alt_count}\t{coverage_prev}\t{coverage_pos}\t{sequence_str}\n")

    # Close files
    bam.close()
    vcf.close()

# Example usage
bam_file = "/scratch/125-emmer/exome_capture/NW21_exome_sorted_dupMarked_RG.bam"
vcf_file = "MODIFIED_all_emmer.biallSNP.hardfiltPASS.nooutl_NATChr.NOhetpositions.rename.phased.Svevomask.Nomiss.vcf.gz"
output_file = "variant_calls_NW21.txt"

variant_calling(bam_file, vcf_file, output_file)

