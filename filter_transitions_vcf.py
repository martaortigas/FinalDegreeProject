import os
import pysam
from cyvcf2 import VCF

#######EXTRACT HAPLOTYPES STEP 1 (EXCLUDING TRANSITIONS)#######
# This script iterates through a VCF and a coordinate file (txt file)
# and for a variant that is not a transition, it retrieves the chromosome, position and genotype
# of each sample if the variant is inside the region of the coordinate file.

# Specify and set the path to the working directory
directory_path = "/scratch/125-emmer/MartaOrtigas/pysam"
os.chdir(directory_path)

input_file = "SHORTER_Dxy_WSL_introgr_to_DSE.txt"
input_vcf = "merged_VCFs.vcf.gz"

vcf = VCF(input_vcf)

# Open the coordinates file
with open(input_file, "r") as inp:
    # Iterate over the coordinates file
    for line in inp:
        columns = line.strip().split("\t")
        chrom, start, end = columns
        start = int(start)
        end = int(end)

        # Create output file name
        output_file = f"NOtrans_{chrom}_{start}-{end}_region.txt"

        # Open the output file
        with open(output_file, "w") as out:
            # Iterate over the VCF file
            for v in vcf(chrom + ":" + str(start) + "-" + str(end)):
                if v.is_snp == True:
                    if v.is_transition == False:
                        row = []
                        chromos = v.CHROM
                        pos = v.POS
                        row.append(chromos)
                        row.append(str(pos))
                        # Iterate over the sample columns (genotypes) of the VCF file
                        for sample in range(0, len(v.genotypes)):
                            if v.genotypes[sample] == [0, 0, True]:
                                ref = v.REF
                                row.append(ref)
                            elif v.genotypes[sample] == [1, 1, True]:
                                alt = v.ALT[0]
                                row.append(alt)
                            elif v.genotypes[sample] == [-1, -1, False]:
                                row.append("N")
                        output_string = "\t".join(row)
                        out.write(f"{output_string}\n")
