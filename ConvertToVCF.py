import os
import pysam

#######CONVERT TO VCF#######
# This script transforms text files into a VCF format with its required structure

# Specify and set the path to the working directory
directory_path = "/scratch/125-emmer/MartaOrtigas/pysam"
os.chdir(directory_path)

def convert_to_vcf(input_file, output_file):
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as vcf:
            vcf.write("##fileformat=VCFv4.3\n")
            vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tIP10\n")
            for line in f:
                fields = line.split()
                chrom = fields[0]
                pos = fields[1]
                ref = fields[2]
                alt = fields[3]
                genotype = fields[4]
                
                geno = genotype.split(",")

                # Convert genotype values to integers
                geno_0 = int(geno[0])
                geno_1 = int(geno[1])

                # Calculate the ratio to determine if one value is at least 70% of the other
                if geno_0 > geno_1:
                    if geno_1 >= 0.3 * geno_0:
                        vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t.|.\n")
                    else:
                        vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t0|0\n")
                elif geno_0 < geno_1:
                    if geno_0 >= 0.3 * geno_1:
                        vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t.|.\n")
                    else:
                        vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t1|1\n")
                else:
                    vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t.|.\n")

# Usage
input_file = "variant_calls_IP10_merged.txt"
output_file = "/scratch/125-emmer/MartaOrtigas/pysam/VCFs/variant_calls_IP10_merged.vcf"
convert_to_vcf(input_file, output_file)

