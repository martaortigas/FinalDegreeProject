import os
import glob

#######MODIFY HEADERS FASTA#######
# This script takes the FASTAs generates by "script3_ExtractHap.py" and
# changes the header of each haplotype using the information of a input text file
# containing the haplotypeID (number of the individual), the name of the country where it belongs, and the ID.

# Specify and set the path to the working directory
directory_path = "/scratch/125-emmer/MartaOrtigas/pysam/SNPs/SNPs_ancient/FastaFiles_ancient/fastas_NOtrans_ancient"
os.chdir(directory_path)

input_file = "haplotype_country.txt"

# Read mapping of haplotypes to countries
haplotype_country = {}
with open(input_file, 'r') as inp:
	for line in inp:
		if line.strip() and not line.startswith("Haplotype_ID"):
			haplotype_id, country, id = line.strip().split()
			haplotype_country[haplotype_id] = (country, id)


for file in glob.glob("*region.fasta*"):
	input_fasta = file
	output_file = f"NewHead_{file}"
	# Read and modify the FASTA file
	with open(input_fasta, 'r') as f:
		with open(output_file, 'w') as out:
			for line in f:
				if line.startswith('>'):
					haplotype_id = line[1:].strip()
					country, id = haplotype_country.get(haplotype_id, ("Unknown", "Unknown"))
					new_header = f">{haplotype_id}_{country}_{id}\n"
					out.write(new_header)
				else:
					out.write(line)

