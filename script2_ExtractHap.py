import os
import csv
import glob

#######EXTRACT HAPLOTYPES STEP 2#######
# This script takes as input the text files generated by script1 and retrieves
# the haplotypes and positions them horizontally. In the input files each column
# corresponds to an haplotype (except column 1 and 2)

# Specify and set the path to the working directory
directory_path = "/scratch/125-emmer/MartaOrtigas/pysam/SNPs/SNPs_ancient/haplotypes_NOtrans_txts"
os.chdir(directory_path)

for file in glob.glob("*region.txt*"):
	input_file = file
	output_file = f"Hap_{file}"

	# Open the input file
	with open(input_file, "r") as inp:
    		# Create a CSV reader
    		reader = csv.reader(inp, delimiter='\t')

    		# Initialize an empty list to store column values
    		column_values = []

    		# Read the rows
    		for row in reader:
        		for i in range(2, len(row)):  # Start from index 2 to iterate over columns
            			if len(column_values) <= i - 2:
                			column_values.append(row[i])
            			else:
                			column_values[i - 2] += row[i]

	# Write the concatenated strings of column values to a new output file
	with open(output_file, "w") as out:
    		for value in column_values:
        		out.write(value + '\n')  # Write each column value followed by a newline
