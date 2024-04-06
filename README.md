# Brief Introduction:
# This script is designed to filter indel variants from a VCF file, specifically targeting those variants 
# that are exclusively present in male samples and meet certain length criteria (insertions > 10bp, deletions > 20bp),
# making them potential markers for male specificity. This script is intended for use in genomic studies 
# where identifying sex-specific markers is required. 
# Before running the script, please replace 'input_vcf_path' with the path to your input VCF file 
# and 'output_vcf_path' with the desired path for the output file containing filtered variants.