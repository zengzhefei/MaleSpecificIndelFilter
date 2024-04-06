from cyvcf2 import VCF, Writer

# Brief Introduction:
# This script is designed to filter indel variants from a VCF file, specifically targeting those variants 
# that are exclusively present in male samples and meet certain length criteria (insertions > 10bp, deletions > 20bp),
# making them potential markers for male specificity. This script is intended for use in genomic studies 
# where identifying sex-specific markers is required. 
# Before running the script, please replace 'input_vcf_path' with the path to your input VCF file 
# and 'output_vcf_path' with the desired path for the output file containing filtered variants.

# Define input and output file paths
input_vcf_path = 'ht.chrall.filtered.indel.recode.vcf'  # Replace with the path to your input VCF file
output_vcf_path = 'male_specific_markers_final.vcf'  # Replace with your desired output file path

# Define lists of male and female sample IDs
male_samples_ids = [
    "male01", "male02", "male03", "male04"
]

female_samples_ids = [
    "female01", "female02", "female03", "female04"
]

# Create VCF reader and writer
vcf_reader = VCF(input_vcf_path)
vcf_writer = Writer(output_vcf_path, vcf_reader)

# Ensure sample IDs are present in the VCF file
actual_male_samples = [sample for sample in male_samples_ids if sample in vcf_reader.samples]
actual_female_samples = [sample for sample in female_samples_ids if sample in vcf_reader.samples]

# Iterate through each variant in the VCF file
for variant in vcf_reader:
    all_male_variant = True  # Flag indicating whether all male samples have the variant
    no_female_variant = True  # Flag indicating whether no female samples have the variant

    # Calculate the difference in length between REF and ALT to determine if the variant meets the length criteria
    ref_len = len(variant.REF)
    alt_len = len(str(variant.ALT[0]))  # Assuming there is only one ALT allele, additional logic may be needed for multiple ALTs
    insert_size = abs(alt_len - ref_len)

    # Skip the variant if it does not meet the length criteria for insertions or deletions
    if not ((alt_len > ref_len and insert_size > 10) or (alt_len < ref_len and insert_size > 20)):
        continue

    # Check male samples
    for sample_id in actual_male_samples:
        sample = variant.genotypes[vcf_reader.samples.index(sample_id)]
        # Mark as False if any male sample does not have the variant
        if not (sample[0] > 0 or sample[1] > 0):
            all_male_variant = False
            break

    # Check female samples
    if all_male_variant:
        for sample_id in actual_female_samples:
            sample = variant.genotypes[vcf_reader.samples.index(sample_id)]
            # Mark as False if any female sample has the variant
            if sample[0] > 0 or sample[1] > 0:
                no_female_variant = False
                break

    # Write the variant to the output file if all male samples have the variant and no female samples have it
    if all_male_variant and no_female_variant:
        vcf_writer.write_record(variant)

# Close the file
vcf_writer.close()
