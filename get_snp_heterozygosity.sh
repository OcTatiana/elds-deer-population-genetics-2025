#!/bin/bash

#vcf="$1"
vcf="./vcf/filtered/gatk_snp.vcf.gz"         # Input VCF file (indexed)
chrom="$1"                # Chromosome to analyze
chr_length="$2"        # Chromosome length (adjust accordingly)
win_size=1000000            # Window size: 1 Mbp
step_size=100000            # Step size: 100 kbp

# Extract positions and genotypes for all samples
# Format: CHROM POS SAMPLE1_GT SAMPLE2_GT ... SAMPLE8_GT
bcftools query -f '%CHROM\t%POS[\t%GT]\n' $vcf -r $chrom > snps_genotypes.txt

# Get number of samples (columns - 2)
num_samples=$(head -1 snps_genotypes.txt | awk '{print NF-2}')

# Sliding windows loop
for ((start=1; start<=chr_length-win_size; start+=step_size)); do
  end=$((start + win_size - 1))

  # Initialize counts array for samples
  declare -a het_counts
  for ((i=0; i<num_samples; i++)); do
    het_counts[i]=0
  done
  

  # Process SNPs in window
	while IFS=$'\t' read -r chr pos rest; do
	  genotypes=($rest)
	  for ((i=0; i<num_samples; i++)); do
		gt=${genotypes[i]}
		if [[ "$gt" =~ ^(0[/|]1|1[/|]0)$ ]]; then
		  (( het_counts[i]++ ))
		fi
	  done
	done < <(awk -v s=$start -v e=$end '($2 >= s && $2 <= e)' snps_genotypes.txt)


	 # Output window and scaled counts (SNPs per kbp)
  # Format: CHROM START END SAMPLE1_SNPs_per_kbp ... SAMPLE8_SNPs_per_kbp
  printf "%s\t%d\t%d" "$chrom" "$start" "$end" >> res_het.txt
  for ((i=0; i<num_samples; i++)); do
    # Scale count by 1000 (window size in kbp)
    scaled=$(echo "scale=4; ${het_counts[i]} / 1000" | bc)
    printf "\t%s" "$scaled" >> res_het.txt
  done
  echo >> res_het.txt
done
