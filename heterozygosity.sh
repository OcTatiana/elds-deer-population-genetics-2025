while IFS=$'\t' read -r chr len rest; do
  if [ "$rest" != "0" ]; then
    bash get_snp_heterozygosity.sh "$chr" "$len"
  fi
done < ./vcf/filtered/snp_len_scaffold_gatk.txt
