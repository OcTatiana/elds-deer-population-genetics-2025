# elds-deer-population-genetics-2025

This repository contains a full pipeline for variant calling, genome annotation, and population genetic analysis.

- **Quality Control:** FastQC, FastQ Screen  
- **Alignment:** BWA-MEM2  
- **Post-processing:** Picard, samtools  
- **Variant Calling:** DeepVariant  
- **Joint Genotyping:** GLnexus
- **Genome Annotation:** AUGUSTUS, GeMoMa
- **Population Structure:**: PLINK, ADMIXTURE
- **Genetic Diversity:** PLINK, vcftools
- **Demographic Histiry:** Stairway-plot-2, PCMS

---
![image](https://github.com/user-attachments/assets/d85b4e60-c3fa-4951-95e7-4eae07bf8302)

---

## 🧪 Input

- Raw paired-end FASTQ files (`sample_R1.fastq.gz`, `sample_R2.fastq.gz`,...)
- Reference genome (`reference.fa`)

---

## 📊 Step 1: Quality Control

```bash
# Run FastQC
fastqc *.fastq.gz -o qc/

# Run FastQ Screen
fastq_screen --outdir qc/ *.fastq.gz
````

---

## 🧬 Step 2: Alignment with BWA-MEM

```bash
# Indexing the reference
bwa-mem2 index -p Rucervus_eldii reference.fna

# Align and sort reads
for file in $ids; do
 	file1="${file}_1.fastq.gz";
 	file2="${file}_2.fastq.gz";
 	bwa-mem2 mem -t 8 -R "@RG\tID:$file\tPL:ILLUMINA\tPU:${file}\tLB:${file}\tSM:$file" .genome_bwa_index/Rucervus_eldii "$file1" "$file2" | samtools view -u - | samtools sort -@8 > "$file.bam";
done
```

---

## 🧹 Step 3: Mark Duplicates with Picard

```bash
for file in $ids; do
  picard MarkDuplicates I="$file.bam" O="${file}_markeddup.bam" M="${file}_metrics.txt";
done
```

---

## 🔍 Step 4: Variant Calling with DeepVariant and Joint Calling with GLnexus
  
```bash
samplesheet="samples.txt"
r1=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samplesheet | awk '{print $1}')
singularity exec -B ./deepvariant/input:/input -B ./deepvariant/output:/output deepvariant.simg \
/opt/deepvariant/bin/run_deepvariant \
--model_type=WGS \
--ref=/input/reference.fna \
--reads="/input/${r1}_markeddup.bam" \
--output_vcf="/output/${r1}.vcf.gz" \
--output_gvcf="/output/${r1}.g.vcf.gz" \
--intermediate_results_dir=/output/intermediate_results \
--num_shards=8
```
---

For multiple samples with gVCFs:

```bash
singularity run -B ./deepvariant/output:/int glnexus_v1.4.1.sif glnexus_cli \
--config DeepVariantWGS \ 
--dir /int/gl \
/int/RE113838.g.vcf.gz /int/RE114379.g.vcf.gz /int/RE114623.g.vcf.gz /int/RE115125.g.vcf.gz \
/int/RE115304.g.vcf.gz /int/RE115445.g.vcf.gz /int/RE115604.g.vcf.gz /int/RE116077.g.vcf.gz \
| bcftools view - | bgzip -c > ./deepvariant/cohort.vcf.gz
```

Filtering:

```bash
bcftools filter -i '(GQ >= 50) && (DP >= 10) && (QUAL >= 30)' cohort.vcf -o filtered_cohort_1.vcf
grep -Fv '^./.' filtered_cohort_1.vcf > filtered_cohort.vcf
```

---

## 🧬 5. Genome Annotation

### AUGUSTUS

```bash
augustus --species=human --UTR=on --outfile=predictions.gff reference.fna
```

### GeMoMa

*Cervus elaphus* (Red deer) genome was used as a reference.

```bash
java -Xmx50G -jar GeMoMa-1.9.jar CLI GeMoMaPipeline GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO restart=true o=true t=reference.fna a=./gemoma/CerEla.gff g=./gemoma/CerEla.fna outdir=output/
```

---

## 🌍 6. Population Genomics

### PCA

```bash
# Сonvert VCF to PLINK
plink --vcf filtered_cohort.vcf --make-bed --out cohort --allow-extra-chr

plink --bfile cohort --pca 8 --out pca_output --allow-extra-chr
```

The same with LD-prunning
```bash
plink --bfile cohort --extract plink.prune.in --make-bed --out prunedData --allow-extra-chr
plink --bfile prunedData --pca --allow-extra-chr --out pruned_pca
```

### ADMIXTURE

```bash
for K in {1..4}; do
  admixture --cv cohort.bed $K | tee log${K}.out
done
```
---

## 🧬 7. Genetic Diversity Analysis

### Heterozygosity

```bash
bash get_snp_heterozygosity.sh
```

### Private Alleles

```bash
bcftools view -i 'COUNT(GT="alt") = 1' filtered_cohort.vcf > private_alleles.vcf

bcftools query -l private_alleles.vcf | while read sample; do
  bcftools view -x -s "$sample" -o "${sample}.private.vcf" private_alleles.vcf
done
```

### Runs of Homozygosity (ROH)

```bash
plink --vcf filtered_cohort.vcf --allow-extra-chr --het --homozyg --homozyg-kb 1 --homozyg-snp 50 --out roh
```

### Nucleotide Diversity (π)

```bash
vcftools --vcf filtered_cohort.vcf --site-pi --out per_site_pi
```

---

## ⏳ 8. Demographic History

### Stairway Plot 2

```bash
# Prepare SFS from VCF
python easySFS/easySFS.py -i filtered_cohort.vcf -p pop.txt -o easy_sfs -a--proj 16

java -cp stairway_plot_es Stairbuilder NZCBIfolded.blueprint
bash NZCBIfolded.blueprint.sh 
```

### PCMS (Pairwise Sequentially Markovian Coalescent model)

```bash
# Create diploid fastq files
for file in $ids; do
 	bamfile="${file}_markeddup.bam";
 	samtools consensus -f FASTQ --ambig "$bamfile" | gzip > "./sfs/${file}.fq.gz"
done

# Run PCMS
for file in $ids; do
	./psmc/utils/fq2psmcfa -q20 "${file}.fq.gz" > "${file}.psmcfa";
	./psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o "${file}.psmc" "${file}.psmcfa"
done

cat *.psmc > combined.psmc
./psmc/utils/psmc_plot.pl -M 'RE113838,RE114379,RE114623,RE115125,RE115304,RE115445,RE115604,RE116077' -g 5 -u 1.5e-08 -P "left top" all ./sfs/combined.psmc 
```
---

## 📚 Tools

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
* [BWA](https://github.com/bwa-mem2/bwa-mem2)
* [Picard](https://broadinstitute.github.io/picard/)
* [DeepVariant](https://github.com/google/deepvariant)
* [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
* [Augustus](https://github.com/Gaius-Augustus/Augustus)
* [GeMoMa](https://www.jstacs.de/index.php/GeMoMa)
* [PLINK](https://www.cog-genomics.org/plink/)
* [Stairway Plot 2](https://github.com/xiaoming-liu/stairway-plot-v2)
* [EasySFS](https://github.com/isaacovercast/easySFS)
* [PCMS](https://github.com/lh3/pcms)
