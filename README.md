# elds-deer-population-genetics-2025

Some intro

- **Quality Control:** FastQC, FastQ Screen  
- **Alignment:** BWA-MEM  
- **Post-processing:** Picard  
- **Variant Calling:** DeepVariant  
- **Joint Genotyping (Optional for multi-sample):** GLnexus  

---

## 🔧 Requirements

- `fastqc`
- `fastq_screen`
- `bwa-mem2`
- `samtools`
- `picard`
- `deepvariant`
- `glnexus`
- `python3`, `bash`, and standard UNIX utilities

We recommend using [conda](https://docs.conda.io/en/latest/) or [Docker](https://www.docker.com/) for environment management.

---

## 🧪 Input

- Raw paired-end FASTQ files (`sample_R1.fastq.gz`, `sample_R2.fastq.gz`)
- Reference genome (`reference.fa`), indexed for BWA and DeepVariant

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

## 🔍 Step 4: Variant Calling with DeepVariant
  
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

## 🧬 Step 5: Joint Calling with GLnexus

For multiple samples with gVCFs:

```bash
# Create a sample list file
echo -e "sample1.g.vcf.gz\nsample2.g.vcf.gz" > gvcf_list.txt

# Run GLnexus
glnexus_cli --config DeepVariant $(cat gvcf_list.txt) > cohort.vcf
```

---

## 📁 Output

* `sample.vcf.gz` — Variant calls (VCF)
* `sample.g.vcf.gz` — gVCF for joint calling
* `cohort.vcf` — Joint genotyped VCF (if multiple samples)

---

## 📚 References

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
* [BWA](https://github.com/bwa-mem2/bwa-mem2)
* [Picard](https://broadinstitute.github.io/picard/)
* [DeepVariant](https://github.com/google/deepvariant)
* [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
