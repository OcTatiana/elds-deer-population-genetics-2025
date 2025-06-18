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
- `bwa`
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
  

# Convert to BAM, sort and index

samtools index sample.sorted.bam
```

---

## 🧹 Step 3: Mark Duplicates with Picard

```bash
picard MarkDuplicates \
    I=sample.sorted.bam \
    O=sample.dedup.bam \
    M=sample.metrics.txt

samtools index sample.dedup.bam
```

---

## 🔍 Step 4: Variant Calling with DeepVariant

```bash
# Run DeepVariant (example with Docker)
docker run -v "${PWD}":"/input" google/deepvariant:1.5.0 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/input/reference.fa \
  --reads=/input/sample.dedup.bam \
  --output_vcf=/input/sample.vcf.gz \
  --output_gvcf=/input/sample.g.vcf.gz \
  --num_shards=8
```

---

## 🧬 Step 5 (Optional): Joint Genotyping with GLnexus

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

## 📜 Notes

* All output paths can be customized.
* Use appropriate `--model_type` in DeepVariant (`WGS`, `WES`, etc.)
* GLnexus config can be customized for better cohort-level calling performance.

---

## 📚 References

* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)
* [BWA](http://bio-bwa.sourceforge.net/)
* [Picard](https://broadinstitute.github.io/picard/)
* [DeepVariant](https://github.com/google/deepvariant)
* [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
