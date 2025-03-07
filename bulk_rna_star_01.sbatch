#!/bin/bash

#SBATCH --job-name=bulk_rna
#SBATCH --nodes=1
#SBATCH --ntasks=60
#SBATCH --ntasks-per-node=60
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --out %x_%a.out
#SBATCH --array=1-6

# adjust --ntasks: number of cores per sample × number of samples. Increasing cores per sample also increasing memo per sample
# keep --ntasks-per-node: Puma (≤ 94 cores), Ocelote (≤ 28 cores), Elgato (≤ 16 cores)  
# Adjust --array=1-n: where n is the total number of samples to run in parallel 


# =========== part 1 claiming softwares input output ===========
# claiming modules from u of arizona hpc repository
# install required: sampatools, R, if using locally
module load samtools/1.19.2
module load R/4.2.2 # installed requeired: R package edgeR, ggplot2, pheatmap, tidyverse, dplyr, optparse, stringr, ggrepel

samples=(1h_1 1h_2 2h_1 2h_2 30m_1 30m_2)
s_indexes=(S3 S4 S1 S2 S5 S6)
index=$((SLURM_ARRAY_TASK_ID-1))
sample=${samples[$index]}
s_index=${s_indexes[$index]}

# don't run. test only
#sample=1h_1
#s_index=S3

# claiming software
trimmomatic=/home/u15/jiachengd/Trimmomatic-0.36
fastqc=/home/u15/jiachengd/anaconda3/envs/python39/bin/fastqc # intall required: fastqc
star=/groups/darrenc/Jiacheng/cellranger-8.0.1/lib/bin/STAR # borrow star from 10x cellranger, avialable for all group members
python=/home/u15/jiachengd/anaconda3/envs/python39/bin/python3.9 # install required: python >3.7
rseqc_dir=/home/u15/jiachengd/anaconda3/envs/python39/bin # install required: rseqc (python >3.7 env prerequisite:numpy scipy matplotlib pysam)
picard=/home/u15/jiachengd/picard.jar # install required: picard
java=/home/u15/jiachengd/jdk-21.0.1/bin/java #install required: jdk
featurecounts=/home/u15/jiachengd/subread-2.0.8-Linux-x86_64/bin/featureCounts # install required: subread

# claiming genome ref
human_hg38_transcript_bed=/groups/darrenc/references/beds/hg38_Gencode_V28.bed
mouse_mm10_transcript_bed=./not_available.bed

human_hg38_fa=/groups/darrenc/references/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa # borrow genome ref from 10x
mouse_mm10_fa=/groups/darrenc/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa # borrow genome ref from 10x

human_hg38_gtf=/groups/darrenc/Jiacheng/10x_hg38_genes.gtf # borrow genome ref from 10x (36601 features)
mouse_mm10_gtf=./not_available.gtf
#human_hg38_gtf=/groups/darrenc/references/gtf/hg38/gencode.v36.annotation.gtf # alternatively, use genecode gtf (59427 features)

human_hg38_star=/groups/darrenc/references/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/star # borrow genome ref from 10x
mouse_mm10_star=/groups/darrenc/references/10x/refdata-cellranger-arc-mm10-2020-A-2.0.0/star # borrow genome ref from 10x

# claiming fastq input
fastq_dir=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/02_fastqs/hk_rna
fastq1=${fastq_dir}/${sample}_${s_index}_R1_001.fastq.gz
fastq2=${fastq_dir}/${sample}_${s_index}_R2_001.fastq.gz

# claiming output dirs / files
#output_dir=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_test
output_dir=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_hk_bulk_rna
bam=${output_dir}/star/${sample}/Aligned.sortedByCoord.out.bam
dedup_bam=${output_dir}/star/${sample}/${sample}_dedup.bam

# create output dirs
mkdir -p ${output_dir}/reports/trimmomatic
mkdir -p ${output_dir}/fastqs_trimmed
mkdir -p ${output_dir}/reports/fastqc
mkdir -p ${output_dir}/star
mkdir -p ${output_dir}/star/${sample}
mkdir -p ${output_dir}/reports/genebody_coverage
mkdir -p ${output_dir}/reports/reads_distribution
mkdir -p ${output_dir}/reports/insertion_size
mkdir -p ${output_dir}/reports/dedup
mkdir -p ${output_dir}/featureCounts

date

# =========== part 2 running trimmomatic fastqc star featurecounts rseqc ===========
echo 'running trimmomatic read adpator (R1 R2) trimming...'
# choose nextera / trueseq read adptors according to your lib
# adaptor clip parameter: 2-seed mismatches, 30-palindrome clip threshold, 10-simple clip threshold, 10-minimum adapter length after clipping, true-keep both reads if only one is trimmed
${java} -jar ${trimmomatic}/trimmomatic-0.36.jar \
PE ${fastq1} ${fastq2} \
${output_dir}/fastqs_trimmed/${sample}_R1.fastq.paired.trimmed.gz ${output_dir}/fastqs_trimmed/${sample}_R1.fastq.unpaired.trimmed.gz \
${output_dir}/fastqs_trimmed/${sample}_R2.fastq.paired.trimmed.gz ${output_dir}/fastqs_trimmed/${sample}_R2.fastq.unpaired.trimmed.gz \
ILLUMINACLIP:${trimmomatic}/adapters/NexteraPE-PE.fa:2:30:10:10:true \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 2> ${output_dir}/reports/trimmomatic/${sample}_trimmomatic.log;

echo 'running fastqc quality check...'
${fastqc} \
    ${output_dir}/fastqs_trimmed/${sample}_R1.fastq.paired.trimmed.gz \
    ${output_dir}/fastqs_trimmed/${sample}_R2.fastq.paired.trimmed.gz \
    -o ${output_dir}/reports/fastqc

echo 'running star mapping...'
# --outFilterMultimapNmax N max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped (default:20, set to 1)
# --quantTranscriptomeBan Singleend prohibits single-end alignment
${star} \
     --runThreadN 10 \
     --genomeDir ${human_hg38_star} \
     --readFilesIn ${output_dir}/fastqs_trimmed/${sample}_R1.fastq.paired.trimmed.gz ${output_dir}/fastqs_trimmed/${sample}_R2.fastq.paired.trimmed.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ${output_dir}/star/${sample}/ \
     --outSAMtype BAM SortedByCoordinate \
     --quantTranscriptomeBan Singleend \
     --quantMode GeneCounts \
     --outSAMattributes NH HI AS NM MD \
     --outFilterMultimapNmax 1
samtools index ${bam}

# echo 'running pciard deduplication...'
${java} -Xmx1G -jar ${picard} MarkDuplicates REMOVE_DUPLICATES=true I=${bam} O=${dedup_bam} M=${output_dir}/reports/dedup/${sample}.dedup.log
samtools index ${dedup_bam}

echo 'running gene-wise quantification...'
# remove -p --countReadPairs if working with single-end reads
# adjust -t options are to exon transcript UTR start_codon stop_codon gene
# adjust -g gene_id if wish to mapped to gene_id
${featurecounts} -T 10 -p --countReadPairs -B \
    -t gene -g gene_name -a ${human_hg38_gtf} \
    -o ${output_dir}/featureCounts/${sample} ${dedup_bam}

echo 'running genebody coverge...'
${python} ${rseqc_dir}/geneBody_coverage.py \
     --input ${dedup_bam} \
     --refgene ${human_hg38_transcript_bed} \
     --minimum_length 100 \
     --out-prefix ${output_dir}/reports/genebody_coverage/${sample}
# move genebody_coverage log file to output folder
mv log.txt ${output_dir}/reports/genebody_coverage/${sample}.txt
     
echo 'running read distribution...'
${python} ${rseqc_dir}/read_distribution.py \
    --input ${dedup_bam} \
    --refgene ${human_hg38_transcript_bed} | sed '1,4d;$d' > ${output_dir}/reports/reads_distribution/${sample}.txt

echo 'running picard insertion size...'
${java} -Xmx1G -jar ${picard} \
    CollectInsertSizeMetrics \
    INPUT=${dedup_bam} \
    OUTPUT=${output_dir}/reports/insertion_size/${sample}.txt \
    H=${output_dir}/reports/insertion_size/${sample}.pdf

echo 'counting reads loss...'
fastq1_counts=$(($(zcat ${fastq1} | wc -l) / 4))
fastq2_counts=$(($(zcat ${fastq2} | wc -l) / 4))
trimmed_fastq_counts=$(($(zcat ${output_dir}/fastqs_trimmed/${sample}_R1.fastq.paired.trimmed.gz | wc -l) / 4))
bam_counts=$(($(samtools view -c ${bam}) / 2 ))
dedup_bam_counts=$(($(samtools view -c ${dedup_bam}) / 2 ))
assigned_counts=$(awk '$1 == "Assigned" {print $2}' ${output_dir}/featureCounts/${sample}.summary)

echo "Fastq1 input read counts: ${fastq1_counts}" > ${output_dir}/reports/${sample}_read_counts_summary.txt
echo "Fastq2 input read counts: ${fastq2_counts}" >> ${output_dir}/reports/${sample}_read_counts_summary.txt
echo "Read pairs passed trimmomatic: ${trimmed_fastq_counts}" >> ${output_dir}/reports/${sample}_read_counts_summary.txt
echo "Uniquely mapped read pairs: ${bam_counts}" >> ${output_dir}/reports/${sample}_read_counts_summary.txt
echo "Uniquely mapped unique read pairs: ${dedup_bam_counts}" >> ${output_dir}/reports/${sample}_read_counts_summary.txt
echo "Uniquely annotated read pairs: ${assigned_counts}" >> ${output_dir}/reports/${sample}_read_counts_summary.txt

date


