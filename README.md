# Bulk RNA-Seq Analysis Pipeline

A two-stage SLURM-based pipeline for processing bulk RNA-seq data, including quality control, alignment, quantification, and differential expression analysis.

## Pipeline Overview

### Stage 1: Per-Sample Processing (`bulk_rna_star_01.sbatch`)
![Stage 1 Workflow](https://via.placeholder.com/800x400.png?text=QC+Alignment+Quantification)

1. **Adapter Trimming**: Trimmomatic 
2. **Quality Control**: FastQC
3. **Alignment**: STAR (hg38 / mm10 reference)
4. **Quantification**: featureCounts 
5. **QC Metrics**:
   - RSeQC (gene body coverage, read distribution)
   - Picard (insert size metrics)
   - Read tracking through pipeline

### Stage 2: Joint Analysis (`bulk_rna_star_02.sbatch`)
![Stage 2 Workflow](https://via.placeholder.com/800x400.png?text=Combine+Counts+DE+Analysis)

1. **Combine Count Tables**
2. **Generate joint QC Reports**
3. **Differential Expression**: edgeR GLM pipeline
4. **Visualization** (heatmap PCA)

## Prerequisites

### Software Dependencies
- **Core Tools**:
  - Trimmomatic v0.36
  - STAR 2.7.2a
  - featureCounts (subread v2.0.8)
  - RSeQC v4.0.0
  - Picard v2.27.x
  - R v4.2.2 (edgeR, ggplot2, pheatmap, tidyverse, dplyr, optparse, stringr, ggrepel)
  - Python > 3.7 (numpy scipy matplotlib pysam)
  - JDK 21.0.1

- **Reference Files**:
  - hg38 / mm10 STAR index
  - GENCODE annotation (GTF/BED12)
  - Nextera adapter sequences

## Job Execution

### 1. Per-Sample Processing (Array Job)
```bash
# Submit parallel array jobs
sbatch bulk_rna_star_01.sbatch

# Submit joint analysis job
sbatch bulk_rna_star_02.sbatch # submit when 01 is done

# Monitor job progress
squeue -u $USER

```
### Output Structure

featureCounts/
└── <sample> # Gene counts
└── combined_counts.txt # Merged gene counts          
star/
└── <sample>/
    ├── Aligned.sortedByCoord.out.bam # Sorted BAM
    └── Log.final.out         # STAR logs
reports/
├── fastqc/           # Quality reports
├── trimmomatic/      # Trimming logs
└── genebody_coverage/ # Coverage plots
└── insertion_size/ # Insertion size
└── reads_distribution/ # Reads distribution on genome
└── <sample>_read_counts_summary.txt # Read tracking
DE_results/
└── heatmap.pdf # Pearson cor heatmap
└── pca.pdf # PCA plot
└── de_test_result.csv # edgeR GLM DE test                
sum_reports/
└── genebody_coverage_sum.csv # genebody coverage data
└── genebody_coverage_sum.pdf # genebody coverage plot
└── insertion_size_sum.csv # insertion size data
└── insertion_size_sum.pdf # insertion size plot
└── read_counts_sum.csv # read tracking data
└── read_counts_sum.pdf # read tracking plot
└── read_distribution_sum.csv # read distribution data
└── read_distribution_sum.pdf # read distribution plot

## Customization guide
### 1. Array Job

```bash
samples=(1h_1 1h_2 2h_1 2h_2 30m_1 30m_2) # enter your sample names
s_indexes=(S3 S4 S1 S2 S5 S6) # according fastq S index (order needs to match ${samples})

# modify to your software path
trimmomatic=/home/u15/jiachengd/Trimmomatic-0.36
fastqc=/home/u15/jiachengd/anaconda3/envs/python39/bin/fastqc 
python=/home/u15/jiachengd/anaconda3/envs/python39/bin/python3.9
rseqc_dir=/home/u15/jiachengd/anaconda3/envs/python39/binrseqc (python >3.7 env prerequisite:numpy scipy matplotlib pysam)
picard=/home/u15/jiachengd/picard.jar
java=/home/u15/jiachengd/jdk-21.0.1/bin/java
featurecounts=/home/u15/jiachengd/subread-2.0.8-Linux-x86_64/bin/featureCounts
```

### 2. Joint Job
```bash
# modify code path
R_figure=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_hk_bulk_rna/bulk_rna_joint_qc_figure.R
R_de=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_hk_bulk_rna/bulk_rna_de.R

# The sample names need to be consistent with ones used in 01 code
samples=(1h_1 1h_2 2h_1 2h_2 30m_1 30m_2)
# defining your controls and treatments here (only one control is allowed; multiple conditions are allowed)
# control group has to have name 'ctrl'
# the order of ${conditions} must match ${samples} order
conditions=(1_hour 1_hour 2_hour 2_hour ctrl ctrl)
# defining one variable you would like to regress out
# the order of ${covar} must match ${samples} order
# don't define covar if not needed
#covar=(tech_1 tech_1 tech_1 tech_1 tech_1 tech_1)
```
