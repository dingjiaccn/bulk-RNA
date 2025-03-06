#!/bin/bash

#SBATCH --job-name=bulk_rna
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=48:00:00
#SBATCH --mail-type=END
#SBATCH --account=cusanovichlab
#SBATCH --partition=standard
#SBATCH --out %x.out

# this code should ONLY be run when all 01 array jobs are done

# =========== part 1 claiming softwares input output ===========
# claiming modules from u of arizona hpc repository
module load R/4.2.2 # installed requeired: R package edgeR, ggplot2, pheatmap, tidyverse, dplyr, optparse, stringr, ggrepel

# claiming R scripts
R_figure=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_hk_bulk_rna/bulk_rna_joint_qc_figure.R
R_de=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_hk_bulk_rna/bulk_rna_de.R

# claiming samples, treatments, covariants, and output folder
# The sample names need to be consistent with ones used in 01 code
samples=(1h_1 1h_2 2h_1 2h_2 30m_1 30m_2)

# defining your controls and treatments here (only one control is allowed; multiple conditions are allowed)
# control group has to have name 'ctrl'
# the order of ${conditions} must match ${samples} order
conditions=(1_hour 1_hour 2_hour 2_hour ctrl ctrl)

# defining one more variable you would like to regress out
# the order of ${covar} must match ${samples} order
# don't define covar if not needed
#covar=(tech_1 tech_1 tech_1 tech_2 tech_2 tech_2)

# The output folder needs to be consistent with one used in array slurm job
output_dir=/xdisk/darrenc/jiachengd/20250304_hyunjin_lydia_radhika/04_hk_bulk_rna

# =========== part 2 running sum-qc DE-test ===========
echo 'running joint qc figures...'
sample_str=${samples[@]}
Rscript ${R_figure} --out "$output_dir" --samples "$sample_str"

echo 'combining count tables across samples...'
# Generate the header with Geneid and sample names
echo -e "Geneid\t$(IFS=$'\t'; echo "${samples[*]}")" > ${output_dir}/featureCounts/combined_counts.txt
# Combine the data columns using paste and process substitution
eval "paste <(tail -n +3 '${output_dir}/featureCounts/${samples[0]}' | cut -f1) \
$(printf "<(tail -n +3 '${output_dir}/featureCounts/%s' | cut -f7) " "${samples[@]}")" >> ${output_dir}/featureCounts/combined_counts.txt


echo 'running edgeR differential expression test...'
condition_str=${conditions[@]}
# Check if covar is defined
if [ -z "${covar+x}" ]; then
  covar_str=""  # If covar is not set, leave covar_str empty
else
  covar_str=${conditions[@]}
fi

# include covar if defined
# Rscript ${R_de} --out ${output_dir} --samples ${sample_str} --conditions ${condition_str} --covar ${covar_str}
# run without covar
Rscript ${R_de} --out "${output_dir}" --samples "${sample_str}" --conditions "${condition_str}"







