#!/usr/bin/env Rscript
# bulk_rna_joint_qc_figure.R

# loading packages
suppressPackageStartupMessages({
  library(optparse)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(edgeR)
})

# ==== defining arguments =====
option_list <- list(
  make_option(c("--out"), 
              type="character", 
              default=".", 
              help="Output directory [default = %default]", 
              metavar="DIR"),
  make_option(c("--samples"), 
              type="character", 
              default=NULL, 
              help="Space-separated list of sample names", 
              metavar="SAMPLES")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

samples <- strsplit(opt$samples, " ")[[1]]
samples <- sort(samples)
output_folder <- opt$out

# ==== finding existing individual qc report ====
reads_sum_files <- list.files(
  paste0(output_folder,'/reports'), 
  pattern = '.txt', 
  full.names = T
  ) %>% sort()

temp_list <- vector('list', length(samples))
for (i in 1:length(samples)){
  temp_list[[i]] <- read.table(reads_sum_files[i],
                               sep = ':')
  colnames(temp_list[[i]]) <- c('Read Pairs', 'Counts')
  temp_list[[i]]$Sample <- samples[i]
}
reads_sum_df <- do.call(rbind, temp_list) %>%
  dplyr::mutate(`Read Pairs` = factor(`Read Pairs`,
                                      c("Fastq1 input read counts",
                                        "Fastq2 input read counts",
                                        "Read pairs passed trimmomatic",
                                        "Uniquely mapped read pairs",
                                        "Uniquely annotated read pairs")
                                      )
                )

p1 <- ggplot(reads_sum_df, aes(x = Sample, y = Counts, fill = `Read Pairs`)) +
  geom_bar(stat = "identity", position = "dodge", color='white') +
  geom_text(aes(label = Counts), 
            position = position_dodge(0.9), 
            vjust = 0.5, hjust = 5 ,size = 3,  angle = 90)  +
  theme_classic() +
  labs(
    x = "Sample", 
    y = "Counts", 
    color = "Read Pairs"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 15),  
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 15), 
    axis.title = element_text(size = 15),  
    axis.ticks = element_line(linewidth = 1), 
    axis.line = element_line(linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 2)
  )

reads_distribution_files <- list.files(
  paste0(output_folder,'/reports/reads_distribution'), 
  pattern = '.txt', 
  full.names = T
) %>% sort()
temp_list <- vector('list', length(samples))
for (i in 1:length(samples)){
  temp_list[[i]] <- read.table(reads_distribution_files[1], header = T) %>%
    dplyr::mutate(Fraction = round(Tag_count / sum(Tag_count) * 100,1))
  temp_list[[i]]$Sample <- samples[i]
}
reads_dist_sum_df <- do.call(rbind, temp_list)

p2 <- ggplot(reads_dist_sum_df, aes(x = "", y = Fraction, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ Sample) +
  theme_classic() +
  labs(
    x = "", 
    y = "Read Distribution", 
    color = "Group"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text = element_text(size = 15),  
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 15), 
    axis.title = element_text(size = 15),  
    axis.ticks = element_line(linewidth = 1), 
    axis.line = element_line(linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 2)
  )

genebody_files <- list.files(
  paste0(output_folder,'/reports/genebody_coverage'), 
  pattern = 'Coverage.txt', 
  full.names = T
) %>% sort()
temp_list <- vector('list', length(samples))
for (i in 1:length(samples)){
  temp_list[[i]] <- read.table(genebody_files[i], row.names = 1) %>% t() %>%
    as.data.frame()
  rownames(temp_list[[i]]) <- NULL
  colnames(temp_list[[i]]) <- c('Percentile of Gene', 'Tags')
  temp_list[[i]] <- temp_list[[i]] %>%
    dplyr::mutate(
      Coverage = (Tags - min(Tags)) / (max(Tags) - min(Tags))
      )
  temp_list[[i]]$Sample <- samples[i]
}
genebody_sum_df <- do.call(rbind, temp_list)
p3 <- ggplot(genebody_sum_df, aes(x = `Percentile of Gene`, y = Coverage, color = Sample)) +
  geom_line(size = 1, alpha = 0.5) +   # Connecting lines
  labs(
    x = "Percentile of Gene", 
    y = "Coverage", 
    color = "Sample"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 15),  
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 15), 
    axis.title = element_text(size = 15),  
    axis.ticks = element_line(linewidth = 1), 
    axis.line = element_line(linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 2)
  )

insertion_files <- list.files(
  paste0(output_folder,'/reports/insertion_size'), 
  pattern = '.txt', 
  full.names = T 
) %>% sort()

temp_list <- vector('list', length(samples))
for (i in 1:length(samples)){
  temp_list[[i]] <- read.table(insertion_files[i], 
                               header = T,
                               skip = 10)
  colnames(temp_list[[i]]) <- c('Insert Size', 'Frequency')
  
  # expand frequency
  temp_list[[i]] <- temp_list[[i]][rep(1:nrow(temp_list[[i]]), temp_list[[i]]$Frequency), ]
  
  temp_list[[i]]$Sample <- samples[i]
  
}
insertion_sum_df <- do.call(rbind, temp_list)
p4 <- ggplot(insertion_sum_df, aes(x = `Insert Size`, fill = Sample, color = Sample)) +
  geom_density(alpha = 0.2) + 
  labs(
    x = "Insert Size", 
    y = "Density", 
    color = "Sample"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 15),  
    legend.text = element_text(size = 8), 
    legend.title = element_text(size = 15), 
    axis.title = element_text(size = 15),  
    axis.ticks = element_line(linewidth = 1), 
    axis.line = element_line(linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 2)
  )

# ==== writing outputs =====
dir.create(paste0(output_folder,'/sum_reports'))

write.csv(reads_sum_df,
          paste0(output_folder,'/sum_reports/read_counts_sum.csv'),
          quote = F,
          row.names = F)

write.csv(reads_dist_sum_df,
          paste0(output_folder,'/sum_reports/read_distribution_sum.csv'),
          quote = F,
          row.names = F)

write.csv(genebody_sum_df,
          paste0(output_folder,'/sum_reports/genebody_coverage_sum.csv'),
          quote = F,
          row.names = F)

write.csv(insertion_sum_df,
          paste0(output_folder,'/sum_reports/insertion_size_sum.csv'),
          quote = F,
          row.names = F)

pdf(paste0(output_folder,'/sum_reports/read_counts_sum.pdf'),
    width=10, height =7)
print(p1)
dev.off()

pdf(paste0(output_folder,'/sum_reports/read_distribution_sum.pdf'),
    width=6, height =4)
print(p2)
dev.off()

pdf(paste0(output_folder,'/sum_reports/genebody_coverage_sum.pdf'),
    width=4, height =3)
print(p3)
dev.off()

pdf(paste0(output_folder,'/sum_reports/insertion_size_sum.pdf'),
    width=4, height =3)
print(p4)
dev.off()



