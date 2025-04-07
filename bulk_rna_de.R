#!/usr/bin/env Rscript
# bulk_rna_de.R

# loading packages
suppressPackageStartupMessages({
  library(optparse)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(edgeR)
  library(ggrepel)
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
              metavar="SAMPLES"),
  make_option(c("--conditions"), 
              type="character", 
              default=NULL, 
              help="Space-separated list of conditions", 
              metavar="CONDITIONS"),
  make_option(c("--covar"), 
              type="character", 
              default=NULL, 
              help="Space-separated list of covariates (optional)", 
              metavar="COVAR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

samples <- strsplit(opt$samples, " ")[[1]]
conditions <- str_split(opt$conditions, ",")[[1]]
conditions <- as.character(conditions)
# Check if covar is provided, otherwise default to NULL
if (!is.null(opt$covar)) {
  covar <- str_split(opt$covar, ",")[[1]]
  covar <- as.character(covar)
} else {
  covar <- NULL
}

output_folder <- opt$out

# samples=c('BJ', '50K_1', '50K_2', 'all9')
# conditions=c('ctrl', 'hspc', 'hspc', 'hspc')
# output_folder <- '/xdisk/darrenc/jiachengd/20250403_ledford_lab_mouse_lung_10x/03_bulk_rna_anlysis'
# covar <- NULL

# ==== running edgeR DE test ====
count_table <- read.table(paste0(output_folder, '/featureCounts/combined_counts.txt'),
                          header =T,
                          check.names = F)
rownames(count_table) <- count_table$Geneid
count_table$Geneid <- NULL

# create metadata
if (length(covar) == 0){
  meta_data <- data.frame(
    sample = samples,
    condition = conditions)
} else {
  meta_data <- data.frame(
    sample = samples,
    condition = conditions,
    covariant = covar)
  meta_data <- meta_data %>%
    dplyr::mutate(covariant = factor(covariant, sort(unique(covariant))))
  
}

meta_data <- meta_data %>%
  dplyr::mutate(condition = factor(condition, 
                                   levels = c('ctrl',
                                              sort(setdiff(unique(condition), 'ctrl')))))
rownames(meta_data) <- meta_data$sample
meta_data <- meta_data[colnames(count_table),]

if (length(covar) == 0){
  design <- model.matrix(~meta_data$condition)
} else {
  design <- model.matrix(~meta_data$condition + meta_data$covar)
}

# create dgelist object
dpb <- DGEList(counts=count_table,group=meta_data$condition)
keep <- filterByExpr(
  dpb,
  group = meta_data$condition
)
dpb <- dpb[keep, , keep.lib.sizes = FALSE]
dpb <- calcNormFactors(dpb)

# plot heatmap cor
logcounts <- cpm(dpb, log=TRUE)
rld_cor <- cor(logcounts)
p1 <- pheatmap(rld_cor,
              annotation = meta_data[,c('condition'), drop=F],
              cluster_cols=T,
              cluster_rows=T,
              cellwidth=10,
              cellheight=10)

# plot pca
mds_result <- plotMDS(cpm(dpb, log=TRUE),gene.selection = 'common',plot = F,top=1000)
mds_coords <- data.frame(
  dim_1 = mds_result$x,
  dim_2 = mds_result$y
)
mds_coords$sample <- colnames(dpb)
mds_coords <- merge(mds_coords, meta_data, by = 'sample')

dim_1_var_exp <- round(mds_result$var.explained[1] * 100, 2)
dim_2_var_exp <- round(mds_result$var.explained[2] * 100, 2)

p2 <- ggplot(mds_coords, aes(x = dim_1, 
                            y = dim_2, 
                            color = condition, 
                            label = sample
)
) +
  geom_point(size = 3, alpha = 1) +
  geom_text_repel(size = 4,            # Adjust label size
                  max.overlaps = 15,   # Increase to allow more labels
                  box.padding = 0.5,   # Increase for more space around labels
                  point.padding = 0.2, # Adjust distance between label and point
                  min.segment.length = 0.1,  # Controls segment visibility
                  segment.color = "grey50",  # Change segment color
                  nudge_y = 0.1) +  # Slight upward shift
  theme_classic() +
  labs(
    title = 'Top 1000 features',
    x = paste0("PC_1 (", dim_1_var_exp, "%)"),
    y = paste0("PC_2 (", dim_2_var_exp, "%)")
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text = element_text(size = 9),  
    legend.text = element_text(size = 9), 
    legend.title = element_text(size = 9), 
    axis.title = element_text(size = 9),  
    axis.ticks = element_line(linewidth = 1), 
    axis.line = element_line(linewidth = 1),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 2)
  )


# perform de test
dpb <- estimateDisp(dpb)
fit <- glmFit(dpb,design = design)

non_ctrl_condition <- length(setdiff(meta_data$condition,'ctrl'))

result_list <- vector('list',non_ctrl_condition)
for (i in 2:(non_ctrl_condition+1)){
  qlf <- glmLRT(fit,coef = i)
  compared_sample <- colnames(design)[i]
  comparison <- paste0(compared_sample, '_vs_ctrl')
  result <- as.data.frame(topTags(qlf, n=nrow(qlf$table)))
  result$comparison <- comparison
  result_list[[i]] <- result
}
result_df <- do.call(rbind, result_list)

# ==== writing outputs =====
dir.create(paste0(output_folder,'/de_test'))
write.csv(result_df,
          paste0(output_folder,'/de_test/de_test_result.csv'),
          quote = F,
          row.names = T)

pdf(paste0(output_folder,'/de_test/cor_heatmap.pdf'),
    width=5, height =5)
print(p1)
dev.off()

pdf(paste0(output_folder,'/de_test/pca.pdf'),
    width=5, height =5)
print(p2)
dev.off()






