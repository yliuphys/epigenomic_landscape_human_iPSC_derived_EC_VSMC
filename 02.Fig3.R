library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(scales)
library(stringr)
library(clusterProfiler)


base_font_size <- 12
theme_set(theme_classic(base_size = base_font_size))

setwd("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data")




################################################################################
### figure 3a-c
scatter_mod = function(input_df, x, sample){
  
  x_axis <- x
  x_label <- names(x)
  
  input_df$sample <- sample
  input_df$x <- input_df[, x]
  input_df$FPKM <- log2(input_df$FPKM + 1)
  plot_title <- sample
  
  lm = coef(lm(FPKM ~ x, data = input_df)); intercept = as.numeric(lm[1]); slope = as.numeric(lm[2])
  cor = cor.test(input_df$FPKM, as.numeric(input_df$x), method = "spearman", exact = FALSE)
  cor_est = cor$estimate; cor_p = cor$p.value
  cor_label = paste0("R = ", round(cor_est, 2), ", p-value = ", formatC(cor_p, format = "e", digits = 2))
  
  p = ggplot(input_df, aes(y = FPKM, x = x)) +
    geom_point(aes(color=sample), alpha = 0.5) +
    stat_density_2d(color="black") + # Add 2D density plot
    theme_classic() + 
    theme(
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 12),
      legend.position = "none"
    ) +
    scale_color_manual(values = c("iECs"="#E1C971", "iVSMCs"="#C0C0C0")) +
    labs(x = x_label, y = "log2(FPKM+1)", title = plot_title) +
    annotate(
      geom = "text",
      x = max(as.numeric(input_df$x), na.rm = TRUE) * 0.6, 
      y = max(input_df$FPKM, na.rm = TRUE) * 0.9,
      label = cor_label, 
      color = "red", 
      size = 4.5
    ) +
    geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed") # Add regression line
  
  print(p)
}

data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3A/EC.narrowPeak.peakAnno.promoter.srt.proc.avgexp.proc", header = T, sep = "\t")
scatter_mod(input_df=data, x= c("Chromatin accessibility"="Nearest_TSS"), sample="iECs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/EC_atac_expr_cor.density.png", width=336/96, height=219/96, dpi=300)


data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3A/VSMC.narrowPeak.peakAnno.promoter.srt.proc.avgexp.proc", header = T, sep = "\t")
scatter_mod(input_df=data, x= c("Chromatin accessibility"="Nearest_TSS"), sample="iVSMCs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/VSMC_atac_expr_cor.density.png", width=336/96, height=219/96, dpi=300)


data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3B/promoter_gene_expression_correlation/IEC.avg.pairwise.txt", header = T, sep = "\t")
data$FPKM <- data$IEC.avg.1
scatter_mod(input_df=data, x= c("5mC in promoter"="IEC.avg"), sample="iECs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/EC_promoter_5mc_expr_cor.density.png", width=336/96, height=219/96, dpi=300)


data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3B/promoter_gene_expression_correlation/IVSMC.avg.pairwise.txt", header = T, sep = "\t")
data$FPKM <- data$IVSMC.avg.1
scatter_mod(input_df=data, x= c("5mC in promoter"="IVSMC.avg"), sample="iVSMCs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/VSMC_promoter_5mc_expr_cor.density.png", width=336/96, height=219/96, dpi=300)


data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3B/FPKM.txt.EC.delXYM.enhancer.20kb", header = T, sep = "\t")
data$FPKM <- data$Expression
scatter_mod(input_df=data, x= c("5mC in enhancer\n(within 20kb upstream of TSS)"="EC.avg.m"), sample="iECs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/EC_enhancer_20kb_5mc_expr_cor.density.png", width=336/96, height=219/96, dpi=300)


data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3B/FPKM.txt.VSMC.delXYM.enhancer.20kb", header = T, sep = "\t")
data$FPKM <- data$Expression
scatter_mod(input_df=data, x= c("5mC in enhancer\n(within 20kb upstream of TSS)"="VSMC.avg.m"), sample="iVSMCs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/VSMC_enhancer_20kb_5mc_expr_cor.density.png", width=336/96, height=219/96, dpi=300)




### figure 3d 
promoter_count = function(input_df, gene_ranges=gene_ranges, expr=expr){
  
  sample <- toupper(strsplit(input_df, "_")[[1]][1])
  resolution <- strsplit(input_df, "_")[[1]][2] %>% gsub("kb", " kb", .)
  
  min_samples <- 2
  mean_threshold <- 0.1
  
  if(sample=="EC"){
    expr_cols <- colnames(expr)[grepl("EC", colnames(expr))]
  }else{
    expr_cols <- colnames(expr)[grepl("VSMC", colnames(expr))]
  }
  
  expr_filter <- expr[rowSums(expr[expr_cols]!=0) >= min_samples &
                        rowMeans(expr[expr_cols]) >= mean_threshold, ]
  rownames(expr_filter) = expr_filter$Gene_ID
  expr_filter = expr_filter[, expr_cols]
  
  expr_avg <- data.frame(gene_id = rownames(expr_filter),
                         expr = rowMeans(log2(expr_filter+1)))
  
  input_df <- standardize_column_names(get(input_df))
  
  mc_BIN1 <- GRanges(
    seqnames = input_df$BIN1_CHR,
    ranges = IRanges(start = input_df$BIN1_START, end = input_df$BIN1_END)
  )
  
  mc_BIN2 <- GRanges(
    seqnames = input_df$BIN2_CHR,
    ranges = IRanges(start = input_df$BIN2_START, end = input_df$BIN2_END)
  )
  
  overlap_BIN1 <- findOverlaps(mc_BIN1, gene_ranges)
  overlap_BIN2 <- findOverlaps(mc_BIN2, gene_ranges)
  
  trans_BIN1 <- data.frame(
    BIN_id = queryHits(overlap_BIN1),
    gene_id = gene_ranges$gene_id[subjectHits(overlap_BIN1)]
  ) %>% unique() %>% group_by(gene_id) %>% 
    summarise(count = n())
  
  trans_BIN2 <- data.frame(
    BIN_id = queryHits(overlap_BIN2),
    gene_id = gene_ranges$gene_id[subjectHits(overlap_BIN2)]
  ) %>% unique() %>% group_by(gene_id) %>% 
    summarise(count = n())
  
  trans_expd = data.frame(gene_id = unique(gene_ranges$gene_id)) %>%
    left_join(trans_BIN1, by = c("gene_id")) %>%
    left_join(trans_BIN2, by = c("gene_id")) %>%
    mutate(across(everything(), ~replace_na(., 0)),
           count = count.x + count.y) %>%
    left_join(expr_avg, by = c("gene_id")) %>%
    mutate(
      sample = sample,
      resolution = resolution
    ) %>% dplyr::select(-count.x, -count.y) %>% unique()
  
  return(trans_expd)
  
}


ec_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iEC3MC_8kb_loops.tsv", header = T, sep = "\t")
vsmc_8kb <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iVSMC7MC_8kb_loops.tsv", header = T, sep = "\t")

ec_8kb <- ec_8kb[ec_8kb$FDR<0.2, ]
vsmc_8kb <- vsmc_8kb[vsmc_8kb$FDR<0.2, ]

expr <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/FPKM.txt.cln", header = T, sep = "\t")

expr <- expr %>%
  mutate(
    Promoter_Start = ifelse(Strand == "+", Start - 2000, Start),
    Promoter_End = ifelse(Strand == "+", End, End + 200)
  ) %>%
  mutate(
    Promoter_Start = ifelse(Promoter_Start < 1, 1, Promoter_Start)
  )

gene_ranges <- GRanges(
  seqnames = Rle(paste0("chr", expr$Reference)),
  ranges = IRanges(start = expr$Promoter_Start, end = expr$Promoter_End),
  strand = Rle(expr$Strand),
  gene_id = expr$Gene_ID
)

ec_8kb_gene <- promoter_count("ec_8kb", gene_ranges, expr_ec)
vsmc_8kb_gene <- promoter_count("vsmc_8kb", gene_ranges, expr_vsmc)


boxplot_mod = function(input_df){
  
  sample <- ifelse(grepl("ec", input_df), "iECs", "iVSMCs")
  resolution <- strsplit(input_df, "_")[[1]][2] %>% gsub("kb", " kb", .)
  plot_title <- paste0(sample, " (", resolution, ")\n")
  
  input_df = get(input_df) %>%
    mutate(resolution = factor(resolution, levels=c("4 kb", "8 kb", "16 kb")),
           count_fct = ifelse(count >= 5, ">=5", as.character(count)),
           count_fct = factor(count_fct, levels = c(as.character(0:4), ">=5")))
  
  lm = coef(lm(expr ~ count, data = input_df)); intercept=as.numeric(lm[1]); slope=as.numeric(lm[2])
  cor=cor.test(input_df$expr,as.numeric(input_df$count), method = "spearman", exact = FALSE); cor_est=cor$estimate; cor_p=cor$p.value 
  cor_label = paste0("R = ", round(cor_est, 2), ", p-value = ", formatC(cor_p, format = "e", digits = 2))
  p = ggplot(input_df[!is.na(input_df$expr),], aes(y=expr, x=count_fct, fill=sample)) +
    geom_boxplot(outlier.shape=NA) +
    theme_classic() + 
    theme(axis.text = element_text(color="black", size=10),
          axis.title = element_text(color="black", size=12),
          legend.position = "none") +
    scale_fill_manual(values = c("EC"="cornsilk1", "VSMC"="grey90")) +
    labs(x = "Gene interaction count", y = "log2(FPKM+1)", title = plot_title) +
    scale_x_discrete(labels = function(x) {
      ifelse(x %in% c("0", "1", "2", "3", "4", ">=5"), x, "")
    }) +
    annotate(geom="text", x=max(as.numeric(input_df$count_fct), na.rm = TRUE)*0.6, y=max(input_df$expr, na.rm = TRUE)*0.8, 
             label = cor_label, color = "red", size = 4.5) +
    geom_abline(intercept = intercept, slope = slope, color = "red", linetype = "dashed")
  
  print(p)
  
}

boxplot_mod("ec_8kb_gene")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/microc-ec_8kb_gene_count_expr_cor.png", width=336/96, height=219/96, dpi=300)

boxplot_mod("vsmc_8kb_gene")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/microc-vsmc_8kb_gene_count_expr_cor.png", width=336/96, height=219/96, dpi=300)

