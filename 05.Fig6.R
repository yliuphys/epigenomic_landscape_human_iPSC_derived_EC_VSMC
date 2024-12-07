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
### figure 6a
library(ComplexUpset)

expand_lengths <- list(expression_ec=5000, expression_vsmc=5000,
                       rrbs_ec=2000, rrbs_vsmc=2000,
                       atac_ec=2000, atac_vsmc=2000,
                       microccontacts_ec=0, microccontacts_vsmc=0,
                       microcloops_ec=0, microcloops_vsmc=0)

target_trait <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gwas_catalog_target_trait.snp.txt", header = T, sep = "\t")
target_trait <- unique(target_trait[, c("SNP", "POS")])

gr_list <- readRDS("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gr_list.rds")

gr_list[['microccontacts_ec']] <- NULL
gr_list[['microccontacts_vsmc']] <- NULL

snp_to_gr <- function(snp_df) {
  
  snp_df <- snp_df %>%  
    dplyr::select(SNP, POS) %>% unique() %>%
    mutate(CHR = gsub("_.+", "", POS),
           POS_START = as.numeric(gsub(".+_", "", POS)),
           POS_END = POS_START)
  
  snp_df <- snp_df[!is.na(snp_df$POS_START) & !is.na(snp_df$POS_END), ]
  GRanges(seqnames = snp_df$CHR, ranges = IRanges(start = snp_df$POS_START, end = snp_df$POS_END), SNP = snp_df$SNP)
}

target_gr <- snp_to_gr(target_trait)

covered_snp_list <- list()
for (assay_name in names(gr_list)) {
  
  assay_gr <- gr_list[[assay_name]]
  expand_length <- ifelse(assay_name %in% names(expand_lengths), expand_lengths[[assay_name]], 0)
  
  expanded_gr <- resize(assay_gr, width(assay_gr) + 2*expand_length, fix = "center")
  
  target_covered_hits <- findOverlaps(target_gr, expanded_gr)
  target_covered_snp <- target_trait$SNP[queryHits(target_covered_hits)]
  
  covered_snp_list[[assay_name]] <- unique(target_covered_snp)
}


all_snps <- unique(unlist(lapply(covered_snp_list, function(x) x)))
snp_matrix <- data.frame(SNP = all_snps)
for (key in names(covered_snp_list)) {
  assay_column <- rep(0, length(all_snps))
  covered_snps <- covered_snp_list[[key]]
  assay_column[which(snp_matrix$SNP %in% covered_snps)] <- 1
  snp_matrix[[key]] <- assay_column
}
rownames(snp_matrix) <- snp_matrix$SNP
snp_matrix_binary <- snp_matrix[,-1]


colnames(snp_matrix_binary) <- colnames(snp_matrix_binary) %>%
  gsub("expression_ec", "RNA-seq (iECs)", .) %>%
  gsub("expression_vsmc", "RNA-seq (iVSMCs)", .) %>%
  gsub("rrbs_ec", "RRBS (iECs)", .) %>%
  gsub("rrbs_vsmc", "RRBS (iVSMCs)", .) %>%
  gsub("atac_ec", "ATAC-seq (iECs)", .) %>%
  gsub("atac_vsmc", "ATAC-seq (iVSMCs)", .) %>%
  gsub("microccontacts_ec", "Micro-C contacts (iECs)", .) %>%
  gsub("microccontacts_vsmc", "Micro-C contacts (iVSMCs)", .) %>%
  gsub("microcloops_ec", "Micro-C loops (iECs)", .) %>%
  gsub("microcloops_vsmc", "Micro-C loops (iVSMCs)", .)


# set_size(8, 3)
size = get_size_mode('exclusive_intersection')

p = upset(snp_matrix_binary,
          intersect = c("RNA-seq (iECs)", "RNA-seq (iVSMCs)", "RRBS (iECs)", "RRBS (iVSMCs)", "ATAC-seq (iECs)", "ATAC-seq (iVSMCs)", "Micro-C loops (iECs)", "Micro-C loops (iVSMCs)"),
          name = "",
          base_annotations=list(
            'Number of SNPs'=intersection_size(
              text_mapping=aes(
                label=paste0(
                  !!size
                  # '\n(', !!upset_text_percentage(), ')'
                ),
                colour='on_background', # ifelse(!!size > 1000, 'on_bar', 'on_background'),
                y=!!size
              )
            ) + ylim(c(0, 1800))
          ),
          width_ratio=0.1, min_size = 40,
          stripes=c('cornsilk1', 'grey90')
) &
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
  )

print(p)

ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/snp_upset.png", width=1000/96, height=388/96, dpi=300)






### figure 6b
calculate_snp_coverage <- function(target_trait, control_trait, gr_list, expand_lengths = list()) {
  
  snp_to_gr <- function(snp_df) {
    
    snp_df <- snp_df %>% 
      dplyr::select(SNP, POS) %>% unique() %>%
      mutate(CHR = gsub("_.+", "", POS),
             POS_START = as.numeric(gsub(".+_", "", POS)),
             POS_END = POS_START)
    
    snp_df <- snp_df[!is.na(snp_df$POS_START) & !is.na(snp_df$POS_END), ]
    GRanges(seqnames = snp_df$CHR, ranges = IRanges(start = snp_df$POS_START, end = snp_df$POS_END))
  }
  
  control_gr <- snp_to_gr(control_trait)
  
  results <- data.frame(
    trait = character(), total_trait_snp = integer(), covered_trait_snp = integer(),
    assay = character(), cell_type = character(),
    control_snp = integer(), covered_control_snp = integer(),
    p_value = numeric(), odds_ratio = numeric(), stringsAsFactors = FALSE
  )
  
  for (trait in target_trait_list) {
    
    target_gr <- target_trait %>% filter(TRAIT == trait) %>%
      snp_to_gr()
    
    for (assay_name in names(gr_list)) {
      
      assay_gr <- gr_list[[assay_name]]
      expand_length <- ifelse(assay_name %in% names(expand_lengths), expand_lengths[[assay_name]], 0)
      
      expanded_gr <- resize(assay_gr, width(assay_gr) + 2*expand_length, fix = "center")
      
      assay_parts <- unlist(strsplit(assay_name, "_"))
      assay <- assay_parts[1]
      cell_type <- assay_parts[2]
      
      target_covered <- countOverlaps(target_gr, expanded_gr)
      num_target_covered <- sum(target_covered > 0)
      num_target_total <- length(target_covered)
      
      control_covered <- countOverlaps(control_gr, expanded_gr)
      num_control_covered <- sum(control_covered > 0)
      num_control_total <- length(control_covered)
      
      fisher_table <- matrix(c(
        num_target_covered, num_target_total - num_target_covered,
        num_control_covered, num_control_total - num_control_covered
      ), nrow = 2)
      
      fisher_test <- fisher.test(fisher_table)
      
      results <- rbind(results, data.frame(
        trait = trait, 
        total_trait_snp = num_target_total,
        covered_trait_snp = num_target_covered,
        covered_trait_ratio = num_target_covered / num_target_total,
        assay = assay,
        cell_type = cell_type,
        control_snp = num_control_total,
        covered_control_snp = num_control_covered,
        covered_control_ratio = num_control_covered / num_control_total,
        p_value = fisher_test$p.value,
        odds_ratio = fisher_test$estimate,
        significance = case_when(
          fisher_test$p.value < 0.001 ~ "***", 
          fisher_test$p.value < 0.01  ~ "**", 
          fisher_test$p.value < 0.05  ~ "*", 
          TRUE                        ~ ""
        )
      ))
    }
  }
  
  return(results)
}

target_trait <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gwas_catalog_target_trait.snp.txt", header = T, sep = "\t")
control_trait <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gwas_catalog_control_trait.snp.txt", header = T, sep = "\t")

target_trait_list <- unique(target_trait$TRAIT)

gr_list <- readRDS("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gr_list.rds")

expand_lengths <- list(expression_ec=5000, expression_vsmc=5000,
                       rrbs_ec=2000, rrbs_vsmc=2000,
                       atac_ec=2000, atac_vsmc=2000,
                       microccontacts_ec=0, microccontacts_vsmc=0,
                       microcloops_ec=0, microcloops_vsmc=0)

result_df <- calculate_snp_coverage(target_trait, control_trait, gr_list, expand_lengths)

result_df$assay <- recode(result_df$assay,
                          "rrbs" = "RRBS",
                          "microccontacts" = "Micro-C (contacts)",
                          "microcloops" = "Micro-C (loops)",
                          "expression" = "RNA-seq",
                          "atac" = "ATAC-seq")
result_df$cell_type <- recode(result_df$cell_type,
                              "ec" = "iECs",
                              "vsmc" = "iVSMCs")
result_df$trait <- factor(result_df$trait, levels = target_trait_list)
result_df$p_adj <- p.adjust(result_df$p_value, method = "fdr")

write.table(result_df, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/SNP_coverage_by_assay.out", col.names = T, row.names = F, sep = "\t", quote = F)



result_df <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/SNP_coverage_by_assay.out", header = T, sep = "\t")
result_df$trait <- factor(result_df$trait, levels = target_trait_list)

p = ggplot(result_df,
           aes(x = trait, # fct_reorder(Trait, Coverage_Rate), 
               y = assay, size=covered_trait_ratio, color=-log10(p_value))) +
  geom_point() + 
  geom_point(data = result_df %>% filter(p_adj < 0.05), 
             aes(x = trait, 
                 y = assay), 
             shape = 21, fill = NA, color = "black", stroke = 1) + 
  theme_classic() + 
  theme(strip.text = element_text(color="black", size = 10),
        axis.text = element_text(color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  scale_color_gradient(low = "white", high = "red") +
  labs(subtitle = "Percent of SNPs covered by assays\n(compared to control traits)", x = "", y="",
       color="-log(p-value)", size="Coverage\npercent") +
  facet_grid(cell_type~.)

print(p)

ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/SNP_coverage_by_assay.v2.png", width=517/96, height=462/96, dpi=300)







### figure 6c

target_trait <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gwas_catalog_target_trait.snp.txt", header = T, sep = "\t")
control_trait <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gwas_catalog_control_trait.snp.txt", header = T, sep = "\t")

snp_to_gr <- function(snp_df) {
  
  snp_df <- snp_df %>% 
    mutate(CHR = gsub("_.+", "", POS),
           POS_START = as.numeric(gsub(".+_", "", POS))-5000,
           POS_END = as.numeric(gsub(".+_", "", POS))+5000)
  
  snp_df <- snp_df[!is.na(snp_df$POS_START) & !is.na(snp_df$POS_END), ]
  GRanges(seqnames = snp_df$CHR, ranges = IRanges(start = snp_df$POS_START, end = snp_df$POS_END))
}

control_gr <- snp_to_gr(control_trait)
target_gr <- snp_to_gr(target_trait)


expr <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/FPKM.txt.cln", header = T, sep = "\t")
gene_ranges <- GRanges(
  seqnames = Rle(paste0("chr", expr$Reference)),
  ranges = IRanges(start = expr$Start, end = expr$End),
  strand = Rle(expr$Strand),
  gene_id = expr$Gene_ID
)

overlap_control <- findOverlaps(control_gr, gene_ranges)
overlap_target <- findOverlaps(target_gr, gene_ranges)

mcols(gene_ranges)$type = "Remaining"
mcols(gene_ranges)[subjectHits(overlap_control),]$type = "Within 5kb of SNPs of control traits"
mcols(gene_ranges)[subjectHits(overlap_target),]$type = "Within 5kb of SNPs of target traits"

gene_type = gene_ranges$type
names(gene_type) = gene_ranges$gene_id

min_samples <- 2
mean_threshold <- 0.1

ec_cols <- colnames(expr)[grepl("EC", colnames(expr))]
expr_ec <- data.frame(FPKM_mean = rowMeans(log2(expr[, ec_cols]+1)),
                      type = gene_type[expr$Gene_ID])

vsmc_cols <- colnames(expr)[grepl("VSMC", colnames(expr))]
expr_vsmc <- data.frame(FPKM_mean = rowMeans(log2(expr[, vsmc_cols]+1)),
                        type = gene_type[expr$Gene_ID])



input_df = expr_ec
input_df$type <- factor(input_df$type, levels = c("Within 5kb of SNPs of target traits",
                                                  "Within 5kb of SNPs of control traits", 
                                                  "Remaining"), 
                        labels = c("Target SNPs", "Control SNPs", "Remaining"))
p_val = kruskal.test(FPKM_mean ~ type, data = input_df)$p.value
p_label = paste0("p-value = ", formatC(p_val, format = "e", digits = 2))
ggplot(input_df, aes(x=FPKM_mean, color=type)) +
  geom_step(aes(y=after_stat(y)),stat="ecdf") +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=10),
        axis.title = element_text(color="black", size=12),
        legend.position = "right") +
  labs(x = "log2(FPKM+1)", y = "Cumulative probability", color = "Gene proximity to SNPs", title = "iECs") +
  annotate(geom="text", x=max(as.numeric(input_df$FPKM_mean), na.rm = TRUE)*0.6, y=0.2, 
           label = p_label, color = "red", size = 4.5)

ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/EC.expr_edfc.png", width=454/96, height=238/96, dpi=300)


input_df = expr_vsmc
input_df$type <- factor(input_df$type, levels = c("Within 5kb of SNPs of target traits",
                                                  "Within 5kb of SNPs of control traits", 
                                                  "Remaining"), 
                        labels = c("Target SNPs", "Control SNPs", "Remaining"))
p_val = kruskal.test(FPKM_mean ~ type, data = input_df)$p.value
p_label = paste0("p-value = ", formatC(p_val, format = "e", digits = 2))
ggplot(input_df, aes(x=FPKM_mean, color=type)) +
  geom_step(aes(y=after_stat(y)),stat="ecdf") +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=10),
        axis.title = element_text(color="black", size=12),
        legend.position = "right") +
  labs(x = "log2(FPKM+1)", y = "Cumulative probability", color = "Gene proximity to SNPs", title = "iVSMCs") +
  annotate(geom="text", x=max(as.numeric(input_df$FPKM_mean), na.rm = TRUE)*0.6, y=0.2, 
           label = p_label, color = "red", size = 4.5)

ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/VSMC.expr_edfc.png", width=454/96, height=238/96, dpi=300)





### figure 6d
deg = read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/iPSC.deseq2.out", header = T, sep = "\t")

deg$type = gene_type[deg$gene_id]

top_genes <- deg %>%
  filter(type == "Within 5kb of SNPs of target traits") %>%
  arrange(padj) %>%
  slice(1:30)

deg %>%
  filter(type=="Within 5kb of SNPs of target traits",
         ! is.na(padj)) %>%
  mutate(Significant = ifelse(padj<0.05, "Yes", "No")) %>%
  ggplot(., 
         aes(x=log2FoldChange, y=-log10(pvalue))) +
  geom_point(aes(color=Significant)) +
  theme_classic() + 
  theme(axis.text = element_text(color="black", size=10),
        axis.title = element_text(color="black", size=12),
        legend.position = "right") +
  scale_color_manual(values = c("grey", "red")) +
  # lims(y=c(0, 250)) +
  labs(x = "log2 fold change", y = "-log10(p-value)", color = "DEG", 
       title = "Expression of genes near target SNPs by cell type\n(iVSMCs vs. iECs)") +
  geom_label_repel(data = top_genes, aes(label = gene_name), size = 3, max.overlaps = Inf)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/deg.volcano.png", width=606/96, height=317/96, dpi=300)





### SNP nomination
tbl <- table(snp_matrix[, c("atac_ec", "microcloops_ec")])
fisher.test(tbl, alternative = "greater")

tbl <- table(snp_matrix[, c("atac_vsmc", "microcloops_vsmc")])
fisher.test(tbl, alternative = "greater")


### calculate nearest gene for each SNP
pc_gene_ranges = gene_ranges[gene_ranges$gene_type=="protein_coding",]
nearest_indices <- nearest(target_gr, pc_gene_ranges)
nearest_genes <- pc_gene_ranges[nearest_indices]

distances <- distance(target_gr, nearest_genes)

results <- data.frame(
  SNP = mcols(target_gr)$SNP,
  Gene_ID = mcols(nearest_genes)$gene_id,
  Gene_Name = mcols(nearest_genes)$gene_name,
  Gene_Type = mcols(nearest_genes)$gene_type,
  Gene_Chr = mcols(nearest_genes)$gene_chr,
  Gene_Start = mcols(nearest_genes)$gene_start,
  Gene_End = mcols(nearest_genes)$gene_end,
  Gene_Strand = mcols(nearest_genes)$gene_strand,
  Distance = distances
)

head(results)
results[results$SNP=="rs9833313", ]
snp_100kb = results[results$Distance>100000, ]$SNP

### extracting SNPs having a protein coding gene in the same loop
# gr_list <- readRDS("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gr_list.rds")
# pc_gene_ranges
# target_gr
loops_ec <- gr_list$microcloops_ec
loops_vsmc <- gr_list$microcloops_vsmc

find_snps_genes_in_loops <- function(snp_gr, gene_gr, loop_gr) {
  
  snp_overlap <- findOverlaps(snp_gr, loop_gr)
  snp_loop_df <- data.frame(
    SNP_ID = snp_gr$SNP[queryHits(snp_overlap)],
    SNP_Location = paste0(seqnames(snp_gr[queryHits(snp_overlap)]), ":", start(snp_gr[queryHits(snp_overlap)]), "-", end(snp_gr[queryHits(snp_overlap)])),
    Loop_Chr = seqnames(loop_gr[subjectHits(snp_overlap)]),
    Loop_Start = start(loop_gr[subjectHits(snp_overlap)]),
    Loop_End = end(loop_gr[subjectHits(snp_overlap)])
  )
  
  gene_overlap <- findOverlaps(gene_gr, loop_gr, type = "within")
  gene_loop_df <-   data.frame(
    Gene_ID = gene_gr$gene_id[queryHits(gene_overlap)],
    Gene_Name = gene_gr$gene_name[queryHits(gene_overlap)],
    Gene_Ranges = paste0(seqnames(gene_gr[queryHits(gene_overlap)]), ":", start(gene_gr[queryHits(gene_overlap)]), "-", end(gene_gr[queryHits(gene_overlap)])),
    Loop_Chr = seqnames(loop_gr[subjectHits(gene_overlap)]),
    Loop_Start = start(loop_gr[subjectHits(gene_overlap)]),
    Loop_End = end(loop_gr[subjectHits(gene_overlap)])
  )
  
  combined_df <- merge(snp_loop_df, gene_loop_df, by = c("Loop_Chr", "Loop_Start", "Loop_End"))
  
  combined_df %>%
    arrange(Loop_Chr, Loop_Start, Loop_End)
  
}


ec_snp_gene_res <- find_snps_genes_in_loops(target_gr, pc_gene_ranges, loops_ec)
vsmc_snp_gene_res <- find_snps_genes_in_loops(target_gr, pc_gene_ranges, loops_vsmc)



snp_cat_df <- enframe(covered_snp_list, name = "Category", value = "SNP") %>%
  unnest(SNP)

SNP_counts <- snp_cat_df %>%
  group_by(SNP) %>%
  summarise(CategoryCount = n_distinct(Category),
            Categories = paste(unique(Category), collapse = ", ")
  )


target_trait <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gwas_catalog_target_trait.snp.txt", header = T, sep = "\t")
BP_SNP_list <- target_trait[target_trait$TRAIT %in% c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension"), ]$SNP
dim(SNP_counts)
# [1] 7574    3
dim(SNP_counts[SNP_counts$SNP %in% BP_SNP_list, ])
# [1] 6050    3
dim(SNP_counts[SNP_counts$Categories=="atac_ec, atac_vsmc, microcloops_ec, microcloops_vsmc", ])
# [1] 240   3

### load LD information
snps_LD = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/SNPS_in_LD.txt", header = T, sep = "\t")
ld_1 = names(table(snps_LD$rsID1)[table(snps_LD$rsID1)==1])

dim(SNP_counts[(SNP_counts$SNP %in% BP_SNP_list) & (SNP_counts$SNP %in% ld_1) & (SNP_counts$SNP %in% snp_100kb) &
                 (SNP_counts$Categories=="atac_ec, atac_vsmc, microcloops_ec, microcloops_vsmc"), ])
# [1] 10  3



### nomination pipeline
library(DiagrammeR)

filter_counts <- list()
snp_filter <- SNP_counts
filter_counts[["Initial"]] <- nrow(snp_filter)

snp_filter <- snp_filter %>%
  filter(SNP %in% BP_SNP_list)
filter_counts[["BP associated"]] <- nrow(snp_filter)

snp_filter <- snp_filter %>%
  filter(Categories=="atac_ec, atac_vsmc, microcloops_ec, microcloops_vsmc")
filter_counts[["Accessible SNPs located within Micro-C loops\n(in both EC and VSMC)"]] <- nrow(snp_filter)

snp_filter <- snp_filter %>%
  filter(SNP %in% ld_1)
filter_counts[["SNPs isolated within LD blocks"]] <- nrow(snp_filter)

snp_filter <- snp_filter %>%
  filter(SNP %in% snp_100kb)
filter_counts[["Distal regulatory element (>100 kb)"]] <- nrow(snp_filter)

snp_filter <- snp_filter %>%
  filter((SNP %in% ec_snp_gene_res$SNP_ID) & (SNP %in% vsmc_snp_gene_res$SNP_ID))
filter_counts[["SNPs sharing loops with protein coding genes"]] <- nrow(snp_filter)

snp_filter[snp_filter$SNP=="rs9833313",]
print(filter_counts)

filter_data = data.frame(Step = names(filter_counts),
                         SNP_Gene_Count = as.numeric(filter_counts))

generate_grViz <- function(filter_data, ranksep = 0.25, nodesep = 0.5) {
  nodes <- paste0(filter_data$Step, ":\n", filter_data$SNP_Gene_Count, " SNPs")
  edges <- paste0(seq(1, nrow(filter_data) - 1), " -> ", seq(2, nrow(filter_data)))
  
  diagram_code <- paste0(
    "digraph G {
       graph [layout = dot, rankdir = TB, ranksep = ", ranksep, ", nodesep = ", nodesep, "]
       node [shape = box, style = filled, fillcolor = lightblue, fontname = Helvetica]
       ",
    paste0(paste0(seq(1, nrow(filter_data)), " [label = \"", nodes, "\"]"), collapse = "\n"),
    "\n",
    paste0(edges, collapse = "\n"),
    "\n}"
  )
  grViz(diagram_code)
}

generate_grViz(filter_data)
