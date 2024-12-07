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
### figure 1b

library(ComplexHeatmap)
library(circlize)
expr <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/FPKM.txt.cln", header = T, sep = "\t")
deg <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/iPSC.deseq2.out", header = T, sep = "\t")

VSMC_gene_list <- c("POSTN", "FBLN5", "ACTA2", "NAV2", "PDGFRB")
EC_gene_list <- c("KDR", "EDN1", "TIE1", "ICAM1", "PECAM1", "VWF", "CDH5", "CD34", "ANGPT2", "LYVE1")
VSMC_gene_list <- VSMC_gene_list[VSMC_gene_list %in% deg[deg$padj<0.05, ]$gene_name]
EC_gene_list <- EC_gene_list[EC_gene_list %in% deg[deg$padj<0.05, ]$gene_name]

select_gene_list <- c(VSMC_gene_list, EC_gene_list)
expr_use <- expr[expr$Gene_Name %in% select_gene_list, ]
rownames(expr_use) <- expr_use$Gene_Name
expr_use <- expr_use[select_gene_list, order(colnames(expr)[8:15])+7]

num_ec_markers <- length(EC_gene_list)
num_vsmc_markers <- length(VSMC_gene_list)
row_groups <- factor(c(rep("VSMC markers", num_vsmc_markers),
                       rep("EC markers", num_ec_markers)))

expr_scaled <- t(scale(t(expr_use)))
col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

col_groups <- factor(c(rep("iECs", 4), rep("iVSMCs", 4)))
col_annotation <- HeatmapAnnotation(Group = col_groups,
                                    col = list(Group = c("iECs" = "#E1C971", "iVSMCs" = "#C0C0C0")),
                                    show_annotation_name = FALSE)

### 5.2x4.7
Heatmap(expr_scaled,
        name = "Expression",
        top_annotation = col_annotation,
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_split = row_groups,
        column_split = col_groups,
        show_row_names = TRUE,
        show_column_names = TRUE,
        heatmap_legend_param = list(title = "Scaled\nexpression"))






### figure 1c
rrbs_count <- read.table("/groups/mliang1/pliu1/Mini-ENCODE-for-human-vascular-cells/data-for-figures/Figure1C/CpG.sites.summary.txt", header = T)
ec_anno <- read.table("/groups/mliang1/pliu1/Mini-ENCODE-for-human-vascular-cells/data-for-figures/Figure1C/iPSC.RRBS.dat.iEC.mincov5.5mC.avg.forAnnot.peakAnno.combined.freq", header = F)
vsmc_anno <- read.table("/groups/mliang1/pliu1/Mini-ENCODE-for-human-vascular-cells/data-for-figures/Figure1C/iPSC.RRBS.dat.iVSMC.mincov5.5mC.avg.forAnnot.peakAnno.combined.freq", header = F)

anno_order <- c("Promoter (1-2kb)", "Promoter (<=1kb)", "5' UTR", "Exon", "Intron", "3' UTR", "Downstream (<=300bp)", "Distal Intergenic")
ec_anno$V1 = factor(gsub("_", " ", ec_anno$V1), levels=anno_order)
vsmc_anno$V1 = factor(gsub("_", " ", vsmc_anno$V1), levels=anno_order)


p1 = rrbs_count %>%
  filter(grepl("EC", Sample.id)) %>%
  ggplot(., aes(x = Sample.id, y = CpG.sites.5x)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9, fill = "#E1C971") +
  labs(x = "", y = "Number of CpG sites", title = "iECs") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")


p2 = ggplot(ec_anno, aes(x = V1, y = V3, fill = V1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  labs(x = "", y = "Proportion", title = "Genomic feature (iECs)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")

p1 + p2 + plot_layout(widths = c(1, 2))
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig1.ec.rrbs.anno.png", width=452/96, height=285/96, dpi=300)


p1 = rrbs_count %>%
  filter(grepl("VSMC", Sample.id)) %>%
  ggplot(., aes(x = Sample.id, y = CpG.sites.5x)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9, fill = "#C0C0C0") +
  labs(x = "", y = "Number of CpG sites", title = "iVSMCs") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")


p2 = ggplot(vsmc_anno, aes(x = V1, y = V3, fill = V1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  labs(x = "", y = "Proportion", title = "Genomic feature (iVSMCs)") +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")

p1 + p2 + plot_layout(widths = c(1, 2))
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig1.vsmc.rrbs.anno.png", width=452/96, height=285/96, dpi=300)






### figure 1d
library(ggvenn)

venn_counts <- c(iECs = 64626, 
                 iVSMCs = 169786,
                 overlap = 129500) 

EC_items <- paste0("EC_Item", 1:venn_counts['iECs'])
VSMC_items <- paste0("VSMC_Item", 1:venn_counts['iVSMCs'])
Overlap_items <- paste0("Overlap_Item", 1:venn_counts['overlap'])

EC_set <- c(EC_items, Overlap_items)
VSMC_set <- c(VSMC_items, Overlap_items)
venn_list <- list(iECs = EC_set, iVSMCs = VSMC_set)

ggvenn(
  venn_list, 
  fill_color = c("#E1C971", "#C0C0C0"),
  stroke_size = 0.5, set_name_size = 4
)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/ATAC.venn.png", width=348/96, height=247/96, dpi=300)


venn_counts <- c(EC = 723, 
                 VSMC = 1062,
                 overlap = 14807) 

EC_items <- paste0("EC_Item", 1:venn_counts['EC'])
VSMC_items <- paste0("VSMC_Item", 1:venn_counts['VSMC'])
Overlap_items <- paste0("Overlap_Item", 1:venn_counts['overlap'])

EC_set <- c(EC_items, Overlap_items)
VSMC_set <- c(VSMC_items, Overlap_items)
venn_list <- list(iECs = EC_set, iVSMCs = VSMC_set)

ggvenn(
  venn_list, 
  fill_color = c("#E1C971", "#C0C0C0"),
  stroke_size = 0.5, set_name_size = 4
)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/micro-c.venn.png", width=348/96, height=247/96, dpi=300)




