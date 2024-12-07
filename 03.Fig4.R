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
### figure 4a-c, f-i
data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure4A_D/EPT.summary.txt", header = T, sep = "\t")

data <- data[grepl("8kb_loops.tsv.proc.BIN", data$Sample), ]

data_long <- data %>%
  pivot_longer(cols = -Sample, names_to = "Category", values_to = "Count") %>%
  mutate(type = ifelse(grepl("EC", Sample), "iECs", "iVSMCs"))

data_long %>%
  filter(Category %in% c("None", "One", "Both")) %>%
  mutate(Category = factor(Category, levels = c("None", "One", "Both"))) %>%
  group_by(type) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ggplot(., aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 1.05), color = "black", size = 3) +
  labs(x = "", y = "Count", title = "Regulatory elements in chromatin contacts") +
  scale_fill_brewer(palette = "YlOrRd") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), 
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~type)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig4a.png", width=405/96, height=250/96, dpi=300)



data_long %>%
  filter(! (Category %in% c("None", "One", "Both"))) %>%
  # mutate(Category = factor(Category, levels = c("None", "One", "Both"))) %>%
  group_by(type) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ggplot(., aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_stack(vjust = 1.05), color = "black", size = 3) +
  labs(x = "", y = "Count", title = "Interaction of regulatory elements") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), 
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~type)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig4b.png", width=534/96, height=250/96, dpi=300)




data_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3C/EC.FDR0.20.None.One.Both.comparison.wilcox.out.proc", header = T, sep = "\t")
data_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure3C/VSMC.FDR0.20.None.One.Both.comparison.wilcox.out.proc", header = T, sep = "\t")
data_ec$cell_type = "iECs"
data_vsmc$cell_type = "iVSMCs"
data = rbind(data_ec, data_vsmc)
data = data[data$resolution=="8kb", ]

data %>%
  mutate(reg = factor(reg, levels = c("None", "One", "Both")),
         resolution = factor(resolution, levels = c("4kb", "8kb", "16kb")),
         x.1 = case_when(
           reg == "None" ~ 1,  # Position of "None" on x-axis
           reg == "One" ~ 1,   # Position of "One" on x-axis
           reg == "Both" ~ 2   # Comparison between "None" and "Both" should start at "None"
         ),
         x.2 = case_when(
           reg == "None" ~ 2,  # Position of "One" on x-axis (comparing None vs. One)
           reg == "One" ~ 3,   # Position of "Both" on x-axis (comparing One vs. Both)
           reg == "Both" ~ 3   # End at "Both" for None vs. Both
         ),
         y_position = case_when(
           reg == "None" ~ max(m + se) + 0.1,  # Offset for the None vs. One comparison
           reg == "One"  ~ max(m + se) + 0.3,  # Offset for the One vs. Both comparison
           reg == "Both" ~ max(m + se) + 0.5   # Offset for the None vs. Both comparison
         )
  ) %>%
  ggplot(., aes(x = reg, y = m, fill = reg)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2, position = position_dodge(0.7)) +
  geom_text(aes(x = (x.1 + x.2) / 2, y = y_position, label = paste0("p=", signif(p, 2))),
            size = 3, color = "black", vjust = -0.5) +
  geom_segment(aes(x = x.1, xend = x.2, y = y_position, yend = y_position), 
               color = "black", linewidth = 0.5) +
  labs(x = "", y = "Gene expression\n(log2(FPKM+1))", title = "", fill="Regulatory\nelement\noccurrence") +
  lims(y = c(0, 2.2)) +
  scale_fill_brewer(palette = "YlOrRd") +
  theme_classic() +
  theme(axis.text = element_text(color="black"), 
        # axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  facet_wrap(~cell_type)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/EC_VSMC_8kb_reg_expr.bar_plot.png", width=405/96, height=250/96, dpi=300)



reg_barplot <- function(data, ylim = c(0, 1.6), y_label = "Methylation level") {
  p = data %>%
    mutate(
      x.1 = 1,
      x.2 = case_when(
        reg == "EE" ~ 2,
        reg == "EP" ~ 3,
        reg == "ET" ~ 4,
        reg == "PP" ~ 5,
        reg == "PT" ~ 6,
        reg == "TT" ~ 7,
        reg == "Ctrl" ~ 1
      ),
      y_position = case_when(
        reg == "EE" ~ 1.8,# 1.8, # 1.2, # max(m + se, na.rm = TRUE) + 0.1,
        reg == "EP" ~ 2,# 2, # 1.5, # max(m + se, na.rm = TRUE) + 0.3,
        reg == "ET" ~ 2.2,# 2.2, # 1.8, # max(m + se, na.rm = TRUE) + 0.5,
        reg == "PP" ~ 2.4,# 2.4, # 2.1, # max(m + se, na.rm = TRUE) + 0.7,
        reg == "PT" ~ 2.6,# 2.6, # 2.4, # max(m + se, na.rm = TRUE) + 0.9,
        reg == "TT" ~ 2.8,# 2.8, # 2.7, # max(m + se, na.rm = TRUE) + 1.1,
        reg == "Ctrl" ~ NA
      )
    ) %>%
    mutate(Cell.type = ifelse(grepl("EC", Cell.type), "iECs", "iVSMCs")) %>%
    ggplot(aes(x = reg, y = m, fill = reg)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2, position = position_dodge(0.9)) +
    geom_text(data = . %>% filter(!is.na(x.1) & !is.na(x.2) & p < 0.05),
              aes(x = (x.1 + x.2) / 2, y = y_position, label = paste0("p=", signif(p, 2))),
              size = 3, color = "black", vjust = -0.5) +
    geom_segment(data = . %>% filter(!is.na(x.1) & !is.na(x.2) & p < 0.05),
                 aes(x = x.1, xend = x.2, y = y_position, yend = y_position), 
                 color = "black", linewidth = 0.5) +
    labs(x = "", y = y_label, title = "", fill = "Regulatory\nelement\noccurrence") +
    lims(y = ylim) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), 
          legend.position = "None") +
    facet_wrap(~Cell.type)
  print(p)
}

data_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5A_5B/iEC.FDR0.20.methylation.new.out.proc", header = T, sep = "\t")
data_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5A_5B/iVSMC.FDR0.20.methylation.new.out.proc", header = T, sep = "\t")
data = rbind(data_ec, data_vsmc)
data = data[data$resolution=="8kb", ]
reg_barplot(data, ylim = c(0, 1.6), y_label = "Methylation level")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig4c.meth_reg.png", width=534/96, height=296/96, dpi=300)


data_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5C_5D/iEC.FDR0.20.signalValue.new.out.proc", header = T, sep = "\t")
data_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5C_5D/iVSMC.FDR0.20.signalValue.new.out.proc", header = T, sep = "\t")
data = rbind(data_ec, data_vsmc)
data = data[data$resolution=="8kb", ]
reg_barplot(data, ylim = c(0, 2.8), y_label = "Chromatin accessibility\n(log2(Signal+1))")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig4g.atac_reg.png", width=534/96, height=296/96, dpi=300)


data_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5G_5H/iEC.FDR0.20.gene_expression.new.out.proc", header = T, sep = "\t")
data_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5G_5H/iVSMC.FDR0.20.gene_expression.new.out.proc", header = T, sep = "\t")
data = rbind(data_ec, data_vsmc)
data = data[data$resolution=="8kb", ]
reg_barplot(data, ylim = c(0, 2.9), y_label = "Gene expression\n(log2(FPKM+1))")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig4i.expr_reg.png", width=534/96, height=296/96, dpi=300)


data <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5E_5F/2kbBPSNPs.freq.summary.complete_list.txt", header = T, sep = "\t")
data_p <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure5E_5F/2kbBPSNPs.full.summary.complete_list.txt.proc.test.out", header = T, sep = "\t")
# data$Cell.type = ifelse(grepl("iEC", data$sample), "EC", "VSMC")
data = data[grepl("8kb_loops.tsv.proc$", data$sample), -1]
data_p = data_p[data_p$Resolution=="8kb" & data_p$FDR==0.2, grepl("\\.p$", colnames(data_p))]
data_merged = data.frame(prop = c(as.numeric(data[1,]), as.numeric(data[2,])),
                         p = c(1, as.numeric(data_p[1,]), 1, as.numeric(data_p[2,])),
                         reg = rep(c("Ctrl", "EE", "PP", "TT", "EP", "ET", "PT"), 2), 
                         Cell.type = c(rep("iECs", 7), rep("iVSMCs", 7)))

data_merged %>%
  mutate(
    x.1 = 1,
    x.2 = case_when(
      reg == "EE" ~ 2,
      reg == "EP" ~ 3,
      reg == "ET" ~ 4,
      reg == "PP" ~ 5,
      reg == "PT" ~ 6,
      reg == "TT" ~ 7,
      reg == "Ctrl" ~ 1
    ),
    y_position = case_when(
      reg == "EE" ~ 0.06, # 1.2, # max(m + se, na.rm = TRUE) + 0.1,
      reg == "EP" ~ 0.07, # 1.5, # max(m + se, na.rm = TRUE) + 0.3,
      reg == "ET" ~ 0.08, # 1.8, # max(m + se, na.rm = TRUE) + 0.5,
      reg == "PP" ~ 0.09, # 2.1, # max(m + se, na.rm = TRUE) + 0.7,
      reg == "PT" ~ 0.1, # 2.4, # max(m + se, na.rm = TRUE) + 0.9,
      reg == "TT" ~ 0.11, # 2.7, # max(m + se, na.rm = TRUE) + 1.1,
      reg == "Ctrl" ~ NA
    )
  ) %>%
  ggplot(aes(x = reg, y = prop, fill = reg)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  geom_text(data = . %>% filter(!is.na(x.1) & !is.na(x.2) & p < 0.05),
            aes(x = (x.1 + x.2) / 2, y = y_position, label = paste0("p=", signif(p, 2))),
            size = 3, color = "black", vjust = -0.5) +
  geom_segment(data = . %>% filter(!is.na(x.1) & !is.na(x.2) & p < 0.05),
               aes(x = x.1, xend = x.2, y = y_position, yend = y_position), 
               color = "black", linewidth = 0.5) +
  labs(x = "", y = "Frequency of vascular-associated SNPs", title = "", fill = "Regulatory\nelement\noccurrence") +
  lims(y = c(0, 0.12)) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"), 
        legend.position = "None") +
  facet_wrap(~Cell.type)
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/fig4h.snp_reg.png", width=534/96, height=296/96, dpi=300)







### figure 4d, e
perform_go_enrichment <- function(gene_ranges, microc_df, expr) {
  
  expr_gene_list <- gsub("\\..*", "", expr$Gene_ID)
  
  microc_df <- microc_df[microc_df$FDR<0.2, ]
  
  microc_ranges <- GRanges(
    seqnames = Rle(c(microc_df$BIN1_CHR, microc_df$BIN2_CHROMOSOME)),
    ranges = IRanges(start = c(microc_df$BIN1_START, microc_df$BIN2_START), 
                     end = c(microc_df$BIN1_END, microc_df$BIN2_END))
  )
  # seqlevels(microc_ranges) = chr_levels
  
  overlaps <- findOverlaps(gene_ranges, microc_ranges)
  overlapped_genes <- gene_ranges[queryHits(overlaps)]
  
  gene_ids <- unique(mcols(overlapped_genes)$gene_id)
  gene_ids <- gene_ids[gene_ids %in% expr_gene_list]
  
  if (length(gene_ids) == 0) {
    stop("No overlapping genes found.")
  }
  
  go_results <- enrichGO(
    gene          = gene_ids,
    OrgDb         = org.Hs.eg.db, 
    keyType       = "ENSEMBL", 
    ont           = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  simplified_results <- simplify(
    go_results,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
  
  return(simplified_results)
  
}


expr <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/FPKM.txt.cln", header = T, sep = "\t")
min_samples <- 2
mean_threshold <- 0.1

ec_cols <- colnames(expr)[grepl("EC", colnames(expr))]
expr_ec <- expr[rowSums(expr[ec_cols]!=0) >= min_samples &
                  rowMeans(expr[ec_cols]) >= mean_threshold, ]
vsmc_cols <- colnames(expr)[grepl("VSMC", colnames(expr))]
expr_vsmc <- expr[rowSums(expr[vsmc_cols]!=0) >= min_samples &
                    rowMeans(expr[vsmc_cols]) >= mean_threshold, ]


microc_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iEC3MC_8kb_loops.tsv", header = T, sep = "\t")
microc_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iVSMC7MC_8kb_loops.tsv", header = T, sep = "\t")

result_8kb_ec <- perform_go_enrichment(gene_ranges, microc_ec, expr_ec)
write.table(result_8kb_ec, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.ec.8kb.out", col.names = T, sep = "\t")

result_8kb_vsmc <- perform_go_enrichment(gene_ranges, microc_vsmc, expr_vsmc)
write.table(result_8kb_vsmc, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.vsmc.8kb.out", col.names = T, sep = "\t")


microc_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iEC3MC_4kb_loops.tsv", header = T, sep = "\t")
microc_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iVSMC7MC_4kb_loops.tsv", header = T, sep = "\t")

result_4kb_ec <- perform_go_enrichment(gene_ranges, microc_ec, expr_ec)
write.table(result_4kb_ec, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.ec.4kb.out", col.names = T, sep = "\t")

result_4kb_vsmc <- perform_go_enrichment(gene_ranges, microc_vsmc, expr_vsmc)
write.table(result_4kb_vsmc, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.vsmc.4kb.out", col.names = T, sep = "\t")


microc_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iEC3MC_16kb_loops.tsv", header = T, sep = "\t")
microc_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iVSMC7MC_16kb_loops.tsv", header = T, sep = "\t")

result_16kb_ec <- perform_go_enrichment(gene_ranges, microc_ec, expr_ec)
write.table(result_16kb_ec, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.ec.16kb.out", col.names = T, sep = "\t")

result_16kb_vsmc <- perform_go_enrichment(gene_ranges, microc_vsmc, expr_vsmc)
write.table(result_16kb_vsmc, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.vsmc.16kb.out", col.names = T, sep = "\t")




visualize_enrichment_res <- function(enrich_df, n=20, title){
  
  enrich_df$GeneRatio_numeric <- sapply(enrich_df$GeneRatio, function(x) {
    eval(parse(text = x))
  })
  enrich_df <- enrich_df[1:n, ]
  p <- ggplot(enrich_df, aes(x = GeneRatio_numeric, y = fct_reorder(Description, -1*log10(qvalue)), fill = -1*log10(qvalue))) + 
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
    labs(x = "Gene ratio", y = "", fill = "-log10(q-value)", title = title) + 
    theme(    axis.text = element_text(colour = 'black'),
              legend.position = "left"
    ) + 
    scale_y_discrete(position = "right") + 
    scale_fill_gradient(low = "white", high="red", limits = c(0, -1*min(log10(enrich_df$qvalue))))
  
  print(p)
  
}

result_8kb_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.ec.8kb.out", sep = "\t")
visualize_enrichment_res(result_8kb_ec, title = "iECs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/clusterprofiler.ec.8kb.png", width=746/96, height=345/96, dpi=300)

result_8kb_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.vsmc.8kb.out", sep = "\t")
visualize_enrichment_res(result_8kb_vsmc, title = "iVSMCs")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/clusterprofiler.vsmc.8kb.png", width=746/96, height=345/96, dpi=300)

result_4kb_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.ec.4kb.out", sep = "\t")
visualize_enrichment_res(result_4kb_ec, title = "iECs (4kb)")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/clusterprofiler.ec.4kb.png", width=746/96, height=345/96, dpi=300)

result_4kb_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.vsmc.4kb.out", sep = "\t")
visualize_enrichment_res(result_4kb_vsmc, title = "iVSMCs (4kb)")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/clusterprofiler.vsmc.4kb.png", width=746/96, height=345/96, dpi=300)

result_16kb_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.ec.16kb.out", sep = "\t")
visualize_enrichment_res(result_16kb_ec, title = "iECs (16kb)")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/clusterprofiler.ec.16kb.png", width=746/96, height=345/96, dpi=300)

result_16kb_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/clusterprofiler/clusterprofiler.vsmc.16kb.out", sep = "\t")
visualize_enrichment_res(result_16kb_vsmc, title = "iVSMCs (16kb)")
ggsave("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/figure/clusterprofiler.vsmc.16kb.png", width=746/96, height=345/96, dpi=300)




