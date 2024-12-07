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
### extract snp information
gwas_folder <- "xdisk/mliang1/qqiu/data/HT-GWAS/"

gwas_list <- c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv", 
               "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv", 
               "gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
               "gwas_catalog_peripheral_arterial_disease.tsv", "gwas_catalog_aortic_aneurysm.tsv")
outfile <- "gwas_catalog_target_trait.snp.txt"

gwas_list <- c("control_traits/gwas_catalog_alzheimer_disease.tsv", "control_traits/gwas_catalog_asthma.tsv",
               "control_traits/gwas_catalog_birth_weight.tsv", "control_traits/gwas_catalog_bone_density.tsv", 
               "control_traits/gwas_catalog_breast_cancer.tsv", "control_traits/gwas_catalog_colorectal_cancer.tsv",
               "control_traits/gwas_catalog_copd.tsv", "control_traits/gwas_catalog_multiple_sclerosis.tsv",
               "control_traits/gwas_catalog_parkinson_disease.tsv", "control_traits/gwas_catalog_psoriasis.tsv",
               "control_traits/gwas_catalog_rheumatoid_arthritis.tsv")
outfile <- "gwas_catalog_control_trait.snp.txt"

SNPS <- c()
for(gwas_file in gwas_list){
  
  output <- gsub("tsv", "snp.txt", gwas_file)
  trait <- gsub("gwas_catalog_|\\.tsv", "", basename(gwas_file))
  
  dat <- read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
  dat <- dat[dat$P.VALUE<5e-8, ]
  
  snp_pos <- data.frame(SNP=dat$SNPS, POS = paste0("chr", dat$CHR_ID, "_", dat$CHR_POS),
                        P.VALUE=dat$P.VALUE, BETA=dat$OR.or.BETA, TRAIT=trait, TYPE=dat$CONTEXT)
  
  snp_pos <- snp_pos[(grepl("^rs", snp_pos$SNP)) & !(grepl("[,|;|x]", snp_pos$SNP) | grepl("NA", snp_pos$POS)), ]
  
  snp_pos_filtered <- snp_pos %>%
    group_by(SNP) %>%
    filter(P.VALUE == min(P.VALUE)) %>%
    ungroup()
  
  SNPS <- rbind(SNPS, snp_pos_filtered)
  
}

write.table(SNPS, outfile, col.names = T, row.names = F, quote = F, sep = '\t')





################################################################################
### load reference files
canonical_chromosomes <- c(1:22, "X", "Y")
chr_levels = paste0("chr", canonical_chromosomes)

ref <- read.table("/xdisk/mliang1/qqiu/project/others/huves/data/Homo_sapiens.GRCh38.112.gene.gtf", sep = '\t', header=F)
ref <- ref[ref$V1 %in% canonical_chromosomes, ]
ref$gene_id <- str_extract(ref$V9, "gene_id [^;]+") %>%
  gsub("gene_id ", "", .)
ref$gene_name <- str_extract(ref$V9, "gene_name [^;]+") %>%
  gsub("gene_name ", "", .)
ref$gene_type <- str_extract(ref$V9, "gene_biotype [^;]+") %>%
  gsub("gene_biotype ", "", .)

gene_ranges <- GRanges(
  seqnames = Rle(paste0("chr", ref$V1)),
  ranges = IRanges(start = ref$V4, end = ref$V5),
  strand = Rle(ref$V7),
  gene_id = ref$gene_id,
  gene_name = ref$gene_name,
  gene_type = ref$gene_type,
  gene_chr = paste0("chr", ref$V1),
  gene_start = ref$V4,
  gene_end = ref$V5,
  gene_strand = ref$V7
)

promoter_ranges <- promoters(gene_ranges, upstream = 2000, downstream = 200)
seqlevels(promoter_ranges) = chr_levels







################################################################################
### expr data process
expr <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1A_1B/FPKM.txt.cln", header = T, sep = "\t")

min_samples <- 2
mean_threshold <- 0.1

ec_cols <- colnames(expr)[grepl("EC", colnames(expr))]
expr_ec <- expr[rowSums(expr[ec_cols]!=0) >= min_samples &
                  rowMeans(expr[ec_cols]) >= mean_threshold, ]
expr_ec_gr <- GRanges(
  seqnames = Rle(paste0("chr", expr_ec$Reference)),
  ranges = IRanges(start = expr_ec$Start, end = expr_ec$End),
  strand = Rle(expr_ec$Strand)
)

vsmc_cols <- colnames(expr)[grepl("VSMC", colnames(expr))]
expr_vsmc <- expr[rowSums(expr[vsmc_cols]!=0) >= min_samples &
                    rowMeans(expr[vsmc_cols]) >= mean_threshold, ]
expr_vsmc_gr <- GRanges(
  seqnames = Rle(paste0("chr", expr_vsmc$Reference)),
  ranges = IRanges(start = expr_vsmc$Start, end = expr_vsmc$End),
  strand = Rle(expr_vsmc$Strand)
)
seqlevels(expr_vsmc_gr) = chr_levels


### rrbs data process
rrbs <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1C/iPSC.RRBS.dat.mincov5.5mC.cln.combat1.proc3.srtbycv", header = T, sep = "\t")
rrbs <- rrbs[rrbs$chr!="M", ]

rrbs_gr <- GRanges(
  seqnames = Rle(paste0("chr", rrbs$chr)),
  ranges = IRanges(start = rrbs$pos, end = rrbs$pos)
)

seqlevels(rrbs_gr) = chr_levels
# common_chromosomes <- intersect(seqlevels(promoter_ranges), seqlevels(rrbs_gr))
# promoter_chr_filter_ranges <- keepSeqlevels(promoter_ranges, common_chromosomes, pruning.mode="coarse")
# rrbs_chr_filter_gr <- keepSeqlevels(rrbs_gr, common_chromosomes, pruning.mode="coarse")
counts <- countOverlaps(promoter_ranges, rrbs_gr)
mcols(promoter_ranges)$RRBS_site_count <- counts
methylated_promoters_gr <- promoter_ranges[counts > 0]


### atac data process
atac_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1F_1G/EC.narrowPeak", sep = "\t")
atac_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1F_1G/VSMC.narrowPeak", sep = "\t")

atac_ec_gr <- GRanges(
  seqnames = Rle(atac_ec$V1),
  ranges = IRanges(start = atac_ec$V2, end = atac_ec$V3)
)
seqlevels(atac_ec_gr) = chr_levels

atac_vsmc_gr <- GRanges(
  seqnames = Rle(atac_vsmc$V1),
  ranges = IRanges(start = atac_vsmc$V2, end = atac_vsmc$V3)
)
seqlevels(atac_vsmc_gr) = chr_levels



### micro-c data process
microc_ec <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iEC3MC_8kb_loops.tsv", header = T, sep = "\t")
microc_vsmc <- read.table("/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/Figure1H/iVSMC7MC_8kb_loops.tsv", header = T, sep = "\t")

# microc_ec <- microc_ec[microc_ec$FDR<0.05, ]
# microc_vsmc <- microc_vsmc[microc_vsmc$FDR<0.05, ]

microc_ec <- microc_ec[microc_ec$FDR<0.2, ]
microc_vsmc <- microc_vsmc[microc_vsmc$FDR<0.2, ]

microc_melt_ec_gr <- GRanges(
  seqnames = Rle(c(microc_ec$BIN1_CHR, microc_ec$BIN2_CHROMOSOME)),
  ranges = IRanges(start = c(microc_ec$BIN1_START, microc_ec$BIN2_START), 
                   end = c(microc_ec$BIN1_END, microc_ec$BIN2_END))
)
seqlevels(microc_melt_ec_gr) = chr_levels

microc_melt_vsmc_gr <- GRanges(
  seqnames = Rle(c(microc_vsmc$BIN1_CHR, microc_vsmc$BIN2_CHROMOSOME)),
  ranges = IRanges(start = c(microc_vsmc$BIN1_START, microc_vsmc$BIN2_START), 
                   end = c(microc_vsmc$BIN1_END, microc_vsmc$BIN2_END))
)
seqlevels(microc_melt_vsmc_gr) = chr_levels


microc_loops_ec <- microc_ec %>%
  filter(BIN1_CHR==BIN2_CHROMOSOME) %>%
  rowwise() %>%
  mutate(
    BIN_END = max(c_across(c(BIN1_START, BIN1_END, BIN2_START, BIN2_END))), 
    BIN_START = min(c_across(c(BIN1_START, BIN1_END, BIN2_START, BIN2_END)))  
  ) %>% ungroup()

microc_loops_ec_gr <- GRanges(
  seqnames = Rle(microc_loops_ec$BIN1_CHR),
  ranges = IRanges(start = microc_loops_ec$BIN_START, 
                   end = microc_loops_ec$BIN_END)
)
seqlevels(microc_loops_ec_gr) = chr_levels


microc_loops_vsmc <- microc_vsmc %>%
  filter(BIN1_CHR==BIN2_CHROMOSOME) %>%
  rowwise() %>%
  mutate(
    BIN_END = max(c_across(c(BIN1_START, BIN1_END, BIN2_START, BIN2_END))), 
    BIN_START = min(c_across(c(BIN1_START, BIN1_END, BIN2_START, BIN2_END)))  
  ) %>% ungroup()

microc_loops_vsmc_gr <- GRanges(
  seqnames = Rle(microc_loops_vsmc$BIN1_CHR),
  ranges = IRanges(start = microc_loops_vsmc$BIN_START, 
                   end = microc_loops_vsmc$BIN_END)
)
seqlevels(microc_loops_vsmc_gr) = chr_levels



gr_list <- list(expression_ec=expr_ec_gr, expression_vsmc=expr_vsmc_gr,
                rrbs_ec=methylated_promoters_gr, rrbs_vsmc=methylated_promoters_gr,
                atac_ec=atac_ec_gr, atac_vsmc=atac_vsmc_gr,
                microccontacts_ec=microc_melt_ec_gr, microccontacts_vsmc=microc_melt_vsmc_gr,
                microcloops_ec=microc_loops_ec_gr, microcloops_vsmc=microc_loops_vsmc_gr)

saveRDS(gr_list, "/xdisk/mliang1/qqiu/project/others/iPSC-induced_EC_and_VSMC/data/gr_list.rds")



