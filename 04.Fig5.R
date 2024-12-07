library(VennDiagram)
library(ggvenn)
library(gridExtra)
library(ggplot2)

setwd("/xdisk/mliang1/pliu1/BP-SNP-PPG/H202SC23091521-human-EukmRNAseq-NVUS2023090840-X202SC23091521-Z01-F001.batch_correction/CUT_Tag")

###Fig 5B
###CTCF and RAD21 genome occupancies in iEC

Contact = 3893+14661
Occupancy = 12128+14661
OV = 14661

  m1 <- paste("a", 1:OV, sep="")
  m2 <- paste("a", 1:OV, sep="")
  m1 <- append(m1, paste("b",1:(Contact-OV), sep=""))
  m2 <- append(m2, paste("c",1:(Occupancy-OV), sep=""))
  
  mm <- list("CTCF.CUT&Tag"=m1, "RAD21.CUT&Tag"=m2)
  
  pdf(file = "iEC_pLKO_CTCF.iPSC-pLKO-Rad21.CUT_Tag.pdf", width = 6, height=6)
  
  mypal<-c("skyblue", "pink")
  p <- ggvenn(mm,fill_color = mypal,fill_alpha = 0.7,stroke_linetype = "solid",set_name_size = 6.5,text_size = 7.5)
  print(p)
  
  dev.off()

###Fig 5C
###CTCF and RAD21 DEGs in iEC

Contact = 255+111
Occupancy = 4073+111
OV = 111

  m1 <- paste("a", 1:OV, sep="")
  m2 <- paste("a", 1:OV, sep="")
  m1 <- append(m1, paste("b",1:(Contact-OV), sep=""))
  m2 <- append(m2, paste("c",1:(Occupancy-OV), sep=""))
  
  mm <- list("CTCF.Knockdown"=m1, "RAD21.Knockdown"=m2)
  
  pdf(file = "iEC_pLKO_CTCF.iPSC-pLKO-Rad21.Knockdown.pdf", width = 6, height=6)
  
  mypal<-c("skyblue", "pink")
  p <- ggvenn(mm,fill_color = mypal,fill_alpha = 0.7,stroke_linetype = "solid",set_name_size = 6.5,text_size = 7.5)
  print(p)
  
  dev.off()
 
###Fig 5E 
###Odd ratios (DNA loops, chromatin contacts, genome occupancies)

library(gridExtra)
library(ggplot2)

setwd("/xdisk/mliang1/pliu1/Mini-ENCODE-for-human-vascular-cells/Figure_knockdown")
fulldat <- read.table("Several.Fisher.testout.proc", header = T, row.names=NULL,sep="\t")


pdf(file = "Several.Fisher.testout.proc.OR.pdf", width = 8, height=4)

subdat <- fulldat[c(1:4),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(5:8),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(9:13),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:5, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(14:18),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:5, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(19:22),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(23:26),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:4, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))


dev.off()

###Fig 5F 
###Odd ratios (chromatin contacts involving interactions vs ctrl)

setwd("/xdisk/mliang1/pliu1/Mini-ENCODE-for-human-vascular-cells/Figure_knockdown")
fulldat <- read.table("summary.interaction_vs_ctrol.combined.testout.proc", header = T, row.names=NULL,sep="\t")

pdf(file = "summary.interaction_vs_ctrol.combined.testout.proc.OR.pdf", width = 6, height=4)

subdat <- fulldat[c(1:6),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:6, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(7:12),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:6, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(13:18),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:6, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(19:24),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:6, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(25:30),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:6, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

subdat <- fulldat[c(31:36),]
dat <- data.frame(Index=subdat$Index,OR=subdat$odd,LL=subdat$down,UL=subdat$up,label=subdat$label)

ggplot(dat, aes(y = Index, x = OR)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:6, labels = dat$label, trans = "reverse") +
  xlab("Odds Ratio (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x.bottom = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"))

dev.off()
