#Set working directory to wherever you have your data
workdir = "C:/Users/a6632/Desktop"
setwd(workdir)
getwd()

# Q30 distributions
QC_out <- read.delim("QC_out.txt", check.names = FALSE)
colnames(QC_out) = c("SampleID","Reads","Avg_Score","Q20","Q28","Q30")
Q30_Avg = paste("Percentage of reads reach Q30 =",round(mean(QC_out$Q30),2),"%")
head(QC_out)

library(ggpubr)

# Basic bar plots with label (Q30 plot)
p = ggbarplot(QC_out, x = "SampleID", y = "Q30", 
              title = "Q30 distribution", font.title = ("bold.italic"),
              subtitle = Q30_Avg, font.subtitle = c(size = 10, "bold"),
              fill="#008fa6", x.text.angle = 90, font.x = c(size=10, "bold"),
              font.y = c(size=10, "bold"), ylab = "Reads percentage", ylim = c(0,100))
# Add horizontal line at y = 75; Change line type and color; Change line size; Change the appearance of titles and labels
Q30_plot = p + geom_hline(yintercept=75, linetype="dashed", color = "#edae49", size=0.8)+
  font("x.text", size = 8)+
  font("y.text", size = 8)

#ggsave("Q30_plot.pdf", width = 10, height = 6)

# Reads information
read_info <- read.delim("read_info.txt", check.names = FALSE)
colnames(read_info) = c("SampleID","Raw_Reads","Assembled_Reads","Effective_Reads","Taxon_Reads","OTUs")
Raw_Avg = paste("Average raw reads =",round(mean(read_info$Raw_Reads),0))
Eff_Avg = paste("Average effective reads =",round(mean(read_info$Effective_Reads),0))
head(read_info)

# Basic bar plots with label (Raw reads plot)
options(scipen = 999)
p = ggbarplot(read_info, x = "SampleID", y = "Raw_Reads", 
              title = "Raw reads distribution", font.title = ("bold.italic"),
              subtitle = Raw_Avg, font.subtitle = c(size = 10, "bold"),
              fill="#008fa6", x.text.angle = 90, font.x = c(size=10, "bold"),
              font.y = c(size=10, "bold"), ylab = "Raw Reads")
# Add horizontal line at y = 75; Change line type and color; Change line size; Change the appearance of titles and labels
RawReads_plot = p + geom_hline(yintercept=50000, linetype="dashed", color = "#edae49", size=0.8)+
  font("x.text", size = 8)+
  font("y.text", size = 8)

#ggsave("RawReads_plot.pdf", width = 10, height = 6)

# Basic bar plots with label (Effective reads plot)
options(scipen = 999)
p = ggbarplot(read_info, x = "SampleID", y = "Effective_Reads", 
              title = "Effective reads distribution", font.title = ("bold.italic"),
              subtitle = Eff_Avg, font.subtitle = c(size = 10, "bold"),
              fill="#008fa6", x.text.angle = 90, font.x = c(size=10, "bold"),
              font.y = c(size=10, "bold"), ylab = "Effective Reads")
# Add horizontal line at y = 75; Change line type and color; Change line size; Change the appearance of titles and labels
EffectiveReads_plot = p + geom_hline(yintercept=20000, linetype="dashed", color = "#edae49", size=0.8)+
  font("x.text", size = 8)+
  font("y.text", size = 8)

#ggsave("EffectiveReads_plot.pdf", width = 10, height = 6)

# for one plot per page
# ggexport(muti_plot, filename="QC_information.pdf", width = 25, height = 12) ?€n = 255?€?
muti_plot = ggarrange(Q30_plot,RawReads_plot,EffectiveReads_plot, nrow=1, ncol=1)
ggexport(muti_plot, filename="QC_information.pdf", width = 10, height = 6)