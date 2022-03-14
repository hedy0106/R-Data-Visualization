#Set working directory to wherever you have your data
workdir = "C:/Users/a6632/Desktop/tmp"
setwd(workdir)
getwd()

rm(list=ls())
library(ggpubr)
library(ggplot2)

#Basic plots (colnames with bacteria full names)
df1 <- read.csv("Relative_Abundance.csv", row.names=1, check.names = FALSE)
#Find the index that matches the k__ prefix and split the taxa name and removes the prefix
colnames(df1)[which(grepl("k__",colnames(df1)) == TRUE)[1]:ncol(df1)] = sapply(strsplit(colnames(df1[which(grepl("k__",colnames(df1)) == TRUE)[1]:(ncol(df1))]),"o__"),"[",2)
#remove Condition with NA
df1 <- subset(df1, !is.na(df1$Condition1))

# as.factor change chr into factors
a <- levels(as.factor(df1$Condition14)) #取Condition14有哪幾項
b <- combn(a,2) # 兩兩配對
my_comparisons <- split(b, col(b)) # 矩陣轉list

xx = length(b[1,])

name_list = names(df1)

#Boxplot
for(i in 9:length(name_list)){
  if (length(a) > 2) {
    abundance_box <- ggboxplot(df1, x="Condition14", y=name_list[i], color = "black", fill = "Condition14", 
                                      add = "jitter", outlier.shape = NA, shape="Condition14") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
      stat_compare_means(label.y = max(df1[,i])*(1+0.1*xx), method = "kruskal.test") + labs(title = paste(name_list[i],"\nStatistic comparisons by kruskal.test & wilcox.test")) + theme(plot.title = element_text(hjust = 0.5))
    
    abundance_box 
    
    ggsave(abundance_box, filename=paste("boxplot.Condition14.",name_list[i],".kruskal.test.wilcox.test.pdf"), dpi=600, width=16,
           height=12, units=c("cm"),colormodel="srgb")
    
  } else {
    abundance_box <- ggboxplot(df1, x="Condition14", y=name_list[i], color = "black", fill = "Condition14", 
                                      add = "jitter", outlier.shape = NA, shape="Condition14") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
      labs(title = paste(name_list[i],"\nStatistic comparisons by wilcox.test",sep='')) + theme(plot.title = element_text(hjust = 0.5))
    
    abundance_box 
    
    ggsave(abundance_box, filename=paste("boxplot.Condition14.",name_list[i],".wilcox.test.pdf",sep=''), dpi=600, width=16,
           height=12, units=c("cm"),colormodel="srgb") 
  } 
}