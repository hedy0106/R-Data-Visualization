#!/usr/bin/env Rscript
library(ggpubr)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 1) {
    stop("The alpha_10_mapping_sort.data file must be supplied.", call.=FALSE)
}

df1 <- read.table(args[1], header=TRUE,sep= "\t", row.names =1)


a <- levels(df1$@@CONDITION@@) #取@@CONDITION@@有哪幾項
b <- combn(a,2) # 兩兩配對
my_comparisons <- split(b, col(b)) # 矩陣轉list

xx = length(b[1,])

# shannon
if (length(a) > 2) {
#多組比較
shannonbox1 <- ggboxplot(df1, x="@@CONDITION@@", y="shannon", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$shannon)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(shannon)\nStatistic comparisons by anova & t.test") + theme(plot.title = element_text(hjust = 0.5))
# shannonbox1 <- shannonbox1 + font("xy.text", size = 8, color = "gray", face = "bold") #改x,y軸字體大小
shannonbox2 <- ggboxplot(df1, x="@@CONDITION@@", y="shannon", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$shannon)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(shannon)\nStatistic comparisons by kruskal.test & t.test") + theme(plot.title = element_text(hjust = 0.5))
shannonbox3 <- ggboxplot(df1, x="@@CONDITION@@", y="shannon", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$shannon)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(shannon)\nStatistic comparisons by anova & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))
shannonbox4 <- ggboxplot(df1, x="@@CONDITION@@", y="shannon", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$shannon)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(shannon)\nStatistic comparisons by kruskal.test & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(shannonbox1, filename="boxplot.@@CONDITION@@.shannon.anova.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(shannonbox2, filename="boxplot.@@CONDITION@@.shannon.kruskal.test.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(shannonbox3, filename="boxplot.@@CONDITION@@.shannon.anova.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(shannonbox4, filename="boxplot.@@CONDITION@@.shannon.kruskal.test.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")

#chao1

chao1box1 <- ggboxplot(df1, x="@@CONDITION@@", y="chao1", color = "black", fill = "@@CONDITION@@", 
                       add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$chao1)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(chao1)\nStatistic comparisons by anova & t.test") + theme(plot.title = element_text(hjust = 0.5))
chao1box2 <- ggboxplot(df1, x="@@CONDITION@@", y="chao1", color = "black", fill = "@@CONDITION@@", 
                       add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$chao1)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(chao1)\nStatistic comparisons by kruskal.test & t.test") + theme(plot.title = element_text(hjust = 0.5))
chao1box3 <- ggboxplot(df1, x="@@CONDITION@@", y="chao1", color = "black", fill = "@@CONDITION@@", 
                       add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$chao1)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(chao1)\nStatistic comparisons by anova & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))
chao1box4 <- ggboxplot(df1, x="@@CONDITION@@", y="chao1", color = "black", fill = "@@CONDITION@@", 
                       add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$chao1)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(chao1)\nStatistic comparisons by kruskal.test & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(chao1box1, filename="boxplot.@@CONDITION@@.chao1.anova.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(chao1box2, filename="boxplot.@@CONDITION@@.chao1.kruskal.test.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(chao1box3, filename="boxplot.@@CONDITION@@.chao1.anova.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(chao1box4, filename="boxplot.@@CONDITION@@.chao1.kruskal.test.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")

#observed_species

observed_speciesbox1 <- ggboxplot(df1, x="@@CONDITION@@", y="observed_species", color = "black", fill = "@@CONDITION@@", 
                                  add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$observed_species)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(observed_species)\nStatistic comparisons by anova & t.test") + theme(plot.title = element_text(hjust = 0.5))
observed_speciesbox2 <- ggboxplot(df1, x="@@CONDITION@@", y="observed_species", color = "black", fill = "@@CONDITION@@", 
                                  add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$observed_species)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(observed_species)\nStatistic comparisons by kruskal.test & t.test") + theme(plot.title = element_text(hjust = 0.5))
observed_speciesbox3 <- ggboxplot(df1, x="@@CONDITION@@", y="observed_species", color = "black", fill = "@@CONDITION@@", 
                                  add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$observed_species)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(observed_species)\nStatistic comparisons by anova & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))
observed_speciesbox4 <- ggboxplot(df1, x="@@CONDITION@@", y="observed_species", color = "black", fill = "@@CONDITION@@", 
                                  add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$observed_species)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(observed_species)\nStatistic comparisons by kruskal.test & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(observed_speciesbox1, filename="boxplot.@@CONDITION@@.observed_species.anova.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(observed_speciesbox2, filename="boxplot.@@CONDITION@@.observed_species.kruskal.test.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(observed_speciesbox3, filename="boxplot.@@CONDITION@@.observed_species.anova.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(observed_speciesbox4, filename="boxplot.@@CONDITION@@.observed_species.kruskal.test.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")

#PD_whole_tree

PD_whole_treebox1 <- ggboxplot(df1, x="@@CONDITION@@", y="PD_whole_tree", color = "black", fill = "@@CONDITION@@", 
                               add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$PD_whole_tree)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(PD_whole_tree)\nStatistic comparisons by anova & t.test") + theme(plot.title = element_text(hjust = 0.5))
PD_whole_treebox2 <- ggboxplot(df1, x="@@CONDITION@@", y="PD_whole_tree", color = "black", fill = "@@CONDITION@@", 
                               add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$PD_whole_tree)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(PD_whole_tree)\nStatistic comparisons by kruskal.test & t.test") + theme(plot.title = element_text(hjust = 0.5))
PD_whole_treebox3 <- ggboxplot(df1, x="@@CONDITION@@", y="PD_whole_tree", color = "black", fill = "@@CONDITION@@", 
                               add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$PD_whole_tree)*(1+0.1*xx), method = "anova") + labs(title = "Alpha Diversity(PD_whole_tree)\nStatistic comparisons by anova & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))
PD_whole_treebox4 <- ggboxplot(df1, x="@@CONDITION@@", y="PD_whole_tree", color = "black", fill = "@@CONDITION@@", 
                               add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+#不同組間的比較
  stat_compare_means(label.y = max(df1$PD_whole_tree)*(1+0.1*xx), method = "kruskal.test") + labs(title = "Alpha Diversity(PD_whole_tree)\nStatistic comparisons by kruskal.test & wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(PD_whole_treebox1, filename="boxplot.@@CONDITION@@.PD_whole_tree.anova.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(PD_whole_treebox2, filename="boxplot.@@CONDITION@@.PD_whole_tree.kruskal.test.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(PD_whole_treebox3, filename="boxplot.@@CONDITION@@.PD_whole_tree.anova.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(PD_whole_treebox4, filename="boxplot.@@CONDITION@@.PD_whole_tree.kruskal.test.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
} else {
#兩組比較
shannonbox1 <- ggboxplot(df1, x="@@CONDITION@@", y="shannon", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+
						 labs(title = "Alpha Diversity(shannon)\nStatistic comparisons by t.test") + theme(plot.title = element_text(hjust = 0.5))
# shannonbox1 <- shannonbox1 + font("xy.text", size = 8, color = "gray", face = "bold") #改x,y軸字體大小
shannonbox3 <- ggboxplot(df1, x="@@CONDITION@@", y="shannon", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
						 labs(title = "Alpha Diversity(shannon)\nStatistic comparisons by wilcox.test") + theme(plot.title = element_text(hjust = 0.5))
ggsave(shannonbox1, filename="boxplot.@@CONDITION@@.shannon.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(shannonbox3, filename="boxplot.@@CONDITION@@.shannon.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")

#chao1

chao1box1 <- ggboxplot(df1, x="@@CONDITION@@", y="chao1", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+
						 labs(title = "Alpha Diversity(chao1)\nStatistic comparisons by t.test") + theme(plot.title = element_text(hjust = 0.5))
# chao1box1 <- chao1box1 + font("xy.text", size = 8, color = "gray", face = "bold") #改x,y軸字體大小
chao1box3 <- ggboxplot(df1, x="@@CONDITION@@", y="chao1", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
						 labs(title = "Alpha Diversity(chao1)\nStatistic comparisons by wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(chao1box1, filename="boxplot.@@CONDITION@@.chao1.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(chao1box3, filename="boxplot.@@CONDITION@@.chao1.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")

#observed_species

observed_speciesbox1 <- ggboxplot(df1, x="@@CONDITION@@", y="observed_species", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+
						 labs(title = "Alpha Diversity(observed_species)\nStatistic comparisons by t.test") + theme(plot.title = element_text(hjust = 0.5))
# observed_speciesbox1 <- observed_speciesbox1 + font("xy.text", size = 8, color = "gray", face = "bold") #改x,y軸字體大小
observed_speciesbox3 <- ggboxplot(df1, x="@@CONDITION@@", y="observed_species", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
						 labs(title = "Alpha Diversity(observed_species)\nStatistic comparisons by wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(observed_speciesbox1, filename="boxplot.@@CONDITION@@.observed_species.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(observed_speciesbox3, filename="boxplot.@@CONDITION@@.observed_species.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")

#PD_whole_tree

PD_whole_treebox1 <- ggboxplot(df1, x="@@CONDITION@@", y="PD_whole_tree", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "t.test")+
						 labs(title = "Alpha Diversity(PD_whole_tree)\nStatistic comparisons by t.test") + theme(plot.title = element_text(hjust = 0.5))
# PD_whole_treebox1 <- PD_whole_treebox1 + font("xy.text", size = 8, color = "gray", face = "bold") #改x,y軸字體大小
PD_whole_treebox3 <- ggboxplot(df1, x="@@CONDITION@@", y="PD_whole_tree", color = "black", fill = "@@CONDITION@@", 
                         add = "jitter", outlier.shape = NA, shape="@@CONDITION@@") + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
						 labs(title = "Alpha Diversity(PD_whole_tree)\nStatistic comparisons by wilcox.test") + theme(plot.title = element_text(hjust = 0.5))

ggsave(PD_whole_treebox1, filename="boxplot.@@CONDITION@@.PD_whole_tree.t.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
ggsave(PD_whole_treebox3, filename="boxplot.@@CONDITION@@.PD_whole_tree.wilcox.test.pdf", dpi=600, width=16,
       height=12, units=c("cm"),colormodel="srgb")
}