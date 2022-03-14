#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 3) {
    stop("mapping.txt, otu_table.txt, and ?_otu_taxonomy_Manhattan_used.txt files must be supplied by order.", call.=FALSE)
}

design = read.table(args[1], header=T, row.names= 1, sep="\t") 

#design <- design[with(design, order(@@CONDITION@@)), ]
design <- design[order(design$@@CONDITION@@), , drop = FALSE]

otu_table = read.delim(args[2], row.names= 1,  header=T, sep="\t")

idx = rownames(design) %in% colnames(otu_table) 
#sub_design = design[idx,] #如要選擇其中幾個樣本畫的畫，請挑出所要樣本組成otu.table
sub_design = subset(design,idx)
count = otu_table[, rownames(sub_design)]


if(as.double(data.frame(colSums(count))[1,1]) > 1) 
{
  norm = t(t(count)/colSums(count,na=T)) # normalization to total 1
}else {
  norm = count
}


# 计算所有样品间相关系数
sim=cor(norm,method="pearson")

# 使用热图可视化，并保存为8x8英寸的PDF
library(gplots)
library(RColorBrewer)


library(edgeR)
# create DGE list
d = DGEList(counts=count, group=sub_design$@@CONDITION@@)
#d = DGEList(counts=count, group=sub_design$genotype)
d = calcNormFactors(d)

# 生成实验设计矩阵
design.mat = model.matrix(~ 0 + d$samples$group)
colnames(design.mat)=levels(as.factor(design$@@CONDITION@@))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)


a <- levels(as.factor(design$@@CONDITION@@)) #取genotype有哪幾項
b <- combn(a,2) # 兩兩配對
j <-1
while( j <= length(b[1,]) ){
  
  # cat("states,condition_1,vs,condition_2,taxonomy_sig", "\n", file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_sigA.txt", sep=""))
  

  BvsA <- makeContrasts(contrasts = paste(b[1,j],"-",b[2,j],sep = ""), levels=design.mat)
  # 组间比较,统计Fold change, Pvalue, to perform likelihood ratio tests:
  lrt = glmLRT(fit,contrast=BvsA)

  de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
  

  pair_group_test = rbind(subset(sub_design, @@CONDITION@@ %in% b[1,j]),subset(sub_design, @@CONDITION@@ %in% b[2,j]))
  #length(subset(sub_design, @@CONDITION@@ %in% b[1,j]))
  #length(subset(sub_design, @@CONDITION@@ %in% b[2,j]))
  aa = as.matrix(norm[, rownames(pair_group_test)])
  
  print(paste(b[1,j],"-",b[2,j],"processing", sep=" "))
  
  for (i in 1:nrow(norm)){
    #print(paste(b[1,j],"-",b[2,j],i, sep=" "))
    g1 <- aa[i,1:length(subset(sub_design, @@CONDITION@@ %in% b[1,j])[,1])]
    g2 <- aa[i,(length(subset(sub_design, @@CONDITION@@ %in% b[1,j])[,1])+1):(length(subset(sub_design, @@CONDITION@@ %in% b[1,j])[,1])+length(subset(sub_design, @@CONDITION@@ %in% b[2,j])[,1]))]
    
    g_12 <- wilcox.test(as.numeric(g1),as.numeric(g2))
    
    lrt$table$wilcoxon_p.value[i] = g_12$p.value
    
    if (g_12$p.value != "NaN")
    {
      
      if(g_12$p.value < 0.05)
        
      {
        if(lrt$table$logFC[i] > 0)
        {
          lrt$table$wilcox_sig[i] = 1
        }
        if(lrt$table$logFC[i] < 0)
        {
          lrt$table$wilcox_sig[i] = -1
        }     
        
      }
      else
      {
        lrt$table$wilcox_sig[i] = 0
      }
      
    }
    
    
    else
    {
      lrt$table$wilcox_sig[i] = 0
    }
    
  }
  
  
  #### # 导出计算结果
  x=lrt$table
  x$sig=de_lrt
  enriched = row.names(subset(x,wilcox_sig==1))
  depleted = row.names(subset(x,wilcox_sig==-1))

  taxonomy = read.delim(args[3], row.names= 1,header=F, sep=";")
  #colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species","evalue")
  colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
  
  # 标记差异OTU类型
  x$level = as.factor(ifelse(x$wilcox_sig==1, "enriched",ifelse(x$wilcox_sig==-1, "depleted","nosig")))
  x$otu = rownames(x)
  # 转换Pvalue为负对数
  x$neglogp = -log10(x$wilcoxon_p.value)
  
  # Taxonomy排序，并筛选OTU表中存在的
  library(dplyr)
  taxonomy$id=rownames(taxonomy)
  taxonomy = arrange(taxonomy, kingdom, phylum, class, order, family, genus, species)
  rownames(taxonomy) = taxonomy$id
  idx = rownames(taxonomy) %in% x$otu
  tax = taxonomy[idx, ] # subset taxonomy from used OTU
  
  # 手动筛选显著的组
  x = x[rownames(tax), ] # reorder according to tax
  x$tax = gsub("p__","",tax$phylum,perl=TRUE) 
  top_phylum=c(" Bacteroidetes"," Firmicutes"," Planctomycetes"," Proteobacteria"," Verrucomicrobia")
  x[!(x$tax %in% top_phylum),]$tax = "Low Abundance" # no level can get value
  # 设置各类的level对应顺序
  x$otu = factor(x$otu, levels=x$otu)   # set x order
  x$level = factor(x$level, levels=c("enriched","depleted","nosig"))
  levels(x$tax)=c(top_phylum,"Low Abundance")

  library(ggplot2)

  if (max(x$logFC)>4){x[x$logFC>4,]$logFC = 4} # norm x axis
  if (min(x$logFC)< -4){x[x$logFC< -4,]$logFC = -4} # norm x axis
  x$level = as.factor(ifelse(x$wilcox_sig==1, "enriched",ifelse(x$wilcox_sig==-1, "depleted","nosig")))
  

  p = ggplot(x, aes(x=logFC, y=neglogp, color=level)) + geom_point()  + 
    scale_colour_manual(values=c("green","red","grey"))+ xlim(-4, 4)+
    labs(x="log2(fold change)",y="Negative log10 p-value", title=paste(b[1,j],"vs",b[2,j], sep=" "))+ theme_bw()+ theme(panel.grid=element_blank())
  p
  ggsave(file=paste("vol_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon.pdf", sep=""), p, width = 8, height = 5)

  
  
  p = ggplot(x, aes(x=logFC, y=neglogp, color=level, size=logCPM, shape=tax)) + geom_point()  + 
    scale_colour_manual(values=c("red","green","grey"))+ xlim(-4, 4)+
    labs(x="log2(fold change)",y="Negative log10 p-value", title=paste(b[1,j],"vs",b[2,j], sep=" "))+ theme_bw()+ theme(panel.grid=element_blank())
  ggsave(file=paste("vol_otu_tax_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon.pdf", sep=""), p, width = 8, height = 5)
  
  j <- j+1
  
}