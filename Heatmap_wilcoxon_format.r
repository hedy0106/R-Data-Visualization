#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2) {
    stop("mapping.txt and otu_table.txt files must be supplied by order.", call.=FALSE)
}

design = read.table(args[1], header=T, row.names= 1, sep="\t") 

#design <- design[with(design, order(@@CONDITION@@)), ]
design <- design[order(design$@@CONDITION@@), , drop = FALSE]

otu_table = read.delim(args[2], row.names= 1,  header=T, sep="\t")


idx = rownames(design) %in% colnames(otu_table) 
#sub_design = design[idx,] 
sub_design = subset(design,idx)
count = otu_table[, rownames(sub_design)]


if(as.double(data.frame(colSums(count))[1,1]) > 1) 
{
  norm = t(t(count)/colSums(count,na=T)) 
}else {
  norm = count
}



sim=cor(norm,method="pearson")


library(gplots)
library(RColorBrewer)
pdf(file=paste("heat_otu_cor_samples.pdf", sep=""), height = 12, width = 12)



heatmap.2(sim, Rowv=TRUE, Colv=TRUE, dendrogram='both', trace='none', margins=c(14,14), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none") 
dev.off()



library(edgeR)

d = DGEList(counts=count, group=sub_design$@@CONDITION@@)

d = calcNormFactors(d)


design.mat = model.matrix(~ 0 + d$samples$group)
colnames(design.mat)=levels(as.factor(design$@@CONDITION@@))
d2 = estimateGLMCommonDisp(d, design.mat)
d2 = estimateGLMTagwiseDisp(d2, design.mat)
fit = glmFit(d2, design.mat)


a <- levels(as.factor(design$@@CONDITION@@)) #取genotype有哪幾項
b <- combn(a,2) # 兩兩配對
j <-1
while( j <= length(b[1,]) ){
  
 cat("states,condition_1,vs,condition_2,taxonomy_sig", "\n", file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_sigA.txt", sep=""))
  

 BvsA <- makeContrasts(contrasts = paste(b[1,j],"-",b[2,j],sep = ""), levels=design.mat)

 lrt = glmLRT(fit,contrast=BvsA)

 de_lrt = decideTestsDGE(lrt, adjust.method="fdr", p.value=0.05)
 
  

   pair_group_test = rbind(subset(sub_design, @@CONDITION@@ %in% b[1,j]),subset(sub_design, @@CONDITION@@ %in% b[2,j]))

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
  cat("enriched", "\n", file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_sigA.txt", sep=""), append=TRUE)
  cat(enriched, "\n", file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_sigA.txt", sep=""), append=TRUE)
  cat("depleted", "\n", file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_sigA.txt", sep=""), append=TRUE)
  cat(depleted, "\n", file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_sigA.txt", sep=""), append=TRUE)
  
  write.table(lrt$table, file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_all_pvalue.txt", sep=""))
  

  pair_group = subset(sub_design, @@CONDITION@@ %in% b[,j])

  DE=c(enriched,depleted)
  sub_norm = as.matrix(norm[DE, rownames(pair_group)])

  pdf(file=paste("heat_otu_@@CONDITION@@_",b[1,j],"-",b[2,j],"_wilcoxon_sigA.pdf", sep=""), height = 14, width = 14)

  if(length(sub_norm[,1]) >= 2 && length(sub_norm[1,]) >= 2 ){
    heatmap.2(sub_norm, scale="row", Colv=FALSE, Rowv=FALSE,dendrogram="none", margins=c(14,15), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)), cexCol=1, keysize=1, density.info="histogram", main=NULL, trace="none")
  }
  dev.off()
  j <- j+1
}