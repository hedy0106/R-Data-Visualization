#!/usr/bin/env Rscript
library(phyloseq)
library(ellipse)
library(ggplot2)
library(plyr)
library(vegan)

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 3) {
    stop("mapping.txt, otu_table_taxonomy.json, and rep_set.tre files must be supplied by order.", call.=FALSE)
}

otutable <- import_biom(BIOMfilename = args[2], treefilename = args[3], parseFunction = parse_taxonomy_greengenes)
mapping <- import_qiime_sample_data(mapfilename = args[1])
phylo_v1 <- merge_phyloseq(otutable, mapping)
phylo_v1 = prune_taxa(taxa_sums(phylo_v1) > 0, phylo_v1)
phylo_v2 = transform_sample_counts(phylo_v1, function(x) x/sum(x))
phylo_selected = subset_samples(phylo_v2)

dist_methods <- list(c("unifrac"),c("wunifrac"),c("bray"))
for( i in dist_methods ){
  method = i

  ordu = ordinate(phylo_selected, "PCoA", method)
  test_selected <- cbind (ordu$vectors[,1:2], sample_data(phylo_selected)[,])
  centroids <- aggregate(cbind(Axis.1,Axis.2)~@@CONDITION@@,test_selected,mean)
  conf.rgn  <- do.call(rbind,lapply(unique(test_selected$@@CONDITION@@),function(t)
    data.frame(@@CONDITION@@=as.character(t),
               ellipse(cov(test_selected[test_selected$@@CONDITION@@==t,1:2]), 
                       centre=as.matrix(centroids[t,2:3]), level=0.95), stringsAsFactors=FALSE)))
  test_selected <- merge(test_selected,centroids,by="@@CONDITION@@",suffixes=c("",".centroid"))

  mapping_table <- read.delim(args[1], row.names=1)
  adonis_table = adonis(phyloseq::distance(phylo_selected, method)~@@CONDITION@@, data=mapping_table, permutations = 10000)
  adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
  adonis_Rvalue = round(as.data.frame(adonis_table$aov.tab[5])[1,1],3)

  p <- ggplot(data=test_selected,(aes(x=Axis.1,y=Axis.2,colour = @@CONDITION@@ )))+
    geom_point(size=3, alpha=0.5)+
    geom_path(data=conf.rgn)+
    ggtitle(paste("PCoA (method=",method,")",sep=''),subtitle=paste("adonis: R=",adonis_Rvalue," P=",format(adonis_pvalue,digits = 3, scientific = TRUE),sep="")) +
    theme(plot.title = element_text(hjust = 0.5, vjust=1))+
    theme_bw()+ theme(panel.grid=element_blank())+
    xlab(paste("PC1 [",round(ordu$values[1,2]*100, digits=1),"%]",sep=""))+
    ylab(paste("PC2 [",round(ordu$values[2,2]*100, digits=1),"%]",sep=""))+
    geom_point(data=centroids, size=2)+
    geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2))
  p
  ggsave(p, filename=paste(method,"_PCoA_ellipse_@@CONDITION@@.pdf", sep=""), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
  
  
  p <- ggplot(data=test_selected,(aes(x=Axis.1,y=Axis.2,colour = @@CONDITION@@ )))+
    geom_point(size=3, alpha=0.5)+
    geom_path(data=conf.rgn)+
    ggtitle(paste("PCoA (method=",method,")",sep=''),subtitle=paste("adonis: R=",adonis_Rvalue," P=",format(adonis_pvalue,digits = 3, scientific = TRUE),sep="")) +
    theme(plot.title = element_text(hjust = 0.5, vjust=1))+
    theme_bw()+ theme(panel.grid=element_blank())+
    xlab(paste("PC1 [",round(ordu$values[1,2]*100, digits=1),"%]",sep=""))+
    ylab(paste("PC2 [",round(ordu$values[2,2]*100, digits=1),"%]",sep=""))+
    geom_point(data=centroids, size=2)+
    geom_segment(aes(x=Axis.1.centroid, y=Axis.2.centroid, xend=Axis.1, yend=Axis.2))+ geom_text(aes(label=X.SampleID),hjust=-0.25, vjust=0, size = 1.5)
  p
  ggsave(p, filename=paste(method,"_PCoA_ellipse_@@CONDITION@@_label.pdf", sep=""), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
  
}


library(phyloseq)
library(ellipse)
library(ggplot2)
library(plyr)
library(vegan)

cat("method condition1~condition2 p-value R", "\n", file="@@CONDITION@@_PCoA_statistics.txt")
otutable <- import_biom(BIOMfilename = args[2], treefilename = args[3], parseFunction = parse_taxonomy_greengenes)
mapping <- import_qiime_sample_data(mapfilename = args[1])
phylo_v1 <- merge_phyloseq(otutable, mapping)
phylo_v1 = prune_taxa(taxa_sums(phylo_v1) > 0, phylo_v1)
phylo_v2 = transform_sample_counts(phylo_v1, function(x) x/sum(x))
phylo_selected = subset_samples(phylo_v2)

dist_methods <- list(c("unifrac"),c("wunifrac"),c("bray"))
for( i in dist_methods ){
  method = i
  ordu = ordinate(phylo_selected, "PCoA", method)
  p <- plot_ordination(phylo_selected, ordu,  color="@@CONDITION@@") +  geom_point(size=3, alpha=0.5)+
    ggtitle(paste("PCoA (method=",method,")",sep=''),subtitle=paste("adonis: R=",adonis_Rvalue," P=",format(adonis_pvalue,digits = 3, scientific = TRUE),sep="")) +
    theme(plot.title = element_text(hjust = 0.5, vjust=1))+
    theme_bw()+ theme(panel.grid=element_blank()) + xlab(paste("PC1 [",round(ordu$values[1,2]*100, digits=1),"%]",sep="")) + ylab(paste("PC2 [",round(ordu$values[2,2]*100, digits=1),"%]",sep=""))
  p
  ggsave(p, filename=paste(method,"_PCoA_@@CONDITION@@.pdf", sep=""), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
  
  p <- plot_ordination(phylo_selected, ordu,  color="@@CONDITION@@") +  geom_point(size=3, alpha=0.5)+
    ggtitle(paste("PCoA (method=",method,")",sep=''),subtitle=paste("adonis: R=",adonis_Rvalue," P=",format(adonis_pvalue,digits = 3, scientific = TRUE),sep="")) +
    theme(plot.title = element_text(hjust = 0.5, vjust=1))+
    theme_bw()+ theme(panel.grid=element_blank()) + xlab(paste("PC1 [",round(ordu$values[1,2]*100, digits=1),"%]",sep="")) + ylab(paste("PC2 [",round(ordu$values[2,2]*100, digits=1),"%]",sep=""))+ geom_text(aes(label=X.SampleID),hjust=-0.25, vjust=0, size = 1.5)
  p
  ggsave(p, filename=paste(method,"_PCoA_@@CONDITION@@_label.pdf", sep=""), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
  
 
  
  mapping_table <- read.delim(args[1], row.names=1)
  a <- levels(mapping_table$@@CONDITION@@) #取genotype有哪幾項
  b <- combn(a,2) # 兩兩配對
  j <-1
  while( j <= length(b[1,]) ){

    
    aa<- as.data.frame(as.matrix(distance(phylo_selected, method)))

    design2_test = subset(mapping_table, @@CONDITION@@ %in% b[,j])
    sub_dis_table_test = aa[rownames(design2_test),rownames(design2_test)]
    sub_dis_table_test <- as.dist(sub_dis_table_test, diag = FALSE, upper = FALSE)
    adonis_table = adonis(sub_dis_table_test~@@CONDITION@@, data=design2_test, permutations = 10000)
    adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
    

    cat(paste(method,paste(b[1,j],"~",b[2,j], sep = ""),format(adonis_pvalue, digits = 3, scientific = TRUE),format(adonis_Rvalue, digits = 3, scientific = TRUE),sep=" "), "\n", file="@@CONDITION@@_PCoA_statistics.txt", append=TRUE)
    j <- j +1
  }

  
  adonis_table = adonis(distance(phylo_selected, method)~@@CONDITION@@, data=mapping_table, permutations = 10000)
  adonis_pvalue = adonis_table$aov.tab$`Pr(>F)`[1]
  #print(paste(method,"All",format(adonis_pvalue, digits = 3, scientific = TRUE)))
  cat(paste(method,"All",format(adonis_pvalue, digits = 3, scientific = TRUE),sep=" "), "\n", file="@@CONDITION@@_PCoA_statistics.txt", append=TRUE)
}

library(phyloseq)
library(ellipse)
library(ggplot2)
library(plyr)
otutable <- import_biom(BIOMfilename = args[2], treefilename = args[3], parseFunction = parse_taxonomy_greengenes)
mapping <- import_qiime_sample_data(mapfilename = args[1])
phylo_v1 <- merge_phyloseq(otutable, mapping)
phylo_v1 = prune_taxa(taxa_sums(phylo_v1) > 0, phylo_v1)
phylo_v2 = transform_sample_counts(phylo_v1, function(x) x/sum(x))

phylo_selected = subset_samples(phylo_v2)

dist_methods <- list(c("unifrac"),c("wunifrac"),c("bray"))
for( i in dist_methods ){
  dist = i
  #dist = "unifrac"
  ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
  plist = llply(as.list(ord_meths), function(i, phylo_selected, dist){
    ordu = ordinate(phylo_selected, method=i, distance=dist)
    plot_ordination(phylo_selected, ordu, color="@@CONDITION@@")
  }, phylo_selected, dist)
  names(plist) <- ord_meths
  pdataframe = ldply(plist, function(x){
    df = x$data[, 1:2]
    colnames(df) = c("Axis_1", "Axis_2")
    return(cbind(df, x$data))
  })
  names(pdataframe)[1] = "method"
  p = ggplot(pdataframe, aes(Axis.1, Axis.2, color=@@CONDITION@@))
  p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=@@CONDITION@@, shape=@@CONDITION@@, fill=@@CONDITION@@))
  p = p + geom_point(size=3, alpha=0.5)
  p = p + facet_wrap(~method, scales="free")
  p = p + ggtitle("Various distance metrics apply on the case") + theme_bw()+ theme(panel.grid=element_blank())
  p

  ggsave(p, filename=paste(dist,"_allplot_@@CONDITION@@.pdf", sep=""), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
  
  
  p = ggplot(pdataframe, aes(Axis.1, Axis.2, color=@@CONDITION@@))
  p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=@@CONDITION@@, shape=@@CONDITION@@, fill=@@CONDITION@@))
  p = p + geom_point(size=3, alpha=0.5)
  p = p + facet_wrap(~method, scales="free")
  p = p + ggtitle("Various distance metrics apply on the case") + theme_bw()+ theme(panel.grid=element_blank())+ geom_text(aes(label=X.SampleID),hjust=-0.25, vjust=0, size = 1.5)
  p

  ggsave(p, filename=paste(dist,"_allplot_@@CONDITION@@_label.pdf", sep=""), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
  
  
}



#==============
