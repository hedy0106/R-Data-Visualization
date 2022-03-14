#Set working directory to wherever you have your data
workdir = "C:/Users/a6632/Desktop/tmp"
setwd(workdir)
getwd()

library(ggpubr)

#Basic plots (skip first row with bacteria full names)
df1 <- read.csv("Relative_Abundance.csv", row.names=1, skip = 1, check.names = FALSE)
names(df1)

#Basic plot
for(ii in 13:ncol(df1)){
  res <- cor.test(df1[,ii], df1$GA,method = "spearman",exact=FALSE)
  aa<-ggscatter(df1, x = colnames(df1)[ii] , y = "GA",
                add = "reg.line",   color = "#2e4057", size = 3,ylab="GA",
                conf.int = TRUE,cor.coef = TRUE,
                cor.method = "spearman",
                add.params = list(color = "#00798c"),cor.coeff.args = list(method = "spearman",ggtheme = theme(panel.background = element_rect(fill = "white", colour = "grey50"),strip.background = element_rect(colour = "grey50", fill = "white"))))
  ggsave(aa, filename=paste("GA.",colnames(df1)[ii],".pdf",sep=''), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
}

#Color by groups
for(ii in 13:ncol(df1)){
  aa<-ggscatter(df1, x = colnames(df1)[ii] , y = "GA",
                add = "reg.line",   color = "Condition4",palette = c("#00798c","#edae49"),shape="Condition4", size = 3,conf.int = TRUE,cor.coef = TRUE,
                cor.method = "spearman",
                add.params = list(color ="Condition4"),cor.coeff.args = list(aes(color =Condition4) ,method = "spearman" ),ggtheme = theme(panel.background = element_rect(fill = "white", colour = "grey50")))
  
  ggsave(aa, filename=paste("GA.",colnames(df1)[ii],".Condition4",".pdf",sep=''), dpi=600, width=16,
         height=12, units=c("cm"),colormodel="srgb")
}
