load('/home/public/data/OA/manuscript/deconvolution/celltypeDE/st_data_vst.rda')
name=colnames(data)
layer=as.numeric(substr(name,nchar(name),nchar(name)))
pheno=c(rep('Normal',20),rep('WHK',19),rep('NK',20))
meta=data.frame(sample=pheno,layer=layer)
vs=apply(data,1,var)
gene_use=names(sort(vs,decreasing=T))[1:2000]
pca<-prcomp(t(data[gene_use,]),center=TRUE,scale=TRUE)
library(dplyr)
pcaImportanceDF <- pca %>%
  broom::tidy(matrix = "pcs") %>%
  dplyr::rename(Individual = percent, Cumulative = cumulative) %>%
  tidyr::gather(key, value, Individual , Cumulative)
pcaDF <- as_data_frame(pca$x,rownames='Sample')
dat=data.frame(pc1=pcaDF$PC1,pc2=pcaDF$PC2,layer=as.factor(layer),sample=pheno)
rownames(dat)=colnames(data) 
library(ggthemes)
p <- ggplot(data = dat, aes(x = pc1, y = pc2, color =layer, shape = pheno)) +
 # geom_hline(yintercept = 0, lty = 2) +
  
#  geom_vline(xintercept = 0, lty = 2) +
  
  geom_point(size=2,alpha = 0.8)+
  xlab(paste0("PC1: ",round(pcaImportanceDF$value[1] * 100),"% variance")) +
  ylab(paste0("PC2: ",round(pcaImportanceDF$value[2] * 100),"% variance")) +
  theme_base()+
theme(plot.margin = margin(1, 1, 1, 1, "cm"),
panel.background = element_blank(),
plot.background = element_blank(),
panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.text = element_text(size = 20,color = "black"),
axis.title.y = element_text(size = 22,face="bold"),
axis.text.x = element_text(size = 22,face="bold"),
axis.title.x=element_blank(),
axis.ticks.x=element_blank(),
#legend.position="bottom",
legend.title=element_text(size = 18,face="bold"),
legend.text=element_text(size = 16),
legend.key = element_rect(colour = "transparent", fill = "white"),
legend.key.size = unit(1.5, 'lines'))+
scale_color_manual(values=c("#F4A582","#80B1D3","#8dd3c7","wheat"))  ###adjust color
#scale_fill_manual(values=c("#F4A582","#80B1D3","#8dd3c7"))+



p

 res.list=list()
 for(i in 0:3){
     norm_counts <- data[,which(meta$sample=='Normal')]
     meta.info <- meta[which(meta$sample=='Normal'),]
     meta.info$condition=0
     meta.info$condition[which(meta.info$layer==i)]=1
     design <- cbind(Grp1 = 1, Grp1vs0 = meta.info$condition)
     fit <- lmFit(norm_counts, design)
     fit <- eBayes(fit)
     res <- topTable(fit, n = Inf, sort.by = "none")
     res$log2FC <- as.matrix(fit$coef)[, 2]
     res$pvalues <- as.matrix(fit$p.value)[, 2]
     
     cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value,method='fdr') < 0.05),"\n"))
     #res <- ConvertGeneID(res)
     order_res <- res[order(abs(res$log2FC), decreasing=TRUE),]
     res.list <- c(res.list, list(order_res))
 }## end for

 names(res.list)=c("L0","L1","L2","L3")
 res.list.normal=res.list
 res.list=list()
 for(i in 0:3){
     norm_counts <- data[,which(meta$sample=='WHK')]
     meta.info <- meta[which(meta$sample=='WHK'),]
     meta.info$condition=0
     meta.info$condition[which(meta.info$layer==i)]=1
     design <- cbind(Grp1 = 1, Grp1vs0 = meta.info$condition)
     fit <- lmFit(norm_counts, design)
     fit <- eBayes(fit)
     res <- topTable(fit, n = Inf, sort.by = "none")
     res$log2FC <- as.matrix(fit$coef)[, 2]
     res$pvalues <- as.matrix(fit$p.value)[, 2]
     
     cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value,method='fdr') < 0.05),"\n"))
     #res <- ConvertGeneID(res)
     order_res <- res[order(abs(res$log2FC), decreasing=TRUE),]
     res.list <- c(res.list, list(order_res))
 }## end for

 names(res.list)=c("L0","L1","L2","L3")
 res.list.whk=res.list
 res.list=list()
 for(i in 0:3){
     norm_counts <- data[,which(meta$sample=='NK')]
     meta.info <- meta[which(meta$sample=='NK'),]
     meta.info$condition=0
     meta.info$condition[which(meta.info$layer==i)]=1
     design <- cbind(Grp1 = 1, Grp1vs0 = meta.info$condition)
     fit <- lmFit(norm_counts, design)
     fit <- eBayes(fit)
     res <- topTable(fit, n = Inf, sort.by = "none")
     res$log2FC <- as.matrix(fit$coef)[, 2]
     res$pvalues <- as.matrix(fit$p.value)[, 2]
     
     cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value,method='fdr') < 0.05),"\n"))
     #res <- ConvertGeneID(res)
     order_res <- res[order(abs(res$log2FC), decreasing=TRUE),]
     res.list <- c(res.list, list(order_res))
 }## end for

 names(res.list)=c("L0","L1","L2","L3")
 res.list.nk=res.list
 
 ######WHOLE#######
 res.list=list()
 for(i in 0:3){
     norm_counts <- data
     meta.info <- meta
     meta.info$condition=0
     meta.info$condition[which(meta.info$layer==i)]=1
     design <- cbind(Grp1 = 1, Grp1vs0 = meta.info$condition)
     fit <- lmFit(norm_counts, design)
     fit <- eBayes(fit)
     res <- topTable(fit, n = Inf, sort.by = "none")
     res$log2FC <- as.matrix(fit$coef)[, 2]
     res$pvalues <- as.matrix(fit$p.value)[, 2]
     
     cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value,method='fdr') < 0.05),"\n"))
     #res <- ConvertGeneID(res)
     order_res <- res[order(abs(res$log2FC), decreasing=TRUE),]
     res.list <- c(res.list, list(order_res))
 }## end for
 
 ####OA#######

 res.list=list()
 for(i in 0:3){
     norm_counts <- data[,which(meta$sample=='NK' | meta$sample=='WHK')]
     meta.info <- meta[which(meta$sample=='NK' | meta$sample=='WHK'),]
     meta.info$condition=0
     meta.info$condition[which(meta.info$layer==i)]=1
     design <- cbind(Grp1 = 1, Grp1vs0 = meta.info$condition)
     fit <- lmFit(norm_counts, design)
     fit <- eBayes(fit)
     res <- topTable(fit, n = Inf, sort.by = "none")
     res$log2FC <- as.matrix(fit$coef)[, 2]
     res$pvalues <- as.matrix(fit$p.value)[, 2]
     
     cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value,method='fdr') < 0.05),"\n"))
     #res <- ConvertGeneID(res)
     order_res <- res[order(abs(res$log2FC), decreasing=TRUE),]
     res.list <- c(res.list, list(order_res))
 }## end for
  names(res.list)=c("L0","L1","L2","L3")
 res.list.OA=res.list
 
library(clusterProfiler)
GO_analysis=function(gene,universe=NULL,go_padj=0.05)
{
id <- bitr(gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
geneList=id[,2]
names(geneList)=as.character(geneList)
if(!is.null(universe))
{
id2 <- bitr(universe,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
universe=id2[,2]
names(universe)=as.character(universe)
}
#go_padj <- 0.05
GO <- enrichGO(OrgDb="org.Hs.eg.db", 
               gene = names(geneList), #a vector of entrez gene id.
               pvalueCutoff = go_padj,
               keyType = "ENTREZID",
               pAdjustMethod = "BH", #选择校正方法one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none
               qvalueCutoff = 0.2,##选择 qvalueCutoff
               ont = "ALL",
               minGSSize = 10, #minimal size of genes annotated by Ontology term for testing
               maxGSSize = 500, #maximal size of genes annotated for testing
               pool=T,
			   universe=universe,
               readable=TRUE)
return(GO)
}	
