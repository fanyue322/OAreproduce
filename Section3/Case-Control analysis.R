library(dplyr)
library(patchwork)

setwd('./workdir/Step3')

###############Pseudobulk Case-control analysis#####################

#############Initialize######################
library(cacoa)
sample.groups=unique(seurat@meta.data$orig.ident)
sample.groups=as.character(sample.groups)
names(sample.groups)=sample.groups
sample.groups[c(1,2,3,4)]='OA'
sample.groups[c(5,6,7)]='Control'
cell.groups=seurat@meta.data$cell.type
names(cell.groups)=colnames(seurat)
sample.per.cell=seurat@meta.data$orig.ident
names(sample.per.cell)=colnames(seurat)
ref.level='Control'
target.level='OA'
cm=seurat@assays$RNA@counts
embedding=seurat@reductions$tsne@cell.embeddings
rownames(embedding)=colnames(cm)

#############remove unshared cell types##############
idx=which(seurat@meta.data$cell.type=='Proliferation-2' | seurat@meta.data$cell.type== 'Regulator-2' |
seurat@meta.data$cell.type== 'Inflammatory')

cao<-Cacoa$new(as.matrix(cm[,-idx]),sample.groups=sample.groups, cell.groups=cell.groups[-idx], sample.per.cell=sample.per.cell[-idx], 
                 ref.level=ref.level, target.level=target.level, embedding=embedding[-idx,])

#######Cell composition############
cao$plotCellGroupSizes(show.significance=TRUE, legend.position=c(1, 1)) 	
			 
############DE analysis##############################
example.celltypes <- names(sort(table(cao$cell.groups),decreasing=T))
example.cell.type <- example.celltypes[1]				 
				 
possible.tests <- c('DESeq2.Wald')		 
possible.resampling <- c('loo')

possible.tests.loo = paste(possible.tests, possible.resampling, sep = '.')		 
		
f.changes = FALSE
for(test in possible.tests){
  # message(test)
  for(resampling in possible.resampling){
    name <- paste(test, resampling, sep='.')
    if(!(name %in% names(cao$test.results))){
      print(name)
      cao$estimateDEPerCellType(test = test, resampling.method = resampling,
                              name = name)
      f.changes = TRUE
    }
  }
}						 

################Smart-seq2 data analysis#####################
load('./seurat_Smartseq2.rda')   ##seurat object  set stage as a meta information S0-S4
###############merge stages#################
###ProC###
idx=which(seurat@meta.data$celltype=='ProC')
counts=seurat@assays$RNA@counts[,idx]
time=seurat@meta.data$stage[idx]
###PreHTC###
idx=which(seurat@meta.data$celltype=='preHTC')
counts=seurat@assays$RNA@counts[,idx]
time=seurat@meta.data$stage[idx]
time[which(time==1)]=0
time[which(time==2)]=0
time[which(time==3)]=1
time[which(time==4)]=2
###HTC###
idx=which(seurat@meta.data$celltype=='HTC')
counts=seurat@assays$RNA@counts[,idx]
time=seurat@meta.data$stage[idx]
time[which(time==3)]=2
time[which(time==4)]=3

###HomC###
idx=which(seurat@meta.data$celltype=='HomC')
counts=seurat@assays$RNA@counts[,idx]
time=seurat@meta.data$stage[idx]

###EC###
idx=which(seurat@meta.data$celltype=='EC')
counts=seurat@assays$RNA@counts[,idx]
time=seurat@meta.data$stage[idx]
time[which(time==4)]=3
###RegC###
idx=which(seurat@meta.data$celltype=='RegC')
counts=seurat@assays$RNA@counts[,idx]
time=seurat@meta.data$stage[idx]
time[which(time==4)]=2
time[which(time==3)]=2

###########Identifying stage-wise DE genes via tradeSeq #####
cellWeights=rep(1,ncol(counts))
sce<-fitGAM(counts=as.matrix(counts),pseudotime=time,cellWeights=cellWeights,nknots=3,verbose=FALSE)
assoRes<-associationTest(sce)
##########gene expression patterns clustering#####
ibrary(clusterExperiment)
nPointsClus<-3
idx=which(assoRes$p<0.05)
clusPat<-clusterExpressionPatterns(sce,nPoints=nPointsClus,rownames(assoRes)[idx],reduceMethod="none")
clusterLabels <- primaryCluster(clusPat$rsec)
gene_sig=rownames(assoRes)[idx]

###############Geo-seq DE analysis  case vs control##############
load('/home/public/data/OA/manuscript/deconvolution/celltypeDE/Geo-seq.vstransformed.rda')
### data: vst transformed GEO-seq data  ###
name=colnames(data)
layer=as.numeric(substr(name,nchar(name),nchar(name))) ##0 AS 1 SZ 2 MZ 3 DZ
pheno=c(rep('OA',20),rep('OA',19),rep('Control',20))
meta=data.frame(sample=pheno,layer=layer)
 res.list=list()
 for(i in 0:3){
     norm_counts <- data[,which(meta$layer==i)]
     meta.info <- meta[which(meta$layer==i),]
     meta.info$condition=0
     meta.info$condition[which(meta.info$pheno=='OA')]=1
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
 
 
###Function for GO and KEGG analysis#####################
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

KEGG_analysis=function(gene,universe=NULL,go_padj=0.05,qvalue=0.2)
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
  GO <- enrichKEGG(gene=geneList,
                   organism='hsa', ###mouse mmu human hsa rat rno
                   keyType='kegg',
                   pvalueCutoff=go_padj,
                   pAdjustMethod='BH',
                   universe=universe,
                   minGSSize=10,
                   maxGSSize=500,
                   qvalueCutoff=qvalue)
  return(GO)
}	