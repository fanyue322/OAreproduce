library(dplyr)
library(patchwork)

setwd('./workdir/Step1')

source('function.R')

###############integrate 7 samples via Liger#####################
library(rliger)
library(Seurat)
library(SeuratWrappers)

counts_OA11<-readRDS("./scRNA-seq/OA11/expr.RDS")
counts_OA12<-readRDS("./scRNA-seq/OA12/expr.RDS")
counts_OA2<-readRDS("./scRNA-seq/OA2/expr.RDS")
counts_OA3<-readRDS("./scRNA-seq/OA3/expr.RDS")
counts_C1<-readRDS("./scRNA-seq/C1/expr.RDS")
counts_C2<-readRDS("./scRNA-seq/C2/expr.RDS")
counts_C3<-readRDS("./scRNA-seq/C3/expr.RDS")


counts_list <- list("OA1-1"=as.matrix(counts_OA11), 
					"OA1-2"=as.matrix(counts_OA12), 
					"OA2"=as.matrix(counts_OA2),
					"OA3"=as.matrix(counts_OA3),
					"C1"=as.matrix(counts_C1),
					"C2"=as.matrix(counts_C2),
					"C3"=as.matrix(counts_C3))
norm_method <- "LogNormalize"
ifnb_liger <- scIA(counts_list, method="liger", res_path=res_path, resolution=0.7, num_dim=40, norm_method = norm_method)	

###############TF identification via SCENIC#########################
library(SCENIC)
library(rliger)
library(SeuratWrappers)		
seurat <- ligerToSeurat(ifnb_liger,use.liger.genes = TRUE,by.dataset = FALSE,renormalize = TRUE)
seurat@meta.data$celltype<-seurat@active.ident
seurat<-FindVariableFeatures(seurat,nfeatures=3000)
exprCount<-as.matrix(seurat@assays$RNA@counts[seurat@assays$RNA@var.features,cell_sel])
scenicOptions <- initializeScenic(org="hgnc", dbDir="./cisTargetDB", nCores=1)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 		
genesKept <- geneFiltering(exprCount, scenicOptions)
exprCount_filtered <- exprCount[genesKept, ]
runCorrelation(exprCount_filtered, scenicOptions)
#exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprCount_filtered, scenicOptions)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
scenicOptions <- SCENIC::runSCENIC_3_scoreCells(scenicOptions, exprCount_filtered)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") 

saveRDS(scenicOptions, file="int/scenicOptions.Rds")	
				
				
norm_method <- "LogNormalize"
res <- scIA(counts_list, method="liger", res_path=res_path, resolution=0.8, num_dim=40, norm_method = norm_method)


############Geo-seq data PCA visualization####################
load('./Geo-seq/ST_raw_counts.rda')  #############counts data
pheno=c(rep('SC1',20),rep('SOA1-1',19),rep('SOA1-2',20))
name=colnames(counts)

############Normalization
vst_out <- sctransform::vst(umi = counts, n_genes = 5000, latent_var = c("log_umi"), method="nb", return_gene_attr = TRUE, return_cell_attr = TRUE, verbosity = 1, return_corrected_umi=TRUE)
data=vst_out$y   

name=colnames(data)
layer=as.numeric(substr(name,nchar(name),nchar(name)))  ##numerize AS SZ MZ DZ as 0,1,2,3
pheno=c(rep('SC1',20),rep('SOA1-1',19),rep('SOA1-2',20))
meta.data=data.frame(sample=pheno,layer=layer)
vs=apply(data,1,var)
###PCA###
gene_use=names(sort(vs,decreasing=T))[1:5000]
pca<-prcomp(t(data[gene_use,]),center=TRUE,scale=TRUE)

###layer marker###
library(limma)

#OA#
idx=grep('OA',meta.data$sample)
dat=data[idx,]
    for(i in 0:3)
        {
		    cond=rep(0,ncol(dat))
		    cond[which(meta.info$layer==i)]=1
			cond=factor(cond,levels=c(1,0))
			design <- cbind(Grp1 = 1, Grp1vs0 = cond)
			fit <- lmFit(dat, design)
			fit <- eBayes(fit)
			res <- topTable(fit, n = Inf, sort.by = "none")
			res$log2FC <- as.matrix(fit$coef)[, 2]
			res$pvalues <- as.matrix(fit$p.value)[, 2]
			
			cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value) < 0.05),"\n"))
			res <- ConvertGeneID(res)
			order_res <- res[order(abs(res$pvalues), decreasing=FALSE),]
			res.list.OA <- c(res.list, list(order_res))
		}## end for
#Control#
idx=grep('OA',meta.data$sample)
dat=data[-idx,]
    for(i in 0:3)
        {
		    cond=rep(0,ncol(dat))
		    cond[which(meta.info$layer==i)]=1
			cond=factor(cond,levels=c(1,0))
			design <- cbind(Grp1 = 1, Grp1vs0 = cond)
			fit <- lmFit(dat, design)
			fit <- eBayes(fit)
			res <- topTable(fit, n = Inf, sort.by = "none")
			res$log2FC <- as.matrix(fit$coef)[, 2]
			res$pvalues <- as.matrix(fit$p.value)[, 2]
			
			cat(paste0("Number DE genes: ", sum(p.adjust(res$P.Value) < 0.05),"\n"))
			res <- ConvertGeneID(res)
			order_res <- res[order(abs(res$pvalues), decreasing=FALSE),]
			res.list.OA <- c(res.list, list(order_res))
		}## end for




