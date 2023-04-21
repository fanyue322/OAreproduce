library(dplyr)
library(patchwork)

setwd('./workdir/Step4')

###############Scissor analysis#####################
############sc data##########
library(Seurat)
sc_dataset<-seurat   ###seurat object scRNA-seq dataset

############bulk dataset#########################
data=readxl::read_xlsx('./GSE114007_raw_counts.xlsx')   ###load GSE114007 bulk RNA-seq dataset
data=as.data.frame(data)
data2=readxl::read_xlsx('./GSE114007_raw_counts.xlsx',sheet='OA')
data2=as.data.frame(data2)
bulk_counts=data[,2:19]
bulk_counts=cbind(bulk_counts,data2[,2:21])
rownames(bulk_counts)=data2$symbol

bulk_dataset=bulk_counts
bulk_state=data.frame(barcode=colnames(bulk_counts),Status=c(rep(0,18),rep(1,20)))
bulk_state$Status=factor(bulk_state$Status)
tag=factor(c(0,1))
phenotype=c(rep(0,18),rep(1,20))
names(phenotype)=colnames(bulk_dataset)
tag<-c("Normal","OA")
infos1 <- Scissor(as.matrix(bulk_dataset), sc_dataset, phenotype,tag=tag, alpha = 0.05, 
                 family = "binomial", Save_file = 'Scissor_OA_stage.RData',cutoff=0.5)
				 
##################bulk data2############
load('./PreprocessedCountsMatrix.RData')  ###E-MTAB-4304 data
bulk_counts=countsMatrix.complete
gene=rownames(bulk_counts)
gene=substr(gene,1,15)
library(clusterProfiler)
symbol=bitr(gene,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
uni_symbol=unique(symbol[,2])
idx=match(uni_symbol,symbol[,2])
symbol=symbol[idx,]
idx=match(symbol[,1],gene)
bulk_counts=bulk_counts[idx,]
rownames(bulk_counts)=symbol[,2]	

bulk_dataset=bulk_counts
bulk_state=data.frame(barcode=colnames(bulk_counts),Status=rep(c(0,1),8))
bulk_state$Status=factor(bulk_state$Status)
tag=factor(c(0,1))
phenotype=rep(c(0,1),8)
names(phenotype)=colnames(bulk_dataset)
tag<-c("Normal","OA")
Idents(sc_dataset)=sc_dataset@meta.data$cell.type
infos2 <- Scissor(as.matrix(bulk_dataset), sc_dataset, phenotype,tag=tag, alpha = 0.05, 
                 family = "binomial", Save_file = 'Scissor_OA2_stage.RData',cutoff=0.5)


###############Subclustering OA patients for GSE114007#####################
res_CibersortX=read.csv('./GSE114007_CibersortX.csv')  ## deconvolution results
res_CibersortX_sub=res_CibersortX[19:38,]
Chondrocyte_prop=res_CibersortX_sub[,2:9]
rownames(Chondrocyte_prop)=res_CibersortX_sub$Mixture
d=dist(mat)
hc <- hclust(d,method="ward.D")

#####DE analysis between GroupA and GroupB via DESeq2#####
library(DESeq2)
metadata=data.frame(ind=colnames(OA_counts),group_id=phenotype,age=age,sex=sex) ##incorporate age and sex as covariates
##phenotype: GroupA or GroupB  age: scaled  
metadata$sex=factor(metadata$sex)  #factor sex
metadata$group_id=factor(metadata$group_id)
dds <- DESeqDataSetFromMatrix((OA_counts), 
                              colData = metadata, 
                              design = ~ group_id)
dds <- DESeq(dds)	
metadata$group_id=factor(metadata$group_id)	
contrast <- c("group_id", levels(metadata$group_id)[1], levels(metadata$group_id)[2])	
res <- results(dds, contrast = contrast,alpha = 0.05)							  

############classification##############
library(reticulate)
library(tidyverse)
source_python('./logistic.py')  ##Perform logistic regression with penalty via sklearn

load('./GSE57218_noralized_data.rda')  ##

sidx_oa=sample(idx_oa)   ##permute OA sample index
sidx_h=sample(idx_h)     ##permute Preserved sample index

idx1=c(sidx_oa[1:8],sidx_h[1:8])
idx2=c(sidx_oa[9:16],sidx_h[9:16])
idx3=c(sidx_oa[17:24],sidx_h[17:24])
idx4=c(sidx_oa[25:33],sidx_h[25:33]) ##we randomly seperated the whole dataset into 4 batches, one for training and 3 for testing 

load('./cacao_Pseudobulk_DE.rda') ###DE results from Pseudobulk analysis
celltype=c('EC','HomC','HTC','PreHTC','ProC','RegC')

auc_tmp=c()
roclist=list()
cout=1
genelist=list()
for(i in 1:6)
{
res=cao$test.results$DESeq2.Wald[[i]]$res
res$absfc=abs(res$log2FoldChange)
res=res[order(res$absfc,decreasing=T),]

feature=res$Gene[which(res$absfc>2)]
gene=intersect(bulk_gene,feature)
x=data[gene,idx4]

x=apply(x,1,scale)
x=t(x) #re-normalized

y=c(rep(0,9),rep(1,9))
res<-LogisticClassification(t(x),y)

odds=as.numeric(res$coef_)

x=data[gene,idx1]
x=apply(x,1,scale)
x=t(x) #re-normalized
pred<-res$predict_proba(t(x))
rocdat=data.frame(truth=c(rep(0,8),rep(1,8)),value=pred[,2])
rocobj1<-roc(rocdat$truth,rocdat[,2])

x=data[gene,idx2]
x=apply(x,1,scale)
x=t(x) #re-normalized
pred<-res$predict_proba(t(x))
rocdat=data.frame(truth=c(rep(0,8),rep(1,8)),value=pred[,2])
rocobj2<-roc(rocdat$truth,rocdat[,2])

x=data[gene,idx3]
x=apply(x,1,scale)
x=t(x) #re-normalized
pred<-res$predict_proba(t(x))
rocdat=data.frame(truth=c(rep(0,8),rep(1,8)),value=pred[,2])
rocobj3<-roc(rocdat$truth,rocdat[,2])
res=mean(c(as.numeric(auc(rocobj1)),as.numeric(auc(rocobj2)),as.numeric(auc(rocobj3))))
auc_tmp=c(auc_tmp,res)

rocobj=rocobj1
rocobj$auc=res
rocobj$sensitivities=(rocobj1$sensitivities+rocobj2$sensitivities+rocobj3$sensitivities)/3
rocobj$specificities=(rocobj1$specificities+rocobj2$specificities+rocobj3$specificities)/3
roclist[[cout]]=rocobj

odds=exp(odds)
names(odds)=gene
genelist[[cout]]=odds
cout=cout+1
}


############Osteoarthritis GWAS analysis####################
library(data.table)
library(dplyr)
library(plyr)
pain_str <- paste0("6159-0.", 0:6)
pain_dat <- fread("./pheno_df_ASA2.txt", select = c("eid", pain_str))             ######phenotype data from UKBiobank
knee_pain <- aaply(pain_dat[, -1], 2, function(a) a == 7) %>%
  t() %>% aaply(1, function(a) sum(a, na.rm = T)) 
knee_pain_dat <- data.frame(eid = pain_dat$eid, pain = ifelse(knee_pain != 0, 1, 0)) %>%
  subset(pain == 1, select = eid) %>%
  # select(eid) |>
  data.frame(pain = 1)

no_pain <- aaply(pain_dat[, -1], 2, function(a) a == -7) %>%
  t() %>% aaply(1, function(a) sum(a, na.rm = T)) 
no_pain_dat <- data.frame(eid = pain_dat$eid, 
                          pain = ifelse(no_pain != 0, 1, 0)) %>%
  subset(pain == 1, select = eid) %>%
  # select(eid) |>
  data.frame(pain = 0)

fin_pain_dat <- rbind(knee_pain_dat, no_pain_dat)

###covariates###
cov_dat<-fread("/home/public/data/OA/manuscript/deconvolution/GWAS/cov_ASA2.txt")     #############covariates
idx=match(fin_pain_dat$eid,cov_dat$eid)
cov_dat=cov_dat[idx,]


fam<-read.table('/home/public/data/OA/manuscript/deconvolution/GWAS/geno/merge.fam')
idx=match(fin_pain_dat$eid,fam$V1)
GRM<-read.table('/home/public/data/OA/manuscript/deconvolution/GWAS/GRM.txt.sXX.txt')
GRM=GRM[idx,idx]
rownames(GRM)=fin_pain_dat$eid
colnames(GRM)=fin_pain_dat$eid                                 #####################Related matrix calculated via GEMMA


######################GMMAT analysis#########################
GRM=as.matrix(GRM)           

library(GMMAT)
geno<-read.table('/home/public/data/OA/manuscript/deconvolution/GWAS/geno/merge.bim')
model0<-glmmkin(fixed=disease~Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+
PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=pheno,kins=GRM,id="eid",family=binomial(link="logit"))

glmm.score(model0,infile=geno.file,outfile="glm.score.testoutfile.txt",infile.nrow.skip=5,
infile.ncol.skip=3,infile.ncol.print=1:3,infile.header.print=c("SNP","Allele1","Allele2"))
                       