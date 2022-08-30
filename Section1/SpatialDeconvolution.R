################read scRNA-seq data #######################
library(rliger)
library(monocle)
library(Seurat)
library(SeuratWrappers)
load("/home/public/data/OA/manuscript/suppl/res.ia.10XWKvs10XNKvsSKS1vsSKS2vs10XS3vs10XS4vs10XS5.rs0.7.PC40.liger.rename.Rdata")
seurat <- ligerToSeurat(ifnb_liger,use.liger.genes = TRUE,by.dataset = FALSE,renormalize = TRUE)
seurat<-subset(seurat,idents = c("Proliferation-1","Proliferation-2","Homeostasis","Effector","Prehypertrophy","Hypertrophy","Regulator-1","Regulator-2","Fibrocartilage-1","Inflammatory"))## 去除Unknown
seurat@meta.data$cell.type<-seurat@active.ident
seurat@meta.data$cell.type<-factor(seurat@meta.data$cell.type,levels = c("Proliferation-1","Proliferation-2","Homeostasis","Effector","Prehypertrophy","Hypertrophy","Regulator-1","Regulator-2","Fibrocartilage-1","Inflammatory"))
seurat@meta.data$cell.type=as.character(seurat@meta.data$cell.type)
idx=grep('Proliferation-',seurat@meta.data$cell.type)
seurat@meta.data$cell.type[idx]='Proliferation'
idx=grep('Regulator',seurat@meta.data$cell.type)
seurat@meta.data$cell.type[idx]='Regulator'

Idents(seurat)=seurat@meta.data$cell.type
 cluster_markers_all <- Seurat::FindAllMarkers(object = seurat, 
                                               # assay = "SCT",
                                                 slot = "data",
                                                verbose = TRUE, 
                                                only.pos = TRUE,
												group.by='cell.type')
gene=unique(cluster_markers_all$gene)

####OA#########
idx=which(seurat@meta.data$orig.ident=='10XNK' | seurat@meta.data$orig.ident=='10XWK')
####normal#####
idx=which(seurat@meta.data$orig.ident=='10XS4')
idx=setdiff(idx,which(seurat@meta.data$cell.type=='Inflammatory'))
idx=setdiff(idx,which(seurat@meta.data$cell.type=='Fibrocartilage-1'))

counts=seurat@assays$RNA@counts[gene,idx]
colnames(counts)=seurat@meta.data$cell.type[idx]
counts=as.matrix(counts)
counts=data.frame(counts)
counts=cbind(Gene=rownames(counts),counts)
write.table(counts,file=paste0('sc_reference_normal_S4.txt'),quote=F,row.names=F,sep='\t')
