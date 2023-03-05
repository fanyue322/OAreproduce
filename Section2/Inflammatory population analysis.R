setwd('./workdir/Step2')

################trajectory analysis#######################
library(rliger)
library(monocle)
library(Seurat)
library(SeuratWrappers)
load("./workdir/Step1/res.ia.10XWKvs10XNKvsSKS1vsSKS2vs10XS3vs10XS4vs10XS5.rs0.7.PC40.liger.rename.Rdata")
seurat <- ligerToSeurat(ifnb_liger,use.liger.genes = TRUE,by.dataset = FALSE,renormalize = TRUE)
seurat<-subset(seurat,idents = c("Proliferation-1","Proliferation-2","Homeostasis","Effector","Prehypertrophy","Hypertrophy","Regulator-1","Regulator-2","Fibrocartilage-1","Inflammatory"))## 去除Unknown
seurat@meta.data$cell.type<-seurat@active.ident
seurat@meta.data$cell.type<-factor(seurat@meta.data$cell.type,levels = c("Proliferation-1","Proliferation-2","Homeostasis","Effector","Prehypertrophy","Hypertrophy","Regulator-1","Regulator-2","Fibrocartilage-1","Inflammatory"))

theme_set(theme_bw())
pd <- data.frame(cell_id = colnames(seurat), 
                 cell_type = seurat@meta.data$cell.type,
				# Cluster = seurat@meta.data$seurat_clusters,
                 row.names = colnames(seurat))
pd <- new("AnnotatedDataFrame", data = pd)
fd <- data.frame(gene_id = rownames(seurat), 
                 gene_short_name = rownames(seurat),
                 row.names = rownames(seurat))			 
fd <- new("AnnotatedDataFrame", data = fd)				 
cds <- newCellDataSet(seurat@assays$RNA@counts, phenoData = pd, featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, clustering_genes$gene_id)
cds <- reduceDimension(cds, num_dim = 40, reduction_method = 'tSNE')
diff_genes <- differentialGeneTest(cds, fullModelFormulaStr = "~ cell_type",
                                   cores = 10)
ordering_genes <- row.names(subset(diff_genes, qval < 1e-3))[order(diff_genes$qval)][1:3000]
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

###############CellChat###################
library(CellChat)
library(patchwork)
library(Seurat)
library(rliger)
options(stringsAsFactors = FALSE)

idx=grep('C',seurat@meta.data$orig.ident)
seurat_h=seurat[,idx]
seurat_oa=seurat[,-idx]
####Performed cell chat analysis for normal samples####
seurat<-seurat_h
idx=which(seurat@meta.data$cell.type=='Inflammatory' | seurat@meta.data$cell.type=='Proliferation-2')
seurat=seurat[,-idx]
###require normalized conunts data###
seurat@meta.data$cell.type<-seurat@active.ident
data.input=seurat@assays$RNA@data
meta=seurat@meta.data

cellchat<-createCellChat(object=data.input,meta=meta,group.by='cell.type')

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
#CellChatDB.use=CellChatDB
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame
#We provide a function `subsetCommunication` to easily access the inferred cell-cell communications of interest. For example, 

#df.net <- subsetCommunication(cellchat)``` returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set `slot.name = "netP"` to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))``` gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5. 

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))``` gives the inferred cell-cell communications mediated by signaling WNT and TGFb. 

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_normal <- cellchat
save(cellchat_normal,file='cellchat_normal_sub.rda')

####Performed cell chat analysis for OA samples####
seurat<-seurat_oa
idx=which(seurat@meta.data$cell.type=='Regulator-2')
seurat=seurat[,-idx]
idx=which(seurat@meta.data$cell.type=='Proliferation-2')
seurat=seurat[,-idx]
###require normalized conunts data###
seurat@meta.data$cell.type<-seurat@active.ident
data.input=seurat@assays$RNA@data
meta=seurat@meta.data

cellchat<-createCellChat(object=data.input,meta=meta,group.by='cell.type')

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
#CellChatDB.use=CellChatDB
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Extract the inferred cellular communication network as a data frame
#We provide a function `subsetCommunication` to easily access the inferred cell-cell communications of interest. For example, 

#df.net <- subsetCommunication(cellchat)``` returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set `slot.name = "netP"` to access the the inferred communications at the level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))``` gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5. 

#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))``` gives the inferred cell-cell communications mediated by signaling WNT and TGFb. 

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat_oa <- cellchat
save(cellchat_oa,file='cellchat_oa.rda')

##########Subclustering for Fibrocartilage via Liger###############
counts_list=list()
cout=1
for(orig in unique(seurat@meta.data$orig.ident))
{
idx=which(object@meta.data$orig.ident==orig)
counts_list[[cout]]=object@assays$RNA@counts[,idx]
cout=cout+1
}
names(counts_list)=unique(object@meta.data$orig.ident)

ifnb_liger<-createLiger(counts_list)
ifnb_liger <- normalize(ifnb_liger)
ifnb_liger <- selectGenes(ifnb_liger)
ifnb_liger <- scaleNotCenter(ifnb_liger)
ifnb_liger <- optimizeALS(ifnb_liger, k = 40)
ifnb_liger <- quantile_norm(ifnb_liger)
ifnb_liger <- louvainCluster(ifnb_liger, resolution = 0.2)
ifnb_liger <- runUMAP(ifnb_liger, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
