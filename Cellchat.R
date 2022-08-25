library(CellChat)
library(patchwork)
library(Seurat)
library(rliger)
library(SeuratWrappers)
options(stringsAsFactors = FALSE)
###OA samples Convert to Seurat object###

library(SeuratWrappers)
load("/home/public/data/OA/manuscript/suppl/res.ia.10XWKvs10XNKvsSKS1vsSKS2vs10XS3vs10XS4vs10XS5.rs0.7.PC40.liger.rename.Rdata")
seurat <- ligerToSeurat(ifnb_liger,use.liger.genes = TRUE,by.dataset = FALSE,renormalize = TRUE)
seurat<-subset(seurat,idents = c("Proliferation-1","Proliferation-2","Homeostasis","Effector","Prehypertrophy","Hypertrophy","Regulator-1","Regulator-2","Fibrocartilage-1","Inflammatory"))## 去除Unknown
seurat@meta.data$cell.type<-seurat@active.ident
seurat@meta.data$cell.type<-factor(seurat@meta.data$cell.type,levels = c("Proliferation-1","Proliferation-2","Homeostasis","Effector","Prehypertrophy","Hypertrophy","Regulator-1","Regulator-2","Fibrocartilage-1","Inflammatory"))
object=seurat
idx=which(seurat@meta.data$orig.ident=='10XS3' | seurat@meta.data$orig.ident=='10XS4' | seurat@meta.data$orig.ident=='10XS5')
seurat_h=object[,idx]
seurat_oa=object[,-idx]
idx=which(seurat@meta.data$orig.ident=='10XWK' | seurat@meta.data$orig.ident=='10XNK')
seurat_oa=object[,idx]
idx=which(seurat@meta.data$cell.type=='Regulator-2' | seurat@meta.data$cell.type=='Proliferation-2')
seurat_oa=seurat_oa[,-idx]
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
save(cellchat_oa,file='cellchat_oa_sub_whknk.rda')

#######plot OA############
library(CellChat)
load('/home/public/data/OA/manuscript/deconvolution/cellchat/cellchat_oa_sub.rda')
cellchat<-cellchat_oa
###interaction num###
groupSize <- as.numeric(table(cellchat@idents))
pdf("OA_interaction_num_whknk.pdf",width=6,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()
pdf("OA_interaction_strength_whknk.pdf",width=6,height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
dev.off()
###Identify signals contributing most to outgoing or incoming signaling of certain cell groups###
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("OA_in&out_num_whknk.pdf",width=10,height=6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()

load('/home/public/data/OA/manuscript/deconvolution/cellchat/cellchat_normal_sub.rda')
cellchat<-cellchat_normal
###interaction num&strength###
groupSize <- as.numeric(table(cellchat@idents))
pdf("Normal_interaction_num.pdf",width=6,height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()
pdf("Normal_interaction_strength.pdf",width=6,height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
dev.off()
###Identify signals contributing most to outgoing or incoming signaling of certain cell groups###
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("Normal_in&out_num.pdf",width=10,height=6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

save(cellchat,file='cellchat_OA.rda')

###OA all###
load('/home/public/data/OA/manuscript/suppl/res.ia.10XWKvs10XNKvsSKS1vsSKS2.rs0.8.PC30.liger.rename.Rdata')
seurat <- ligerToSeurat(ifnb_liger,use.liger.genes = TRUE,by.dataset = FALSE,renormalize = TRUE)
seurat@meta.data$cell.type<-seurat@active.ident

###require normalized conunts data###
data.input=seurat@assays$RNA@data
meta=seurat@meta.data

cellchat<-createCellChat(object=data.input,meta=meta,group.by='cell.type')

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$count
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

save(cellchat,file='cellchat_OA_all.rda')
###OA samples Convert to Seurat object###
load('/home/public/data/OA/manuscript/suppl/res.ia.NKvsWK.PC30.liger.rename.Rdata')
seurat <- ligerToSeurat(ifnb_liger,use.liger.genes = TRUE,by.dataset = FALSE,renormalize = TRUE)
seurat@meta.data$cell.type<-seurat@active.ident

###require normalized conunts data###
data.input=seurat@assays$RNA@data
meta=seurat@meta.data

cellchat<-createCellChat(object=data.input,meta=meta,group.by='cell.type')

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
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

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

save(cellchat,file='cellchat_OA.rda')
