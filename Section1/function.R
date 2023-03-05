scIA <- function(counts_list, method="liger", res_path="./", gsea_method="clusterProfiler", de_method="wilcox", norm_method = "LogNormalize", species="hsa",resolution=0.5, thr_fc=1.0, num_dim=25){
	
	if(method=="liger"){# 
		system("which R")
		library(rliger)
		gc()
		## http://htmlpreview.github.io/?https://github.com/MacoskoLab/liger/blob/master/vignettes/Integrating_scRNA_and_scATAC_data.html
		if(class(counts_list)=="liger"){
			ifnb_liger <- counts_list
			library(ggplot2)
			gg <- rliger::plotByDatasetAndCluster(ifnb_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
			ggplot2::ggsave(paste0(res_path,"/res.ia.clustering.byData.umap.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.png"),plot=gg[[1]] + gg[[2]], dpi=320, height=5, width=14)
		}else if(!file.exists(paste0(res_path,"/res.ia.",paste0(names(counts_list),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.Rdata"))){
			cat(paste0("##Performing Integrative Analysis using LIGER ...\n"))
			ifnb_liger <- rliger::createLiger(counts_list)
			ifnb_liger <- rliger::normalize(ifnb_liger) ## store in ifnb_liger@norm.data
			ifnb_liger <- rliger::selectGenes(ifnb_liger)
			ifnb_liger <- rliger::scaleNotCenter(ifnb_liger) ## store in ifnb_liger@scale.data
			cat(paste0("### optimizeALS ...\n"))
			gc()
			ifnb_liger <- rliger::optimizeALS(ifnb_liger, k = num_dim)
			ifnb_liger <- rliger::quantile_norm(ifnb_liger)
			cat(paste0("### louvainCluster ...\n"))
			ifnb_liger <- rliger::louvainCluster(ifnb_liger, resolution = resolution)
			cat(paste0("### runUMAP ...\n"))
			ifnb_liger <- rliger::runUMAP(ifnb_liger, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
			
			save(ifnb_liger,file=paste0(res_path,"/res.ia.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.Rdata"))
			
			## plot clutering
			library(ggplot2)
			gg <- rliger::plotByDatasetAndCluster(ifnb_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
			ggplot2::ggsave(paste0(res_path,"/res.ia.clustering.byData.umap.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.png"),plot=gg[[1]] + gg[[2]], dpi=320, height=5, width=14)
		}else{
			cat(paste0("/res.ia.",paste0(names(counts_list),collapse="vs"),".PC",num_dim,".liger.Rdata file exists! Loading \n"))
			load(paste0(res_path,"/res.ia.",paste0(names(counts_list),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.Rdata"))
		}# end fi
		## loading specific genes
		# create a folder if does not exist
		library(ggplot2)
		
		cat(paste0("##Plotting the figures ...\n"))
		if(!dir.exists(paste0(res_path,"/de.pc"))){
			dir.create(paste0(res_path,"/de.pc"))
		}# end funcs
		res_hvg_path <- paste0(res_path,"/de.pc")
		gene_loadings <- rliger::plotGeneLoadings(ifnb_liger, do.spec.plot = FALSE, return.plots = TRUE)
		#for(id in 1:num_dim){
		for(id in 1:length(gene_loadings)){
			#ggl <- gene_loadings[[id]]
			if(file.exists(paste0(res_hvg_path,"/res.ia.loading.byData.umap.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".PC",id,".liger.png"))){next;}
			ggplot2::ggsave(paste0(res_hvg_path,"/res.ia.loading.byData.umap.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".PC",id,".liger.png"),plot=gene_loadings[[id]], dpi=320, height=12, width=8)
		}# end funcs
		
		p_wordClouds <- rliger::plotWordClouds(ifnb_liger, num.genes = 10, return.plots = TRUE)
		#p_wordClouds[[1]]
		#ggsave(paste0(res_path,"/res.ia.wordclouds.byData.",paste0(names(ifnb_liger@H),collapse="vs"),".liger.png"))
		for(id in 1:length(p_wordClouds)){
			if(file.exists(paste0(res_hvg_path,"/res.ia.wordclouds.byData.umap.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".PC",id,".liger.png"))){next;}
			ggplot2::ggsave(paste0(res_hvg_path,"/res.ia.wordclouds.byData.umap.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".PC",id,".liger.png"),plot=p_wordClouds[[id]],dpi=320, height=12, width=10)
		}# end funcs
		## save the liger objects
		
		
		############# cluster specific DE
		# create a folder if does not exist
		cat(paste0("##Performing DE analysis ...\n"))
		if(!dir.exists(paste0(res_path,"/de.ct"))){
			dir.create(paste0(res_path,"/de.ct"))
		}# end funcs
		res_de_path <- paste0(res_path,"/de.ct")
		library(openxlsx)
		if(!file.exists(paste0(res_path,"/res.ia.DE.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.wilcox.xlsx"))){
			cluster.results <- rliger::runWilcoxon(ifnb_liger, compare.method = "clusters")
			source("/home/sqsunsph/code/funcs.R")
			#res_annot <- GetGeneAnnotation(hgnc_symbol=cluster.results$feature, species=species)
			#combined_results <- merge(cluster.results,res_annot, by.x="feature", by.y="hgnc_symbol")
			if(species=="hsa"){
				load("/home/sqsunsph/data/gene.sets/gene.annotation.RData")
				combined_results <- merge(cluster.results,gene_annot, by.x="feature", by.y="hgnc_symbol")
			}else{
				combined_results <- cluster.results
			}
			openxlsx::write.xlsx(combined_results, file = paste0(res_path,"/res.ia.DE.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.wilcox.xlsx"), sheetName = "DEA", append = FALSE, row.names=FALSE)
		
			#The number of significant genes identified by runWilcoxon, to filter the output by taking markers which have padj (Benjamini-Hochberg adjusted p-value) less than 0.05 and logFC (log fold change between observations in group versus out) larger than 3:
			#cluster.results <- cluster.results[cluster.results$padj < 0.05,]
			#cluster.results <- cluster.results[cluster.results$logFC > 1,]
			# Cluster 3 and take the top 20 markers by typing these commands:
			for(ict in unique(cluster.results$group)){
				if(file.exists(paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".liger.wilcox.xlsx"))){next;}
				wilcoxon.cluster <- cluster.results[cluster.results$group == ict, ]
				wilcoxon.cluster <- wilcoxon.cluster[order(wilcoxon.cluster$padj), ]
				#markers <- wilcoxon.cluster[1:20, ]
				openxlsx::write.xlsx(wilcoxon.cluster, file = paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".liger.wilcox.xlsx"), sheetName = "DEA.ct", append = FALSE, row.names=FALSE)
				head(wilcoxon.cluster)
			}
		}
		if(!file.exists(paste0(res_path,"/res.ia.DE.byData.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.wilcox.xlsx"))){
			datasets.results <- rliger::runWilcoxon(ifnb_liger, compare.method = "datasets")
			source("/home/sqsunsph/code/funcs.R")
			#res_annot <- GetGeneAnnotation(hgnc_symbol=datasets.results$feature, species=species)
			#combined_results <- merge(datasets.results,res_annot, by.x="feature", by.y="hgnc_symbol")
			#load("/home/sqsunsph/data/gene.sets/gene.annotation.RData")
			#combined_results <- merge(datasets.results, gene_annot, by.x="feature", by.y="hgnc_symbol")
			
			openxlsx::write.xlsx(datasets.results, file = paste0(res_path,"/res.ia.DE.byData.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".liger.wilcox.xlsx"), sheetName = "DEA", append = FALSE, row.names=FALSE)
		}else{
			cat(paste0("/res.ia.DE.byData.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".wilcox.xlsx, file exists! \n"))
		}# end fi
		# return the result if want to plot the figures
		
		############### GSEA
		cat(paste0("##Performing Gene Set Enrichment Analysis ...\n"))
		#gsea_method <- "clusterProfiler"
		de_method <- "wilcox"
		#species <- "hsa"
		#thr_fc <- 1.0
		code_path <- "/home/sqsunsph/projects/code"
		source(paste0(code_path,"/sca_funcs.R"))
		if(!dir.exists(paste0(res_path,"/gsea.ct"))){
			dir.create(paste0(res_path,"/gsea.ct"))
		}# end funcs
		res_gsea_path <- paste0(res_path,"/gsea.ct")
		## does not distinct up and down regularted genes
		for(idb in c("GO", "KEGG")){
			for(ict in unique(ifnb_liger@clusters)){
				if(file.exists(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx"))){
					cat(paste0("The file: res.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".xlsx, exists! \n"))
					next;
				}# end fi
				cat(paste0(idb,": Cell-type ", ict, " ....\n"))
				markers <- openxlsx::read.xlsx(paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".liger.wilcox.xlsx"), rowNames=TRUE)
				summ_stats <- data.frame(symbol=rownames(markers), pvalue=markers$pval, adjusted_pvalue=markers$padj, stats=markers$logFC)
				summ_stats <- na.omit(summ_stats)
				if(any(is.na(summ_stats))){next;}
				summ_stats$g <- ifelse(summ_stats$adjusted_pvalue < 0.05, 'UP', 'stable' )
				table(summ_stats$g)
				#print("run scGSEA before")
				res <- scGSEA(summ_stats, species=species, database=idb, method=gsea_method)
				#print("run scGSEA after")
				if(is.null(res)){cat(paste0("## There is no DE genes (from cell-type ",ict , ") enrich on ", idb, " database! \n")); next;}
				ggrd <- clusterProfiler::dotplot(res)
				print("run scGSEA ggsave before")
				ggplot2::ggsave(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".liger.png"), plot=ggrd, dpi=320, height=5, width=10)
				print("run scGSEA ggsave after")
				openxlsx::write.xlsx(res, file = paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx"), sheetName = idb, append = FALSE)
				rm(res)
			}# end for	
		}# end for
		
		####### upreg and downreg
		for(idb in c("GO", "KEGG")){
			for(ict in unique(ifnb_liger@clusters)){
				if(file.exists(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".upreg.",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx"))){
					cat(paste0("The file: res.upreg.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx, exists! \n"))
					next;
				}# end fi
				cat(paste0(idb,": Cell-type ", ict, " ....\n"))
				markers <- openxlsx::read.xlsx(paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".liger.wilcox.xlsx"), rowNames=TRUE)
				summ_stats <- data.frame(symbol=rownames(markers), pvalue=markers$pval, adjusted_pvalue=markers$padj, stats=markers$logFC)
				summ_stats <- na.omit(summ_stats)
				# consider the log2FC>0
				summp_statas <- summ_stats[summ_stats$stats>thr_fc,]
				if(any(is.na(summp_statas))){next;}
				summp_statas$g <- ifelse(summp_statas$adjusted_pvalue < 0.05, 'UP', 'stable' )
				table(summp_statas$g)
				res <- scGSEA(summp_statas, species=species, database=idb, method=gsea_method)
				if(is.null(res)){cat(paste0("## There is no DE genes (from cell-type ",ict , ") enrich on ", idb, " database! \n")); next;}
				ggrd <- clusterProfiler::dotplot(res)
				ggplot2::ggsave( paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".upreg.",gsea_method,".",de_method,".ct",ict,".",idb,".liger.png"), plot=ggrd, dpi=320, height=5, width=10)
				openxlsx::write.xlsx(res, file = paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".upreg.",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx"), sheetName = idb, append = FALSE)
				rm(res)
				
				if(file.exists(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".downreg.",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx"))){
					cat(paste0("The file: res.downreg.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx, exists! \n"))
					next;
				}# end fi
				# consider the log2FC<0
				summn_stats <- summ_stats[summ_stats$stats< (-thr_fc),]
				if(any(is.na(summn_stats))){next;}
				summn_stats$g <- ifelse(summn_stats$adjusted_pvalue < 0.05, 'UP', 'stable')
				table(summn_stats$g)
				res <- scGSEA(summn_stats, species=species, database=idb, method=gsea_method)
				if(is.null(res)){cat(paste0("## There is no DE genes (from cell-type ",ict , ") enrich on ", idb, " database! \n")); next;}
				ggru <- clusterProfiler::dotplot(res)
				ggplot2::ggsave( paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".downreg.",gsea_method,".",de_method,".ct",ict,".",idb,".liger.png"), plot=ggru, dpi=320, height=5, width=10)
				openxlsx::write.xlsx(res, file = paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ifnb_liger@H),collapse="vs"),".rs",resolution,".PC",num_dim,".downreg.",gsea_method,".",de_method,".ct",ict,".",idb,".liger.xlsx"), sheetName = idb, append = FALSE)
				rm(res)
			}# end for	
		}# end for
		return(ifnb_liger)
	}else if(method=="Seurat"){
		library(Seurat)
		##lapply(ifnb.list, FUN=function(x){names(ncol(x))})
		ID <- unlist( lapply(X = counts_list, FUN = function(x) { return(ncol(x)) }) )
		if(file.exists(paste0(res_path,"/res.ia.",paste0(names(ID),collapse="vs"),".",norm_method,".F4q.seurat.Rdata"))){
			load(paste0(res_path,"/res.ia.",paste0(names(ID),collapse="vs"),".",norm_method,".F4q.seurat.Rdata"))
		}else{
			## normalize and identify variable features for each dataset independently
			counts_list <- lapply(X = counts_list, FUN = function(x) {
				x <- NormalizeData(x)
				x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
			})
		
			## select features that are repeatedly variable across datasets for integration
			features <- SelectIntegrationFeatures(object.list = counts_list)
			
			if(norm_method == "SCT"){
				anchors <- FindIntegrationAnchors(object.list = counts_list, normalization.method = "SCT", anchor.features = features)
				combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
			}else{
				anchors <- FindIntegrationAnchors(object.list = counts_list, anchor.features = features)
				## this command creates an 'integrated' data assay
				combined <- IntegrateData(anchorset = anchors)
			}## end fi
			combined[["ID"]] <- rep(names(ID), times=ID)
			save(combined, file=paste0(res_path,"/res.ia.",paste0(names(ID),collapse="vs"),".",norm_method,".F2.5q.seurat.Rdata"))
		}## end fi
			
		combined[["ID"]] <- rep(names(ID), times=ID)
		
		rm(counts_list)
		gc()
		
		## specify that we will perform downstream analysis on the corrected data note that the
		## original unmodified data still resides in the 'RNA' assay
		DefaultAssay(combined) <- "integrated"
		if(!file.exists(paste0(res_path,"/res.ia.split.byData.umap.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",norm_method,".seurat.png"))){
			# Run the standard workflow for visualization and clustering
			combined <- ScaleData(combined, verbose = FALSE)
			combined <- RunPCA(combined, npcs = num_dim, verbose = FALSE)
			combined <- RunUMAP(combined, reduction = "pca", dims = 1:num_dim)
			combined <- FindNeighbors(combined, reduction = "pca", dims = 1:num_dim)
			combined <- FindClusters(combined, resolution = 0.8)
		
			p1 <- DimPlot(combined, reduction = "umap", group.by = "ID")
			p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
			#p1 + p2
			## plot clutering
			library(ggplot2)
			ggplot2::ggsave(paste0(res_path,"/res.ia.clustering.byData.umap.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",norm_method,".seurat.png"),plot=p1+p2, dpi=320, height=5, width=14)
			
			ggdiv <- DimPlot(combined, reduction = "umap", split.by = "ID")
			ggplot2::ggsave(paste0(res_path,"/res.ia.split.byData.umap.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",norm_method,".seurat.png"),plot=ggdiv, dpi=320, height=5, width=14)
		}## end fi
			
		
		
		## For performing differential expression after integration, we switch back to the original
		## data
		cat(paste0("##Performing DE analysis ...\n"))
		if(!dir.exists(paste0(res_path,"/de.ct"))){
			dir.create(paste0(res_path,"/de.ct"))
		}# end funcs
		res_de_path <- paste0(res_path,"/de.ct")
		library(openxlsx)
		for(ict in levels(Idents(combined))){
			if(file.exists(paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".",norm_method,".seurat.wilcox.xlsx"))){next;}
			DefaultAssay(combined) <- "RNA"
			nk.markers <- FindConservedMarkers(combined, ident.1 = ict, grouping.var = "ID", logfc.threshold = 0.1, test.use = 'wilcox', min.pct = 0,verbose = TRUE)
			head(nk.markers)
			openxlsx::write.xlsx(nk.markers, file = paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".",norm_method,".seurat.wilcox.xlsx"), sheetName = "DEA.ct", append = FALSE, row.names=TRUE)
		}# end for 
		
		############### GSEA
		cat(paste0("##Performing Gene Set Enrichment Analysis ...\n"))
		#gsea_method <- "clusterProfiler"
		de_method <- "wilcox"
		#species <- "hsa"
		#thr_fc <- 1.0
		code_path <- "/home/sqsunsph/projects/code"
		source(paste0(code_path,"/sca_funcs.R"))
		if(!dir.exists(paste0(res_path,"/gsea.ct"))){
			dir.create(paste0(res_path,"/gsea.ct"))
		}# end funcs
		res_gsea_path <- paste0(res_path,"/gsea.ct")
		## does not distinct up and down regularted genes
		for(idb in c("GO", "KEGG")){
			for(ict in levels(Idents(combined))){
				if(file.exists(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".",norm_method,".seurat.xlsx"))){
					#cat(paste0("The file: res.",paste0(unique(combined[["ID"]]$ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".xlsx, exists! \n"))
					cat(paste0("The file: res.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".xlsx, exists! \n"))
					next;
				}# end fi
				cat(paste0(idb,": Cell-type ", ict, " ....\n"))
				markers <- openxlsx::read.xlsx(paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".",norm_method,".seurat.wilcox.xlsx"), rowNames=TRUE)
				summ_stats <- data.frame(symbol=rownames(markers), pvalue=markers$minimump_p_val, adjusted_pvalue=p.adjust(markers$minimump_p_val))
				summ_stats <- na.omit(summ_stats)
				if(any(is.na(summ_stats))){next;}
				summ_stats$g <- ifelse(summ_stats$adjusted_pvalue < 0.05, 'UP', 'stable' )
				table(summ_stats$g)
				#print("run scGSEA before")
				res <- scGSEA(summ_stats, species=species, database=idb, method=gsea_method)
				#print("run scGSEA after")
				if(is.null(res)){cat(paste0("## There is no DE genes (from cell-type ",ict , ") enrich on ", idb, " database! \n")); next;}
				ggrd <- clusterProfiler::dotplot(res)
				print("run scGSEA ggsave before")
				ggplot2::ggsave(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".",norm_method,".seurat.png"), plot=ggrd, dpi=320, height=5, width=10)
				print("run scGSEA ggsave after")
				openxlsx::write.xlsx(res, file = paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".",norm_method,".seurat.xlsx"), sheetName = idb, append = FALSE)
				rm(res)
			}# end for	
		}# end for
		
		####### upreg and downreg
		if(FALSE){for(idb in c("GO", "KEGG")){
			for(ict in levels(Idents(combined))){
				if(file.exists(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".upreg.",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.xlsx"))){
					cat(paste0("The file: res.upreg.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.xlsx, exists! \n"))
					next;
				}# end fi
				cat(paste0(idb,": Cell-type ", ict, " ....\n"))
				markers <- openxlsx::read.xlsx(paste0(res_de_path,"/res.ia.DE.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".CellType",ict,".seurat.wilcox.xlsx"), rowNames=TRUE)
				summ_stats <- data.frame(symbol=rownames(markers), pvalue=markers$pval, adjusted_pvalue=markers$padj, stats=markers$logFC)
				summ_stats <- na.omit(summ_stats)
				# consider the log2FC>0
				summp_statas <- summ_stats[summ_stats$stats>thr_fc,]
				if(any(is.na(summp_statas))){next;}
				summp_statas$g <- ifelse(summp_statas$adjusted_pvalue < 0.05, 'UP', 'stable' )
				table(summp_statas$g)
				res <- scGSEA(summp_statas, species=species, database=idb, method=gsea_method)
				if(is.null(res)){cat(paste0("## There is no DE genes (from cell-type ",ict , ") enrich on ", idb, " database! \n")); next;}
				ggrd <- clusterProfiler::dotplot(res)
				ggplot2::ggsave( paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".upreg.",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.png"), plot=ggrd, dpi=320, height=5, width=10)
				openxlsx::write.xlsx(res, file = paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".upreg.",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.xlsx"), sheetName = idb, append = FALSE)
				rm(res)
				
				if(file.exists(paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".downreg.",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.xlsx"))){
					cat(paste0("The file: res.downreg.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.xlsx, exists! \n"))
					next;
				}# end fi
				# consider the log2FC<0
				summn_stats <- summ_stats[summ_stats$stats< (-thr_fc),]
				if(any(is.na(summn_stats))){next;}
				summn_stats$g <- ifelse(summn_stats$adjusted_pvalue < 0.05, 'UP', 'stable')
				table(summn_stats$g)
				res <- scGSEA(summn_stats, species=species, database=idb, method=gsea_method)
				if(is.null(res)){cat(paste0("## There is no DE genes (from cell-type ",ict , ") enrich on ", idb, " database! \n")); next;}
				ggru <- clusterProfiler::dotplot(res)
				ggplot2::ggsave( paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".downreg.",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.png"), plot=ggru, dpi=320, height=5, width=10)
				openxlsx::write.xlsx(res, file = paste0(res_gsea_path,"/res.ia.GSEA.byCluster.",paste0(names(ID),collapse="vs"),".rs",resolution,".PC",num_dim,".downreg.",gsea_method,".",de_method,".ct",ict,".",idb,".seurat.xlsx"), sheetName = idb, append = FALSE)
				rm(res)
			}# end for	
		}}# end for
		return(combined)
	}# end fi
}# end func
