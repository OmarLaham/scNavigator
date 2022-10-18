library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

root.dir = "/app/"
app.dir = paste0(root.dir, "planet/")
scripts.dir = paste0(root.dir, "scripts/")
runs.dir = paste0(app.dir, "media/runs/")

#load data depending on upload type
load.data <-function(dir, runID, upload.name) {


	rds.file.name <- paste0(dir, upload.name, ".rds")
	rds.gz.file.name <- paste0(dir, upload.name, ".rds")

	data <- NULL

	if(dir.exists(rds.file.name) == TRUE) {
		print(">> Loading .rds..")
		data <- readRDS(rds.file.name)
	} else if (dir.exists(rds.gz.file.name) == TRUE) {
		print(">> Loading .rds.gz..")
		data <- readRDS(rds.gz.file.name)
	} else { #10Genomex
		print(">> Loading 10 Genomex..")
		data <- Read10X(data.dir = dir)
	}

	return (data)

}

# filter, normalized, find variable features, scale and run PCA
process.sample <- function(data, nFeature.RNA.min, nFeature.RNA.max, percentMT) {

	# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
	data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

	data <- subset(data, subset = nFeature_RNA > nFeature.RNA.min & nFeature_RNA < nFeature.RNA.max & percent.mt < percentMT)

	#Normalizing the data
	data <- NormalizeData(data)
	print("> Normalized data")

	#Identification of highly variable features (feature selection)
	data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

	# Scaling the data
	all.genes <- rownames(data)
	data <- ScaleData(data, features = all.genes)
	print("> Scaled data")

	#Perform linear dimensional reduction
	print("> Running PCA..")
	data <- RunPCA(data, features = VariableFeatures(object = data))
	print("> PCA done!")

	return(data)

}

#plot elbow plot
elbow.plot <- function(data, runs.dir, n.dims, timestamp) {

	# Plot the elbow plot to determine number of clusters
	p <- ElbowPlot(object = data, ndims = n.dims)
	p <- AugmentPlot(p)
	ggsave(p, filename=paste0(runs.dir, runID, "/tmp/elbow_plot_", timestamp, ".png"))


	#return processed seurat object
	return(data)

}

#clustering
find.clusters <- function(data, clustering.res, nDims) {


	#Cluster the cells
	data <- FindNeighbors(data, dims = 1:nDims)
	data <- FindClusters(data, resolution = clustering.res)

}
#dimensionality reduction UMAP and t-SNE
reduce.dimensions <- function(data, nDims, runs.dir, expTitle, timestamp) {

	#Run non-linear dimensional reduction (UMAP/tSNE) and save plots
	data <- RunUMAP(data, dims = 1:nDims)
	data <- RunTSNE(data, dims = 1:nDims)

	#plot UMAP
	p <- DimPlot(data, reduction = "umap")
	ggsave(path = paste0(runs.dir, runID, "/data/experiments/", expTitle, "/"), device = "png", filename = paste0("umap_", timestamp, ".png"), plot = p)
	print("Saved UMAP Plot")

	#export UMAP cell embeddings (coords) for better visualization than Seurat plots
	umap.cell.embeddings <- as.data.frame(data[["umap"]]@cell.embeddings)
	write.csv(umap.cell.embeddings, paste0(runs.dir, runID, "/data/experiments/", expTitle, "/umap_cell_embeddings.csv"))


	#plot t-SNE
	p <- DimPlot(data, reduction = "tsne")
	ggsave(path = paste0(runs.dir, runID, "/data/experiments/", expTitle, "/"), device = "png", filename = paste0("tsne_", timestamp, ".png"), plot = p)
	print("Saved t-SNE Plot")

	#export t-SNE cell embeddings (coords) for better visualization than Seurat plots
	tsne.cell.embeddings <- as.data.frame(data[["tsne"]]@cell.embeddings)
	write.csv(tsne.cell.embeddings, paste0(runs.dir, runID, "/data/experiments/", expTitle, "/tsne_cell_embeddings.csv"))

	#save .rds
	saveRDS(data, file=paste0(runs.dir, runID, "/data/experiments/", expTitle, "/data.rds"))
	print("Saved data as .rds")

	#return processed seurat object
	return(data)

}

# a helper function to identify the root principal points for trajectory analysis:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- rownames(colData(cds))#which(colData(cds)[, "embryo.time.bin"] == time_bin)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}


#data is your seurat object
run.trajectory <- function(data, root.cell.ids = c(), runs.dir, expTitle, timestamp) {

	#We can convert the Seurat object to a CellDataSet object using the as.cell_data_set() function from SeuratWrappers

	#cds means cell data set
	data.cds <- as.cell_data_set(data)
	data.cds <- cluster_cells(cds = data.cds, reduction_method = "UMAP")
	data.cds <- learn_graph(data.cds, use_partition = TRUE)

	# order cells
	if (length(root.cell.ids) == 0) {
		data.cds <- order_cells(data.cds, root_pr_nodes=get_earliest_principal_node(data.cds))

	} else {
		data.cds <- order_cells(data.cds, reduction_method = "UMAP", root_cells = root.cell.ids)
	}

	#plot trajectories colored by pseudotime and dave
	p <- plot_cells(
	  cds = data.cds,
	  color_cells_by = "seurat_clusters",
	  group_cells_by = "cluster",
	  label_groups_by_cluster = T,
	  show_trajectory_graph = TRUE,
	  label_cell_groups = F,
	  group_label_size = 3
	)

	ggsave(path = paste0(runs.dir, runID, "/data/experiments/", expTitle, "/"), device = "png", filename = paste0("trajectory_pseudotime_plot_", timestamp, ".png"), plot = p)
	print("> Saved pseudotime plot.")

}
