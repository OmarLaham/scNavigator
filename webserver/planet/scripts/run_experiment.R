library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadType <- args[[2]];
uploadName <- args[[3]];
expID <- args[[4]];
nDims <- args[[5]];
clusteringResolution <- as.double(args[[6]])

#very important to use same seed so we dont have different results for different runs
set.seed(1234)


process.sample <- function(data, clustering.res, nDims) {

	# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
	data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")

	data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)

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

	#Cluster the cells
	data <- FindNeighbors(data, dims = 1:nDims)
	data <- FindClusters(data, resolution = clustering.res)
	
	#Run non-linear dimensional reduction (UMAP/tSNE) and save plots
	data <- RunUMAP(data, dims = 1:nDims)
	data <- RunTSNE(data, dims = 1:nDims)
	
	#plot UMAP
	p <- DimPlot(data, reduction = "umap")
	ggsave(path = paste0("../media/runs/", runID, "/data/experiments/", expID, "/"), device = "png", filename = "umap.png", plot = p)
	print("Saved UMAP Plot")
	
	#plot t-SNE
	p <- DimPlot(data, reduction = "tsne")
	ggsave(path = paste0("../media/runs/", runID, "/data/experiments/", expID, "/"), device = "png", filename = "t-sne.png", plot = p)
	print("Saved t-SNE Plot")
	
	#save .rds
	saveRDS(data, file=paste0("../media/runs/", runID, "/data/experiments/", expID, "/data.rds"))
	print("Saved data as .rds")	
	
	#return processed seurat object
	return(data)
			
}

# a helper function to identify the root principal points:
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
run.trajectory <- function(data, root.cell.ids = c()) {

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
	
	ggsave(path = paste0("../media/runs/", runID, "/data/experiments/", expID, "/"), device = "png", filename = "trajectory_pseudotime_plot.png", plot = p)
	print("> Saved pseudotime plot.")
	
}

#create exp dir if not exist
mainDir <- paste0("../media/runs/", runID, "/data/experiments/")
subDir <- expID
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}

print("> Processing sample..")

dir <- paste0("../media/runs/", runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- NULL
if(uploadType == "rds") {
	data <- readRDS(paste0(dir, uploadName, ".rds"))
} else if (uploadType == "rds.gz") {
	data <- readRDS(paste0(dir, uploadName, ".rds.gz"))
} else { #10 Genomex
	data <- Read10X(data.dir = dir)
}

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data <- process.sample(data, clustering.res = clusteringResolution, nDims=nDims)

#run trajectory to plot pseudotime trajectory
print("> Trajectory analysis..")
run.trajectory(data)

print("Done with clustering and trajectory. UMAP and t-SNE generated.")


