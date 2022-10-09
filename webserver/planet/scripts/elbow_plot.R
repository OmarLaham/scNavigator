library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadType <- args[[2]];
uploadName <- args[[3]];
nDims <- args[[4]];

#very important to use same seed so we dont have different results for different runs
set.seed(1234)

process.sample <- function(data) {

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

	# Plot the elbow plot to determine number of clusters
	p <- ElbowPlot(object = data, ndims = 50)
	p <- AugmentPlot(p)
	ggsave(p, filename=paste0("../media/runs/", runID, "/tmp/elbow_plot.png"))

	
	#return processed seurat object
	return(data)
			
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

data <- process.sample(data)


print("Done and elbow plot saved!")




