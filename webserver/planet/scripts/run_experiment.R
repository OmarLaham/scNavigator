library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

source('helper_functions.R')

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadType <- args[[2]];
uploadName <- args[[3]];
expID <- args[[4]];
nFeature.RNA.min <- as.double(args[[5]]);
nFeature.RNA.max <- as.double(args[[6]]);
percent.mt <- as.double(args[[7]]);
nDims <- args[[8]];
clusteringResolution <- as.double(args[[9]])


#very important to use same seed so we dont have different results for different runs
set.seed(1234)



#create exp dir if not exist
mainDir <- paste0("../media/runs/", runID, "/data/experiments/")
subDir <- expID
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}

print("> Processing sample..")

dir <- paste0("../media/runs/", runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, uploadType, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data <- process.sample(data)

print("> Clustering..")
data <- find.clusters(data, clusteringResolution, nDims)

print("> Dimensionality Reduction..")
data <- reduce.dimensions(data, nDims)

#run trajectory to plot pseudotime trajectory
print("> Trajectory analysis..")
run.trajectory(data)

print("Done with clustering and trajectory. UMAP and t-SNE generated.")


