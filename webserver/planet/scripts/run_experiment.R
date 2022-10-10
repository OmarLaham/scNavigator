library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

root.dir = "/app/"
app.dir = paste0(root.dir, "planet/")
scripts.dir = paste0(root.dir, "scripts/")
runs.dir = paste0(app.dir, "media/runs/")

source('helper_functions.R')

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
expTitle <- args[[2]];
uploadName <- args[[3]];
nFeature.RNA.min <- as.double(args[[4]]);
nFeature.RNA.max <- as.double(args[[5]]);
percentMT <- as.double(args[[6]]);
nDims <- args[[7]];
clusteringResolution <- as.double(args[[8]])


#very important to use same seed so we dont have different results for different runs
set.seed(1234)



#create exp dir if not exist
mainDir <- paste0("../media/runs/", runID, "/data/experiments/")
subDir <- expTitle
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}

print("> Processing sample..")

dir <- paste0("../media/runs/", runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, uploadType, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data <- process.sample(data, nFeature.RNA.min, nFeature.RNA.max, percentMT)

print("> Clustering..")
data <- find.clusters(data, clusteringResolution, nDims)

print("> Dimensionality Reduction..")
data <- reduce.dimensions(data, nDims, runs.dir)

#run trajectory to plot pseudotime trajectory
print("> Trajectory analysis..")
run.trajectory(data, c(), runs.dir)

print("Done with clustering and trajectory. UMAP and t-SNE generated.")


