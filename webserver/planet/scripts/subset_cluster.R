library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratWrappers)
library(monocle3)

root.dir = "/app/"
app.dir = paste0(root.dir, "planet/")
scripts.dir = paste0(app.dir, "scripts/")
runs.dir = paste0(app.dir, "media/runs/")

source(paste0(scripts.dir, 'helper_functions.R'))

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadName <- args[[2]];
expTitle <- args[[3]];
cluster <- args[[4]];
nFeature.RNA.min <- as.integer(args[[5]]);
nFeature.RNA.max <- as.integer(args[[6]]);
percentMT <- as.integer(args[[7]]);
nDims <- args[[8]];
clusteringResolution <- as.double(args[[9]]);

#very important to use same seed so we dont have different results for different runs
set.seed(1234)

subset.folder.name <- paste(uploadName, expTitle, cluster, sep="_")

#create exp dir if not exist
mainDir <- paste0(runs.dir, runID, "/data/raw_uploads/")
subDir <- subset.folder.name
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}

print("> subsetting sample..")

dir <- paste0(runs.dir, runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, uploadType, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data <- process.sample(data, nFeature.RNA.min, nFeature.RNA.max, percentMT)

print("> Clustering..")
data <- find.clusters(data, clusteringResolution, nDims)

data <- subset(data, seurat.clusters = cluster)

subset.export.path <- paste0(runs.dir, runID, "/data/raw_uploads/", subset.folder.name)
saveRDS(data, subset.export.path)

print("Done with clustering and trajectory. UMAP and t-SNE generated.")


