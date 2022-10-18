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
expTitle <- args[[2]];
uploadName <- args[[3]];
nFeature.RNA.min <- as.integer(args[[4]]);
nFeature.RNA.max <- as.integer(args[[5]]);
percentMT <- as.integer(args[[6]]);
nDims <- args[[7]];
clusteringResolution <- as.double(args[[8]]);
timestamp <- args[[9]];#will be added to img name to avoid web page caching

#very important to use same seed so we dont have different results for different runs
set.seed(1234)



#create exp dir if not exist
mainDir <- paste0(runs.dir, runID, "/data/experiments/")
subDir <- expTitle
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}

print("> Processing sample..")

dir <- paste0(runs.dir, runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, uploadType, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data <- process.sample(data, nFeature.RNA.min, nFeature.RNA.max, percentMT)

print("> Clustering..")
data <- find.clusters(data, clusteringResolution, nDims)

print("> Dimensionality Reduction..")
data <- reduce.dimensions(data, nDims, runs.dir, expTitle, timestamp)

#run trajectory to plot pseudotime trajectory
print("> Trajectory analysis..")
run.trajectory(data, c(), runs.dir, expTitle, timestamp)

print("Done with clustering and trajectory. UMAP and t-SNE generated.")

#DEA
#create dea dir if not exist

expDir <- paste0(mainDir, expTitle, "/");
mainDir <- expDir;
subDir <- "dea"
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}

print("> Running DEA for all clusters")

for(cluster in levels(data)) {

    print(paste("- DEA cluster: ", cluster))

    dea.res <- as.data.frame(FindMarkers(data, ident.1 = cluster, max.cells.per.ident = 50))
    dea.export.path <- paste0(expDir, "dea/cluster_", cluster, "_dea.csv")
    write.csv(dea.res, dea.export.path)

}

print("Done with DEA for all clusters.")

