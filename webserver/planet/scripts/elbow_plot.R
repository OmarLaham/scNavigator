library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

root.dir = "/app/"
app.dir = paste0(root.dir, "planet/")
scripts.dir = paste0(root.dir, "scripts/")
runs.dir = paste0(app.dir, "media/runs/")

source('helper_functions.R')

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadName <- args[[2]];
nFeature.RNA.min <- args[[3]];
nFeature.RNA.max <-  args[[4]];
percentMT <-  args[[5]];
nDims <- args[[6]];

#very important to use same seed so we dont have different results for different runs
set.seed(1234)

print("> Processing sample..")


dir <- paste0("../media/runs/", runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, runID, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data

data <- process.sample(data, nFeature.RNA.min, nFeature.RNA.max, percentMT)

#plot elbow plot
elbow.plot(data, runs.dir)


print("Done and elbow plot saved!")




