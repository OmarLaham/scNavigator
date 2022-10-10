library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

source('helper_functions.R')

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadName <- args[[2]];
nDims <- args[[3]];

#very important to use same seed so we dont have different results for different runs
set.seed(1234)

print("> Processing sample..")


dir <- paste0("../media/runs/", runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, runID, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)

data

data <- process.sample(data)

#plot elbow plot
elbow.plot(data)


print("Done and elbow plot saved!")




