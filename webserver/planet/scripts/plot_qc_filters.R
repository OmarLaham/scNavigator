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
nFeature.RNA.min <- as.double(args[[6]]);
nFeature.RNA.max <- as.double(args[[7]]);
percent.mt <- as.double(args[[8]]);

#very important to use same seed so we dont have different results for different runs
set.seed(1234)


print("> Processing sample..")

dir <- paste0("../media/runs/", runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, uploadType, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)



data <- process.sample(data)

print("Plotting QC metrics..")
p <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path = paste0("../media/runs/", runID, "/tmp/"), device = "png", filename = "qc_metrics.png", plot = p)
	print("Saved UMAP Plot")

print("Done and plotted QC metrics..")


