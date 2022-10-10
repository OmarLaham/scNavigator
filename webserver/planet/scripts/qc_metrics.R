library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

root.dir = "/app/"
app.dir = paste0(root.dir, "planet/")
scripts.dir = paste0(app.dir, "scripts/")
runs.dir = paste0(app.dir, "media/runs/")

source(paste0(scripts.dir, 'helper_functions.R'))

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
uploadName <- args[[2]];
nFeature.RNA.min <- as.integer(args[[3]]);
nFeature.RNA.max <-  as.integer(args[[4]]);
percentMT <-  as.integer(args[[5]]);


#very important to use same seed so we dont have different results for different runs
set.seed(1234)


print("> Processing sample..")

dir <- paste0(runs.dir, runID, "/data/raw_uploads/", uploadName, "/")

# Load the dataset
data <- load.data(dir, uploadType, uploadName)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data, project = "scNavigator", min.cells = 3, min.features = 200)



data <- process.sample(data, nFeature.RNA.min, nFeature.RNA.max, percentMT)

print("Plotting QC metrics..")
p <- VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path = paste0(runs.dir, runID, "/tmp/"), device = "png", filename = "qc_metrics.png", plot = p)
	print("Saved UMAP Plot")

print("Done and plotted QC metrics..")


