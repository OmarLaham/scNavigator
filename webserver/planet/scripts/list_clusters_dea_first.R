library(dplyr)
library(Seurat)
library(patchwork)

root.dir = "/app/"
app.dir = paste0(root.dir, "planet/")
scripts.dir = paste0(app.dir, "scripts/")
runs.dir = paste0(app.dir, "media/runs/")


source(paste0(scripts.dir, 'helper_functions.R'))

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE);
runID <- args[[1]];
expTitle <- args[[2]];

exp.dir <- paste0(runs.dir, runID, "/data/experiments/", expTitle, "/")


#very important to use same seed so we dont have different results for different runs
set.seed(1234)


#create dea dir if not exist
mainDir <- paste0(exp.dir)
subDir <- "dea"
if(dir.exists(file.path(mainDir, subDir)) == FALSE) {
	dir.create(file.path(mainDir, subDir))
}


print("> Running DEA on all clusters..")

data.rds.path <- paste0(exp.dir, "data.rds")

# Load the dataset
data <- readRDS(data.rds.path)


#export available clusters
clusters <- as.data.frame(levels(data))
clusters.export.path <- paste0(exp.dir, "dea/", "clusters.csv")
write.csv(clusters, clusters.export.path)

#run DEA for first cluster only


cluster.num <- 0
print(paste("Cluster:", cluster.num, ".."))

dea.res <- as.data.frame(FindMarkers(data, ident.1 = cluster.num, max.cells.per.ident = 50))
dea.export.path <- paste0(exp.dir, "dea/cluster_", cluster.num, "_dea.csv")
write.csv(dea.res, dea.export.path)

print("Done with DEA for all clusters.")


