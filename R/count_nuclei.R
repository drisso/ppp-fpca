suppressPackageStartupMessages(library(anndataR))

file_path <- "/projects/shared/TCGA_data/h5ad"
files <- list.files(file_path)
nuclei <- numeric(length(files))

for(i in seq_along(files)) {
  ann <- read_h5ad(files[i])
  spat <- ann$obsm$spatial
  nuclei[i] <- nrow(spat)
}
