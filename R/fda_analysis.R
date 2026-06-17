suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(fda))
source("plot_funs.R")
source("ppp_fun.R")

# ==========================================
# 2. Define Paths
# ==========================================
data_dir <- "../data/"
input_file <- paste0(data_dir, "combined_Lcross_KDE.csv")

# ==========================================
# 3. Read Data and Apply Function
# ==========================================
cat("Reading data from:", input_file, "\n")
dt <- fread(input_file)
cat("Successfully loaded data with", nrow(dt), "rows and", ncol(dt), "columns.\n")

# Read useful annotation
library(readr)
anno <- read_csv( paste0(data_dir, "TCGA_filtered_metadata_VR.csv"))
dim(anno)

anno <- anno |>
  mutate(
    time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    event = ifelse(!is.na(days_to_death), 1, 0)
  )

patient_ids_from_dt <- substr(colnames(dt)[-1], 1, 12)
dt_ids <- data.frame(id=colnames(dt)[-1], submitter_id=patient_ids_from_dt)
head(dt_ids)

dt_ids |>
  left_join(anno, by = "submitter_id") -> dt_anno

head(dt_anno)

cat("Applying fsmooth function...\n")
smoothed_Lfun <- fsmooth_safe(dt, M = 6, genLfun.fd = 4, centrata = FALSE)

anno_vec <- dt_anno$disease_code
names(anno_vec) <- dt_anno$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored

# fPCA 
#' number of basis function: length(r) = 70,  M = 6
#' K = 70 + 6 − 2 = 74 
fpca <- pca.fd(smoothed_Lfun$fd, nharm = 10, smoothed_Lfun$fdPar)
ncol(dt)
rownames(fpca$scores) <- colnames(dt)[-1] # sample names without r_values

round(fpca$varprop, 5)

# clustering
set.seed(123)
fpca_km <- kmeans(fpca$scores[,1:2], centers = 6, nstart = 10)
fpca_hc <- hclust(dist(fpca$scores[,1:2]), method = "ward.D2")
memb_hc <- cutree(fpca_hc, k = 6)

table(fpca_km$cluster)
table(memb_hc)
table(fpca_km$cluster, memb_hc)

# plot of the L funs
plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun, 
  clusters_vec = memb_hc, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored


# PC scores
par(mfrow=c(1,2))
plot(fpca$scores[,1], fpca$scores[,2], col = fpca_km$cluster, pch = 3, main = "kmeans")
abline(v = 0, lty = "dashed", col = "grey")
abline(h = 0, lty = "dashed", col = "grey")
plot(fpca$scores[,1], fpca$scores[,2], col = memb_hc, pch = 3, main = "Ward.d2")
abline(v = 0, lty = "dashed", col = "grey")
abline(h = 0, lty = "dashed", col = "grey")
par(mfrow=c(1,1))

