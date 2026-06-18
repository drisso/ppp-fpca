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

table(anno_vec)

library(Polychrome)
pal <- unname(palette36.colors(n = 31))

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored +   scale_color_manual(values = pal)

# fPCA 
#' number of basis function: length(r) = 70,  M = 6
#' K = 70 + 6 − 2 = 74 
fpca <- pca.fd(smoothed_Lfun$fd, nharm = 10, smoothed_Lfun$fdPar)
rownames(fpca$scores) <- colnames(dt)[-1] # sample names without r_values

round(fpca$varprop, 5)

df_pca <- data.frame(fpca$scores, anno=anno_vec)
head(df_pca)

df_pca |>
  ggplot(aes(X1, X2, color = anno)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = -50)
  
plot(mean.fd(smoothed_Lfun$fd))  
plot(fpca$harmonics[1], main="PC1")
plot(fpca$harmonics[2], main="PC2")

hist(fpca$scores[,1])

## Identify and inspect outliers
samples <- union(names(which(fpca$scores[,1]>500)), names(which(fpca$scores[,2]< (-50))))
samples

library(imageFeatureTCGA)
library(HistoImagePlot)

hn_spe_list <- list()
thumb_paths <- character(length(samples))

for(i in seq_along(samples)) {
hov_files <- paste0(
    "https://store.cancerdatasci.org/hovernet/h5ad/",
    samples[i],
    ".h5ad.gz")

thumb_path[i] <- paste0(
    "https://store.cancerdatasci.org/hovernet/thumb/",
    samples[i],
    ".png")

hn_spe_list[[i]] <- HoverNet(hov_files, outClass = "SpatialExperiment") |>
    import()

}

for(i in seq_along(samples)) {
print(plotHoverNetH5ADOverlay(hn_spe_list[[i]], thumb_path[i], 
  title=substr(samples[i], 1, 23),
  point_size = 0.02,
  legend_point_size = 3))
}

## Clearly, the outliers are weird samples with very few nuclei.
lapply(hn_spe_list, ncol)

## Let's look more generally at the relation between PCs and some general stats
nuclei_stats <- read_csv("../data/nuclei_stats.csv")

dt_anno |>
  left_join(nuclei_stats, by = c("id" = "sample_id")) -> dt_anno

scores <- fpca$scores
colnames(scores) <- paste0("PC", seq_len(ncol(scores)))

dt_anno <- cbind(dt_anno, scores)

dt_anno |>
  ggplot(aes(PC1, PC2, color = disease_code)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = pal) +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = -50)

dt_anno |>
  ggplot(aes(PC1, PC2, color = log(total_nuclei))) +
  geom_point(size = 3) +
  theme_bw() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = -50) +
  scale_color_viridis_c()

dt_anno |>
  ggplot(aes(log(total_nuclei), log(PC1), color = disease_code)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno |>
  ggplot(aes(PC1, PC2, color = purity)) +
  geom_point(size = 3) +
  theme_bw() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = -50) +
  scale_color_viridis_c()

dt_anno |>
  ggplot(aes(total_nuclei)) +
  geom_histogram()

summary(dt_anno$total_nuclei)
table(dt_anno$total_nuclei > 40000)

## Since we are computing the L curves from a sketch of 40k nuclei, let's remove all samples with less than that and repeat PCA.

## Check where these samples are
dt_anno |>
  mutate(low_quality = ifelse(total_nuclei > 40000, FALSE, TRUE)) -> dt_anno
table(dt_anno$low_quality)

dt_anno |>
  ggplot(aes(PC1, PC2, color = low_quality)) +
  geom_point(size = 3) +
  theme_bw() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = -50)

## remove samples
dt_anno_pass <- filter(dt_anno, !low_quality)
dim(dt_anno_pass)
dt <- as.data.frame(dt)
dt_pass <- cbind(dt[,1], dt[,which(colnames(dt) %in% dt_anno_pass$id)])
dim(dt_pass)

## Perform fPCA
dt_pass_tb <- as_tibble(dt_pass)
smoothed_Lfun_pass <- fsmooth_safe(dt_pass_tb, M = 6, genLfun.fd = 4, centrata = FALSE)
fpca_pass <- pca.fd(smoothed_Lfun_pass$fd, nharm = 10, smoothed_Lfun_pass$fdPar)
rownames(fpca_pass$scores) <- colnames(dt_pass_tb)[-1] # sample names without r_values

round(fpca_pass$varprop, 5)

dt_anno_pass[,paste0("PC", 1:10)] <- fpca_pass$scores

## Check results
dt_anno_pass |>
  ggplot(aes(PC1, PC2, color = disease_code)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_pass |>
  ggplot(aes(PC1, PC2, color = log(total_nuclei))) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

dt_anno_pass |>
  ggplot(aes(log(total_nuclei), PC1)) +
  geom_point(size = 2, aes(color = disease_code)) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_pass |>
  ggplot(aes(log(total_nuclei), PC1)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal) +
  facet_wrap(~disease_code)

dt_anno_pass |>
  ggplot(aes(PC1, PC2, color = purity)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

dt_anno_pass |>
  ggplot(aes(purity, PC2)) +
  geom_point(size = 2, aes(color = disease_code)) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_pass |>
  ggplot(aes(purity, PC2)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal) +
  facet_wrap(~disease_code)

dt_anno_pass |>
  ggplot(aes(disease_code, PC1)) +
  geom_boxplot() +
  theme_bw()

dt_anno_pass |>
  ggplot(aes(disease_code, PC2)) +
  geom_boxplot() +
  theme_bw()

dt_anno_pass |>
  ggplot(aes(disease_code, purity)) +
  geom_boxplot() +
  theme_bw()

fit1 <- lm(PC1 ~ log(total_nuclei), data=dt_anno_pass)
summary(fit1)
fit2 <- lm(PC1 ~ disease_code, data=dt_anno_pass)
summary(fit2)
fit3 <- lm(PC1 ~ disease_code + log(total_nuclei) + purity + inflammatory + hospital_code, data=dt_anno_pass)
summary(fit3)
fit4 <- lm(PC2 ~ disease_code + log(total_nuclei) + purity + inflammatory + hospital_code, data=dt_anno_pass)
summary(fit4)

# average curve per disease

# centroid in fPCA



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

