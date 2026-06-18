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

## Smooth functions
dt_pass_tb <- as_tibble(dt_pass)
smoothed_Lfun_pass <- fsmooth_safe(dt_pass_tb, M = 6, genLfun.fd = 4, centrata = FALSE)

anno_vec <- dt_anno_pass$disease_code
names(anno_vec) <- dt_anno_pass$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_pass, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored +   scale_color_manual(values = pal)

## Perform fPCA
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

## Explore some examples
idx <- names(which.max(fpca_pass$scores[,1]))
plot_sample(idx)

idx <- names(which.min(fpca_pass$scores[,1]))
plot_sample(idx)

idx <- names(which.max(fpca_pass$scores[,2]))
plot_sample(idx)

idx <- names(which.min(fpca_pass$scores[,2]))
plot_sample(idx)


## Let's do another round of removing outliers
## Check where these samples are
dt_anno_pass |>
  mutate(remove = ifelse({PC1>100 | abs(PC2)>10}, TRUE, FALSE)) -> dt_anno_pass
table(dt_anno_pass$remove)

dt_anno_pass |>
  ggplot(aes(PC1, PC2, color = remove)) +
  geom_point(size = 3) +
  theme_bw() +
  geom_vline(xintercept = 100) +
  geom_hline(yintercept = -10) +
  geom_hline(yintercept = 10)

## remove samples
dt_anno_ext <- filter(dt_anno_pass, !remove)
dim(dt_anno_ext)
dt_ext <- cbind(dt[,1], dt[,which(colnames(dt) %in% dt_anno_ext$id)])
dim(dt_ext)

## smooth functions
dt_ext_tb <- as_tibble(dt_ext)
smoothed_Lfun_ext <- fsmooth_safe(dt_ext_tb, M = 6, genLfun.fd = 4, centrata = FALSE)

anno_vec <- dt_anno_ext$disease_code
names(anno_vec) <- dt_anno_ext$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_ext, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored +   scale_color_manual(values = pal)

## Perform fPCA
fpca_ext <- pca.fd(smoothed_Lfun_ext$fd, nharm = 10, smoothed_Lfun_ext$fdPar)
rownames(fpca_ext$scores) <- colnames(dt_ext_tb)[-1] # sample names without r_values

round(fpca_ext$varprop, 5)

dt_anno_ext[,paste0("PC", 1:10)] <- fpca_ext$scores

## Check results
dt_anno_ext |>
  ggplot(aes(PC1, PC2, color = disease_code)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_ext |>
  ggplot(aes(PC1, PC2, color = log(total_nuclei))) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

dt_anno_ext |>
  ggplot(aes(log(total_nuclei), PC1)) +
  geom_point(size = 2, aes(color = disease_code)) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_ext |>
  ggplot(aes(log(total_nuclei), PC1)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal) +
  facet_wrap(~disease_code)

dt_anno_ext |>
  ggplot(aes(PC1, PC2, color = purity)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

dt_anno_ext |>
  ggplot(aes(purity, PC2)) +
  geom_point(size = 2, aes(color = disease_code)) +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_ext |>
  ggplot(aes(purity, PC2)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw() +
  scale_color_manual(values = pal) +
  facet_wrap(~disease_code)

dt_anno_ext |>
  ggplot(aes(disease_code, PC1)) +
  geom_boxplot() +
  theme_bw()

dt_anno_ext |>
  ggplot(aes(disease_code, PC2)) +
  geom_boxplot() +
  theme_bw()

fit1 <- lm(PC1 ~ log(total_nuclei), data=dt_anno_ext)
summary(fit1)
fit2 <- lm(PC1 ~ disease_code, data=dt_anno_ext)
summary(fit2)
fit3 <- lm(PC1 ~ disease_code + log(total_nuclei) + purity + inflammatory + hospital_code, data=dt_anno_ext)
summary(fit3)
fit4 <- lm(PC2 ~ disease_code + log(total_nuclei) + purity + inflammatory + hospital_code, data=dt_anno_ext)
summary(fit4)

## Explore some examples
idx <- names(which.max(fpca_ext$scores[,1]))
plot_sample(idx)

idx <- names(which.min(fpca_ext$scores[,1]))
plot_sample(idx)

idx <- names(which.max(fpca_ext$scores[,2]))
plot_sample(idx)

idx <- names(which.min(fpca_ext$scores[,2]))
plot_sample(idx)

# average curve per disease
plot_group_avg_Lfun(
  smoothed_obj = smoothed_Lfun_ext,
  group_vec = dt_anno_ext$disease_code, 
  mark_i = "neoplastic", 
  mark_j = "stromal", 
  center_plot = TRUE, 
  se = FALSE
) + scale_color_manual(values = pal)

# centroid in fPCA
fpc1 <- tapply(fpca_ext$scores[,1], dt_anno_ext$disease_code, mean)
fpc2 <- tapply(fpca_ext$scores[,2], dt_anno_ext$disease_code, mean)

centroids <- data.frame(PC1=fpc1, PC2=fpc2, disease_code = names(fpc1)) 
centroids |>
  ggplot(aes(PC1, PC2, color = disease_code)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  ggrepel::geom_label_repel(aes(label=disease_code))

dt_anno_ext |>
  ggplot(aes(PC1, PC2, color = disease_code)) +
  geom_point(size = 1) +
  theme_bw() +
  ggrepel::geom_label_repel(data=centroids, aes(label=disease_code), max.overlaps = Inf) +
  theme(legend.position = "none") +
  scale_color_manual(values = pal)

dt_anno_ext <- mutate(dt_anno_ext, class = NA)
dt_anno_ext[dt_anno_ext$disease_code %in% c("OV", "UCEC", "UCS"), "class"] <- "Gynecological"
dt_anno_ext[dt_anno_ext$disease_code %in% c("KIRC", "KIRP", "KICH"), "class"] <- "Renal"
dt_anno_ext[dt_anno_ext$disease_code %in% c("GBM", "LGG", "UVM"), "class"] <- "CNS and eye"
dt_anno_ext[dt_anno_ext$disease_code %in% c("HNSC", "LUSC", "CESC", "ESCA", "BLCA"), "class"] <- "Epithelial"
dt_anno_ext[dt_anno_ext$disease_code %in% c("COAD", "READ", "STAD", "PAAD"), "class"] <- "GI"
dt_anno_ext[dt_anno_ext$disease_code %in% c("DLBC"), "class"] <- "Lymphoma"
dt_anno_ext[dt_anno_ext$disease_code %in% c("SKCM", "SARC"), "class"] <- "Mesenchymal"


dt_anno_ext |>
  dplyr::select(PC1, PC2, class) |>
  na.omit() |>
  ggplot(aes(PC1, PC2, color = class)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_ext |>
  dplyr::select(PC1, PC2, class) |>
  na.omit() |>
  ggplot(aes(class, PC1)) +
  geom_boxplot() +
  theme_bw()

dt_anno_ext |>
  dplyr::select(PC1, PC2, class) |>
  na.omit() |>
  ggplot(aes(class, PC2)) +
  geom_boxplot() +
  theme_bw()

# curves color-coded by PC1 and PC2
anno_vec <- dt_anno_ext$PC1
names(anno_vec) <- dt_anno_ext$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_ext, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored + scale_color_viridis_c()

anno_vec <- dt_anno_ext$PC2
names(anno_vec) <- dt_anno_ext$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_ext, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored + scale_color_viridis_c()

# clustering
set.seed(123)
fpca_km <- kmeans(fpca_ext$scores, centers = 6, nstart = 10)
fpca_hc <- hclust(dist(fpca_ext$scores), method = "ward.D2")
memb_hc <- cutree(fpca_hc, k = 6)

table(fpca_km$cluster)
table(memb_hc)
table(fpca_km$cluster, memb_hc)

# plot of the L funs
plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_ext, 
  clusters_vec = as.factor(memb_hc), 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored

dt_anno_ext <- mutate(dt_anno_ext, kmeans = as.factor(fpca_km$cluster))
dt_anno_ext <- mutate(dt_anno_ext, hclust = as.factor(memb_hc))

dt_anno_ext |>
  ggplot(aes(PC1, PC2, color = kmeans)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = pal)

dt_anno_ext |>
  ggplot(aes(PC1, PC2, color = hclust)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = pal)

## Focus on Lung, Stomach and esophageal, 
dt_anno_brca <- filter(dt_anno_ext, disease_code %in% c("LUAD", "LUSC"))
dt_anno_brca <- filter(dt_anno_ext, disease_code %in% c("STAD", "ESCA"))
dt_anno_brca <- filter(dt_anno_ext, disease_code %in% c("BLCA"))

dt_brca <- cbind(dt[,1], dt[,which(colnames(dt) %in% dt_anno_brca$id)])
dim(dt_brca)

## smooth functions
dt_brca_tb <- as_tibble(dt_brca)
smoothed_Lfun_brca <- fsmooth_safe(dt_brca_tb, M = 6, genLfun.fd = 4, centrata = FALSE)

anno_vec <- dt_anno_brca$disease_code
names(anno_vec) <- dt_anno_brca$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_brca, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored + scale_color_manual(values = pal)

## Perform fPCA
fpca_brca <- pca.fd(smoothed_Lfun_brca$fd, nharm = 10, smoothed_Lfun_brca$fdPar)
rownames(fpca_brca$scores) <- colnames(dt_brca_tb)[-1] # sample names without r_values

round(fpca_brca$varprop, 5)

dt_anno_brca[,paste0("PC", 1:10)] <- fpca_brca$scores

## Check results
dt_anno_brca |>
  ggplot(aes(PC1, PC2, color = log(total_nuclei))) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

dt_anno_brca |>
  ggplot(aes(PC1, PC2, color = purity)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_viridis_c()

anno_vec <- dt_anno_brca$PC1
names(anno_vec) <- dt_anno_brca$id

plot_colored <- plot_smoothed_Lfun(
  smoothed_obj = smoothed_Lfun_brca, 
  clusters_vec = anno_vec, 
  mark_i = "neoplastic",  # Define marks manually here for the title
  mark_j = "stromal",
  center_plot = TRUE      # Set to TRUE to subtract 'r' since centrata=FALSE in fsmooth
)
plot_colored + scale_color_viridis_c()

## Check association with survival
set.seed(123)
fpca_km <- kmeans(fpca_brca$scores, centers = 4, nstart = 10)
table(fpca_km$cluster)
table(fpca_km$cluster, dt_anno_brca$disease_code)

dt_anno_brca <- mutate(dt_anno_brca, kmeans = as.factor(fpca_km$cluster))

dt_anno_brca |>
  ggplot(aes(PC1, PC2, color = kmeans)) +
  geom_point(size = 1) +
  theme_bw() +
  scale_color_manual(values = pal)


# --- Survival analysis by cluster (Kaplan-Meier + log-rank test) ---------------------------------
# Requires: dt_anno_brca with columns 'time' (follow-up), 'event' (0=censored,1=event), and 'kmeans' (cluster factor)
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

# Basic input checks
if (!all(c("time", "event", "kmeans") %in% colnames(dt_anno_brca))) {
  stop("dt_anno_brca must contain 'time', 'event', and 'kmeans' columns for survival analysis.")
}

# Build survival object
surv_obj <- with(dt_anno_brca, Surv(time = time, event = event))

# Fit KM curves by cluster
fit_km <- survfit(surv_obj ~ kmeans, data = dt_anno_brca)

# Plot Kaplan-Meier curves with risk table and p-value
# Use existing palette if available; otherwise let survminer choose
palette_for_clusters <- if (exists("pal") && length(pal) >= length(unique(dt_anno_brca$kmeans))) pal else NULL

km_plot <- ggsurvplot(
  fit_km,
  data = dt_anno_brca,
  risk.table = FALSE,
  pval = TRUE,
  conf.int = FALSE,
  palette = palette_for_clusters,
  ggtheme = theme_minimal(),
  legend.title = "Cluster",
  risk.table.height = 0.2
)

# Print the plot (in interactive sessions this renders; in scripts it will be written if wrapped)
print(km_plot)

# Log-rank (survdiff) test
lr <- survdiff(surv_obj ~ kmeans, data = dt_anno_brca)
# Compute p-value from chi-square distribution
pval_lr <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
cat(sprintf("Log-rank test: chi-square = %.3f, df = %d, p = %g\n", lr$chisq, length(lr$n) - 1, pval_lr))

# Optionally, show pairwise comparisons (Benjamini-Hochberg adjusted) using pairwise_survdiff if survminer available
if (requireNamespace("survminer", quietly = TRUE)) {
  if (length(unique(dt_anno_brca$kmeans)) > 1) {
    pw <- survminer::pairwise_survdiff(Surv(time, event) ~ kmeans, data = dt_anno_brca)
    cat("Pairwise log-rank p-values (raw):\n")
    print(pw$p.value)
  }
}

# -----------------------------------------------------------------------------------------------

# average curve per cluster
plot_group_avg_Lfun(
  smoothed_obj = smoothed_Lfun_brca,
  group_vec = dt_anno_brca$kmeans, 
  mark_i = "neoplastic", 
  mark_j = "stromal", 
  center_plot = TRUE, 
  se = FALSE
) + scale_color_manual(values = pal)

idx <- names(fpca_brca$scores[dt_anno_brca$kmeans==1,1])[1]
plot_sample(idx)

idx <- names(fpca_brca$scores[dt_anno_brca$kmeans==2,1])[1]
plot_sample(idx)

idx <- names(fpca_brca$scores[dt_anno_brca$kmeans==3,1])[1]
plot_sample(idx)

idx <- names(fpca_brca$scores[dt_anno_brca$kmeans==4,1])[1]
plot_sample(idx)
