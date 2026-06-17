suppressPackageStartupMessages(library(anndataR))
suppressPackageStartupMessages(library(fda))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(spatstat))
#suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(tidyr))


# Function to Convert h5ad to ppp object ----
ann2ppp <- function(data_path, mask = NULL, magnif, spat_coords_name = c("x", "y"), marks_col = "type") {
  
  #' magnif is the value of the magnification parameter in the metadata (can be 20 or 40 (20x or 40x respectively)) 
  #' it's probably an error if it's noted as 20x.
  
  # 1. Read AnnData
  ann <- read_h5ad(data_path)
  
  # 2. Extract spatial coordinates
  spat <- ann$obsm$spatial
  x_coords <- spat[, 1]
  y_coords <- spat[, 2]
  
  # 3. Extract marks (e.g., cell types)
  marks_data <- NULL
  if (marks_col %in% colnames(ann$obs)) {
    marks_data <- as.factor(ann$obs[[marks_col]])
  }
  
  # 4. Define the window (owin)
  if (!is.null(mask)) {
    # Extract the 2D matrix from the imager cimg object (dropping z and c dimensions).
    # imager keeps 'x' as rows and 'y' as columns, but spatstat's owin mask 
    # expects 'y' as rows and 'x' as columns, so we transpose it.
    # We also convert it to a logical matrix (TRUE for foreground > 0)
    mask_mat <- t(mask[, , 1, 1] > 0.5) 

    # Create an owin using the image mask and the ranges of the spatial coordinates
    win_thumb <- spatstat.geom::owin(
      xrange = c(0, dim(mask)[1]), 
      yrange = c(0, dim(mask)[2]), 
      mask = mask_mat
    )
    # C. Calculate exact scale
    exact_scale <- magnif / 1.25
    
    # D. Transform the window in the full resolution scale of the coordinates
    final_win <- spatstat.geom::affine(win_thumb, mat = diag(c(exact_scale, exact_scale)))
    
  } else {
    # Create a basic bounding box window if no mask is provided
    final_win <- spatstat.geom::owin(
      xrange = range(x_coords, na.rm = TRUE), 
      yrange = range(y_coords, na.rm = TRUE)
    )
  }
  
  # 5. Build the PPP object
  p <- spatstat.geom::ppp(
    x = x_coords, 
    y = y_coords, 
    window = final_win,
    marks = marks_data
  )
  
  return(p)
}

# Lcross estimate with down-sampling ----
Lcross_ds <- function(ppp_obj, mark_i, mark_j, mark_column, downsample = TRUE, n_down = 10000, ...) {
  
  # 0. If marks are a dataframe, extract the specific column we want to use!
  if (is.data.frame(marks(ppp_obj))) {
    marks(ppp_obj) <- marks(ppp_obj)[[mark_column]]
  }
  # 1. Get the indices of the specific cell types dynamically
  idx_i <- which(marks(ppp_obj) == mark_i)
  idx_j <- which(marks(ppp_obj) == mark_j)
  
  # Safety check: ensure both cell types actually exist in the image!
  if(length(idx_i) < 2 | length(idx_j) < 2) {
    warning(paste("Not enough points for", mark_i, "and/or", mark_j, "to compute Lcross. Returning NULL."))
    return(NULL)
  }
  
  # 2. Handle Downsampling
  if (downsample) {
    # Randomly sample a maximum of n_down cells of each type
    keep_i <- sample(idx_i, min(length(idx_i), n_down))
    keep_j <- sample(idx_j, min(length(idx_j), n_down))
    
    # Subset the ppp object
    ppp_sub <- ppp_obj[c(keep_i, keep_j)]
    
  } else {
    # Even if we don't downsample, subsetting to JUST the two marks of interest 
    # saves a massive amount of memory for spatstat calculation
    ppp_sub <- ppp_obj[c(idx_i, idx_j)]
  }
  
  # 3. Calculate Lcross dynamically using the provided mark names
  L <- Lcross(ppp_sub, 
              i = mark_i, 
              j = mark_j, ...)
  
  return(L)
}

# The spatial sketching helper function
sketch_spatial_R <- function(ppp_obj, idx, n_target, n_threads = 1) {
  if (length(idx) <= n_target) return(idx)
  
  coords <- cbind(x = ppp_obj$x[idx], y = ppp_obj$y[idx])
  # fastKeamns from Morpho
  km <- Morpho::fastKmeans(coords, k = n_target, iter.max = 15, threads = n_threads)
  
  # Fast Vectorized Sampling
  # Randomly shuffle all the indices
  rand_order <- sample(seq_along(idx))
  # Match the clusters to the shuffled order
  shuffled_clusters <- km$class[rand_order]
  # Pick the first occurrence
  first_occurrences <- !duplicated(shuffled_clusters)
  # Extract the selected relative indices
  sketched_rel_idx <- rand_order[first_occurrences]
  
  return(idx[sketched_rel_idx])
}

# sketching with Kernel Density Estimate
downsample_KDE <- function(ppp_obj, marks_name, n_target, ...) {
  
  # 1. Subset to the specified marks using safe spatstat API
  pp_neo <- subset(ppp_obj, marks %in% marks_name)
  spatstat.geom::marks(pp_neo) <- factor(spatstat.geom::marks(pp_neo))
  
  # 2. Split into a named list based on the marks
  pp_split <- spatstat.geom::split.ppp(pp_neo)
  
  # 3. Process each split in parallel
  sampled_list <- mclapply(pp_split, function(pp) {
    
    n_available <- spatstat.geom::npoints(pp)
    
    # Safely handle empty patterns
    if (n_available == 0) return(pp[integer(0)])
    
    # Skip KDE if we don't have enough points to downsample
    if (n_available <= n_target) return(pp)
    
    # Calculate KDE safely with error handling for bandwidth selection
    kern <- tryCatch({
      spatstat.explore::density.ppp(pp, at="pixels", ...)
    }, error = function(e) {
      # Fallback: if auto-bandwidth fails, try a simple rule or fixed sigma
      # Adjust or expose this fallback based on your data scale
      spatstat.explore::density.ppp(pp, at="pixels", sigma = spatstat.explore::bw.ppl(pp))
    })
    
    xy <- spatstat.geom::coords(pp)
    pint <- spatstat.geom::lookup.im(kern, xy$x, xy$y)
    
    good <- which(is.finite(pint) & !is.na(pint) & pint > 0)
    if (length(good) == 0) return(pp[integer(0)])
    if (length(good) <= n_target) return(pp[good])
    
    w <- 1 / pint[good]
    w <- w / sum(w)
    
    # Defensive programming: use a robust sampling wrapper to avoid the sample(x) gotcha
    sampled_idx <- if (length(good) == 1) good else sample(good, n_target, replace = FALSE, prob = w)
    
    return(pp[sampled_idx])
    
  }, mc.cores = min(length(pp_split), parallel::detectCores()))
  
  # 4. Recombine using the explicit geometry package namespace
  pp_final <- do.call(spatstat.geom::superimpose, sampled_list)
  
  return(pp_final)
}

# Lcross function using geographic sketching
Lcross_sketch <- function(ppp_obj, mark_i, mark_j, mark_column = NULL, n_threads = 1, n_down = 10000, ...) {
  
  # 0. Handle marks safely
  if (is.data.frame(marks(ppp_obj))) {
    if(is.null(mark_column)) {
      stop("marks(ppp_obj) is a dataframe but 'mark_column' is missing. Please provide it (e.g., mark_column = 'merged_type')")
    }
    marks_vec <- marks(ppp_obj)[[mark_column]]
  } else {
    marks_vec <- marks(ppp_obj)
  }

  # 1. Get the indices of the specific cell types dynamically
  idx_i <- which(marks(ppp_obj) == mark_i)
  idx_j <- which(marks(ppp_obj) == mark_j)
  
  # Safety check
  if(length(idx_i) < 2 | length(idx_j) < 2) {
    warning(paste("Not enough points for", mark_i, "and/or", mark_j, "to compute Lcross. Returning NULL."))
    return(NULL)
  }
  
  # 2. Handle Downsampling using the Spatial Sketching approach
  keep_i <- sketch_spatial_R(ppp_obj, idx_i, n_down, n_threads)
  keep_j <- sketch_spatial_R(ppp_obj, idx_j, n_down, n_threads)
  
  # Subset the ppp object
  ppp_sub <- ppp_obj[c(keep_i, keep_j)]
  marks(ppp_sub) <- factor(marks(ppp_sub), levels = c(mark_i, mark_j))
  
  # 4. Calculate Lcross dynamically 
  L <- spatstat.explore::Lcross(ppp_sub, 
                                i = mark_i, 
                                j = mark_j, ...)
  
  return(L)
}

# smoothing function ----
fsmooth <- function(dt, M = 6, genLfun.fd = 4, centrata = F, ...)
{
  # 1. Extract 'r' and the function matrix from the data.table
  # Assuming the first column is the distance vector ('r_value')
  r <- dt$r_value
  
  # Drop the first column (r_value) and convert the rest to a matrix
  mat_fun <- as.matrix(dt[, -1, with = FALSE])
  
  # 2. Apply centering ONLY if requested
  if(centrata == TRUE){
    mat_fun <- sweep(mat_fun, MARGIN = 1, STATS = r, FUN = "-")
  }
  
  # 3. Create Basis and Smooth (requires library(fda))
  K <- length(r) + M - 2
  genbasis <- create.bspline.basis(range(r), K, M, r)
  genfdPar <- fdPar(genbasis, genLfun.fd, lambda = 1e-11)
  
  genfdSmooth <- smooth.basis(r, mat_fun, genfdPar)
  Genfun.fd <- genfdSmooth$fd
  
  # 4. Set Names
  fdn_dist <- paste("Distanza (da 0 a", max(r), "micron)")
  fdnames <- list(fdn_dist,
                  "Sample" = colnames(mat_fun), # Automatically uses the TCGA names from dt
                  "Funzione")
  Genfun.fd$fdnames <- fdnames
  
  return(list(fd = Genfun.fd, fdPar = genfdPar, 
              basis = genbasis, mat = mat_fun, r = r))
}

