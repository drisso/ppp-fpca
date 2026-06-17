library(fda)

fsmooth <- function(fcross_list, fun, M = 6, genLfun.fd = 4, centrata = F, ...)
{
  
  # 1. Calculate raw L-functions (Common for both)
  r <- fcross_list[[1]]$r
  mat_fun <- sapply(fcross_list, function(x) x$border)
  
  # 2. Apply centering ONLY if requested
  if(centrata == TRUE){
    mat_fun <- sweep(mat_fun, MARGIN = 1, STATS = r, FUN = "-")
  }
  
  # 3. Create Basis and Smooth
  K <- length(r) + M - 2
  genbasis <- create.bspline.basis(range(r), K, M, r)
  genfdPar <- fdPar(genbasis, genLfun.fd, lambda = 1e-11)
  
  genfdSmooth <- smooth.basis(r, mat_fun, genfdPar)
  Genfun.fd <- genfdSmooth$fd
  
  # 4. Set Names
  fdn_dist <- paste("Distanza (da 0 a", max(r), "micron)")
  fdnames <- list(fdn_dist,
                  "Gene" = colnames(mat_fun), # Crucial: use names of surviving columns
                  "Funzione")
  Genfun.fd$fdnames <- fdnames
  
  return(list(fd = Genfun.fd, fdPar = genfdPar, 
              basis = genbasis, mat = mat_fun, r = r))
}