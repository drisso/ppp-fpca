library(ggplot2)
library(dplyr)
library(tidyr)

# Lfuns of selected samples
plot_fv_list <- function(list_of_funcs, correction = "border") {
  
  # 1. Extract unique mapping of unique_id to sample_type from the SPE object
  # meta_map <- as.data.frame(colData(spe)) |> 
  #   dplyr::select(unique_id, sample_type) |> 
  #   dplyr::distinct()
  # 
  # 2. Run the map_dfr block dynamically on whatever list you provide
  df_lcross <- map_dfr(names(list_of_funcs), function(slide_name) {
    
    # Extract the specific L-function object from the provided list
    L_obj <- list_of_funcs[[slide_name]]
    
    as.data.frame(L_obj) |> 
      mutate(
        Slide = slide_name,
        # Use .data[[correction]] to dynamically grab the 'border' column (or 'iso', etc.)
        L_r_minus_r = .data[[correction]] - r,  
        
        # Extract the marks names exactly as you designed
        mark_i = trimws(strsplit(gsub("list\\(|\\)", "", attr(L_obj, "fname")[2]), ",")[[1]][1]),
        mark_j = trimws(strsplit(gsub("list\\(|\\)", "", attr(L_obj, "fname")[2]), ",")[[1]][2])
      )
  })
  
  # 3. Join the sample_type metadata to our plot dataframe
  # df_lcross <- df_lcross |> 
  #   left_join(meta_map, by = c("Slide" = "unique_id"))
  # 
  # 4. Plot with color mapped to sample_type
  p <- ggplot(df_lcross, aes(x = r, y = L_r_minus_r)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(alpha = 0.7, linewidth = 1.5) + 
    theme_minimal() +
    labs(
      title = "",
      subtitle = "",
      x = "Radius (r)",
      y = "L(r) - r",
      color = "Sample Type"
    ) +
    theme(legend.position = "right",
          legend.title = element_text(size = 36, face = "bold"),
          legend.text  = element_text(size = 26),
          axis.title = element_text(size = 26),
          axis.text  = element_text(size = 18),
    ) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1, linewidth = 4), ncol = 1))
  
  return(p)
}

# plot points

# spatial plot of points
points_plot <- \(ppp) {
  ggplot() +
    geom_raster(
      data = as.data.frame(Window(ppp)),
      aes(x = x, y = y),
      fill = "grey92",
      alpha = 1
    ) +
    geom_point(
      data = as.data.frame(ppp),
      aes(x = x, y = y, color = marks),
      size = 0.05,
      alpha = 0.75
    ) +
    coord_fixed(expand = FALSE) +
    scale_x_reverse() +
    scale_y_reverse() +
    #scale_colour_paletteer_d("colorblindr::OkabeIto", name = "Cell types") +
    # scale_color_manual(values = my_colors,
    #                    name = "Cell types") +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      #plot.title = element_text(size = 22, face = "bold"),
      #plot.subtitle = element_text(size = 14),
      #axis.title = element_text(size = 16),
      #axis.text  = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 36, face = "bold"),
      legend.text  = element_text(size = 26),
      legend.key.height = unit(10, "mm"),
      legend.key.width  = unit(10, "mm"),
      plot.margin = margin(8, 12, 8, 8, "mm"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    labs(
      title = "",
      subtitle = "",
      x = "X coordinate",
      y = "Y coordinate"
    )
}


# plot functional Lfun
plot_smoothed_Lfun <- function(smoothed_obj, clusters_vec, mark_i = "Mark 1", mark_j = "Mark 2", center_plot = TRUE) {
  
  # 1. Extract the distance vector (r)
  r_vec <- smoothed_obj$r
  
  # 2. Evaluate the functional data object (fd) to get smoothed values
  # This creates a matrix where rows are 'r' distances and columns are samples
  smooth_mat <- eval.fd(r_vec, smoothed_obj$fd)
  
  # 3. Convert the matrix to a data.frame and prepare for ggplot
  df_smooth <- as.data.frame(smooth_mat)
  df_smooth$r <- r_vec
  
  # Pivot to long format
  df_long <- pivot_longer(
    df_smooth,
    cols = -r,
    names_to = "Slide",
    values_to = "L_value"
  )
  
  # 4. Center the L-function (L(r) - r) if requested 
  # (Since you ran fsmooth with centrata = FALSE, we can center it visually here)
  if (center_plot) {
    df_long <- df_long |> mutate(L_plot_value = L_value - r)
  } else {
    df_long <- df_long |> mutate(L_plot_value = L_value)
  }
  
  # 5. Prepare the cluster metadata from the named vector (memb_hc)
  cluster_df <- data.frame(
    Slide = names(clusters_vec),
    cluster = as.factor(clusters_vec)
  )
  
  # 6. Join the functional data with the cluster assignments
  df_final <- df_long |> 
    left_join(cluster_df, by = "Slide")
  
  # 7. Build the ggplot
  p <- ggplot(df_final, aes(x = r, y = L_plot_value, group = Slide, color = cluster)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(alpha = 0.7, linewidth = 0.6) + 
    theme_minimal() +
    labs(
      title = paste0("Smoothed L-cross Function: ", mark_i, " vs ", mark_j),
      subtitle = ifelse(center_plot, 
                        "Centered L-function (L(r) - r). Positive = Attraction, Negative = Repulsion", 
                        "L-function L(r)"),
      x = "Radius (r)",
      y = ifelse(center_plot, "L(r) - r", "L(r)"),
      color = "Cluster"
    ) +
    theme(legend.position = "right")
  
  return(p)
}
