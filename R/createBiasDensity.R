#' Create a probability density map based on the occurence points of all species.
#'
#' The creation of a probablity density map (bias file) helps avoid sampling bias in models. Restricting background point selection
#' to areas around occurence points decreases the potential of unbalanced sampling between occurence points and background points. This
#' may occur if not all areas in the training extent were uniformly randomly sampled.
#'
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param occurrence_points_list A list of occurence point dataframes for each species sampled.
#' @param scale A 2 element vector containing the minimum and maximum range for the bias probablity scale.
#' @param shape A shapefile to overlay on the plotted bias density map.
#' @param output_folder Location to save data and plots.
#' @return A RasterLayer of selection probablities ranging from 0 to 100 throughout the training extent.
#' @export
createBiasDensity <- function(predictors, occurrence_points_list, scale, shape, output_folder) {
  combined_occurrence_points <- data.table::rbindlist(occurrence_points_list)[, 1:2]
  occurrence_raster <- raster::rasterize(combined_occurrence_points, predictors, 1)
  occurrence_coordinates <- raster::coordinates(occurrence_raster)[which(occurrence_raster@data@values == 1), ]

  kernel_density <- MASS::kde2d(occurrence_coordinates[, 1], occurrence_coordinates[, 2],
    n = c(nrow(occurrence_raster), ncol(occurrence_raster)),
    lims = c(
      raster::extent(predictors)@xmin, raster::extent(predictors)@xmax,
      raster::extent(predictors)@ymin, raster::extent(predictors)@ymax
    )
  )

  kernel_density$z <- .scaleValues(kernel_density$z, scale[1], scale[2])
  bias_density <- raster::resample(raster::raster(kernel_density), predictors)
  bias_density@data@values[which(rowSums(is.na(raster::values(predictors))) > 0)] <- NA

  raster::writeRaster(bias_density, paste0(output_folder, "biasfile.asc"), overwrite = TRUE)
  .create_bias_map(bias_density, shape, output_folder)

  return(bias_density)
}

.scaleValues <- function(values, new_min, new_max) {
  value_max <- max(values)
  value_min <- min(values)
  value_range <- value_max - value_min

  new_range <- new_max - new_min

  scaled_values <- apply(values, c(1, 2), function(x) {
    new_range * ((x - value_min) / value_range) + new_min
  })

  return(scaled_values)
}

.create_bias_map <- function(bias_density, shape, output_folder) {
  png(paste0(output_folder, "bias_density.png"), width = 1000, height = 800)

  raster::plot(bias_density, col = viridis::viridis(256, direction = -1), axes = FALSE, frame.plot = FALSE)
  box(col = "white")
  raster::plot(shape, border = "black", add = TRUE)
  prettymapr::addscalebar()
  prettymapr::addnortharrow("bottomleft", padin = c(1.1, 0.6), scale = 1.25)

  dev.off()
}
