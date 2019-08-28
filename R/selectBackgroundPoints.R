#' Randomly select background points for a given extent.
#'
#' If provided, background points can be selected from a probablity density layer, such as a bias file. Otherwise,
#' background points are uniformly selected across the layer extent. Selected background points will be paired with
#' environmental data for each environmental predictor provided.
#'
#' @param sample_layer A layer with the extent from which to randomly select background points.
#' @param occurrence_data_list A list of occurrence point dataframes for each species sampled. This ensures background points will not be selected in the same positions.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param categoricals A vector of indecies specifying the environmental predictors that contain catagorical or discrete data.
#' @param num_background_points The number of background points to select.
#' @param use_probablity Set TRUE if the background points should be sampled using the bias probablity density. FALSE for uniform sampling.
#' @param shape A shapefile to overlay on the plotted data points map.
#' @param output_folder Location to save data and plots.
#' @return Dataframe of selected background coordinates and extracted environmental predictor data.
#' @export
selectBackgroundPoints <- function(sample_layer, occurrence_data_list, predictors, categoricals, num_background_points, use_probabilty, shape, output_folder) {
  combined_occurrence_data <- data.table::rbindlist(occurrence_data_list)[, 1:2]
  background_points <- dismo::randomPoints(sample_layer, num_background_points, combined_occurrence_data, prob = use_probabilty)
  background_data <- na.omit(cbind(background_points, as.data.frame(raster::extract(predictors, background_points))))
  for (i in categoricals) {
    background_data[, i + 2] <- as.factor(background_data[, i + 2])
  }

  colnames(background_data)[1:2] <- c("longitude", "latitude")

  .create_points_map(combined_occurrence_data, background_data[, 1:2], shape, output_folder)
  write.csv(background_data, paste0(output_folder, "background_data.csv"), row.names = FALSE)

  return(background_data)
}

.create_points_map <- function(occurence_data, background_data, shape, output_folder) {
  png(paste0(output_folder, "data_points.png"), width = 1000, height = 800)

  raster::plot(shape, border = "black")
  points(background_data, col = "navy", pch = 19, cex = 0.5)
  points(occurence_data, col = "orangered", pch = 15, cex = 0.5)
  prettymapr::addscalebar()
  prettymapr::addnortharrow("bottomleft", padin = c(1.1, 0.6), scale = 1.25)

  dev.off()
}
