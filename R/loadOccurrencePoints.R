#' Loads the occurrence data from the specified folder.
#'
#' Occurrence coordinates should be in a CSV format with the file extension '.csv'. Files should contain two columns,
#' longitude and latitude. Coordinates will be paired with environmental data for each environmental predictor provided.
#'
#' @param occurrence_folder A path to a directory that contains the occurrence files.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param categoricals A vector of indicies specifying the environmental predictors that contain catagorical or discrete data.
#' @return List of dataframes of occurrence coordinates and extracted environmental predictor data.
#' @export
loadOccurrencePoints <- function(occurrence_folder, predictors, categoricals) {
  species_files <- list.files(occurrence_folder, pattern = "\\.csv$", full.names = TRUE)
  occurrence_points_list <- lapply(species_files, read.csv)
  occurrence_data_list <- lapply(occurrence_points_list, function(occurrence_points) {
    occurrence_data <- na.omit(cbind(occurrence_points, as.data.frame(raster::extract(predictors, occurrence_points))))
    for (i in categoricals) {
      occurrence_data[, i + 2] <- as.factor(occurrence_data[, i + 2])
    }

    return(occurrence_data)
  })

  names(occurrence_data_list) <- lapply(species_files, function(file) {
    name <- sub("\\.csv$", "", basename(file))
    return(name)
  })

  return(occurrence_data_list)
}
