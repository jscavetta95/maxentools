#' Generate subsample values for niche overlap given two species occurrence data.
#'
#' @param occurrence_data1 An occurrence dataframe containing cooridnates and environmental data at each coordinate.
#' @param args1 Arguments, including hyperparameters, to use when training the first pseudo model.
#' @param occurrence_data2 An occurrence dataframe containing cooridnates and environmental data at each coordinate.
#' @param args2 Arguments, including hyperparameters, to use when training the second pseudo model.
#' @param background_data Dataframe of selected background coordinates and extracted environmental predictor data.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param iterations Number of subsample values to generate.
#' @return A list of subsampled niche overlap values.
#' @export
subsampleNicheOverlap <- function(occurrence_data1, args1, occurrence_data2, args2, background_data, predictors, iterations) {
  combined_occurrence_data <- rbind(occurrence_data1, occurrence_data2)

  progress_bar <- dplyr::progress_estimated(iterations)

  require(doSNOW)
  cluster <- parallel::makeCluster(parallel::detectCores())
  doSNOW::registerDoSNOW(cluster)
  opts <- list(progress = function() progress_bar$tick()$print())
  scores <- foreach(i = 1:iterations, .options.snow = opts) %dopar% {
    rows <- sample(nrow(combined_occurrence_data), size = nrow(occurrence_data1), replace = FALSE)

    pseudo1_occurrence_data <- combined_occurrence_data[rows, ]
    pseudo2_occurrence_data <- combined_occurrence_data[-rows, ]

    predictive_map1 <- .runPseudoModel(pseudo1_occurrence_data, background_data, predictors, args1)
    predictive_map2 <- .runPseudoModel(pseudo2_occurrence_data, background_data, predictors, args2)
    overlap <- dismo::nicheOverlap(predictive_map1, predictive_map2)

    return(overlap)
  }

  parallel::stopCluster(cluster)
  progress_bar$stop()

  return(scores)
}

.runPseudoModel <- function(pseudo_occurrence_data, background_data, predictors, args) {
  predictor_data <- rbind(pseudo_occurrence_data, background_data)
  point_class <- c(rep(1, nrow(pseudo_occurrence_data)), rep(0, nrow(background_data)))

  full_model <- dismo::maxent(predictor_data, point_class, args = args, removeDuplicates = TRUE)
  predictive_map <- dismo::predict(full_model, predictors, args = c("outputformat=cloglog"))

  return(predictive_map)
}
