#' Generate subsample Moran's I values for niche overlap given two species occurrence data.
#'
#' @param species_models A list of species lists where each contains a path, predictive_map, occurrence_data, and args element.
#' @param background_data Dataframe of selected background coordinates and extracted environmental predictor data.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param iterations Number of random replicates to run for each species combination.
#' @param cores Number of cores to use. Defaults to a single core; specify more for multiprocessing.
#' @return A dataframe containing the species combination, the real Moran's I, all replicate Moran's I values, and a P value.
#' @export
calculateNicheOverlapProbability <- function(species_models, background_data, predictors, iterations, cores=1) {
  comparisons <- combn(species_models, 2)
  results <- setNames(data.frame(matrix(ncol = 4+iterations, nrow = ncol(comparisons))),
                      c("species_1", "species_2", "real_I", paste0("pseudo_I_", seq(1,iterations)), "p_value"))

  foreach(i = 1:ncol(comparisons)) %do% {
    species1 <- comparisons[[1, i]]
    species2 <- comparisons[[2, i]]

    results$species_1[i] <- species1$path
    results$species_2[i] <- species2$path

    results$real_I[i] <- dismo::nicheOverlap(species1$predictive_map, species2$predictive_map)

    pseudo_scores <- subsampleNicheOverlap(species1$occurrence_data[, -(1:2)], species1$args,
                                           species2$occurrence_data[, -(1:2)], species2$args,
                                           background_data[, -(1:2)], predictors, iterations,
                                           cores)
    j <- 4
    for(pseudo in pseudo_scores)
    {
      results[i,j] <- pseudo
      j = j + 1
    }

    results$p_value[i] <- 1 - mean(pseudo_scores > results$real_I[i])
  }

  return(results)
}

#' Generate subsample values for niche overlap given two species occurrence data.
#'
#' @param occurrence_data1 An occurrence dataframe containing cooridnates and environmental data at each coordinate.
#' @param args1 Arguments, including hyperparameters, to use when training the first pseudo model.
#' @param occurrence_data2 An occurrence dataframe containing cooridnates and environmental data at each coordinate.
#' @param args2 Arguments, including hyperparameters, to use when training the second pseudo model.
#' @param background_data Dataframe of selected background coordinates and extracted environmental predictor data.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param iterations Number of subsample values to generate.
#' @param cores Number of cores to use.
#' @return A list of subsampled niche overlap values.
#' @export
subsampleNicheOverlap <- function(occurrence_data1, args1, occurrence_data2, args2, background_data, predictors, iterations, cores) {
  combined_occurrence_data <- rbind(occurrence_data1, occurrence_data2)

  progress_bar <- dplyr::progress_estimated(iterations)

  require(doSNOW)
  cluster <- parallel::makeCluster(cores)
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
