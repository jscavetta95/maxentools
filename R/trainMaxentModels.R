#' Train maxent models with the provided occurence data, background data, environmental predictors, and hyperparameter ranges.
#'
#' Train models with varying factor combinations and beta regularization values. Each factor combination provided will be
#' separately trained the desired amount of beta regularization values. Beta regularization values are found using the LIPO
#' search technique, which will generally tend towards the best hyperparameter value for the test case. Many metrics are calculated
#' and provided, however, currently the symmetric extremal dependence index (SEDI) metric is used to determine the best model.
#' Total number of models trained = factors * regularization iterations.
#'
#' @param occurrence_data An occurrence dataframe containing cooridnates and environmental data at each coordinate.
#' @param background_data Dataframe of selected background coordinates and extracted environmental predictor data.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param factors The combination of factors to use. Should be structured like "LQHPT".
#' L = linear, Q = quadratic, H = hinge, P = product, T = threshold.
#' @param regularization_iterations The number of regularization values to try for each factor value.
#' Essentially the number of models to run for each factor.
#' @param regularization_range A 2 element vector containing the start and end of the regularization values
#' to consider in the search.
#' @param shape A shapefile to overlay on produced predictive maps.
#' @param output_folder The path where data and plots will be saved.
#' @return A list containing the best model, the occurrence data used, the arguments (hyperparameters) used,
#' the predictive map generated from the best model, and the path where data was saved.
#' @export
trainMaxentModels <- function(occurrence_data, background_data, predictors, factors,
                              regularization_iterations, regularization_range, shape, output_folder) {
  group_data <- ENMeval::get.block(occurrence_data[, 1:2], background_data[, 1:2])

  require(doSNOW)
  cluster <- parallel::makeCluster(parallel::detectCores())
  doSNOW::registerDoSNOW(cluster)
  out <- foreach(i = 1:length(factors)) %dopar% {
    tuneRegularization(
      occurrence_data[, -(1:2)], background_data[, -(1:2)],
      predictors, group_data, factors[i], regularization_range, regularization_iterations
    )
  }

  parallel::stopCluster(cluster)

  full_models <- unlist(lapply(out, function(x) lapply(x, function(x) x$full_model)), recursive = FALSE)
  stats_list <- unlist(lapply(out, function(x) lapply(x, function(x) x$stats)), recursive = FALSE)
  factor_list <- unlist(lapply(out, function(x) lapply(x, function(x) x$parameters$factor)), recursive = FALSE)
  regularization_list <- unlist(lapply(out, function(x) lapply(x, function(x) x$parameters$regularization)), recursive = FALSE)

  stats <- data.table::rbindlist(stats_list)
  stats <- cbind(factor = unlist(factor_list), regularization = unlist(regularization_list), stats)
  stats <- .rankModels(stats)

  best_index <- order(stats$mean_symmetric_extremal_dependence_index, decreasing = FALSE)[1]
  best_model <- full_models[[best_index]]
  args <- ENMeval::make.args(stats$regularization[best_index], stats$factor[best_index])[[1]]

  predictive_map <- dismo::predict(best_model, predictors, args = c("outputformat=cloglog"))
  .plotModel(predictive_map, shape, output_folder)

  predictor_out <- .getTestGainAndResponse(occurrence_data[, -(1:2)], background_data[, -(1:2)], args)

  save(full_models, file = paste0(output_folder, "/models.RData"))
  write.csv(stats, paste0(output_folder, "/stats.csv"), row.names = FALSE)
  write.csv(predictor_out$test_gain, paste0(output_folder, "/test_gain.csv"))
  write.csv(predictor_out$response_curve, paste0(output_folder, "/response_curve.csv"), row.names = FALSE)

  return(list(
    best_model = best_model, occurrence_data = occurrence_data,
    args = args, predictive_map = predictive_map, path = output_folder
  ))
}

.getTestGainAndResponse <- function(occurrence_data, background_data, args) {
  predictor_data <- rbind(occurrence_data, background_data)
  point_class <- c(rep(1, nrow(occurrence_data)), rep(0, nrow(background_data)))
  full_model <- dismo::maxent(predictor_data, point_class,
    args = c("jackknife=true", "responsecurves=true", "replicates=10", "writeplotdata=true", args),
    removeDuplicates = TRUE
  )

  test_gain <- full_model@results[grepl("Test.gain(.with.only.*)?$", rownames(full_model@results)), "species (average)"]
  most_important_predictor <- tail(strsplit(names(which.max(test_gain[-1])), "\\.")[[1]], n = 1)

  path <- gsub("maxent.html", "", full_model@html)
  response_curve <- read.csv(paste0(path, "/plots/species_", most_important_predictor, "_only.dat"))

  return(list(test_gain = test_gain, response_curve = response_curve))
}

.rankModels <- function(stats) {
  ranks <- matrix(nrow = nrow(stats), ncol = 8)
  ranks[order(stats$mean_test_auc, decreasing = TRUE), 1] <- 1:nrow(stats)
  ranks[order(stats$mean_omission_rate, decreasing = FALSE), 2] <- 1:nrow(stats)
  ranks[order(stats$mean_kappa, decreasing = TRUE), 3] <- 1:nrow(stats)
  ranks[order(stats$mean_true_skill_statistic, decreasing = TRUE), 4] <- 1:nrow(stats)
  ranks[order(stats$mean_odds_ratio_skill_score, decreasing = TRUE), 5] <- 1:nrow(stats)
  ranks[order(stats$mean_symmetric_extremal_dependence_index, decreasing = TRUE), 6] <- 1:nrow(stats)
  ranks[order(stats$mean_f1, decreasing = TRUE), 7] <- 1:nrow(stats)
  ranks[order(stats$mean_f2, decreasing = TRUE), 8] <- 1:nrow(stats)

  stats$rank <- apply(ranks, 1, mean)

  return(stats)
}

.plotModel <- function(predictive_map, shape, output_folder) {
  png(paste0(output_folder, "/best_model.png"), width = 171, height = 150, units = "mm", res = 600)

  raster::plot(predictive_map, col = viridis::viridis(256, direction = -1), axes = FALSE, frame.plot = FALSE)
  box(col = "white")
  raster::plot(shape, border = "black", add = TRUE)
  prettymapr::addscalebar()
  prettymapr::addnortharrow("bottomleft", padin = c(0.4, 0.4), scale = 1)

  dev.off()
}
