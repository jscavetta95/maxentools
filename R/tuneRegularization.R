#' Use the AdaLipo technique to search for the best value for the regularization hyperparameter.
#'
#' @param occurrence_data Dataframe of species occurrence coordinates and extracted environmental predictor data.
#' @param background_data Dataframe of selected background coordinates and extracted environmental predictor data.
#' @param predictors A RasterStack of RasterLayer objects representing the environmental predictors.
#' @param group_data A named list of bin designation for occurrence points and background points.
#' @param factor The factor hyperparameter value to use in model creation.
#' @param range A 2 element vector containing the start and end point of the range of regularization values to search.
#' @param num_iterations Number of models to create.
#' @return A list of results from all created models.
#' @export
tuneRegularization <- function(occurrence_data, background_data, predictors, group_data, factor, range, num_iterations) {
  combined_results <- vector("list", num_iterations)

  candidate <- runif(1, range[1], range[2])
  results <- .createModel(occurrence_data, background_data, predictors, group_data, factor, candidate)
  combined_results[[1]] <- results
  score <- .scoreResults(results)
  selected_values <- matrix(c(candidate, score), ncol = 2)

  iteration <- 1
  lipschitz_constant <- 0
  while (iteration < num_iterations) {
    if (sample(c(0, 1), 1)) {
      candidate <- runif(1, range[1], range[2])
    } else {
      candidate_not_found <- TRUE
      tries <- 0
      while (candidate_not_found & tries < 1000) {
        candidate <- runif(1, range[1], range[2])

        min_upper_bound <- min(apply(selected_values, 1, function(param) {
          param[2] + lipschitz_constant * sqrt(sum((candidate - param[1])^2))
        }))

        if (min_upper_bound >= max(selected_values[, 2])) {
          candidate_not_found <- FALSE
        } else {
          tries <- tries + 1
        }
      }
    }

    results <- .createModel(occurrence_data, background_data, predictors, group_data, factor, candidate)
    combined_results[[iteration + 1]] <- results
    score <- .scoreResults(results)
    selected_values <- rbind(selected_values, matrix(c(candidate, score), ncol = 2))
    iteration <- iteration + 1

    slopes <- combn(
      nrow(selected_values), 2,
      function(combo) {
        abs(selected_values[combo[1], 2] - selected_values[combo[2], 2]) /
          sqrt(sum((selected_values[combo[1], 1] - selected_values[combo[2], 1])^2))
      }
    )

    lipschitz_constant <- ceiling(max(slopes))
  }

  return(combined_results)
}

.scoreResults <- function(results) {
  return(mean(results$stats$mean_symmetric_extremal_dependence_index))
}

.createModel <- function(occurrence_data, background_data, predictors, group_data, factor, regularization) {
  predictor_data <- rbind(occurrence_data, background_data)
  point_class <- c(rep(1, nrow(occurrence_data)), rep(0, nrow(background_data)))

  args <- ENMeval::make.args(regularization, factor)[[1]]

  full_model <- dismo::maxent(predictor_data, point_class, args = args, removeDuplicates = TRUE)
  # predictive_map <- dismo::predict(full_model, predictors, args = c("outputformat=cloglog"))

  results <- setNames(
    data.frame(matrix(ncol = 10, nrow = 4)),
    c(
      "train_auc", "test_auc", "auc_diff", "omission_rate", "kappa",
      "true_skill_statistic", "odds_ratio_skill_score",
      "symmetric_extremal_dependence_index", "f1", "f2"
    )
  )

  for (bin in 1:4) {
    train_occurrence_data <- occurrence_data[group_data$occ.grp != bin, ]
    test_occurrence_data <- occurrence_data[group_data$occ.grp == bin, ]
    bin_background_data <- background_data[group_data$bg.grp != bin, ]

    predictor_data <- rbind(train_occurrence_data, bin_background_data)
    point_class <- c(rep(1, nrow(train_occurrence_data)), rep(0, nrow(bin_background_data)))

    model <- dismo::maxent(predictor_data, point_class, args = args, removeDuplicates = TRUE)

    train_eval <- dismo::evaluate(train_occurrence_data, bin_background_data, model)
    test_eval <- dismo::evaluate(test_occurrence_data, bin_background_data, model)
    train_prediction <- dismo::predict(model, train_occurrence_data)
    test_prediction <- dismo::predict(model, test_occurrence_data)

    results <- .calculateMetrics(train_eval, test_eval, train_prediction, test_prediction, results, bin)
  }

  stats <- data.frame(
    mean_train_auc = mean(results$train_auc), var_train_auc = var(results$train_auc),
    mean_test_auc = mean(results$test_auc), var_test_auc = var(results$test_auc),
    mean_auc_diff = mean(results$auc_diff), var_auc_diff = var(results$auc_diff),
    mean_omission_rate = mean(results$omission_rate), var_omission_rate = var(results$omission_rate),
    mean_kappa = mean(results$kappa), var_kappa = var(results$kappa),
    mean_true_skill_statistic = mean(results$true_skill_statistic), var_true_skill_statistic = var(results$true_skill_statistic),
    mean_odds_ratio_skill_score = mean(results$odds_ratio_skill_score), var_odds_ratio_skill_score = var(results$odds_ratio_skill_score),
    mean_symmetric_extremal_dependence_index = mean(results$symmetric_extremal_dependence_index), var_symmetric_extremal_dependence_index = var(results$symmetric_extremal_dependence_index),
    mean_f1 = mean(results$f1), var_f1 = var(results$f1),
    mean_f2 = mean(results$f2), var_f2 = var(results$f2)
  )

  return(list(full_model = full_model, stats = stats, parameters = list(factor = factor, regularization = regularization)))
}

.calculateMetrics <- function(train_eval, test_eval, train_prediction, test_prediction, results, bin) {
  results$train_auc[bin] <- train_eval@auc
  results$test_auc[bin] <- test_eval@auc
  results$auc_diff[bin] <- train_eval@auc - test_eval@auc
  results$omission_rate[bin] <- mean(test_prediction < min(train_prediction))
  results$kappa[bin] <- max(test_eval@kappa, na.rm = TRUE)
  results$true_skill_statistic[bin] <- max(test_eval@TPR - test_eval@FPR, na.rm = TRUE)

  true_positives <- test_eval@confusion[, 1]
  false_positives <- test_eval@confusion[, 2]
  false_negatives <- test_eval@confusion[, 3]
  true_negatives <- test_eval@confusion[, 4]

  results$odds_ratio_skill_score[bin] <- max(0,
    (true_positives * true_negatives - false_positives * false_negatives) /
      (true_positives * true_negatives + false_positives * false_negatives),
    na.rm = TRUE
  )

  results$symmetric_extremal_dependence_index[bin] <- max(apply(test_eval@confusion, 1, .calculateSEDI), na.rm = TRUE)

  results$f1[bin] <- .calculateFScore(true_positives, false_negatives, false_positives, 1)
  results$f2[bin] <- .calculateFScore(true_positives, false_negatives, false_positives, 2)

  return(results)
}

.calculateSEDI <- function(confusion) {
  true_positive <- confusion[1]
  false_positive <- confusion[2]
  false_negative <- confusion[3]
  true_negative <- confusion[4]

  if (false_positive + false_negative == 0) {
    return(1)
  }

  if (true_positive + true_negative == 0) {
    return(-1)
  }

  if (true_positive + false_positive == 0 || true_negative + false_negative == 0) {
    return(0)
  }

  sedi.tp <- .replaceZero(true_positive)
  sedi.fp <- .replaceZero(false_positive)
  sedi.fn <- .replaceZero(false_negative)
  sedi.tn <- .replaceZero(true_negative)

  sedi.fpr <- sedi.fp / (sedi.fp + sedi.tn)
  sedi.tpr <- sedi.tp / (sedi.tp + sedi.fn)

  return((log(sedi.fpr) - log(sedi.tpr) - log(1 - sedi.fpr) + log(1 - sedi.tpr)) /
    (log(sedi.fpr) + log(sedi.tpr) + log(1 - sedi.fpr) + log(1 - sedi.tpr)))
}

.replaceZero <- function(value) {
  return(ifelse(value > 0, value, 1e-05))
}

.calculateFScore <- function(true_positive, false_negative, false_positive, beta) {
  f_beta <- (1 + beta^2) * true_positive / ((1 + beta^2) * true_positive + beta^2 * false_negative + false_positive)
  return(max(f_beta, na.rm = TRUE))
}
