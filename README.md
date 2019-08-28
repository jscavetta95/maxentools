# maxentools
A set of high level functions to provide easy training and evaluation of Maxent species distribution models.

Install with ```devtools::install_github("jscavetta95/maxentools")```

## Example
### Using maxentools to train, evaluate, and produce a species distribtuion model. Also includes niche overlap evaluation.

The data for this example is available in the ./data/ folder. Please be sure to unzip the .asc files in the layers folder.

First, load in the required data. This includes a shape file and environmental layers.
```{r load, eval = FALSE}
library(maxentools)

shape <- rgdal::readOGR(dsn = "./data/states.gdb")
shape = raster::subset(shape, STATE_NAME != "District of Columbia" & STATE_NAME != "Alaska" & STATE_NAME != "Hawaii")
shape <- sp::spTransform(shape, sp::CRS('+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'))

predictors <- raster::stack(list.files("./data/layers", pattern = '\\.asc$', full.names = TRUE))
```

Next, we want to load in our occurrence coordinates for each species. We also want the coordinates to be matched with the environmental data at that point. To do this, we can use the loadOccurrencePoints function. Note, we must specify which if any of our environmental layers are categorical.
```{r occurrence, eval = FALSE}
occurrence_data_list <- loadOccurrencePoints(occurrence_folder = "./data/occurrences/", predictors = predictors, categoricals = c(1))
```

Now, we need to have background (pseudo-absence) points for our models. However, we wish to control for potential sampling bias in our data, so we generate a sampling bias map. This will be a layer of probablities for sampling a background point, such that background points have a higher chance of being sampled near occurrence points. We can use all of our species occurrences to generate a the bias map.
```{r bias, eval = FALSE}
bias_density <- createBiasDensity(predictors = predictors, occurrence_points_list = occurrence_data_list, 
                                  scale = c(1,20), shape = shape, output_folder = "./data/")
```

With the bias layer, we can now sample background points. Let's sample 10000 background points using the bias probablities.
```{r background, eval = FALSE}
background_data <- selectBackgroundPoints(sample_layer = bias_density, occurrence_data_list = occurrence_data_list, predictors = predictors,
                                          categoricals = c(1), num_background_points = 10000, use_probabilty = TRUE, 
                                          shape = shape, output_folder = "./data/")
```

We must choose our hyperparameter values we wish to test. There are two to consider, the allowed factor combinations (linear, quadratic, etc.) and the beta regularization multiplier. The total models that will be ran is regularization iterations * length of factors.
```{r hyperparameters, eval = FALSE}
factors <- c("L", "LQ", "LQH", "LQHP", "LQHPT")
regularization_iterations <- 10
```

Before we run the models, let's create the folders we wish to store data in.
```{r folders, eval = FALSE}
lapply(names(occurrence_data_list), function(name) if(!dir.exists(paste0("./data/", name))) dir.create(paste0("./data/", name)))
```

Let's run our models. We will use parallel processing with a progress bar to monitor progress. This may take a while.
```{r run, eval = FALSE}
progress_bar <- dplyr::progress_estimated(length(occurrence_data_list))

require(doSNOW)
cluster <- parallel::makeCluster(parallel::detectCores() / 2)
doSNOW::registerDoSNOW(cluster)
opts <- list(progress = function() progress_bar$tick()$print())

best_species_models <- foreach(i = seq_along(occurrence_data_list), .packages = c("maxentools"), .options.snow = opts) %dopar% {
  return(trainMaxentModels(occurrence_data = occurrence_data_list[[i]], background_data = background_data, predictors = predictors,
                           factors = factors, regularization_iterations = regularization_iterations, regularization_range = c(0, 2.5), 
                           shape = shape, output_folder = paste0("./data/", names(occurrence_data_list)[[i]])))
}

parallel::stopCluster(cluster)
progress_bar$stop()
```

Our data should be saved in folders for each species. This includes the predictive map output and all the metrics for each model.

As an extra evaluation, we can use a subsampling technique to determine if two species significantly overlap. We will reduce our environmental resolution by a factor of 10 to save time.
```{r overlap, eval = FALSE}
predictors_aggregate <- raster::aggregate(predictors, fact=10)

results <- calculateNicheOverlapProbability(species_models = best_species_models, background_data = background_data, 
                                            predictors = predictors, iterations = 100)

write.csv(results, "./data/niche_overlap.csv", row.names = FALSE)
```
