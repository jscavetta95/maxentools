% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectBackgroundPoints.R
\name{selectBackgroundPoints}
\alias{selectBackgroundPoints}
\title{Randomly select background points for a given extent.}
\usage{
selectBackgroundPoints(sample_layer, occurrence_data_list, predictors,
  categoricals, num_background_points, use_probabilty, shape,
  output_folder)
}
\arguments{
\item{sample_layer}{A layer with the extent from which to randomly select background points.}

\item{occurrence_data_list}{A list of occurrence point dataframes for each species sampled. Only required for plotting.}

\item{predictors}{A RasterStack of RasterLayer objects representing the environmental predictors.}

\item{categoricals}{A vector of indecies specifying the environmental predictors that contain catagorical or discrete data.}

\item{num_background_points}{The number of background points to select.}

\item{shape}{A shapefile to overlay on the plotted data points map.}

\item{output_folder}{Location to save data and plots.}

\item{use_probablity}{Set TRUE if the background points should be sampled using the bias probablity density. FALSE for uniform sampling.}
}
\value{
Dataframe of selected background coordinates and extracted environmental predictor data.
}
\description{
If provided, background points can be selected from a probablity density layer, such as a bias file. Otherwise,
background points are uniformly selected across the layer extent. Selected background points will be paired with
environmental data for each environmental predictor provided.
}
