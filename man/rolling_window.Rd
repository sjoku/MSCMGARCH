\name{rolling_window}
\alias{rolling_window}
\title{Perform Rolling Window Analysis for Forecasting}
\description{
This function performs rolling window analysis for forecasting using the Markov-Switching Multivariate GARCH model with copula-distributed innovations.
}
\usage{
rolling_window(series, type, copula_type, asymmetric, window_length, portfolio_weights, signs, nc)
}
\arguments{
  \item{series}{A matrix of return series.}
  \item{type}{A character string specifying the type of model to be used.}
  \item{copula_type}{A numeric vector specifying the types of copulas to be used.}
  \item{asymmetric}{A logical value indicating whether to use an asymmetric model.}
  \item{window_length}{An integer specifying the length of the rolling window.}
  \item{portfolio_weights}{A numeric vector specifying the portfolio weights.}
  \item{signs}{A matrix of signs for the asymmetric model (optional).}
  \item{nc}{An integer specifying the number of cores to be used for parallel computation.}
}
\value{
A list containing the forecasted values.
}
\examples{
# Load necessary libraries
library(MSCMGARCH)

# Simulate data
sim_data <- Sim_MSCMGARCH(type = "Normal Gumbel", n = 1000, amount = 10, true_par = c(0.9, 0.2, 12), lower_bound = c(0, 0, 1), upper_bound = c(1, 1, 17), nc = 1, seed = 123)

# Perform rolling window analysis
forecast <- rolling_window(series = sim_data, type = "MS_CMGARCH", copula_type = c(4, 4, 4), asymmetric = FALSE, window_length = 250, portfolio_weights = c(0.5, 0.5), signs = NULL, nc = 1)
}
