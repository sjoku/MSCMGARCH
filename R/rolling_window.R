#' Perform Rolling Window Analysis for Forecasting
#'
#' This function performs rolling window analysis for forecasting using the Markov-Switching Multivariate GARCH model with copula-distributed innovations.
#'
#' @param series A matrix of return series.
#' @param type A character string specifying the type of model to be used.
#' @param copula_type A numeric vector specifying the types of copulas to be used.
#' @param asymmetric A logical value indicating whether to use an asymmetric model.
#' @param window_length An integer specifying the length of the rolling window.
#' @param portfolio_weights A numeric vector specifying the portfolio weights.
#' @param signs A matrix of signs for the asymmetric model (optional).
#' @param nc An integer specifying the number of cores to be used for parallel computation.
#' @return A list containing the forecasted values.
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(MSCMGARCH)
#'
#' # Simulate data
#' sim_data <- Sim_MSCMGARCH(type = "Normal Gumbel", n = 1000, amount = 10, true_par = c(0.9, 0.2, 12), lower_bound = c(0, 0, 1), upper_bound = c(1, 1, 17), nc = 1, seed = 123)
#'
#' # Perform rolling window analysis
#' forecast <- rolling_window(series = sim_data, type = "MS_CMGARCH", copula_type = c(4, 4, 4), asymmetric = FALSE, window_length = 250, portfolio_weights = c(0.5, 0.5), signs = NULL, nc = 1)
#' }
#' @export
rolling_window <- function(series, type, copula_type, asymmetric, window_length, portfolio_weights, signs, nc) {
  # Initialize an empty list to store the forecasted values
  forecasted_values <- list()
  
  # Loop through the series using a rolling window approach
  for (i in 1:(nrow(series) - window_length + 1)) {
    # Define the rolling window
    window <- series[i:(i + window_length - 1), ]
    
    # Estimate the model parameters using the MLE_MSCMGARCH function
    if (asymmetric) {
      estimation <- MLE_MSCMGARCH(r = window, type = type, copula_type = copula_type, asymmetric = asymmetric, signs = signs, nc = nc)
    } else {
      estimation <- MLE_MSCMGARCH(r = window, type = type, copula_type = copula_type, nc = nc)
    }
    
    # Perform the forecasting and store the forecasted values in the list
    forecasted_values[[i]] <- estimation
  }
  
  return(forecasted_values)
}
