\name{MSCMGARCH-package}
\alias{MSCMGARCH-package}
\alias{MSCMGARCH}
\docType{package}
\title{
\packageTitle{MSCMGARCH}
}
\description{
This package provides a Markov-switching multivariate GARCH model with copula-distributed innovations. It allows for the dynamic assessment of financial market interdependencies, capturing potential regime shifts in asymmetric tail dependence structures among volatility innovations in speculative return series.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MSCMGARCH}
\packageIndices{MSCMGARCH}
This package includes functions for simulating data, estimating model parameters, and forecasting. The most important functions are:
\itemize{
  \item \code{Sim_MSCMGARCH}: Simulates data from the Markov-switching multivariate GARCH model.
  \item \code{MLE_MSCMGARCH}: Estimates the model parameters using maximum likelihood estimation.
  \item \code{rolling_window}: Performs rolling window analysis for forecasting.
}
}
\author{
\packageAuthor{MSCMGARCH}

Maintainer: \packageMaintainer{MSCMGARCH}
}
\references{
Markov-Switching Multivariate GARCH Model with Copula-Distributed Innovations
}
\examples{
# Simulate data
sim_data <- Sim_MSCMGARCH(type = "Normal Gumbel", n = 1000, amount = 10, true_par = c(0.9, 0.2, 12), lower_bound = c(0, 0, 1), upper_bound = c(1, 1, 17), nc = 1, seed = 123)

# Estimate model parameters
estimation <- MLE_MSCMGARCH(r = sim_data, type = "Normal Gumbel", nc = 1)

# Perform rolling window analysis
forecast <- rolling_window(series = sim_data, type = "MS_CMGARCH", copula_type = c(4, 4, 4), asymmetric = FALSE, window_length = 250, portfolio_weights = c(0.5, 0.5), signs = NULL, nc = 1)
}
\keyword{ package }
\seealso{
\code{\link[MSCMGARCH:MSCMGARCH-package]{MSCMGARCH}}
}
