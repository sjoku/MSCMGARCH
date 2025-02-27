\name{MLE_MSCMGARCH}
\alias{MLE_MSCMGARCH}
\title{Estimate Markov-Switching Multivariate GARCH Model with Copula-Distributed Innovations}
\description{
This function estimates the parameters of a Markov-Switching Multivariate GARCH model with copula-distributed innovations using maximum likelihood estimation.
}
\usage{
MLE_MSCMGARCH(r, type, start_val = NULL, BEKK1 = NULL, nc)
}
\arguments{
  \item{r}{A matrix of return series.}
  \item{type}{A character string specifying the type of copula to be used.}
  \item{start_val}{Initial values for the parameters (optional).}
  \item{BEKK1}{Initial values for the BEKK parameters (optional).}
  \item{nc}{An integer specifying the number of cores to be used for parallel computation.}
}
\value{
A list containing the BEKK parameters and the estimated parameters.
}
\examples{
# Load necessary libraries
library(MSCMGARCH)

# Simulate data
sim_data <- Sim_MSCMGARCH(type = "Normal Gumbel", n = 1000, amount = 10, true_par = c(0.9, 0.2, 12), lower_bound = c(0, 0, 1), upper_bound = c(1, 1, 17), nc = 1, seed = 123)

# Estimate model parameters
estimation <- MLE_MSCMGARCH(r = sim_data, type = "Normal Gumbel", nc = 1)
}
