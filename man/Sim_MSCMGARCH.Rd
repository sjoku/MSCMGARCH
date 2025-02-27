\name{Sim_MSCMGARCH}
\alias{Sim_MSCMGARCH}
\title{Simulate Markov-Switching Multivariate GARCH Model with Copula-Distributed Innovations}
\description{
This function simulates a Markov-Switching Multivariate GARCH model with copula-distributed innovations.
}
\usage{
Sim_MSCMGARCH(type, n, amount, true_par, lower_bound, upper_bound, nc, cl = NULL, seed)
}
\arguments{
  \item{type}{A character string specifying the type of copula to be used.}
  \item{n}{An integer specifying the number of observations to be simulated.}
  \item{amount}{An integer specifying the number of simulations to be performed.}
  \item{true_par}{A numeric vector specifying the true parameters of the model.}
  \item{lower_bound}{A numeric vector specifying the lower bounds for the parameters.}
  \item{upper_bound}{A numeric vector specifying the upper bounds for the parameters.}
  \item{nc}{An integer specifying the number of cores to be used for parallel computation.}
  \item{cl}{A cluster object for parallel computation (optional).}
  \item{seed}{An integer specifying the seed for random number generation.}
}
\value{
A list containing the simulated data.
}
\examples{
\dontrun{
true_par <- c(0.5, 0.5, 0.5, 0.5, 0.5)
lower_bound <- c(0, 0, 0, 0, 0)
upper_bound <- c(1, 1, 1, 1, 1)
result <- Sim_MSCMGARCH("Clayton Gumbel", 100, 10, true_par, lower_bound, upper_bound, 2, seed = 123)
}
}
