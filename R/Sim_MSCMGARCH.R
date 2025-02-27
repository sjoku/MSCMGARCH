#' Simulate Markov-Switching Multivariate GARCH Model with Copula-Distributed Innovations
#'
#' This function simulates a Markov-Switching Multivariate GARCH model with copula-distributed innovations.
#'
#' @param type A character string specifying the type of copula to be used.
#' @param n An integer specifying the number of observations to be simulated.
#' @param amount An integer specifying the number of simulations to be performed.
#' @param true_par A numeric vector specifying the true parameters of the model.
#' @param lower_bound A numeric vector specifying the lower bounds for the parameters.
#' @param upper_bound A numeric vector specifying the upper bounds for the parameters.
#' @param nc An integer specifying the number of cores to be used for parallel computation.
#' @param cl A cluster object for parallel computation (optional).
#' @param seed An integer specifying the seed for random number generation.
#' @return A list containing the simulated data.
#' @examples
#' \dontrun{
#' true_par <- c(0.5, 0.5, 0.5, 0.5, 0.5)
#' lower_bound <- c(0, 0, 0, 0, 0)
#' upper_bound <- c(1, 1, 1, 1, 1)
#' result <- Sim_MSCMGARCH("Clayton Gumbel", 100, 10, true_par, lower_bound, upper_bound, 2, seed = 123)
#' }
#' @export
Sim_MSCMGARCH <- function(type, n, amount, true_par, lower_bound, upper_bound, nc, cl = NULL, seed) {
  amount_parallel <- rep(1, amount)
  
  f <- function(x) {
    return(Test(type, n, x, true_par, lower_bound, upper_bound, seed, 1))
  }
  
  # for replicability
  theta_list <- future_lapply(X = amount_parallel, FUN = f, future.seed = seed)
  
  list_final <- unlist(theta_list[[1]])
  for (i in 2:length(theta_list)) {
    list_final <- rbind(list_final, unlist(theta_list[[i]]))
  }
  return(list_final)
}
