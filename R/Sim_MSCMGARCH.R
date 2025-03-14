#' Simulate Markov-Switching Copula Multivariate GARCH Models
#'
#' Performs Monte Carlo simulations for evaluating parameter estimation of Markov-Switching Copula Multivariate GARCH (MSCMGARCH) models. The function leverages parallel processing to efficiently run multiple simulation replications, each providing a set of estimated parameters and their respective accuracy metrics.
#'
#' @param type Character vector specifying the copula types to use in the simulation. Implemented copula types are:
#'   - "Clayton Gumbel"
#'   - "Clayton Gumbel BEKK"
#'   - "Normal Gumbel"
#'   - "Normal Gumbel BEKK"
#'   - "Clayton Gumbel Survival"
#'   - "Clayton Gumbel Survival BEKK"
#' @param n Integer. The sample size for each simulation replication.
#' @param amount Integer. The number of simulation replications to perform.
#' @param true_par Numeric vector. The true parameter values used for generating simulated data.
#' @param lower_bound Numeric vector. Lower bounds for the parameters used in estimation.
#' @param upper_bound Numeric vector. Upper bounds for the parameters used in estimation.
#' @param nc Integer. Number of cores to use for parallel processing.
#' @param cl Optional cluster object created by \code{makeCluster} for parallel computation. Defaults to \code{NULL}.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A matrix containing simulation results. Each row corresponds to one simulation replication and includes parameter estimates along with an accuracy measure (e.g., hit rate or goodness-of-fit).
#'
#' @examples
#' \dontrun{
#' # Example simulation with the "Clayton Gumbel BEKK" copula
#' type <- "Clayton Gumbel BEKK"
#' n <- 1000
#' amount <- 10
#' true_par <- c(0.05, 0.01, 0.95, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.8, 0.2, 2, 3)
#' lower_bound <- rep(0, length(true_par))
#' upper_bound <- c(rep(1, 11), 1, 1, 10, 10)
#' nc <- 4
#' seed <- 12345
#'
#' # Run the simulation
#' sim_results <- Sim_MSCMGARCH(type, n, amount, true_par, lower_bound, upper_bound, nc, seed=seed)
#'
#' # View the results
#' head(sim_results)
#' }
#'
#' @export
Sim_MSCMGARCH<-function(type,n,amount,true_par,lower_bound,upper_bound,nc,cl=NULL,seed){
  
  
  amount_parallel=rep(1,amount)
 
  
  
 
  f=function(x){return(Test(type,n,x,true_par,lower_bound,upper_bound,seed,1))}
  #for replicability
  theta_list<-future_lapply(X=amount_parallel,FUN=f,future.seed=seed)
  
  list_final=unlist(theta_list[[1]])
   for(i in 2:length(theta_list)){
     list_final=rbind(list_final,unlist(theta_list[[i]]))
   }
  return(list_final)
}
