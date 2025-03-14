# library(NMOF)
# library(readxl)
# library(stats)
# library(stats4)
# library(copula)
# library(VineCopula)
# library(parallel)
# library(mvtnorm)
# library(optimParallel)
# library(foreach)
# library(doParallel)
# library(maxLik)
# library(mgarchBEKK)
# library(pbapply)
# library(Rcpp)
# library(RcppArmadillo)
# library(BEKKs)
# 


#' @title Compute BEKK Conditional Covariance Matrices
#' @description Computes a sequence of conditional covariance matrices using the BEKK model specification.
#' @param C Lower triangular matrix representing the constant term in the BEKK model.
#' @param A Matrix of ARCH parameters in the BEKK model.
#' @param G Matrix of GARCH parameters in the BEKK model.
#' @param data Matrix of time series data where rows are observations and columns are variables.
#' @return A list containing the conditional covariance matrices for each time point.
#' @details This function implements the BEKK(1,1) recursion formula:
#'          H_t = CC' + A'(r_{t-1}r_{t-1}')A + G'H_{t-1}G
#'          where the Gt term appears to be missing from the current implementation.
#' @examples
#' # Generate some random bivariate data
#' set.seed(123)
#' data <- matrix(rnorm(200), ncol = 2)
#' 
#' # Define BEKK parameters
#' C <- matrix(c(0.2, 0, 0.1, 0.3), ncol = 2)  # Lower triangular
#' A <- matrix(c(0.3, 0.05, 0.05, 0.3), ncol = 2)
#' G <- matrix(c(0.8, 0.02, 0.02, 0.8), ncol = 2)
#' 
#' # Compute conditional covariance matrices
#' H_list <- comph_bekk(C, A, G, data)
comph_bekk <- function(C, A, G, data) {
  # dimensions
  n      = nrow(data)
  
  # redefine for convenience
  CC     = t(C) %*% C
  At     = t(A)
  Gt     = t(G)
  
  # compute uncond. covariance matrix
  Hu     = (t(data) %*% data) / nrow(data)
  
  # compute H trajectory
  H      = vector(mode = "list", n)
  H[[1]] = Hu
  for (i in 2:n) {
    H[[i]] <- (CC + At %*% (data[i - 1,] %*% t(data[i - 1,])) %*% A
               + Gt %*% H[[i - 1]] %*% G)
  }
  return(H)
}









#' Generate Simulated Time Series Data for Various Markov-Switching Copula Models
#'
#' This function generates simulated bivariate time series data based on specified Markov-switching copula models combined optionally with a BEKK multivariate GARCH process. Supported types include "Clayton Gumbel BEKK", "Normal Gumbel BEKK", "Clayton Gumbel Survival BEKK", "Clayton Gumbel", "Normal Gumbel", and "Clayton Gumbel Survival".
#'
#' @param type A character string specifying the type of Markov-switching copula model to simulate. Implemented options:
#'   - "Clayton Gumbel BEKK"
#'   - "Normal Gumbel BEKK"
#'   - "Clayton Gumbel Survival BEKK"
#'   - "Clayton Gumbel"
#'   - "Normal Gumbel"
#'   - "Clayton Gumbel Survival"
#' @param n Integer specifying the desired length of the time series to simulate.
#' @param par Numeric vector containing all necessary parameters for simulation:
#'   - For BEKK models ("Clayton Gumbel BEKK", "Normal Gumbel BEKK", "Clayton Gumbel Survival BEKK"), `par` includes:
#'     - 11 BEKK parameters (elements for matrices C, A, and G)
#'     - 2 transition probabilities (p and q)
#'     - Copula parameters (Clayton, Gumbel, or Survival, depending on type)
#'   - For non-BEKK models ("Clayton Gumbel", "Normal Gumbel", "Clayton Gumbel Survival"), `par` includes:
#'     - 2 transition probabilities (p and q)
#'     - Copula parameters (Clayton, Gumbel, or Survival)
#'
#' @return A list containing two unnamed elements:
#'   - First element: Numeric matrix of simulated time series data (dimension: n x 2)
#'   - Second element: Numeric vector indicating the Markov-switching regime (0 or 1) at each time point
#'
#' @examples
#' \dontrun{
#' type <- "Clayton Gumbel BEKK"
#' n <- 1000
#' true_par <- c(0.05, 0.01, 0.03, 0.1, 0.05, 0.05, 0.1, 0.8, 0.05, 0.05, 0.8,
#'               0.8, 0.2, 2, 3)
#'
#' simulated_data <- example_Creator(type, n, true_par)
#'
#' # Visualizing the first few rows of simulated series
#' head(simulated_data[[1]])
#'
#' # Visualizing state regimes
#' plot(simulated_data[[2]], type='l', main='Regime states', ylab='State',
#'      xlab='Time', col='blue')
#' }
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom VineCopula BiCopSim
#' @export
example_Creator<-function(type,n, par){
   if(type == "Clayton Gumbel BEKK"){
  series=rmvnorm(n+200,mean=rep(0,2))
  bekk_par=par[1:11]
  p=par[12]
  q=par[13]
  param1=par[14]
  param2=par[15]

  C=matrix(c(bekk_par[1],0,bekk_par[2:3]),2)
  Ct=t(C)
  CC=Ct%*%C
  A=t(matrix(bekk_par[4:7],2))
  At=t(A)
  G=t(matrix(bekk_par[8:11],2))
  Gt=t(G)
  Uncond_var=matrix(solve(diag(4) - t(kronecker(A, A)) - t(kronecker(G, G))) %*% c(CC),2)



  covarianz_matrix=cor(BiCopSim(100000,4,param2))
  copula_chol_inv_gumb=solve(t(chol(covarianz_matrix)))
  covarianz_matrix_clay=cor(BiCopSim(100000,3,param1))
  copula_chol_inv_clay=solve(t(chol(covarianz_matrix_clay)))
  state=numeric(n+200)
  H=Uncond_var
  H_decomp=eigen_value_decomposition(H)
  state[1]=rbinom(1,1,0.5)
  if(state[1]==0)
  {
    series[1,]=H_decomp%*%copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
  } else
  {
    series[1,]=H_decomp%*%copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param2)))
  }


  for(i in 2:(n+200))
  {
    state[i]=state[i-1]*rbinom(1,1,q)+(1-state[i-1])*rbinom(1,1,1-p)
    H=CC+At%*%series[i-1,]%*%t(series[i-1,])%*%A+Gt%*%H%*%G
    H_decomp=eigen_value_decomposition(H)
    if(state[i]==0)
    {
      series[i,]=H_decomp%*%copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
    } else
    {
      series[i,]=H_decomp%*%copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param2)))
    }

  }

  return(list(series[201:(n+200),],state[201:(n+200)]))
   }
  if(type == "Normal Gumbel BEKK"){
    series=rmvnorm(n+200,mean=rep(0,2))
    bekk_par=par[1:11]
    p=par[12]
    q=par[13]
    param1=par[14]
   
    
    C=matrix(c(bekk_par[1],0,bekk_par[2:3]),2)
    Ct=t(C)
    CC=Ct%*%C
    A=t(matrix(bekk_par[4:7],2))
    At=t(A)
    G=t(matrix(bekk_par[8:11],2))
    Gt=t(G)
    Uncond_var=matrix(solve(diag(4) - t(kronecker(A, A)) - t(kronecker(G, G))) %*% c(CC),2)
    
    
    
    covarianz_matrix=cor(BiCopSim(100000,4,param1))
    copula_chol_inv_gumb=solve(t(chol(covarianz_matrix)))
    
    state=numeric(n+200)
    H=Uncond_var
    H_decomp=eigen_value_decomposition(H)
    state[1]=rbinom(1,1,0.5)
    if(state[1]==0)
    {
      series[1,]=H_decomp%*%t(rmvnorm(1,sigma=diag(2)))
    } else
    {
      series[1,]=H_decomp%*%copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param1)))
    }
    
    
    for(i in 2:(n+200))
    {
      state[i]=state[i-1]*rbinom(1,1,q)+(1-state[i-1])*rbinom(1,1,1-p)
      H=CC+At%*%series[i-1,]%*%t(series[i-1,])%*%A+Gt%*%H%*%G
      H_decomp=eigen_value_decomposition(H)
      if(state[i]==0)
      {
        series[i,]=H_decomp%*%t(rmvnorm(1,sigma=diag(2)))
      } else 
      {
        series[i,]=H_decomp%*%copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param1)))
      }
      
    }
    
    return(list(series[201:(n+200),],state[201:(n+200)]))
  }
  if(type == "Clayton Gumbel Survival BEKK"){
    series=rmvnorm(n+200,mean=rep(0,2))
    bekk_par=par[1:11]
    p=par[12]
    q=par[13]
    param1=par[14]
    param2=par[15]
    
    C=matrix(c(bekk_par[1],0,bekk_par[2:3]),2)
    Ct=t(C)
    CC=Ct%*%C
    A=t(matrix(bekk_par[4:7],2))
    At=t(A)
    G=t(matrix(bekk_par[8:11],2))
    Gt=t(G)
    Uncond_var=matrix(solve(diag(4) - t(kronecker(A, A)) - t(kronecker(G, G))) %*% c(CC),2)
    
    
    
    covarianz_matrix=cor(qnorm(BiCopSim(100000,24,param2)))
    copula_chol_inv_gumb=solve(t(chol(covarianz_matrix)))
    covarianz_matrix_clay=cor(qnorm(BiCopSim(100000,3,param1)))
    copula_chol_inv_clay=solve(t(chol(covarianz_matrix_clay)))
    state=numeric(n+200)
    H=Uncond_var
    H_decomp=eigen_value_decomposition(H)
    state[1]=rbinom(1,1,0.5)
    if(state[1]==0)
    {
      series[1,]=H_decomp%*%copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
    } else
    {
      series[1,]=H_decomp%*%copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,24,param2)))
    }
    
    
    for(i in 2:(n+200))
    {
      state[i]=state[i-1]*rbinom(1,1,q)+(1-state[i-1])*rbinom(1,1,1-p)
      H=CC+At%*%series[i-1,]%*%t(series[i-1,])%*%A+Gt%*%H%*%G
      H_decomp=eigen_value_decomposition(H)
      if(state[i]==0)
      {
        series[i,]=H_decomp%*%copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
      } else
      {
        series[i,]=H_decomp%*%copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,24,param2)))
      }
      
    }
    
    return(list(series[201:(n+200),],state[201:(n+200)]))
  }
  if(type == "Clayton Gumbel"){
    series=rmvnorm(n,mean=rep(0,2))

    p=par[1]
    q=par[2]
    param1=par[3]
    param2=par[4]

    covarianz_matrix=cor(qnorm(BiCopSim(100000,4,param2)))
    copula_chol_inv_gumb=solve(t(chol(covarianz_matrix)))
    covarianz_matrix_clay=cor(qnorm(BiCopSim(100000,3,param1)))
    copula_chol_inv_clay=solve(t(chol(covarianz_matrix_clay)))
    state=numeric(n)

    state[1]=rbinom(1,1,0.5)
    if(state[1]==0)
    {
      series[1,]=copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
    } else
    {
      series[1,]= copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param2)))
    }


    for(i in 2:n)
    {
      state[i]=state[i-1]*rbinom(1,1,q)+(1-state[i-1])*rbinom(1,1,1-p)
      if(state[i]==0)
      {
        series[i,]=copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
      } else
      {
        series[i,]= copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param2)))
      }

    }

    return(list(series,state))
  }

  if(type == "Normal Gumbel"){
    p=par[1]
    q=par[2]
    param1=par[3]
    series=rmvnorm(n,mean=rep(0,2))

    copula_gumb=cov(qnorm(BiCopSim(100000,4,param1)))
    copula_chol_inv_gumb=solve(t(chol(copula_gumb)))
    state=numeric(n)
    state[1]=0
    for(i in 2:n)
    {
      state[i]=state[i-1]*rbinom(1,1,q)+(1-state[i-1])*rbinom(1,1,1-p)
      if(state[i]==0)
      {
        series[i,]=rmvnorm(1,mean=rep(0,2))
      } else
      {
        series[i,]= copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,4,param1)))
      }


      # compute uncond. covariance matrix



    }

    return(list(series,state))

  }
  if(type== "Clayton Gumbel Survival"){

    series=rmvnorm(n,mean=rep(0,2))
    p=par[1]
    q=par[2]
    param1=par[3]
    param2=par[4]



    covarianz_matrix_clay=cov(qnorm(BiCopSim(100000,3,param1)))
    copula_chol_inv_clay=solve(t(chol(covarianz_matrix_clay)))
    covarianz_matrix=cov(qnorm(BiCopSim(100000,24,param2)))
    copula_chol_inv_gumb=solve(t(chol(covarianz_matrix)))
    state=numeric(n)
    state[1]=0
    series[1,]=copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
    for(i in 2:n)
    {
      state[i]=state[i-1]*rbinom(1,1,q)+(1-state[i-1])*rbinom(1,1,1-p)
      if(state[i]==0)
      {
        series[i,]=copula_chol_inv_clay%*%t(qnorm(BiCopSim(1,3,param1)))
      } else
      {
        series[i,]= copula_chol_inv_gumb%*%t(qnorm(BiCopSim(1,24,param2)))
      }


      # compute uncond. covariance matrix



    }

    return(list(series,state))
  }
}


Filtered_Probabilities<-function(type,par,n,series,j){
  if(type== "Clayton Gumbel"){
  #generate vector for results
  Prob=numeric(n)
  #Filter Probabilities
    p=par[1]
    q=par[2]
    param1=par[3]
    par1=par[4]
    set.seed(j)
      covarianz_matrix=cov(qnorm(BiCopSim(100000,4,par1)))
      copula_chol=t(chol(covarianz_matrix))
      covarianz_gumb=cov(qnorm(BiCopSim(100000,3,param1)))
      gumbel=t(chol(covarianz_gumb))

      #Probability for beeing in state 1 in time t=1
      p1=1
      #filtered_prob=p1
      #LogLikelihood for first time t=0
      #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
      #H_inv=pblapply(H,function(f){solve(t(chol(f)))})
      #series[1,]=H_inv[[1]]%*%series[1,]
      llv=numeric(length(series)/2)
      llv[1]=log(dnorm((gumbel%*%series[1,])[1])*dnorm((gumbel%*%series[1,])[2])*BiCopPDF(pnorm((gumbel%*%series[1,]))[1], pnorm((gumbel%*%series[1,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))
      )
      Prob[1]=1

      for(i in 2:(length(series)/2))
      {
        #series[i,]=H_inv[[i]]%*%series[i,]
        likelihood=(p1*p+(1-p1)*(1-q))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 4, par1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))
        llv[i]=log(likelihood)
        p1t=(p*p1+(1-q)*(1-p1))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1,0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))/likelihood
        #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
        #Filtered probabilities for state 1
        Prob[i]=p1t
        #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
        #filtered_prob=cbind(filtered_prob,p_filtered)
        #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
        #print()
        p1<-p1t

      }
      return(Prob)

  }
  if(type== "Clayton Gumbel BEKK"){
    #generate vector for results
    Prob=numeric(n)
    #Filter Probabilities
    bekk_par=par[1:11]
    p=par[12]
    q=par[13]
    param1=par[14]
    par1=par[15]
    covarianz_matrix=cov(qnorm(BiCopSim(100000,4,par1)))
    copula_chol=t(chol(covarianz_matrix))
    covarianz_gumb=cov(qnorm(BiCopSim(100000,3,param1)))
    gumbel=t(chol(covarianz_gumb))

    #Probability for beeing in state 1 in time t=1
    p1=1
    C=t(matrix(c(bekk_par[1],0,bekk_par[2:3]),2))

    A=matrix(bekk_par[4:7],2)

    G=matrix(bekk_par[8:11],2)

    H=comph_bekk(C,A,G,series)


    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    det_H=pblapply(H_inv,det)
    series[1,]=H_inv[[1]]%*%series[1,]
    llv=numeric(length(series)/2)
    llv[1]=log(dnorm((gumbel%*%series[1,])[1])*dnorm((gumbel%*%series[1,])[2])*BiCopPDF(pnorm((gumbel%*%series[1,]))[1], pnorm((gumbel%*%series[1,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel)*det_H[[1]])
    )
    Prob[1]=1

    for(i in 2:(length(series)/2))
    {
      series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))*det_H[[i]]+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 4, par1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))*det_H[[i]]
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1,0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))*det_H[[i]]/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t

    }
    return(Prob)

  }
  if(type=="Normal Gumbel"){
    Prob=numeric(n)
    #Filter Probabilities
    p=par[1]
    q=par[2]
    param1=par[3]
    set.seed(j)
    covarianz_matrix=cov(qnorm(BiCopSim(100000,4,param1)))
    copula_chol=t(chol(covarianz_matrix))


    #Probability for beeing in state 1 in time t=1
    p1=1
    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    #H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    #series[1,]=H_inv[[1]]%*%series[1,]
    n=length(series)/2
    llv=numeric(n)
    llv[1]=log(dmvnorm(series[1,]))
    Prob[1]=1

    for(i in 2:n)
    {
      #series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dmvnorm(series[i,])+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 4, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dmvnorm(series[i,])/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t

    }
    return(Prob)
  }
  if(type== "Clayton Gumbel Survival BEKK"){
    #generate vector for results
    Prob=numeric(n)
    #Filter Probabilities
    bekk_par=par[1:11]
    p=par[12]
    q=par[13]
    param1=par[14]
    par1=par[15]
    covarianz_matrix=cov(qnorm(BiCopSim(100000,24,par1)))
    copula_chol=t(chol(covarianz_matrix))
    covarianz_gumb=cov(qnorm(BiCopSim(100000,3,param1)))
    gumbel=t(chol(covarianz_gumb))
    
    #Probability for beeing in state 1 in time t=1
    p1=1
    C=t(matrix(c(bekk_par[1],0,bekk_par[2:3]),2))
    
    A=matrix(bekk_par[4:7],2)
    
    G=matrix(bekk_par[8:11],2)
    
    H=comph_bekk(C,A,G,series)
    
    
    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    det_H=pblapply(H_inv,det)
    series[1,]=H_inv[[1]]%*%series[1,]
    llv=numeric(length(series)/2)
    llv[1]=log(dnorm((gumbel%*%series[1,])[1])*dnorm((gumbel%*%series[1,])[2])*BiCopPDF(pnorm((gumbel%*%series[1,]))[1], pnorm((gumbel%*%series[1,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel)*det_H[[1]])
    )
    Prob[1]=1
    
    for(i in 2:(length(series)/2))
    {
      series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))*det_H[[i]]+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 24, par1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))*det_H[[i]]
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1,0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))*det_H[[i]]/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t
      
    }
    return(Prob)
    
  }
  if(type== "Normal Gumbel BEKK"){
    #generate vector for results
    Prob=numeric(n)
    #Filter Probabilities
    bekk_par=par[1:11]
    p=par[12]
    q=par[13]
    param1=par[14]
    
    covarianz_gumb=cov(qnorm(BiCopSim(100000,4,param1)))
    gumbel=t(chol(covarianz_gumb))
    
    #Probability for beeing in state 1 in time t=1
    p1=1
    C=t(matrix(c(bekk_par[1],0,bekk_par[2:3]),2))
    
    A=matrix(bekk_par[4:7],2)
    
    G=matrix(bekk_par[8:11],2)
    
    H=comph_bekk(C,A,G,series)
    
    
    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    det_H=pblapply(H_inv,det)
    series[1,]=H_inv[[1]]%*%series[1,]
    llv=numeric(length(series)/2)
    llv[1]=log(prod(dnorm((series[1,])))*det_H[[1]])
    
    Prob[1]=1
    
    for(i in 2:(length(series)/2))
    {
      series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dnorm((series[i,])[1])*dnorm((series[i,])[2])*det_H[[i]]+(p1*(1-p)+(1-p1)*(q))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 4, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))*det_H[[i]]
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dnorm((series[i,])[1])*dnorm((series[i,])[2])*det_H[[i]]/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t
      
    }
    return(Prob)
    
  }
  if(type=="Normal Gumbel"){
    Prob=numeric(n)
    #Filter Probabilities
    p=par[1]
    q=par[2]
    param1=par[3]
    set.seed(j)
    covarianz_matrix=cov(qnorm(BiCopSim(100000,4,param1)))
    copula_chol=t(chol(covarianz_matrix))
    
    
    #Probability for beeing in state 1 in time t=1
    p1=1
    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    #H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    #series[1,]=H_inv[[1]]%*%series[1,]
    n=length(series)/2
    llv=numeric(n)
    llv[1]=log(dmvnorm(series[1,]))
    Prob[1]=1
    
    for(i in 2:n)
    {
      #series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dmvnorm(series[i,])+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 4, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dmvnorm(series[i,])/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t
      
    }
    return(Prob)
  }
  if(type=="Normal Gumbel"){
    Prob=numeric(n)
    #Filter Probabilities
    p=par[1]
    q=par[2]
    param1=par[3]
    set.seed(j)
    covarianz_matrix=cov(qnorm(BiCopSim(100000,4,param1)))
    copula_chol=t(chol(covarianz_matrix))
    
    
    #Probability for beeing in state 1 in time t=1
    p1=1
    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    #H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    #series[1,]=H_inv[[1]]%*%series[1,]
    n=length(series)/2
    llv=numeric(n)
    llv[1]=log(dmvnorm(series[1,]))
    Prob[1]=1
    
    for(i in 2:n)
    {
      #series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dmvnorm(series[i,])+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 4, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dmvnorm(series[i,])/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t
      
    }
    return(Prob)
  }
  if(type == "Clayton Gumbel Survival"){
    Prob=numeric(n)
    #Filter Probabilities
    p=par[1]
    q=par[2]
    param1=par[3]
    par1=par[4]
    set.seed(j)
    covarianz_matrix=cov(qnorm(BiCopSim(100000,24,par1)))
    copula_chol=t(chol(covarianz_matrix))
    covarianz_gumb=cov(qnorm(BiCopSim(100000,3,param1)))
    gumbel=t(chol(covarianz_gumb))

    #Probability for beeing in state 1 in time t=1
    p1=1
    #filtered_prob=p1
    #LogLikelihood for first time t=0
    #H=comph_bekk(matrix(c(BEKK[1],0,BEKK[2],BEKK[3]),2),matrix(c(BEKK[4], BEKK[5], BEKK[6], BEKK[7]),2),matrix(BEKK[8:11],2),series)
    #H_inv=pblapply(H,function(f){solve(t(chol(f)))})
    #series[1,]=H_inv[[1]]%*%series[1,]
    llv=numeric(length(series)/2)
    llv[1]=log(dnorm((gumbel%*%series[1,])[1])*dnorm((gumbel%*%series[1,])[2])*BiCopPDF(pnorm((gumbel%*%series[1,]))[1], pnorm((gumbel%*%series[1,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))
    )
    Prob[1]=1

    for(i in 2:(length(series)/2))
    {
      #series[i,]=H_inv[[i]]%*%series[i,]
      likelihood=(p1*p+(1-p1)*(1-q))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1, 0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))+(p1*(1-p)+(1-p1)*(q))*dnorm((copula_chol%*%series[i,])[1])*dnorm((copula_chol%*%series[i,])[2])*BiCopPDF(pnorm((copula_chol%*%series[i,]))[1], pnorm((copula_chol%*%series[i,]))[2], 24, par1, 0, obj = NULL,check.pars = TRUE)*abs(det(copula_chol))
      llv[i]=log(likelihood)
      p1t=(p*p1+(1-q)*(1-p1))*dnorm((gumbel%*%series[i,])[1])*dnorm((gumbel%*%series[i,])[2])*BiCopPDF(pnorm((gumbel%*%series[i,]))[1], pnorm((gumbel%*%series[i,]))[2], 3, param1,0, obj = NULL,check.pars = TRUE)*abs(det(gumbel))/likelihood
      #p1t=(1-p11)*(1-p1)*dmvnorm(t(chol_covMatrix%*%series[i,]),sigma=covMatrix)*det(chol_covMatrix)/likelihood
      #Filtered probabilities for state 1
      Prob[i]=p1t
      #p_filtered=(p1t*p1*dmvnorm(series[i,]))/likelihood
      #filtered_prob=cbind(filtered_prob,p_filtered)
      #print((p1t*p1*dmvnorm(series[i,]))/likelihood)
      #print()
      p1<-p1t

    }
    return(Prob)
  }
}

#' BEKK Parameter Estimation Function
#'
#' Estimates BEKK parameters using a combination of random grid search for initial values and the Berndt-Hall-Hall-Hausman (BHHH) optimization algorithm.
#'
#' @param r Numeric matrix. The input data (time series) for estimating the BEKK model.
#' @param init_values Either "random" or a numeric vector of initial parameter values.
#' @param nc Integer. Number of cores for parallel processing (currently unused, placeholder for future use).
#' @param max_iter Integer. Maximum number of iterations for the optimization algorithm.
#' @param crit Numeric. Convergence criterion for stopping optimization.
#'
#' @return A list containing:
#'   - `C0`: Numeric matrix representing the estimated C matrix.
#'   - `A`: Numeric matrix of estimated A parameters.
#'   - `G`: Numeric matrix of estimated G parameters.
#'   - `full_result`: Full output from the optimization function (`bhh_bekk`).
#'
#' @examples
#' \dontrun{
#' simulated_series <- example_Creator("Clayton Gumbel BEKK", 1000, true_par)
#' estimated_bekk <- bekk(r=simulated_series[[1]], init_values="random", nc=4)
#' estimated_C0 <- estimated_bekk$C0
#' estimated_A <- estimated_bekk$A
#' estimated_G <- estimated_bekk$G
#' }
#'
#' @export
bekk <- function(r, init_values = "random", nc = 1, max_iter = 100, crit = 1e-6){

  if (identical(init_values, "random")) {
    initial_search <- random_grid_search_BEKK(r)
    theta_init <- initial_search$theta
  } else if(is.numeric(init_values)){
    theta_init <- init_values
  } else {
    stop("init_values must be 'random' or numeric vector of initial parameters")
  }

  estimated_params <- bhh_bekk(r = r, theta = theta_init, max_iter = max_iter, crit = crit)
  estimated_theta <- estimated_params$theta

  d <- ncol(r)
  num_C_params <- d*(d+1)/2
  num_A_params <- d^2
  num_G_params <- d^2

  C0 <- matrix(estimated_theta[1:num_C_params], d, d)
  A <- matrix(estimated_theta[(num_C_params+1):(num_C_params+num_A_params)], d, d)
  G <- matrix(estimated_theta[(num_C_params+num_A_params+1):(num_C_params+num_A_params+num_G_params)], d, d)

  return(list(C0=C0, A=A, G=G, full_result=estimated_params))
}



Test<-function(type,n,amount,true_par,lower_bound,upper_bound,seed,nc){
  
  
  if(type=="Clayton Gumbel"){
    par_est=matrix(0,nrow = amount,ncol=5)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      
      BEKK1=true_par[1:11]
      ex=example_Creator("Clayton Gumbel BEKK",n,true_par)
      #dimension noch einmal testen
      
      
      
      start_val=random_grid_search_clayton_gumbel(BEKK1,true_par[12:15],ex[[1]],1)$theta_optim
      
      result<-optim(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,method="BFGS",control=list(fnscale=-1),bekk=BEKK1,r=ex[[1]])
      
      correct_states=0
      Filtered_P=Filtered_Probabilities("Clayton Gumbel BEKK",c(BEKK1,result$par),n,ex[[1]],i)
      for(j in 1:n){
        if(ex[[2]][j]==0 && Filtered_P[j]>= 0.5 || ex[[2]][j]==1 && Filtered_P[j]<= 0.5 ){correct_states=correct_states+1}
      }
      
      hit_rate[i]=correct_states/n
      par_est[i,1:4]=result$par
      par_est[i,5]= hit_rate[i]
      
    }
    return(par_est)
    
  }
  if(type=="Clayton Gumbel BEKK"){
    par_est=matrix(0,nrow = amount,ncol=16)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      
      ex=example_Creator(type,n,true_par)
      #dimension noch einmal testen
      
      BEKK1=bekk(ex[[1]])
      
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
      par_est[i,1:11]=BEKK1
      
      #nc1=rep(nc,nc)
      # paralell_grid_search<-function(bekk, theta, r, nc){
      #   return(random_grid_search_clayton_gumbel(bekk, theta, r, nc))
      #   }
      # theta_list <- parLapply(X=nc1,fun=paralell_grid_search,bekk=BEKK1,theta=true_par[12:15], r = ex[[1]])
      # max_index <- which.max(sapply(X=theta_list, '[[', 'best_val'))
      # theta <- theta_list[[max_index]]
      # start_val <- theta[[1]]
      start_val=random_grid_search_clayton_gumbel(BEKK1,true_par[12:15],ex[[1]],1)$theta_optim
      #ML<-function(par){return(loglike_Clayton_Gumbel(BEKK1,par, ex[[1]]))}
      #result<-optimParallel(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,lower=c(0,0,0,1),upper=c(1,1,17,17),control=list(fnscale=-1),cl=cl, bekk=BEKK1 , r=ex[[1]] )
      result<-optim(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,method="BFGS",control=list(fnscale=-1),bekk=BEKK1,r=ex[[1]])
      correct_states=0
      Filtered_P=Filtered_Probabilities(type,c(BEKK1,result$par),n,ex[[1]],i)
      for(j in 1:n){
        if(ex[[2]][j]==0 && Filtered_P[j]>= 0.5 || ex[[2]][j]==1 && Filtered_P[j]<= 0.5 ){correct_states=correct_states+1}
      }
      hit_rate[i]=correct_states/n
      par_est[i,12:15]=result$par
      par_est[i,16]= hit_rate[i]
      
    }
    return(par_est)

  }
  if(type=="BEKK"){
    par_est=matrix(0,nrow = amount,ncol=11)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      
      ex=example_Creator("Clayton Gumbel BEKK",n,true_par)
      #dimension noch einmal testen
      
      BEKK1=bekk(ex[[1]],init_values="random",nc=1)
      print(summary(BEKK1))
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
      par_est[i,]=BEKK1
      
      #nc1=rep(nc,nc)
      # paralell_grid_search<-function(bekk, theta, r, nc){
      #   return(random_grid_search_clayton_gumbel(bekk, theta, r, nc))
      #   }
      # theta_list <- parLapply(X=nc1,fun=paralell_grid_search,bekk=BEKK1,theta=true_par[12:15], r = ex[[1]])
      # max_index <- which.max(sapply(X=theta_list, '[[', 'best_val'))
      # theta <- theta_list[[max_index]]
      # start_val <- theta[[1]]
      
      
    }
    return(par_est)
    
  }
  if(type=="Clayton Gumbel Survival BEKK"){
    par_est=matrix(0,nrow = amount,ncol=16)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      
      ex=example_Creator(type,n,true_par)
     
      BEKK1=bekk(ex[[1]])
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
      par_est[i,1:11]=BEKK1
      
      
      start_val=random_grid_search_clayton_gumbel90(BEKK1,true_par[12:15],ex[[1]],1)$theta_optim
      #ML<-function(par){return(loglike_Clayton_Gumbel(BEKK1,par, ex[[1]]))}
      #result<-optimParallel(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,lower=c(0,0,0,1),upper=c(1,1,17,17),control=list(fnscale=-1),cl=cl, bekk=BEKK1 , r=ex[[1]] )
      result<-optim(par=as.vector(start_val),fn=loglike_Clayton_Gumbel90,method="BFGS",control=list(fnscale=-1),bekk=BEKK1,r=ex[[1]])
      correct_states=0
      Filtered_P=as.vector(FilterProbs_Clayton_Gumbel90(BEKK1,result$par,ex[[1]]))
      for(j in 1:n){
        if(ex[[2]][j]==0 && Filtered_P[j]>= 0.5 || ex[[2]][j]==1 && Filtered_P[j]<= 0.5 ){correct_states=correct_states+1}
      }
      
      par_est[i,12:15]=result$par
      par_est[i,16]= correct_states/n
    }
    return(par_est)
    
  }
  if(type=="Clayton Gumbel Survival"){
    par_est=matrix(0,nrow = amount,ncol=5)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      BEKK1=true_par[1:11]
      
      ex=example_Creator("Clayton Gumbel Survival BEKK",n,true_par)
      
      
      
      start_val=random_grid_search_clayton_gumbel90(BEKK1,true_par[12:15],ex[[1]],1)$theta_optim
      #ML<-function(par){return(loglike_Clayton_Gumbel(BEKK1,par, ex[[1]]))}
      #result<-optimParallel(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,lower=c(0,0,0,1),upper=c(1,1,17,17),control=list(fnscale=-1),cl=cl, bekk=BEKK1 , r=ex[[1]] )
      result<-optim(par=as.vector(start_val),fn=loglike_Clayton_Gumbel90,method="BFGS",control=list(fnscale=-1),bekk=BEKK1,r=ex[[1]])
      correct_states=0
      Filtered_P=as.vector(FilterProbs_Clayton_Gumbel90(BEKK1,result$par,ex[[1]]))
      for(j in 1:n){
        if(ex[[2]][j]==0 && Filtered_P[j]>= 0.5 || ex[[2]][j]==1 && Filtered_P[j]<= 0.5 ){correct_states=correct_states+1}
      }
      
      par_est[i,1:4]=result$par
      par_est[i,5]= correct_states/n
    }
    return(par_est)
    
  }
  if(type=="Normal Gumbel BEKK"){
    par_est=matrix(0,nrow = amount,ncol=15)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      
      ex=example_Creator(type,n,true_par)
      
      
      BEKK1=bekk(ex[[1]])
      
      BEKK1=c(c(BEKK1$C0)[c(1,3,4)],c(BEKK1$A),c(BEKK1$G))
      par_est[i,1:11]=BEKK1
      
      #nc1=rep(nc,nc)
      # paralell_grid_search<-function(bekk, theta, r, nc){
      #   return(random_grid_search_clayton_gumbel(bekk, theta, r, nc))
      #   }
      # theta_list <- parLapply(X=nc1,fun=paralell_grid_search,bekk=BEKK1,theta=true_par[12:15], r = ex[[1]])
      # max_index <- which.max(sapply(X=theta_list, '[[', 'best_val'))
      # theta <- theta_list[[max_index]]
      # start_val <- theta[[1]]
      start_val=random_grid_search_normal_gumbel1(BEKK1,true_par[12:14],ex[[1]],1)$theta_optim
      #ML<-function(par){return(loglike_Clayton_Gumbel(BEKK1,par, ex[[1]]))}
      #result<-optimParallel(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,lower=c(0,0,0,1),upper=c(1,1,17,17),control=list(fnscale=-1),cl=cl, bekk=BEKK1 , r=ex[[1]] )
      result<-optim(par=start_val,fn=loglike_Normal_Gumbel1,method="BFGS",control=list(fnscale=-1),bekk=BEKK1,r=ex[[1]])
      correct_states=0
      Filtered_P=Filtered_Probabilities(type,c(BEKK1,result$par),n,ex[[1]],i)
      for(j in 1:n){
        if(ex[[2]][j]==0 && Filtered_P[j]>= 0.5 || ex[[2]][j]==1 && Filtered_P[j]<= 0.5 ){correct_states=correct_states+1}
      }
      hit_rate[i]=correct_states/n
      par_est[i,12:14]=result$par
      par_est[i,15]= hit_rate[i]
      
    }
    return(par_est)
    
  }
  if(type=="Normal Gumbel"){
    par_est=matrix(0,nrow = amount,ncol=4)
    hit_rate=numeric(amount)
    for(i in 1: amount)  {
      BEKK1=true_par[1:11]
      
      ex=example_Creator("Normal Gumbel BEKK",n,true_par)
      #dimension noch einmal testen
      
    
      
      #nc1=rep(nc,nc)
      # paralell_grid_search<-function(bekk, theta, r, nc){
      #   return(random_grid_search_clayton_gumbel(bekk, theta, r, nc))
      #   }
      # theta_list <- parLapply(X=nc1,fun=paralell_grid_search,bekk=BEKK1,theta=true_par[12:15], r = ex[[1]])
      # max_index <- which.max(sapply(X=theta_list, '[[', 'best_val'))
      # theta <- theta_list[[max_index]]
      # start_val <- theta[[1]]
      start_val=random_grid_search_normal_gumbel1(BEKK1,true_par[12:14],ex[[1]],1)$theta_optim
      #ML<-function(par){return(loglike_Clayton_Gumbel(BEKK1,par, ex[[1]]))}
      #result<-optimParallel(par=as.vector(start_val),fn=loglike_Clayton_Gumbel,lower=c(0,0,0,1),upper=c(1,1,17,17),control=list(fnscale=-1),cl=cl, bekk=BEKK1 , r=ex[[1]] )
      result<-optim(par=start_val,fn=loglike_Normal_Gumbel1,method="BFGS",control=list(fnscale=-1),bekk=BEKK1,r=ex[[1]])
      correct_states=0
      Filtered_P=Filtered_Probabilities("Normal Gumbel BEKK",c(BEKK1,result$par),n,ex[[1]],i)
      for(j in 1:n){
        if(ex[[2]][j]==0 && Filtered_P[j]>= 0.5 || ex[[2]][j]==1 && Filtered_P[j]<= 0.5 ){correct_states=correct_states+1}
      }
      hit_rate[i]=correct_states/n
      par_est[i,1:3]=result$par
      par_est[i,4]= hit_rate[i]
      
    }
    return(par_est)
    
  }
  
}

