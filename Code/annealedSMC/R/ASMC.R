#' Annealed sequential Monte Carlo
#'
#' Annealed SMC for estimating DE parameters & likelihood parameters
#'
#' @param K Number of particles
#' @param smc_tuning_param List of parameters for the SMC algorithm. Must contain eps and phi.
#' @param data Observed data
#' @param is_unknowPar Boolean vector indicating which parameters are being estimated
#' @param ODEmodel Function defining the ODE model. Must be compatible with deSolve functions.
#' @param ODEparameters Named vector of DE parameters
#' @param likelihood List containing density/mass of observations
#' @param likeliParam_fixed Additional parameters to be passed to the likelihood which are known
#' @param likeliParam_unknown Additional parameters of the likelihood which are to be estimated. NULL if none applicable.
#' @param reference_prior_list List containing reference & prior distributions
#' @param hyperpar List of hyperparameters. Must contain sigma_theta1, sigma_theta2, and ksi
#' @param n_core Number of cores to be used for parallel processing. Default value is 1.
#'
#' @return List of results
#'
#' @export
ASMC <- function(K, smc_tuning_param, data, is_unknownPar, ODEmodel, ODEparameters,
                     likelihood, likeliParam_fixed, likeliParam_unknown=NULL, reference_prior_list, hyperpar, n_core=1){
  # if there are no likelihood parameters to be estimated, add dummy components to prior & reference distributions
  if( is.null(likeliParam_unknown) )
  {
    reference_prior_list$prior_par$likeliParam = NA
    reference_prior_list$logPriorsRatio$likeliParam = function(x,y,z){0}
    reference_prior_list$reference_par$likeliParam = NA
    reference_prior_list$logReferenceRatio$likeliParam = function(x,y,z){0}
    reference_prior_list$logPriorReference$likeliParam = function(x,y,z){0}
    reference_prior_list$reference_sim$likeliParam = function(x){return(NULL)}
  }
  
  with(as.list(c(smc_tuning_param, data, ODEmodel,
                 likelihood, likeliParam_fixed, reference_prior_list, hyperpar)), {
                   # initialize storage list
                   alpha = list()         # a list for tempering parameters
                   alphaDiff = list()
                   ESS = list()           # a list for ESS
                   # initialize values
                   r = 1                  # SMC iteration
                   alpha[[r]] = 0         # tempering parameters
                   alphaDiff[[r]] = 0     # difference b/w two successive tempering parameters
                   logZ = 0
                   n_param = length(ODEparameters)      # number of ODE parameters
                   n_unknownPar = sum(is_unknownPar)  # number of unknown ODE parameters
                   index_unknownPar = which(is_unknownPar)
                   
                   # initialization for i-th particle
                   init <- function(i){
                     # initialize DE parameters from reference distribution
                     theta  = reference_sim$theta(reference_par$theta)
                     names(theta) = names(ODEparameters)
                     
                     # initialize likelihood parameters from reference distribution
                     likeliParam = reference_sim$likeliParam(reference_par$likeliParam)
                     names(likeliParam) = names(likeliParam_unknown)
                     
                     return(list(theta = theta,
                                 accept_theta = rep(F, sum(is_unknownPar)),
                                 select_theta = rep(F, sum(is_unknownPar)),
                                 likeliParam=likeliParam,
                                 accept_LP=rep(F, length(likeliParam)),
                                 select_LP = rep(F, length(likeliParam))))
                   }
                   
                   # Initialization of the K particles
                   particle_list = lapply(1:K, init)
                   
                   # initialize weights
                   # normalized particle weights.
                   # Particles were sampled from the same reference distribution. They have the equal weights.
                   W = rep(1/K,K)
                   # log normalized particle weights
                   logW = log(W)
                   # initialization of the unnormalized weights
                   w = rep(1,K)
                   # log unnormalized particle weights
                   logw = rep(0,K)
                   
                   # initialize cores for parallel processing
                   if( n_core > 1){
                     cl = makeCluster(n_core)
                     clusterExport(cl, varlist=ls(), envir=environment())
                   }
                   
                   while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
                   {
                     cat("iteration:",r,"\n")
                     r = r+1  # SMC iteration
                     
                     # evaluate the log-likelihood for updating alpha
                     u <- rep(0, K)   # incremental log-importance weights
                     
                     if(n_core>1){
                       u = (parSapply(cl, 1:K, function(k){
                         LogL = logLik(ODEmodel,
                                       particle_list[[k]]$theta,
                                       data,
                                       likelihood,
                                       likeliParam_fixed,
                                       particle_list[[k]]$likeliParam,
                                       alpha=1)
                         return (LogL + logPriorReference$theta(particle_list[[k]]$theta, prior_par$theta, reference_par$theta) +
                                   logPriorReference$likeliParam(particle_list[[k]]$likeliParam, prior_par$likeliParam, reference_par$likeliParam))
                       }))
                     }
                     else{
                       u = (sapply(1:K, function(k){
                         LogL = logLik(ODEmodel,
                                       particle_list[[k]]$theta,
                                       data,
                                       likelihood,
                                       likeliParam_fixed,
                                       particle_list[[k]]$likeliParam,
                                       alpha=1)
                         return (LogL + logPriorReference$theta(particle_list[[k]]$theta, prior_par$theta, reference_par$theta) +
                                   logPriorReference$likeliParam(particle_list[[k]]$likeliParam, prior_par$likeliParam, reference_par$likeliParam))
                       }))
                     }
                     
                     # update alpha with bisection
                     alphaDiff[[r]] = bisection( 0, 1,  W, u, phi )
                     alpha[[r]] = alpha[[r-1]] + alphaDiff[[r]]
                     cat("annealing parameter:",alpha[[r]],"\n")
                     
                     # if alpha is set greater than 1, fix by setting to 1
                     if( alpha[[r]]>1 ){
                       alpha[[r]] = 1
                       alphaDiff[[r]] = 1-alpha[[r-1]]
                     }
                     
                     MCMCmove <- function(i, particle_list){
                       with( particle_list[[i]],{
                         # determine which theta to move
                         rand_num <- sample(1:n_unknownPar, 1)
                         which_theta <- index_unknownPar[rand_num]
                         select_theta <- rep(F, n_unknownPar)
                         select_theta[rand_num] <- T
                         
                         # perform n_move MCMC moves
                         for( i in 1:n_move )
                         {
                           # reset acceptance rates
                           accept_theta <- rep(F, n_unknownPar)
                           
                           # sample new theta
                           new_theta <- theta
                           sig = ifelse(sample(1:2, size=1, prob=c(ksi, 1-ksi))==1,
                                        sigma_theta1[which_theta],
                                        sigma_theta2[which_theta])
                           new_theta[which_theta] <- rnorm(1, theta[which_theta], sd=sig)
                           names(new_theta) <- names(theta)
                           
                           # perform MCMC move
                           if(MH_theta(new_theta,
                                       ODEmodel,
                                       theta,
                                       data,
                                       likelihood,
                                       likeliParam_fixed,
                                       likeliParam,
                                       alpha[[r]],
                                       prior_par,
                                       reference_par,
                                       logPriorsRatio,
                                       logReferenceRatio,
                                       is_MCMC = F))
                           {
                             theta  <- new_theta
                             accept_theta[rand_num] <- TRUE
                           }
                         }

                         
                         ## sample likelihood parameters
                         if( is.null(likeliParam_unknown) ){
                           likeliParam = NULL
                         }else{
                           which_LP <- sample(1:length(likeliParam_unknown), size=1)
                           select_LP <- rep(F, length(likeliParam))
                           select_LP[which_LP] <- T
                           
                           for(i in 1:n_move)
                           {
                             accept_LP <- rep(F, length(likeliParam))
                             LP_new <- likeliParam
                             LP_new[which_LP] <- rnorm(1, likeliParam[which_LP], sigma_likeliParam[which_LP])
                             names(LP_new) <- names(likeliParam_unknown)
                             
                             # compute MH acceptance
                             if(MH_likeliParam(LP_new,
                                               likeliParam,
                                               ODEmodel,
                                               theta,
                                               data,
                                               likelihood,
                                               likeliParam_fixed,
                                               alpha[[r]],
                                               prior_par,
                                               reference_par,
                                               logPriorsRatio,
                                               logReferenceRatio,
                                               F)){
                               likeliParam = LP_new
                               accept_LP[which_LP] = T
                             }
                           }
                         }
                         
                         return(list(theta=theta,
                                     accept_theta=accept_theta,
                                     select_theta = select_theta,
                                     likeliParam=likeliParam,
                                     accept_LP=accept_LP,
                                     select_LP = select_LP))
                       })
                     }
                     
                     if(n_core>1){
                       particle_list = parLapply(cl, 1:K,  MCMCmove, particle_list)
                     }else
                       particle_list = lapply(1:K,  MCMCmove, particle_list)
                     
                     n_select_theta <- apply(matrix(unlist(map(particle_list, "select_theta")),nrow=K, byrow=T), 2, sum)
                     n_select_LP <- apply(matrix(unlist(map(particle_list, "select_LP")),nrow=K, byrow=T), 2, sum)
                     
                     accept_rate_theta = apply(matrix(unlist(map(particle_list,"accept_theta")),nrow=K, byrow=T), 2, sum)/n_select_theta
                     cat("theta acceptance rate: ", accept_rate_theta, "\n")
                     accept_rate_LP = apply(matrix(unlist(map(particle_list,"accept_LP")),nrow=K, byrow=T), 2, sum)/n_select_LP
                     cat("likeliParam acceptance rate: ", accept_rate_LP, "\n")
                     
                     # compute the ESS
                     log_incremental_w = alphaDiff[[r]]*u
                     logw = log_incremental_w + logw  # log unnormalized weights
                     logmax = max(logw)
                     W = exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
                     logW = log(W)                                # log normalized weights
                     logZ = logZ + logsum(log_incremental_w + logW)
                     
                     ESS[[r]] = ESS(logW)
                     
                     # resample if ESS below threshold
                     if( ESS[[r]]<eps )
                     {
                       cat("Resample: ESS=", ESS[[r]], '\n')
                       ancestors = systematic_resample( W )
                       particle_list =  particle_list[ancestors]
                       W = rep(1/K,K)
                       logW = log(W)
                       w = rep(1,K)
                       logw = rep(0,K)
                     }
                   }
                   if( n_core > 1 ) stopCluster(cl)
                   
                   return( list(particle_list=particle_list,W=W,logZ=logZ,alpha=alpha))
                 })
  
}

logLik <- function(ODEmodel, ODEparameters, data, likelihood, likeliParam_fixed, likeliParam_unknown, alpha) {
  with(c(data), {
    # combine the likelihood parameters
    likeliParam = unlist(c(likeliParam_fixed, likeliParam_unknown))
    
    # solve the DE trajectory
    out = deSolve::ode(y=initial_state, times=t, func=ODEmodel, parms=unlist(ODEparameters))
    
    # compute the tempered likelihood
    if( partial_likelihood ){
      temp = alpha * sum(sapply(1:length(y), function(i){
        if (is_obs[i])
        {
          logLiks = likelihood$dFun(y[[i]], out[,i+1], likeliParam, log=TRUE)
          return(sum(logLiks[logLiks!=-Inf], na.rm = TRUE))
        }
        else return(0)
      }))
    } else{
      temp = likelihood$dFun(data, out, ODEparameters, likeliParam, log=T)
    }
    
    return(temp)
  })
  
}

# bisection function
# - recursive implementation of bisection algorithm
bisection <- function(low, high, W, u, phi ){
  
  mid <- (low+high)/2
  f.low <- rCESS( W, u, low, phi )
  f.mid <- rCESS( W, u, mid, phi )
  f.high <- rCESS( W, u, high, phi )
  
  # browser()
  if( f.low*f.high>0 )
    stop('Invalid endpoint for bisection.')
  
  try({if( low>=high )
    stop('bisection overlap')
    
    if( (abs(f.mid)<1e-10)||((high-low)/2<1e-10) )
      return( mid )
    if( (f.low*f.mid)<0 )
      return( bisection( low, mid, W, u, phi ) )
    if( (f.high*f.mid)<0 )
      return( bisection( mid, high, W, u, phi ) )
  })
  
  stop('bisection flawed')
}

# log-sum-exponential evaluation of log(sum(w))
logsum <- function(logw)
{
  logmax = max(logw)
  log(sum(exp(logw-logmax)))+logmax
}

# relative conditional effective sample size
rCESS <- function(W, u, a, phi) {
  logw <- a*u          # weight update
  exp(2*logsum(log(W)+logw) - logsum(log(W)+2*logw)) - phi
}

# effective sample size
ESS <- function(logW){
  K <- length(logW)
  logWmax <- max(logW)
  logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(K)
  return(exp(logRESS))
}

# systematic resampling algorithm
systematic_resample <- function( W ){
  K <- length(W)
  U <- runif(1,0,1/K) + 0:(K-1)/K
  W.sum <- cumsum(W)
  N <- rep(NA,K)
  j <- 1
  for( i in 1:K )
  {
    found = F
    while( !found )
    {
      if( U[i]>W.sum[j] )
        j <- j+1
      else
        found = T
    }
    N[i] <- j
  }
  
  return( N )
}