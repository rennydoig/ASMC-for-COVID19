MH_theta <- function(theta_new, ODEmodel, theta_old, data, likelihood,
                          likeliParam_fixed, likeliParam_unknown, alpha, prior_par, reference_par,
                          logPriorsRatio, logReferenceRatio,  is_MCMC = F){
  # # check whether being used in MCMC or aSMC
  if( missing(alpha) )
  {
    alpha = 1
    reference_par = logReferenceRatio = NULL
    is_MCMC = T
  }

  # compute the numerator and denominator of the MH ratio

  num <- logLik(ODEmodel, theta_new, data, likelihood, likeliParam_fixed, likeliParam_unknown, alpha)
  den <- logLik(ODEmodel, theta_old, data, likelihood, likeliParam_fixed, likeliParam_unknown, alpha)

  # compute the MH ratio
  ratio = ifelse(is_MCMC,
                 min(1, exp(num-den+
                              logPriorsRatio$theta(theta_new, theta_old, prior_par$theta))),
                 min(1, exp(num-den+
                              alpha*logPriorsRatio$theta(theta_new, theta_old, prior_par$theta)+
                              (1-alpha)*logReferenceRatio$theta(theta_new, theta_old, reference_par$theta))))

  # accept-reject step
  U <- runif(1)
  return( U<ratio )
}


MH_likeliParam <- function(LP_new, LP_old, ODEmodel, theta, data, likelihood,
                         likeliParam_fixed, alpha, prior_par, reference_par,
                         logPriorsRatio, logReferenceRatio,  is_MCMC = F){
  # # check whether being used in MCMC or aSMC
  if( missing(alpha) )
  {
    alpha = 1
    reference_par = logReferenceRatio = NULL
    is_MCMC = T
  }

  # compute the numerator and denominator of the MH ratio
  num <- logLik2_LP(ODEmodel, theta, data, likelihood, likeliParam_fixed, LP_new, alpha)
  den <- logLik2_LP(ODEmodel, theta, data, likelihood, likeliParam_fixed, LP_old, alpha)

  # compute the MH ratio
  ratio = ifelse(is_MCMC,
                 min(1, exp(num-den+
                              logPriorsRatio$likeliParam(LP_new, LP_old, prior_par$likeliParam))),
                 min(1, exp(num-den+
                              alpha*logPriorsRatio$likeliParam(LP_new, LP_old, prior_par$likeliParam)+
                              (1-alpha)*logReferenceRatio$likeliParam(LP_new, LP_old, reference_par$likeliParam))))

  # accept-reject step
  U <- runif(1)
  return( U<ratio )
}
