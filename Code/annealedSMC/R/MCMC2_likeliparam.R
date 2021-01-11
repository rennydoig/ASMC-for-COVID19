#' MCMC for SIR model
#'
#' Markov chain Monte Carlo for estimating parameters of susceptible-infected-recovered model with or without additional likelihood parameters.
#'
#' @param K Number of samples to generate
#' @param data List of observed data
#' @param is_unknownPar Boolean vector indicating which parameters are to be estimated
#' @param ODEmodel Function defining the DE system. Must be compatible with deSolve.
#' @param init_DEpar Initial estimates for the DE parameters
#' @param likelihood List containing funcitons for simulating from and evaluating the likelihood of the data.
#' @param likeliParam_fixed Vector of fixed likelihood parameters.
#' @param init_likeliParam_unknown Vector of initial values for the unknown likeliood parameters.
#' @param prior_list List of prior distributions and parameters
#' @param hyperpar List of hyperparameters
#'
#' @return List containing posterior distribution of DE parameters and basis coefficients
#'
#' @export
MCMC2_LP <- function(K, data, is_unknownPar, ODEmodel, init_DEpar,
                     likelihood, likeliParam_fixed, init_likeliParam_unknown=NULL,
                     prior_list, hyperpar) {
  # if there are no likelihood parameters to be estimated, add dummy components to prior & reference distributions
  if( is.null(init_likeliParam_unknown) )
  {
    prior_list$prior_par$likeliParam = NA
    prior_list$logPriorsRatio$likeliParam = function(x,y,z){0}
  }

  with(as.list(c(hyperpar, data, prior_list)),{

    # initialize values
    n_param = length(init_DEpar)
    sample_list = list()

    # initialize sample_list
    sample_list[[1]] = list(theta = init_DEpar,
                            likeliParam = init_likeliParam_unknown)

    # used for keeping track of which and how many DE parameters are being estimated
    index_unknownPar = which(is_unknownPar==TRUE)
    n_unknownPar = length(index_unknownPar)

    # initialize counters for theta and c acceptance
    n_accept_theta = rep(0, n_unknownPar)
    n_sample_theta = rep(0, n_unknownPar)
    n_accept_LP = 0

    for( r in 2:K )
    {
      if( (r%%100)==0 )
        cat(r," ")

      sample_list[[r]] = sample_list[[r-1]]

      # sample the DE parameters
      new_theta = sample_list[[r-1]]$theta
      rand_num = sample(1:n_unknownPar,1)
      n_sample_theta[rand_num] = n_sample_theta[rand_num] + 1
      new_theta[index_unknownPar[rand_num]] = rnorm(1, sample_list[[r-1]]$theta[index_unknownPar[rand_num]], sd=sigma_theta[index_unknownPar[rand_num]])

      if(MH_theta2_LP(new_theta,
                      ODEmodel,
                      sample_list[[r-1]]$theta,
                      data,
                      likelihood,
                      likeliParam_fixed,
                      sample_list[[r-1]]$likeliParam,
                      1,
                      prior_par,
                      reference_par,
                      logPriorsRatio,
                      logReferenceRatio,
                      is_MCMC = T))
      {
        sample_list[[r]]$theta = new_theta
        n_accept_theta[rand_num] = n_accept_theta[rand_num]+1
      }

      if( is.null(init_likeliParam_unknown)){
        sample_list[[r]]$likeliParam = NULL
      }else{
        ## Sample the likelihood parameter(s)
        new_LP = rnorm(1, sample_list[[r-1]]$likeliParam, sigma_likeliParam)
        names(new_LP) = names(init_likeliParam_unknown)

        if( MH_likeliParam2(new_LP,
                            sample_list[[r-1]]$likeliParam,
                            ODEmodel,
                            sample_list[[r]]$theta,
                            data,
                            likelihood,
                            likeliParam_fixed,
                            1,
                            prior_par,
                            reference_par,
                            logPriorsRatio,
                            logReferenceRatio,
                            T))
        {
          sample_list[[r]]$likeliParam = new_LP
          n_accept_LP = n_accept_LP + 1
        }
      }
    }

    return( list(sample_list=sample_list,n_accept_theta=n_accept_theta/n_sample_theta, n_accept_likeliParam=n_accept_LP/K) )
  })
}
