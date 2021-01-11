#' Estimate marginal likelihood from MCMC samples
#'
#' Uses bridge sampling to estimate the marginal likelihood of the data from the MCMC samples. Uses the bridgesampling package.
#'
#' @param result_list List of same format as the sample_list compoenent of the MCMC_SIR output
#' @param burn_in Number of samples to use as burn_in
#' @param is_unknownPar Boolean vector indicating which DE parameters are to be estimated
#' @param data List of data
#' @param basis_funs List of basis function expansion components matching basisFuns output
#' @param ODEmodel ODE model
#' @param likelihood List containing likelihood density/mass
#' @param likeliParam_fixed Fixed parameters to be passed to the likelihood. Leave as NULL if all additional parameters are estimated.
#' @param hyperpar List of hyperparameters; same as was used to generate the MCMC output
#' @param prior_list List of prior densities for theta and basis function coefficients
#' @param lower_bound Named vector of lower bounds on estimated parameters to be passed to bridge_sampler
#' @param upper_bound Named vector of upper bounds on estimated paramters to be passed to bridge_sampler
#' @param ... Additional parameters to be passed to bridge_sampler
#'
#' @return Estimated marginal likelihood of the data
#'
#' @export
getMarginalLikelihood2_LP = function(result_list, burn_in, is_unknownPar, data, ODEmodel, likelihood, likeliParam_fixed=NULL, hyperpar, prior_list, lower_bound, upper_bound, ...){
  with( c(data), {
    # extract the number of DE parameters, number of unknown parameters, DE param names, fixed DE parameters
    n_DEpar = length(result_list[[1]]$theta)
    n_unknownPar = sum(is_unknownPar)
    DE_par_names = names(result_list[[1]]$theta)
    theta_known = result_list[[1]]$theta[!is_unknownPar]

    n_LP = length(result_list[[1]]$likeliParam)
    LP_names = names(result_list[[1]]$likeliParam)

    # convert result_list into a matrix with the ith row being the ith sample
    theta = matrix(unlist(lapply(result_list, `[[`, 'theta')), ncol=n_DEpar, byrow=T)[,is_unknownPar]
    likeliParam = matrix(unlist(lapply(result_list, `[[`, "likeliParam")), ncol=n_LP, byrow=T)

    samples = cbind(theta, likeliParam)[-(1:burn_in),]
    colnames(samples) = c(DE_par_names[is_unknownPar], LP_names)

    # function for evaluating the log-posterior
    # takes in a row from the samples matrix
    posterior = function(samples, data, theta_known, is_unknownPar, ODEmodel, likelihood, likeliParam_fixed, hyperpar, prior_list){
      with( c(data, hyperpar, prior_list), {
        # extract DE parameters from the sample matrix
        n_unknownPar = sum(is_unknownPar)
        theta = rep(NA, length(is_unknownPar))
        theta[is_unknownPar] = samples[1:n_unknownPar]
        names(theta)[is_unknownPar] = names(samples)[1:n_unknownPar]
        theta[!is_unknownPar] = theta_known
        names(theta)[!is_unknownPar] = names(theta_known)

        # extract likelihood parameters
        likeliParam = samples[-(1:n_unknownPar)]
        names(likeliParam) = names(samples)[-(1:n_unknownPar)]

        # function for evaluating the likelihood of a single sample
        eval_likelihood = logLik2_LP(ODEmodel, theta, data, likelihood, likeliParam_fixed, likeliParam, 1)

        # evaluate the prior of theta
        ind = prior_par$theta$index_unknown
        theta_prior = sum(dprior$theta(theta, prior_par$theta))
        LP_prior = sum(dprior$likeliParam(likeliParam, prior_par$likeliParam))

        return( eval_likelihood+theta_prior+LP_prior)
      })
    }

    # use the bridge_sampler function to estimate the marginal likelihood
    marginal_likelihood = bridge_sampler( samples = samples,
                                          log_posterior = posterior,
                                          theta_known=theta_known, is_unknownPar=is_unknownPar,  ODEmodel=ODEmodel, likelihood=likelihood, likeliParam_fixed=likeliParam_fixed, hyperpar=hyperpar, prior_list=prior_list,
                                          data = data,
                                          lb = lb,
                                          ub = ub,
                                          ... )

    return( marginal_likelihood$logml)
  })
}
