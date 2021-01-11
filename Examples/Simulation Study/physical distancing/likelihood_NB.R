####################################################
# # Likelihood model for the observations: Poisson distribution
# # No parameters in the likelihood functions
likelihood_NB <- list()
likelihood_NB$dFun <- function(data, fit, ODEpars, likeliParam, log=T)
{
  with(as.list(likeliParam), {
    psi = exp(logpsi)
    
    # t=0 corresponds to 2020-03-03
    # 2020-03-14 (t=11) is the change point in the sampling proportion
    p = ifelse(fit[,"time"] < 11, p1, p2)

    mu = p * (fit[,"I"] + fit[,"Id"])

    # this likelihood assumes only I is observed, so there is only one element of is_obs=T
    temp = dnbinom(data$y, size=psi, mu=mu, log=log)
    return(sum(temp[!is.nan(temp)]))
  })
}
likelihood_NB$getMean <- function(fit, pars)
{
  with(as.list(pars), {
    # t=0 corresponds to 2020-03-03
    # 2020-03-14 (t=11) is the change point in the sampling proportion
    p = ifelse(fit[,"time"] < 11, p1, p2)
    mu = p * (fit[,"I"] + fit[,"Id"])
    return(data.frame(t=fit[,"time"], mu=mu))
  })
}
