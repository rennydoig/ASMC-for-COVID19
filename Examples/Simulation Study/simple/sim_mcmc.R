# Annealed SMC for COVID-19
# Renny Doig
# November 5, 2020
#
# Repeat ASMC simulation but with MCMC
# - run targeting the same run-time as the SMC simulation

# Preliminaries ----------------------------------------------------------------
rm(list=ls())
library(annealedSMC)
library(purrr)
library(deSolve)
library(furrr)

source("likelihood_NB.R")

sim_data <- readRDS("../sim_data.rds")
true <- sim_data$true
data <- sim_data$sim

# Data simulation --------------------------------------------------------------

# physical distancing transmission model
Simplemodel = function(t, state, parms){
  with(as.list(c(state, parms)), {
    R0 <- exp(logR0)
    f1 <- exp(logitf1) / (1 + exp(logitf1))
    f2 <- exp(logitf2) / (1 + exp(logitf2))
    
    beta <- R0 * gamma
    f <- ifelse(t<tau1, 1, ifelse(t<tau2, f1, f2))
    
    dSdt <- -f * beta * I * S / N
    dEdt <- f * beta * I * S / N - delta*E
    dIdt <- delta*E - gamma*I
    
    return( list(c(dSdt, dEdt, dIdt)) )
  })
}

# true pars
N <- 5.1e6
pars <- list(N = N,
             logR0 = log(2.5),
             delta = 1/3,
             tau1 = 16,
             tau2 = 74,
             logitf1 = log(0.4/0.6),
             logitf2 = log(0.7/0.3),
             gamma = 1/5)


# initial states
# assume 8 active cases on February 8th
i0 = 8
initial_states <- c(S = N - i0,
                    E = 0.4 * i0,
                    I = 0.6 * i0)
data <- map(data, function(.x){
  .x$initial_state <- initial_states
  .x
})


# Prior distributions ----------------------------------------------------------

prior_par <- list()
logPriorsRatio <- list()


## DE parameters
is_unknownPar <- c(F, T, F, F, F, T, T, F)

mean_theta <- c(log(2.5), 0.4, 0.7)
spread_theta <- c(0.1, 20, 20)

prior_par$theta <- list(mean=mean_theta, sd=spread_theta)
logPriorsRatio$theta = function(xstar, x, par){
  with(par, {
    is.logit <- 2:3
    # extract the estimated parameters & means
    thstar <- xstar
    th <- x
    
    thstar[is.logit] = exp(thstar[is.logit])/(1+exp(thstar[is.logit]))
    th[is.logit] = exp(th[is.logit])/(1+exp(th[is.logit]))
    
    dxstar = sum(dnorm(thstar[!is.logit], mean[!is.logit], sd[!is.logit], log=T)) +
      sum(dbeta(thstar[is.logit], mean[is.logit]*sd[is.logit], (1-mean[is.logit])*sd[is.logit], log=T))
    
    dx = sum(dnorm(th[!is.logit], mean[!is.logit], sd[!is.logit], log=T)) +
      sum(dbeta(th[is.logit], mean[is.logit]*sd[is.logit], (1-mean[is.logit])*sd[is.logit], log=T))
    
    logJ = sum(xstar[is.logit]) - 2*sum(log(1+exp(xstar[is.logit]))) - sum(x[is.logit]) + 2*sum(log(1+exp(x[is.logit])))
    
    return( dxstar - dx + logJ )
  })
}



## likelihood parameters
Lpars <- list(p1 = 0.1,
              p2 = 0.3,
              logpsi = log(1))
is_unknownLPar <- c(F, F, T)

mean_eta <- 1
spread_eta <- 3
prior_par$likeliParam <- list(mean=mean_eta, sd=spread_eta)
logPriorsRatio$likeliParam = function(xstar, x, par){
  # transform parameters
  psistar = exp(xstar)
  psi = exp(x)
  
  # evaluate prior densities
  dxstar = dgamma(psistar, shape=par$mean^2/par$sd, rate=par$mean/par$sd, log=T)
  dx = dgamma(psi, shape=par$mean^2/par$sd, rate=par$mean/par$sd, log=T)
  
  # log Jacobian
  logJ = xstar - x
  
  return( dxstar - dx + logJ )
}

dprior <- list()
dprior$theta <- function(x, par){
  with(par, {
    is.logit <- 2:3
    # extract the estimated parameters & means
    th <- x
    
    th[is.logit] = exp(th[is.logit])/(1+exp(th[is.logit]))
    
    dx = sum(dnorm(th[!is.logit], mean[!is.logit], sd[!is.logit], log=T)) +
      sum(dbeta(th[is.logit], mean[is.logit]*sd[is.logit], (1-mean[is.logit])*sd[is.logit], log=T))
    
    logJ = sum(x[is.logit]) - 2*sum(log(1+exp(x[is.logit])))
    
    return(dx + logJ)
  })
}
dprior$likeliParam <- function(x, par){
  with(par, {
    psi <- exp(x)
    dx <- dgamma(psi, shape=mean^2/sd, rate=mean/sd) + x
    return(dx)
  })
}

prior_list <- list(prior_par = prior_par,
                   logPriorsRatio = logPriorsRatio,
                   dprior = dprior)
                        


# Running ASMC -----------------------------------------------------------------

## algorithm settings
K <- 23000
hyperpar <- list(sigma_theta = c(0, 0.025, 0, 0, 0, 0.1, 0.6, 0),
                 sigma_likeliParam = 0.5)

if(file.exists("fit_si_mcmc.rds")){fits <- readRDS("fit_si_mcmc.rds")}else{
  set.seed(51)
  plan(multisession(workers=5))
  fits <- future_map(data[1:30], ~ MCMC(K, ., is_unknownPar, Simplemodel, unlist(pars), likelihood_NB,
                                        unlist(Lpars[!is_unknownLPar]), unlist(Lpars[is_unknownLPar]), prior_list, hyperpar),
                     .options=furrr_options(seed=T))
  plan(sequential)
  saveRDS(fits, "fit_si_mcmc.rds")
}


# Estimate marginal likelihood -------------------------------------------------
lb <- rep(-Inf, sum(is_unknownPar) + 1)
ub <- rep(Inf, sum(is_unknownPar) + 1)
names(ub) <- names(lb) <- c(names(pars)[is_unknownPar], names(Lpars)[is_unknownLPar])

plan(multisesstion(workers=5))
set.seed(52)
ml_mcmc <- future_map(1:30, ~ getMarginalLikelihood(fits[[.]]$sample_list, 23000/2, is_unknownPar, data[[.]], Simplemodel,
                                                    likelihood_NB, Lpars[!is_unknownLPar], hyperpar, prior_list, lb, ub, silent=T),
                      .options=furrr_options(seed=T))
plan(sequential)

saveRDS(ml_mcmc, "ml_mcmc.rds")