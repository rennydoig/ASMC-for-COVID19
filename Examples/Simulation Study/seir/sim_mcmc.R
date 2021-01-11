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

sim_data <- readRDS("../physdist/sim_data.rds")
true <- sim_data$true
data <- sim_data$sim

# Data simulation --------------------------------------------------------------

# physical distancing transmission model
SEIRmodel = function(t, state, parms){
  with(as.list(c(state, parms)), {
    R0 <- exp(logR0)
    beta <- R0 * gamma
    
    dSdt <- - beta * I * S / N
    dEdt <- beta * I * S / N - delta*E
    dIdt <- delta*E - gamma*I
    
    return( list(c(dSdt, dEdt, dIdt)) )
  })
}

# true pars
N <- 5.1e6
pars <- list(N = N,
             logR0 = log(2.5),
             delta = 1/3,
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
is_unknownPar <- c(F, T, F, F)

mean_theta <- c(log(2.5))
spread_theta <- c(0.1)

prior_par$theta <- list(mean=mean_theta, sd=spread_theta)
logPriorsRatio$theta = function(xstar, x, par){
  with(par, {
    dxstar = sum(dnorm(xstar, mean, sd, log=T))
    dx = sum(dnorm(x, mean, sd, log=T))
    return( dxstar - dx )
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
    dx = sum(dnorm(x, mean, sd, log=T))
    return(dx)
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
K <- 18000
hyperpar <- list(sigma_theta = c(0, 0.05, 0, 0),
                 sigma_likeliParam = 0.5)

if(file.exists("fit_se_mcmc.rds")){fits <- readRDS("fit_se_mcmc.rds")}else{
  set.seed(51)
  plan(multisession(workers=5))
  fits <- future_map(data[1:30], ~ MCMC2_LP(K, ., is_unknownPar, SEIRmodel, unlist(pars), likelihood_NB,
                                            unlist(Lpars[!is_unknownLPar]), unlist(Lpars[is_unknownLPar]), prior_list, hyperpar),
                     .options=furrr_options(seed=T))
  plan(sequential)
  saveRDS(fits, "fit_se_mcmc.rds")
}


# Estimate marginal likelihood -------------------------------------------------
lb <- rep(-Inf, sum(is_unknownPar) + 1)
ub <- rep(Inf, sum(is_unknownPar) + 1)
names(ub) <- names(lb) <- c(names(pars)[is_unknownPar], names(Lpars)[is_unknownLPar])

plan(multisesstion(workers=5))
set.seed(52)
ml_mcmc <- future_map(1:30, ~ getMarginalLikelihood2_LP(fits[[.]]$sample_list, 17000, is_unknownPar, data[[.]], SEIRmodel,
                                                       likelihood_NB, Lpars[!is_unknownLPar], hyperpar, prior_list, lb, ub, silent=T),
                      .progress=T, .options=furrr_options(seed=T))
plan(sequential)

saveRDS(ml_mcmc, "ml_mcmc.rds")