# Annealed SMC for COVID-19
# Renny Doig
# November 14, 2020
#
# - Generate 50 datasets from physical distancing transmission model with
#   negative binomial likelihood
# - Run ASMC on assuming SEIR model

# Preliminaries ----------------------------------------------------------------
rm(list=ls())
library(annealedSMC)
library(purrr)
library(deSolve)
library(tidyr)
library(tictoc)

source("likelihood_NB.R")

sim_data <- readRDS("../sim_data.rds")
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


## reference distributions
reference_sim <- list()
reference_sim$theta <- function(par){
  with(par, {
    param[index] <- rnorm(1, mean, sd)
    return(param)
  })
}
reference_sim$likeliParam <- function(par){
  log(rgamma(1, shape=par$mean^2/par$sd, rate=par$mean/par$sd))
}

reference_par <- prior_par
reference_par$theta$index <- is_unknownPar
reference_par$theta$param <- unlist(pars)

logReferenceRatio <- logPriorsRatio
logPriorReference <- list()
logPriorReference$theta <- function(x, y, z){ 0 }
logPriorReference$likeliParam <- function(x, y, z){ 0 }

# list containing all objects necessary for prior and reference distribution evaluation
reference_prior_list = list(prior_par = prior_par,
                            reference_par = reference_par,
                            reference_sim = reference_sim,
                            logPriorsRatio = logPriorsRatio,
                            logReferenceRatio = logReferenceRatio,
                            logPriorReference = logPriorReference)



# Running ASMC -----------------------------------------------------------------

## algorithm settings
K <- 1000
ncore <- 10
smc_param <- list(eps=0.5, phi=0.99, n_move=1)
hyperpar <- list(sigma_theta1 = c(0, 0.05, 0, 0),
                 sigma_theta2 = NA,
                 ksi = 1,
                 sigma_likeliParam = 0.5)

set.seed(51)

if(file.exists("fit_si_asmc.rds")){fits <- readRDS("fit_si_asmc.rds")}else{
  fits <- map(data, ~ ASMC(K, smc_param, ., is_unknownPar, SEIRmodel, pars, likelihood_NB,
                   Lpars[!is_unknownLPar], Lpars[is_unknownLPar], reference_prior_list, hyperpar, ncore))
  saveRDS(fits, "fit_si_asmc.rds")
}
