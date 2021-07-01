# Annealed SMC for COVID-19
# Renny Doig
# November 14, 2020
#
# - Generate 50 datasets from physical distancing transmission model with
#   negative binomial likelihood
# - Run ASMC on assuming simple physical distancing model

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


## reference distributions
reference_sim <- list()
reference_sim$theta <- function(par){
  with(par, {
    temp <- c(rnorm(1, mean[1], sd[1]),
              rbeta(2, shape1=mean[2:3]*sd[2:3], shape2=(1-mean[2:3])*sd[2:3]))
    temp[2:3] <- log(temp[2:3]/(1-temp[2:3]))
    
    param[index] <- temp
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
hyperpar <- list(sigma_theta1 = c(0, 0.025, 0, 0, 0, 0.1, 0.6, 0),
                 sigam_theta2 = NA,
                 ksi = 1,
                 sigma_likeliParam = 0.5)

set.seed(51)

if(file.exists("fit_si_asmc.rds")){fits <- readRDS("fit_si_asmc.rds")}else{
  fits <- map(data, ~ ASMC(K, smc_param, ., is_unknownPar, Simplemodel, pars, likelihood_NB,
                   Lpars[!is_unknownLPar], Lpars[is_unknownLPar], reference_prior_list, hyperpar, ncore))
  saveRDS(fits, "fit_si_asmc.rds")
}
