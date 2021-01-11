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

# Data simulation --------------------------------------------------------------

# physical distancing transmission model
PhysDistmodel = function(t, state, parms){
  with(as.list(c(state, parms)), {
    # inverse-transform estimated parameters
    R0 = exp(logR0)
    q = exp(logq)
    ud = exp(logud)
    f1 <- exp(logitf1)/(1+exp(logitf1))
    f2 <- exp(logitf2)/(1+exp(logitf2))
    
    beta = R0/(1/(q+1/D) + 1/k2)
    f = ifelse(t<tau1, 1, ifelse(t<tau2, f1, f2))
    
    dSdt = -beta*(I+E2+f*(Id+E2d))*S/N - ud*S + ur*Sd
    dE1dt = beta*(I+E2+f*(Id+E2d))*S/N - k1*E1 - ud*E1 + ur*E1d
    dE2dt = k1*E1 - k2*E2 - ud*E2 + ur*E2d
    dIdt = k2*E2 - q*I - I/D - ud*I + ur*Id
    dQdt = q*I - Q/D - ud*Q + ur*Qd
    dRdt = I/D + Q/D - ud*R + ur*Rd
    
    dSddt = -f*beta*(I+E2+f*(Id+E2d))*Sd/N + ud*S - ur*Sd
    dE1ddt = f*beta*(I+E2+f*(Id+E2d))*Sd/N - k1*E1d + ud*E1 - ur*E1d
    dE2ddt = k1*E1d - k2*E2d + ud*E2 - ur*E2d
    dIddt = k2*E2d - q*Id - Id/D + ud*I - ur*Id
    dQddt = q*Id - Qd/D + ud*Q - ur*Qd
    dRddt = Id/D + Qd/D + ud*R - ur*Rd
    
    return( list(c(dSdt, dE1dt, dE2dt, dIdt, dQdt, dRdt, dSddt, dE1ddt, dE2ddt, dIddt, dQddt, dRddt)) )
  })
}

# true pars
N <- 5.1e6
pars <- list(N = N,
             logR0 = log(2.57),
             k1 = 1/5,
             k2 = 1,
             logq = log(0.05),
             D = 5,
             logud = log(0.1),
             ur = 0.02,
             tau1 = 16,
             tau2 = 75,
             logitf1 = log(0.35/0.65),
             logitf2 = log(0.65/0.35))

# initial states
# assume 8 active cases on February 8th
i0 <- 8
e <- with(pars, {ur/(exp(logud) + ur)})
initial_states <- c(S = (1-e) * (N - i0),
                    E1 = 0.4 * (1-e) * i0,
                    E2 = 0.1 * (1-e) * i0,
                    I = 0.5 * (1-e) * i0,
                    Q = 0,
                    R = 0,
                    Sd = e * (N - i0),
                    E1d = 0.4 * e * i0,
                    E2d = 0.1 * e * i0,
                    Id = 0.5 * e * i0,
                    Qd = 0,
                    Rd = 0)

# likelihood parameters
Lpars <- list(p1 = 0.1,
              p2 = 0.3,
              logpsi = log(5))

# simulation settings
n_samples <- 100
times <- seq(-30, 100, by=0.1)
M <- 50

# get simulated data
set.seed(50)

if(file.exists("sim_data.rds")){sim <- readRDS("sim_data.rds")}else{
  sim <- multiSim(M, n_samples, likelihood_NB, PhysDistmodel, initial_states, times, c(pars, Lpars))
  saveRDS(sim, "sim_data.rds")
}

truth <- as.data.frame(sim$true)
data <- sim$sim


# Set prior and reference distributions ----------------------------------------

prior_par <- list()
logPriorsRatio <- list()

## DE parameters
is_unknownPar <- rep(F, length(pars))
is_unknownPar[c(2,5,7,11,12)] <- T
mean_theta <- c(log(2.5), log(0.075), log(0.15), 0.4, 0.7)
spread_theta <- c(0.5, 0.5, 0.5, 5, 5)
prior_par$theta <- list(mean=mean_theta, sd=spread_theta)
logPriorsRatio$theta = function(xstar, x, par){
  with(par, {
    # extract the estimated parameters & means
    thstar <- xstar
    th <- x
    
    thstar[4:5] = exp(thstar[4:5])/(1+exp(thstar[4:5]))
    th[4:5] = exp(th[4:5])/(1+exp(th[4:5]))
    
    dxstar = sum(dnorm(thstar[1:3], mean[1:3], sd[1:3], log=T)) +
      sum(dbeta(thstar[4:5], mean[4:5]*sd[4:5], (1-mean[4:5])*sd[4:5], log=T))
    
    dx = sum(dnorm(th[1:3], mean[1:3], sd[1:3], log=T)) +
      sum(dbeta(th[4:5], mean[4:5]*sd[4:5], (1-mean[4:5])*sd[4:5], log=T))
    
    logJ = sum(xstar[4:5]) - 2*sum(log(1+exp(xstar[4:5]))) - sum(x[4:5]) + 2*sum(log(1+exp(x[4:5])))
    
    return( dxstar - dx + logJ )
  })
}


## likelihood parameters
is_unknownLpar <- c(F, F, T)
mean_eta <- 2
spread_eta <- 5
prior_par$likeliParam <- list(mean=mean_eta, sd=spread_eta)
logPriorsRatio$likeliParam = function(xstar, x, par){
  with(par, {
    # transform parameters
    psistar = exp(xstar)
    psi = exp(x)
    
    # evaluate prior densities
    dxstar = dgamma(psistar, shape=mean^2/sd, rate=mean/sd, log=T)
    dx = dgamma(psi, shape=mean^2/sd, rate=mean/sd, log=T)
    
    # log Jacobian
    logJ = xstar - x
    
    return( dxstar - dx + logJ )
  })
}

dprior <- list()
dprior$theta <- function(x, par){
  with(par, {
    is.logit <- 4:5
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
K <- 55000

hyperpar <- list(sigma_theta = c(0, 0.025, 0, 0, 0.45, 0, 0.1, 0, 0, 0, 0.15, 0.75),
                 sigma_likeliParam = 0.5)

set.seed(51)

if(file.exists("fit_pd_mcmc.rds")){ fits <- readRDS("fit_pd_mcmc.rds")}else{
  plan(multisession(workers=15))
  fits <- future_map(data[1:30], ~ MCMC2_LP(K, ., is_unknownPar, PhysDistmodel, unlist(pars), likelihood_NB,
                                            unlist(Lpars[!is_unknownLpar]), unlist(Lpars[is_unknownLpar]), prior_list, hyperpar),
                     .options=furrr_options(seed=T))
  plan(sequential)
  saveRDS(fits, "fit_pd_mcmc.rds")
}


# Estimate marginal likelihood -------------------------------------------------
lb <- rep(-Inf, sum(is_unknownPar) + 1)
ub <- rep(Inf, sum(is_unknownPar) + 1)
names(ub) <- names(lb) <- c(names(pars)[is_unknownPar], names(Lpars)[is_unknownLpar])

plan(multisession(workers=5))
set.seed(52)
ml_mcmc <- future_map(1:30, ~ getMarginalLikelihood2_LP(fits[[.]]$sample_list, 54000, is_unknownPar, data[[.]], PhysDistmodel,
                                                       likelihood_NB, Lpars[!is_unknownLpar], hyperpar, prior_list, lb, ub, silent=T),
                      .progress=T, .options=future_options(seed=T))
plan(sequential)

saveRDS(ml_mcmc, "ml_mcmc.rds")