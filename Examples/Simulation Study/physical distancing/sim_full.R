# Annealed SMC for COVID-19
# Renny Doig
# November 14, 2020
#
# - Generate 50 datasets from physical distancing transmission model with
#   negative binomial likelihood
# - Run ASMC on all 10 datasets using the "high" prior variance from simstudy2

# Preliminaries ----------------------------------------------------------------
rm(list=ls())
library(annealedSMC)
library(purrr)
library(deSolve)
library(tidyr)
library(tictoc)

library(tidyverse)
library(cowplot)

source("functions.R")
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

if(file.exists("../sim_data.rds")){sim <- readRDS("../sim_data.rds")}else{
  sim <- multiSim(M, n_samples, likelihood_NB, PhysDistmodel, initial_states, times, c(pars, Lpars))
  saveRDS(sim, "../sim_data.rds")
}

truth <- as.data.frame(sim$true)
data <- sim$sim


# Simulated data ---------------------------------------------------------------

## check the true and simulated trajectories to make sure the results are reasonable

as.data.frame(truth) %>%
  mutate(active = I + Id) %>%
  ggplot(aes(x=time, y=active)) + geom_line() +
  geom_line(data=map_dfr(data, ~ data.frame(time=.$t, active=.$y), .id="dataset"),
            aes(group=dataset), size=0.1)


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

## Reference distribution settings
reference_par <- prior_par
reference_par$theta$index <- is_unknownPar
reference_par$theta$param <- unlist(pars)
reference_sim <- list()
reference_sim$theta <- function(par){
  with(par, {
    temp <- c(rnorm(3, mean[1:3], sd[1:3]),
              rbeta(2, shape1=mean[4:5]*sd[4:5], shape2=(1-mean[4:5])*sd[4:5]))
    temp[4:5] <- log(temp[4:5]/(1-temp[4:5]))
    
    param[index] <- temp
    return(param)
  })
}
reference_sim$likeliParam <- function(par){
  log(rgamma(1, shape=par$mean^2/par$sd, rate=par$mean/par$sd))
}
logPriorReference <- list()
logPriorReference$theta <- function(x, y, z){ 0 }
logPriorReference$likeliParam <- function(x, y, z){ 0 }

prior_reference_list <- list(prior_par = prior_par,
                             logPriorsRatio = logPriorsRatio,
                             reference_par = reference_par,
                             reference_sim = reference_sim,
                             logReferenceRatio = logPriorsRatio,
                             logPriorReference = logPriorReference)



# Running ASMC -----------------------------------------------------------------

## algorithm settings
K <- 1000
ncore <- 10
smc_param <- list(eps=0.5, phi=0.99, n_move=1)
hyperpar <- list(sigma_theta1 = c(0, 0.025, 0, 0, 0.45, 0, 0.1, 0, 0, 0, 0.15, 0.75),
                 sigam_theta2 = NA,
                 ksi = 1,
                 sigma_likeliParam = 0.5)

set.seed(51)

if(file.exists("fit_pd_asmc.rds")){ fits_asmc <- readRDS("fit_pd_asmc.rds")}else{
  fits <- map(data[1:10], ~ ASMC(K, smc_param, ., is_unknownPar, PhysDistmodel, pars, likelihood_NB,
                                Lpars[!is_unknownLpar], Lpars[is_unknownLpar], prior_reference_list, hyperpar, ncore))
  saveRDS(fits, "fit_pd_asmc.rds")
}

fits_mcmc <- readRDS("fit_pd_mcmc.rds")


# Parameter point estimates ----------------------------------------------------

resample_seeds <- 1996 + 1:M

true_pars <- data.frame(Parameter = c("R0", "q", "ud", "f1", "f2", "psi"),
                        truth = c(2.57, 0.05, 0.1, 0.35, 0.65, 5))

theta <- map(fits_asmc, function(.x){
  map_dfr(.x$particle_list, "theta") %>%
    bind_cols(map_dfr(.x$particle_list, "likeliParam")) %>%
    mutate(R0 = exp(logR0),
           q = exp(logq),
           ud = exp(logud),
           f1 = exp(logitf1) / (1 + exp(logitf1)),
           f2 = exp(logitf2) / (1 + exp(logitf2)), 
           psi = exp(logpsi), .keep="none") %>%
    bind_cols(W = .x$W)
})

CIs <- map_dfr(c("R0", "q", "ud", "f1", "f2", "psi"), function(.par){
  map_dfr(theta, ~ get_CI(., .par), .id="dataset")
}) %>%
  full_join(true_pars) %>%
  mutate(dataset=as.numeric(dataset))

filter(CIs, dataset<=10) %>%
  ggplot(aes(x=mean, y=NA, group=dataset)) + geom_point(position=position_dodge(0.6)) +
  geom_linerange(aes(xmin=lower.bound, xmax=upper.bound), position=position_dodge(0.6), size=1) +
  facet_wrap(~ Parameter, scales="free_x") +
  geom_vline(aes(xintercept=truth), colour="red") +
  theme(panel.background = element_rect(fill="white", colour="black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  background_grid(major="x", colour.major="grey90") +
  labs(x="", y="")
ggsave("credible_bands.pdf", width=6, height=4, units="in", dpi=900)

# table of means and their standard errors
group_by(CIs, Parameter) %>%
  summarise(p.mean=mean(mean), se=sd(mean)/sqrt(50), .groups="keep")

# coverage probabilities
mutate(CIs, in95=(truth>=lower.bound & truth<=upper.bound)) %>%
  group_by(Parameter) %>%
  summarise(lower=mean(lower.bound), upper=mean(upper.bound), CP=mean(in95))


# histogram
theta_long.1 <- get_theta_long(fits_asmc[[1]], resample_seed)

ggplot(theta_long.1, aes(x=value)) + geom_histogram(bins=50) +
  facet_wrap(~Parameter, scales="free")


## Compare against MCMC
K_mcmc <- 1:(length(fits_mcmc[[1]]$sample_list) / 2)
theta_mcmc <- map_dfr(fits_mcmc, function(.x){
  map_dfr(.x$sample_list[-K_mcmc], "theta") %>%
    bind_cols(map_dfr(.x$sample_list[-K_mcmc], "likeliParam")) %>%
    mutate(R0 = exp(logR0),
           q = exp(logq),
           ud = exp(logud),
           f1 = exp(logitf1) / (1 + exp(logitf1)),
           f2 = exp(logitf2) / (1 + exp(logitf2)), 
           psi = exp(logpsi), .keep="none")
}, .id="dataset")

CIs_mcmc <- pivot_longer(theta_mcmc, cols=-dataset, names_to="Parameter", values_to="value") %>%
  mutate(dataset=as.numeric(dataset)) %>%
  group_by(dataset, Parameter) %>%
  summarise(mean=mean(value),
            lower.bound=quantile(value, probs=0.025),
            upper.bound=quantile(value, probs=0.975), .groups="keep") %>%
  full_join(true_pars)

bind_rows(cbind(filter(CIs, dataset<=10), Method="ASMC"),
          cbind(filter(CIs_mcmc,dataset<=10), Method="MCMC")) %>%
  ggplot(aes(x=mean, y=Method, group=dataset)) + 
  geom_point(position = position_dodge(0.6)) +
  geom_linerange(aes(xmin=lower.bound, xmax=upper.bound), position=position_dodge(0.6)) +
  facet_wrap(~ Parameter, scales="free_x") +
  geom_vline(aes(xintercept=truth), colour="red") +
  theme(panel.background = element_rect(fill="white", colour="black"),
        axis.ticks.y = element_blank()) +
  background_grid(major="x", colour.major="grey90") +
  labs(x="")
ggsave("credible_bands_comparison.pdf", width=6, height=3.5, units="in", dpi=900)



# Coverage probabilities -------------------------------------------------------

cred_levels <- c(seq(90, 10, by=-10)/100, 0.05, 0.01)

CP <- map_dfr(cred_levels, function(.level){
  map_dfr(c("R0", "q", "ud", "f1", "f2", "psi"), function(.par){
    map_dfr(theta, ~ get_CI(., .par, .level))
  })  %>%
    full_join(true_pars) %>%
    mutate(inCI = ((truth>=lower.bound) & (truth<=upper.bound))) %>%
    group_by(Parameter) %>%
    summarise(CP = mean(inCI), .groups="keep") %>%
    bind_cols(Level=1-.level)
})

CP_mcmc <- map_dfr(cred_levels, function(.level){
  pivot_longer(theta_mcmc, cols=-dataset, names_to="Parameter", values_to="value") %>%
    mutate(dataset=as.numeric(dataset)) %>%
    group_by(dataset, Parameter) %>%
    summarise(mean=mean(value),
              lower.bound=quantile(value, probs=.level/2),
              upper.bound=quantile(value, probs=1-.level/2), .groups="keep") %>%
    full_join(true_pars) %>%
    mutate(inCI = ((truth>=lower.bound) & (truth<=upper.bound))) %>%
    group_by(Parameter) %>%
    summarise(CP = mean(inCI), .groups="keep") %>%
    bind_cols(Level=1-.level)
})

bind_rows(cbind(CP, Method="ASMC"), cbind(CP_mcmc, Method="MCMC")) %>%
  ggplot(aes(x=Level, y=CP, shape=Method)) + geom_point() + geom_abline(slope=1) +
  facet_wrap(~ Parameter) +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="xy", colour.major="grey90") +
  labs(x="Credible level", y="Estimated coverage probability")
ggsave("coverage_probability.pdf", width=6, height=4, units="in", dpi=900)

# Active cases -----------------------------------------------------------------

# first 10 only
if(file.exists("active_cases.rds")){active_cases <- readRDS("active_cases.rds")}else{
  active_cases <- map_dfr(1:10, ~ get_active_cases(fits_asmc[[.]], times, initial_states, PhysDistmodel, resample_seeds[.]), .id="dataset")
  saveRDS(active_cases, "active_cases.rds")
}

group_by(active_cases, dataset, time) %>%
  summarise(mean=mean(active_cases), L95=quantile(active_cases, probs=0.025), U95 = quantile(active_cases,probs=0.975), .groups="keep") %>%
  mutate(time = time + 30) %>%
  ggplot(aes(x=time, y=mean, group=dataset)) + geom_line() + geom_ribbon(aes(ymin=L95,ymax=U95),alpha=0.1) +
  geom_line(data=data.frame(time=truth[,"time"]+30, active=truth[,"I"]+truth[,"Id"]), inherit.aes=F, aes(x=time, y=active), colour="red") +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey90") +
  labs(x="Time", y="Active cases")
ggsave("active_credbands.pdf", width=6, height=3, units="in", dpi=900)

# Model comparison -------------------------------------------------------------

fits_asmc.simple <- readRDS("../simple/fit_si_asmc.rds")[1:30]
fits_asmc.seir <- readRDS("../seir/fit_se_asmc.rds")[1:30]

BF <- data.frame(PD = map_dbl(fits_asmc, "logZ"),
                 SI = map_dbl(fits_asmc.simple, "logZ"),
                 SE = map_dbl(fits_asmc.seir, "logZ")) %>%
  mutate(PDvSI = exp(PD - SI),
         PDvSE = exp(PD - SE),
         SIvSE = exp(SI - SE), .keep="none")

summarise_BF <- function(bf){
  data.frame(Against = mean(bf < 1),
             Weak = mean(between(bf, 1, 3)),
             Positive = mean(between(bf, 3, 20)),
             Strong = mean(between(bf, 20, 150)),
             V.Strong = mean(bf > 150))
}

apply(BF, 2, summarise_BF) 

## compare to MCMC results
ml_mcmc.full <- readRDS("ml_mcmc.rds")
ml_mcmc.simple <- readRDS("../simple/ml_mcmc.rds")
ml_mcmc.seir <- readRDS("../seir/ml_mcmc.rds")

BF_mcmc <- data.frame(PD = map_dbl(ml_mcmc.full, ~ .),
                      SI = map_dbl(ml_mcmc.simple, ~ .),
                      SE = map_dbl(ml_mcmc.seir, ~ .)) %>%
  mutate(PDvSI = exp(PD - SI),
         PDvSE = exp(PD - SE),
         SIvSE = exp(SI - SE), .keep="none") 

apply(BF_mcmc, 2, summarise_BF)


# MCMC traceplots --------------------------------------------------------------

filter(theta_mcmc, dataset=="1") %>%
  rowid_to_column(var="Sample") %>%
  pivot_longer(cols=-c(dataset,Sample), names_to="Parameter", values_to="value") %>%
  ggplot(aes(x=Sample, y=value)) + geom_line() + facet_wrap(~Parameter, scales="free_y")
