# Annealed SMC for COVID-19
# November 2, 2020
# Renny Doig
#
# Analysis of BC case counts using only phystical distancing model

# Preliminaries ----------------------------------------------------------------
rm(list=ls())
library(annealedSMC)
library(tidyverse)
library(lubridate)
library(cowplot)
library(gridExtra)

source("functions.R")
source('likelihood_NB.R')
source("likelihood_delay.R")


# Load data --------------------------------------------------------------------

dat <- read.csv("../BCCDC_COVID19_Dashboard_Case_Details.csv") %>%
  mutate(I = 1, Date=ymd(Reported_Date)) %>%
  filter(Date>"2020-02-02") %>%
  group_by(Date) %>%
  summarise(I=sum(I)) %>%
  mutate(t = as.numeric(Date-min(Date)))

# reformat data to be compatible with our methods
data <- list()
data$y <- dat$I
data$t <- dat$t
data$is_obs <- T
data$partial_likelihood <- F

start_date <- ymd(dat$Date[1])


# Initialize model & parameters ------------------------------------------------

PhysDistmodel = function(t, state, parms){
  with(as.list(c(state, parms)), {
    R0 = exp(logR0)
    q <- exp(logq)
    ud <- exp(logud)
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

# initialize population parameters
# split initial number of cases into exposed and infectious
# a 40:60 split is what was used in the "Long time frames..." paper with Caroline & I trust their judgement
# assume 40 cases; maybe estimate this parameter as well
N <- 5.1e6
# f = 1 until Mar 31 (28 days from start)
# f = 0.35 Mar 18 until May 16 (74 days from start)
# f = 0.65 May 17 onwards 
# so: tau1 = 29 and tau2 = 75
# named list of DE parameters
pars <- list(N = N,
             logR0=log(2.57),
             k1 = 1/5,
             k2 = 1,
             logq = log(0.05),
             D = 5,
             logud = log(0.1),
             ur = 0.02,
             tau1 = ymd("2020-03-18") - start_date,
             tau2 = ymd("2020-05-17") - start_date,
             logitf1 = 0.35,
             logitf2 = 0.65)

# initial states
e <- with(pars, {ur/(ur+exp(logud))})

i0 = 8
initial_states <- c(S = (1-e) * (N-i0),
                     E1 = (1-e) * 0.4 * i0,
                     E2 = (1-e) * 0.1 * i0,
                     I = (1-e) * 0.5 * i0,
                     Q = 0,
                     R = 0,
                     Sd = e * (N-i0),
                     E1d = e * 0.4 * i0,
                     E2d = e * 0.1 * i0,
                     Id = e * 0.5 * i0,
                     Qd = 0,
                     Rd = 0)

data$n_DE <- length(initial_states)
data$DE_val_names <- names(initial_states)
data$initial_state <- initial_states


# Prior distributions ----------------------------------------------------------

prior_par <- list()
logPriorsRatio <- list()


## DE parameters
is_unknownPar <- rep(F, length(pars))
is_unknownPar[c(2,5,7,11,12)] <- T

mean_theta <- c(log(2.57), log(0.05), log(0.1), 0.35, 0.65)
spread_theta <- c(0.05, 0.15, 0.1, 20, 20)

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
Lpars <- list(delay_scale = 9.85,
              delay_shape = 1.73,
              p1 = 0.35,
              p2 = 0.68,
              logpsi = log(5))
is_unknownLPar <- c(F, F, F, F, T)

mean_eta <- 5
spread_eta <- 1
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


prior_plot(c(mean_theta, mean_eta), c(spread_theta, spread_eta))

# ASMC -------------------------------------------------------------------------

hyperpar <- list(sigma_theta1 = c(0, 0.01, 0, 0, 1.5, 0, 0.1, 0, 0, 0, 0.05, 0.1),
                 sigam_theta2 = NA,
                 ksi = 1,
                 sigma_likeliParam = 0.65)

# hyperpar_mcmc <- list(sigma_theta = c(0, 0.02, 0, 0), sigma_likeliParam=0.5)

smc_tuning_param <- list(eps=0.5, phi=0.95, n_move=1)

K <- 1000
ncore <- 8
set.seed(501)

if(file.exists("fit_asmc.rds")){ fit_asmc <- readRDS("fit_asmc.rds")}else{
  fit_asmc <- ASMC2_LP(K,
                       smc_tuning_param,
                       data,
                       is_unknownPar,
                       PhysDistmodel,
                       pars,
                       likelihood_delay,
                       unlist(Lpars[!is_unknownLPar]),
                       unlist(Lpars[is_unknownLPar]),
                       reference_prior_list,
                       hyperpar,
                       ncore)
  saveRDS(fit_asmc, "fit_asmc.rds")
}


# Results ----------------------------------------------------------------------

resample_seed <- 502

theta_long <- get_theta_long(fit_asmc, resample_seed)

# histogram of results
ggplot(theta_long, aes(x=value)) +
  geom_histogram(bins=50) +
  facet_wrap( ~ Parameter, scales="free") +
  labs(x="", y="") +
  theme(panel.background = element_rect(fill="white", colour="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  background_grid(major="y", colour.major="grey90")


# data.frame of weighted particles
theta <- map_dfr(fit_asmc$particle_list, "theta") %>% 
  bind_cols(map_dfr(fit_asmc$particle_list, "likeliParam")) %>%
  mutate(R0 = exp(logR0),
         q = exp(logq),
         ud = exp(logud),
         f1 = exp(logitf1) / (1 + exp(logitf1)),
         f2 = exp(logitf2) / (1 + exp(logitf2)),
         psi = exp(logpsi),
         .keep = "none") %>%
  bind_cols(W = fit_asmc$W)

# 95% credible intervals and posterior means from IS estimates
CIs <- map_dfr(unique(theta_long$Parameter), ~ get_CI(theta, .)) %>%
  mutate(Mean = round(mean, 4),
         "Credible Interval" = paste0("(",round(lower.bound,4),",",round(upper.bound,4),")"),
         .keep="unused")
CIs
# ggsave("CI_table.pdf", tableGrob(CIs), width=6, height=5, units="in", dpi=900)

active_cases <-  get_active_cases(fit_asmc, seq(0, 146, by=1), initial_states, PhysDistmodel, resample_seed)

obs_data <- data.frame(time = data$t, mean = data$y) %>%
  mutate(Date = ymd(start_date + time)) %>%
  ggplot(aes(x=Date, y=mean)) + geom_point(size=1) +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey85") +
  labs(y="Observed counts")

ac <- group_by(active_cases, time) %>%
  summarise(mean = mean(active_cases),
            L95 = quantile(active_cases, probs=0.025),
            U95 = quantile(active_cases, probs=0.975)) %>%
  mutate(Date = ymd(time + start_date)) %>%
  ggplot(aes(x=Date, y=mean)) + geom_line() +
  geom_ribbon(aes(ymin=L95, ymax=U95), alpha=0.3) +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey85") +
  labs(x="", y="Active cases")
plot_grid(ac, obs_data, align="v", nrow=2, rel_heights=c(2/3,1/3))
ggsave("active_cases.pdf", width=6, height=5, units="in", dpi=900)  


# Projections ------------------------------------------------------------------

full_traj <- get_full_trajectory(fit_asmc, seq(0, 146, by=1), initial_states, PhysDistmodel, resample_seed)

proj_times <- seq(146, 200, by=1)
proj_states <- map_dfr(1:K, ~ filter(full_traj, time==max(time))[.,-1])

proj_cases <- get_active_casesV2(fit_asmc, proj_times, proj_states, PhysDistmodel, resample_seed)

projection <- rbind(cbind(active_cases, Type="Obs"), cbind(proj_cases, Type="Pred"))
group_by(projection, time, Type) %>%
  summarise(mean = mean(active_cases),
            L95 = quantile(active_cases, probs=0.025),
            U95 = quantile(active_cases, probs=0.975), .groups="keep") %>%
  mutate(Date = ymd(time + start_date)) %>%
  ggplot(aes(x=Date, y=mean, colour=Type, fill=Type)) + geom_line() +
  geom_ribbon(aes(ymin=L95, ymax=U95), alpha=0.3) +
  scale_fill_manual(values = c("grey60", "steelblue")) +
  scale_colour_manual(values = c("grey60", "steelblue")) +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey85") +
  labs(y="Active cases")
ggsave("predictions.pdf", width=6, height=4, units="in", dpi=900)


## Project case counts 

proj_times_obs <- seq(129, 200, by=1)
proj_states_obs <- map_dfr(1:K, ~ filter(full_traj, time==129)[.,-1])
proj_counts <- projection_sim(fit_asmc, 54, likelihood_delay, Lpars[!is_unknownLPar], PhysDistmodel, proj_states_obs, proj_times_obs, resample_seed)

# read in updated data list (as of November 14)
dat_update <- read.csv("../data_1114.csv") %>%
  mutate(I = 1, Date=ymd(Reported_Date)) %>%
  filter(Date>"2020-02-02") %>%
  group_by(Date) %>%
  summarise(I=sum(I)) %>%
  mutate(t = as.numeric(Date-min(Date))) %>%
  filter(t <= 200)


rbind(data.frame(time=data$t, cases=data$y, Type="Obs"), cbind(proj_counts, Type="Pred")) %>%
  mutate(Date = ymd(time + start_date)) %>%
  ggplot(aes(x=Date, y=cases)) + geom_point()

group_by(proj_counts, time) %>%
  summarise(mean = mean(cases),
            L95 = quantile(cases, probs=0.025), U95 = quantile(cases, probs=0.975),
            L50 = quantile(cases, probs=0.25), U50 = quantile(cases, probs=0.75)) %>%
  mutate(Date = ymd(time + start_date)) %>%
  ggplot(aes(x=Date, y=mean)) + geom_line() +
  geom_ribbon(aes(ymin=L95, ymax=U95), alpha=0.3) +
  geom_ribbon(aes(ymin=L50, ymax=U50), alpha=0.3) +
  geom_point(data=data.frame(Date=dat$Date, mean=dat$I))

group_by(proj_counts, time) %>%
  summarise(mean = mean(cases),
            L95 = quantile(cases, probs=0.025), U95 = quantile(cases, probs=0.975),
            L50 = quantile(cases, probs=0.25), U50 = quantile(cases, probs=0.75)) %>%
  mutate(Date = ymd(time + start_date)) %>%
  ggplot(aes(x=Date, y=mean)) + geom_point(size=1, colour="steelblue") +
  geom_linerange(aes(ymin=L95, ymax=U95), size=0.3, colour="steelblue") +
  geom_linerange(aes(ymin=L50, ymax=U50), size=0.5, colour="steelblue") +
  geom_point(data=data.frame(Date=dat_update$Date, mean=dat_update$I), size=0.75) +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey85")
ggsave("proj_counts.pdf", width=6, height=4, units="in", dpi=900)

saveRDS(proj_counts, "proj_counts_physdist.rds")


# Model comparison -------------------------------------------------------------

# load in other fits
fit_contact <- readRDS("../1117_bcanalysis_contactV2/fit_asmc.rds")
fit_simple <- readRDS("../1116_bcanalysis_simpleV2/fit_asmc.rds")

(logZ.1 <- fit_simple$logZ)
(logZ.2 <- fit_asmc$logZ)
(logZ.3 <- fit_contact$logZ)

exp(logZ.2 - logZ.1) # physdist v. simple
exp(logZ.2 - logZ.3) # physdist v. contact
exp(logZ.1 - logZ.3) # simple v. contact


