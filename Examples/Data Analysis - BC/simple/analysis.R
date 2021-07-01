# Annealed SMC for COVID-19
# Renny Doig
# November 15, 2020
#
# Examining results from simple physical distancing model fit to confirm model
#  fit was reasonable

# Preliminaries ----------------------------------------------------------------
rm(list=ls())
library(deSolve)
library(tidyverse)
library(lubridate)
library(cowplot)
library(annealedSMC)

source('likelihood_delay.R')
source("functions.R")


# Load data --------------------------------------------------------------------

dat <- read.csv("../../../Data/BCCDC_COVID19_Dashboard_Case_Details.csv") %>%
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

# initialize population parameters
# split initial number of cases into exposed and infectious
# a 40:60 split is what was used in the "Long time frames..." paper with Caroline & I trust their judgement
# assume 40 cases; maybe estimate this parameter as well
N <- 5.1e6
# f = 1 until Mar 31 (28 days from start)
# f = 0.35 Apl 1 until May 16 (74 days from start)
# f = 0.65 May 17 onwards 
# so: tau1 = 29 and tau2 = 75
# named list of DE parameters
pars <- list(N = N,
             logR0 = log(2.57),
             delta = 1/3,
             tau1 = ymd("2020-03-18") - start_date,
             tau2 = ymd("2020-05-17") - start_date,
             logitf1 = log(0.35/0.65),
             logitf2 = log(0.65/0.35),
             gamma = 1/5)

# initial states
i0 = 8
initial_states <- c(S = N - i0,
                    E = 0.4 * i0,
                    I = 0.6 * i0)

data$n_DE <- length(initial_states)
data$DE_val_names <- names(initial_states)
data$initial_state <- initial_states


# Prior distributions ----------------------------------------------------------

prior_par <- list()
logPriorsRatio <- list()


## DE parameters
is_unknownPar <- c(F, T, F, F, F, T, T, F)

mean_theta <- c(log(2.57), 0.35, 0.65)
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
Lpars <- list(delay_scale = 9.85,
              delay_shape = 1.73,
              p1 = 0.35,
              p2 = 0.68,
              logpsi = log(5))
is_unknownLPar <- c(F, F, F, F, T)

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


# ASMC -------------------------------------------------------------------------

hyperpar <- list(sigma_theta1 = c(0, 0.01, 0, 0, 0, 0.05, 0.25, 0),
                 sigma_theta2 = NA,
                 ksi = 1,
                 sigma_likeliParam = 0.5)

smc_tuning_param <- list(eps=0.5, phi=0.95, n_move=1)

K <- 1000
ncore <- 8
set.seed(501)

if(file.exists("fit_asmc.rds")){ fit_asmc <- readRDS("fit_asmc.rds")}else{
  fit_asmc <- ASMC(K,
                   smc_tuning_param,
                   data,
                   is_unknownPar,
                   Simplemodel,
                   pars,
                   likelihood_delay,
                   unlist(Lpars[!is_unknownLPar]),
                   unlist(Lpars[is_unknownLPar]),
                   reference_prior_list,
                   hyperpar,
                   ncore)
  saveRDS(fit_asmc, "fit_asmc.rds")
}


# Trajectories -----------------------------------------------------------------

resample_seed <- 1996

# full trajectory
trajectory <- get_full_trajectory(fit_asmc, seq(0, 146, by=1), initial_states, Simplemodel, resample_seed)

ac.plot <- mutate(trajectory, active=I) %>%
  group_by(time) %>%
  summarise(mean=mean(active), L95=quantile(active, probs=0.025), U95=quantile(active, probs=0.975), .groups="keep") %>%
  mutate(Date = ymd(start_date + time)) %>%
  ggplot(aes(x=Date, y=mean)) + geom_line() +
  geom_ribbon(aes(ymin=L95, ymax=U95), alpha=0.3) + 
  geom_vline(xintercept = ymd(start_date + c(pars$tau1, pars$tau2, 69)), linetype="dotted") +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey90") +
  labs(x="", y="Active cases")

ob <- data.frame(Date = dat$Date, active=data$y) %>%
  ggplot(aes(x=Date, y=active)) + geom_point(size=0.9) +
  geom_vline(xintercept = ymd(start_date + c(pars$tau1, pars$tau2, 69)), linetype="dotted") +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey90") +
  labs(y="Case counts")
plot_grid(ac.plot, ob, align="v", nrow=2, rel_heights=c(2/3, 1/3))
ggsave("active_cases.pdf", width=6, height=4, units="in", dpi=900)


# Parameter summaries ----------------------------------------------------------

theta_long <- get_theta_long(fit_asmc, resample_seed)

# data.frame of weighted particles
theta <- map_dfr(fit_asmc$particle_list, "theta") %>% 
  bind_cols(map_dfr(fit_asmc$particle_list, "likeliParam")) %>%
  mutate(R0 = exp(logR0),
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


# Compare trajectories ---------------------------------------------------------

traj_full <- readRDS("../physical_distancing/trajectory.rds") %>%
  mutate(active = I + Id, model="Full") %>%
  select(c(time, active, model))
traj_contact <- readRDS("../contact_tracing/trajectory.rds") %>%
  mutate(active=I+H, model="Contact") %>%
  select(c(time, active, model))

ac <- mutate(trajectory, active=I, model="Simple") %>%
  select(c(time, active, model)) %>%
  bind_rows(traj_full) %>%
  bind_rows(traj_contact) %>%
  group_by(time, model) %>%
  summarise(mean=mean(active), L95=quantile(active, probs=0.025), U95=quantile(active, probs=0.975), .groups="keep") %>%
  mutate(Date = ymd(start_date + time)) %>%
  ggplot(aes(x=Date, y=mean)) + geom_line(aes(colour=model)) +
  geom_ribbon(aes(ymin=L95, ymax=U95, fill=model), alpha=0.25) + 
  geom_vline(xintercept = ymd(start_date + c(pars$tau1, pars$tau2)), linetype="dotted") +
  scale_fill_manual(name="Model", values=c("tomato", "darkgreen", "dodgerblue")) +
  scale_colour_manual(name="Model", values=c("tomato", "darkgreen", "dodgerblue")) +
  theme(panel.background=element_rect(fill="white", colour="black"),
        legend.position=c(0.75, 0.65),
        legend.background = element_rect(colour="black")) +
  background_grid(major="y", colour.major="grey90") +
  labs(y="Active cases", x="")
ac
ggsave("trajectories.pdf", width=6, height=3, units="in", dpi=900)


ob <- data.frame(Date = dat$Date, active=data$y) %>%
  ggplot(aes(x=Date, y=active)) + geom_point() +
  theme(panel.background = element_rect(fill="white", colour="black")) +
  background_grid(major="y", colour.major="grey90") +
  labs(y="Case counts")

plot_grid(ac, ob, align="v", nrow=2, rel_heights=c(2/3, 1/3))
ggsave("trajectories.pdf", width=6, height=4, units="in", dpi=900)


# Bayes factors ----------------------------------------------------------------

data.frame(PD = readRDS("../1102_bcanalysis/fit_asmc.rds")$logZ,
                 SI = fit_asmc$logZ,
                 SE = readRDS("../1117_bcanalysis_contactV2/fit_asmc_EqE.rds")$logZ) %>%
  mutate(PDvSI = exp(PD - SI),
         PDvSE = exp(PD - SE),
         SIvSE = exp(SI - SE))
