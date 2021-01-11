# Annealed SMC for COVID-19
# November 17, 2020
# Renny Doig
#
# Analysis of BC case counts using contact tracing model with delay likelihood

# Preliminaries ----------------------------------------------------------------
rm(list=ls())
library(annealedSMC)
library(purrr)
library(lubridate)

source('likelihood_delay.R')


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

Contactmodel <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    beta <- exp(logitbeta) / (1 + exp(logitbeta))
    q <- exp(logitq)/(1 + exp(logitq))
    rho <- exp(logitrho) / (1 + exp(logitrho))
    deltaI <- exp(logdeltaI)
    deltaq <- exp(logdeltaq)
    gammaI <- exp(loggammaI)
    theta <- exp(logittheta) / (1 + exp(logittheta))
    c <- ifelse(t <= tau1, exp(logc1), ifelse(t<=tau2, exp(logc2), exp(logc3)))

    dSdt <- -(beta*c + c*q*(1-beta)) * S * (I + theta*A)/N + lambda*Sq
    dEdt <- beta*c*(1-q) * S * (I + theta*A)/N - sigma * E
    dIdt <- sigma * rho * E - (deltaI + alpha + gammaI) * I
    dAdt <- sigma * (1-rho) * E - gammaA * A
    dSqdt <- (1-beta) * c * q * S * (I + theta * A)/N - lambda * Sq
    dEqdt <- beta * c * q * S * (I + theta*A)/N - deltaq * Eq
    dHdt <- deltaI * I + deltaq * Eq - (alpha + gammaH) * H
    dRdt <- gammaI * I + gammaA * A + gammaH * H
    dNdt <- -alpha * (I + H)
    
    return(list(c(dNdt, dSdt, dEdt, dIdt, dAdt, dSqdt, dEqdt, dHdt, dRdt)))
  })
}


N <- 5.1e6
pars <- list(logitbeta = log(.145/.855),
             logitq = log(0.1/.9),
             sigma = 1/5,
             lambda = 1/14,
             logitrho = log(.6/.4),
             logdeltaI = log(.1),
             logdeltaq = log(.1),
             loggammaI = log(.1),
             gammaA = 0.139,
             gammaH = 0.2,
             alpha = 0.008,
             logittheta = log(.05/.95),
             tau1 = ymd("2020-03-18") - start_date,
             tau2 = ymd("2020-05-17") - start_date,
             logc1 = log(11),
             logc2 = log(5),
             logc3 = log(8))

# initial states
i0 <- 8

initial_states <- c(N = N,
                    S = N - i0,
                    E = 0.4 * i0,
                    I = 0.5 * i0,
                    A = 0.1 * i0,
                    Sq = 0,
                    Eq = 0,
                    H = 0,
                    R = 0)

data$n_DE <- length(initial_states)
data$DE_val_names <- names(initial_states)
data$initial_state <- initial_states


# Prior distributions ----------------------------------------------------------

prior_par <- list()
logPriorsRatio <- list()


## DE parameters
is_unknownPar <- rep(F, length(pars))
is_unknownPar[c(1,2,5,6,7,8,12,15,16,17)] <- T

mean_theta <- c(0.145, 0.1, 0.6, log(0.1), log(0.1), log(0.2), 0.05, log(11), log(5), log(8))
spread_theta <- c(20, 20, 20, 0.1, 0.1, 0.1, 50, 0.1, 0.1, 0.1)

prior_par$theta <- list(mean=mean_theta, sd=spread_theta)
logPriorsRatio$theta = function(xstar, x, par){
  with(par, {
    # extract the estimated parameters & means
    thstar <- xstar
    th <- x
    is.logit <- c(1,2,3,7)
    
    thstar[is.logit] = exp(thstar[is.logit])/(1+exp(thstar[is.logit]))
    th[is.logit] = exp(th[is.logit])/(1+exp(th[is.logit]))
    
    dxstar = sum(dnorm(thstar[-is.logit], mean[-is.logit], sd[-is.logit], log=T)) +
      sum(dbeta(thstar[is.logit], mean[is.logit]*sd[is.logit], (1-mean[is.logit])*sd[is.logit], log=T))
    
    dx = sum(dnorm(th[-is.logit], mean[-is.logit], sd[-is.logit], log=T)) +
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

mean_eta <- 5
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
    is.logit <- c(1,2,3,7)
    temp <- rep(NA, length(mean))
    temp[-is.logit] <- rnorm(6, mean[-is.logit], sd[-is.logit])
    temp[is.logit] <- rbeta(4, shape1=mean[is.logit]*sd[is.logit], shape2=(1-mean[is.logit])*sd[is.logit])
    temp[is.logit] <- log(temp[is.logit]/(1-temp[is.logit]))
    
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

hyperpar <- list(sigma_theta1 = c(0.015, 0.25, 0, 0, 0.05, 0.1, 1.5, 0.055, 0, 0, 0, 0.45, 0, 0, 0.05, 0.03, 0.25),
                 sigam_theta2 = NA,
                 ksi = 1,
                 sigma_likeliParam = 0.45)

smc_tuning_param <- list(eps=0.5, phi=0.95, n_move=1)

K <- 1000
ncore <- 8
set.seed(501)

fit_asmc <- ASMC2_LP(K,
                     smc_tuning_param,
                     data,
                     is_unknownPar,
                     Contactmodel,
                     pars,
                     likelihood_delay,
                     unlist(Lpars[!is_unknownLPar]),
                     unlist(Lpars[is_unknownLPar]),
                     reference_prior_list,
                     hyperpar,
                     ncore)
saveRDS(fit_asmc, "fit_asmc.rds")
