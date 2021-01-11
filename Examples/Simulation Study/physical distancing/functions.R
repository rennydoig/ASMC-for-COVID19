# Helper functions for plotting, processing, etc.

# Single simulation function ---------------------------------------------------
sim_data <- function(n, likelihood, ODEmodel, initial_states, times, pars, ...){
  # solve the DE
  fit <- ode(y=initial_states, times=times, func=ODEmodel, parms=pars, ...)
  # get the mean of the likelihood
  means <- likelihood$getMean(fit, pars)
  # sampling indices
  sampled_indices <- seq(1, nrow(means), length.out=n)
  
  psi <- exp(pars$logpsi)
  
  list( y = rnbinom(n, size=psi, mu=means$mu[sampled_indices]),
        t = means$t[sampled_indices],
        n_DE = length(initial_states),
        DE_val_names = names(initial_states),
        initial_state = initial_states,
        partial_likelihood=F )
}


# Multiple simulations ---------------------------------------------------------

multiSim <- function(M, n, likelihood, ODEmodel, initial_states, times, pars, ...){
  # solve the DE
  fit <- ode(y=initial_states, times=times, func=ODEmodel, parms=pars, ...)
  # get the mean of the likelihood
  means <- likelihood$getMean(fit, pars)
  # sampling indices
  sampled_indices <- seq(1, nrow(means), length.out=n)
  
  psi <- exp(pars$logpsi)
  
  obs <- map(1:M, ~ list( y = rnbinom(n, size=psi, mu=means$mu[sampled_indices]),
                          t = means$t[sampled_indices],
                          n_DE = length(initial_states),
                          DE_val_names = names(initial_states),
                          initial_state = initial_states,
                          partial_likelihood=F ))
  return(list(true=fit, sim=obs))
}


# Plot prior densities ---------------------------------------------------------

prior_plot <- function(par1, par2){
  names_par <- c("R0", "q", "ud", "f1", "f2", "psi")
  n <- 500000
  densities <- list(function(p1, p2) rnorm(n, p1, p2),
                    function(p1, p2) rnorm(n, p1, p2),
                    function(p1, p2) rnorm(n, p1, p2),
                    function(p1, p2) rbeta(n, p1*p2, (1-p1)*p2),
                    function(p1, p2) rbeta(n, p1*p2, (1-p1)*p2),
                    function(p1, p2) rgamma(n, p1^2/p2, p1/p2))
  trans <- list(function(x) exp(x),
                function(x) exp(x),
                function(x) exp(x),
                function(x) x,
                function(x) x,
                function(x) x)
  map_dfr(seq_along(names_par), function(.p){
    data.frame(value=densities[[.p]](par1[.p], par2[.p]),
               param = names_par[.p]) %>%
      mutate(value = trans[[.p]](value))
  }) %>%
    ggplot(aes(x=value)) + geom_histogram(bins=50) +
    facet_wrap(~ param, nrow=3, scales="free")
}


# Summary functions ------------------------------------------------------------

get_theta_long <- function(results, seed){
  set.seed(seed)
  
  indices <- sample(1:K, size=K, replace=T, prob=results$W)
  particles <- results$particle_list[indices]
  
  temp <- map_dfr(particles, "theta") %>%
    bind_cols(map_dfr(particles, "likeliParam")) %>%
    mutate(R0 = exp(logR0),
           q = exp(logq),
           ud = exp(logud),
           f1 = exp(logitf1) / (1 + exp(logitf1)),
           f2 = exp(logitf2) / (1 + exp(logitf2)),
           psi = exp(logpsi),
           .keep = "none") %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "value")
}

get_CI <- function(results, param, level=.05){
  lb <- level / 2
  ub <- 1 - level / 2
  
  theta <- select(results, all_of(c(param, "W")))
  theta_bound <- theta %>%
    arrange(.data[[param]]) %>%
    mutate(cum.W = cumsum(W)) %>%
    filter((cum.W > lb) & (cum.W < ub))
  
  data.frame(Parameter = param,
             mean = sum(theta[param] * theta$W),
             lower.bound = min(theta_bound[param]),
             upper.bound = max(theta_bound[param]))
}

get_CIV2 <- function(results, param, seed, level=.05){
  set.seed(seed)
  indices <- sample(1:K, K, T, results$W)
  particles <- results[indices,param] %>% set_names(nm="value")
  
  lb <- level / 2
  ub <- 1 - level / 2
  
  summarise(particles, mean=mean(value), lower.bound=quantile(value, probs=lb), upper.bound=quantile(value, probs=ub)) %>%
    bind_cols(Parameter=param)
}


get_active_cases <- function(results, times, initial_state, trans_model, seed){
  set.seed(seed)
  indices <- sample(1:K, K, T, results$W)
  particles <- results$particle_list[indices]
  theta <- map_dfr(particles, "theta")
  
  map_dfr(1:nrow(theta), function(.i){
    as.data.frame(ode(y=initial_state, times=times, func=trans_model, parms=theta[.i,])) %>%
      mutate(active_cases = I + Id) %>%
      select(time, active_cases)
  })
}
