# Helper functions for plotting, processing, etc.

# Plot prior densities ---------------------------------------------------------

prior_plot <- function(par1, par2){
  names_par <- c("R0", "f1", "f2", "psi")
  n <- 500000
  densities <- list(function(p1, p2) rnorm(n, p1, p2),
                    function(p1, p2) rbeta(n, p1*p2, (1-p1)*p2),
                    function(p1, p2) rbeta(n, p1*p2, (1-p1)*p2),
                    function(p1, p2) rgamma(n, p1^2/p2, p1/p2))
  trans <- list(function(x) exp(x),
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

get_full_trajectory <- function(results, times, initial_state, trans_model, seed){
  set.seed(seed)
  indices <- sample(1:K, K, T, results$W)
  particles <- results$particle_list[indices]
  theta <- map_dfr(particles, "theta")
  
  map_dfr(1:nrow(theta), function(.i){
    as.data.frame(ode(y=initial_state, times=times, func=trans_model, parms=theta[.i,]))
  })
}

get_active_casesV2 <- function(results, times, initial_state, trans_model, seed){
  set.seed(seed)
  indices <- sample(1:K, K, T, results$W)
  particles <- results$particle_list[indices]
  theta <- map_dfr(particles, "theta")
  
  map_dfr(1:nrow(theta), function(.i){
    init_states <- as.numeric(initial_state[.i,])
    names(init_states) <- names(initial_state)
    
    as.data.frame(ode(y=init_states, times=times, func=trans_model, parms=theta[.i,])) %>%
      mutate(active_cases = I + Id) %>%
      select(time, active_cases)
  })
}


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

projection_sim <- function(results, n, likelihood, Lpar, trans_model, initial_states, times, seed){
  set.seed(seed)
  indices <- sample(1:K, K, T, prob=results$W)
  particles <- results$particle_list[indices]
  theta <- map_dfr(particles, "theta") %>% bind_cols(map_dfr(particles, "likeliParam"))
  
  map_dfr(1:K, function(.k){
    init_states <- as.numeric(initial_states[.k,])
    names(init_states) <- names(initial_states)
    
    fit <- ode(y=init_states, times=times, func=trans_model, parms=theta[.k,])
    means <- likelihood$getMean(fit, c(theta[.k,], Lpar))
    sampled_indices <- seq(1, nrow(means), length.out=n)
    psi <- exp(theta$logpsi[.k])
    
    data.frame(time=means$t[sampled_indices],
               cases = rnbinom(n, size=psi, mu=means$mu[sampled_indices]))
  })
}
