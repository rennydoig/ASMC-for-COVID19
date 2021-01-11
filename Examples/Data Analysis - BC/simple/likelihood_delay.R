# Likelihood functions
#
# Negative binomial with observation delay
likelihood_delay <- list()

likelihood_delay$dFun <- function(data, fit, ODEpars, likeliParam, log=T)
{
  require(dplyr)
  with(as.list(c(likeliParam, data, ODEpars)), {
    mean_delay = delay_scale*gamma(1 + 1/delay_shape)
    
    get_mu = function(fit, t){
      # need to be able to check days as early as 2*mean_delay days before t
      if(min(fit[,"time"]) > t - 2*mean_delay) stop("Need earlier start time for data.")
      
      # find the time points that fall between t and 2 x (mean delay time)
      indices = which(between(fit[,"time"], t - 2*mean_delay, t))
      
      # time spacing (assuming constant spacing)
      dt = fit[indices[2], "time"] - fit[indices[1], "time"]
      
      # cases arising from each time
      new_cases = delta*fit[indices, "E"]
      
      # pull correct sampling ratio
      # t=0 corresponds to 2020-02-01, so t=11 is 2020-03-14, the change point of the sampling proportion
      change_time = lubridate::ymd("2020-04-14") - lubridate::ymd("2020-02-05")
      samp_frac = ifelse( t<=change_time, p1, p2 )
      
      # evaluate the inside of the integral
      ft = samp_frac * new_cases * dweibull(max(fit[indices,"time"])-fit[indices,"time"],
                                            delay_shape, delay_scale)
      
      # approximate integral via trapezoid rule
      return( 0.5*dt * (ft[1] + 2*sum(ft[2:(length(ft)-1)]) + ft[length(ft)]) )
    }
    
    mu = vapply(data$t[-(1:(2*mean_delay+1))],
                function(x){get_mu(fit, x)},
                FUN.VALUE = 1)
    
    psi = exp(logpsi)
    
    temp <- dnbinom(y[-(1:(2*mean_delay+1))], size=psi, mu=mu, log=log)

    return( sum(temp[!is.nan(temp)]) )
  })
}
