##--------------------------------------------------------------------------------------------------------
## SCRIPT : Simuler des donn?es de croissance
##
## Authors : Matthieu Authier
## Last update : 2018-11-08
##
## R version 3.5.1 (2018-07-02) -- "Feather Spray"
## Copyright (C) 2016 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

lapply(c("rgdal", "ggplot2", "ggthemes", "dplyr", "broom", "rgeos", "maptools", "lubridate", "rstan", 'loo', "mvtnorm", 'coda'), 
       library, character.only = TRUE)
# rstan (Version 2.14.1, packaged: 2017-04-19 05:03:57 UTC, GitRev: 2e1f913d3ca3)
# For execution on a local, multicore CPU with excess RAM we recommend calling
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

rm(list = ls())

### compilation
massage <- stan_model(file = paste("src/stan_files", "Age_masse.stan", sep = "/"),
                      model_name = "Mass-Age model (lognormal likelihood) with monotonic splines")

### data simulation
simdata <- function(seed = 20181108, n_knot = 16, n_obs = 1e3, age = NULL, mu = log(30), sigma_spline = 1, sigma_res = 0.1, slope = log(2)/3) {
  set.seed(seed)
  ### prepare splines: B-splines Eilers & Marx (2010)
  design_matrix <- function(x, xl, xr, ndx, deg = 3, d = 2) {
    # x = covariate
    # xl = lower bound
    # xr = upper bound
    # ndx = number of knots
    # deg = order of splines
    # d = order difference (1 or 2) for smoothness
    ### from Eilers & Marx
    bbase <- function(x, xl, xr, ndx, deg){
      tpower <- function(x, t, p){
        # Truncated p-th power function
        (x - t) ^ p * (x > t)
      }
      # Construct a B-spline basis of degree 'deg'
      dx <- (xr - xl) / ndx
      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      P <- outer(x, knots, tpower, deg)
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
      B <- (-1) ^ (deg + 1) * P %*% t(D)
      return(B)
    }
    B <- bbase(x, xl, xr, ndx, deg)
    # matrix of d-order differences
    D <- diff(diag(dim(B)[2]), diff = d)
    Q <- solve(D %*% t(D), D)
    X <- outer(x, 1:(d-1), "^") # fixed part without the intercept
    Z <- B %*% t(Q) # random part for mixed model formulation
    return(list(B = B, X = X, Z = Z))
  }
  ### rescaling
  rescale <- function(x) { (x - mean(x))/sd(x) }
  ### spline coef.
  b <- sigma_spline * rescale(1 - exp(-0.3 * 0:(n_knot + 2)))
  ### spline basis
  if(is.null(age)) {
    p = abs(rnorm(30))
    p = p[rev(order(p))]/sum(p)
    age = (sample.int(30, size = n_obs, replace = TRUE, prob = p) - 1) / 2
    n_obs = length(age)
  }
  B <- design_matrix(x = age, xl = 0, xr = 15, ndx = n_knot)$B
  ### linear predictor
  linpred <- mu + drop(B %*% b)
  ### std. deviation
  sigma <- sigma_res * exp(slope * rescale(age))
  ### log mean
  # logmean <- linpred - 0.5 * sigma * sigma
  ### data
  return(list(n_obs = n_obs, 
              n_knot = n_knot, 
              n_deg = 3,
              MASS = rlnorm(length(linpred), meanlog = linpred - 0.5 * sigma * sigma, sdlog = sigma),
              AGE = age,
              B = B,
              prior_scale_sigma = 1,
              prior_scale_slope = log(2)/2,
              prior_scale_mu = 5,
              prior_location_mu = 0
              )
         )
}
standata <- simdata()
with(standata, plot(AGE, MASS,las = 1, bty = 'n', xlab = "Age", ylab = "Mass (kg)"))

### sampling ...
fit <- sampling(massage, 
                data = standata, 
                pars=c('mu', 'slope', 'b', 'sigma_res', 'sigma_spline', 'log_lik', 'y_rep'), 
                chains = 4, 
                iter = 400, 
                warmup = 150,
                thin = 1,
                control = list(adapt_delta = 0.9, max_treedepth = 15)
                )
stan_rhat(fit)
loo_fit <- loo(extract_log_lik(fit))
print(loo_fit)

print(fit, digits_summary = 2)

stanpredict <- function(age = seq(0, 15, 0.5), fit, standata) {
  ### prepare splines: B-splines Eilers & Marx (2010)
  design_matrix <- function(x, xl, xr, ndx, deg = 3, d = 2) {
    # x = covariate
    # xl = lower bound
    # xr = upper bound
    # ndx = number of knots
    # deg = order of splines
    # d = order difference (1 or 2) for smoothness
    ### from Eilers & Marx
    bbase <- function(x, xl, xr, ndx, deg){
      tpower <- function(x, t, p){
        # Truncated p-th power function
        (x - t) ^ p * (x > t)
      }
      # Construct a B-spline basis of degree 'deg'
      dx <- (xr - xl) / ndx
      knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
      P <- outer(x, knots, tpower, deg)
      n <- dim(P)[2]
      D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
      B <- (-1) ^ (deg + 1) * P %*% t(D)
      return(B)
    }
    B <- bbase(x, xl, xr, ndx, deg)
    # matrix of d-order differences
    D <- diff(diag(dim(B)[2]), diff = d)
    Q <- solve(D %*% t(D), D)
    X <- outer(x, 1:(d-1), "^") # fixed part without the intercept
    Z <- B %*% t(Q) # random part for mixed model formulation
    return(list(B = B, X = X, Z = Z))
  }
  ### summary function
  lower <- function(x, alpha = 0.80) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[1]) }
  upper <- function(x, alpha = 0.80) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[2]) }
  get_summary <- function(x) { return(c(lower(x), mean(x), upper(x))) }
  
  stdage <- (age - mean(standata$AGE))/sd(standata$AGE)
  linpred <- matrix(rep(rstan::extract(fit, 'mu')$mu, each = length(age)), ncol = length(age), byrow = TRUE) +
    rstan::extract(fit, 'b')$b %*% t(design_matrix(x = age, xl = 0, xr = 15, ndx = standata$n_knot)$B)
  sigma <- matrix(rep(rstan::extract(fit, 'sigma_res')$sigma_res, each = length(age)), ncol = length(age), byrow = TRUE) *
    exp(rstan::extract(fit, 'slope')$slope %*% t(stdage))
  return(list(linpred = linpred, sigma = sigma, 
              mass = data.frame(age = age,
                                lower_bound = apply(exp(linpred - 0.5 * sigma * sigma), 2, lower),
                                mean_mass = apply(exp(linpred - 0.5 * sigma * sigma), 2, mean),
                                upper_bound = apply(exp(linpred - 0.5 * sigma * sigma), 2, upper),
                                se_mass = apply(exp(linpred - 0.5 * sigma * sigma), 2, sd)
                                )
              )
         )
}

my_pred <- stanpredict(fit = fit, standata = standata)

theme_set(theme_bw(base_size = 20))
ggplot() +
  geom_ribbon(data = my_pred$mass,
              aes(x = age, ymin = lower_bound, ymax = upper_bound), 
              fill = "lightblue", color = "lightblue"
              ) +
  geom_line(data = my_pred$mass,
            aes(x = age, y = mean_mass),
            color = "midnightblue"
            ) +
  geom_point(data = with(standata, data.frame(x = AGE, y = MASS)),
             aes(x = x, y = y), alpha = 0.1
             ) +
  ylab("Mass") + xlab("Age") +
  theme(legend.position = "top", 
        plot.title = element_text(lineheight = 0.8, face = "bold"), 
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        )
