data {
 int<lower = 1> n_obs;
 int<lower = 1> n_knot;
 int<lower = 1> n_deg;
 vector<lower=0>[n_obs] MASS;
 vector[n_obs] AGE;
 matrix[n_obs, n_knot + n_deg] B; // Beziers-splines basis
 real<lower = 0.0> prior_scale_sigma; // pb dependent
 real<lower = 0.0> prior_scale_slope; // log(2)/3 for example
 real<lower = 0.0> prior_scale_mu;
 real prior_location_mu;
}

transformed data {
 vector[n_obs] STDAGE;
 STDAGE = (AGE - mean(AGE)) / sd(AGE);
}

parameters {
 real u_slope;
 real u_mu;
 ordered[n_knot + n_deg] u_b;
 real<lower = 0> tau;
 real<lower = 0> unscaled_sigma2;
 real  logit_prop;
}

transformed parameters {
 real prop;
 real slope;
 real mu;
 vector[n_obs] sigma;
 vector[n_obs] logmean;
 vector[n_obs] res;
 vector[n_knot + n_deg] b;
 real sigma_tot;
 real sigma_spline;
 real sigma_res;
 slope = u_slope * prior_scale_slope;
 mu = prior_location_mu + u_mu * prior_scale_mu;
 prop = inv_logit(1.5 * logit_prop);
 sigma_tot = prior_scale_sigma * sqrt(unscaled_sigma2/tau);
 sigma_res = sqrt(prop) * sigma_tot;
 sigma_spline = sqrt(1 - prop) * sigma_tot;
 // heteroskedasticity
 sigma = sigma_res * exp(slope * STDAGE);
 b = sigma_spline * u_b;
 // mean
 logmean = rep_vector(mu, n_obs) + B * b - 0.5 * square(sigma);
 res = (log(MASS) - logmean) ./ sigma;
}

model {
 // variance partitioning
 logit_prop ~ normal(0.0, 1.0);

 // scaled beta-2 prior for variance
 tau ~ gamma(1.0, 1.0);
 unscaled_sigma2 ~ gamma(0.5, 1.0);

 // intercept and slope parameters
 u_slope ~ normal(0.0, 1.0);
 u_mu ~ normal(0.0, 1.0);

 // splines
 u_b ~ normal(0.0, 1.0);

 // lognormal likelihood because a mass can't be negative
 MASS ~ lognormal(logmean, sigma);
}

generated quantities {
 vector[n_obs] log_lik;
 vector[n_obs] y_rep;
 for(i in 1:n_obs) {
  log_lik[i] = lognormal_lpdf(MASS[i] | logmean[i], sigma[i]);
  y_rep[i] = lognormal_rng(logmean[i], sigma[i]);
 }
}
