# Fitting
library(rstan)

# Data mining
library(dplyr) 
library(tidyr)

# Plotting
library(ggplot2)

stanmodel <- stan_model(file = paste("src/stan_files", "TK_model.stan", sep = "/"),
                      model_name = "stupid TK model")

standata <- list(n_adu = 4, # 4 adults sampled
                 Contam_internal_adu = c(0.1,0.2,0.08,0.15),
                 n_prey_species = 5,
                 mean_B_sp = 10*c(1,2,5,3,2),
                 mean_contam_item = c(0.01,0.02,0.3,0.4,0.2)
)
fit <- sampling(stanmodel, 
                data = standata, 
                #pars=c('prop_item', 'Contam_internal_adu_yrep'), 
                pars=c('prop_item','Contam_internal_adu_yrep'), 
                chains = 1, 
                iter = 400, 
                warmup = 150,
                thin = 1)

# --- ATTENTION, work only with 1 MCMC chain
df_conctam_internal_adu <- data.frame(
  ind_1 = fit@sim$samples[[1]]$`Contam_internal_adu_yrep[1]`,
  ind_2 = fit@sim$samples[[1]]$`Contam_internal_adu_yrep[2]`,
  ind_3 = fit@sim$samples[[1]]$`Contam_internal_adu_yrep[3]`,
  ind_4 = fit@sim$samples[[1]]$`Contam_internal_adu_yrep[4]`
) %>% gather(key_ind, value_ind)

df_obs <- data.frame(key_ind = c("ind_1", "ind_2", "ind_3", "ind_4"),
                     value_ind = standata$Contam_internal_adu)


ggplot() + theme_minimal() +
  scale_y_log10() +
  geom_violin(data = df_conctam_internal_adu,
              aes(x = key_ind, y = value_ind)) +
  geom_point(data = df_obs,
             aes(x = key_ind, y = value_ind))
  
