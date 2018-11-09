data {
 //int<lower = 0> n_juv;
 int<lower = 0> n_adu;
 int<lower = 0> n_prey_species;

 //vector<lower=0>[n_juv] C_juv;
 vector[n_adu] Contam_internal_adu;
 
 vector[n_prey_species] mean_B_sp;
 vector[n_prey_species] mean_contam_item;
}

parameters {
 // Trophic links
 simplex[n_prey_species] prop_item;

 // Contam in prey items
 real prec_contam_item;
 vector[n_prey_species] contam_prey_sp;
 real<lower = 0> k_out_adu;
 real assimil_rate;
 real prec_contam_int;
}
transformed parameters {
 row_vector[n_prey_species] Ing_prey_sp;
 real mean_Ing_contam;
 real mean_contam_int;

  // trophic ingestion
  for(i in 1:n_prey_species){
    Ing_prey_sp[i] = prop_item[i] * mean_B_sp[i];
  }

 // // trophic contamination
 mean_Ing_contam = assimil_rate * (Ing_prey_sp * contam_prey_sp); // sum: row_vector * vector !
 mean_contam_int = mean_Ing_contam / k_out_adu; //* (1-exp(- k_out_adu * time));

}
model {
 // Trophic links
 prop_item ~ dirichlet(mean_B_sp); // K-simplex distribution, PARAMETERIZATION ???

 // Contam in prey items
 prec_contam_item ~ gamma(.01,.01); // NO IDEA !!!
 mean_contam_item ~ normal(contam_prey_sp, prec_contam_item);
 assimil_rate ~ gamma(.01,.01);
 
 // // Internal contam
 prec_contam_int ~ gamma(1,1); // NO IDEA !!!
 Contam_internal_adu ~ lognormal(mean_contam_int, prec_contam_int);
 k_out_adu ~ gamma(1,1); // NO IDEA !!!
 
}
generated quantities {
   vector[n_adu] Contam_internal_adu_yrep;
   for(i in 1:n_adu){
     Contam_internal_adu_yrep[i] = lognormal_rng(mean_contam_int, prec_contam_int);
   }
}
