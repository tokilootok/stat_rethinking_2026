data {
  int<lower=1> N;
  int<lower=1> NW;
  int<lower=1> NJ;
  vector[N] score_z;
  array[N] int<lower=1, upper=NW> wine_id;
  array[N] int<lower=1, upper=NJ> judge_id;
  vector[NW] wine_origin_c;
  vector[NJ] judge_country_c;
}

parameters {
  vector[NW] quality;

  real alpha_harshness;
  real delta_harshness;
  real<lower=0> sigma_harshness;
  vector[NJ] z_harshness;

  real alpha_log_discrimination;
  real delta_log_discrimination;
  real<lower=0> sigma_log_discrimination;
  vector[NJ] z_log_discrimination;

  real delta_origin;
  real<lower=0> sigma_score;
}

transformed parameters {
  vector[NJ] harshness = alpha_harshness
                          + delta_harshness * judge_country_c
                          + sigma_harshness * z_harshness;
  vector[NJ] log_discrimination = alpha_log_discrimination
                                  + delta_log_discrimination * judge_country_c
                                  + sigma_log_discrimination * z_log_discrimination;
  vector<lower=0>[NJ] discrimination = exp(log_discrimination);
  vector[N] mu;

  for (n in 1:N) {
    mu[n] = discrimination[judge_id[n]]
            * (quality[wine_id[n]] - harshness[judge_id[n]]);
  }
}

model {
  // Fix the latent quality scale to one. Estimating this scale together with
  // every discrimination parameter would create an avoidable scale ridge.
  quality ~ normal(delta_origin * wine_origin_c, 1);

  alpha_harshness ~ normal(0, 0.5);
  delta_harshness ~ normal(0, 0.5);
  sigma_harshness ~ exponential(2);
  z_harshness ~ std_normal();

  alpha_log_discrimination ~ normal(0, 0.35);
  delta_log_discrimination ~ normal(0, 0.35);
  sigma_log_discrimination ~ exponential(2);
  z_log_discrimination ~ std_normal();

  delta_origin ~ normal(0, 0.5);
  sigma_score ~ exponential(1);

  score_z ~ normal(mu, sigma_score);
}

generated quantities {
  real discrimination_ratio_us_fr = exp(delta_log_discrimination);
  vector[N] log_lik;
  vector[N] score_z_rep;

  for (n in 1:N) {
    log_lik[n] = normal_lpdf(score_z[n] | mu[n], sigma_score);
    score_z_rep[n] = normal_rng(mu[n], sigma_score);
  }
}

