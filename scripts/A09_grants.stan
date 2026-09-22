data {
  int<lower=1> K;
  int<lower=1> ND;
  array[K] int<lower=0> awards;
  array[K] int<lower=0> applications;
  array[K] int<lower=1, upper=2> gender_id;
  array[K] int<lower=1, upper=ND> discipline_id;
  simplex[ND] weight_pooled;
  simplex[ND] weight_female;
  simplex[ND] weight_male;
}

parameters {
  matrix[2, ND] logit_award;
}

model {
  to_vector(logit_award) ~ normal(0, 1);
  for (k in 1:K) {
    awards[k] ~ binomial_logit(
      applications[k],
      logit_award[gender_id[k], discipline_id[k]]
    );
  }
}

generated quantities {
  matrix[2, ND] p_award;
  vector[ND] cell_risk_difference;
  real total_female;
  real total_male;
  real total_risk_difference;
  real standardized_female;
  real standardized_male;
  real standardized_direct_risk_difference;
  array[K] int awards_rep;
  vector[K] log_lik;

  for (g in 1:2) {
    for (d in 1:ND) {
      p_award[g, d] = inv_logit(logit_award[g, d]);
    }
  }
  for (d in 1:ND) {
    cell_risk_difference[d] = p_award[1, d] - p_award[2, d];
  }

  // Total contrast preserves each gender's observed discipline distribution.
  total_female = dot_product(weight_female, to_vector(p_award[1]));
  total_male = dot_product(weight_male, to_vector(p_award[2]));
  total_risk_difference = total_female - total_male;

  // Standardized controlled contrast uses one common target distribution.
  standardized_female = dot_product(weight_pooled, to_vector(p_award[1]));
  standardized_male = dot_product(weight_pooled, to_vector(p_award[2]));
  standardized_direct_risk_difference = standardized_female - standardized_male;

  for (k in 1:K) {
    real p = p_award[gender_id[k], discipline_id[k]];
    awards_rep[k] = binomial_rng(applications[k], p);
    log_lik[k] = binomial_lpmf(awards[k] | applications[k], p);
  }
}

