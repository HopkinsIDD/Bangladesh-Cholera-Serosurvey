//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
    int<lower=1> N_cohort_pos;
    int<lower=1> N_cohort_pos_ids;
    array[N_cohort_pos] int<lower=1, upper=N_cohort_pos_ids> id_cohort_pos;
    array[N_cohort_pos] int<lower=0,upper=1> z_cohort_pos;
    int<lower=1> T_max;
    int<lower=1> K;
    matrix[N_cohort_pos,K] t_ort;
    matrix[T_max,K] t_new_ort;
    int<lower=1> N_cohort_neg;
    int<lower=1> N_cohort_neg_ids;
    array[N_cohort_neg] int<lower=1, upper=N_cohort_neg_ids> id_cohort_neg;
    array[N_cohort_neg] int<lower=0,upper=1> z_cohort_neg;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
    real alpha_sens;
    vector[K] beta;
    real alpha_spec;
    simplex[T_max] d;
    real<lower=0> sigma_sens_id;
    vector[N_cohort_pos_ids] eta_sens_id;
    real<lower=0> sigma_spec_id;
    vector[N_cohort_neg_ids] eta_spec_id;
}

// 'sens_tv' is the time-varying sensitivity at each 'd' from 1:T_max. 'sens' is
// the overall sensitivity integrated across the estimated days back 'd'.
transformed parameters{
    // vector[N_cohort_pos_ids] alpha_sens_j;
    vector[N_cohort_pos] sens_id;
    // vector[N_cohort_neg_ids] alpha_spec_j;
    vector<lower=0, upper=1>[N_cohort_neg] spec_id;

    // alpha_sens_j=alpha_sens + sigma_sens_id * eta_sens_id;
    // sens_id=alpha_sens_j[id_cohort_pos] + t_ort * beta;
    sens_id=alpha_sens + sigma_sens_id * eta_sens_id[id_cohort_pos] + t_ort * beta;
    // alpha_spec_j=alpha_spec + sigma_spec_id * eta_spec_id;
    // spec_id=inv_logit(alpha_spec_j[id_cohort_neg]);
    spec_id=inv_logit(alpha_spec + sigma_spec_id * eta_spec_id[id_cohort_neg]);
}

model {
    target += bernoulli_logit_lpmf(z_cohort_pos | sens_id);
    target += bernoulli_lpmf(z_cohort_neg | 1-spec_id);
    target += normal_lpdf(beta | 0,10);
    target += normal_lpdf(alpha_sens | 0,1);
    target += normal_lpdf(alpha_spec | 0,1);
    target += normal_lpdf(eta_sens_id | 0,1);
    target += student_t_lpdf(sigma_sens_id | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
    target += normal_lpdf(eta_spec_id | 0,1);
    target += student_t_lpdf(sigma_spec_id | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
}

generated quantities {
    vector<lower=0, upper=1>[T_max] sens_tv;
    real<lower=0, upper=1> sens;
    real<lower=0, upper=1> spec;
    vector[N_cohort_pos] log_lik;

    for(t in 1:T_max){
        // sens_tv[t]=mean(inv_logit(alpha_sens_j + t_new_ort[t,] * beta));
        sens_tv[t]=mean(inv_logit(alpha_sens + sigma_sens_id * eta_sens_id + t_new_ort[t,] * beta));
    }

    sens=sens_tv' * d;
    // spec=mean(inv_logit(alpha_spec_j));
    spec=mean(inv_logit(alpha_spec + sigma_spec_id * eta_spec_id));


    for(i in 1:N_cohort_pos){
        log_lik[i] = bernoulli_logit_lpmf(z_cohort_pos[i] | sens_id[i]);
    }
}
