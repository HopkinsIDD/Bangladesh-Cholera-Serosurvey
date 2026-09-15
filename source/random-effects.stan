//
// This Stan program defines a model for adjusting a predicted
// seroincidence by the sensitivity and specificity of the diagnostic.
// We assume that the sensitivity of the diagnostic decays with time.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// We have input data from a cohort (which determines the diagnostic) and
// the survey (which we'd like to estimate). 'z_survey' are the predictions
// of the survey serostatus per individual from a random forest model that has 'N_survey'
// observations. 'N_cohort_pos' are the number of seroincident observations
// in the cohort data. 't' is the time since infection for each seropositive
// cohort observation and 'z_cohort_pos' are the predictions of these
// observations. There were 'N_cohort_neg' observations of seronegatives in the
// cohort study, which were predicted as 'z_cohort_neg'.
data {
    int<lower=1> N_survey;
    int<lower=1> N_comm;
    int<lower=1> N_hh;
    array[N_survey] int<lower=1, upper=N_comm> comm;
    array[N_survey] int<lower=1, upper=N_hh> hh;
    array[N_comm] int<lower=1> comm_pop;
    array[N_survey] int<lower=0, upper=1> z_survey;
    real sens_logit_mean;
    real<lower=0> sens_logit_sd;
    real spec_logit_mean;
    real<lower=0> spec_logit_sd;
    int<lower=1> K;
    matrix[N_survey,K] X; // model matrix
}

// 'p' is our parameter of interest, the true seroincidence. 'spec' is the
// specificity of the RF predictions. The 'beta' values are the coefficients
// of the logistic regression model for finding the time-varying sensitivity
// (below). 'd' is the probability of infection 1:T_max days ago.
parameters {
    vector[K] beta;
    real sens_logit;
    real spec_logit;
    real gamma_0;
    real alpha_0;
    real<lower=0> sigma_comm;
    vector[N_comm] eta_comm;
    real<lower=0> sigma_hh;
    vector[N_hh] eta_hh;
}

// 'sens_tv' is the time-varying sensitivity at each 'd' from 1:T_max. 'sens' is
// the overall sensitivity integrated across the estimated days back 'd'.
transformed parameters{
    vector<lower=0, upper=1>[N_survey] p;
    vector[N_comm] alpha_comm;
    real<lower=0, upper=1> sens=inv_logit(sens_logit);
    real<lower=0, upper=1> spec=inv_logit(spec_logit);

    for(j in 1:N_comm){
        alpha_comm[j] = alpha_0 + gamma_0*log(comm_pop[j]) + sigma_comm * eta_comm[j];
    }

    p=inv_logit(alpha_comm[comm] + X * beta + sigma_hh * eta_hh[hh]);
}

//  We observe 'z_survey' cases as a result of a bernoulli process
//  with the probabilty of testing positive being the sum of the true
//  positive rate p*sens and false positive rate (1-p)*(1-spec). Sensitivity
//  is modeled as a logistic regression with a cubic polynomial for log-time.
//  We assume that the specificity is distributed as a binomial.
model {
    target+= bernoulli_lpmf(z_survey | p*sens + (1-p)*(1-spec));
    target += normal_lpdf(alpha_0 | 0,1);
    target += normal_lpdf(gamma_0 | 0,1);
    target += normal_lpdf(beta | 0,1);
    target += normal_lpdf(eta_comm | 0,1);
    target += student_t_lpdf(sigma_comm | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
    target += normal_lpdf(eta_hh | 0,1);
    target += student_t_lpdf(sigma_hh | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);

    // sensitivity and specificity modeling
    // target+= bernoulli_logit_lpmf(z_cohort_pos | sens_id);
    // target+= bernoulli_lpmf(z_cohort_neg | 1-spec_id);
    // target += normal_lpdf(alpha_sens | 0,1);
    // target += normal_lpdf(alpha_spec | 0,1);
    // target += normal_lpdf(eta_sens_id | 0,1);
    // target += student_t_lpdf(sigma_sens_id | 3, 0, 10)
    // - 1 * student_t_lccdf(0 | 3, 0, 10);
    // target += normal_lpdf(eta_spec_id | 0,1);
    // target += student_t_lpdf(sigma_spec_id | 3, 0, 10)
    // - 1 * student_t_lccdf(0 | 3, 0, 10);
    target += normal_lpdf(sens_logit | sens_logit_mean, sens_logit_sd);
    target += normal_lpdf(spec_logit | spec_logit_mean, spec_logit_sd);

}

generated quantities{
    vector<lower=0, upper=1>[N_survey] y_survey;
    vector[N_survey] log_lik;

    for(i in 1:N_survey)
        y_survey[i] = bernoulli_rng(p[i]);

    for(i in 1:N_survey){
        log_lik[i] = bernoulli_lpmf(z_survey[i] | p[i]*sens + (1-p[i])*(1-spec));
    }
}
