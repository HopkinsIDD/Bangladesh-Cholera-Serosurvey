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

// We have input data from a cohort (which determines the diagnostic) and the survey (which we'd like to estimate).
// 'N_cohort_pos' are the number of seroincident observations in the cohort data.
// 't' is the time since infection for each seroincident cohort observation ('T_max' is the maximum for 't').
// 'z_cohort_pos' are the predictions of the seroincident observations.
// There were 'N_cohort_neg' observations of seronegatives in the cohort study, which were predicted as 'z_cohort_neg'.
//
// The serosurvey consists of 'N_survey' households from 'N_cell' grid cells, indexed by 'cell' with populations 'cell_pop'.
// Each household has 'hh_obs' members.
// 'z_survey' are the predictions of the seroincidence for each household from a random forest model.
data {
    int<lower=1> N_survey;
    int<lower=1> N_cell;
    int<lower=1, upper=N_cell> cell[N_survey];
    int<lower=1> cell_pop[N_cell];
    int<lower=1> hh_obs[N_survey];
    int<lower=0> z_survey[N_survey];
    int<lower=0> N_cohort_pos;
    int<lower=1> T_max;
    vector<lower=1,upper=T_max>[N_cohort_pos] t;
    int<lower=0,upper=1> z_cohort_pos[N_cohort_pos];
    int<lower=0> N_cohort_neg;
    int<lower=0,upper=1> z_cohort_neg[N_cohort_neg];
}

// We also find the log-time since infection for each cohort observation 't_log', with its orthogonal squared and cubic terms ('t_log2' and 't_log3').
// We will later estimate the sensitivity for all values from t=1:T_max, for which we need 't_log_survey', 'tls_2', and 'tls_3'.
transformed data {
    vector[T_max] t_log_raw;
    real t_log_mean;
    real t_log_sd;
    vector[N_cohort_pos] t_log;
    vector[N_cohort_pos] t_log2;
    vector[N_cohort_pos] t_log3;
    vector[T_max] t_log_survey;
    vector[T_max] tls_2;
    vector[T_max] tls_3;

    for(i in 1:T_max){
        t_log_raw[i] = log(i);
    }
    t_log_mean = sum(t_log_raw)/T_max;
    t_log_sd = sd(t_log_raw);
    for(i in 1:T_max){
        t_log_survey[i]=(log(i)-t_log_mean)/t_log_sd;
        tls_2[i]=t_log_survey[i]^2;
        tls_3[i]=t_log_survey[i]^3;
    }

    for(i in 1:N_cohort_pos){
        t_log[i]=(log(t[i])-t_log_mean)/t_log_sd;
        t_log2[i]=t_log[i]^2;
        t_log3[i]=t_log[i]^3;
    }

}

// 'spec' is the specificity of the RF predictions from the cohort model.
// The 'beta' values are the coefficients of the logistic regression model for finding the time-varying sensitivity.
// 'd' is the probability of infection 1:T_max days ago.
// 'alpha_0', 'gamma_0', 'sigma_0', and 'alpha_cell' are hierarchical parameters for estimating 'pi_h' (below).
parameters {
    real<lower=0, upper=1> spec;
    real beta_0;
    real beta_1;
    real beta_2;
    real beta_3;
    real gamma_0;
    simplex[T_max] d;
    real alpha_0;
    real<lower=0> sigma_0;
    real alpha_cell[N_cell];
}

// 'pi_h' is the underlying household level seroincidence rate.
// 'sens_tv' is the time-varying sensitivity at each 'd' from 1:T_max.
// 'sens' is the overall sensitivity integrated across the estimated days back 'd'.
transformed parameters{
    vector<lower=0, upper=1>[N_survey] pi_h;
    vector<lower=0, upper=1>[T_max] sens_tv;
    real<lower=0, upper=1> sens;

    for(i in 1:N_survey){
        pi_h[i]=inv_logit(alpha_cell[cell[i]]);
    }
    sens_tv=inv_logit(beta_0+beta_1*t_log_survey+beta_2*tls_2+beta_3*tls_3);
    sens=sens_tv' * d;
}

//  We observe 'z_survey' seroincident cases as a binomial distribution based on the number of observations in each household with the proportion equal to the sum of the true positive rate pi_h*sens and false negative rate (1-pi_h)*(1-spec).
//  The predictions of the seroincident observations in the cohort data are a logistic regression with a cubic function for log-time since infection.
//  The predictions of seronegative observations in the cohort data are Bernoulli random variables with probability (1-spec).
// The alpha_cell for each cell is distributed normally with the mean equal to alpha_0 + gamma_0 * the log population of the cell and standard deviation equal to sigma_0.
// alpha_cell determines pi_h for each household.
model {
    z_survey ~ binomial(hh_obs, pi_h*sens + (1-pi_h)*(1-spec));
    z_cohort_pos ~ bernoulli_logit(beta_0+beta_1*t_log+beta_2*t_log2+beta_3*t_log3);
    z_cohort_neg ~ bernoulli(1-spec);

    for(j in 1:N_cell){
        alpha_cell[j] ~ normal(alpha_0+gamma_0*log(cell_pop[j]), sigma_0);
    }

    alpha_0 ~ normal(0,1);
    gamma_0 ~ normal(0,1);
    sigma_0 ~ normal(0,1) T[0,];
}

//  We generate predictions for the sampled grid cells 'y_cell' and proportion 'pi_cell'.
generated quantities{
    vector<lower=0>[N_cell] y_cell;
    vector<lower=0, upper=1>[N_cell] pi_cell;

    for(i in 1:N_cell){
        y_cell[i] = binomial_rng(cell_pop[i], inv_logit(alpha_cell[i]));
        pi_cell[i] = y_cell[i]/cell_pop[i];
    }
}
