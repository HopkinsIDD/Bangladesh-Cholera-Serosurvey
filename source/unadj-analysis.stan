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


// In the unadjusted analysis, we do not adjust our predictions for sensitivity and specificity and thus only use data from the serosurvey.
// The serosurvey consists of 'N_survey' households from 'N_comm' communities, indexed by 'comm' with populations 'comm_pop'.
// Each household has 'hh_obs' members.
// 'z_survey' are the predictions of the seroincidence for each household from a random forest model.
// We extrapolate the results to the unsampled population within sampled communities, 'pop_unobs', as well as to 'N_unsamp' communities not within the sample with populations 'pop_unsamp'.
data {
    int<lower=1> N_survey;
    int<lower=1> N_comm;
    int<lower=1, upper=N_comm> comm[N_survey];
    int<lower=1> comm_pop[N_comm];
    int<lower=1> hh_obs[N_survey];
    int<lower=0> z_survey[N_survey];
    int<lower=0> pop_unobs[N_comm];
    int<lower=0> N_unsamp;
    int<lower=0> pop_unsamp[N_unsamp];
}

// 'alpha_0', 'gamma_0', 'sigma_0', and 'alpha_comm' are hierarchical parameters for estimating 'pi_h' (below).
parameters {
    real gamma_0;
    real alpha_0;
    real<lower=0> sigma_0;
    real alpha_comm[N_comm];
}

// 'pi_h' is the underlying household level seroincidence rate.
transformed parameters{
    vector<lower=0, upper=1>[N_survey] pi_h;

    for(i in 1:N_survey){
        pi_h[i]=inv_logit(alpha_comm[comm[i]]);
    }
}

//  We observe 'z_survey' seroincident cases as a binomial distribution based on the number of observations in each household and the proportion pi_h.
// The alpha_comm for each community is distributed normally with the mean equal to alpha_0 + gamma_0 * the log population of the community and standard deviation equal to sigma.
// alpha_comm determines pi_h for each household.
model {
    z_survey ~ binomial(hh_obs, pi_h);

    for(j in 1:N_comm){
        alpha_comm[j] ~ normal(alpha_0+gamma_0*log(comm_pop[j]), sigma_0);
    }

    alpha_0 ~ normal(0,1);
    gamma_0 ~ normal(0,1);
    sigma_0 ~ normal(0,1) T[0,];
}

//  We generate predictions for the sampled communities 'y_survey', the unsampled people within sampled communities 'y_unobs', and the unsampled communities 'y_unsamp'.
//  The predictions for the sampled communities are based on the parameters for each sampled community.
//  The predictions for the unsampled communities are based on the randomly drawn alpha_unsamp, which comes from the country-wide random effects and the community size.
// To generate the overall seroincidence 'pi_overall' within the code, it can be uncommented here (otherwise it is calculated in R).
// To generate the log-likelihood, uncomment the lines with 'log_lik'.
generated quantities{
    vector<lower=0>[N_survey] y_survey;
    vector<lower=0>[N_comm] y_unobs;
    vector[N_unsamp] alpha_unsamp;
    vector<lower=0>[N_unsamp] y_unsamp;
    real<lower=0, upper=1> pi_overall;
    // vector[N_survey] log_lik;

    for(i in 1:N_survey){
        y_survey[i] = binomial_rng(hh_obs[i],pi_h[i]);
    }

    for(j in 1:N_comm){
        y_unobs[j] = binomial_rng(pop_unobs[j], inv_logit(alpha_comm[j]));
    }

    for(j in 1:N_unsamp){
        alpha_unsamp[j] = normal_rng(alpha_0+gamma_0*log(pop_unsamp[j]), sigma_0);
        y_unsamp[j] = binomial_rng(pop_unsamp[j], inv_logit(alpha_unsamp[j]));
    }

    pi_overall=(sum(y_survey)+sum(y_unobs)+sum(y_unsamp))/(sum(hh_obs)+sum(pop_unobs)+sum(pop_unsamp));

    // for(i in 1:N_survey){
    //     log_lik[i] = binomial_lpmf(z_survey[i] | hh_obs[i], pi_h[i]*sens + (1-pi_h[i])*(1-spec));
    // }
}
