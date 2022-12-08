functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real sum_alpha = sum(alpha);
    return lgamma(sum_alpha) - lgamma(sum(y) + sum_alpha)
           // + lgamma(sum(y)+1) - sum(lgamma(to_vector(y)+1) // constant, may omit
           + sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
  }
  real atbswitch(real t,
                 real t1,
                 real t2,
                 real delta,
                 real iota1,
                 real iota2,
                 real xi
  ) {
    real p_t = (iota1/(1+exp(-delta*(t-t1+xi))) + (1-iota1)/(1+exp(-delta*(t-t1))) + iota2/(1+exp(delta*(t-t2-xi))) + (1-iota2)/(1+exp(delta*(t-t2)))) - 1;
    return(p_t);
  }
  real[] ode(real t, 
             real[] y, // [K+1] S, I_k
             real[] theta, // [5] beta, mu, nu, epsilon, tau
             real[] x_r, // [7] t1, t2, delta, iota1, iota2, xi, zeta
             int[] x_i // [2] K, t0
  ) {
    int K = x_i[1];
    real dydt[K+1];
    real beta_t = theta[1] * (1+x_r[7]*(t-x_i[2]));
    real p = atbswitch(t,x_r[1],x_r[2],x_r[3],x_r[4],x_r[5],x_r[6]);
    real epsilon[K];
    for(k in 1:K) epsilon[k] =  1.0-(k-1.0)/(K-1.0)*(1.0-theta[4]);
  
    dydt[1] = - beta_t*y[1]*sum(y[2:(K+1)]) + 
        p*theta[5]*(1-theta[2])*sum(to_vector(epsilon[1:(K-1)]) .* to_vector(y[2:K])) + 
        p*theta[5]*epsilon[K]*y[K+1] +
        (1-p)*theta[5]*sum(y[2:(K+1)]) +
        theta[3]*sum(y[2:(K+1)]);
    dydt[2] = beta_t*y[1]*y[2] - 
    p*theta[5]*theta[2]*y[2] - 
    p*theta[5]*(1-theta[2])*epsilon[1]*y[2] - 
    (1-p)*theta[5]*y[2] - 
    theta[3]*y[2];
    for(k in 3:(K)) {
      dydt[k] = beta_t*y[1]*y[k] + 
        p*theta[5]*theta[2]*y[k-1] - 
        p*theta[5]*theta[2]*y[k] - 
        p*theta[5]*(1-theta[2])*epsilon[k-1]*y[k] - 
        (1-p)*theta[5]*y[k] - 
        theta[3]*y[k];
    }
    dydt[K+1] = beta_t*y[1]*y[K+1] + 
      p*theta[5]*theta[2]*y[K] - 
      p*theta[5]*epsilon[K]*y[K+1] - 
      (1-p)*theta[5]*y[K+1] - 
      theta[3]*y[K+1];
    return(dydt);
  }
}
data {
  // simulation
  int t0; // start_year
  int S; // duration of simulation
  real ts[S]; // time bins
  // structure
  int C; // number of compartments
  int L; // length of data
  int G; // number of groups
  int K; // number of MIC classes
  int K_init; // number of MIC classes circulating at t0
  // data
  int mic_k[G,L,K]; // number of individuals with resistance
  int mic_n[G,L]; // sample size
  int res_k[G,L]; // number of individuals with resistance
  int res_n[G,L]; // sample size
  real treatment_use[G,7]; // parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu) and time-dependent increase in beta (zeta)
  // priors
  real p_nu[2];
  real p_tau;
  real p_mu[2];
  real p_epsilon[2];
  real p_prevalence[G,2]; 
  vector[K_init] p_init_mic;
  // prediction
  int pred_sample_size[G];
  int max_threshold; // maximal threshold
  int inference;
}
transformed data {
  int xis[2] = {K,t0}; 
  real xrs[G,7] = treatment_use;
}
parameters {
  real<lower=0> nu; // recovery rate  (unique)
  ordered[G] tau_raw; // treatment rate (by group)
  real<lower=0,upper=1> mu; // mutation probability (unique)
  real<lower=0,upper=1> prevalence[G]; // prevalence (by group)
  real<lower=0,upper=1> epsilon; // treatment efficacy (unique)
  simplex[K_init] init_mic[G]; // initial resistance level (by group)
  real<lower=0> phi; // overdispersion
}
transformed parameters {
  real<lower=0> tau[G];
  real<lower=0> beta[G]; // transmission rate (by group, computed from the other parameters to get a stable prevalence)
  real init[G,C]; // initial values
  real theta[G,5]; // vector of parameters
  real y[G,S,C]; // ODE output
  simplex[K] output_mic[G,S];
  
  for(g in 1:G) {
    tau[g] = exp(tau_raw[g]);
    beta[g] = (tau[g]+nu)/(1-prevalence[g]);
    init[g,1] =  1-prevalence[g];
    for(k in 2:(1+K_init)) init[g,k] = prevalence[g] * init_mic[g,k-1];
    for(k in (2+K_init):C) init[g,k] = 0.0;
    theta[g] =  {beta[g], mu, nu, epsilon, tau[g]};
    y[g] = integrate_ode_rk45(
      ode, 
      init[g], // initial states
      t0, // t0
      ts, // evaluation dates (ts)
      theta[g], // parameters
      xrs[g], // real data
      xis, // integer data
      1.0E-6, 1.0E-6, 1.0E3);
    for(i in 1:S) for(j in 1:C) y[g,i,j] =  y[g,i,j] + 2*1.0E-6; // avoid any negative values du to approximation of small values by adding some value greater than the ODE solver's precision
    for(i in 1:S) output_mic[g,i] = ( to_vector(y[g,i,2:C]) ) ./ sum( to_vector(y[g,i,2:C]) );
  }
}
model {
  // priors
  epsilon ~ beta(p_epsilon[1],p_epsilon[2]);
  mu ~ beta(p_mu[1],p_mu[2]);
  nu ~ gamma(p_nu[1],p_nu[2]);
  phi ~ exponential(0.01);
  for(g in 1:G) {
    tau_raw[g] ~ normal(0,p_tau);
    prevalence[g] ~ normal(p_prevalence[g,1],p_prevalence[g,2]);
    init_mic[g] ~ dirichlet(p_init_mic);
  }
  // likelihood
  if(inference==1) {
    for(g in 1:G) for(i in 1:L) mic_k[g,i] ~ dirichlet_multinomial(phi*output_mic[g,i] );
  }
}
generated quantities {
  simplex[K] pred_trans_output_mic[G,S];
  int pred_mic_n[G,S,K];
  int pred_res_n[G,S];
  vector[K] pred_mic_prop[G,S];
  vector[S] pred_resistance[G];
  vector[max_threshold] time_to_threshold[G];
  real p[G,S];
  
  for(g in 1:G) {
    for(i in 1:S) {      
      pred_trans_output_mic[g,i] = dirichlet_rng(phi * output_mic[g,i]);
      pred_mic_n[g,i] = multinomial_rng( pred_trans_output_mic[g,i] , i<=L ? mic_n[g,i] : pred_sample_size[g] );
      pred_res_n[g,i] = pred_mic_n[g,i,K];
      pred_mic_prop[g,i] = (to_vector(pred_mic_n[g,i]) + 0.0) ./ sum(pred_mic_n[g,i]);
      pred_resistance[g,i] = (pred_res_n[g,i] + 0.0) / sum(pred_mic_n[g,i]);
    }
    for(i in 1:max_threshold) {
      for(j in 1:S) {
        time_to_threshold[g,i] = 0.0;
        time_to_threshold[g,i] = (pred_resistance[g,j]<=i/100.0) ? j : time_to_threshold[g,i] ;
      }
    }
    for(i in 1:S) p[g,i] = atbswitch(ts[i],treatment_use[g,1],treatment_use[g,2],treatment_use[g,3],treatment_use[g,4],treatment_use[g,5],treatment_use[g,6]);
  }
}
