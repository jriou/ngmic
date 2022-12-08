functions {
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
             real[] y, // [3] S, I, J
             real[] theta, // [5] beta, mu, nu, epsilon, tau
             real[] x_r, // [6] t1, t2, delta, iota1, iota2, xi
             int[] x_i // [0] dummy
  ) {
    real dydt[3];
    real p = atbswitch(t,x_r[1],x_r[2],x_r[3],x_r[4],x_r[5],x_r[6]);
    dydt[1] = - theta[1]*y[1]*(y[2]+y[3]) + ((1-p)*theta[5]+p*theta[5]*(1-theta[2])+theta[3])*y[2] + 
              ((1-p)*theta[5]+p*theta[5]*theta[4]+theta[3])*y[3];
    dydt[2] = theta[1]*y[1]*y[2] - ((1-p)*theta[5]+p*theta[5]*(1-theta[2])+theta[3])*y[2] - p*theta[5]*theta[2]*y[2];
    dydt[3] = theta[1]*y[1]*y[3] + p*theta[5]*theta[2]*y[2] - ((1-p)*theta[5]+p*theta[5]*theta[4]+theta[3])*y[3];
    return(dydt);
  }
}
data {
  int inference;
  // simulation
  int t0; // start_year
  int S; // duration of simulation
  real ts[S]; // time bins
  // structure
  int C; // number of compartments
  int L; // length of data
  int G; // number of groups
  // data
  int res_k[G,L]; // number of individuals with resistance
  int res_n[G,L]; // sample size
  real target_prevalence[G,2]; 
  real treatment_use[G,6]; // parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu)
  // priors
  real p_nu;
  real p_tau;
  real p_mu[2];
  real p_epsilon[2];
  real p_init_res[2];
  // prediction
  int pred_sample_size[G];
  int max_threshold; // maximal threshold
}
transformed data {
  int xis[0]; 
  real xrs[G,6] = treatment_use;
}
parameters {
  real<lower=0> nu; // recovery rate  (unique)
  ordered[G] tau_raw; // treatment rate (by group)
  real<lower=0,upper=1> mu; // mutation probability (unique)
  real<lower=0,upper=1> prevalence[G]; // prevalence (by group)
  real<lower=0,upper=1> epsilon; // treatment efficacy (unique)
  real<lower=0,upper=1> init_res[G]; // initial resistance level (by group)
  real<lower=0> phi; // overdispersion
}
transformed parameters {
  real<lower=0> tau[G];
  real<lower=0> beta[G]; // transmission rate (by group, computed from the other parameters to get a stable prevalence)
  real kappa = phi + 2;
  real init[G,3]; // initial values
  real theta[G,5]; // vector of parameters
  real y[G,S,C]; // ODE output
  vector[S] output_res[G]; 
  
  for(g in 1:G) {
    tau[g] = exp(tau_raw[g]);
    beta[g] = (tau[g]+nu)/(1-prevalence[g]);
    init[g] =  {1-prevalence[g],prevalence[g]*(1-init_res[g]),prevalence[g]*init_res[g]};
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
    for(i in 1:S) output_res[g,i] = y[g,i,3] / sum(y[g,i,2:3]);
  }
}
model {
  // priors
  epsilon ~ beta(p_epsilon[1],p_epsilon[2]);
  mu ~ beta(p_mu[1],p_mu[2]);
  nu ~ exponential(p_nu);
  phi ~ exponential(0.01);
  for(g in 1:G) {
    tau_raw[g] ~ normal(0,p_tau);
    prevalence[g] ~ normal(target_prevalence[g,1],target_prevalence[g,2]);
    init_res[g] ~ beta(p_init_res[1],p_init_res[2]);
  }
  // likelihood
  if(inference==1) for(g in 1:G) for(i in 1:L) target += beta_binomial_lpmf( res_k[g,i] | res_n[g,i] , output_res[g,i] * kappa, (1-output_res[g,i]) * kappa );
}
generated quantities {
  int pred_res_k[G,S];
  vector[S] pred_res_p[G];
  vector[max_threshold] time_to_threshold[G];
  real p[G,S];
  
  for(g in 1:G) {
    for(i in 1:S) {
      pred_res_k[g,i] = beta_binomial_rng( i<=L ? res_n[g,i] : pred_sample_size[g] , output_res[g,i] * kappa, (1-output_res[g,i]) * kappa );
      pred_res_p[g,i] = pred_res_k[g,i] / ((i<=L ? res_n[g,i] : pred_sample_size[g]) + 0.0); // add 0.0 to avoid numerical problems with division by integer
    }
    for(i in 1:max_threshold) {
      for(j in 1:S) {
        time_to_threshold[g,i] = 0.0;
        time_to_threshold[g,i] = (pred_res_p[g,j]<=i/100.0) ? j : time_to_threshold[g,i] ;
      }
    }
    for(i in 1:S) p[g,i] = atbswitch(ts[i],treatment_use[g,1],treatment_use[g,2],treatment_use[g,3],treatment_use[g,4],treatment_use[g,5],treatment_use[g,6]);
  }
}
