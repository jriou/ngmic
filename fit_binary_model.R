fit_binary_model = function(res_data, groups, period_simulated, target_prevalence, treatment_use, inference ) {
  require(rstan)
  rstan_options(auto_write = TRUE,javascript=FALSE)
  data_selected = filter(res_data,gender_or %in% groups) %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  t0=period_simulated[1]-1
  S=length(period_simulated)
  C = 3
  L = length(unlist(data_selected[data_selected$gender_or == groups[[1]],"n"]))
  G = length(groups)
  res_k = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_k[i,] = unlist(data_selected[data_selected$gender_or == groups[[i]],"n"])
  res_n = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_n[i,] = unlist(data_selected[data_selected$gender_or == groups[[i]],"sample_size"])
  
  data_list = list(
    inference=inference,
    # simulation
    t0=t0, # start_year
    S=S, # duration of simulation
    ts=period_simulated, # time bins
    # structure
    C=C, # number of compartments
    L=L, # length of data
    G=G, # number of groups
    # data
    res_k=res_k, # number of individuals with resistance
    res_n=res_n,
    target_prevalence=target_prevalence, 
    treatment_use=treatment_use, # parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu)
    # priors
    p_nu=1,
    p_tau=1,
    p_mu=c(1,1000),
    p_epsilon=c(1,1),
    p_init_res=c(1,1),
    # prediction
    pred_sample_size=res_n[,L],
    max_threshold=20 # maximal threshold
  )
  samples = stan("models/binary_model.stan",data=data_list,chains=4,iter=1000,init=0,control = list(adapt_delta=.9))
  print(samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"))
  return(list(samples=samples,data_list=data_list,data_selected=data_selected))
}

fit_binary_model_flat = function(res_data, groups, period_simulated, target_prevalence, treatment_use, inference ) {
  require(rstan)
  rstan_options(auto_write = TRUE,javascript=FALSE)
  data_selected = filter(res_data,gender_or %in% groups) %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  t0=period_simulated[1]-1
  S=length(period_simulated)
  C = 3
  L = length(unlist(data_selected[data_selected$gender_or == groups[[1]],"n"]))
  G = length(groups)
  res_k = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_k[i,] = unlist(data_selected[data_selected$gender_or == groups[[i]],"n"])
  res_n = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_n[i,] = unlist(data_selected[data_selected$gender_or == groups[[i]],"sample_size"])
  
  data_list = list(
    inference=inference,
    # simulation
    t0=t0, # start_year
    S=S, # duration of simulation
    ts=period_simulated, # time bins
    # structure
    C=C, # number of compartments
    L=L, # length of data
    G=G, # number of groups
    # data
    res_k=res_k, # number of individuals with resistance
    res_n=res_n,
    target_prevalence=target_prevalence, 
    treatment_use=treatment_use, # parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu)
    # priors
    p_nu=1,
    p_tau=1,
    p_mu=c(1,1),
    p_epsilon=c(1,1),
    p_init_res=c(1,1),
    # prediction
    pred_sample_size=res_n[,L],
    max_threshold=20 # maximal threshold
  )
  samples = stan("models/binary_model.stan",data=data_list,chains=4,iter=1000,init=0,control = list(adapt_delta=.9))
  print(samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"))
  return(list(samples=samples,data_list=data_list,data_selected=data_selected))
}

fit_binary_model_incprev = function(res_data, groups, period_simulated, target_prevalence, treatment_use, zeta, inference ) {
  require(rstan)
  rstan_options(auto_write = TRUE,javascript=FALSE)
  data_selected = filter(res_data,gender_or %in% groups) %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  t0=period_simulated[1]-1
  S=length(period_simulated)
  C = 3
  L = length(unlist(data_selected[data_selected$gender_or == groups[[1]],"n"]))
  G = length(groups)
  res_k = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_k[i,] = unlist(data_selected[data_selected$gender_or == groups[[i]],"n"])
  res_n = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_n[i,] = unlist(data_selected[data_selected$gender_or == groups[[i]],"sample_size"])
  
  data_list = list(
    inference=inference,
    # simulation
    t0=t0, # start_year
    S=S, # duration of simulation
    ts=period_simulated, # time bins
    # structure
    C=C, # number of compartments
    L=L, # length of data
    G=G, # number of groups
    # data
    res_k=res_k, # number of individuals with resistance
    res_n=res_n,
    target_prevalence=target_prevalence, 
    treatment_use= cbind(treatment_use,zeta=zeta), # parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu) and time-dependent increase in beta (zeta)
    # priors
    p_nu=1,
    p_tau=1,
    p_mu=c(1,1),
    p_epsilon=c(1,1),
    p_init_res=c(1,1),
    # prediction
    pred_sample_size=res_n[,L],
    max_threshold=20 # maximal threshold
  )
  samples = stan("models/binary_modelB.stan",data=data_list,chains=4,iter=1000,init=0,control = list(adapt_delta=.9))
  print(samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"))
  return(list(samples=samples,data_list=data_list,data_selected=data_selected))
}