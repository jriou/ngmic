fit_multistep_model = function(mic_data, groups, period_simulated, target_prevalence, treatment_use, inference=1 ) {
  rstan_options(auto_write = TRUE,javascript=FALSE)
  data_selected = filter(mic_data,gender_or %in% groups) %>%
    arrange(gender_or,year,mic) %>%
    group_by(gender_or,year,mic,sample_size,atb,res) %>%
    summarise(n=sum(n)) 
  datares = filter(mic_data,gender_or %in% groups) %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  t0=period_simulated[1]-1
  S=length(period_simulated)
  mics = unique(data_selected$mic)
  K = length(mics)
  C = K+1
  years = unique(data_selected$year)
  L = length(years)
  G = length(groups)
  mic_classes = unique(data_selected$mic)
  res_classes = unique(filter(data_selected,res==1)$mic)
  which_res = which(mic_classes %in% res_classes)
  eucast_breakpoint = min(which_res)
  mic_k = array(dim=c(G,L,K))
  for(g in 1:G) {
    for(l in 1:L) {
      mic_k[g,l,] = unlist(data_selected[data_selected$gender_or==groups[[g]] & data_selected$year==years[l],"n"])
    }
  }
  mic_n = matrix(nrow=G,ncol=L)
  for(g in 1:G) mic_n[g,] = unlist(data_selected[data_selected$gender_or == groups[[g]] & data_selected$mic == mics[[1]],"sample_size"])
  res_k = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_k[i,] = unlist(datares[datares$gender_or == groups[[i]],"n"])
  res_n = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_n[i,] = unlist(datares[datares$gender_or == groups[[i]],"sample_size"])
  K_init = max(which(apply(mic_k[,1,],2,sum)>0))
  
  data_list = list(
    # simulation
    t0=t0, # start_year
    S=S, # duration of simulation
    ts=period_simulated, # time bins
    # structure
    C=C, # number of compartments
    L=L, # length of data
    G=G, # number of groups
    K=K, # number of MIC classes
    K_init=K_init, # number of MIC classes circulating at t0
    # data
    mic_k=mic_k, # number of individuals with resistance
    mic_n=mic_n,
    res_k=res_k,
    res_n=res_n,
    treatment_use=treatment_use, # parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu)
    eucast_breakpoint=eucast_breakpoint,
    # priors
    p_prevalence=target_prevalence, 
    p_tau=1,
    p_nu=c(2,8),
    p_mu=c(1,100),
    p_epsilon=c(1,1),
    p_init_mic=rep(1,K_init),
    # prediction
    pred_sample_size=mic_n[,L],
    max_threshold=20, # maximal threshold
    inference=inference
  )
  samples = stan("models/multistep_model.stan",data=data_list,chains=4,iter=2000,init=0.5,control = list(adapt_delta=.9,max_treedepth=10))
  print(samples)
  check_hmc_diagnostics(samples)
  return(list(samples=samples,data_list=data_list,data_selected=data_selected))
}


fit_multistep_model_omit = function(mic_data, groups, period_simulated, omit, target_prevalence, treatment_use, inference=1 ) {
  rstan_options(auto_write = TRUE,javascript=FALSE)
  data_selected = filter(mic_data,gender_or %in% groups) %>%
    arrange(gender_or,year,mic) %>%
    group_by(gender_or,year,mic,sample_size,atb,res) %>%
    summarise(n=sum(n)) %>%
    filter(!(year %in% omit) )
  datares = filter(mic_data,gender_or %in% groups) %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1) %>%
    filter(!(year %in% omit) )
  period_simulated = period_simulated[-which(period_simulated %in% omit)]
  t0=period_simulated[1]-1
  S=length(period_simulated)
  mics = unique(data_selected$mic)
  K = length(mics)
  C = K+1
  years = unique(data_selected$year)
  L = length(years)
  G = length(groups)
  mic_classes = unique(data_selected$mic)
  res_classes = unique(filter(data_selected,res==1)$mic)
  which_res = which(mic_classes %in% res_classes)
  eucast_breakpoint = min(which_res)
  mic_k = array(dim=c(G,L,K))
  for(g in 1:G) {
    for(l in 1:L) {
      mic_k[g,l,] = unlist(data_selected[data_selected$gender_or==groups[[g]] & data_selected$year==years[l],"n"])
    }
  }
  mic_n = matrix(nrow=G,ncol=L)
  for(g in 1:G) mic_n[g,] = unlist(data_selected[data_selected$gender_or == groups[[g]] & data_selected$mic == mics[[1]],"sample_size"])
  res_k = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_k[i,] = unlist(datares[datares$gender_or == groups[[i]],"n"])
  res_n = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_n[i,] = unlist(datares[datares$gender_or == groups[[i]],"sample_size"])
  K_init = max(which(apply(mic_k[,1,],2,sum)>0))
  
  data_list = list(
    # simulation
    t0=t0, # start_year
    S=S, # duration of simulation
    ts=period_simulated, # time bins
    # structure
    C=C, # number of compartments
    L=L, # length of data
    G=G, # number of groups
    K=K, # number of MIC classes
    K_init=K_init, # number of MIC classes circulating at t0
    # data
    mic_k=mic_k, # number of individuals with resistance
    mic_n=mic_n,
    res_k=res_k,
    res_n=res_n,
    treatment_use=treatment_use, # parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu)
    eucast_breakpoint=eucast_breakpoint,
    # priors
    p_prevalence=target_prevalence, 
    p_tau=1,
    p_nu=c(2,8),
    p_mu=c(1,100),
    p_epsilon=c(1,1),
    p_init_mic=rep(1,K_init),
    # prediction
    pred_sample_size=mic_n[,L],
    max_threshold=20, # maximal threshold
    inference=inference
  )
  samples = stan("models/multistep_model.stan",data=data_list,chains=4,iter=2000,init=0.5,control = list(adapt_delta=.9,max_treedepth=10))
  print(samples)
  check_hmc_diagnostics(samples)
  return(list(samples=samples,data_list=data_list,data_selected=data_selected))
}


fit_multistep_model_incprev = function(mic_data, groups, period_simulated, target_prevalence, treatment_use, zeta, inference=1 ) {
  rstan_options(auto_write = TRUE,javascript=FALSE)
  data_selected = filter(mic_data,gender_or %in% groups) %>%
    arrange(gender_or,year,mic) %>%
    group_by(gender_or,year,mic,sample_size,atb,res) %>%
    summarise(n=sum(n)) 
  datares = filter(mic_data,gender_or %in% groups) %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  t0=period_simulated[1]-1
  S=length(period_simulated)
  mics = unique(data_selected$mic)
  K = length(mics)
  C = K+1
  years = unique(data_selected$year)
  L = length(years)
  G = length(groups)
  mic_classes = unique(data_selected$mic)
  res_classes = unique(filter(data_selected,res==1)$mic)
  which_res = which(mic_classes %in% res_classes)
  eucast_breakpoint = min(which_res)
  mic_k = array(dim=c(G,L,K))
  for(g in 1:G) {
    for(l in 1:L) {
      mic_k[g,l,] = unlist(data_selected[data_selected$gender_or==groups[[g]] & data_selected$year==years[l],"n"])
    }
  }
  mic_n = matrix(nrow=G,ncol=L)
  for(g in 1:G) mic_n[g,] = unlist(data_selected[data_selected$gender_or == groups[[g]] & data_selected$mic == mics[[1]],"sample_size"])
  res_k = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_k[i,] = unlist(datares[datares$gender_or == groups[[i]],"n"])
  res_n = matrix(nrow=G,ncol=L)
  for(i in 1:G) res_n[i,] = unlist(datares[datares$gender_or == groups[[i]],"sample_size"])
  K_init = max(which(apply(mic_k[,1,],2,sum)>0))
  
  data_list = list(
    # simulation
    t0=t0, # start_year
    S=S, # duration of simulation
    ts=period_simulated, # time bins
    # structure
    C=C, # number of compartments
    L=L, # length of data
    G=G, # number of groups
    K=K, # number of MIC classes
    K_init=K_init, # number of MIC classes circulating at t0
    # data
    mic_k=mic_k, # number of individuals with resistance
    mic_n=mic_n,
    res_k=res_k,
    res_n=res_n,
    treatment_use=cbind(treatment_use,zeta), # parameters of the functional form of antibiotic use (t1,t2,delta,iota1,iota2,nu) and zeta
    eucast_breakpoint=eucast_breakpoint,
    # priors
    p_prevalence=target_prevalence, 
    p_tau=1,
    p_nu=c(2,8),
    p_mu=c(1,100),
    p_epsilon=c(1,1),
    p_init_mic=rep(1,K_init),
    # prediction
    pred_sample_size=mic_n[,L],
    max_threshold=20, # maximal threshold
    inference=inference
  )
  samples = stan("models/multistep_modelB.stan",data=data_list,chains=4,iter=2000,init=0.5,control = list(adapt_delta=.9,max_treedepth=10))
  print(samples)
  check_hmc_diagnostics(samples)
  return(list(samples=samples,data_list=data_list,data_selected=data_selected))
}