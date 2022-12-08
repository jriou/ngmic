
# Comparison R and Stan
source("R/R_multistep.R")

## get posterior samples
(load("models/samples_2022-02-09/S_multistep_grasp_ceftriaxone_2022-02-09.Rdata"))
samples = S_multistep_grasp_ceftriaxone$samples
test_sample = rstan::extract(samples)
test_nsample = 1
test_group = 1

## get fixed data
data_list = S_multistep_grasp_ceftriaxone$data_list
period_simulated= data_list$ts
treatment_use =  data_list$treatment_use[test_group,]
ts = period_simulated

## select one sample
parameters = list(
  beta=test_sample$beta[test_nsample,test_group],
  mu=test_sample$mu[test_nsample], 
  nu=test_sample$nu[test_nsample], 
  epsilon=test_sample$epsilon[test_nsample], 
  tau=test_sample$tau[test_nsample,test_group], 
  t1=treatment_use[1], 
  t2=treatment_use[2],
  delta=treatment_use[3],
  iota1=treatment_use[4],
  iota2=treatment_use[5],
  xi=treatment_use[6])
init_mic = test_sample$init_mic[test_nsample,test_group,]
target_prevalence = test_sample$prevalence[test_nsample,test_group]
init = c(S=1-target_prevalence,
         I1=target_prevalence*init_mic[1],
         I2=target_prevalence*init_mic[2],
         I3=target_prevalence*init_mic[3],
         I4=target_prevalence*init_mic[4],
         I5=target_prevalence*init_mic[5],
         I6=target_prevalence*init_mic[6],
         I7=target_prevalence*init_mic[7],
         I8=target_prevalence*init_mic[8])

out = ode(y=init,times=ts,multistep_model,parms=parameters)

out %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  pivot_longer(3:10) %>% 
  mutate(value=value/(1-S)+1e-7) %>% 
  filter(time>0) %>% 
  ggplot() +
  geom_line(aes(x=time,y=value,colour=name)) +
  facet_wrap(~name,nrow=1) +
  labs(y="proportion",colour="MIC class")

# check with stan output
test_sample$y[test_nsample,test_group,,] %>% 
  as_tibble() %>% 
  mutate(time=ts) %>% 
  pivot_longer(2:9) %>% 
  mutate(value=value/(1-V1)) %>% 
  ggplot() +
  geom_line(aes(x=time,y=value,colour=name)) +
  facet_wrap(~name,nrow=1) +
  labs(y="proportion",colour="MIC class")




