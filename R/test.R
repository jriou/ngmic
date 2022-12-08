source("setup.R")

load("models/samples_2022-02-22/S_multistep_grasp_ceftriaxone_2022-02-22.Rdata")

S_multistep_grasp_ceftriaxone2 = S_multistep_grasp_ceftriaxone
SIM_multistep_grasp_ceftriaxone2 = SIM_multistep_grasp_ceftriaxone


(load("models/samples_2022-02-10/S_multistep_grasp_ceftriaxone_2022-02-10.Rdata"))



plot_summary2(SIM_multistep_grasp_ceftriaxone)
plot_summary2(SIM_multistep_grasp_ceftriaxone2)


plot_summary2(S_multistep_grasp_ceftriaxone)
plot_summary2(S_multistep_grasp_ceftriaxone2)

print(S_multistep_grasp_ceftriaxone$samples,pars=c("prevalence","mu","beta","epsilon","tau","nu","phi"))
print(S_multistep_grasp_ceftriaxone2$samples,pars=c("prevalence","mu","beta","epsilon","tau","nu","phi"))

S_multistep_grasp_ceftriaxone$data_list
S_multistep_grasp_ceftriaxone2$data_list

check_hmc_diagnostics(S_multistep_grasp_ceftriaxone$samples)
check_hmc_diagnostics(S_multistep_grasp_ceftriaxone2$samples)

pairs(SIM_multistep_grasp_ceftriaxone2$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))

pairs(S_multistep_grasp_ceftriaxone$samples,pars=c("mu","epsilon","tau","nu","phi"))
pairs(S_multistep_grasp_ceftriaxone2$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))


print(S_multistep_grasp_ceftriaxone$samples,pars=c("prevalence","mu","beta","epsilon","tau","nu","phi"))
print(S_multistep_grasp_ceftriaxone2$samples,pars=c("prevalence","mu","beta","epsilon","tau","nu","phi"))


stan_trace(S_multistep_grasp_ceftriaxone2$samples,pars=c("tau_raw","nu","mu","epsilon","nu"),inc_warmup = TRUE)
pairs(S_multistep_grasp_ceftriaxone2$samples,pars=c("tau","nu","mu","epsilon","nu"))
pairs(S_multistep_grasp_ceftriaxone2$samples,pars=c("init_mic"))
pairs(S_multistep_grasp_ceftriaxone2$samples,pars=c("init_mic"))
print(S_multistep_grasp_ceftriaxone2$samples,pars=c("tau_raw","nu","mu","epsilon","nu"))

check_hmc_diagnostics(S_multistep_grasp_ceftriaxone$samples)
check_hmc_diagnostics(S_multistep_grasp_ceftriaxone2$samples)

plot_summary2(S_multistep_grasp_ceftriaxone,lim=2050,legend.pos=c(.8,.7),colmic = "Blues")
pairs(SIM_multistep_grasp_ceftriaxone$samples,pars=c("nu","tau_raw","mu","epsilon","phi"))
stan_trace(SIM_multistep_grasp_ceftriaxone$samples,pars=c("mu","beta","epsilon","tau","nu","phi"),inc_warmup = TRUE)
stan_dens(SIM_multistep_grasp_ceftriaxone$samples,pars=c("mu","epsilon","tau","nu","phi"))

shinystan::launch_shinystan(S_multistep_grasp_ceftriaxone2$samples)


