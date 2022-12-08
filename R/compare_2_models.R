
(load("models/samples_2022-02-10/S_multistep_grasp_ceftriaxone_2022-02-10.Rdata"))

S_multistep_grasp_ceftriaxone2 = S_multistep_grasp_ceftriaxone
SIM_multistep_grasp_ceftriaxone2 = SIM_multistep_grasp_ceftriaxone


(load("models/samples_2022-02-09/S_multistep_grasp_ceftriaxone_2022-02-09.Rdata"))



plot_summary2(SIM_multistep_grasp_ceftriaxone)
plot_summary2(SIM_multistep_grasp_ceftriaxone2)


plot_summary2(S_multistep_grasp_ceftriaxone)
plot_summary2(S_multistep_grasp_ceftriaxone2)

check_hmc_diagnostics(S_multistep_grasp_ceftriaxone$samples)
check_hmc_diagnostics(S_multistep_grasp_ceftriaxone2$samples)

pairs(SIM_multistep_grasp_ceftriaxone$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))

pairs(S_multistep_grasp_ceftriaxone$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))
pairs(S_multistep_grasp_ceftriaxone2$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))

print(S_multistep_grasp_ceftriaxone$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))
print(S_multistep_grasp_ceftriaxone2$samples,pars=c("mu","beta","epsilon","tau","nu","phi"))
