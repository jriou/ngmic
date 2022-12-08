source("setup.R")
source("data_management.R") 
source("fit_multistep_model.R")
load("data/pred_doublelog_coef.Rdata")

# scp ngres/main-multistep3.R  UBELIX:projects/ngres/.main-multistep3.R
# scp ngres/models/sb  UBELIX:projects/ngres/.main-multistep3.R

# GRASP #####################################################################


# Ceftriaxone

SIM_multistep_grasp_ceftriaxone = fit_multistep_model(
  mic_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8], 
  inference=0
)

# plot_summary2(SIM_multistep_grasp_ceftriaxone,lim=2050,legend.pos=c(.8,.7),colmic = "Blues")

# check_hmc_diagnostics(SIM_multistep_grasp_ceftriaxone$samples)
# print(SIM_multistep_grasp_ceftriaxone$samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"),digits_summary = 4)
# print(SIM_multistep_grasp_ceftriaxone$samples,pars=c("output_res"),digits_summary = 4)
# 
# plot_summary(SIM_multistep_grasp_ceftriaxone)
# 

S_multistep_grasp_ceftriaxone = fit_multistep_model(
  mic_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8], 
  inference=1
)
# check_hmc_diagnostics(S_multistep_grasp_ceftriaxone$samples)
# print(S_multistep_grasp_ceftriaxone$samples,pars=c("mu","epsilon","prevalence","init_mic","nu","beta","tau"),digits_summary = 4)
# print(S_multistep_grasp_ceftriaxone$samples,pars=c("output_res"),digits_summary = 4)
# 
# plot_summary(S_multistep_grasp_ceftriaxone,lim=2020)
save(SIM_multistep_grasp_ceftriaxone,S_multistep_grasp_ceftriaxone,file=paste0("models/S_multistep_grasp_ceftriaxone",datetag,".Rdata"))



# Azithromycin

SIM_multistep_grasp_azithro = fit_multistep_model(
  mic_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8], 
  inference=0
)
# check_hmc_diagnostics(SIM_multistep_grasp_azithro$samples)
# print(SIM_multistep_grasp_azithro$samples,pars=c("mu","epsilon","prevalence","nu","beta","tau"),digits_summary = 4)
# 
# plot_summary3(SIM_multistep_grasp_azithro$samples,SIM_multistep_grasp_azithro$data_list,data_selected=grasp_azithromycin_mic,lim=2050,colmic = "Reds")

S_multistep_grasp_azithro = fit_multistep_model(
  mic_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8], 
  inference=1
)
# check_hmc_diagnostics(S_multistep_grasp_azithro$samples)
# print(S_multistep_grasp_azithro$samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"),digits_summary = 4)
# print(S_multistep_grasp_azithro$samples,pars=c("output_res"),digits_summary = 4)
# 
# plot_summary(S_multistep_grasp_azithro,lim=2050)
save(SIM_multistep_grasp_azithro,S_multistep_grasp_azithro,file=paste0("models/S_multistep_grasp_azithro",datetag,".Rdata"))



# Cefixime

SIM_multistep_grasp_cefixime = fit_multistep_model(
  mic_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8], 
  inference=0
)
# check_hmc_diagnostics(SIM_multistep_grasp_cefixime$samples)
# print(SIM_multistep_grasp_cefixime$samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"),digits_summary = 4)
# print(SIM_multistep_grasp_cefixime$samples,pars=c("output_res"),digits_summary = 4)
# 
# plot_summary(SIM_multistep_grasp_cefixime)


S_multistep_grasp_cefixime = fit_multistep_model(
  mic_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8], 
  inference=1
)
# check_hmc_diagnostics(S_multistep_grasp_cefixime$samples)
# print(S_multistep_grasp_cefixime$samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"),digits_summary = 4)
# print(S_multistep_grasp_cefixime$samples,pars=c("output_res"),digits_summary = 4)
# 
# plot_summary(S_multistep_grasp_cefixime,lim=2050)
save(SIM_multistep_grasp_cefixime,S_multistep_grasp_cefixime,file=paste0("models/S_multistep_grasp_cefixime",datetag,".Rdata"))




