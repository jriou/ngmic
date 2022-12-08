source("setup.R")
source("data_management.R") 
source("fit_binary_model.R")
load("data/pred_doublelog_coef.Rdata")


# GRASP #####################################################################

# Azithromycin


SIM_binary_grasp_azithro_flat = fit_binary_model_flat(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8], 
  inference=0
)

save(SIM_binary_grasp_azithro_flat,file=paste0("models/S_binary_grasp_azithro_flat",datetag,".Rdata"))



SIM_binary_grasp_azithro = fit_binary_model(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8], 
  inference=0
)


S_binary_grasp_azithro = fit_binary_model(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8], 
  inference=1
)

save(SIM_binary_grasp_azithro,S_binary_grasp_azithro,file=paste0("models/S_binary_grasp_azithro",datetag,".Rdata"))


# Ciprofloxacin

SIM_binary_grasp_ciprofloxacin = fit_binary_model(
  res_data=grasp_ciprofloxacin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2000:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(2,6),3:8], 
  inference=0
)

S_binary_grasp_ciprofloxacin = fit_binary_model(
  res_data=grasp_ciprofloxacin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2000:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(2,6),3:8], 
  inference=1
)

save(SIM_binary_grasp_ciprofloxacin,S_binary_grasp_ciprofloxacin,file=paste0("models/S_binary_grasp_ciprofloxacin",datetag,".Rdata"))


# cefixime

SIM_binary_grasp_cefixime = fit_binary_model(
  res_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8], 
  inference=0
)

S_binary_grasp_cefixime = fit_binary_model(
  res_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8], 
  inference=1
)

save(SIM_binary_grasp_cefixime,S_binary_grasp_cefixime,file=paste0("models/S_binary_grasp_cefixime",datetag,".Rdata"))




# ceftriaxone

SIM_binary_grasp_ceftriaxone = fit_binary_model(
  res_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8], 
  inference=0
)

S_binary_grasp_ceftriaxone = fit_binary_model(
  res_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8], 
  inference=1
)
save(SIM_binary_grasp_ceftriaxone,S_binary_grasp_ceftriaxone,file=paste0("models/S_binary_grasp_ceftriaxone",datetag,".Rdata"))




