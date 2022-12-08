source("setup.R")
source("data_management.R") 
source("fit_binary_model.R")
load("data/pred_doublelog_coef.Rdata")


# GRASP #####################################################################

# Azithromycin

S_binary_grasp_azithro_incprev001 = fit_binary_model_incprev(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8],
  zeta = 0.001,
  inference=1
)

S_binary_grasp_azithro_incprev005 = fit_binary_model_incprev(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8],
  zeta = 0.005,
  inference=1
)

S_binary_grasp_azithro_incprev01 = fit_binary_model_incprev(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8],
  zeta = 0.01,
  inference=1
)

save(S_binary_grasp_azithro_incprev001,S_binary_grasp_azithro_incprev005,S_binary_grasp_azithro_incprev01,file=paste0("models/S_binary_grasp_azithro_incprev",datetag,".Rdata"))


# Ciprofloxacin

S_binary_grasp_ciprofloxacin_incprev001 = fit_binary_model_incprev(
  res_data=grasp_ciprofloxacin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2000:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(2,6),3:8],
  zeta = 0.001,
  inference=1
)
S_binary_grasp_ciprofloxacin_incprev005 = fit_binary_model_incprev(
  res_data=grasp_ciprofloxacin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2000:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(2,6),3:8],
  zeta = 0.005,
  inference=1
)

S_binary_grasp_ciprofloxacin_incprev01 = fit_binary_model_incprev(
  res_data=grasp_ciprofloxacin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2000:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(2,6),3:8],
  zeta = 0.01,
  inference=1
)
save(S_binary_grasp_ciprofloxacin_incprev001,S_binary_grasp_ciprofloxacin_incprev005,S_binary_grasp_ciprofloxacin_incprev01,file=paste0("models/S_binary_grasp_ciprofloxacin_incprev",datetag,".Rdata"))


# cefixime
S_binary_grasp_cefixime_incprev001 = fit_binary_model_incprev(
  res_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8],
  zeta = 0.001,
  inference=1
)
S_binary_grasp_cefixime_incprev005 = fit_binary_model_incprev(
  res_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8],
  zeta = 0.005,
  inference=1
)
S_binary_grasp_cefixime_incprev01 = fit_binary_model_incprev(
  res_data=grasp_cefixime_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(3,7),3:8],
  zeta = 0.01,
  inference=1
)
save(S_binary_grasp_cefixime_incprev001,S_binary_grasp_cefixime_incprev005,S_binary_grasp_cefixime_incprev01,file=paste0("models/S_binary_grasp_cefixime_incprev",datetag,".Rdata"))




# ceftriaxone
S_binary_grasp_ceftriaxone_incprev001 = fit_binary_model_incprev(
  res_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8],
  zeta = 0.001,
  inference=1
)
S_binary_grasp_ceftriaxone_incprev005 = fit_binary_model_incprev(
  res_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8],
  zeta = 0.005,
  inference=1
)
S_binary_grasp_ceftriaxone_incprev01 = fit_binary_model_incprev(
  res_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8],
  zeta = 0.01,
  inference=1
)
save(S_binary_grasp_ceftriaxone_incprev001,S_binary_grasp_ceftriaxone_incprev005,S_binary_grasp_ceftriaxone_incprev01,file=paste0("models/S_binary_grasp_ceftriaxone_incprev",datetag,".Rdata"))



