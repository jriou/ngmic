source("setup.R")
source("data_management.R") 
source("fit_multistep_model.R")
load("data/pred_doublelog_coef.Rdata")

# scp ngres/supp-multistep2.R  UBELIX:projects/ngres/.
# scp ngres/models/sb_multistep.sh  UBELIX:projects/ngres/.
# scp ngres/fit_multistep_model.R  UBELIX:projects/ngres/.

S_multistep_grasp_ceftriaxone_ss1 = fit_multistep_model_omit(
  mic_data=grasp_ceftriaxone_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2004:2050,
  omit=2009:2010,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                              ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(4,8),3:8], 
  inference=1
)

save(S_multistep_grasp_ceftriaxone_ss1,file=paste0("models/S_multistep_grasp_ceftriaxone_ss1_",datetag,".Rdata"))




