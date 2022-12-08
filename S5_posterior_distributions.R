#' ---
#' title: "Additional file 5: posterior distributions"
#' author: ""
#' date: ""
#' header-includes: 
#'  \usepackage{tikz}
#'  \usepackage{pgfplots}
#' output:
#'    html_document:
#'      theme: "journal"
#'      toc: true
#'      toc_depth: 3
#'      toc_float: true
#'      toc_collapsed: false
#'      code_folding : hide
#' ---



#+ include=FALSE
source("setup.R")
source("data_management.R") 
source("fit_multistep_model.R")
load("data/pred_doublelog_coef.Rdata")


#' We present here the posterior distributions of all parameters, first separately for 
#' each combination of model and antibiotic, then all together.
#' 
#' ## Ciprofloxacin 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_ciprofloxacin_2022-05-17.Rdata"))

#' ### Single-step model
post_summary(S_binary_grasp_ciprofloxacin)

#' ### Multi-step model
#' 
#' Not done because of the lack of detail on MIC after 2009 in the data
#' 

#' ## Cefixime 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_cefixime_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_cefixime_2022-05-17.Rdata"))

#' ### Single-step model
post_summary(S_binary_grasp_cefixime)

#' ### Multi-step model
post_summary(S_multistep_grasp_cefixime)


#' ## Azithromycin 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_azithro_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_azithro_2022-05-17.Rdata"))

#' ### Single-step model
post_summary(S_binary_grasp_azithro)

#' ### Multi-step model
post_summary(S_multistep_grasp_azithro)



#' ## Ceftriaxone 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_ceftriaxone_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_2022-05-17.Rdata"))

#' ### Single-step model
post_summary(S_binary_grasp_ceftriaxone)

#' ### Multi-step model
post_summary(S_multistep_grasp_ceftriaxone)



#' ## Summary of posterior distributions
#' 
#' 
#' ### Single-step models
#+ fig.width=8, fig.height=8

atbcols = RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)]
bind_rows(
  cbind(post_summary(S_binary_grasp_azithro,FALSE),atb="Azithromycin"),
  cbind(post_summary(S_binary_grasp_ciprofloxacin,FALSE),atb="Ciprofloxacin"),
  cbind(post_summary(S_binary_grasp_cefixime,FALSE),atb="Cefixime"),
  cbind(post_summary(S_binary_grasp_ceftriaxone,FALSE),atb="Ceftriaxone")
) %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone"))) %>%
  ggplot() +
  geom_pointrange(aes(x=atb,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=atb),position=position_dodge(.7)) +
  theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_colour_manual(values=atbcols,guide="none") +
  facet_wrap(~comment,scale="free",ncol=4) +
  labs(x="Antibiotics",y="Posterior value")

#' ### Multi-step models
#+ fig.width=8, fig.height=8
atbcols = atbcols[-1]
bind_rows(cbind(post_summary(S_multistep_grasp_azithro,FALSE),atb="Azithromycin"),
          cbind(post_summary(S_multistep_grasp_cefixime,FALSE),atb="Cefixime"),
          cbind(post_summary(S_multistep_grasp_ceftriaxone,FALSE),atb="Ceftriaxone")) %>%
  mutate(atb=factor(atb,levels=c("Cefixime","Azithromycin","Ceftriaxone"))) %>% 
  ggplot() +
  geom_pointrange(aes(x=atb,y=`50%`,ymin=`2.5%`,ymax=`97.5%`,colour=atb),position=position_dodge(.7)) +
  theme(axis.text = element_text(angle=45,hjust=1)) +
  scale_colour_manual(values=atbcols,guide="none") +
  facet_wrap(~comment,scale="free",ncol=4) +
  labs(x="Antibiotics",y="Posterior value")


# rmarkdown::render("S5_posterior_distributions.R")
