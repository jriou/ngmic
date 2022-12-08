#' ---
#' title: "Additional file 4: model fit and predictions"
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


#' We present here the full results from the application of the single-step and multi-step models to data from GRASP.
#' All the posterior samples used for this work are available in `models`.
#' Files bigger than 50 Mb were split into several `.7z` archives for storing on github.
#' 
#' For each model and each antibiotic, we present a figure summarizing different aspects of the model results.
#' Panel (A) shows the model fit; panel (B) is the evolution of treatment usage with time; panel (C) is the model-based projection 
#' of NG antimicrobial resistance (using the appropriate EUCAST threshold) until 2030; panel (D) is the projection of NG prevalence if influenced only by resistance levels; and (E) is the evolution of the 
#' probability of reaching the 5%, 10% and 20% thresholds. We present results from the single-step model, and when 
#' applicable (that is, for all antibiotics except ciprofloxacin), results from the multi-step model.
#' 
#' ## Ciprofloxacin 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_ciprofloxacin_2022-05-17.Rdata"))

#' ### Single-step model
#+ fig.width=8, fig.height=8
plot_summary(S_binary_grasp_ciprofloxacin,lim=2031,legend.pos=c(.8,.3),colmic = "Purples")

#' ### Multi-step model
#' 
#' Not done because of the lack of detail on MIC after 2009 in the data
#' 

#' ## Cefixime 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_cefixime_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_cefixime_2022-05-17.Rdata"))

#' ### Single-step model
#+ fig.width=8, fig.height=8
plot_summary(S_binary_grasp_cefixime,lim=2031,legend.pos=c(.8,.7),colmic = "Greens")

#' ### Multi-step model
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_cefixime,lim=2031,legend.pos=c(.8,.7),colmic = "Greens")


#' ## Azithromycin 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_azithro_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_azithro_2022-05-17.Rdata"))

#' ### Single-step model
#+ fig.width=8, fig.height=8
plot_summary(S_binary_grasp_azithro,lim=2031,legend.pos=c(.4,.3),colmic = "Reds")

#' ### Multi-step model
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_azithro,lim=2031,legend.pos=c(.4,.3),colmic = "Reds")



#' ## Ceftriaxone 
#+ include=FALSE
(load("models/samples_2022-05-17/S_binary_grasp_ceftriaxone_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_ss1__2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_incprev_2022-05-17.Rdata"))

#' ### Single-step model
#+ fig.width=8, fig.height=8
plot_summary(S_binary_grasp_ceftriaxone,lim=2031,legend.pos=c(.8,.7),colmic = "Blues")

#' ### Multi-step model
#' 
#' #### With all data from 2004 to 2018
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_ceftriaxone,lim=2031,legend.pos=c(.8,.7),colmic = "Blues")

#' #### Sensitivity analysis removing data from 2009-2010
#' 
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_ceftriaxone_ss1,lim=2031,legend.pos=c(.8,.7),colmic = "Blues")

#' #### Sensitivity analysis with increasing transmission
#' ##### Yearly increase of $\beta$ by 1.001
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_ceftriaxone_incprev001,lim=2031,legend.pos=c(.8,.7),colmic = "Blues")

#' ##### Yearly increase of $\beta$ by 1.005
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_ceftriaxone_incprev005,lim=2031,legend.pos=c(.8,.7),colmic = "Blues")

#' ##### Yearly increase of $\beta$ by 1.01
#+ fig.width=8, fig.height=8
plot_summary2(S_multistep_grasp_ceftriaxone_incprev01,lim=2031,legend.pos=c(.8,.7),colmic = "Blues")



# rmarkdown::render("S4_model_fits_predictions.R")
