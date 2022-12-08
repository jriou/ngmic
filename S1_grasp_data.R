#' ---
#' title: "Additional file 1: GRASP data"
#' author: ""
#' date: ""
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
library(kableExtra)

#' ## Resistance data
#' 
#' Data from the GRASP network in the United Kingdom from 2000 to 2018. 
#' Susceptibility to antibiotics can be dichotomized into sensitive or resistant using the EUCAST
#' breakpoints (>1 mg/L for azithromycin, >0.06 mg/L for ciprofloxacin, >0.125 mg/L for cefixime, 
#' and >0.125 mg/L for ceftriaxone). This data is stratified by sex and sexual
#' orientation in two groups: heterosexual men and women (HMW) and men who have sex with men (MSM).
#' 
#+ allres, fig.width=8, fig.height=5, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE

res_data = bind_rows(grasp_azithromycin_mic,grasp_cefixime_mic,grasp_ciprofloxacin_mic,grasp_ceftriaxone_mic) %>%
  group_by(atb,EUCAST_breakpoint,year,gender_or,res) %>%
  summarise(n=sum(n)) %>%
  mutate(sample_size=sum(n),prev=n/sample_size) %>%
  filter(res==1) %>%
  ungroup()  %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")))
ggplot(res_data) +
  geom_line(data=res_data,aes(x=year,y=prev),size=.5) +
  geom_point(data=res_data,aes(x=year,y=prev,fill=atb),colour="black",shape=21,size=1.8) +
  facet_grid(gender_or ~ atb) +
  scale_y_continuous(labels=scales::percent,limits=c(0,.55)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide=FALSE) +
  coord_cartesian(xlim=c(2000,2019)) +
  labs(x="Year",y="Proportion resistant")

# #+ table1
# res_data %>%
#   dplyr::select(atb,year,gender_or,n,sample_size) %>%
#   pivot_wider(names_from=gender_or,values_from=c(n,sample_size)) %>%
#   dplyr::select(atb,year,n_HMW,sample_size_HMW,n_MSM,sample_size_MSM) %>%
#   arrange(atb,year) %>%
#   kbl(caption="Table. Number of isolates with resistance according to EUCAST thresholds by antibiotic and by gender and sexual orientation (GRASP, 2000-2018).",
#       col.names = c("Antibiotic","Year","N with resistance (HMW)","Sample size (HMW)","N with resistance (MSM)","Sample size (MSM)")) %>%
#   kable_styling()

#' ## MIC data
#' 
#' For each antibiotic, we also present the evolution of the MIC in time as a stacked area graph,
#' more adapted to observe time trends (the "MIC drift"). For visualization purposes, we group together 
#' MIC classes that correspond to resistance according to the EUCAST thresholds. 
#' 
#' ### Ciprofloxacin
#' 
#' The interpretation of ciprofloxacin MIC data is hindered by a technical change in 2010.
#' 
#+ ciprofloxacin_mic, fig.width=8, fig.height=3.5, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE

col1 = RColorBrewer::brewer.pal(7,"Purples")[7]
grasp_ciprofloxacin_mic %>%
  group_by(atb,year,gender_or,micres) %>%
  summarise(n=sum(n)) %>%
  mutate(sample_size=sum(n),type="MIC") %>%
  ungroup() %>%
  ggplot() +
  geom_area(aes(x=year,y=n/sample_size,group=micres,fill=micres),colour="black",size=0.3,alpha=.7) +
  scale_fill_brewer(type="seq",palette="Purples") +
  scale_x_continuous(expand=c(0,0),limits=c(2000,2018)) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1),breaks=c(.1,.3,.5,.7,.9)) +
  facet_grid( ~ gender_or) +
  labs(x=NULL,y="Proportion",fill="Ciprofloxacine MIC:") +
  theme(panel.spacing = unit(0.3, "lines"),legend.position="bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

#' ### Cefixime
#' 
#' We observe a slight general MIC drift for cefixime for HMW and MSM. This general trend is perturbed 
#' by an unusual downward shift of MIC levels in 2006, that remains unexplained. We also observe a temporary
#' upward shift of MIC levels around 2009-2010, mostly among MSM, that can be related to the circulation of an international
#' clone with a novel penA mosaic allele (Unemo et al. 2012).
#' 
#+ cefixime_mic, fig.width=8, fig.height=3.5, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE

col1 = RColorBrewer::brewer.pal(8,"Greens")[8]
grasp_cefixime_mic %>%
  group_by(atb,year,gender_or,micres) %>%
  summarise(n=sum(n)) %>%
  mutate(sample_size=sum(n),type="MIC") %>%
  ungroup() %>%
  ggplot() +
  geom_area(aes(x=year,y=n/sample_size,group=micres,fill=micres),colour="black",size=0.3,alpha=.7) +
  scale_fill_brewer(type="seq",palette="Greens") +
  scale_x_continuous(expand=c(0,0),limits=c(2000,2018)) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1),breaks=c(.1,.3,.5,.7,.9)) +
  facet_grid( ~ gender_or) +
  labs(x=NULL,y="Proportion",fill="Cefixime MIC:") +
  theme(panel.spacing = unit(0.3, "lines"),legend.position="bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))


#' ### Azithromycin
#' 
#' MIC levels for azithromycin appear relatively stable throughout the study period, with a sudden increase
#' around 2015. 
#' 
#+ azithromycin_mic, fig.width=8, fig.height=3.5, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE

col1 = RColorBrewer::brewer.pal(6,"Reds")[6]
grasp_azithromycin_mic %>%
  group_by(atb,year,gender_or,micres) %>%
  summarise(n=sum(n)) %>%
  mutate(sample_size=sum(n),type="MIC") %>%
  ungroup() %>%
  ggplot() +
  geom_area(aes(x=year,y=n/sample_size,group=micres,fill=micres),colour="black",size=0.3,alpha=.7) +
  scale_fill_brewer(type="seq",palette="Reds") +
  scale_x_continuous(expand=c(0,0),limits=c(2000,2018)) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1),breaks=c(.1,.3,.5,.7,.9)) +
  facet_grid( ~ gender_or) +
  labs(x=NULL,y="Proportion",fill="Azithromycin MIC:") +
  theme(panel.spacing = unit(0.3, "lines"),legend.position="bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

#' ### Ceftriaxone
#' 
#' We observe a continuous MIC drift for ceftriaxone in both populations from 2005 on. As for cefixime, 
#' we also observe a temporary upward shift of MIC levels around 2009-2010 among MSM, that could be related to the 
#' circulation of the international clone.
#' 
#+ ceftriaxone_mic, fig.width=8, fig.height=3.5, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE

grasp_ceftriaxone_mic %>%
  group_by(atb,year,gender_or,micres) %>%
  summarise(n=sum(n)) %>%
  mutate(sample_size=sum(n),type="MIC") %>%
  ungroup() %>%
  ggplot() +
  geom_area(aes(x=year,y=n/sample_size,group=micres,fill=micres),colour="black",size=0.3,alpha=.7) +
  scale_fill_brewer(type="seq",palette="Blues") +
  scale_x_continuous(expand=c(0,0),limits=c(2000,2018)) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1),breaks=c(.1,.3,.5,.7,.9)) +
  facet_grid( ~ gender_or) +
  labs(x=NULL,y="Proportion",fill="Ceftriaxone MIC:") +
  theme(panel.spacing = unit(0.3, "lines"),legend.position="bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))


# rmarkdown::render("S1_grasp_data.R")
