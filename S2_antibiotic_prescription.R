#' ---
#' title: "Additional file 2: antibiotic prescription"
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

#' ## Antibiotic usage by prescription
#' 
#' Six main combinations of antibiotics prescribed between 2000 and 2018.

#+ prescription_data, fig.width=8, fig.height=3, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE
ggplot(grasp_prescriptions) +
  geom_bar(aes(x=year,y=n/sample_size,group=treatment,fill=treatment),stat="identity",colour="black",width=.8) +
  facet_wrap( ~ gender_or) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent) +
  labs(x="Year",y="Proportion",fill="Therapy")

#' ## Probability of usage for each of the antibiotics
#' 
#' Data on antibiotic prescription can be transformed into the probability of usage of each of the four antibiotics individually.
#' The inverse triangles at the top indicate the periods of time during which each antibiotic was recommended
#' as a first-line treatment against *N. gonorrhoeae* infection in the UK.
#' 
#+ ind_prescription_data, fig.width=8, fig.height=5, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE
grasp_pred = grasp_prescriptions_by_atb %>%
  mutate(atb=factor(atb,levels=c(c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone"))))
recommendations = data.frame(atb=factor(c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")),
                             atbin=c(1900,2005,2010,2010),
                             atbout=c(2005,2010,2019,2100))
ggplot() +
  geom_bar(data=grasp_pred,aes(x=year,y=p,fill=atb),colour="black",stat="identity",width=1) +
  geom_segment(data=recommendations,aes(x=atbin,xend=atbout,y=1,yend=1),size=.5,colour="black") +
  geom_point(data=recommendations,aes(x=atbin,y=1),shape=25,fill="white",colour="black",size=2) +
  geom_point(data=recommendations,aes(x=atbout,y=1),shape=25,fill="white",colour="black",size=2) +
  facet_grid(gender_or ~ atb) +
  scale_y_continuous(labels=scales::percent,limits=c(0,1),expand=expansion(mult=c(0,0.05))) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide=FALSE) +
  coord_cartesian(xlim=c(2000,2019)) +
  labs(x="Year",y="Proportion")


#' ## Functional form
#'
#' One key aspect of modelling the dynamics of resistance in *N. gonorrhoeae* is to account
#' for the selection pressure exerted by antibiotics, and thus to include this information
#' about the probabity of usage of each antibiotic. However, this data in discrete by year,
#' which does not adapt well to ordinary differential equation-based models, that are continuous
#' in nature.
#' 
#' To solve this issue, we chose to use "forcing functions", that is continuous approximations
#' of the data using a closed-form function. We considered using splines or polynomials, but
#' decided on a "two-waves logistic growth and degrowth function" of the following form:
#' $$ p(t) = \frac{\iota_1}{1+\exp(-\delta*(t-t_1+\xi))} + \frac{1-\iota_1}{1+\exp(-\delta*(t-t_1))} + \frac{\iota_2}{1+\exp(\delta*(t-t_2-\xi))} + \frac{1-\iota_2}{1+\exp(\delta*(t-t_2))} - 1$$
#' 
#' This seemingly complicated formula is simply a weighted sum of four logistic functions,
#' well-suited to approximate proportions (range 0-1) that increase and decrease over time in
#' response to changing recommendations, sometimes in two waves as for ceftriaxone.
#' 
#' Within this formulation:
#' 
#' * $\delta$ is the slope of increase or decrease;
#' 
#' * $t_1$ and $t_2$ are median points of introduction or replacement in antibiotic usage; 
#' 
#' * (if applicable) if a two-wave introduction is needed, $\xi$ is the delay between the two waves and $\iota$ is the intermediate point.
#' 
#' For azithromycin, we fixed $t_2$ to 2022 to account for the change in recommendations in 2019 and
#' the expected diminution of usage from this date (we assume that the slope of diminution will mirror
#' the slope of increase around 2005).
#' For ciprofloxacin, we fixed $t_1$ at a very early date as ciprofloxacin was in use at the start of the study period.
#' For the opposite reason, we fixed  $t_2$ at a very late date for ceftriaxone.
#' The other parameters are estimated from data with least squares.
#' The R code for this procedure is available in `model-dev.R`.
#' This procedure results in the following functions, that are included in the models.

#+ forcing_functions, fig.width=8, fig.height=6, warning=FALSE, echo=FALSE, message=FALSE, results=FALSE
load("data/pred_doublelog_coef.Rdata")
atbswitch4 = function(x,t1,t2,delta,iota1,iota2,xi)  (iota1/(1+exp(-delta*(x-t1+xi))) + (1-iota1)/(1+exp(-delta*(x-t1))) + 
                                                        iota2/(1+exp(delta*(x-t2-xi))) + (1-iota2)/(1+exp(delta*(x-t2)))) - 1
expand.grid(atb=unique(grasp_prescriptions_by_atb$atb),gender_or=unique(grasp_prescriptions_by_atb$gender_or),year=seq(2000,2018,by=0.01)) %>%
  tbl_df() %>%
  left_join(pred_doublelog_coef, by = c("atb", "gender_or")) %>%
  mutate(pred_p=atbswitch4(x=year,t1=t1,t2=t2,delta=delta,iota1=iota1,iota2=iota2,xi=xi)) %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone"))) %>%
  ggplot() +
  geom_point(data=grasp_pred,aes(x=year,y=p),shape=21) +
  geom_line(aes(x=year,y=pred_p,colour=atb),size=1) +
  scale_colour_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide=FALSE) +
  scale_y_continuous(labels=scales::percent) +
  facet_grid(gender_or~atb) +
  labs(x="Year",y="Proportion",colour=NULL)

# rmarkdown::render("S2_antibiotic_prescription.R")
