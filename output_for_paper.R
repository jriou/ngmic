# Setup

source("setup.R")
source("data_management.R") 
source("fit_multistep_model.R")
load("data/pred_doublelog_coef.Rdata")

(load("models/samples_2022-05-17/S_binary_grasp_azithro_2022-05-17.Rdata")) 
(load("models/samples_2022-05-17/S_binary_grasp_ciprofloxacin_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_binary_grasp_cefixime_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_binary_grasp_ceftriaxone_2022-05-17.Rdata"))

(load("models/samples_2022-05-17/S_multistep_grasp_azithro_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_cefixime_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_ss1__2022-05-17.Rdata"))
(load("models/samples_2022-05-17/S_multistep_grasp_ceftriaxone_incprev_2022-05-17.Rdata"))

# Text ----

## totals
grasp_prescriptions %>%
  group_by(gender_or) %>%
  summarise(n=sum(n))

bind_rows(grasp_azithromycin_mic,grasp_cefixime_mic,grasp_ceftriaxone_mic,grasp_ciprofloxacin_mic) %>%
  group_by(atb) %>%
  summarise(n=sum(n)) 

## ciprofloxacin

extract_pred(S_binary_grasp_ciprofloxacin,"Ciprofloxacin",2050) %>%
  filter(year %in% c(2000,2018)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

## Cefixime
grasp_cefixime_mic %>%
  group_by(year,gender_or,res) %>%
  summarise(n=sum(n)) %>%
  mutate(p=n/sum(n)) %>%
  filter(res==1) %>%
  as.data.frame()


extract_pred(S_binary_grasp_cefixime,"Cefixime",2050) %>%
  filter(year %in% c(2004,2010)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

extract_mic(S_multistep_grasp_cefixime)$res %>%
  filter(year %in% c(2004,2010)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

## Azithromycin

extract_pred(S_binary_grasp_azithro,"Azithromycin",2050) %>%
  filter(year %in% c(2001,2015,2018,2049)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

extract_mic(S_multistep_grasp_azithro)$res %>%
  filter(year %in% c(2001,2015,2018, 2049)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))


## Ceftriaxone

grasp_ceftriaxone_mic %>%
  filter(res==1) %>%
  filter(n>0)

extract_pred(S_binary_grasp_ceftriaxone,"Ceftriaxone",2050) %>%
  filter(year %in% c(2004,2018,2030)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

extract_thres(S_binary_grasp_ceftriaxone,"Ceftriaxone",2050) %>%
  filter(year%in%c(2030,2050))

extract_mic(S_multistep_grasp_ceftriaxone)$res %>%
  filter(year %in% c(2004,2018,2030)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

extract_mic(S_multistep_grasp_ceftriaxone)$thres %>%
  filter(year %in% c(2030)) 


## Ceftriaxone sensitivity analysis without 2009-2010

extract_mic(S_multistep_grasp_ceftriaxone_ss1)$res %>%
  filter(year %in% c(2004,2018,2030, 2049)) %>%
  transmute(year=year,
            gender_or=gender_or,
            post=format_post(`50%`,`2.5%`,`97.5%`))

extract_mic(S_multistep_grasp_ceftriaxone_ss1)$thres %>%
  filter(year %in% c(2030)) 

# Fig 1 ----

grasp_pred = grasp_prescriptions_by_atb %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")))
res_data = bind_rows(grasp_azithromycin_mic,grasp_cefixime_mic,grasp_ciprofloxacin_mic,grasp_ceftriaxone_mic) %>%
  group_by(atb,EUCAST_breakpoint,year,gender_or,res) %>%
  summarise(n=sum(n)) %>%
  mutate(sample_size=sum(n),prev=n/sample_size) %>%
  filter(res==1) %>%
  ungroup() %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")))
recommendations = data.frame(atb=factor(c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")),
                  atbin=c(1900,2005,2010,2010),
                  atbout=c(2005,2010,2019,2100))
ggplot() +
  geom_bar(data=grasp_pred,aes(x=year,y=p,fill=atb),colour="grey70",stat="identity",alpha=0.2,width=1) +
  geom_line(data=res_data,aes(x=year,y=prev),size=.5) +
  geom_point(data=res_data,aes(x=year,y=prev,fill=atb),colour="black",shape=21,size=1.8) +
  geom_segment(data=recommendations,aes(x=atbin,xend=atbout,y=1,yend=1),size=.5,colour="black") +
  geom_point(data=recommendations,aes(x=atbin,y=1),shape=25,fill="white",colour="black",size=2) +
  geom_point(data=recommendations,aes(x=atbout,y=1),shape=25,fill="white",colour="black",size=2) +
  facet_grid(gender_or ~ atb) +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide=FALSE) +
  scale_colour_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide=FALSE) +
  coord_cartesian(xlim=c(2000,2019)) +
  labs(x="Year",y="Proportion")
ggsave(file="figures/fig1.pdf",height=4,width=8)


# Fig 2 (using tex in /figures) ----

# Fig 3 ----


single_step_pred = bind_rows(
  extract_pred(S_binary_grasp_azithro,"Azithromycin",2019),
  extract_pred(S_binary_grasp_ciprofloxacin,"Ciprofloxacin",2019),
  extract_pred(S_binary_grasp_cefixime,"Cefixime",2019),
  extract_pred(S_binary_grasp_ceftriaxone,"Ceftriaxone",2031)) %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone"))) 

single_step_pred %>% 
  filter(year<2019) %>% 
  left_join(res_data) %>%
  mutate(gut=ifelse(prev<=`97.5%` & prev>=`2.5%`,1,0)) %>% 
  group_by(atb) %>% 
  summarise(gut=sum(gut,na.rm=TRUE),n=n())
  

multi_step_pred = bind_rows(
  cbind(extract_mic(S_multistep_grasp_azithro)$res,atb="Azithromycin"),
  cbind(extract_mic(S_multistep_grasp_cefixime)$res,atb="Cefixime"),
  cbind(extract_mic(S_multistep_grasp_ceftriaxone)$res,atb="Ceftriaxone")) %>%
  filter(!(year>=2019 & atb!="Ceftriaxone")) %>%
  filter(year<2031) %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")))

multi_step_pred %>% 
  filter(year<2019) %>% 
  left_join(res_data) %>%
  mutate(gut=ifelse(prev<=`97.5%` & prev>=`2.5%`,1,0)) %>% 
  group_by(atb) %>% 
  summarise(gut=sum(gut,na.rm=TRUE),n=n())

recommendations2 = recommendations %>% 
  mutate(atbin=ifelse(atbin==1900,2000,atbin)) %>% 
  mutate(atbout=ifelse(atbout==2100,2030,atbout))

predrect = expand_grid(year=2018:2032,atb=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone")) %>%
  mutate(atb=factor(atb,levels=c("Ciprofloxacin","Cefixime","Azithromycin","Ceftriaxone"))) %>% 
  filter(!(year>2022 & atb!="Ceftriaxone"))

ggplot(single_step_pred,aes(x=year)) +
  geom_hline(yintercept=.05,linetype=1,size=.2,col="orange") +
  geom_line(data=res_data,aes(x=year,y=prev),size=.3) +
  geom_point(data=res_data,aes(x=year,y=prev),fill="white",colour="black",shape=21,size=1.3) +
  geom_ribbon(data=predrect,aes(ymin=0,ymax=1),alpha=.1) +
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=atb),alpha=.6) +
  geom_ribbon(data=multi_step_pred,aes(ymin=`2.5%`,ymax=`97.5%`,fill=atb),alpha=.4) +
  geom_line(aes(y=`50%`)) +
  geom_line(data=multi_step_pred,aes(y=`50%`),linetype=2) +
  geom_segment(data=recommendations2,aes(x=atbin,xend=atbout,y=.55,yend=.55),size=.5,colour="black") +
  geom_point(data=recommendations2,aes(x=atbin,y=.55),shape=25,fill="white",colour="black",size=2) +
  geom_point(data=recommendations2,aes(x=atbout,y=.55),shape=25,fill="white",colour="black",size=2) +
  facet_grid(gender_or ~ atb,scales="free_x", space='free_x') +
  annotate(geom="point",x=2000,y=0,size=0.01,colour="transparent") + 
  annotate(geom="point",x=2020,y=0,size=0.01,colour="transparent") + 
  scale_y_continuous(labels=scales::percent,breaks=c(0,0.05,.2,.4,.6),minor_breaks = NULL) +
  scale_x_continuous(expand=expansion(add=c(1.5,0)),breaks=seq(2000,2030,by=10)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide="none") +
  scale_colour_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide="none") +
  coord_cartesian(ylim=c(0,0.55)) +
  labs(x="Year",y="Proportion")
ggsave(file="figures/fig3_v3.pdf",height=4,width=8)



g1 = ggplot(single_step_pred,aes(x=year)) +
  geom_hline(yintercept=.05,linetype=2,size=.2) +
  geom_line(data=res_data,aes(x=year,y=prev),size=.3) +
  geom_point(data=res_data,aes(x=year,y=prev),fill="white",colour="black",shape=21,size=1.3) +
  geom_ribbon(data=predrect,aes(ymin=0,ymax=1),alpha=.1) +
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=atb),alpha=.6) +
  geom_line(aes(y=`50%`)) +
  geom_segment(data=recommendations2,aes(x=atbin,xend=atbout,y=.55,yend=.55),size=.5,colour="black") +
  geom_point(data=recommendations2,aes(x=atbin,y=.55),shape=25,fill="white",colour="black",size=2) +
  geom_point(data=recommendations2,aes(x=atbout,y=.55),shape=25,fill="white",colour="black",size=2) +
  facet_grid(gender_or ~ atb,scales="free_x", space='free_x') +
  annotate(geom="point",x=2000,y=0,size=0.01,colour="transparent") + 
  annotate(geom="point",x=2020,y=0,size=0.01,colour="transparent") + 
  scale_y_continuous(labels=scales::percent,breaks=c(0,0.05,.2,.4,.6),minor_breaks = NULL) +
  scale_x_continuous(expand=expansion(add=c(1.5,0)),breaks=seq(2000,2030,by=10)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide="none") +
  scale_colour_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(4,3,1,2)],guide="none") +
  coord_cartesian(ylim=c(0,0.55)) +
  labs(x="Year",y="Proportion")
g2 = ggplot(multi_step_pred,aes(x=year)) +
  geom_hline(yintercept=.05,linetype=2,size=.2) +
  geom_line(data=res_data,aes(x=year,y=prev),size=.3) +
  geom_point(data=res_data,aes(x=year,y=prev),fill="white",colour="black",shape=21,size=1.3) +
  geom_ribbon(data=predrect,aes(ymin=0,ymax=1),alpha=.1) +
  geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=atb),alpha=.6) +
  # geom_ribbon(data=multi_step_pred,aes(ymin=`2.5%`,ymax=`97.5%`,fill=atb),alpha=.4) +
  geom_line(aes(y=`50%`)) +
  # geom_line(data=multi_step_pred,aes(y=`50%`),linetype=2) +
  geom_segment(data=recommendations2,aes(x=atbin,xend=atbout,y=.55,yend=.55),size=.5,colour="black") +
  geom_point(data=recommendations2,aes(x=atbin,y=.55),shape=25,fill="white",colour="black",size=2) +
  geom_point(data=recommendations2,aes(x=atbout,y=.55),shape=25,fill="white",colour="black",size=2) +
  facet_grid(gender_or ~ atb,scales="free_x", space='free_x') +
  annotate(geom="point",x=2000,y=0,size=0.01,colour="transparent") + 
  annotate(geom="point",x=2020,y=0,size=0.01,colour="transparent") + 
  scale_y_continuous(labels=scales::percent,breaks=c(0,0.05,.2,.4,.6),minor_breaks = NULL) +
  scale_x_continuous(expand=expansion(add=c(1.5,0)),breaks=seq(2000,2030,by=10)) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(3,1,2)],guide="none") +
  scale_colour_manual(values=RColorBrewer::brewer.pal(4,"Set1")[c(3,1,2)],guide="none") +
  coord_cartesian(ylim=c(0,0.55)) +
  labs(x="Year",y="Proportion")
ggarrange(g1,g2,nrow=2,labels=c("A","B"))
ggsave(file="figures/fig3_v4.pdf",height=6,width=8)



# Fig 4 ----


multi_step_pred = extract_mic(S_multistep_grasp_ceftriaxone)
tmp_threshold = data.frame(mic=0.25,
                           thres=0.05)
gA = ggplot() +
  geom_area(data=grasp_ceftriaxone_mic,aes(x=year,y=n/sample_size,group=micres,fill=mic),colour="black",size=0.3,alpha=.7) +
  scale_fill_gradient(low = RColorBrewer::brewer.pal(5,"Blues")[1],
                      high = RColorBrewer::brewer.pal(7,"Blues")[7],
                      trans="log2",
                      breaks=unique(grasp_ceftriaxone_mic$mic)) +
  scale_x_continuous(expand=c(0,0),limits=c(2004,2018)) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1)) +
  facet_grid(. ~ gender_or) +
  labs(x=NULL,y="Proportion",fill="MIC (mg/L)") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

gB = ggplot() +
  geom_ribbon(data=filter(single_step_pred,year>2018),aes(x=year,ymin=0,ymax=1),alpha=.1) +
  geom_ribbon(data=multi_step_pred$mic,aes(x=year,ymax=`97.5%`,ymin=`2.5%`,fill=mic),alpha=.8) + 
  geom_line(data=multi_step_pred$mic,aes(x=year,y=`50%`) ) +
  geom_point(data=grasp_ceftriaxone_mic,aes(x=year,y=n/sample_size,fill=mic),colour="black",shape=21,size=1.5) +
  scale_fill_gradient(low = RColorBrewer::brewer.pal(5,"Blues")[1],
                      high = RColorBrewer::brewer.pal(7,"Blues")[7],
                      trans="log2",
                      breaks=unique(multi_step_pred$mic$mic),guide=FALSE) +
  facet_grid(gender_or ~ mic) +
  scale_y_continuous(labels=scales::percent) +
  coord_cartesian(xlim=c(2004,2031),ylim=c(0,0.9)) +
  geom_hline(data=tmp_threshold,aes(yintercept=thres),linetype=2,size=.2) +
  labs(x="Year",y="Proportion",fill="MIC (mg/L)")

threspred = extract_thres(S_binary_grasp_ceftriaxone,"Ceftriaxone",2050) %>%
  mutate(model="Single-step")
threspred = filter(multi_step_pred$thres,threshold=="5%") %>%
  mutate(model="Multi-step") %>%
  bind_rows(threspred) %>%
  mutate(model=factor(model,levels=c("Single-step","Multi-step")))
  
gC = ggplot(threspred) +
  geom_line(aes(x=year,y=value,linetype=gender_or,colour=model)) +
  labs(x="Year",y="Probability",colour="Model type",linetype="Subpopulation") +
  coord_cartesian(xlim=c(2004,2030)) +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  scale_colour_manual(values=c("black","orange")) +
  theme(legend.position=c(0.26,.7),
        legend.margin=margin(1,2,2,2,"pt"),
        legend.background = element_blank(),
        legend.text=element_text(size=6),
        legend.key.height=unit(7,"pt"))
  
ggarrange(ggarrange(gC,gA,widths = c(1,2),labels=c("A","B")),gB,nrow=2,labels=c("","C"))
ggsave(file="figures/fig4.pdf",height=5.5,width=8.3)


