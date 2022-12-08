# Loading libraries

library(tidyverse)
library(cowplot)
library(rstan)
options(mc.cores = parallel::detectCores())
library(loo)
library(readxl)
library(deSolve)
library(ggpubr)
library(MCMCpack)
library(RColorBrewer)

# tag
datetag = paste0("_",Sys.Date())

# ggplot2 theme
theme_custom = theme_light() + theme(text = element_text(size = 8),
                                     axis.text=element_text(size=7),
                                     strip.text=element_text(size=8))
theme_set(theme_custom)

# Functions
distsum = function(x,dist) {
  if(dist=="gamma") r = paste0(round(x[1]/x[2],2)," (",round(sqrt(x[1]/(x[2]^2)),2),")")
  if(dist=="beta") r = paste0(round(x[1]/(x[1]+x[2]),2)," (",round(sqrt(x[1]*x[2]/( (x[1] + x[2])^2*(x[1]+x[2]+1))),2),")")
  return(r)
}
get_gammapar = function(mu,sigma) {
  alpha = mu^2/sigma^2
  beta = mu/sigma^2
  return(c(alpha=alpha,beta=beta))
}
get_betapar = function(mu,sigma) {
  alpha = -(mu^3 + mu*sigma^2 - mu^2)/sigma^2
  beta = (mu^3 + (sigma^2 + 1)*mu - 2*mu^2 - sigma^2)/sigma^2
  return(c(alpha=alpha,beta=beta))
}
get_lnormpar = function(mu,sigma) {
  logmu = log(mu/sqrt(1+sigma^2/mu^2)); 
  logsigma = sqrt(log(1+sigma^2/mu^2));
  return(c(logmu=logmu,logsigma=logsigma))
}

check_gamma = function(shape,rate) {
  xx = rgamma(1000,shape,rate)
  plot(density(xx),main=paste0("~Gamma(",shape,",",rate,")"))
  cat(paste0("Distribution of rate: ",round(mean(xx),2)," [",round(quantile(xx,0.025),2),"-",round(quantile(xx,0.975),2),"] years^-1\n"))
  cat(paste0("Distribution of duration: ",round(mean(1/xx),2)," [",round(quantile(1/xx,0.025),2),"-",round(quantile(1/xx,0.975),2),"] years"))
}
lognormal_repar = function(raw,p) {
  Lp=c(NA,NA)
  Lp[1] = log(p[1]/sqrt(1+p[2]^2/p[1]^2))
  Lp[2] = sqrt(log(1+p[2]^2/p[1]^2))
  repar = exp(Lp[1] + raw * Lp[2])
  return(repar)
}
lognormal_getpar = function(p) {
  Lp=c(NA,NA)
  Lp[1] = log(p[1]/sqrt(1+p[2]^2/p[1]^2))
  Lp[2] = sqrt(log(1+p[2]^2/p[1]^2))
  return(Lp)
}
check_lognorm = function(m,s) {
  xx = rnorm(1000,0,1)
  xx = lognormal_repar(xx,c(m,s))
  plot(density(xx),main=paste0("~Lognormal with mean ",m," and sd ",s,""))
  cat(paste0("Distribution of rate: ",round(mean(xx),2)," [",round(quantile(xx,0.025),2),"-",round(quantile(xx,0.975),2),"] years^-1\n"))
  cat(paste0("Distribution of duration: ",round(mean(1/xx),2)," [",round(quantile(1/xx,0.025),2),"-",round(quantile(1/xx,0.975),2),"] years"))
}
cortocov = function(corr,sd) sd %*% t(sd) * corr

logit = function(x) log(x/(1-x))
inv.logit = function(x) exp(x)/(1+exp(x))

g_legend = function(a.gplot) { 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 
log_growth = function(t,tau,delay,slope) tau / (1+exp(-slope*(t+1-delay)))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

softmax = function(x) {
  return(exp(x)/sum(exp(x)))
}

log_softmax = function(x) {
  return(x-log(sum(exp(x))))
}

mic_wide = function(x,subgroup,start) {
  xx = filter(x,gender_or==subgroup) %>%
    dplyr::select(year,mic,n) %>%
    spread(mic,n) %>%
    filter(year>=start+1) %>%
    mutate(year_sim=year-start)
  xx[is.na(xx)] = 0
  return(xx)
}
range_to_normal = function(x,quant) {
  m = (x[1] + x[2])/2
  q = x
  s = mean((q[1]-m)/qnorm(quant[1]),(q[2]-m)/qnorm(quant[2]))
  return(c(m,s))
}
qsum = function(x) c(reached=1-mean(is.na(x)),mean=mean(x,na.rm=T),quantile(x,0.025,na.rm=T),quantile(x,0.975,na.rm=T))
qsum2 = function(x) c(`50%`=median(x),quantile(x,c(0.025,0.975)))
DIC = function(x) {
  require(loo)
  y = extract_log_lik(x)
  pointwise = apply(y,2,mean)
  deviance = -2*sum(pointwise)
  p_D = 1/2*var(pointwise)
  return(DIC=deviance +2*p_D)
}
atbswitch = function(x,t1,t2,delta)  (1/(1+exp(-delta*(x-t1))) + 1/(1+exp(delta*(x-t2)))) - 1
atbswitch4 = function(x,t1,t2,delta,iota1,iota2,xi)  (iota1/(1+exp(-delta*(x-t1+xi))) + (1-iota1)/(1+exp(-delta*(x-t1))) + iota2/(1+exp(delta*(x-t2-xi))) + (1-iota2)/(1+exp(delta*(x-t2)))) - 1

format_post = function(median,lb,ub) {
  paste0(sprintf("%.1f",median*100),
         "% (95%CrI: ",
         sprintf("%.1f",lb*100),
         "-",
         sprintf("%.1f",ub*100),
         "%)"
         )
}
format_post2 = function(median,lb,ub) {
  paste0(sprintf("%.2f",median),
         " (",
         sprintf("%.2f",lb),
         "-",
         sprintf("%.2f",ub),
         ")"
  )
}
format_post3 = function(median,lb,ub) {
  paste0(sprintf("%.3f",median),
         " (",
         sprintf("%.3f",lb),
         "-",
         sprintf("%.3f",ub),
         ")"
  )
}
format_post4 = function(median,lb,ub) {
  paste0(sprintf("%.4f",median),
         " (",
         sprintf("%.4f",lb),
         "-",
         sprintf("%.4f",ub),
         ")"
  )
}
format_post5 = function(median,lb,ub) {
  paste0(sprintf("%.5f",median),
         " (",
         sprintf("%.5f",lb),
         "-",
         sprintf("%.5f",ub),
         ")"
  )
}

# Plot functions

extract_pred = function(sim,atb,lim=NULL) {
  S = sim$data_list$S
  L = sim$data_list$L
  ts = sim$data_list$ts
  dat = sim$data_selected
  if(is.null(lim)) lim = max(ts)
  dd_res = summary(sim$samples,pars=c("pred_res_p"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),each=S),
           year=rep(ts,2),
           atb=atb) %>%
    left_join(dat,by=c("gender_or","year")) %>%
    filter(year<lim)
  return(dd_res)
}

extract_thres = function(sim,atb,lim=NULL) {
  S = sim$data_list$S
  L = sim$data_list$L
  ts = sim$data_list$ts
  dat = sim$data_selected
  if(is.null(lim)) lim = max(ts)
  thres = rstan::extract(sim$samples,pars=c("pred_res_p"))[[1]] 
  thres5 = t(apply(thres,2:3,function(x) mean(x>0.05))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=sim$data_list$ts,threshold="5%") %>%
    gather("gender_or","value",1:2)
  return(thres5)
}


plot_summary = function(sim,lim=NULL,save=FALSE,legend.pos=NULL,colmic="Accent") {
  S = sim$data_list$S
  L = sim$data_list$L
  ts = sim$data_list$ts
  dat = sim$data_selected
  colres = RColorBrewer::brewer.pal(5,colmic)[5]
  if(is.null(lim)) lim = max(ts)
  dd_res = summary(sim$samples,pars=c("pred_res_p"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),each=S),
           year=rep(ts,2)) %>%
    left_join(dat,by=c("gender_or","year")) %>%
    filter(year<lim)
  dd_p = summary(sim$samples,pars=c("p"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),each=S),
           year=rep(ts[1:S],2)) %>%
    left_join(dat,by=c("gender_or","year"))%>%
    filter(year<lim)
  dd_prev = summary(sim$samples,pars=c("y"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),each=S*3),
           year=rep(rep(ts,each=3),2),
           comp=rep(c("S","I_1","I_2"),S*2)) %>%
    filter(comp=="S") %>%
    filter(year<lim)
  thres = rstan::extract(sim$samples,pars=c("pred_res_p"))[[1]] 
  thres5 = t(apply(thres,2:3,function(x) mean(x>0.05))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=sim$data_list$ts,threshold="5%") %>%
    gather("gender_or","value",1:2)
  thres10 = t(apply(thres,2:3,function(x) mean(x>0.1))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=sim$data_list$ts,threshold="10%") %>%
    gather("gender_or","value",1:2)
  thres20 = t(apply(thres,2:3,function(x) mean(x>0.2))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=sim$data_list$ts,threshold="20%") %>%
    gather("gender_or","value",1:2)
  dd_thres = bind_rows(thres5,thres10,thres20) %>%
    mutate(threshold=factor(threshold,levels=c("5%","10%","20%"))) %>%
    filter(year<lim)
  g0 = ggplot() +
    geom_ribbon(data=filter(dd_res,!is.na(n)),aes(x=year,ymin=`2.5%`,ymax=`97.5%`),alpha=.5,fill=colres) +
    geom_line(data=filter(dd_res,!is.na(n)),aes(x=year,y=`50%`)) +
    geom_point(data=dat,aes(x=year,y=n/sample_size)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Resistance") +
    scale_y_continuous(labels=scales::percent)
  g1 = ggplot() +
    geom_ribbon(data=dd_res,aes(x=year,ymin=`2.5%`,ymax=`97.5%`),alpha=.5,fill=colres) +
    geom_line(data=dd_res,aes(x=year,y=`50%`)) +
    geom_point(data=dat,aes(x=year,y=n/sample_size)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Resistance") +
    scale_y_continuous(labels=scales::percent)
  g2 = ggplot() +
    geom_line(data=dd_p,aes(x=year,y=`50%`),colour=colres) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Treatment usage") +
    scale_y_continuous(labels=scales::percent,limits=c(0,1))
  g3 = ggplot() +
    geom_ribbon(data=dd_prev,aes(x=year,ymin=1-`2.5%`,ymax=1-`97.5%`),alpha=.5) +
    geom_line(data=dd_prev,aes(x=year,y=1-`50%`)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Prevalence") +
    scale_y_continuous(labels=scales::percent)
  g4 = ggplot() +
    geom_line(data=dd_thres,aes(x=year,y=value,colour=threshold)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Probability",colour="Threshold") +
    scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
    theme(legend.position=c(0.1,.7),
                  legend.margin=margin(1,2,2,2,"pt"),
                  legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
                  legend.text=element_text(size=6),
                  legend.key.height=unit(7,"pt"))
  if(!is.null(legend.pos)) g4 = g4 +theme(legend.position=legend.pos,
                                            legend.margin=margin(1,2,2,2,"pt"),
                                            legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
                                            legend.text=element_text(size=6),
                                            legend.key.height=unit(7,"pt"))
  if(!save) print(
    ggarrange(g0,
              ggarrange(g2,g1,g3,g4,ncol=2,nrow=2,labels=LETTERS[2:5]),
              ncol=1,labels="A",heights=c(1,2))
  )
  if(save) return(list(fitmic=g0,fitres=g1,treatment=g2,prevalence=g3,thresholds=g4))
}




plot_summary_prior = function(sim,lim=NULL,colmic="Accent") {
  S = sim$data_list$S
  L = sim$data_list$L
  ts = sim$data_list$ts
  dat = sim$data_selected
  colres = RColorBrewer::brewer.pal(5,colmic)[5]
  if(is.null(lim)) lim = max(ts)
  dd_res = rstan::extract(sim$samples,pars=c("pred_res_p")) %>%
    unlist() %>%
    as.data.frame() %>%
    mutate(it=rep(1:2000,S*2),
           gender_or=rep(rep(c("HMW","MSM"),each=2000),50),
           year=rep(ts,each=2000*2)) %>%
  filter(gender_or=="MSM")
  ggplot() +
    geom_line(data=dd_res,aes(x=year,y=.,group=it),alpha=.1) +
    labs(x="Year",y="Resistance") +
    scale_y_continuous(labels=scales::percent)
}


plot_summary2 = function(sim,lim=NULL,save=FALSE,legend.pos=NULL,colmic="Accent") {
  samples = sim$samples
  data_list = sim$data_list
  data_selected = sim$data_selected
  S = data_list$S
  L = data_list$L
  K = data_list$K
  ts = data_list$ts
  dat = data_selected  %>%
    arrange(gender_or,year,mic)
  datres = data_selected %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n),.groups="drop_last") %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  mic_classes = unique(dat$mic)
  res_classes = unique(filter(dat,res==1)$mic)
  which_res = which(mic_classes %in% res_classes)
  colres = RColorBrewer::brewer.pal(min(9,length(mic_classes)),colmic)[min(9,length(mic_classes))]
  if(is.null(lim)) lim = max(ts)
  dd_mic = summary(samples,pars=c("pred_mic_prop"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>% 
    mutate(gender_or=rep(c("HMW","MSM"),each=S*K),
           year=rep(rep(ts,each=K),2),
           mic=rep(mic_classes,S*2)) %>%
    left_join(dat,by=c("gender_or","year","mic")) %>%
    filter(year<=lim)
  pred_mic_n = rstan::extract(samples,pars="pred_mic_n")[[1]] 
  pred_res = apply(pred_mic_n[,,,which_res],c(1,2,3),sum) / apply(pred_mic_n,c(1,2,3),sum)
  pred_res = apply(pred_res,c(2,3),qsum2)
  dd_res = apply(pred_res,1,rbind) %>%
    as.data.frame() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),S),
           year=rep(ts,each=2)) %>%
    left_join(datres,by = c("gender_or", "year"))%>%
    filter(year<=lim)
  dd_p = summary(samples,pars=c("p"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),each=S),
           year=rep(ts[1:S],2)) %>%
    left_join(dat,by=c("gender_or","year")) %>%
    filter(year<=lim)
  dd_prev = summary(samples,pars=c("y"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>% 
    mutate(gender_or=rep(c("HMW","MSM"),each=S*(K+1)),
           year=rep(rep(ts,each=(K+1)),2),
           comp=rep(c("S",paste0("K_",1:K)),S*2)) %>%
    filter(comp=="S") %>%
    filter(year<=lim)
  thres = apply(rstan::extract(samples,pars=c("pred_mic_prop"))[[1]][,,,which_res],1:3,sum)
  thres5 = t(apply(thres,2:3,function(x) mean(x>0.05))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=data_list$ts,threshold="5%") %>%
    gather("gender_or","value",1:2)
  thres10 = t(apply(thres,2:3,function(x) mean(x>0.1))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=data_list$ts,threshold="10%") %>%
    gather("gender_or","value",1:2)
  thres20 = t(apply(thres,2:3,function(x) mean(x>0.2))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=data_list$ts,threshold="20%") %>%
    gather("gender_or","value",1:2)
  dd_thres = bind_rows(thres5,thres10,thres20) %>%
    mutate(threshold=factor(threshold,levels=c("5%","10%","20%"))) %>%
    filter(year<lim)
  g0 = ggplot() +
    geom_ribbon(data=filter(dd_mic,!is.na(n)),aes(x=year,ymax=`97.5%`,ymin=`2.5%`,fill=as.factor(mic)),alpha=.5) + 
    geom_point(data=dat,aes(x=year,y=n/sample_size),shape=20) +
    geom_line(data=filter(dd_mic,!is.na(n)),aes(x=year,y=`50%`) ) +
    scale_fill_brewer(type="seq",palette="Accent") +
    facet_grid(gender_or ~ mic) +
    labs(x="Year",y="Proportion") +
    theme(legend.position="none")
  g1 = ggplot() +
    geom_ribbon(data=dd_res,aes(x=year,ymin=`2.5%`,ymax=`97.5%`),fill=colres,alpha=.5) +
    geom_line(data=dd_res,aes(x=year,y=`50%`)) +
    geom_point(data=datres,aes(x=year,y=n/sample_size)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Resistance") +
    scale_y_continuous(labels=scales::percent)
  g2 = ggplot() +
    geom_line(data=dd_p,aes(x=year,y=`50%`),colour=colres) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Treatment usage") +
    scale_y_continuous(labels=scales::percent,limits=c(0,1))
  g3 = ggplot() +
    geom_ribbon(data=dd_prev,aes(x=year,ymin=1-`2.5%`,ymax=1-`97.5%`),alpha=.5) +
    geom_line(data=dd_prev,aes(x=year,y=1-`50%`)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Prevalence") +
    scale_y_continuous(labels=scales::percent)
  g4 = ggplot() +
    geom_line(data=dd_thres,aes(x=year,y=value,colour=threshold)) +
    facet_wrap(~gender_or) +
    labs(x="Year",y="Probability",colour="Threshold") +
    scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
    theme(legend.position=c(0.1,.7),
          legend.margin=margin(1,2,2,2,"pt"),
          legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
          legend.text=element_text(size=6),
          legend.key.height=unit(7,"pt"))
  if(!is.null(legend.pos)) g4 = g4 +theme(legend.position=legend.pos,
                                          legend.margin=margin(1,2,2,2,"pt"),
                                          legend.box.background = element_rect(fill="white",colour = "grey60",size=.4),
                                          legend.text=element_text(size=6),
                                          legend.key.height=unit(7,"pt"))
  if(!save) print(
    suppressWarnings(ggarrange(g0,
                               ggarrange(g2,g1,g3,g4,ncol=2,nrow=2,labels=LETTERS[2:5]),
                               ncol=1,labels="A",heights=c(1,1.5)))
  )
  if(save) return(list(fitmic=g0,fitres=g1,treatment=g2,prevalence=g3,thresholds=g4))
}



plot_summary_prior2 = function(sim,lim=NULL,colmic="Accent") {
  samples = sim$samples
  data_list = sim$data_list
  data_selected = sim$data_selected
  S = data_list$S
  L = data_list$L
  K = data_list$K
  ts = data_list$ts
  dat = data_selected  %>%
    arrange(gender_or,year,mic)
  datres = data_selected %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n),.groups="drop_last") %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  mic_classes = unique(dat$mic)
  res_classes = unique(filter(dat,res==1)$mic)
  which_res = which(mic_classes %in% res_classes)
  colres = RColorBrewer::brewer.pal(min(9,length(mic_classes)),colmic)[min(9,length(mic_classes))]
  if(is.null(lim)) lim = max(ts)
  dd_mic = summary(samples,pars=c("pred_mic_prop"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>% 
    mutate(gender_or=rep(c("HMW","MSM"),each=S*K),
           year=rep(rep(ts,each=K),2),
           mic=rep(mic_classes,S*2)) %>%
    left_join(dat,by=c("gender_or","year","mic")) %>%
    filter(gender_or=="MSM") %>%
    filter(year<=lim)
  pred_mic_n = rstan::extract(samples,pars="pred_mic_n")[[1]] 
  pred_res = apply(pred_mic_n[,,,which_res],c(1,2,3),sum) / apply(pred_mic_n,c(1,2,3),sum)
  dd_res = structure(pred_res,dim=400000) %>%
    as.data.frame() %>%
    mutate(it=rep(1:4000,S*2),
           gender_or=rep(rep(c("HMW","MSM"),each=4000),50),
           year=rep(ts,each=4000*2)) %>%
    filter(gender_or=="MSM")
  ggplot() +
    geom_line(data=dd_res,aes(x=year,y=.,group=it),alpha=.1) +
    labs(x="Year",y="Resistance") +
    scale_y_continuous(labels=scales::percent)
}


extract_mic = function(sim) {
  samples = sim$samples
  data_list = sim$data_list
  data_selected = sim$data_selected
  S = data_list$S
  L = data_list$L
  K = data_list$K
  ts = data_list$ts
  dat = data_selected  %>%
    arrange(gender_or,year,mic)
  datres = data_selected %>%
    group_by(year,gender_or,res) %>%
    summarise(n=sum(n)) %>%
    mutate(sample_size=sum(n)) %>%
    filter(res==1)
  mic_classes = unique(dat$mic)
  res_classes = unique(filter(dat,res==1)$mic)
  which_res = which(mic_classes %in% res_classes)
  dd_mic = summary(samples,pars=c("pred_mic_prop"))[[1]] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    tbl_df() %>% 
    mutate(gender_or=rep(c("HMW","MSM"),each=S*K),
           year=rep(rep(ts,each=K),2),
           mic=rep(mic_classes,S*2)) %>%
    left_join(dat,by=c("gender_or","year","mic")) 
  pred_mic_n = rstan::extract(samples,pars="pred_mic_n")[[1]] 
  pred_res = apply(pred_mic_n[,,,which_res],c(1,2,3),sum) / apply(pred_mic_n,c(1,2,3),sum)
  pred_res = apply(pred_res,c(2,3),qsum2)
  dd_res = apply(pred_res,1,rbind) %>%
    as.data.frame() %>%
    tbl_df() %>%
    mutate(gender_or=rep(c("HMW","MSM"),S),
           year=rep(ts,each=2)) %>%
    left_join(datres,by = c("gender_or", "year"))
  thres = apply(rstan::extract(samples,pars=c("pred_mic_prop"))[[1]][,,,which_res],1:3,sum)
  thres5 = t(apply(thres,2:3,function(x) mean(x>0.05))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=data_list$ts,threshold="5%") %>%
    gather("gender_or","value",1:2)
  thres10 = t(apply(thres,2:3,function(x) mean(x>0.1))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=data_list$ts,threshold="10%") %>%
    gather("gender_or","value",1:2)
  thres20 = t(apply(thres,2:3,function(x) mean(x>0.2))) %>%
    as.data.frame() %>%
    tbl_df() %>%
    rename(HMW=V1,MSM=V2) %>%
    mutate(year=data_list$ts,threshold="20%") %>%
    gather("gender_or","value",1:2)
  dd_thres = bind_rows(thres5,thres10,thres20) %>%
    mutate(threshold=factor(threshold,levels=c("5%","10%","20%")))
  return(list(mic=dd_mic,res=dd_res,thres=dd_thres))
}






post_summary = function(sim, format.table=TRUE) {
  pars = c("beta[1]","beta[2]","tau[1]","tau[2]","nu","mu","epsilon")
  parnames = c("$\\beta$ (HMW)","$\\beta$ (MSM)","$\\tau$ (HMW)","$\\tau$ (MSM)","$\\nu$","$\\mu$","$\\epsilon$")
  parcom = c("Transmission rate (HMW)","Transmission rate (MSM)","Treatment rate (HMW)","Treatment rate (MSM)",
             "Spontaneous recovery rate","Mutation rate","Reduced treatment efficacy")
  parunit = c("Per year","Per year","Per year","Per year","Per year","Probability","Probability")
  ss = summary(sim$samples,pars=c("beta", "mu", "nu", "epsilon", "tau"))[[1]] %>% 
    as.data.frame() %>% 
    rownames_to_column("pars") %>% 
    as_tibble() 
  tt = tibble(pars=pars,
              name=parnames,
              comment=parcom,
              unit=parunit) %>% 
    left_join(ss,by="pars") %>% 
    mutate(post=if_else(`50%`<.01,
                        format_post5(`50%`,`2.5%`,`97.5%`),
                        format_post2(`50%`,`2.5%`,`97.5%`)))
  if(format.table) {
    tt = tt %>% 
      dplyr::select(name,comment,post,unit) %>%
      knitr::kable("html",
                   col.names =  c("Parameter",
                                  "Interpretation",
                                  "Posterior median (95% credible interval)",
                                  "Unit"))
  }
  return(tt)
}

