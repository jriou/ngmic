### Setup

source("setup.R")
source("data_management.R")


### Introduction and replacement of ATB and thus of selective pressure -----------------------------------------------------

# Double logistic function to model the progressive introduction of an antibiotic (p(t) = proportion of a particular 
# antibiotic) then its replacement by another one

atbswitch = function(x,t1,t2,delta) (1/(1+exp(-delta*(x-t1))) + 1/(1+exp(delta*(x-t2)))) - 1

ts = 2000:2018
t1 = 2005
t2 = 2014
delta = 3

ggplot(data.frame(x = seq(min(ts), max(ts), by=.1)), aes(x)) + 
  stat_function(fun = atbswitch,args=c(t1=t1,t2=t2,delta=delta),geom="line",colour="seagreen") +
  stat_function(fun = atbswitch,args=c(t1=t1,t2=t2,delta=delta),geom="area",fill="seagreen",alpha=.5) +
  geom_vline(xintercept=t1,linetype=2) +
  geom_vline(xintercept=t2,linetype=2) +
  scale_y_continuous(labels=scales::percent,expand=c(0,0),limits=c(0,1)) +
  annotate(geom="label",label="ATB of interest",x=(t1+t2)/2,y=.95) +
  annotate(geom="label",label="Previous ATB",x=(min(ts)+t1)/2,y=.95) +
  annotate(geom="label",label="Following ATB",x=(t2+max(ts))/2,y=.95) +
  labs(x="Time",y="p(t)")

# Double logistic function
atbswitch = function(x,t1,t2,delta)  (1/(1+exp(-delta*(x-t1))) + 1/(1+exp(delta*(x-t2)))) - 1
plot(2000:2018,atbswitch(x=2000:2018,t1=2002.5,t2=2010.5,delta=2),type="l")

# Triple logistic function
atbswitch_double_wave = function(x,t1,t2,delta,iota,nu,delta_prime)  (iota/(1+exp(-delta_prime*(x-t1+nu))) + (1-iota)/(1+exp(-delta*(x-t1))) + 1/(1+exp(delta*(x-t2)))) - 1
# simple growth if iota=0, nu=0, delta_prime=0 and t2=large number
plot(2000:2018,atbswitch_double_wave(x=2000:2018,t1=2010,t2=3000,delta=.7,iota=0,nu=0,delta_prime=0),type="l",ylim=c(0,1))
# double growth if t2=far future
plot(2000:2018,atbswitch_double_wave(x=2000:2018,t1=2010,t2=3000,delta=3,iota=.3,nu=5,delta_prime=0.8),type="l",ylim=c(0,1))
# simple degrowth if t1=far past
plot(2000:2018,atbswitch_double_wave(x=2000:2018,t1=0,t2=2015,delta=.7,iota=0,nu=0,delta_prime=0),type="l",ylim=c(0,1))
# growth and degrowth if iota=0, nu=0, delta_prime=0 
plot(2000:2018,atbswitch_double_wave(x=2000:2018,t1=2010,t2=2017,delta=2,iota=.1,nu=5,delta_prime=1),type="l",ylim=c(0,1))

# Quadruple logistic function
atbswitch4 = function(x,t1,t2,delta,iota1,iota2,xi)  (iota1/(1+exp(-delta*(x-t1+xi))) + (1-iota1)/(1+exp(-delta*(x-t1))) + iota2/(1+exp(delta*(x-t2-xi))) + (1-iota2)/(1+exp(delta*(x-t2)))) - 1
# simple growth if iota1=0, xi=0, delta_prime=0 and t2=large ximber
plot(2000:2018,atbswitch4(x=2000:2018,t1=2010,t2=3000,delta=.7,iota1=0,iota2=0,xi=0),type="l",ylim=c(0,1))
# double growth if t2=far future 
plot(2000:2018,atbswitch4(x=2000:2018,t1=2010,t2=3000,delta=3,iota1=.3,iota2=0,xi=5),type="l",ylim=c(0,1))
# simple degrowth if t1=far past
plot(2000:2018,atbswitch4(x=2000:2018,t1=0,t2=2010,delta=1.5,iota1=0,iota2=0,xi=0),type="l",ylim=c(0,1))
# double degrowth if t1=far past
plot(2000:2018,atbswitch4(x=2000:2018,t1=0,t2=2010,delta=1.5,iota1=0,iota2=.3,xi=5),type="l",ylim=c(0,1))
# growth and degrowth if iota=0, xi=0, delta_prime=0 
plot(2000:2018,atbswitch4(x=2000:2018,t1=2010,t2=2017,delta=2,iota1=0,iota2=0,xi=5),type="l",ylim=c(0,1))


pred_doublelog_coef = expand.grid(atb=unique(grasp_prescriptions_by_atb$atb),gender_or=c("HMW","MSM"),
                                  t1=0,t2=3000,delta=NA,iota1=0,iota2=0,xi=0)

par(mfrow=c(4,2))
for(i in 1:8) {
  dd = filter(grasp_prescriptions_by_atb,atb==pred_doublelog_coef[i,"atb"],gender_or==pred_doublelog_coef[i,"gender_or"])
  x = dd$year
  I = dd$p
  if(i %in% c(1,5)) { # simple growth and degrowth for azithromycin
    mod = nls(I ~ atbswitch4(x,t1,t2=2022,delta,iota1=0,iota2=0,xi=0), start=list(delta=1,t1=2000))
    pred_doublelog_coef[i,"t2"] = 2022
    pred_doublelog_coef[i,c("delta","t1")] = coefficients(mod)
    plot(predict(mod)~dd$year,type="l",main=i,ylim=c(0,1))
    points(dd$p ~ dd$year)
  } else if(i %in% c(2,6)) { # double degrowth for ciprofloxacin
    mod = nls(I ~ atbswitch4(x,t1=0,t2,delta,iota1=0,iota2,xi=7), start=list(delta=1,t2=2000,iota2=.01))
    pred_doublelog_coef[i,c("delta","t2","iota2")] = coefficients(mod)
    pred_doublelog_coef[i,"xi"] = 7
    plot(predict(mod)~dd$year,type="l",main=i)
    points(dd$p ~ dd$year)
  } else if(i %in% c(3,7)) { # simple growth and degrowth for cefixime
    mod = nls(I ~ atbswitch4(x,t1,t2,delta,iota1=0,iota2=0,xi=0), start=list(delta=1,t1=2000,t2=2010))
    pred_doublelog_coef[i,c("delta","t1","t2")] = coefficients(mod)
    plot(predict(mod)~dd$year,type="l",main=i)
    points(dd$p ~ dd$year)
  } else if(i %in% c(4,8)) { # double growth for ceftriaxone
    mod = nls(I ~ atbswitch4(x,t1=2010.5,t2=3000,delta,iota1,iota2=0,xi), start=list(iota1=.5,xi=5,delta=2))
    pred_doublelog_coef[i,c("iota1","xi","delta")] = coefficients(mod)
    pred_doublelog_coef[i,"t1"] = 2010.5
    plot(predict(mod)~dd$year,type="l",main=i,ylim=c(0,1))
    points(dd$p ~ dd$year)
  }
}
par(mfrow=c(1,1))

save(pred_doublelog_coef,file="data/pred_doublelog_coef.Rdata")


# Polynomial model fitted to data

pred_poly = function(x,c) {
  c[is.na(c)] = 0
  return(inv.logit(c[[1]]+c[[2]]*x+c[[3]]*x^2+c[[4]]*x^3+c[[5]]*x^4+c[[6]]*x^5))
}

pred_poly_coef = expand.grid(atb=unique(grasp_prescriptions_by_atb$atb),gender_or=c("HMW","MSM"),c0=NA,c1=NA,c2=NA,c3=NA,c4=NA,c5=NA)

par(mfrow=c(4,2))
for(i in 1:length(pred_poly_coef)) {
  dd = filter(grasp_prescriptions_by_atb,atb==pred_poly_coef[i,"atb"],gender_or==pred_poly_coef[i,"gender_or"])
  mm = lm(logit(p+1e-7) ~ poly(year,degree = 6,raw=TRUE),data=dd)
  plot(pred_poly(dd$year,mm$coefficients)~dd$year,type="l",main=i)
  points(dd$p ~ dd$year,ylim=c(0,1))
  pred_poly_coef[i,3:8] = mm$coefficients
}
par(mfrow=c(1,1))


# Splines
require(splines)

pred_splines = function(x,c) {
  c[is.na(c)] = 0
  return(inv.logit(c[[1]]+c[[2]]*x+c[[3]]*x^2+c[[4]]*x^3+c[[5]]*x^4+c[[6]]*x^5))
}

pred_splines_coef = expand.grid(atb=unique(grasp_prescriptions_by_atb$atb),gender_or=c("HMW","MSM"),intercept=NA,b1=NA,b2=NA,b3=NA,b4=NA,b5=NA,b6=NA)

par(mfrow=c(4,2))
for(i in 1:dim(pred_splines_coef)[1]) {
  dd = filter(grasp_prescriptions_by_atb,atb==pred_splines_coef[i,"atb"],gender_or==pred_splines_coef[i,"gender_or"])
  mm = lm(logit(p+1e-7) ~ ns(year,df=6),data=dd)
  plot(inv.logit(predict(mm))~dd$year,type="l",main=i)
  points(dd$p ~ dd$year,ylim=c(0,1))
  # pred_splines_coef[i,3:9] = mm$coefficients
}
par(mfrow=c(1,1))


plot(inv.logit(predict(mm)))
knots = attr(mm$model$`ns(year, df = 6)`, "knots")
betas = mm$coefficients

pred_ns = function(x,betas,knots) {
  return(inv.logit(betas[1] + 
    betas[2] * x + 
    betas[3]*(x-knots[1])^3 + 
    betas[4]*(x-knots[2])^3 +
    betas[5]*(x-knots[3])^3 +
    betas[6]*(x-knots[4])^3 +
    betas[7]*(x-knots[5])^3))
}
pred_ns(2000:2018,betas,knots)
inv.logit(predict(mm))

# Splines with rms
require(rms)

pred_splines_coef = expand.grid(atb=unique(grasp_prescriptions_by_atb$atb),gender_or=c("HMW","MSM"),intercept=NA,b1=NA,b2=NA,b3=NA,b4=NA,b5=NA,b6=NA)

par(mfrow=c(4,2))
for(i in 1:dim(pred_splines_coef)[1]) {
  dd = filter(grasp_prescriptions_by_atb,atb==pred_splines_coef[i,"atb"],gender_or==pred_splines_coef[i,"gender_or"])
  mm = lm(logit(p+1e-7) ~ ns(year,df=6),data=dd)
  plot(inv.logit(predict(mm))~dd$year,type="l",main=i)
  points(dd$p ~ dd$year,ylim=c(0,1))
  # pred_splines_coef[i,3:9] = mm$coefficients
}
par(mfrow=c(1,1))


inv.logit(predict(mm))



### Binary model ------------------------------------------------------------------------------------------------------------

### Chosen parameters

ts = 2000:2030
t1 = 2005
t2 = 2010
target_prevalence = 0.02
init_res = 0.001
beta = 3.5
mu = 0.2
nu = 1
epsilon = .5
tau = 1.5
delta = 4

# System of ODEs ------------------

binary_model = function(t,state,parameters){
  with(
    as.list(c(state,parameters)),
      {
        # time-dependent p(t)
        p = atbswitch(t,t1,t2,delta);
        # susceptible
        dS = - beta*S*(I+J) + ((1-p)*tau+p*tau*(1-mu)+nu)*I + ((1-p)*tau+p*tau*epsilon+nu)*J;
        # infected with wild-type strain
        dI = beta*S*I - ((1-p)*tau+p*tau*(1-mu)+nu)*I - p*tau*mu*I;
        # infected with resistant strain
        dJ = beta*S*J + p*tau*mu*I - ((1-p)*tau+p*tau*epsilon+nu)*J;
        return(list(c(dS,dI,dJ)))
      }
    )
}

# Solve
init = c(S=1-target_prevalence,I=target_prevalence*(1-init_res),J=target_prevalence*init_res)
parameters = c(beta=beta, mu=mu, nu=nu, epsilon=epsilon, tau=tau, delta=delta, t1=t1, t2=t2)
out = ode(y=init,times=ts,binary_model,parms=parameters)

tbl_df(as.data.frame(out)) %>%
  mutate(prev=1-S,prop_res=J/(I+J)) %>%
  dplyr::select(time,prev,prop_res) %>%
  gather("var","value",2:3) %>%
  ggplot(.) +
  geom_line(aes(x=time,y=value,colour=var)) +
  geom_vline(xintercept=t1,linetype=2) +
  geom_vline(xintercept=t2,linetype=2) +
  annotate("point",x=t1+2,y=0,colour="white") +
  scale_y_continuous(labels=scales::percent,expand=expand_scale(mult=c(0,.05))) +
  scale_colour_discrete(labels=c("Prevalence","Proportion of resistance")) +
  labs(x="Time",y="Proportion",colour=NULL) 

# Next generation matrix to find R0 to find beta so that prevalence stays stable at target_prevalence ----------

# if p=1
#>>>>>>>>>>>>>>>>>>>>>>>> In SageMath:
beta,mu,nu,epsilon,tau= var("beta,mu,nu,epsilon,tau")            # declarations
F = matrix([[beta,0],[0,beta]])                                  # infection matrix
V = matrix([[tau+nu,0],[-mu*tau,tau*epsilon+nu]])                # migration matrix
K = F*V.inverse()                                                # multiply F and V^-1 to get the NGM
r = K.eigenvalues()                                              # eigenvalues of the NGM
r                                                                # R0 = max(beta/(epsilon*tau + nu), beta/(nu + tau)) = beta/(nu + tau)
R0,target_prevalence=var("R0,target_prevalence")
R0=beta/(tau + nu)
solve(target_prevalence== 1-1/R0,beta)                           # beta =  (nu + tau)/(1-target_prevalence)


# Solve with beta constrained to maintain a stable prevalence at equilibrium if p=0 or p=1
init = c(S=1-target_prevalence,I=target_prevalence*(1-init_res),J=target_prevalence*init_res)
parameters = c(beta= (nu + tau)/(1-target_prevalence),
               mu=mu, nu=nu, epsilon=epsilon, tau=tau, delta=delta, t1=t1, t2=t2)
out = ode(y=init,times=2000:2050,binary_model,parms=parameters)

tbl_df(as.data.frame(out)) %>%
  mutate(prev=1-S,prop_res=J/(I+J)) %>%
  dplyr::select(time,prev,prop_res) %>%
  gather("var","value",2:3) %>%
  ggplot(.) +
  geom_line(aes(x=time,y=value,colour=var)) +
  geom_vline(xintercept=t1,linetype=2) +
  geom_vline(xintercept=t2,linetype=2) +
  annotate("point",x=t1+2,y=0,colour="white") +
  scale_y_continuous(labels=scales::percent,expand=expand_scale(mult=c(0,.05))) +
  scale_colour_discrete(labels=c("Prevalence","Proportion of resistance")) +
  labs(x="Time",y="Proportion",colour=NULL) 

# Solve with increasing prevalence

binary_model_incprev = function(t,state,parameters){
  with(
    as.list(c(state,parameters)),
    {
      # time-dependent p(t)
      p = atbswitch(t,t1,t2,delta);
      # susceptible
      dS = - beta*(1+zeta*(t-2000))*S*(I+J) + ((1-p)*tau+p*tau*(1-mu)+nu)*I + ((1-p)*tau+p*tau*epsilon+nu)*J;
      # infected with wild-type strain
      dI = beta*(1+zeta*(t-2000))*S*I - ((1-p)*tau+p*tau*(1-mu)+nu)*I - p*tau*mu*I;
      # infected with resistant strain
      dJ = beta*(1+zeta*(t-2000))*S*J + p*tau*mu*I - ((1-p)*tau+p*tau*epsilon+nu)*J;
      return(list(c(dS,dI,dJ)))
    }
  )
}

# Solve
zeta = 0.005
init = c(S=1-target_prevalence,I=target_prevalence*(1-init_res),J=target_prevalence*init_res)
parameters = c(beta= (nu + tau)/(1-target_prevalence), zeta=zeta,
               mu=mu, nu=nu, epsilon=epsilon, tau=tau, delta=delta, t1=t1, t2=t2)
out = ode(y=init,times=2000:2050,binary_model_incprev,parms=parameters)

tbl_df(as.data.frame(out)) %>%
  mutate(prev=1-S,prop_res=J/(I+J)) %>%
  dplyr::select(time,prev,prop_res) %>%
  gather("var","value",2:3) %>%
  ggplot(.) +
  geom_line(aes(x=time,y=value,colour=var)) +
  geom_vline(xintercept=t1,linetype=2) +
  geom_vline(xintercept=t2,linetype=2) +
  annotate("point",x=t1+2,y=0,colour="white") +
  scale_y_continuous(labels=scales::percent,expand=expand_scale(mult=c(0,.05))) +
  scale_colour_discrete(labels=c("Prevalence","Proportion of resistance")) +
  labs(x="Time",y="Proportion",colour=NULL) 






# Introduce actual treatment use with a functional form (forcing function)
pred_doublelog_coef
chosen_scenario = 6
parameters = list(beta= (nu + tau)/(1-target_prevalence),
               mu=mu, 
               nu=nu, 
               epsilon=epsilon, 
               tau=tau, 
               t1=pred_doublelog_coef[chosen_scenario,"t1"], 
               t2=pred_doublelog_coef[chosen_scenario,"t2"],
               delta=pred_doublelog_coef[chosen_scenario,"delta"],
               iota1=pred_doublelog_coef[chosen_scenario,"iota1"],
               iota2=pred_doublelog_coef[chosen_scenario,"iota2"],
               xi=pred_doublelog_coef[chosen_scenario,"xi"])

binary_model2 = function(t,state,parameters){
  with(
    as.list(c(state,parameters)),
    {
      # time-dependent p(t)
      p = atbswitch4(t,t1,t2,delta,iota1,iota2,xi);
      # susceptible
      dS = - beta*S*(I+J) + ((1-p)*tau+p*tau*(1-mu)+nu)*I + ((1-p)*tau+p*tau*epsilon+nu)*J;
      # infected with wild-type strain
      dI = beta*S*I - ((1-p)*tau+p*tau*(1-mu)+nu)*I - p*tau*mu*I;
      # infected with resistant strain
      dJ = beta*S*J + p*tau*mu*I - ((1-p)*tau+p*tau*epsilon+nu)*J;
      return(list(c(dS,dI,dJ)))
    }
  )
}

# Solve
init = c(S=1-target_prevalence,I=target_prevalence*(1-init_res),J=target_prevalence*init_res)
out = ode(y=init,times=ts,binary_model2,parms=parameters)

tbl_df(as.data.frame(out)) %>%
  mutate(prev=1-S,prop_res=J/(I+J),
         p=atbswitch4(time,parameters$t1,parameters$t2,parameters$delta,parameters$iota1,parameters$iota2,parameters$xi)) %>%
  dplyr::select(time,prev,prop_res,p) %>%
  gather("var","value",2:4) %>%
  ggplot(.) +
  geom_line(aes(x=time,y=value,colour=var)) +
  annotate("point",x=t1+2,y=0,colour="white") +
  scale_y_continuous(labels=scales::percent,expand=expand_scale(mult=c(0,.05))) +
  scale_colour_discrete(labels=c("ATB use","Prevalence","Proportion of resistance")) +
  labs(x="Time",y="Proportion",colour=NULL) 


### Implement in Stan ---------------------------------------------------
load("data/pred_doublelog_coef.Rdata")
source("fit_binary_model.R")

# prior predictive check (inference=0)
SIM_grasp_azithro = fit_binary_model(
  res_data=grasp_azithromycin_mic, 
  groups=c("HMW","MSM"), 
  period_simulated=2001:2050,
  target_prevalence =  matrix(c(range_to_normal(c(0.0016,0.0038),quant=c(0.025,0.975)),
                                range_to_normal(c(0.0119,0.0279),quant=c(0.025,0.975))),
                                ncol=2,nrow=2,byrow=TRUE), 
  treatment_use = pred_doublelog_coef[c(1,5),3:8], 
  inference=0
  )
check_hmc_diagnostics(SIM_grasp_azithro$samples)
print(SIM_grasp_azithro$samples,pars=c("mu","epsilon","prevalence","init_res","nu","beta","tau"),digits_summary = 4)
print(SIM_grasp_azithro$samples,pars=c("output_res"),digits_summary = 4)

plot_summary(SIM_grasp_azithro)



### multistep model ------------------------------------------------------------------------------------------------------------

### Chosen parameters

ts = 2000:2050
t1 = 2010
t2 = 2050
target_prevalence = 0.02
init_res = rdirichlet(1,c(1000,76,30,10,0,0,0))
mu = 0.02
nu = 0.68
epsilon = seq(1,0.7,length.out = 7)
tau = 1
delta = 1.6
K = 8
beta = (tau+nu)/(1-target_prevalence)

# System of ODEs ------------------

multistep_model = function(t,state,parameters){
  with(
    as.list(c(state,parameters)),
    {
      # time-dependent p(t)
      p = atbswitch(t,t1,t2,delta);
      # susceptible
      dS = - beta*S*(I1+I2+I3+I4+I5+I6+I7) + 
        p*tau*(1-mu)*(epsilon[1]*I1+epsilon[2]*I2+epsilon[3]*I3+epsilon[4]*I4+epsilon[5]*I5+epsilon[6]*I6) + 
        p*tau*epsilon[7]*I7 +
        (1-p)*tau*(I1+I2+I3+I4+I5+I6+I7) +
        nu*(I1+I2+I3+I4+I5+I6+I7);
      # K MIC classes
      dI1 = beta*S*I1 - p*tau*mu*I1 - p*tau*(1-mu)*epsilon[1]*I1 - (1-p)*tau*I1 - nu*I1;
      dI2 = beta*S*I2 + p*tau*mu*I1 - p*tau*mu*I2 - p*tau*(1-mu)*epsilon[2]*I2 - (1-p)*tau*I2 - nu*I2;
      dI3 = beta*S*I3 + p*tau*mu*I2 - p*tau*mu*I3 - p*tau*(1-mu)*epsilon[3]*I3 - (1-p)*tau*I3 - nu*I3;
      dI4 = beta*S*I4 + p*tau*mu*I3 - p*tau*mu*I4 - p*tau*(1-mu)*epsilon[4]*I4 - (1-p)*tau*I4 - nu*I4;
      dI5 = beta*S*I5 + p*tau*mu*I4 - p*tau*mu*I5 - p*tau*(1-mu)*epsilon[5]*I5 - (1-p)*tau*I5 - nu*I5;
      dI6 = beta*S*I6 + p*tau*mu*I5 - p*tau*mu*I6 - p*tau*(1-mu)*epsilon[6]*I6 - (1-p)*tau*I6 - nu*I6;
      dI7 = beta*S*I7 + p*tau*mu*I6 - p*tau*epsilon[7]*I7 - (1-p)*tau*I7 - nu*I7;
      return(list(c(dS,dI1,dI2,dI3,dI4,dI5,dI6,dI7)))
    }
  )
}

# Solve
init = c(S=1-target_prevalence,I1=target_prevalence*init_res[1],I2=target_prevalence*init_res[2],I3=target_prevalence*init_res[3],I4=target_prevalence*init_res[4],I5=target_prevalence*init_res[5],I6=target_prevalence*init_res[6],I7=target_prevalence*init_res[7])
parameters = c(beta=beta, mu=mu, nu=nu, epsilon=epsilon, tau=tau, delta=delta, t1=t1, t2=t2)
out = ode(y=init,times=ts,multistep_model,parms=parameters)

ggarrange(
tbl_df(as.data.frame(out)) %>%
  mutate(prev=1-S,res1=I1/prev,res2=I2/prev,res3=I3/prev,res4=I4/prev,res5=I5/prev,res6=I6/prev,res7=I7/prev) %>%
  dplyr::select(time,prev,res1,res2,res3,res4,res5,res6,res7) %>%
  gather("var","value",3:9) %>%
  ggplot(.) +
  geom_area(aes(x=time,y=value,group=var,fill=as.factor(var)),colour="black") +
  scale_fill_brewer(type="seq",palette="Accent",guide=FALSE) +
  labs(x="Time",y="Proportion",fill=NULL) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  geom_vline(xintercept=c(t1,t2),linetype=2)
,
tbl_df(as.data.frame(out)) %>%
  mutate(prev=1-S) %>%
  ggplot(.) +
  geom_line(aes(x=time,y=prev),colour="black") +
  labs(x="Time",y="Proportion",fill=NULL) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  annotate("point",x=t1+2,y=0,colour="white") +
  geom_vline(xintercept=c(t1,t2),linetype=2)
,ncol=1)


# Solve with beta constrained to maintain a stable prevalence at equilibrium if p=0 or p=1
init = c(S=1-target_prevalence,I1=target_prevalence*init_res[1],I2=target_prevalence*init_res[2],I3=target_prevalence*init_res[3],I4=target_prevalence*init_res[4],I5=target_prevalence*init_res[5],I6=target_prevalence*init_res[6],I7=target_prevalence*init_res[7])
parameters = c(beta=(nu + tau)/(1-target_prevalence), mu=mu, nu=nu, epsilon=epsilon, tau=tau, delta=delta, t1=t1, t2=t2)
out = ode(y=init,times=ts,multistep_model,parms=parameters)

ggarrange(
  tbl_df(as.data.frame(out)) %>%
    mutate(prev=1-S,res1=I1/prev,res2=I2/prev,res3=I3/prev,res4=I4/prev,res5=I5/prev,res6=I6/prev,res7=I7/prev) %>%
    dplyr::select(time,prev,res1,res2,res3,res4,res5,res6,res7) %>%
    gather("var","value",3:9) %>%
    ggplot(.) +
    geom_area(aes(x=time,y=value,group=var,fill=as.factor(var)),colour="black") +
    scale_fill_brewer(type="seq",palette="Accent",guide=FALSE) +
    labs(x="Time",y="Proportion",fill=NULL) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    geom_vline(xintercept=c(t1,t2),linetype=2)
  ,
  tbl_df(as.data.frame(out)) %>%
    mutate(prev=1-S) %>%
    ggplot(.) +
    geom_line(aes(x=time,y=prev),colour="black") +
    labs(x="Time",y="Proportion",fill=NULL) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    annotate("point",x=t1+2,y=0,colour="white") +
    geom_vline(xintercept=c(t1,t2),linetype=2)
  ,ncol=1)





# Introduce actual treatment use with a functional form (forcing function)
pred_doublelog_coef
chosen_scenario = 6
parameters = list(beta= (nu + tau)/(1-target_prevalence),
                  mu=mu, 
                  nu=nu, 
                  epsilon=epsilon, 
                  tau=tau, 
                  t1=pred_doublelog_coef[chosen_scenario,"t1"], 
                  t2=pred_doublelog_coef[chosen_scenario,"t2"],
                  delta=pred_doublelog_coef[chosen_scenario,"delta"],
                  iota1=pred_doublelog_coef[chosen_scenario,"iota1"],
                  iota2=pred_doublelog_coef[chosen_scenario,"iota2"],
                  xi=pred_doublelog_coef[chosen_scenario,"xi"])


multistep_model2 = function(t,state,parameters){
  with(
    as.list(c(state,parameters)),
    {
      # time-dependent p(t)
      p = atbswitch4(t,t1,t2,delta,iota1,iota2,xi);
      # susceptible
      dS = - beta*S*(I1+I2+I3+I4+I5+I6+I7) + 
        p*tau*(1-mu)*(epsilon[1]*I1+epsilon[2]*I2+epsilon[3]*I3+epsilon[4]*I4+epsilon[5]*I5+epsilon[6]*I6) + 
        p*tau*epsilon[7]*I7 +
        (1-p)*tau*(I1+I2+I3+I4+I5+I6+I7) +
        nu*(I1+I2+I3+I4+I5+I6+I7);
      # K MIC classes
      dI1 = beta*S*I1 - p*tau*mu*I1 - p*tau*(1-mu)*epsilon[1]*I1 - (1-p)*tau*I1 - nu*I1;
      dI2 = beta*S*I2 + p*tau*mu*I1 - p*tau*mu*I2 - p*tau*(1-mu)*epsilon[2]*I2 - (1-p)*tau*I2 - nu*I2;
      dI3 = beta*S*I3 + p*tau*mu*I2 - p*tau*mu*I3 - p*tau*(1-mu)*epsilon[3]*I3 - (1-p)*tau*I3 - nu*I3;
      dI4 = beta*S*I4 + p*tau*mu*I3 - p*tau*mu*I4 - p*tau*(1-mu)*epsilon[4]*I4 - (1-p)*tau*I4 - nu*I4;
      dI5 = beta*S*I5 + p*tau*mu*I4 - p*tau*mu*I5 - p*tau*(1-mu)*epsilon[5]*I5 - (1-p)*tau*I5 - nu*I5;
      dI6 = beta*S*I6 + p*tau*mu*I5 - p*tau*mu*I6 - p*tau*(1-mu)*epsilon[6]*I6 - (1-p)*tau*I6 - nu*I6;
      dI7 = beta*S*I7 + p*tau*mu*I6 - p*tau*epsilon[7]*I7 - (1-p)*tau*I7 - nu*I7;
      return(list(c(dS,dI1,dI2,dI3,dI4,dI5,dI6,dI7)))
    }
  )
}


# Solve
init = c(S=1-target_prevalence,I1=target_prevalence*init_res[1],I2=target_prevalence*init_res[2],I3=target_prevalence*init_res[3],I4=target_prevalence*init_res[4],I5=target_prevalence*init_res[5],I6=target_prevalence*init_res[6],I7=target_prevalence*init_res[7])
out = ode(y=init,times=ts,multistep_model2,parms=parameters)

ggarrange(
  tbl_df(as.data.frame(out)) %>%
    mutate(prev=1-S,res1=I1/prev,res2=I2/prev,res3=I3/prev,res4=I4/prev,res5=I5/prev,res6=I6/prev,res7=I7/prev) %>%
    dplyr::select(time,prev,res1,res2,res3,res4,res5,res6,res7) %>%
    gather("var","value",3:9) %>%
    ggplot(.) +
    geom_area(aes(x=time,y=value,group=var,fill=as.factor(var)),colour="black") +
    scale_fill_brewer(type="seq",palette="Accent",guide=FALSE) +
    labs(x="Time",y="Proportion",fill=NULL) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) 
  ,
  tbl_df(as.data.frame(out)) %>%
    mutate(prev=1-S,
           p=atbswitch4(time,parameters$t1,parameters$t2,parameters$delta,parameters$iota1,parameters$iota2,parameters$xi)) %>%
    ggplot(.) +
    geom_line(aes(x=time,y=prev),colour="black") +
    # geom_line(aes(x=time,y=p),colour="red",linetype=2) +
    labs(x="Time",y="Proportion",fill=NULL) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    annotate("point",x=t1+2,y=0,colour="white") 
  ,ncol=1)




