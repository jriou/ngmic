atbswitch = function(x,t1,t2,delta,iota1,iota2,xi)  (iota1/(1+exp(-delta*(x-t1+xi))) + (1-iota1)/(1+exp(-delta*(x-t1))) + iota2/(1+exp(delta*(x-t2-xi))) + (1-iota2)/(1+exp(delta*(x-t2)))) - 1


multistep_model = function(t,state,parameters){
  with(
    as.list(c(state,parameters)),
    {
      # time-dependent p(t)
      p = atbswitch(t,t1,t2,delta,iota1,iota2,xi)
      # epsilon for each class
      epsilon_k = array(8)
      for(k in 1:8) epsilon_k[k] = 1.0-(k-1.0)/(8-1.0)*(1.0-epsilon);
      # susceptible
      dS = - beta*S*(I1+I2+I3+I4+I5+I6+I7+I8) + 
        p*tau*(1-mu)*(epsilon_k[1]*I1+epsilon_k[2]*I2+epsilon_k[3]*I3+epsilon_k[4]*I4+epsilon_k[5]*I5+epsilon_k[6]*I6+epsilon_k[7]*I7) + 
        p*tau*epsilon_k[8]*I8 +
        (1-p)*tau*(I1+I2+I3+I4+I5+I6+I7+I8) +
        nu*(I1+I2+I3+I4+I5+I6+I7+I8);
      # K MIC classes
      dI1 = beta*S*I1 - p*tau*mu*I1 - p*tau*(1-mu)*epsilon_k[1]*I1 - (1-p)*tau*I1 - nu*I1;
      dI2 = beta*S*I2 + p*tau*mu*I1 - p*tau*mu*I2 - p*tau*(1-mu)*epsilon_k[2]*I2 - (1-p)*tau*I2 - nu*I2;
      dI3 = beta*S*I3 + p*tau*mu*I2 - p*tau*mu*I3 - p*tau*(1-mu)*epsilon_k[3]*I3 - (1-p)*tau*I3 - nu*I3;
      dI4 = beta*S*I4 + p*tau*mu*I3 - p*tau*mu*I4 - p*tau*(1-mu)*epsilon_k[4]*I4 - (1-p)*tau*I4 - nu*I4;
      dI5 = beta*S*I5 + p*tau*mu*I4 - p*tau*mu*I5 - p*tau*(1-mu)*epsilon_k[5]*I5 - (1-p)*tau*I5 - nu*I5;
      dI6 = beta*S*I6 + p*tau*mu*I5 - p*tau*mu*I6 - p*tau*(1-mu)*epsilon_k[6]*I6 - (1-p)*tau*I6 - nu*I6;
      dI7 = beta*S*I7 + p*tau*mu*I6 - p*tau*mu*I7 - p*tau*(1-mu)*epsilon_k[7]*I7 - (1-p)*tau*I7 - nu*I7;
      dI8 = beta*S*I8 + p*tau*mu*I7 - p*tau*epsilon_k[8]*I8 - (1-p)*tau*I8 - nu*I8;
      return(list(c(dS,dI1,dI2,dI3,dI4,dI5,dI6,dI7,dI8)))
    }
  )
}
