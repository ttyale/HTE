###Yang/Li formula
#Sample size calculation
sample_size_cal<- function(eff, rhox, rhoy, varx, vary, beta, alpha, m, p ){
  kappa = 1/(1-rhoy)*(1 + (m-2)*rhoy - (m-1)*rhox*rhoy)/(1+(m-1)*rhoy)
  n = (qnorm(1-alpha/2) + qnorm(1-beta))^2*vary/(eff^2*m*p*(1-p)*kappa*varx)
  return (n)
}
#sample_size_cal(eff, rhox, rhoy, varx, vary, beta, alpha, m,p)

#power calculation
power_cal<- function(eff, rhox, rhoy, varx, vary, k, alpha, m, p ){
  # k - number of clusters
  kappa = 1/(1-rhoy)*(1 + (m-2)*rhoy - (m-1)*rhox*rhoy)/(1+(m-1)*rhoy)
  power =pnorm( sqrt(k*eff^2*m*p*(1-p)*kappa*varx/vary)-qnorm(1-alpha/2))
  return (power)
}
#power_cal(eff, rhox, rhoy, varx, vary, k=50, alpha, m, p)

####Tong/Li formula
#New sample size calculation
sample_size_cal_new<- function(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv){
  varw=p*(1-p)
  kappa<-vary*(1-rhoy)*(1+(m-1)*rhoy)^3/(varw*varx*m*((1+(m-2)*rhoy-(m-1)*rhox*rhoy)*(1+(m-1)*rhoy)^2+cv^2*m*rhoy*(1-rhoy)*(rhoy-rhox)))
  n = (qnorm(1-alpha/2) + qnorm(1-beta))^2*kappa/(eff^2)
  return (n)
}
#sample_size_cal_new(eff, rhox, rhoy, varx, vary, beta, alpha, m,p,cv=0)

#New power calculation
power_cal_new<- function(eff, rhox, rhoy, varx, vary, k, alpha, m, p, cv){
  # k - number of clusters
  varw=p*(1-p)
  kappa<-vary*(1-rhoy)*(1+(m-1)*rhoy)^3/(varw*varx*m*((1+(m-2)*rhoy-(m-1)*rhox*rhoy)*(1+(m-1)*rhoy)^2+cv^2*m*rhoy*(1-rhoy)*(rhoy-rhox)))
  power<-pnorm( sqrt(k*eff^2/kappa)-qnorm(1-alpha/2))
  return (power)
}
#power_cal_new(eff, rhox, rhoy, varx, vary, k=50, alpha, m, p, cv=0)


#########################################################################################
####Data generating process
#ss -number of clusters for both arms
datagen<-function(mbar, rhox, rhoy, p, cv, b4, ss, tx="continuous", met="new"){
  #Fixed parameters
  alpha = 0.05  #type-1 error
  beta = 0.20   #type-2 errpr
  p=0.5         #treatment proportion
  pcov=0.3      #binary X rate
  vary = 1      #variance of Y condition on X
 
  b1=0          #intercept
  b2=0.25       #treatment effect
  b3=0.1        #covariate effect
  eff=b4        #interaction effect to be powered
 
  if (tx=="binary") varx = pcov*(1-pcov)   #variance of X
  if (tx=="continuous") varx= 1
  
  #draw cluster size based on cv
  if (cv==0){
    cl_size=rep(mbar,ss)
  }else{
    #generate cl size from gamma
    a0=1/cv^2
    b0=1/mbar/cv^2
    cl_size=round(rgamma(n=ss,shape=a0,rate=b0))
    
    for (k in length(cl_size)){
      if (cl_size[k]<2) cl_size[k]=2
    }
  }
  
  #total sample size
  n <- sum(cl_size)
  cluster = rep(1:ss, cl_size)
  
  # Generate binary covariate with ICC from beta binomial distribution  
  if (tx=="binary"){
    # Generate from beta binomial
    ab_sum <- 1/rhox - 1
    a = pcov *ab_sum
    b = ab_sum - a
    
    p <- rbeta(ss, a, b)
    X <-c()
    for (i in 1:ss){
      X<-c(X,rbinom(cl_size[i],1,p[i]))
    }
    
  }
  if (tx=="continuous"){
    # Generate from multilevel normal
    sigmaxc <- sqrt(rhox*varx) #cluster level sd
    sigmaxe <- sqrt((1-rhox)*varx) # individual level sd
    
    miu<-rep(rnorm(ss, 0, sigmaxc), cl_size)
    X <-rnorm(n, 0, sigmaxe)+miu+1/2
  }
 
  #randomize treatment assignment
  trt <- rep(0, ss)
  t <- sample(1:ss, size=ss/2)
  trt[t] <- 1
  
  #individual treatment variable Z
  Z<- rep(trt, cl_size)
  
  #generate id
  tid <- seq(1:n)

  #beta <- effect_size_cal(k, rhox, rhoy , varx, vary= sigma, beta= 0.2, alpha = 0.05, m, p = 0.5 )
  sigmac <- sqrt(rhoy*vary)  #cluster level sd; sigma= vary
  sigmae <- sqrt((1-rhoy)*vary)  # individual level sd
  
  #individual random effect
  epsilon <- rnorm(n, 0, sigmae)
  
  #cluster random effect
  g <- rnorm(ss, 0, sigmac)
  gamma <- rep(g, cl_size)
  
  # with positive interaction  
  Y <- b1+b2*Z+b3*X+b4*X*Z+gamma+epsilon
  # with null interaction
  Y0 <- b1+b2*Z+b3*X+gamma+epsilon
  
  df <- data.frame(tid, cluster, Z,X,Y,Y0)

  return(df)
}



############################
#Simulation
require("nlme")
sim<-function(mbar=10, rhox=0.1, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new",nsim=1000){
  #Fixed parameters
  alpha = 0.05  #type-1 error
  beta = 0.20   #type-2 error
  p = 0.5       #treatment proportion
  pcov = 0.3    #binary X rate
  vary = 1      #variance of Y condition on X
  
  eff = b4
  m = mbar
  
  if (tx=="binary") varx = pcov*(1-pcov)   #variance of X
  if (tx=="continuous") varx= 1
  
  #Calculate sample sizes
    ss<-ceiling(sample_size_cal_new(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv))
    tpower<-power_cal_new(eff, rhox, rhoy, varx, vary, k=ss, alpha, m, p,cv=cv)
    
    #no CV sample size
    ss0<-ceiling(sample_size_cal(eff, rhox, rhoy, varx, vary, beta, alpha, m, p))
    #tpower<-power_cal(eff, rhox, rhoy, varx, vary, k=ss, alpha, m, p)
  
  if (ss %% 2 ==1) ss=ss+1
  if (ss0 %% 2 ==1) ss0=ss0+1
  #print(ss)
  #print(tpower)
  
  #Store results
  output0=array(NA,dim=c(4,5,nsim))
  output1=array(NA,dim=c(4,5,nsim))
  
  for (s in 1:nsim){
    #Generate data
    df<-datagen(mbar=mbar, rhox=rhox, rhoy=rhoy, cv=cv, b4=b4, ss=ss, tx=tx, met=met)
    
    #Fit model
    fit0 = try(lme(Y0 ~  Z*X, random = ~ 1|cluster, data = df),silent=T)
    if(class(fit0)=="try-error"){ s<- s-1; break}
    output0[,,s]<-summary(fit0)$tTable
    
    fit1 = try(lme(Y  ~  Z*X, random = ~ 1|cluster, data = df),silent=T)
    if(class(fit1)=="try-error"){s<- s-1; break}
    output1[,,s]<-summary(fit1)$tTable
  }  
  tp1<-round(mean(output0[4,5,]<0.05,na.rm=T) , 3)
  empower<-round(mean(output1[4,5,]<0.05,na.rm=T) , 3)
  
  #print(tp1)
  #print(empower)
  print(c(mbar, rhox, rhoy, cv))
  return(data.frame("s0"=ss0, "sample_size"=ss, "type1_error"=round(tp1*100,1), "true_power"=round(tpower*100,1), "em_power"=round(empower*100,1)))
}

######calculate sample size
sample_size<-function(type="binary"){
  
  delta_c<-c(0.1,0.15,0.25)
  delta_b<-c(0.25,0.35,0.45)
  rhoy_v<-c(0.01,0.05,0.10)
  rhox_v<-c(0.10,0.25,0.50)
  m<-c(20,50,100)
  cv<-c(0,0.3,0.6,0.9)

  S<-length(delta_c)*length(rhoy_v)*length(rhox_v)*length(m)*length(cv)
  NC<-matrix(NA,S,7)
  colnames(NC)<-c("delta","rhox","rhoy","m","cv","n0","n1")
  
  if (type=="binary"){
    delta=delta_b
    varx=0.21
  }
  if (type=="continuous"){
    delta=delta_c
    varx=1
  }  
  
  count=1
  for (i in 1:length(delta)){
    for (j in 1:length(rhox_v)){
      for (k in 1:length(rhoy_v)){
        for (l in 1:length(m)){
          for (q in 1:length(cv)){
            ss1<-ceiling(sample_size_cal_new(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5, cv=cv[q]))
            ss0<-ceiling(sample_size_cal(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5))
            
            #ss1<-sample_size_cal_new(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5, cv=cv[q])
            #ss0<-sample_size_cal(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5)
            if (ss0 %% 2 ==1) ss0=ss0+1
            if (ss1 %% 2 ==1) ss1=ss1+1
            
            NC[count,]<-c(delta[i],rhox_v[j],rhoy_v[k],m[l],cv[q],ss0,ss1)
            count=count+1
          }
        }
      }
    }
  }
  
  return(NC)
}
  






