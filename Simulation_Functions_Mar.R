#source("/Users/tonytong/Dropbox/HTE Variable Sizes/Simulation/Simulation_Functions.R")

#########################################
# Marginal
# Sample size mar calculation (t test)
sample_size_cal_mar<- function(eff, rhox, rhoy, varx, vary, beta, alpha, m, p){
  varw=p*(1-p)
  kappa = vary*(1+(m-1)*rhoy)/varw/m
  
  n0=3
  n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  
  while(n0<n){
    n0=n0+1
    n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  }
  
  if (n0%%2==1) ss=n0+1 else ss=n0
  return(ss)
}

#marginal power (t test)
power_cal_mar<- function(eff, rhox, rhoy, varx, vary, k, alpha, m, p, cv){
  # k - number of clusters
  varw=p*(1-p)
  
  kappa<-vary*(1+(m-1)*rhoy)/varw/m
  power<-pt(sqrt(k*eff^2/kappa)-qt(1-alpha/2,df=k-2),df=k-2)
  return (power)
}

#New sample size calculation (t test)
sample_size_cal_new_mar<- function(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv){
  varw=p*(1-p)
  kappa<-vary*(1+(m-1)*rhoy)/varw/m*(1-cv^2*m*rhoy*(1-rhoy)/(1+(m-1)*rhoy)^2)^{-1}
  n0=3
  n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  
  while(n0<n){
    n0=n0+1
    n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  }
  
  if (n0%%2==1) ss=n0+1 else ss=n0
  return(ss)
}

#marginal power new (t test)
power_cal_new_mar<- function(eff, rhox, rhoy, varx, vary, k, alpha, m, p, cv){
  # k - number of clusters
  varw=p*(1-p)
  kappa<-vary*(1+(m-1)*rhoy)/varw/m*(1-cv^2*m*rhoy*(1-rhoy)/(1+(m-1)*rhoy)^2)^{-1}
  power<-pt(sqrt(k*eff^2/kappa)-qt(1-alpha/2,df=k-2),df=k-2)
  return (power)
}

# marginal Y with no X
sample_size_cal_Y<- function(eff, rhox, rhoy, varx, vary, beta, alpha, m, p){
  
  varw=p*(1-p)
  B=b3^2+(b4^2+2*b3*b4)*p
  omega=vary/(vary+B*varx)
  rhoym=omega*rhoy+(1-omega)*rhox
  
  varym=vary+B*varx
  
  kappa<-varym*(1+(m-1)*rhoym)/varw/m
  
  n0=3
  n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  
  while(n0<n){
    n0=n0+1
    n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  }
  
  if (n0%%2==1) ss=n0+1 else ss=n0
  return(ss)
}

#marginal power (t test)
power_cal_Y<- function(eff, rhox, rhoy, varx, vary, k, alpha, m, p, cv){
  # k - number of clusters
  varw=p*(1-p)
  B=b3^2+(b4^2+2*b3*b4)*p
  omega=vary/(vary+B*varx)
  rhoym=omega*rhoy+(1-omega)*rhox
  
  varym=vary+B*varx
  
  kappa<-varym*(1+(m-1)*rhoym)/varw/m
  power<-pt(sqrt(k*eff^2/kappa)-qt(1-alpha/2,df=k-2),df=k-2)
  return (power)
}

#New sample size calculation (t test)
sample_size_cal_new_Y<- function(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv){
  
  varw=p*(1-p)
  B=b3^2+(b4^2+2*b3*b4)*p
  omega=vary/(vary+B*varx)
  rhoym=omega*rhoy+(1-omega)*rhox
  
  varym=vary+B*varx
  
  kappa<-varym*(1+(m-1)*rhoym)/varw/m*(1-cv^2*m*rhoym*(1-rhoym)/(1+(m-1)*rhoym)^2)^{-1}
  n0=3
  n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  
  while(n0<n){
    n0=n0+1
    n=ceiling((qt(alpha/2,df=n0-2)+qt(beta,df=n0-2))^2*kappa/eff^2)
  }
  
  if (n0%%2==1) ss=n0+1 else ss=n0
  return(ss)
}

#marginal power new (t test)
power_cal_new_Y<- function(eff, rhox, rhoy, varx, vary, k, alpha, m, p, cv){
  # k - number of clusters
  varw=p*(1-p)
  B=b3^2+(b4^2+2*b3*b4)*p
  omega=vary/(vary+B*varx)
  rhoym=omega*rhoy+(1-omega)*rhox
  
  varym=vary+B*varx
  
  kappa<-varym*(1+(m-1)*rhoym)/varw/m*(1-cv^2*m*rhoym*(1-rhoym)/(1+(m-1)*rhoym)^2)^{-1}
  power<-pt(sqrt(k*eff^2/kappa)-qt(1-alpha/2,df=k-2),df=k-2)
  return (power)
}


#########################################################################################
####Data generating process
#ss -number of clusters for both arms
datagen<-function(mbar, rhox, rhoy, p, cv, b4, ss, tx="continuous", met="new",mar="M1"){
  #M1 is conditional on covariate M2 is marginal over covariate
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
  # with null effect
  #M1 is conditional on covariate M2 is marginal over covariate
  if (mar=="M1"){
    Y0 <- b1+b3*X+gamma+epsilon
  }
  if (mar=="M2"){
    Y0 <- b1+gamma+epsilon
  }
  df <- data.frame(tid, cluster, Z,X,Y,Y0)
  
  return(df)
}


############################
#Simulation
require("nlme")
sim<-function(mbar=10, rhox=0.1, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", mar="M1", nsim=1000){
  #Fixed parameters
  alpha = 0.05  #type-1 error
  beta = 0.20   #type-2 error
  p = 0.5       #treatment proportion
  pcov = 0.3    #binary X rate
  vary = 1      #variance of Y condition on X
 
  b2<<-0.25
  b3<<-0.1
  b4<<-b4
  
  #eff is now marginal effect
  eff = b2+0.5*b4
  m = mbar
  
  if (tx=="binary") varx = pcov*(1-pcov)   #variance of X
  if (tx=="continuous") varx= 1
  
  ##HTE: Calculate sample sizes
  #ss<-ceiling(sample_size_cal_new(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv))
  #tpower<-power_cal_new(eff, rhox, rhoy, varx, vary, k=ss, alpha, m, p,cv=cv)
  #no CV cluster size
  #ss0<-ceiling(sample_size_cal(eff, rhox, rhoy, varx, vary, beta, alpha, m, p))
  #tpower<-power_cal(eff, rhox, rhoy, varx, vary, k=ss, alpha, m, p)
  
  if (mar=="M1"){
    #marginal sample size and power calculation based on HTE (b4)
    ss<-sample_size_cal_new_mar(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv)
    tpower<-power_cal_new_mar(eff, rhox, rhoy, varx, vary, k=ss, alpha, m, p,cv=cv)
    #no CV cluster size
    ss0<-sample_size_cal_mar(eff, rhox, rhoy, varx, vary, beta, alpha, m, p)
    #mtpower0<-power_cal_mar
 
    #Store results
    output0=array(NA,dim=c(4,5,nsim))
    output1=array(NA,dim=c(4,5,nsim))
    erhoy<-evary<-rep(NA,nsim)
    
    for (s in 1:nsim){
      df<-datagen(mbar=mbar, rhox=rhox, rhoy=rhoy, cv=cv, b4=b4, ss=ss, tx=tx, met=met,mar="M1")  
      df$X<-df$X-mean(df$X)
      
      fit0 = try(lme(Y0 ~  Z*X, random = ~ 1|cluster, data = df),silent=T)
      if(class(fit0)=="try-error"){ s<- s-1; break}
      output0[,,s]<-summary(fit0)$tTable   
      
      fit1 = try(lme(Y  ~  Z*X, random = ~ 1|cluster, data = df),silent=T)
      if(class(fit1)=="try-error"){s<- s-1; break}
      output1[,,s]<-summary(fit1)$tTable
      
      var_comp<-as.numeric(VarCorr(fit1)[,"Variance"])
      erhoy[nsim]=var_comp[1]/sum(var_comp)
      evary[nsim]=var_comp[2]
    }
  }  
    
  if (mar=="M2"){
    #marginal sample size and power calculation based on marginalized X (b3&b4W)
    ss<- sample_size_cal_new_Y(eff, rhox, rhoy, varx, vary, beta, alpha, m, p, cv)
    tpower<-power_cal_new_Y(eff, rhox, rhoy, varx, vary, k=ss, alpha, m, p,cv=cv)
    #no CV cluster size
    ss0<-sample_size_cal_Y(eff, rhox, rhoy, varx, vary, beta, alpha, m, p)
    #mmtpower0<-power_cal_Y(eff, rhox, rhoy, varx, vary, k=mmss, alpha, m, p,cv=cv, b3=b3, b4=b4)
 
    #Store results
    output0=array(NA,dim=c(2,5,nsim))
    output1=array(NA,dim=c(2,5,nsim))
    
    B=b3^2+(b4^2+2*b3*b4)*p
    omega=vary/(vary+B*varx)
    rhoym=omega*rhoy+(1-omega)*rhox
    varym=vary+B*varx
    
    erhoy<-evary<-rep(NA,nsim)
    
    for (s in 1:nsim){
      df<-datagen(mbar=mbar, rhox=rhox, rhoy=rhoy, cv=cv, b4=b4, ss=ss, tx=tx, met=met,mar="M2")  
    
      fit0 = try(lme(Y0 ~  Z, random = ~ 1|cluster, data = df),silent=T)
      if(class(fit0)=="try-error"){ s<- s-1; break}
      output0[,,s]<-summary(fit0)$tTable  
    
      fit1 = try(lme(Y  ~  Z, random = ~ 1|cluster, data = df),silent=T)
      if(class(fit1)=="try-error"){s<- s-1; break}
      output1[,,s]<-summary(fit1)$tTable
      
      var_comp<-as.numeric(VarCorr(fit1)[,"Variance"])
      erhoy[nsim]=var_comp[1]/sum(var_comp)
      evary[nsim]=var_comp[2]
    }
    
    rhoy=rhoym
    vary=varym
  }  
  
  tp1<-round(mean(output0[2,5,]<0.05,na.rm=T) , 3)
  empower<-round(mean(output1[2,5,]<0.05,na.rm=T) , 3)
  erhoym<-round(mean(erhoy,na.rm=T),3)
  evarym<-round(mean(evary,na.rm=T),3)
  
  #srho<-round(sd(erhoyx,na.rm=T),3)
  #print(tp1)
  #print(empower)
  print(c(mbar, rhox, rhoy, cv))
  print(c(erhoym,rhoy))
  print(c(evarym,vary))
  return(data.frame("s0"=ss0, "sample_size"=ss, "rhoy"=rhoy, "erhoy"=erhoym, "vary"=vary, "evary"=evarym, "type1_error"=round(tp1*100,1), "true_power"=round(tpower*100,1), "em_power"=round(empower*100,1)))
}


######calculate sample size
sample_size_mar<-function(mar="M1"){
  
  b2 <<- 0.25
  b3 <<- 0.1
  
  delta_c<-b2+0.5*c(0.1,0.15,0.25)
  #delta_b<-c(0.25,0.35,0.45)
  rhoy_v<-c(0.01,0.05,0.10)
  rhox_v<-c(0.10,0.25,0.50)
  m<-c(20,50,100)
  cv<-c(0,0.3,0.6,0.9)

  
  S<-length(delta_c)*length(rhoy_v)*length(rhox_v)*length(m)*length(cv)
  NC<-matrix(NA,S,7)
  colnames(NC)<-c("delta","rhox","rhoy","m","cv","n0","n1")
  
  #if (type=="binary"){
  #  delta=delta_b
  #  varx=0.21
  #}
  #if (type=="continuous"){
    delta=delta_c
    varx=1
  #}  
  
  count=1
  for (i in 1:length(delta)){
    for (j in 1:length(rhox_v)){
      for (k in 1:length(rhoy_v)){
        for (l in 1:length(m)){
          for (q in 1:length(cv)){
            if (mar=="M1"){
              b4<<-delta[i]
              ss1<-sample_size_cal_new_mar(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5, cv=cv[q])
              ss0<-sample_size_cal_mar(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5)           
            }
            if (mar=="M2"){
              b4<<-delta[i]
              ss1<-sample_size_cal_new_Y(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5, cv=cv[q])
              ss0<-sample_size_cal_Y(eff=delta[i], rhox=rhox_v[j], rhoy=rhoy_v[k], varx=varx, vary=1, beta=0.2, alpha=0.05, m=m[l], p=0.5)           
            }
            
            NC[count,]<-c(delta[i],rhox_v[j],rhoy_v[k],m[l],cv[q],ss0,ss1)
            count=count+1
          }
        }
      }
    }
  }
  
  return(NC)
}


