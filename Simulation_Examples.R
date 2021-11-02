
#########################################################################################################################
######Simulation examples for empirical power and type 1 error rate (as well as sample size and predicted power).########
#########################################################################################################################

source(".../Simulation_Functions_HTE.R")
set.seed(234)

#Example 1: 
#Continuous endpoint with delta=0.10 and cv=0.0 for heterogeneous treatment effect (HTE) 
start_time <- Sys.time()
con_delta010_cv00<-rbind(
  sim(mbar=20, rhox=0.10, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  
  sim(mbar=50, rhox=0.10, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  
  sim(mbar=100, rhox=0.10, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.01, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.05, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.10, cv=0.0, b4=0.10, tx="continuous", met="new", nsim=5000)
)
end_time <- Sys.time()
end_time - start_time
write.csv(con_delta010_cv00,"m1con_delta010_cv00.csv")


#Example 2: Binary endpoint with delta=0.45 and cv=0.6 for heterogeneous treatment effect (HTE) 
start_time <- Sys.time()
bin_delta045_cv60<-rbind(
  sim(mbar=20, rhox=0.10, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  
  sim(mbar=50, rhox=0.10, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  
  sim(mbar=100, rhox=0.10, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.01, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.05, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.10, cv=0.6, b4=0.45, tx="binary", met="new", nsim=5000)
)
end_time <- Sys.time()
end_time - start_time
write.csv(bin_delta045_cv60,"bin_delta045_cv60.csv")


source(".../Simulation_Functions_Mar.R")
set.seed(234)

#Example 3: Continuous endpoint with delta=0.25 and cv=0.3 for conditional overall effect (M1) 
start_time <- Sys.time()
con_delta025_cv30<-rbind(
  sim(mbar=20, rhox=0.10, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  
  sim(mbar=50, rhox=0.10, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  
  sim(mbar=100, rhox=0.10, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.01, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.05, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.10, cv=0.3, b4=0.25, tx="continuous", met="new", mar="M1", nsim=5000)
)
end_time <- Sys.time()
end_time - start_time
write.csv(con_delta025_cv30,"M1con_delta025_cv30.csv")




#Example 4: Continuous endpoint with delta=0.15 and cv=0.0 for marginal overall effect (M2) 
start_time <- Sys.time()
con_delta015_cv00<-rbind(
  sim(mbar=20, rhox=0.10, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.10, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.25, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=20, rhox=0.50, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  
  sim(mbar=50, rhox=0.10, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.10, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.25, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=50, rhox=0.50, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  
  sim(mbar=100, rhox=0.10, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.10, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.25, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.01, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.05, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000),
  sim(mbar=100, rhox=0.50, rhoy=0.10, cv=0.0, b4=0.15, tx="continuous", met="new", mar="M2", nsim=5000)
)
end_time <- Sys.time()
end_time - start_time
write.csv(con_delta015_cv00,"M2con_delta015_cv00.csv")



