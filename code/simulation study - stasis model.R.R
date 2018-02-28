#### Simulations testing for type I error in adequacy tests ###
library(adePEM)
library(paleoTS)
### 

replicates<-500
nr_datasets<-1000
nr_parameters<-4

##########
sequence.length<-10
omega.true<-0.01
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n10_omega0.01_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n10_omega0.01_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n10_omega0.01_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n10_omega0.01_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n10_omega0.01_theta1_trait.var0.05))


for (i in 1:length(out_n10_omega0.01_theta1_trait.var0.05[,1])){
  
  if (out_n10_omega0.01_theta1_trait.var0.05[i,1]>0.975 | out_n10_omega0.01_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n10_omega0.01_theta1_trait.var0.05[i,2]>0.975 | out_n10_omega0.01_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n10_omega0.01_theta1_trait.var0.05[i,3]>0.975 | out_n10_omega0.01_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n10_omega0.01_theta1_trait.var0.05[i,4]>0.975 | out_n10_omega0.01_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n10_omega0.01_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n10_omega0.01_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.01_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.01_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.01_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n10_omega0.01_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 


########
sequence.length<-20
omega.true<-0.01
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n20_omega0.01_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n20_omega0.01_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n20_omega0.01_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n20_omega0.01_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n20_omega0.01_theta1_trait.var0.05))


for (i in 1:length(out_n20_omega0.01_theta1_trait.var0.05[,1])){
  
  if (out_n20_omega0.01_theta1_trait.var0.05[i,1]>0.975 | out_n20_omega0.01_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n20_omega0.01_theta1_trait.var0.05[i,2]>0.975 | out_n20_omega0.01_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n20_omega0.01_theta1_trait.var0.05[i,3]>0.975 | out_n20_omega0.01_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n20_omega0.01_theta1_trait.var0.05[i,4]>0.975 | out_n20_omega0.01_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n20_omega0.01_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n20_omega0.01_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.01_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.01_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.01_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n20_omega0.01_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 


########
sequence.length<-40
omega.true<-0.01
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n40_omega0.01_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n40_omega0.01_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n40_omega0.01_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n40_omega0.01_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n40_omega0.01_theta1_trait.var0.05))


for (i in 1:length(out_n40_omega0.01_theta1_trait.var0.05[,1])){
  
  if (out_n40_omega0.01_theta1_trait.var0.05[i,1]>0.975 | out_n40_omega0.01_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n40_omega0.01_theta1_trait.var0.05[i,2]>0.975 | out_n40_omega0.01_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n40_omega0.01_theta1_trait.var0.05[i,3]>0.975 | out_n40_omega0.01_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n40_omega0.01_theta1_trait.var0.05[i,4]>0.975 | out_n40_omega0.01_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n40_omega0.01_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n40_omega0.01_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.01_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.01_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.01_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n40_omega0.01_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 


########
sequence.length<-80
omega.true<-0.01
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n80_omega0.01_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n80_omega0.01_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n80_omega0.01_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n80_omega0.01_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n80_omega0.01_theta1_trait.var0.05))


for (i in 1:length(out_n80_omega0.01_theta1_trait.var0.05[,1])){
  
  if (out_n80_omega0.01_theta1_trait.var0.05[i,1]>0.975 | out_n80_omega0.01_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n80_omega0.01_theta1_trait.var0.05[i,2]>0.975 | out_n80_omega0.01_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n80_omega0.01_theta1_trait.var0.05[i,3]>0.975 | out_n80_omega0.01_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n80_omega0.01_theta1_trait.var0.05[i,4]>0.975 | out_n80_omega0.01_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n80_omega0.01_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n80_omega0.01_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.01_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.01_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.01_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n80_omega0.01_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

 

########
sequence.length<-10
omega.true<-0.02
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n10_omega0.02_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n10_omega0.02_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n10_omega0.02_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n10_omega0.02_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n10_omega0.02_theta1_trait.var0.05))


for (i in 1:length(out_n10_omega0.02_theta1_trait.var0.05[,1])){
  
  if (out_n10_omega0.02_theta1_trait.var0.05[i,1]>0.975 | out_n10_omega0.02_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n10_omega0.02_theta1_trait.var0.05[i,2]>0.975 | out_n10_omega0.02_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n10_omega0.02_theta1_trait.var0.05[i,3]>0.975 | out_n10_omega0.02_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n10_omega0.02_theta1_trait.var0.05[i,4]>0.975 | out_n10_omega0.02_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n10_omega0.02_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n10_omega0.02_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.02_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.02_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.02_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n10_omega0.02_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets


save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-20
omega.true<-0.02
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n20_omega0.02_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n20_omega0.02_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n20_omega0.02_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n20_omega0.02_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n20_omega0.02_theta1_trait.var0.05))


for (i in 1:length(out_n20_omega0.02_theta1_trait.var0.05[,1])){
  
  if (out_n20_omega0.02_theta1_trait.var0.05[i,1]>0.975 | out_n20_omega0.02_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n20_omega0.02_theta1_trait.var0.05[i,2]>0.975 | out_n20_omega0.02_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n20_omega0.02_theta1_trait.var0.05[i,3]>0.975 | out_n20_omega0.02_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n20_omega0.02_theta1_trait.var0.05[i,4]>0.975 | out_n20_omega0.02_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n20_omega0.02_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n20_omega0.02_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.02_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.02_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.02_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n20_omega0.02_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-40
omega.true<-0.02
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n40_omega0.02_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n40_omega0.02_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n40_omega0.02_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n40_omega0.02_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n40_omega0.02_theta1_trait.var0.05))


for (i in 1:length(out_n40_omega0.02_theta1_trait.var0.05[,1])){
  
  if (out_n40_omega0.02_theta1_trait.var0.05[i,1]>0.975 | out_n40_omega0.02_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n40_omega0.02_theta1_trait.var0.05[i,2]>0.975 | out_n40_omega0.02_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n40_omega0.02_theta1_trait.var0.05[i,3]>0.975 | out_n40_omega0.02_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n40_omega0.02_theta1_trait.var0.05[i,4]>0.975 | out_n40_omega0.02_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n40_omega0.02_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n40_omega0.02_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.02_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.02_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.02_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n40_omega0.02_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-80
omega.true<-0.02
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
    y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
    y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
    y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
    y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
    y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n80_omega0.02_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n80_omega0.02_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n80_omega0.02_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n80_omega0.02_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n80_omega0.02_theta1_trait.var0.05))


for (i in 1:length(out_n80_omega0.02_theta1_trait.var0.05[,1])){
  
  if (out_n80_omega0.02_theta1_trait.var0.05[i,1]>0.975 | out_n80_omega0.02_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n80_omega0.02_theta1_trait.var0.05[i,2]>0.975 | out_n80_omega0.02_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n80_omega0.02_theta1_trait.var0.05[i,3]>0.975 | out_n80_omega0.02_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n80_omega0.02_theta1_trait.var0.05[i,4]>0.975 | out_n80_omega0.02_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n80_omega0.02_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n80_omega0.02_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.02_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.02_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.02_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n80_omega0.02_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-10
omega.true<-0.04
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
    y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
    y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
    y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
    y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
    y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n10_omega0.04_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n10_omega0.04_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n10_omega0.04_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n10_omega0.04_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n10_omega0.04_theta1_trait.var0.05))


for (i in 1:length(out_n10_omega0.04_theta1_trait.var0.05[,1])){
  
  if (out_n10_omega0.04_theta1_trait.var0.05[i,1]>0.975 | out_n10_omega0.04_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n10_omega0.04_theta1_trait.var0.05[i,2]>0.975 | out_n10_omega0.04_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n10_omega0.04_theta1_trait.var0.05[i,3]>0.975 | out_n10_omega0.04_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n10_omega0.04_theta1_trait.var0.05[i,4]>0.975 | out_n10_omega0.04_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n10_omega0.04_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n10_omega0.04_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.04_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.04_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.04_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n10_omega0.04_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-20
omega.true<-0.04
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n20_omega0.04_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n20_omega0.04_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n20_omega0.04_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n20_omega0.04_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n20_omega0.04_theta1_trait.var0.05))


for (i in 1:length(out_n20_omega0.04_theta1_trait.var0.05[,1])){
  
  if (out_n20_omega0.04_theta1_trait.var0.05[i,1]>0.975 | out_n20_omega0.04_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n20_omega0.04_theta1_trait.var0.05[i,2]>0.975 | out_n20_omega0.04_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n20_omega0.04_theta1_trait.var0.05[i,3]>0.975 | out_n20_omega0.04_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n20_omega0.04_theta1_trait.var0.05[i,4]>0.975 | out_n20_omega0.04_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n20_omega0.04_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n20_omega0.04_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.04_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.04_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.04_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n20_omega0.04_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-40
omega.true<-0.04
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
    y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
    y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
    y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
    y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
    y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n40_omega0.04_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n40_omega0.04_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n40_omega0.04_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n40_omega0.04_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n40_omega0.04_theta1_trait.var0.05))


for (i in 1:length(out_n40_omega0.04_theta1_trait.var0.05[,1])){
  
  if (out_n40_omega0.04_theta1_trait.var0.05[i,1]>0.975 | out_n40_omega0.04_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n40_omega0.04_theta1_trait.var0.05[i,2]>0.975 | out_n40_omega0.04_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n40_omega0.04_theta1_trait.var0.05[i,3]>0.975 | out_n40_omega0.04_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n40_omega0.04_theta1_trait.var0.05[i,4]>0.975 | out_n40_omega0.04_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n40_omega0.04_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n40_omega0.04_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.04_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.04_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.04_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n40_omega0.04_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-80
omega.true<-0.04
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
    y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
    y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
    y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
    y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
    y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n80_omega0.04_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n80_omega0.04_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n80_omega0.04_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n80_omega0.04_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n80_omega0.04_theta1_trait.var0.05))


for (i in 1:length(out_n80_omega0.04_theta1_trait.var0.05[,1])){
  
  if (out_n80_omega0.04_theta1_trait.var0.05[i,1]>0.975 | out_n80_omega0.04_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n80_omega0.04_theta1_trait.var0.05[i,2]>0.975 | out_n80_omega0.04_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n80_omega0.04_theta1_trait.var0.05[i,3]>0.975 | out_n80_omega0.04_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n80_omega0.04_theta1_trait.var0.05[i,4]>0.975 | out_n80_omega0.04_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n80_omega0.04_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n80_omega0.04_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.04_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.04_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.04_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n80_omega0.04_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-10
omega.true<-0.06
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
    y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
    y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
    y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
    y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
    y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n10_omega0.06_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n10_omega0.06_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n10_omega0.06_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n10_omega0.06_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n10_omega0.06_theta1_trait.var0.05))


for (i in 1:length(out_n10_omega0.06_theta1_trait.var0.05[,1])){
  
  if (out_n10_omega0.06_theta1_trait.var0.05[i,1]>0.975 | out_n10_omega0.06_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n10_omega0.06_theta1_trait.var0.05[i,2]>0.975 | out_n10_omega0.06_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n10_omega0.06_theta1_trait.var0.05[i,3]>0.975 | out_n10_omega0.06_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n10_omega0.06_theta1_trait.var0.05[i,4]>0.975 | out_n10_omega0.06_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n10_omega0.06_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n10_omega0.06_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.06_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.06_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n10_omega0.06_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n10_omega0.06_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-20
omega.true<-0.06
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n20_omega0.06_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n20_omega0.06_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n20_omega0.06_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n20_omega0.06_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n20_omega0.06_theta1_trait.var0.05))


for (i in 1:length(out_n20_omega0.06_theta1_trait.var0.05[,1])){
  
  if (out_n20_omega0.06_theta1_trait.var0.05[i,1]>0.975 | out_n20_omega0.06_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n20_omega0.06_theta1_trait.var0.05[i,2]>0.975 | out_n20_omega0.06_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n20_omega0.06_theta1_trait.var0.05[i,3]>0.975 | out_n20_omega0.06_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n20_omega0.06_theta1_trait.var0.05[i,4]>0.975 | out_n20_omega0.06_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n20_omega0.06_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n20_omega0.06_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.06_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.06_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n20_omega0.06_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n20_omega0.06_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-40
omega.true<-0.06
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
    y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
    y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
    y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
    y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
    y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n40_omega0.06_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n40_omega0.06_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n40_omega0.06_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n40_omega0.06_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n40_omega0.06_theta1_trait.var0.05))


for (i in 1:length(out_n40_omega0.06_theta1_trait.var0.05[,1])){
  
  if (out_n40_omega0.06_theta1_trait.var0.05[i,1]>0.975 | out_n40_omega0.06_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n40_omega0.06_theta1_trait.var0.05[i,2]>0.975 | out_n40_omega0.06_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n40_omega0.06_theta1_trait.var0.05[i,3]>0.975 | out_n40_omega0.06_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n40_omega0.06_theta1_trait.var0.05[i,4]>0.975 | out_n40_omega0.06_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n40_omega0.06_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n40_omega0.06_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.06_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.06_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n40_omega0.06_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n40_omega0.06_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets

save.image(file="simulations_Stasis.RData") 

#############
sequence.length<-80
omega.true<-0.06
theta.true<-1
trait.variance<-0.05


out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("slope.larger_than", "net.change_larger_than", "runs_test", "auto_corr")


for (i in 1:nr_datasets){
  
  
  y<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
  
  y_auto.corr<-adePEM::auto.corr(y$mm, model="stasis")
  y_runs.test<-adePEM::runs.test(y$mm, model="stasis", y$tt, theta = theta.true)
  y_slope.test<-adePEM::slope.test(y$mm,model="stasis", y$tt, theta=theta.true)
  y_net.change.test<-adePEM::net.change.test(y$mm, model="stasis")
  
  slope<-rep(NA, replicates)
  net_change<-rep(NA, replicates)
  runs_test<-rep(NA, replicates)
  auto_corr<-rep(NA, replicates)
  
  
  
  for (j in 1:replicates){
    
    x<-sim.Stasis(ns=sequence.length, theta=theta.true, omega=omega.true, vp=(theta.true*trait.variance), nn = sequence.length, tt = seq(1,sequence.length,1))
    
    slope[j]<-adePEM::slope.test(x$mm, model="stasis", y$tt, theta = theta.true)
    net_change[j]<-adePEM::net.change.test(x$mm, model="stasis") 
    runs_test[j]<-adePEM::runs.test(x$mm, model="stasis", y$tt, theta = theta.true)
    auto_corr[j]<-adePEM::auto.corr(x$mm, model="stasis")
    
  }
  
  print(i)
  
  
  out[i,1]<-length(slope[slope>y_slope.test])/replicates
  
  out[i,2]<-length(net_change[net_change>y_net.change.test])/replicates
  
  out[i,3]<-length(runs_test[runs_test>y_runs.test])/replicates
  
  out[i,4]<-length(auto_corr[auto_corr>y_auto.corr])/replicates
  
  
  
}


out_n80_omega0.06_theta1_trait.var0.05<-out


adequate_models.slope<-rep(NA, nrow(out_n80_omega0.06_theta1_trait.var0.05))
adequate_models.net.evolution<-rep(NA, nrow(out_n80_omega0.06_theta1_trait.var0.05))
adequate_models.runstest<-rep(NA, nrow(out_n80_omega0.06_theta1_trait.var0.05))
adequate_models.autocorr<-rep(NA, nrow(out_n80_omega0.06_theta1_trait.var0.05))


for (i in 1:length(out_n80_omega0.06_theta1_trait.var0.05[,1])){
  
  if (out_n80_omega0.06_theta1_trait.var0.05[i,1]>0.975 | out_n80_omega0.06_theta1_trait.var0.05[i,1]<0.025) adequate_models.slope[i]<-FALSE else adequate_models.slope[i]<-TRUE
  if (out_n80_omega0.06_theta1_trait.var0.05[i,2]>0.975 | out_n80_omega0.06_theta1_trait.var0.05[i,2]<0.025) adequate_models.net.evolution[i]<-FALSE else adequate_models.net.evolution[i]<-TRUE
  if (out_n80_omega0.06_theta1_trait.var0.05[i,3]>0.975 | out_n80_omega0.06_theta1_trait.var0.05[i,3]<0.025) adequate_models.runstest[i]<-FALSE else adequate_models.runstest[i]<-TRUE
  if (out_n80_omega0.06_theta1_trait.var0.05[i,4]>0.975 | out_n80_omega0.06_theta1_trait.var0.05[i,4]<0.025) adequate_models.autocorr[i]<-FALSE else adequate_models.autocorr[i]<-TRUE
  
}

# Combine vectors of model adequacy with data frame:
temp_data.frame<-as.data.frame(cbind(adequate_models.slope, adequate_models.net.evolution,adequate_models.runstest, adequate_models.autocorr))

results_out_n80_omega0.06_theta1_trait.var0.05<-rep(NA,nr_parameters+1)

# Number of traits not fulfilling the different criteria:

results_out_n80_omega0.06_theta1_trait.var0.05[2]<-(length(adequate_models.slope)- sum(adequate_models.slope, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.06_theta1_trait.var0.05[3]<-(length(adequate_models.net.evolution)- sum(adequate_models.net.evolution, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.06_theta1_trait.var0.05[4]<-(length(adequate_models.runstest)- sum(adequate_models.runstest, na.rm=TRUE))/nr_datasets

results_out_n80_omega0.06_theta1_trait.var0.05[5]<-(length(adequate_models.autocorr)- sum(adequate_models.autocorr, na.rm=TRUE))/nr_datasets


all_failed<-rep(NA,nr_datasets)
for (i in 1:length(adequate_models.autocorr)){
  
  if (adequate_models.slope[i]==TRUE && adequate_models.net.evolution[i] == TRUE && adequate_models.runstest[i]==TRUE && adequate_models.autocorr[i]==TRUE) all_failed[i]<-0 else all_failed[i]<-1
  
  
}

results_out_n80_omega0.06_theta1_trait.var0.05[1]<-sum(all_failed)/nr_datasets


#######################
#######################
#######################
#######################
#######################
#######################


par(mar=c(3.5,3.5,1,0.5)+1, mgp=c(2,1,0), mfrow=c(2,2))

plot(results_out_n80_omega0.06_theta1_trait.var0.05,  ylim=c(0,0.25), pch=2, col="black", xaxt = "n", ylab = "Type-I error",  xlab = "")
axis(1, at = 1:5, labels = paste(c("All", "Slope", "Net evol.", "Runs", "Auto"), sep = " "), cex.axis = 0.9)
abline(h=0.05, lty=2)
points(results_out_n40_omega0.06_theta1_trait.var0.05, pch=3, col="blue")
points(results_out_n20_omega0.06_theta1_trait.var0.05, pch=0, col="red")
points(results_out_n10_omega0.06_theta1_trait.var0.05, pch=5, col="orange")
legend("top", "omega = 0.06", bty="n") 

plot(results_out_n80_omega0.04_theta1_trait.var0.05,  ylim=c(0,0.25), pch=2, col="black", xaxt = "n", ylab = "Type-I error",  xlab = "")
axis(1, at = 1:5, labels = paste(c("All", "Slope", "Net evol.", "Runs", "Auto"), sep = " "), cex.axis = 0.9)
abline(h=0.05, lty=2)
points(results_out_n40_omega0.04_theta1_trait.var0.05, pch=3, col="blue")
points(results_out_n20_omega0.04_theta1_trait.var0.05, pch=0, col="red")
points(results_out_n10_omega0.04_theta1_trait.var0.05, pch=5, col="orange")
legend("top", "omega = 0.04", bty="n") 

plot(results_out_n80_omega0.02_theta1_trait.var0.05,  ylim=c(0,0.25), pch=2, col="black", xaxt = "n", ylab = "Type-I error",  xlab = "")
axis(1, at = 1:5, labels = paste(c("All", "Slope", "Net evol.", "Runs", "Auto"), sep = " "), cex.axis = 0.9)
abline(h=0.05, lty=2)
points(results_out_n40_omega0.02_theta1_trait.var0.05, pch=3, col="blue")
points(results_out_n20_omega0.02_theta1_trait.var0.05, pch=0, col="red")
points(results_out_n10_omega0.02_theta1_trait.var0.05, pch=5, col="orange")
legend("top", "omega = 0.02", bty="n") 

plot(results_out_n80_omega0.01_theta1_trait.var0.05,  ylim=c(0,0.25), pch=2, col="black", xaxt = "n", ylab = "Type-I error",  xlab = "")
axis(1, at = 1:5, labels = paste(c("All", "Slope", "Net evol.", "Runs", "Auto"), sep = " "), cex.axis = 0.9)
abline(h=0.05, lty=2)
points(results_out_n40_omega0.01_theta1_trait.var0.05, pch=3, col="blue")
points(results_out_n20_omega0.01_theta1_trait.var0.05, pch=0, col="red")
points(results_out_n10_omega0.01_theta1_trait.var0.05, pch=5, col="orange")
legend("top", "omega = 0.01", bty="n") 

