library(adePEM)
library(paleoTS)

#####
##### Calculate model statistics on all time series in uploaded folders #####
#####

# set directory to folder with fossil time series. 
setwd("~/fossil time series")

temp = list.files(pattern="*.txt")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

nr_datasets<-length(temp)
nr_parameters<-66
replicates<-1000

out<-matrix(nrow=nr_datasets, ncol=nr_parameters, data=NA)
colnames(out)<-c("AICc.Directional.change","AICc.Random.walk", "AICc.Stasis", 
                 "Akaike.wt.Directional.change","Akaike.wt.Random.walk","Akaike.wt.Stasis",
                 "theta","omega", 
                 "auto.corr_stasis_est", "auto.corr_stasis_min","auto.corr_stasis_max","auto.corr_stasis_p", "auto.corr_stasis_result",
                 "runs.test_stasis_est", "runs.test_stasis_min","runs.test_stasis_max","runs.test_stasis_p", "runs.test_stasis_result",
                 "slope.test_stasis_est", "slope.test_stasis_min","slope.test_stasis_max","slope.test_stasis_p", "slope.test_stasis_result",
                 "net.change.test_stasis_est", "net.change.test_stasis_min","net.change.test_stasis_max","net.change.test_stasis_p", "net.change.test_stasis_result",
                 "vstep_RW",
                 "auto.corr_RW_est", "auto.corr_RW_min","auto.corr_RW_max","auto.corr_RW_p", "auto.corr_RW_result",
                 "runs.test_RW_est", "runs.test_RW_min","runs.test_RW_max","runs.test_RW_p", "runs.test_RW_result",
                 "slope.test_RW_est", "slope.test_RW_min","slope.test_RW_max","slope.test_RW_p", "slope.test_RW_result",
                 "mstep_trend", "vstep_trend",
                 "auto.corr_trend_est", "auto.corr_trend_min","auto.corr_trend_max","auto.corr_trend_p", "auto.corr_trend_result",
                 "runs.test_trend_est", "runs.test_trend_min","runs.test_trend_max","runs.test_trend_p", "runs.test_trend_result",
                 "slope.test_trend_est", "slope.test_trend_min","slope.test_trend_max","slope.test_trend_p", "slope.test_trend_result",
                 "interval.length.MY", "generation.time", "taxa", "dimensions", "species")

out<-as.data.frame(out)

for (i in 1:nr_datasets){
  indata<-read.table(temp[i], header=T)
  y<-as.paleoTS(indata$mm, indata$vv, indata$N, indata$age.in.MY)
  
  if (indata$log[1] =="no")   y<-ln.paleoTS(y);
  if (indata$dimension[1] =="area") y$mm<-y$mm/2; y$vv<-y$vv/4;  
  
  
  out[i,1:3]<-fit3models(y, method = "Joint")[1:3,3] 
  out[i,4:6]<-fit3models(y, method = "Joint")[1:3,4]
  
  theta<-opt.joint.Stasis(y)$par[1]
  omega<-opt.joint.Stasis(y)$par[2]
  out[i,7]<-theta
  out[i,8]<-omega
  
  # adequacy test of stasis
  temp_Stasis<-fit4adequacy.stasis(y, nrep=replicates, plot = FALSE)
  out[i,9:13]<-as.matrix(temp_Stasis$summary)[1,1:5]
  out[i,14:18]<-as.matrix(temp_Stasis$summary)[2,1:5]
  out[i,19:23]<-as.matrix(temp_Stasis$summary)[3,1:5]
  out[i,24:28]<-as.matrix(temp_Stasis$summary)[4,1:5]
  
  
  vstep_RW<-opt.joint.URW(y)$par[2]
  out[i,29]<-vstep_RW
  
  # adequacy test of RW
  temp_RW<-fit3adequacy.RW(y, nrep=replicates, plot = FALSE)
  out[i,30:34]<-as.matrix(temp_RW$summary)[1,1:5]
  out[i,35:39]<-as.matrix(temp_RW$summary)[2,1:5]
  out[i,40:44]<-as.matrix(temp_RW$summary)[3,1:5]
  
  mstep_trend<-opt.joint.GRW(y)$par[2]
  vstep_trend<-opt.joint.GRW(y)$par[3]
  out[i,45]<-mstep_trend
  out[i,46]<-vstep_trend
  
  # adequacy tests of trend
  temp_trend<-fit3adequacy.trend(y,nrep=replicates, plot = FALSE)
  out[i,47:51]<-as.matrix(temp_trend$summary)[1,1:5]
  out[i,52:56]<-as.matrix(temp_trend$summary)[2,1:5]
  out[i,57:61]<-as.matrix(temp_trend$summary)[3,1:5]
  
  
  out[i,62]<-indata$interval.length.MY[1]
  out[i,63]<-indata$generation.time[1]
  out[i,64]<-levels(droplevels(indata$taxa[1]))
  out[i,65]<-levels(droplevels(indata$dimension[1]))
  out[i,66]<-levels(droplevels(indata$species[1]))
  
  
  
  print(i)
  
  
}


############################
############################
############################




# Pick out the models that show a certain fit to the different modes:

mode<-rep(NA, length(out[,1]))

for (i in 1:length(mode)){
  
  if ((max(c(out[i,4], out[i,5],out[i,6])) == out[i,4]))  mode[i]<-"Directional"
  if ((max(c(out[i,4], out[i,5],out[i,6])) == out[i,5]))  mode[i]<-"Random walk"
  if ((max(c(out[i,4], out[i,5],out[i,6])) == out[i,6])) mode[i]<-"Stasis"
}



#Frequency of the different modes:
table(mode)

# Combine mode and data matrix:
indata_3_modes<-cbind(out,mode)
#write.csv(indata_stasis, "out_file.csv")

stasis_only<-indata_3_modes[indata_3_modes[, "mode"] == "Stasis",]
BM_only<-indata_3_modes[indata_3_modes[, "mode"] == "Random walk",]
DT_only<-indata_3_modes[indata_3_modes[, "mode"] == "Directional",]


# Check models that showed both a realative and abolsute fit to the data:
adequate_stasis.models<-rep(NA, nrow(stasis_only))
adequate_BM.models<-rep(NA, nrow(BM_only))
adequate_DT.models<-rep(NA, nrow(DT_only))


for (i in 1:length(stasis_only[,1])){
  
  if (stasis_only[i,13] == "FAILED" | stasis_only[i,18]=="FAILED" | stasis_only[i,23]=="FAILED" | stasis_only[i,28]=="FAILED") adequate_stasis.models[i]<-FALSE else adequate_stasis.models[i]<-TRUE
  
}

for (i in 1:length(BM_only[,1])){
  
  if (BM_only[i,34] == "FAILED" | BM_only[i,39]=="FAILED" | BM_only[i,44]=="FAILED") adequate_BM.models[i]<-FALSE else adequate_BM.models[i]<-TRUE
  
}

for (i in 1:length(DT_only[,1])){
  
  if (DT_only[i,51] == "FAILED" | DT_only[i,56]=="FAILED" | DT_only[i,61]=="FAILED") adequate_DT.models[i]<-FALSE else adequate_DT.models[i]<-TRUE
  
}


# Combine data frames:
temp_data.frame<-as.data.frame(rbind(stasis_only, BM_only ,DT_only))

new_3_modes<-cbind(temp_data.frame,c(adequate_stasis.models, adequate_BM.models, adequate_DT.models))
dim(new_3_modes)
names(new_3_modes)[68]<-"adequate"

stasis_models<-cbind(stasis_only, adequate_stasis.models)
BM_models<-cbind(BM_only, adequate_BM.models)
DT_models<-cbind(DT_only, adequate_DT.models)

# Remove lineages that failed to fit the model adequacy criteria:
new_3_modes_adequate<-new_3_modes[new_3_modes[, "adequate"] == "TRUE",]
dim(new_3_modes_adequate)
table(new_3_modes_adequate$mode)

### Done preparing data frames ###



######### 
######### Analyze time series
#########

########    STASIS    ########    

#Filter out the data sets that fitted the stasis model in an adequate way
stasis_adequate<-stasis_models[stasis_models[, "adequate_stasis.models"] == "TRUE",]
dim(stasis_adequate)


########    BM    ########   

#Filter out the data sets that fitted the BM model in an adequate way
BM_adequate<-BM_models[BM_models[, "adequate_BM.models"] == "TRUE",]
dim(BM_adequate)

########    Directional change    ######## 

#Filter out the data sets that fitted the DT model in an adequate way
DT_adequate<-DT_models[DT_models[, "adequate_DT.models"] == "TRUE",]
dim(DT_adequate)


####### Case studies

### Comparing the rate of evolution in coccoliths

# Success BM

indata<-read.table(temp[6], header=T)
y<-as.paleoTS(indata$mm, indata$vv, indata$N, indata$age.in.MY)

if (indata$log[1] =="no")   y<-ln.paleoTS(y);
if (indata$dimension[1] =="area") y$mm<-y$mm/2; y$vv<-y$vv/4;  
par(mfrow=c(1,1))
plot.paleoTS(y)

fit3models(y, method = "Joint")
fit3adequacy.RW(y, nrep=1000, plot = TRUE)



# Failure Stasis
indata<-read.table(temp[65], header=T)
y<-as.paleoTS(indata$mm, indata$vv, indata$N, indata$age.in.MY)

if (indata$log[1] =="no")   y<-ln.paleoTS(y);
if (indata$dimension[1] =="area") y$mm<-y$mm/2; y$vv<-y$vv/4;  
plot.paleoTS(y)

fit3models(y, method = "Joint", pool=FALSE)
fit4adequacy.stasis(y,nrep=1000, plot = TRUE)

fit3adequacy.RW(y, nrep=1000, plot = TRUE)

fit3adequacy.trend(y, nrep=1000, plot = TRUE)

# Failure DT
indata<-read.table(temp[48], header=T)
y<-as.paleoTS(indata$mm, indata$vv, indata$N, indata$age.in.MY)

if (indata$log[1] =="no")   y<-ln.paleoTS(y);
if (indata$dimension[1] =="area") y$mm<-y$mm/2; y$vv<-y$vv/4;  
layout(1:1)
plot.paleoTS(y)

fit3models(y, method = "Joint")
fit3adequacy.trend(y, nrep=1000, plot = TRUE)

fit3adequacy.RW(y, nrep=1000, plot = TRUE)

###################
################### Model adequacy when not the best model according to AICc

### BM adequate model when stasis shows best relative fit

adequate_BM.stasis.best<-rep(NA, length(stasis_adequate[,1]))


for (i in 1:length(adequate_BM.stasis.best)){
  
  if (stasis_adequate[i,34] == "PASSED" & stasis_adequate[i,39]=="PASSED" & stasis_adequate[i,44]=="PASSED") adequate_BM.stasis.best[i]<-TRUE else adequate_BM.stasis.best[i]<-FALSE
  
}

table(adequate_BM.stasis.best)


### DT adequate model when stasis shows best relative fit

adequate_DT.stasis.best<-rep(NA, length(stasis_adequate[,1]))


for (i in 1:length(adequate_DT.stasis.best)){
  
  if (stasis_adequate[i,51] == "PASSED" & stasis_adequate[i,56]=="PASSED" & stasis_adequate[i,61]=="PASSED" ) adequate_DT.stasis.best[i]<-TRUE else adequate_DT.stasis.best[i]<-FALSE
  
}

table(adequate_DT.stasis.best)



### DT adequate model when BM shows best relative fit

adequate_DT.BM.best<-rep(NA, length(BM_adequate[,1]))


for (i in 1:length(adequate_DT.BM.best)){
  
  if (BM_adequate[i,51] == "PASSED" & BM_adequate[i,56]=="PASSED" & BM_adequate[i,61]=="PASSED" ) adequate_DT.BM.best[i]<-TRUE else adequate_DT.BM.best[i]<-FALSE
  
}

table(adequate_DT.BM.best)


### stasis adequate model when BM shows best relative fit

adequate_stasis.BM.best<-rep(NA, length(BM_adequate[,1]))


for (i in 1:length(adequate_DT.BM.best)){
  
  if (BM_adequate[i,13] == "PASSED" & BM_adequate[i,18]=="PASSED" & BM_adequate[i,23]=="PASSED" & BM_adequate[i,28]=="PASSED") adequate_stasis.BM.best[i]<-TRUE else adequate_stasis.BM.best[i]<-FALSE
  
}

table(adequate_stasis.BM.best)



### BM adequate model when DT shows best relative fit

adequate_BM.DT.best<-rep(NA, length(DT_adequate[,1]))


for (i in 1:length(adequate_BM.DT.best)){
  
  if (DT_adequate[i,34] == "PASSED" & DT_adequate[i,39]=="PASSED" & DT_adequate[i,44]=="PASSED") adequate_BM.DT.best[i]<-TRUE else adequate_BM.DT.best[i]<-FALSE
  
}

table(adequate_BM.DT.best)


### stasis adequate model when DT shows best relative fit

adequate_stasis.DT.best<-rep(NA, length(DT_adequate[,1]))


for (i in 1:length(adequate_stasis.DT.best)){
  
  if (DT_adequate[i,13] == "PASSED" & DT_adequate[i,18]=="PASSED" & DT_adequate[i,23]=="PASSED" & DT_adequate[i,28]=="PASSED") adequate_stasis.DT.best[i]<-TRUE else adequate_stasis.DT.best[i]<-FALSE
  
}

table(adequate_stasis.DT.best)

