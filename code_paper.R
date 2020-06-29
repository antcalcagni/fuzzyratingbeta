rm(list=ls());graphics.off()
setwd("/home/antonio/MEGA/Lavoro_sync/My papers/Submitted/Modeling random and non-random uncertainty in ratings data_A fuzzy beta model/code_paper/")
source("fuzzybetareg.R")

##### Simulation study #####
# X: estimates (list)
# Y: simulated dataset, matrix of dimension 5 algorithms x J+H parameters x B replicates

## Compute summary statistics
res = list()
for(k in 1:16){
  print(paste0("Processing design: ",k))
  load(paste0("simdata/datasim_",k,".rda"))
  B=dim(Y)[3]
  
  for(j in 1:5){ #remove abnormal estimates in betareg
    for(jh in 1:(X$J+X$H)){Y[j,jh,abs(Y[j,jh,])>=10]=NA}
  }
  
  NAestimates = matrix(0,5,1); rownames(NAestimates) = rownames(Y); NAestimates_iid = list()
  for(j in 1:5){
    iid = which(apply(apply(Y[j,1:(X$J+X$H),],2,is.na),2,sum)>0)
    NAestimates[j,] = length(iid)  
    NAestimates_iid[[j]] = iid
  }
  
  ratioBias = matrix(NA,5,4); rownames(ratioBias) = rownames(Y); colnames(ratioBias) = c("r_beta","over_beta","r_gamma","over_gamma");
  for(j in 1:5){
    ratioBias[j,1] = sum((Y[j,1:X$J,setdiff(1:B,NAestimates_iid[[j]])]-X$b)>0)/sum((Y[j,1:X$J,setdiff(1:B,NAestimates_iid[[j]])]-X$b)<0)
    ratioBias[j,2] = sum((Y[j,1:X$J,setdiff(1:B,NAestimates_iid[[j]])]-X$b)>0)/length((Y[j,1:X$J,setdiff(1:B,NAestimates_iid[[j]])]))
    ratioBias[j,3] = sum((Y[j,(X$J+1):(X$J+X$H),setdiff(1:B,NAestimates_iid[[j]])]-X$g)>0)/sum((Y[j,(X$J+1):(X$J+X$H),setdiff(1:B,NAestimates_iid[[j]])]-X$g)<0)
    ratioBias[j,4] = sum((Y[j,(X$J+1):(X$J+X$H),setdiff(1:B,NAestimates_iid[[j]])]-X$g)<0)/length((Y[j,(X$J+1):(X$J+X$H),setdiff(1:B,NAestimates_iid[[j]])]))
  }
  
  
  theta_mean=theta_var=theta_bias=theta_rmse=matrix(NA,5,(X$J+X$H))
  colnames(theta_mean) = colnames(theta_var) = colnames(theta_bias) = colnames(theta_rmse) = c(paste0("beta",1:X$J),paste0("gamma",1:X$H)) 
  rownames(theta_mean) = rownames(theta_var) = rownames(theta_bias) = rownames(theta_rmse) = rownames(Y)
  theta_armse=theta_pa=matrix(NA,5,2); colnames(theta_armse) = colnames(theta_pa) = c("beta","gamma"); rownames(theta_armse) = rownames(theta_pa) = rownames(Y)
  for(j in 1:5){
    theta_mean[j,] = apply(Y[j,1:(X$J+X$H),setdiff(1:B,NAestimates_iid[[j]])],1,function(x)mean(x))
    theta_var[j,] = apply(Y[j,1:(X$J+X$H),setdiff(1:B,NAestimates_iid[[j]])],1,function(x)var(x))
    theta_bias[j,] = theta_mean[j,]-c(X$b,X$g)
    theta_rmse[j,] = sqrt(apply(mapply(function(b){(Y[j,1:(X$J+X$H),b]-c(X$b,X$g))^2},setdiff(1:B,NAestimates_iid[[j]])),1,function(x)mean(x)))
    theta_armse[j,1] = mean(sqrt(apply(as.matrix(mapply(function(b){((Y[j,1:X$J,b]-X$b)/X$b)^2},setdiff(1:B,NAestimates_iid[[j]]))),2,function(x)mean(x)))) #over beta
    theta_armse[j,2] = mean(sqrt(apply(matrix(mapply(function(b){((Y[j,(X$J+1):(X$J+X$H),b]-X$g)/X$g)^2},setdiff(1:B,NAestimates_iid[[j]])),nrow=X$H),1,function(x)mean(x)))) #over gamma
  }
  
  Ystar=list()
  
  res[[k]] = list(n=X$n,J=X$J,H=X$H,theta_mean=theta_mean,theta_var=theta_var,theta_bias=theta_bias,theta_rmse=theta_rmse,theta_armse=theta_armse,
                  NAestimates=NAestimates,NAestimates_iid=NAestimates_iid,ratioBias=ratioBias)
} 

## NA's 
NA_overall = cbind(X$design,t(mapply(function(k)res[[k]]$NAestimates,1:length(res))))
colnames(NA_overall)[4:8] = rownames(res[[1]]$NAestimates)
print(NA_overall)


## Table 1: Beta
iid = c(which(X$design[,1]==50),which(X$design[,1]==100),which(X$design[,1]==250),which(X$design[,1]==500))
data_table = data.frame(armse = t(mapply(function(k)as.matrix(res[[k]]$theta_armse[c(4,2,3,5),1]),iid)),
                        bias = t(mapply(function(k)apply(as.matrix(res[[k]]$theta_bias[c(4,2,3,5),1:res[[k]]$J]),1,mean),iid)),
                        matrix(mapply(function(k)c(res[[k]]$n,res[[k]]$J,res[[k]]$H),iid),ncol=3,byrow = TRUE))
names(data_table)[c(1:4,9:11)] = c("armse.fuzzy","armse.defuzz_mean","armse.defuzz_mode","armse.crisp_re","n","J","H");

print(round(data_table[data_table$J==2,c(9,11,5,1,6,2,7,3,8,4)],3),row.names = FALSE)

rmse_ov_beta = round(apply(data_table[,c(1:4)],2,mean),3)


## Table 1: Gamma
data_table = data.frame(armse = t(mapply(function(k)as.matrix(res[[k]]$theta_armse[c(4,2,3,5),2]),iid)),
                        bias = t(mapply(function(k)apply(as.matrix(res[[k]]$theta_bias[c(4,2,3,5),(res[[k]]$J+1):(res[[k]]$J+res[[k]]$H)]),1,mean),iid)),
                        matrix(mapply(function(k)c(res[[k]]$n,res[[k]]$J,res[[k]]$H),iid),ncol=3,byrow = TRUE))
names(data_table)[c(1:4,9:11)] = c("armse.fuzzy","armse.defuzz_mean","armse.defuzz_mode","armse.crisp_re","n","J","H");

print(round(data_table[data_table$J==2,c(9,11,5,1,6,2,7,3,8,4)],3),row.names = FALSE)

rmse_ov_gamma = round(apply(data_table[,c(1:4)],2,mean),3)


## Table 2: Beta
iid = c(which(X$design[,1]==50),which(X$design[,1]==100),which(X$design[,1]==250),which(X$design[,1]==500))
data_table = data.frame(armse = t(mapply(function(k)as.matrix(res[[k]]$theta_armse[c(4,2,3,5),1]),iid)),
                        bias = t(mapply(function(k)apply(as.matrix(res[[k]]$theta_bias[c(4,2,3,5),1:res[[k]]$J]),1,mean),iid)),
                        matrix(mapply(function(k)c(res[[k]]$n,res[[k]]$J,res[[k]]$H),iid),ncol=3,byrow = TRUE))
names(data_table)[c(1:4,9:11)] = c("armse.fuzzy","armse.defuzz_mean","armse.defuzz_mode","armse.crisp_re","n","J","H");

print(round(data_table[data_table$J==4,c(9,11,5,1,6,2,7,3,8,4)],3),row.names = FALSE)


## Table 2: Gamma
data_table = data.frame(armse = t(mapply(function(k)as.matrix(res[[k]]$theta_armse[c(4,2,3,5),2]),iid)),
                        bias = t(mapply(function(k)apply(as.matrix(res[[k]]$theta_bias[c(4,2,3,5),(res[[k]]$J+1):(res[[k]]$J+res[[k]]$H)]),1,mean),iid)),
                        matrix(mapply(function(k)c(res[[k]]$n,res[[k]]$J,res[[k]]$H),iid),ncol=3,byrow = TRUE))
names(data_table)[c(1:4,9:11)] = c("armse.fuzzy","armse.defuzz_mean","armse.defuzz_mode","armse.crisp_re","n","J","H");

print(round(data_table[data_table$J==4,c(9,11,5,1,6,2,7,3,8,4)],3),row.names = FALSE)


## Table 3
rmse_ov_beta
rmse_ov_gamma
fem = apply(mapply(function(k)res[[k]]$ratioBias[4,],1:length(res)),1,mean)
dml_mean = apply(mapply(function(k)res[[k]]$ratioBias[2,],1:length(res)),1,mean)
dml_mode = apply(mapply(function(k)res[[k]]$ratioBias[3,],1:length(res)),1,mean)
dreml = apply(mapply(function(k)res[[k]]$ratioBias[5,],1:length(res)),1,mean)


##### Case study 1 #####

## data
load("casestudy1_data.rda")

## models
fm1 = "RDB_extremepos_m,RDB_extremepos_s ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | 1"
mod1 = fuzzybetareg_fit(formula = fm1,data = data_rd,verbose=TRUE,init = "fixed")
summary_fuzzybetareg(mod1)

fm2 = "RDB_extremepos_m,RDB_extremepos_s ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | sex"
LRtest_fuzzybetareg(fm2,fm1,data_rd,init="fixed")
mod2 = fuzzybetareg_fit(formula = fm2,data = data_rd,verbose=TRUE,init = "fixed")
summary_fuzzybetareg(mod2)

fm3 = "RDB_extremepos_m,RDB_extremepos_s ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | sex+frequency_driving"
LRtest_fuzzybetareg(fm3,fm2,data_rd,init="fixed")
mod3 = fuzzybetareg_fit(formula = fm3,data = data_rd,verbose=TRUE,init = "fixed")
summary_fuzzybetareg(mod3)


##### Case study 2 #####

## data
load("casestudy2_data.rda")

## models
fm1 = "QR7_m,QR7_s ~  quality_food+quality_employees  | restaurant_type+pricebool"
mod1 = fuzzybetareg_fit(formula = fm1,data = data_cs,verbose=FALSE,init="fixed")
summary_fuzzybetareg(mod1)



