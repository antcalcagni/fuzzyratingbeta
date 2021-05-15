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

## figure 2
tikzDevice::tikz(file='../AStA/fig1.tex',width=4.5,height=3)
set.seed(121);sbjs = c(sample(which(data_rd$RDB_extremepos_m<0.2),5),which(data_rd$RDB_extremepos_m>0.3))
i=sbjs[1];curve(beta_fn(x,data_rd$RDB_extremepos_m[i],data_rd$RDB_extremepos_s[i]),0,1,bty="n",xlab="",ylab="",lwd=1.5)
for(i in 2:length(sbjs)){curve(beta_fn(x,data_rd$RDB_extremepos_m[sbjs[i]],data_rd$RDB_extremepos_s[sbjs[i]]),0,1,add=TRUE,lwd=1.5)}
dev.off()


## models
data_rd$RDB_defuzz = mapply(function(i)defuzzify_centroid(data_rd$RDB_extremepos_m[i],data_rd$RDB_extremepos_s[i],what = "mean"),1:NROW(data_rd))

fm1 = "RDB_extremepos_m,RDB_extremepos_s ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | 1"
mod1 = fuzzybetareg_fit(formula = fm1,data = data_rd,verbose=TRUE,init = "fixed")
summary_fuzzybetareg(mod1)
PseudoR2(formula = fm1,data = data_rd,init = "fixed")

fm1_d = "RDB_defuzz ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | 1"
mod1_defuzz = betareg::betareg(formula = fm1_d,data=data_rd)
summary(mod1_defuzz)
PseudoR2(formula = fm1_d,data = data_rd,init = "fixed")

fm2 = "RDB_extremepos_m,RDB_extremepos_s ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | sex"
LRtest_fuzzybetareg(fm2,fm1,data_rd,init="fixed")
mod2 = fuzzybetareg_fit(formula = fm2,data = data_rd,verbose=TRUE,init = "fixed")
summary_fuzzybetareg(mod2)
PseudoR2(formula = fm2,data = data_rd,init = "fixed")
predictive_check(model=mod2,Ydata=cbind(data_rd$RDB_extremepos_m,data_rd$RDB_extremepos_s))$fit_ov

fm2_d = "RDB_defuzz ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | sex"
mod2_defuzz = betareg::betareg(formula = fm2_d,data=data_rd)
summary(mod2_defuzz)
PseudoR2(formula = fm2_d,data = data_rd,init = "fixed")

fm3 = "RDB_extremepos_m,RDB_extremepos_s ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | sex+frequency_driving"
LRtest_fuzzybetareg(fm3,fm2,data_rd,init="fixed")
mod3 = fuzzybetareg_fit(formula = fm3,data = data_rd,verbose=TRUE,init = "fixed")
summary_fuzzybetareg(mod3)
#PseudoR2(formula = fm3,data = data_rd,init = "fixed")

fm3_d = "RDB_defuzz ~ DAS_slowd_m+DAS_slowd_s+FCRS_cat+RDB_drugs_categ | sex+frequency_driving"
mod3_defuzz = betareg::betareg(formula = fm3_d,data=data_rd)
summary(mod3_defuzz)
#PseudoR2(formula = fm3_d,data = data_rd,init = "fixed")

## figure 2
modx=mod2 #choosing final model 
cols = c("orangered3","orchid4","palegreen4","peru")

ym = data_rd$RDB_extremepos_m; ys = data_rd$RDB_extremepos_s
x1 = data_rd$DAS_slowd_m; x2 = data_rd$DAS_slowd_s
c1 = as.numeric(data_rd$FCRS_cat); c2 = as.numeric(data_rd$RDB_drugs_categ);  

n = NROW(data_rd); c = rep(NA,n)
c[c1==1&c2==1]=cols[1]; c[c1==1&c2==2]=cols[2];
c[c1==2&c2==1]=cols[3]; c[c1==2&c2==2]=cols[4];

Yobs_a0 = predictive_check(model = mod1,Ydata = cbind(ym,ys),alphax = 0.5)$Yobs_alpha

tikzDevice::tikz(file='fig2.tex',width=6.5,height=4)
par(mfrow=c(1,2),mai=c(1.15, 0.85, 0.45, 0.15))

## fuzzy observations (alpha-cuts=0.5) - over anger (m)
I=n;jit=0.05
plot(x1[1:I],ym[1:I],bty="n",xlim=c(0,1),col="white",ylim=c(0,0.6),lwd=0,ylab="rdb",xlab="anger (m)",main="(A)")
rect(ytop = Yobs_a0[1:I,1],ybottom = Yobs_a0[1:I,2],xleft = x1[1:I],xright = x1[1:I]+jit,border=c,lty=2,lwd=2)

# line c1=1 & c2=1 (over x1)
bx = modx$coefficients$mean[1:2]
lines(seq(0,1,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[1],lwd=2.3)

# line c1=1 & c2=2 (over x1)
bx = c(modx$coefficients$mean[1]+modx$coefficients$mean[5],modx$coefficients$mean[2]); 
lines(seq(0,1,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[2],lwd=2.3)

# line c1=2 & c2=1 (over x1)
bx = c(modx$coefficients$mean[1]+modx$coefficients$mean[4],modx$coefficients$mean[2]); 
lines(seq(0,1,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[3],lwd=2.3)

# line c1=2 & c2=2 (over x1)
bx = c(modx$coefficients$mean[1]+modx$coefficients$mean[4]+modx$coefficients$mean[5],modx$coefficients$mean[2]); 
lines(seq(0,1,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[4],lwd=2.3)

# fuzzy observations (alpha-cuts=0.5) - over anger (s)
I=n;jit=15.25
plot(x2[1:I],ym[1:I],bty="n",xlim=c(0,250),col="white",ylim=c(0,0.6),lwd=0,ylab="rdb",xlab="anger (s)",main="(B)")
rect(ytop = Yobs_a0[1:I,1],ybottom = Yobs_a0[1:I,2],xleft = x2[1:I],xright = x2[1:I]+jit,border=c,lty=2,lwd=2)

# line c1=1 & c2=1 (over x2)
bx = c(modx$coefficients$mean[1],modx$coefficients$mean[3])
lines(seq(0,250,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[1],lwd=2.3)

# line c1=1 & c2=2 (over x2)
bx = c(modx$coefficients$mean[1]+modx$coefficients$mean[5],modx$coefficients$mean[3]); 
lines(seq(0,250,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[2],lwd=2.3)

# line c1=2 & c2=1 (over x2)
bx = c(modx$coefficients$mean[1]+modx$coefficients$mean[4],modx$coefficients$mean[3]); 
lines(seq(0,250,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[3],lwd=2.3)

# line c1=2 & c2=2 (over x2)
bx = c(modx$coefficients$mean[1]+modx$coefficients$mean[4]+modx$coefficients$mean[5],modx$coefficients$mean[3]); 
lines(seq(0,250,length.out = n),plogis(cbind(1,seq(0,1,length.out = n))%*%bx),col=cols[4],lwd=2.3)

add_legend("bottom",fill = cols,legend = c("fcrs:bad,drugs:non-use","fcrs:bad,drugs:use","fcrs:satisf,drugs:non-use","fcrs:satisf,drugs:use"),border = FALSE,bty = "n",ncol = 2)

dev.off()




##### Case study 2 #####

## data
load("casestudy2_data.rda")

## figure 3
tikzDevice::tikz(file='fig3.tex',width=4.5,height=3)
load("data_fig3.rda")
curve(beta_fn(x,datax$m[1],datax$s[1]),0,1,bty="n",xlab="",ylab="",lwd=1.5); curve(tpzfun(x,a=datax$inf0[1],d=datax$sup0[1],b=datax$inf1[1],c=datax$sup1[1]),0,1,bty="n",xlab="",ylab="",lwd=1.25,add=TRUE,lty=2)
for(i in 2:NROW(datax)){curve(beta_fn(x,datax$m[i],datax$s[i]),0,1,bty="n",xlab="",ylab="",lwd=1.5,add=TRUE); curve(tpzfun(x,a=datax$inf0[i],d=datax$sup0[i],b=datax$inf1[i],c=datax$sup1[i]),0,1,bty="n",xlab="",ylab="",lwd=1.25,add=TRUE,lty=2)}
dev.off()

## models
data_cs$QR7_defuzz = mapply(function(i)defuzzify_centroid(data_cs$QR7_m[i],data_cs$QR7_s[i],what = "mean"),1:NROW(data_cs))

fm1 = "QR7_m,QR7_s ~  quality_food + quality_employees  | restaurant_type + pricebool"
mod1 = fuzzybetareg_fit(formula = fm1,data = data_cs,verbose=FALSE,init="fixed")
summary_fuzzybetareg(mod1)
PseudoR2(formula = fm1,data = data_cs,init = "fixed")

fm1_d = "QR7_defuzz ~  quality_food + quality_employees  | restaurant_type + pricebool"
mod1_defuzz = betareg::betareg(formula = fm1_d,data = data_cs)
summary(mod1_defuzz)
PseudoR2(formula = fm1_d,data = data_cs,init = "fixed")
predictive_check(model=mod1,Ydata=cbind(data_cs$QR7_m,data_cs$QR7_s))$fit_ov

## figure 4
modx=mod1 #choosing final model 

ym = data_cs$QR7_m; ys = data_cs$QR7_s
x1 = data_cs$quality_food; x2 = data_cs$quality_employees
n = NROW(data_cs); 

Yobs_a0 = predictive_check(model = mod1,Ydata = cbind(ym,ys),alphax = 0.5)$Yobs_alpha

tikzDevice::tikz(file='fig4.tex',width=6.5,height=4)
par(mfrow=c(1,2),mai=c(1.15, 0.85, 0.45, 0.15))

## fuzzy observations (alpha-cuts=0.5) - over quality_food
I=n;jit=0.05
plot(x1[1:I],ym[1:I],bty="n",xlim=c(0.1,1.1),col="white",ylim=c(0.2,1.1),lwd=0,ylab="service quality",xlab="food",main="(A)")
rect(ytop = Yobs_a0[1:I,1],ybottom = Yobs_a0[1:I,2],xleft = x1[1:I],xright = x1[1:I]+jit,border="snow4",lty=2,lwd=2)

bx = modx$coefficients$mean[1:2];
lines(seq(0.2,1,length.out = n),plogis(cbind(1,seq(0.2,1,length.out = n))%*%bx),col="skyblue4",lwd=2.3)

# fuzzy observations (alpha-cuts=0.5) - over employees
I=n;jit=0.05
plot(x2[1:I],ym[1:I],bty="n",xlim=c(0.1,1.1),col="white",ylim=c(0.2,1.1),lwd=0,ylab="service quality",xlab="employees",main="(B)")
rect(ytop = Yobs_a0[1:I,1],ybottom = Yobs_a0[1:I,2],xleft = x2[1:I],xright = x2[1:I]+jit,border="snow4",lty=2,lwd=2)

bx = c(modx$coefficients$mean[1],modx$coefficients$mean[3]);
lines(seq(0.2,1,length.out = n),plogis(cbind(1,seq(0.2,1,length.out = n))%*%bx),col="skyblue4",lwd=2.3)

dev.off()
