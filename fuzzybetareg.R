beta_fn = function(x,mi,gi){
  a = 1+gi*mi; b= 1+gi*(1-mi)  
  C = (((a-1)/(a+b-2))^(a-1) * (1-((a-1)/(a+b-2)))^(b-1))
  fx = 1/C * x^(a-1) * (1-x)^(b-1)
  return(fx)
}

tpzfun = function(x,a,b,c,d){
  mux=trapezoid::dtrapezoid(x=x,min=a,max=d,mode1=b,mode2=c)
  mux=mux/max(mux)
  return(mux)
}

defuzzify_centroid = function(mi,si,what=c("mean","mode")){
  a1 = 1+si*mi; b1 = 1+si*(1-mi) 
  if(what=="mean"){(a1)/(a1+b1)}
  else if(what=="mode"){(a1-1)/(a1+b1-2)}
}


compute_weights = function(x,m,s){
  fx = beta_fn(x,m,s) #compute weight as membership of x in fuzzy set (m,s)
  if(is.infinite(fx)|is.nan(fx)){fx=ifelse(abs(m-x)<1e-9,1,0)} #check whether fuzzy set is degenerate in m
  return(fx)
}

compute_residuals = function(m,s,mu_est,phi_est){
  # Pearson like residual for fuzzy observations:
  # difference between mode and estimated mu weighted by membership function computed for mu.
  # Formally: (m - mu_hat)*1-xi_ytilde{mu_hat} where xi_ytilde{mu_hat} weights how much mu_hat belongs to fuzzy set y_tilde
  ws = mapply(function(i)compute_weights(mu_est[i],m[i],s[i]),1:length(m))
  var_yhat = (mu_est*(1-mu_est)/(1+phi_est))
  resx = ((m-mu_est)*(1-ws))/var_yhat
  return(resx)
}

ystar1_compute = function(mu,phi,m,s){
  a2 = mu*phi; b2 = phi*(1-mu) #parameters for beta(y;a2,b2)
  a1 = 1+s*m; b1 = 1+s*(1-m) 
  eta = a1+a2-1; nu = b1+b2-1
  log(eta/(eta+nu)) - (nu/(2*eta*(1+eta+nu)))
  
  
}

ystar2_compute = function(mu,phi,m,s){
  a2 = mu*phi; b2 = phi*(1-mu) #parameters for beta(y;a2,b2)
  a1 = 1+s*m; b1 = 1+s*(1-m) 
  eta = a1+a2-1; nu = b1+b2-1
  log(nu/(eta+nu)) - (eta/(2*nu*(1+eta+nu)))
}

loglikel_taylor_dispVar = function(X,b,Z,g,ystar1,ystar2){
  #c(print(exp(-X%*%b)),print(exp(Z%*%g)))
  sum(log(gamma(exp(Z%*%g))) - log(gamma((1/(1+exp(-X%*%b)))*exp(Z%*%g))) - log(gamma(exp(Z%*%g)-(1/(1+exp(-X%*%b)))*exp(Z%*%g))) + 
        ystar1*(exp(Z%*%g)/(exp(-X%*%b) + 1) - 1) - ystar2*(exp(Z%*%g)*(1/(exp(-X%*%b) + 1) - 1) + 1))
} #loglikelihood with plugin of Taylor expanded expectations [expected complete-data loglikelihood]

diff_likel_beta_dispVar = function(X,b,Z,g,ystar1,ystar2){
  t(X)%*%((exp(X%*%b + Z%*%g)*(ystar1 - ystar2 - 
                                 pracma::psi(exp(X%*%b + Z%*%g)/(exp(X%*%b) + 1)) + 
                                 pracma::psi(exp(Z%*%g)/(exp(X%*%b) + 1))))/(exp(X%*%b) + 1)^2)  
}

diff_likel_gamma_dispVar = function(X,b,Z,g,ystar1,ystar2){
  t(Z)%*%( (pracma::psi(exp(Z%*%g))*exp(Z%*%g) - 
              pracma::psi(exp(Z%*%g) - exp(Z%*%g)/(exp(-X%*%b) + 1))*exp(Z%*%g) + 
              pracma::psi(exp(Z%*%g) - exp(Z%*%g)/(exp(-X%*%b) + 1))*(exp(Z%*%g))/(exp(-X%*%b) + 1) + 
              (ystar1*exp(Z%*%g))/(exp(-X%*%b) + 1) - 
              (exp(Z%*%g)*pracma::psi(exp(Z%*%g)/(exp(-X%*%b) + 1)))/(exp(-X%*%b) + 1) - ystar2*exp(Z%*%g)*(1/(exp(-X%*%b) + 1) - 1) ) )
}

compute_SEs = function(X,Z,m,s,beta_est,gamma_est){
  ## See: McLachlan & Peel (2000). Finite Mixture Models. Wiley (pp. 64-67)
  J=length(beta_est); H=length(gamma_est);n=NROW(X)
  
  # compute expectations
  mux_est = plogis(X%*%beta_est); phix_est = exp(Z%*%gamma_est)
  ystar1_est = mapply(function(i)ystar1_compute(mux_est[i],phix_est[i],m[i],s[i]),1:n)
  ystar2_est = mapply(function(i)ystar2_compute(mux_est[i],phix_est[i],m[i],s[i]),1:n)
  
  # compute empirical information matrix as sum of (1,..i,..n) squared expected score statistics
  I_empirical = matrix(0,J+H,J+H)
  for(i in 1:n){
    sc = rbind(diff_likel_beta_dispVar(X = matrix(X[i,],1,J),b = beta_est,Z = matrix(Z[i,],1,H),g = gamma_est,ystar1 = ystar1_est[i],ystar2 = ystar2_est[i]),
               diff_likel_gamma_dispVar(X = matrix(X[i,],1,J),b = beta_est,Z = matrix(Z[i,],1,H),g = gamma_est,ystar1 = ystar1_est[i],ystar2 = ystar2_est[i]))
    I_empirical = I_empirical + sc%*%t(sc)
  }
  
  # compute asymptotic standard errors
  ses = matrix(sqrt(diag(solve(I_empirical))),J+H,1)
  return(list(SEs=ses,ICov=solve(I_empirical)))
}

fuzzybetareg_fit = function(formula=NULL,data=NULL,K=250,init=c("random","fixed"),eps=1e-7,verbose=TRUE,use.optim=FALSE,bfgs.lower=NULL,bfgs.upper=NULL){
  
  ## Extract data from formula object
  out = modelframe_fuzzybetareg(formula,data)
  X = out$X; Z = out$Z
  s = out$s; m = out$m
  J=NCOL(X); H=NCOL(Z); n=NROW(X)
  
  ## Starting values for beta and gamma
  type_init = match.arg(init)
  if(type_init=="random"){beta0=runif(J,-0.5,0.5); gamma0=runif(H,-0.5,0.5)}
  else if(type_init=="fixed"){
    ydefuzz1 = defuzzify_centroid(m,s,"mean"); if(!(min(ydefuzz1)>0&max(ydefuzz1)<1)){ydefuzz1[ydefuzz1>=1]=0.999;ydefuzz1[ydefuzz1<=0]=1-0.999}
    out0=betareg::betareg(ydefuzz1~-1+X|-1+Z); beta0=out0$coefficients$mean;gamma0=out0$coefficients$precision
  }
  if(use.optim){if(is.null(bfgs.lower)){bfgs.lower = c(rep(-5,J),rep(-5,H))}; if(is.null(bfgs.upper)){bfgs.upper = c(rep(5,J),rep(5,H))}}
  
  ## EM routine
  Thetax = matrix(NA,K,J+H); loglikx = matrix(NA,K,1)
  conv=1; k=2; 
  Thetax[1,] = c(beta0,gamma0); loglikx[1] = -100000; 
  while(conv>eps){
    if(verbose){message(paste0("@ Iter no.: ",k))}
    ## Compute expectations
    mux_k = plogis(X%*%Thetax[k-1,1:J])
    phix_k = exp(Z%*%Thetax[k-1,(J+1):(J+H)])
    ystar1 = mapply(function(i)ystar1_compute(mux_k[i],phix_k[i],m[i],s[i]),1:n) # E[ln Y|y_tilde] = ln E[Y|y_tilde] - V[Y|y_tilde] * (2*E[Y|y_tilde]^2)^-1
    ystar2 = mapply(function(i)ystar2_compute(mux_k[i],phix_k[i],m[i],s[i]),1:n)# E[ln Y|y_tilde] = ln (1-E[Y|y_tilde]) - V[Y|y_tilde] * (2*(1-E[Y|y_tilde])^2)^-1
    ## Maximization 
    if(use.optim==FALSE){
    if(J>1){b_k = pracma::fsolve(function(x){diff_likel_beta_dispVar(X,x,Z,Thetax[k-1,(J+1):(J+H)],ystar1,ystar2)},x0 = Thetax[k-1,1:J])$x;}
    else{b_k = uniroot(function(x){diff_likel_beta_dispVar(X,x,Z,Thetax[k-1,(J+1):(J+H)],ystar1,ystar2)},lower = -5,upper = 5)$root}
    if(H>1){g_k = pracma::fsolve(function(x){diff_likel_gamma_dispVar(X,b_k,Z,x,ystar1,ystar2)},x0 = Thetax[k-1,(J+1):(J+H)])$x}
    else{g_k = uniroot(function(x){diff_likel_gamma_dispVar(X,b_k,Z,x,ystar1,ystar2)},lower = -5,upper = 5)$root}
    Thetax[k,] = c(b_k,g_k)
    }else{
    Thetax[k,] = optim(par=rep(0.1,J+H),fn=function(x){loglikel_taylor_dispVar(X,x[1:J],Z,x[(J+1):(J+H)],ystar1,ystar2)},method = "L-BFGS-B",lower = bfgs.lower,upper = bfgs.upper,control=list(trace=0,fnscale=-1))$par
    }
    if(verbose){message(paste0("Current params: ", paste(round(Thetax[k,],4),collapse = ", ") ))}
    ## Compute loglik and check for convergence
    loglikx[k] = loglikel_taylor_dispVar(X,Thetax[k,1:J],Z,Thetax[k,(J+1):(J+H)],ystar1,ystar2)
    if(verbose){message(paste0("Current loglik: ", round(loglikx[k],4)),"\n")}
    conv = abs(loglikx[k]-loglikx[k-1])/abs(loglikx[k])
    if(k>=K){break};
    k=k+1
  }
  ## Prepare output
  outlist = list(init_par=list(beta=beta0,gamma=gamma0),Theta=Thetax[1:(k-1),],loglikel=loglikx[1:(k-1)],beta_EM=Thetax[k-1,1:J],gamma_EM=Thetax[k-1,(J+1):(J+H)],conv=ifelse(conv<=eps,1,-1))
  ses = compute_SEs(X,Z,m,s,Thetax[k-1,1:J],Thetax[k-1,(J+1):(J+H)])
  df = NCOL(X)+NCOL(Z)
  
  if(outlist$conv==1|k<K){xconv=1}else{xconv=0}
  
  ## Output
  out=list()
  out$coefficients$mean = Thetax[k-1,1:J]; attributes(out$coefficients$mean) = list(names=colnames(X))
  out$coefficients$precision = Thetax[k-1,(J+1):(J+H)]; attributes(out$coefficients$precision) = list(names=colnames(Z))
  out$residuals = compute_residuals(m,s,mux_k,phix_k);
  out$type.residuals = "Pearson"
  out$fitted.values = list(mean = mux_k, precision = phix_k)
  out$type = "ML via EM-fuzzy"
  out$use.optim = use.optim 
  out$EM = outlist
  out$loglikel = outlist$loglikel[length(outlist$loglikel)]
  out$AIC = -2*out$loglikel+2*df
  out$control = list(method="EM-fuzzy",maxit=K,epstol=eps,init=type_init)
  out$start = c(beta0,gamma0); attributes(out$start) = list(names=c(colnames(X),colnames(Z)))
  out$n = n
  out$df = df
  out$vcov = ses$ICov
  out$standard.errors = ses$SEs; attributes(out$standard.errors) = list(names=c(colnames(X),colnames(Z)))
  out$converged = xconv
  out$data = list(m=m,s=s,X=X,Z=Z)
  out$formula = formula
  return(out)
}


LRtest_fuzzybetareg = function(formula1,formula0,data,use.optim=FALSE,init="random"){
  ## m0
  out0 = fuzzybetareg_fit(formula0,data,verbose = FALSE,use.optim = use.optim,init=init)
  ll0 = out0$loglikel
  df0 = out0$df
  ## m1
  out1 = fuzzybetareg_fit(formula1,data,verbose = FALSE,use.optim = use.optim,init=init)
  ll1 = out1$loglikel
  df1 = out1$df
  ## lrtest
  diffll = 2*abs(ll1-ll0)
  diffdf = df1-df0
  pdiffll = pchisq(diffll,diffdf,lower.tail = FALSE)
  ## anova() style table
  res = matrix(NA,2,6); colnames(res)=c("#Df","LogLik","AIC","Df","Chisq","Pr(>Chisq)"); rownames(res)=c("m0","m1")
  res[1,] = c(df0,ll0,out0$AIC,NA,NA,NA)
  res[2,] = c(df1,ll1,out1$AIC,diffdf,diffll,pdiffll)
  
  formulas = c(paste0("Model0: ~ ",strsplit(x = formula0,split = "~")[[1]][2]),paste0("Model1: ~ ",strsplit(x = formula1,split = "~")[[1]][2],"\n"))
  structure(as.data.frame(res), heading = c("Likelihood ratio test\n", formulas),class = c("anova", "data.frame"))
  
}

PseudoR2 = function(formula=NULL,data=NULL,use.optim=FALSE,init="random"){
  
  fu=unlist(strsplit(x=unlist(strsplit(x=formula,split="~"))[1],split=","))
  if(length(fu)>1){ #fuzzy-betareg
    ## m1
    out1 = fuzzybetareg_fit(formula,data,verbose = FALSE,use.optim = use.optim,init=init)
    ll1 = out1$loglikel
    ## m0 (model with a constant only)
    fm0 = paste0(unlist(strsplit(x=formula,split="~"))[1],"~1|1")
    out0 = fuzzybetareg_fit(fm0,data,verbose = FALSE,use.optim = use.optim,init=init)
    ll0 = out0$loglikel
  }else{ #crisp-betareg
    ## m1
    out1 = betareg::betareg(formula,data)
    ll1 = out1$loglik
    ## m0 (model with a constant only)
    fm0 = paste0(unlist(strsplit(x=formula,split="~"))[1],"~1|1")
    out0 = betareg::betareg(fm0,data)
    ll0 = out0$loglik
  }
  
  ## LRT-based pseudoR2 (e.g., see: Veall, M. R., & Zimmermann, K. F. (1994). Evaluating Pseudo-R 2's for binary probit models. Quality and Quantity, 28(2), 151-164.)
  n = NROW(data)
  lrt = 2*(ll1-ll0)
  rs = 2*ll0/n
  pr2 = (-lrt*(1-rs))/(rs*(lrt+n))
  
  return(pr2)
}

predictive_check = function(B=5000,model=NULL,Ydata=NULL,seedx=NULL,alphax=0.01){
  if(!is.null(seedx)){set.seed(seedx)}
  
  n=NROW(Ydata)
  mux=model$fitted.values$mean
  phix=model$fitted.values$precision
  Yhat = t(mapply(function(i)rbeta(B,mux[i]*phix[i],phix[i]*(1-mux[i])),1:n))
  
  xsup=seq(0,1,length.out = 500)
  Fy_0 = t(mapply(function(i){
    fy=beta_fn(x = xsup,mi = Ydata[i,1],gi = Ydata[i,2])
    return(c(min(xsup[fy>alphax]),max(xsup[fy>alphax])))
    },1:n))
  
  yhat_fy = mapply(function(i)sum(Yhat[i,]>=Fy_0[i,1]&Yhat[i,]<=Fy_0[i,2])/B,1:n)
  return(list(fit_i=yhat_fy,fit_ov=mean(yhat_fy),Ypredict=Yhat,Yobs_alpha=Fy_0))
}


AIC_fuzzybetareg = function(X,Z,m,s,beta_par,gamma_par){
  mux = plogis(X%*%beta_par); phix = exp(Z%*%gamma_par)
  ll = sum(mapply(function(i)log(integrate(function(x)beta_fn(x,m[i],s[i])*dbeta(x,mux[i]*phix[i],phix[i]*(1-mux[i])),0,1)$value),1:length(m)))
  df = NCOL(X)+NCOL(Z)
  aic = -2*ll+2*df
  return(aic)
}


summary_fuzzybetareg = function(x,digits=3){
  
  ## extend coefficient table
  J = length(x$coefficients$mean)
  H = length(x$coefficients$precision)
  cf = as.vector(do.call("c", x$coefficients))
  se = x$standard.errors
  cf = cbind(cf, se, cf/se, 2*pnorm(-abs(cf/se))); colnames(cf) = c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf = list(mean = cf[seq.int(length.out = J), , drop = FALSE], precision = cf[seq.int(length.out = H) + J, , drop = FALSE])
  rownames(cf$mean) = names(x$coefficients$mean)
  rownames(cf$precision) <- names(x$coefficients$precision)
  x$coefficients = cf
  
  ## table to print
  if(!x$converged) {
    cat("Model did not converge\n")
  } else {
    cat("Beta variable dispersion model fitted on fuzzy observations via fuzzy-EM algorithm\n")
    cat("\nFormula:", x$formula,"\n")
    cat(sprintf("\nResiduals (%s):\n", x$type.residuals))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
                    .Names = c("Min", "1Q", "Median", "3Q", "Max")))
    
    if(NROW(x$coefficients$mean)) {
      cat(paste("\nMu Coefficients (mean model):\n", sep = ""))
      printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in mean model)\n")
    
    if(NROW(x$coefficients$precision)) {
      cat(paste("\nPhi coefficients (precision model):\n", sep = ""))
      printCoefmat(x$coefficients$precision, digits = digits, signif.legend = FALSE)
    } else cat("\nNo coefficients (in precision model)\n")
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")
    
    cat("\nType of estimator:", x$type)
    cat("\nLog-likelihood:", formatC(x$loglik, digits = digits),
        "on", x$df, "Df", sprintf("(Num. of observations: %s)", x$n))
    cat("\nAIC:",x$AIC)
    if(x$use.optim){
      cat(paste("\nNumber of iterations:", NROW(x$EM$Theta)+1),"(Likelihood equations solved via BFGS)")  
    }else{
      cat(paste("\nNumber of iterations:", NROW(x$EM$Theta)+1),"(Likelihood equations solved via Newton method)")
    }
    
  }
  invisible(x)
}


modelframe_fuzzybetareg = function(formula,data){
  ## formula like m,s ~ x1+..+xJ | z1+..+zH
  
  formula_lhs = strsplit(x = formula,split = "~")[[1]][1]
  formula_rhs = strsplit(x = formula,split = "~")[[1]][2]
  
  yvar=strsplit(x = formula_lhs,split = ",")[[1]]
  m=data[,gsub(x = yvar[1],pattern = " ",replacement = "")==names(data)]
  s=data[,gsub(x = yvar[2],pattern = " ",replacement = "")==names(data)]
  
  fx = gsub(x = strsplit(x = formula_rhs,split = "\\|")[[1]][1],pattern = " ",replacement = "")
  if(fx!=""){
    xvar=strsplit(x = fx,split = "\\+")[[1]]
    X = model.matrix(as.formula(paste0("~",paste(xvar,collapse="+"))),data=data)
  }else{
    ones = rep(1,NROW(X)) #mu is just scalar
    X = model.matrix(as.formula(paste0("~",paste(ones,collapse="+"))),data=data)
  }
  
  fz = gsub(x = strsplit(x = formula_rhs,split = "\\|")[[1]][2],pattern = " ",replacement = "")
  if(fz!=""){
    zvar=strsplit(x = fz,split = "\\+")[[1]]
    Z = model.matrix(as.formula(paste0("~",paste(zvar,collapse="+"))),data=data)
  }else{
    ones = rep(1,NROW(X)) #phi is just scalar
    Z = model.matrix(as.formula(paste0("~",paste(ones,collapse="+"))),data=data)
  }
  return(list(X=X,Z=Z,s=s,m=m))
}

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}
