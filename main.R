rm(list=ls())

# set.seed(1)
## Code used by UNC computing cluster to pull in inputs from Bash shell;
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);
node.idx <- Sys.getenv("SLURM_ARRAY_TASK_ID")

hypothesis='alternative' # or "null"

## set node.idx for running in windows environment (i.e., not on UNC computing cluster);
if (.Platform$OS.type == "windows") { node.idx = 1 }
# if (!is.na(node.idx)) node.idx = formatC(node.idx, width = 4, format = "d", flag = "0")


node.idx = as.numeric(node.idx)
node.idx = ifelse(is.na(node.idx), 3, node.idx)

# nchain = node.idx %% 5
targetEvents = ifelse(node.idx%%11==0,100,ifelse(node.idx%%11==1,110,ifelse(node.idx%%11==2,120,ifelse(node.idx%%11==3,130,ifelse(node.idx%%11==4,140,ifelse(node.idx%%11==5,150,ifelse(node.idx%%11==6,160,ifelse(node.idx%%11==7,170,ifelse(node.idx%%11==8,180,ifelse(node.idx%%11==9,190,ifelse(node.idx%%11==10,200)))))))))))

filename = paste0("CRDesign_realisticDatGen_B",node.idx,"_targetEvents",targetEvents,"_H",ifelse(hypothesis=='null',0,1))


library(mvtnorm)
library(lattice)
library(data.table) ; library(tidyr) ; library(tidyverse)
library(Rcpp) ; library(RcppNumerical)
library(survival) ; library(survminer) ; library(eha) # Survival

## ----------------------------------------------------------------------------------------- ##
## ----------------------------------------------------------------------------------------- ##

setwd('prespecified location')
sourceCpp('Gcopula_jointMAR.cpp')
sourceCpp('CureRateMCMC_cpp.cpp')

### For different scenarios, p, parameters(corr_r,lambda,gamma,beta), and Z changed. ###
## DATA SIZE
n = 500 ; # sample size
k = 2 ;   # number of promotion times
p = 2 ;   # number of covariates for Poisson mean

## PARAMETERS
corr_r       = -2.325 ;               # correlation for the Gaussian copula (off-diagonal)
invisible(ifelse(k==1,GammaCorrmat <- matrix(1,1,1),GammaCorrmat <- Corrmat(k,corr_r)$Gamma))
lambda        = c(1.596,1.742)#c(3,2)[1:k] ;            # Weibull shape
gamma       = c(0.653,0.098)#c(2,3)[1:k]#(0.8,0.5) ;            # promotion time Exponential mean
beta         = c(0.753,-0.494,1.233,-0.451) ; # covariates regression coefficients
if (hypothesis == "null") {
  beta[seq(p,length(beta),by=p)] = numeric(k)
} 

## study end
censored.time = 5 ; targetFollowUp = 3

## Supplementary Data ##
Z = matrix(1,nr=n,nc=p*k) ; # intercept
Z[,seq(2,ncol(Z),by=p)] = rbinom(n*2,1,0.5) # TREATMENT



theta = matrix(NA,n,k)
for (i in 1:k) {
  theta[,i] = exp(Z[,(1+p*(i-1)):(i*p)] %*% beta[(1+p*(i-1)):(i*p)])
}
upper = 1-exp(-theta)

## DATA
X = rmvnorm(n,mean=numeric(k),sigma=GammaCorrmat)
pX = pnorm(X)
Y=matrix(NA,n,k) ; for (i in 1:k) Y[,i] = (-1/gamma[i] * log(1+(log(1-pX[,i]))/theta[,i]))^(1/lambda[i])
Y[pX>upper] = 1e10 


accrual = matrix(runif(n*k,0,1),nrow=n)
studytime = 0.001 
c = matrix(rexp(n*k,studytime),nrow=n) # censored time
delta = 1*(c>=Y) # observation indicator
YY = pmin(c,Y) ; # observed survival time

cat("Event proportion: ",mean(delta),"\n")
# When delta==1, N cannot be zero. always N>0

### find the total elapsed time when the target event number has been reached (this may not happen)
if (colSums(delta)[1]>=(targetEvents)) {
  tottime = accrual + YY
  studyendpoint = sort(tottime[delta[,1]==1,1])[targetEvents]
} else {
  studyendpoint = 3
}

YY[YY>studyendpoint] = studyendpoint
delta = 1*(YY>=Y)

cat("Event proportion: ",mean(delta),"\n")

dat        = data.frame(Y=YY,delta,Z) ; # Final Data
names(dat) = c(paste0(rep(c("Y","delta"),each=k),1:k),paste0("Z",rep(1:k,each=p),1:p))
true=list(lambdat=lambda,betat=beta,thetat=gamma,rt=corr_r) ; dat = as.matrix(dat)

## data
n = nrow(dat) ; k = sum(colnames(dat)%like%"Y") ; p = sum(colnames(dat) %like%"Z")/k
Y = matrix(dat[,colnames(dat) %like%"Y"],nc=k) ; delta = matrix(dat[,colnames(dat) %like%"delta"],nc=k) ; Z = matrix(dat[,colnames(dat) %like%"Z"],nc=p*k) ;  

## hyper-parameters
beta0 = numeric(length(beta)) ; betaSig0 = solve((diag(length(beta)))*100)
lambda_s0 = theta_s0 = 1; lambda_r0 = theta_r0 = 0.1 ; corr_var0 = 0.7*0.7
priors = list(beta0=beta0,betaSig0=betaSig0, lambda_s0=lambda_s0,lambda_r0=lambda_r0, theta_s0=theta_s0,theta_r0=theta_r0, corr_var0=corr_var0)
## initial values
beta = numeric(ncol(Z)) ; lambda = theta = rep(1,k) ; corr_r = numeric(k*(k-1)/2);

invisible(ifelse(k==1,GammaCorrmat <- matrix(1,1,1),GammaCorrmat <- Corrmat(k,corr_r)$Gamma))
inits = list(beta=beta, lambda=lambda, theta=theta, corr_r=corr_r, GammaCorrmat=GammaCorrmat)

### SLICE SAMPLING PARAMETER ###
wlong = ifelse(sd(Y[,1], na.rm=TRUE) > 1,1,sd(Y[,1], na.rm=TRUE))
### Poisson integration related parameters
zz = 3

MM = 100000 ; burnin = 5000 ; thin = 1 ; saving = 10000 ; MHmcmc=0 ; accepttol = 0.25

start = Sys.time()

setwd('./results')
postsampling = MCMC_curerate(priors=priors, Y=Y,X=Z, delta=delta, dat=dat, inits=inits, slice_param=list(m=1,w=wlong), M=MM, burnin=burnin,  MHmcmc=MHmcmc, thin=thin, saving=saving, accepttol=accepttol, filename=filename)

end = Sys.time()

end - start

postsampling$time = end-start ; postsampling$true = true

saveRDS(postsampling, file=paste0(format(Sys.time(),"%y%m%d_%H%M%S"),"-",filename,".RData"))

