######## Bioequivalence tests based on estimating NLMEMs ########
### Author: Kathrin Möllenhoff
### code for applying the tests to two example datasets with parallel design
### N patients outcome is recorded at n points in time 
### Test is performed on the log-parameters (log AUC,log Cmax) 
### the threshold is chosen according to the 80/125-rule (log 1.25)
### example dataset 1 was simulated with rich design, high BSV and under H0
### example dataset 2 was simulated with sparse design, low BSV and under H1


library(saemix)
library(VGAM)
library(car)

alpha <- 0.05 #significance level of the test
epsilon <- log(1.25) #threshold of the test


##################################################
#Definition of the models used for the analysis###
##################################################
nb_chains=10
nb_iter=c(300, 100)
seed=23456
###defining the model for saemix / cf page 4 doc package saemix
model1cpt<-function(psi,id,xidep) {
  dose<-xidep[,1]
  tim<-xidep[,2]
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
#concentration in mg/L given by a one-compartment model
#model parameter given by (ka,V,Cl)
pk_modell <- function(D,v){
  Ka <- v[1]
  V <- v[2]
  Cl <- v[3]
  f=function(x){D*Ka/(Cl-Ka*V)*(exp(-Ka*x)-exp((-Cl/V)*x))}
  return(f)
}

#converting to (log-)parameters
#concentration in mg/L given by a one-compartment model
#model parameter given by (lK,lKa,lCl)
log_pk_modell <- function(D,v){
  lK <- v[1]
  lKa <- v[2]
  lCl <- v[3]
  f=function(x){(D*exp(lKa+lK-lCl)/(exp(lKa)-exp(lK)))*(exp(-exp(lK)*x)-exp(-exp(lKa)*x))}
  return(f)
}

# functions giving the AUC and the Cmax on the log scale, evaluated at the log-parameters 
#(log k,log ka,log Cl) and the corresponding true value of the test statistic
AUC_log <- function(D,v) {log(D/exp(v[3]))} #function evaluating the log(AUC)=log(D/Cl)
tmax <- function(v) {(v[2]-v[1])/(exp(v[2])-exp(v[1]))}
Cmax_log <- function(D,v) {log(D/(exp(v[3]-v[1]))*exp(-exp(v[1])*tmax(v)))} #function evaluating the log(Cmax)

#(help)function for calculating the euclidean norm of a vector
norm_vec <- function(x) sqrt(sum(x^2))

i <- 1 # or 2! number of dataset

  ####Loading of dataset
  tab <- read.table(paste(paste("dataset",i, sep="_"),".txt",sep=""),header=T,sep=" ",dec=".",na=".")
  assign(paste("data", i, sep="_"),tab)
  ####Appropriate format to fit with saemix
  tab <- subset(tab, is.na(tab$Dose))
  tab$Dose=NULL
  tab$Dose=4
  tab <- tab[,c(1,5,2,3,4)]
  
  ####Create treatment indicator: 0 for Reference treatment, 1 for Test treatment
  for (j in 1:nrow(tab)){
    if(tab$Tr[j]=="R"){tab$Treat[j]<-0} else{tab$Treat[j]<-1}
  }
  tab$Tr=NULL
  tab<-tab[,c(1,2,3,5,4)]
  tab$Treat <- as.numeric(tab$Treat)
  
  ####Number of sampling times
  nb_t<-sum(tab$Id[tab$Id==1])
  ####Number of subjects
  n<-nrow(tab)/nb_t
  
  #preparing everything for using saemix
  saemix.data<- saemixData(name.data=tab,header=TRUE,sep=" ",na=NA,
                           name.group=c("Id"),name.predictors=c("Dose","Time"),name.covariates = c("Treat"),name.response=c("Concentration"),
                           name.X="Time",units=list(x="hr",y="mg/L"))
  
  saemix.modelb<-saemixModel(model=model1cpt,
                             description="One-compartment model with first-order absorption",
                             psi0=matrix(c(1.5,0.5,0.04,0,0,0),ncol=3,nrow=2, byrow=TRUE,
                                         dimnames=list(NULL, c("ka","V","CL"))),
                             transform.par=c(1,1,1), 
                             covariate.model=matrix(c(1,1,1),ncol=3), # fixed.estim=c(1,1,1),
                             covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
                             omega.init=matrix(c(0.05,0,0,0,0.0125,0,0,0,0.05),ncol=3,byrow=TRUE),
                             error.model="combined",error.init=c(0.1,0.1))
  
  ####fitting the model (parameters -> ka, Cl, V)
  objet=paste("fit", i, sep="_")
  mod1 <- saemix(saemix.modelb,saemix.data,list(seed=seed,nb.chains=nb_chains,nbiter.saemix = nb_iter))
  assign(objet,mod1)
  
  Ka_estim <- mod1@results@fixed.effects[1]
  beta_Ka_estim<- mod1@results@fixed.effects[2]
  V_estim <- mod1@results@fixed.effects[3]
  beta_V_estim <- mod1@results@fixed.effects[4]
  Cl_estim <- mod1@results@fixed.effects[5]
  beta_Cl_estim <- mod1@results@fixed.effects[6]
  
  a_estim <- mod1@results@respar[1]
  b_estim <- mod1@results@respar[2]
  
  se_beta_Cl <- mod1@results@se.fixed[6]
  

  # Test on AUC
  
  ##############################################
  ##########Calculating the test statistic######
  ##############################################
  
  t.stat <- abs(beta_Cl_estim)
  se.beta <- se_beta_Cl
  
  ###################################
  #####new Test for Equivalence######
  ###################################
  
  BOT_stat <- t.stat
  BOT_limit <- qfoldnorm(alpha,epsilon,se.beta)
  BOT_rejectH0 <- ifelse(BOT_stat <= BOT_limit,1,0)
  #decision
  if(BOT_rejectH0==1){print("AUC: Reject H0")}else{print("AUC: Keep H0")}
  
  ###################################
  ###########model-based TOST#######
  ###################################
  
  stat_minusepsilon <- (beta_Cl_estim+epsilon)/se.beta
  stat_plusepsilon <- (beta_Cl_estim-epsilon)/se.beta
  #decision
  TOST_rejectH0_minus <- ifelse((stat_minusepsilon> qnorm(1-alpha)),1,0)
  TOST_rejectH0_plus <- ifelse((stat_plusepsilon< -qnorm(1-alpha)),1,0)
  TOST_rejectH0 <- ifelse((stat_minusepsilon>qnorm(1-alpha)) & (stat_plusepsilon< -qnorm(1-alpha)),1,0)
  if(TOST_rejectH0==1){print("AUC-TOST: Reject H0")}else{print("AUC-TOST: Keep H0")}
  
  # Test on Cmax
  
  ##############################################
  ##########Calculating the test statistic######
  ##############################################
  
  
    ####Delta method to compute beta_Cmax and to get asymptotic SE(beta_Cmax)
    estimates <- c(mod1@results@fixed.effects[1],mod1@results@fixed.effects[2],mod1@results@fixed.effects[3],mod1@results@fixed.effects[4],mod1@results@fixed.effects[5],mod1@results@fixed.effects[6])
    names(estimates) <- c("Ka","beta_Ka","V","beta_V","Cl","beta_Cl")
    FIM <- (mod1@results@fim)
    fim_fixed <- FIM[1:6,1:6]
    colnames(fim_fixed) <- c("Ka","beta_Ka","V","beta_V","Cl","beta_Cl")
    rownames(fim_fixed) <- c("Ka","beta_Ka","V","beta_V","Cl","beta_Cl")
    varcov<-solve(FIM)[1:6,1:6] #we only include the columns of the fix effects
    
    colnames(varcov) <- c("Ka","beta_Ka","V","beta_V","Cl","beta_Cl")
    rownames(varcov) <- c("Ka","beta_Ka","V","beta_V","Cl","beta_Cl")
    
    SE_deltam <- deltaMethod(estimates,"-beta_V-(log((Ka*V)/Cl)+beta_Ka+beta_V-beta_Cl)*Cl*exp(beta_Cl)/(Ka*V*exp(beta_Ka+beta_V)-Cl*exp(beta_Cl))+(Cl/(Ka*V-Cl))*log(Ka*V/Cl)",vcov.=varcov)
    
    beta_cmax_estim <- SE_deltam[,1]
    secmax <- SE_deltam[,2]
    t.stat <- abs(beta_cmax_estim)
    se.beta <- secmax
  
  ###################################
  ######new Test for Equivalence#####
  ###################################
  
 BOT_stat <- t.stat
 BOT_limit <- qfoldnorm(alpha,epsilon,se.beta)
 BOT_rejectH0 <- ifelse(BOT_stat <= BOT_limit,1,0)
 if(BOT_rejectH0==1){print("Cmax: Reject H0")}else{print("Cmax: Keep H0")}

 ###################################
 ###########model-based TOST#######
 ###################################
 
 stat_minusepsilon <- (beta_cmax_estim+epsilon)/se.beta
 stat_plusepsilon <- (beta_cmax_estim-epsilon)/se.beta
 #decision
 TOST_rejectH0_minus <- ifelse((stat_minusepsilon> qnorm(1-alpha)),1,0)
 TOST_rejectH0_plus <- ifelse((stat_plusepsilon< -qnorm(1-alpha)),1,0)
 TOST_rejectH0 <- ifelse((stat_minusepsilon>qnorm(1-alpha)) & (stat_plusepsilon< -qnorm(1-alpha)),1,0)
 if(TOST_rejectH0==1){print("Cmax-TOST: Reject H0")}else{print("Cmax-TOST: Keep H0")}
  


