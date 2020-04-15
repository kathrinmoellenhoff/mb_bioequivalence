######## Bioequivalence tests based on NCA ########
### Author: Kathrin Möllenhoff
### code for applying the tests to two example datasets with parallel design
### N patients outcome is recorded at n points in time 
### Test is performed on the log-parameters (log AUC,log Cmax) 
### the threshold is chosen according to the 80/125-rule (log 1.25)
### N=40 subjects, 20 for reference, 20 for test treatment
### example dataset 1 was simulated with rich design, high BSV and under H0
### example dataset 2 was simulated with sparse design, low BSV and under H1

library(VGAM)
library(MESS) 

alpha <- 0.05
epsilon <- log(1.25)

###########################################
###########################################

i <- 1 # or 2! number of dataset

  ####Loading of dataset
  tab<-read.table(paste(paste("dataset",i, sep="_"),".txt",sep=""),header=T,sep=" ",dec=".",na=".")
  assign(paste("data", i, sep="_"),tab)
  tab$Dose=NULL
  tab$Dose=4
  
  ####Create treatment indicator: 0 for Reference treatment, 1 for Test treatment
  for (j in 1:nrow(tab)){
    if(tab$Tr[j]=="R"){tab$Treat[j]<-0} else{tab$Treat[j]<-1}
  }
  tab$Tr=NULL
  tab$Treat <- as.numeric(tab$Treat)
  
  ####Number of sampling times
  nb_t<-sum(tab$Id[tab$Id==1])
  ####Number of subjects
  n<-nrow(tab)/nb_t
  
  time1<-unique(tab$Time)
  id_vect<-c(unique(tab$Id))
  
  auc_ref_vect<-c()
  log_auc_ref_vect<-c()
  cmax_ref_vect<-c()
  log_cmax_ref_vect<-c()
  auc_treat_vect<-c()
  log_auc_treat_vect<-c()
  cmax_treat_vect<-c()
  log_cmax_treat_vect<-c()
  
  #calculate NCA estimates#
  
  for (j in 1:(n/2)){
    # Id
    tabbis<-tab[tab$Id==j,]
    
    #ref group
    conc1<-tabbis$Concentration[tabbis$Treat==0]
    auc_ref<-auc(time1, conc1, type='linear')
    
    log_auc_ref=log(auc_ref)
    auc_ref_vect<-c(auc_ref_vect,auc_ref)
    log_auc_ref_vect<-c(log_auc_ref_vect,log_auc_ref)
    
    cmax_ref<-max(conc1)
    cmax_ref_vect<-c(cmax_ref_vect,cmax_ref)
    log_cmax_ref<-log(cmax_ref)
    log_cmax_ref_vect<-c(log_cmax_ref_vect,log_cmax_ref)
    
    #treat group
    tabbis<-tab[tab$Id==(20+j),]
    conc2<-tabbis$Concentration[tabbis$Treat==1]
    auc_treat<-auc(time1, conc2, type='linear')
    
    log_auc_treat=log(auc_treat)
    auc_treat_vect<-c(auc_treat_vect,auc_treat)
    log_auc_treat_vect<-c(log_auc_treat_vect,log_auc_treat)
    
    cmax_treat<-max(conc2)
    cmax_treat_vect<-c(cmax_treat_vect,cmax_treat)
    log_cmax_treat<-log(cmax_treat)
    log_cmax_treat_vect<-c(log_cmax_treat_vect,log_cmax_treat)
  }
  
  AUC<-c(auc_ref_vect,auc_treat_vect)
  log_AUC<-c(log_auc_ref_vect,log_auc_treat_vect)
  Cmax<-c(cmax_ref_vect,cmax_treat_vect)
  log_Cmax<-c(log_cmax_ref_vect,log_cmax_treat_vect)
  
  
  sd_log_AUC<-sd(log_AUC[1:20]-log_AUC[21:40])/sqrt(20)
  sd_log_Cmax<-sd(log_Cmax[1:20]-log_Cmax[21:40])/sqrt(20)
  
  mean_diff_log_auc<-mean(log_AUC[1:20]-log_AUC[21:40])
  mean_diff_log_cmax<-mean(log_Cmax[1:20]-log_Cmax[21:40])
  
  ###################################
  #######Tests for Equivalence########
  ###################################
  
  #######################################
  ########NCA-based optimal BE test#####
  ######################################
  
  #AUC
  t.stat <- abs(mean_diff_log_auc)
  se.beta <- sd_log_AUC

  BOT_stat<-t.stat
  BOT_limit<-qfoldnorm(alpha,epsilon,se.beta)
  BOT_rejectH0<-ifelse(BOT_stat <= BOT_limit,1,0)
  #decision
  if(BOT_rejectH0==1){print("AUC: Reject H0")}else{print("AUC: Keep H0")}
  
  #same for Cmax
  t.stat <- abs(mean_diff_log_cmax)
  se.beta <- sd_log_Cmax
  
  BOT_stat<-t.stat
  BOT_limit<-qfoldnorm(alpha,epsilon,se.beta)
  BOT_rejectH0<-ifelse(BOT_stat <= BOT_limit,1,0)
  #decision
  if(BOT_rejectH0==1){print("Cmax: Reject H0")}else{print("Cmax: Keep H0")}
  
  
  ###################################
  ########### NCA-based TOST#######
  ###################################  
  
  #AUC
  stat_minusepsilon <- (mean_diff_log_auc+epsilon)/sd_log_AUC
  stat_plusepsilon <- (mean_diff_log_auc-epsilon)/sd_log_AUC
  #decision
  TOST_rejectH0_minus <- ifelse((stat_minusepsilon> qnorm(1-alpha)),1,0)
  TOST_rejectH0_plus <- ifelse((stat_plusepsilon< -qnorm(1-alpha)),1,0)
  TOST_rejectH0 <- ifelse((stat_minusepsilon>qnorm(1-alpha)) & (stat_plusepsilon< -qnorm(1-alpha)),1,0)
  if(TOST_rejectH0==1){print("AUC-TOST: Reject H0")}else{print("AUC-TOST: Keep H0")}
  
  #Cmax
  stat_minusepsilon <- (mean_diff_log_cmax+epsilon)/sd_log_Cmax
  stat_plusepsilon <- (mean_diff_log_cmax-epsilon)/sd_log_Cmax
  #decision
  TOST_rejectH0_minus <- ifelse((stat_minusepsilon> qnorm(1-alpha)),1,0)
  TOST_rejectH0_plus <- ifelse((stat_plusepsilon< -qnorm(1-alpha)),1,0)
  TOST_rejectH0 <- ifelse((stat_minusepsilon>qnorm(1-alpha)) & (stat_plusepsilon< -qnorm(1-alpha)),1,0)
  if(TOST_rejectH0==1){print("Cmax-TOST: Reject H0")}else{print("Cmax-TOST: Keep H0")}
  
  
