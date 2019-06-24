#this code expects a data file that was the input to the "data_prep_for_addm" script and so has those required variables; it also expects a simulated dataset from the "sim_addm_fun" script and saved as "sims".  It produces some goodness of fit statistics called CS (chi-square) and ML (log likelihood)
#(ADDM)
#note that there is some code below that does some specific trimming of the data based on absolute value difference. comment this out for your datasets.


fit_addm_fun_v5<-function(data,sims){
  ####chi-square and MLE fit with correct vs. incorrect by valdiff
  pquantsp<-c(0,0.1,0.3,0.5,0.7,0.9,1) #correct choice, last gaze to correct item quantiles
  pquantsn<-c(0,0.1,0.3,0.5,0.7,0.9,1) #correct choice, last gaze to incorrect item quantiles
  nquantsn<-c(0,0.5,1) #incorrect choice, last gaze to incorrect item quantiles
  nquantsp<-c(0,0.5,1) #incorrect choice, last gaze to correct item quantiles
  zquants<-c(0,0.1,0.3,0.5,0.7,0.9,1) #zero quantiles 
  unbound<-1 #set to 1 if you want the lower bound to be zero and the upper bound to be 100secs.
  maxvaldiff<-5 #use this to set a maximum value difference 
  roundnum <-0.1 # round the val diffs to the nearest x (can be set to any decimal or integer)
  #############################################
  
  CS<-0#this is where we store the chi-square statistic
  ML<-0#this is where we store the MLE statistic
  if (length(data$revfixnum)>0) {data<-data[data$revfixnum==1,]}#only want one row per trial
  if (length(data$Lastfix)>0) {data<-data[data$Lastfix==1,]}#only want one row per trial
  if (length(data$lastfix)>0) {data<-data[data$lastfix==1,]}#only want one row per trial
  sims<-sims[sims$revfixnum==1,] #only want one row per trial
  
  # midround function
  midround <- function(x,base){ 
    base*round(x/base) 
  }
  #define value difference and absolute value difference
  data$valdiff<-midround(data$leftval-data$rightval,roundnum) 
  data$absvaldiff<-abs(data$valdiff) 
  sims$valdiff<-midround(sims$leftval-sims$rightval,roundnum) 
  sims$absvaldiff<-abs(sims$valdiff) 
  
  #remove extreme value differences
  data<-data[data$absvaldiff<=maxvaldiff,]
  sims<-sims[sims$absvaldiff<=maxvaldiff,]
  
  absvaldiffs<-sort(unique(data$absvaldiff))
  
  #preallocate matrices to store the quantiles and the counts/probabilities
  dataprobs<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=((length(pquantsp)+length(pquantsn)+length(nquantsn)+length(nquantsp)-4)))
  posquantsp<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=length(pquantsp))
  posquantsn<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=length(pquantsn))
  negquantsn<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=length(nquantsn))
  negquantsp<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=length(nquantsp))
  
  #check if there are zero value difference trials and if so, skip them in the loop
  start<-1
  if (absvaldiffs[1]==0) (start<-2)
  
  for (i in start:length(absvaldiffs)){
    temp<-data[data$absvaldiff==absvaldiffs[i],]
    posdatp<-temp[((temp$valdiff>0) & temp$choice==1 & temp$roi==1) | ((temp$valdiff<0) & temp$choice==0 & temp$roi==2),]
    posdatn<-temp[((temp$valdiff>0) & temp$choice==1 & temp$roi==2) | ((temp$valdiff<0) & temp$choice==0 & temp$roi==1),]
    negdatn<-temp[((temp$valdiff>0) & temp$choice==0 & temp$roi==2) | ((temp$valdiff<0) & temp$choice==1 & temp$roi==1),]
    negdatp<-temp[((temp$valdiff>0) & temp$choice==0 & temp$roi==1) | ((temp$valdiff<0) & temp$choice==1 & temp$roi==2),]
    posquantsp[(i-start+1),]<-quantile(posdatp$rt,probs=pquantsp,na.rm=TRUE,names=FALSE)
    posquantsn[(i-start+1),]<-quantile(posdatn$rt,probs=pquantsn,na.rm=TRUE,names=FALSE)
    negquantsn[(i-start+1),]<-quantile(negdatn$rt,probs=nquantsn,na.rm=TRUE,names=FALSE)
    negquantsp[(i-start+1),]<-quantile(negdatp$rt,probs=nquantsp,na.rm=TRUE,names=FALSE)
    
    #set lower and upper bounds at 0 and 100s if "unbound" is turned on
    if (unbound==1){
      posquantsp[(i-start+1),1]<-0
      posquantsp[(i-start+1),(length(pquantsp))]<-100000 #some really big number
      posquantsn[(i-start+1),1]<-0
      posquantsn[(i-start+1),(length(pquantsn))]<-100000 #some really big number
      negquantsn[(i-start+1),1]<-0
      negquantsn[(i-start+1),(length(nquantsn))]<-100000 #some really big number
      negquantsp[(i-start+1),1]<-0
      negquantsp[(i-start+1),(length(nquantsp))]<-100000 #some really big number		
    }
    
    posdatahistp<-hist(posdatp$rt,breaks=posquantsp[(i-start+1),],plot=FALSE)
    posdataprobsp<-posdatahistp$counts
    posdatahistn<-hist(posdatn$rt,breaks=posquantsn[(i-start+1),],plot=FALSE)
    posdataprobsn<-posdatahistn$counts
    negdatahistn<-hist(negdatn$rt,breaks=negquantsn[(i-start+1),],plot=FALSE)
    negdataprobsn<-negdatahistn$counts
    negdatahistp<-hist(negdatp$rt,breaks=negquantsp[(i-start+1),],plot=FALSE)
    negdataprobsp<-negdatahistp$counts
    dataprobs[(i-start+1),]<-c(posdataprobsp,posdataprobsn,negdataprobsn,negdataprobsp)
  }
  
  #is this right?
  posquantsn[is.na(posquantsn)] <- 0
  posquantsp[is.na(posquantsp)] <- 0
  negquantsn[is.na(negquantsn)] <- 0
  negquantsp[is.na(negquantsp)] <- 0
  
  #case where there is zero value difference
  if (start==2){
    temp<-data[data$absvaldiff==0,]
    ldatl<-temp[(temp$choice==1 & temp$roi==1) | (temp$choice==0 & temp$roi==2),]
    ldatr<-temp[(temp$choice==1 & temp$roi==2) | (temp$choice==0 & temp$roi==1),]
    lzeroquantsl<-quantile(ldatl$rt,probs=zquants,na.rm=TRUE,names=FALSE)
    lzeroquantsr<-quantile(ldatr$rt,probs=zquants,na.rm=TRUE,names=FALSE)
    
    #set lower and upper bounds at 0 and 100s if "unbound" is turned on
    if (unbound==1){
      lzeroquantsl[1]<-0
      lzeroquantsl[length(zquants)]<-100000 #some really big number
      lzeroquantsr[1]<-0
      lzeroquantsr[length(zquants)]<-100000 #some really big number
    }
    
    lzerodatahistl<-hist(ldatl$rt,breaks=lzeroquantsl,plot=FALSE)
    lzerodatahistr<-hist(ldatr$rt,breaks=lzeroquantsr,plot=FALSE)
    zerodataprobs<-c(lzerodatahistl$counts,lzerodatahistr$counts)
  }
  lzeroquantsl[is.na(lzeroquantsl)] <- 0
  lzeroquantsr[is.na(lzeroquantsr)] <- 0
  
  ###now handle the simulations
  simsprobs<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=((length(pquantsp)+length(pquantsn)+length(nquantsn)+length(nquantsp)-4)))
  
  for (i in start:length(absvaldiffs)){
    temp<-sims[sims$absvaldiff==absvaldiffs[i],]	
    possimp<-temp[((temp$valdiff>0) & temp$choice==1 & temp$roi==1) | ((temp$valdiff<0) & temp$choice==0 & temp$roi==2),]
    posrthistp<-hist(possimp$rt,breaks=posquantsp[(i-start+1),],plot=FALSE)
    posrtprobsp<-posrthistp$counts
    possimn<-temp[((temp$valdiff>0) & temp$choice==1 & temp$roi==2) | ((temp$valdiff<0) & temp$choice==0 & temp$roi==1),]
    posrthistn<-hist(possimn$rt,breaks=posquantsn[(i-start+1),],plot=FALSE)
    posrtprobsn<-posrthistn$counts
    negsimn<-temp[((temp$valdiff>0) & temp$choice==0 & temp$roi==2) | ((temp$valdiff<0) & temp$choice==1 & temp$roi==1),]
    negrthistn<-hist(negsimn$rt,breaks=negquantsn[(i-start+1),],plot=FALSE)
    negrtprobsn<-negrthistn$counts
    negsimp<-temp[((temp$valdiff>0) & temp$choice==0 & temp$roi==1) | ((temp$valdiff<0) & temp$choice==1 & temp$roi==2),]
    negrthistp<-hist(negsimp$rt,breaks=negquantsp[(i-start+1),],plot=FALSE)
    negrtprobsp<-negrthistp$counts
    
    simsprobs[(i-start+1),]<-c(posrtprobsp,posrtprobsn,negrtprobsn,negrtprobsp)
    
    simsprobs[dataprobs == 0] <- 1
    
    
    #tempCS<-chisq.test(dataprobs[(i-start+1),],p=simsprobs[(i-start+1),],rescale.p=TRUE)$statistic
    tempML<-sum(dataprobs[(i-start+1),]*log(simsprobs[(i-start+1),]/(sum(simsprobs[(i-start+1),]))))
    #if (is.na(tempCS)) (tempCS <- 0)
    if (is.na(tempML)) (tempML <- 0)
    #CS<-CS+tempCS
    ML<-ML+tempML
  }
  
  #case where there is zero value difference
  if (start==2){
    temp<-sims[sims$absvaldiff==0,]	
    lsiml<-temp[(temp$choice==1 & temp$roi==1) | (temp$choice==0 & temp$roi==2),]
    lsimr<-temp[(temp$choice==1 & temp$roi==2) | (temp$choice==0 & temp$roi==1),]
    lzerorthistl<-hist(lsiml$rt,breaks=lzeroquantsl,plot=FALSE)
    lzerorthistr<-hist(lsimr$rt,breaks=lzeroquantsr,plot=FALSE)
    zerosimsprobs<-c(lzerorthistl$counts,lzerorthistr$counts)
    
    zerosimsprobs[zerodataprobs == 0] <- 1
    
    #tempCS<-chisq.test(zerodataprobs,p=zerosimsprobs,rescale.p=TRUE)$statistic
    tempML<-sum(zerodataprobs*log(zerosimsprobs/(sum(zerosimsprobs))))
    #if (is.na(tempCS)) (tempCS <- 0)
    if (is.na(tempML)) (tempML <- 0)
    #CS<-CS+tempCS
    ML<-ML+tempML
  }	
  gof<-c(CS,ML)
  return(gof)
}







