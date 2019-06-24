# This code expects a data file that was the input to the "data_prep_for_addm" script and so has those required variables.
# It also expects a simulated dataset from the "sim_addm_fun" script and saved as "sims"  
# It produces some goodness of fit statistics called CS (chi-square) and ML (log likelihood)

# Note that there is some code below that does some specific trimming of the data based on absolute value difference. 
# Comment this out for your datasets if you don't need it.

# 10.10: changed correct percentage tadv sep to absolute dwell time tadv sep
# 11.17: fixed the zero valdiff section in sims... still didn't fix the data
# 12.1: fixed the midround issue with absvaldiff (was getting multiple valdiffs for the same number)


fit_addm_tadvsep_fun_v3 <- function(data,sims) {
  ####chi-square and MLE fit with correct vs. incorrect by valdiff and tadvsep
  pquants<-c(0,0.5,1) #desired correct and zero value-difference quantiles (INCLUDE 0 and 1)
  nquants<-c(0,0.5,1) #desired incorrect quantiles (INCLUDE 0 and 1)
  tadvseps <- c(0.333,0.667) #desired separations (DON'T INCLUDE 0 or 1)
  unbound<-1 #set to 1 if you want the lower bound to be zero and the upper bound to be 100secs.
  maxvaldiff<-100000 #use this to set a maximum value difference 
  fourmode<-3 #use this to set how you want to bin trials with 4 items. #1: only use the top items, #2:only use the bottom items, #3: use the average of the top and bottom
  #############################################
  
  ntadvseps <- length(tadvseps)+1 # number of time advantage "chunks" to slice the data into for binning
  CS<-0 # this is where we store the chi-square statistic
  ML<-0 # this is where we store the MLE statistic
  sims<-fixsumsFN(sims)
  
  
  #2 item case
  if (max(data$roi)<3){
    data$valdiff<-(data$leftval-data$rightval)
    data$absvaldiff<-midround(abs(data$valdiff),0.001)
    sims$valdiff<-(sims$leftval-sims$rightval)
    sims$absvaldiff<-midround(abs(sims$valdiff),0.001)
  }
  
  #4 item case
  if (max(data$roi)>2 & fourmode==1){
    data$valdiff<-data$upleftval-data$uprightval
    data$absvaldiff<-midround(abs(data$valdiff),0.001)
    sims$valdiff<-sims$upleftval-sims$uprightval
    sims$absvaldiff<-midround(abs(sims$valdiff),0.001)
  }
  if (max(data$roi)>2 & fourmode==2){
    data$valdiff<-data$downleftval-data$downrightval
    data$absvaldiff<-midround(abs(data$valdiff),0.001)
    sims$valdiff<-sims$downleftval-sims$downrightval
    sims$absvaldiff<-midround(abs(sims$valdiff),0.001)
  }
  if (max(data$roi)>2 & fourmode==3){
    data$valdiff<-(data$upleftval+data$downleftval-data$uprightval-data$downrightval)/2
    data$absvaldiff<-midround(abs(data$valdiff),0.001)
    sims$valdiff<-(sims$upleftval+sims$downleftval-sims$uprightval-sims$downrightval)/2
    sims$absvaldiff<-midround(abs(sims$valdiff),0.001)
  }
  
  #remove extreme value differences
  data<-data[data$absvaldiff<=maxvaldiff,]
  sims<-sims[sims$absvaldiff<=maxvaldiff,]
  
  # edited to help with missing 0 valdiff sims trials
  absvaldiffs<-sort(intersect(unique(data$absvaldiff),unique(sims$absvaldiff)))
  
  if (median(data$tleft) > 50) { # if median total dwell time left is over 50, we know it's in ms
    data$tadvl <- with(data, tleft-tright)
  } else { # if not, we need to convert to ms
    data$tadvl <- with(data, tleft-tright)*1000
  }
  data$tadvc <- with(data, ((valdiff) >= 0)*tadvl + ((valdiff) < 0)*(-1*tadvl))
  sims$tadvl <- with(sims, tleft-tright)
  sims$tadvc <- with(sims, ((valdiff) >= 0)*tadvl + ((valdiff) < 0)*(-1*tadvl))
  # here, all left choices in 0-valuedifference trials are coded as correct
  # all right choices in 0-valuedifference trials are coded as incorrect
  
  #preallocate matrices to store the quantiles and the counts/probabilities
  dataprobs<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=(ntadvseps*(length(pquants)+length(nquants)-2)))
  posquants<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=ntadvseps*length(pquants))
  negquants<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=ntadvseps*length(nquants))
  
  #check if there are zero value difference trials and if so, skip them in the loop
  start<-1
  if (absvaldiffs[1] == 0) {start<-2}
  # added min(sims$absvaldiff) == 0 
  # sometimes the sims are not including 0 valdiff trials (because they don't finish and get excluded?)
  
  postemp <- NULL
  negtemp <- NULL
  if (length(absvaldiffs)>1) {
    for (i in start:length(absvaldiffs)){
      postemp2 <- NULL
      negtemp2 <- NULL
      temp<-data[data$absvaldiff==absvaldiffs[i],]
      #seps <- c(-1000000,quantile(temp$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      posdat<-temp[((temp$valdiff>0) & temp$choice==1) | ((temp$valdiff<0) & temp$choice==0),]
      negdat<-temp[((temp$valdiff>0) & temp$choice==0) | ((temp$valdiff<0) & temp$choice==1),]
      posseps <- c(-1000000,quantile(posdat$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      negseps <- c(-1000000,quantile(negdat$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      for (sepnum in 1:ntadvseps) {
        posdattemptadv <- posdat[posdat$tadvc >= posseps[sepnum] & posdat$tadvc < posseps[sepnum+1],]
        negdattemptadv <- negdat[negdat$tadvc >= negseps[sepnum] & negdat$tadvc < negseps[sepnum+1],]
        posquants[(i-start+1),((sepnum-1)*length(pquants)+1):(sepnum*length(pquants))]<-quantile(posdattemptadv$rt,probs=pquants,na.rm=TRUE,names=FALSE)
        negquants[(i-start+1),((sepnum-1)*length(nquants)+1):(sepnum*length(nquants))]<-quantile(negdattemptadv$rt,probs=nquants,na.rm=TRUE,names=FALSE)
      }
      posquants[is.na(posquants)] <- 0
      negquants[is.na(negquants)] <- 0
      
      #set lower and upper bounds at 0 and 100s if "unbound" is turned on
      if (unbound==1){
        for (sepnum in 1:ntadvseps) {
          posquants[(i-start+1),((sepnum-1)*length(pquants)+1)]<-0
          posquants[(i-start+1),(sepnum*length(pquants))]<-10000000 #some really big number
          negquants[(i-start+1),((sepnum-1)*length(nquants)+1)]<-0
          negquants[(i-start+1),(sepnum*length(nquants))]<-10000000 #some really big number	
        }
      }
      
      for (sepnum in 1:ntadvseps) {
        posdattemptadv <- posdat[posdat$tadvc >= posseps[sepnum] & posdat$tadvc < posseps[sepnum+1],]
        negdattemptadv <- negdat[negdat$tadvc >= negseps[sepnum] & negdat$tadvc < negseps[sepnum+1],]
        posdatahist<-hist(posdattemptadv$rt,breaks=posquants[(i-start+1),((sepnum-1)*length(pquants)+1):(sepnum*length(pquants))],plot=FALSE)
        posdataprobs<-posdatahist$counts
        negdatahist<-hist(negdattemptadv$rt,breaks=negquants[(i-start+1),((sepnum-1)*length(nquants)+1):(sepnum*length(nquants))],plot=FALSE)
        negdataprobs<-negdatahist$counts
        dataprobs[(i-start+1),((sepnum-1)*(length(posdataprobs)+length(negdataprobs))+1):(sepnum*(length(posdataprobs)+length(negdataprobs)))]<-c(posdataprobs,negdataprobs)
        postemp2 <- c(postemp2, posdataprobs)
        negtemp2 <- c(negtemp2, negdataprobs)
      }
      postemp <- rbind(postemp, postemp2)
      negtemp <- rbind(negtemp, negtemp2)
    }
    
    #case where there is zero value difference
    if (start==2){
      temp <- data[data$absvaldiff==0,]
      zeroseps <- c(-1000000,quantile(temp$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      zeroquants<-matrix(0,nrow=1,ncol=ntadvseps*length(pquants))
      for (sepnum in 1:ntadvseps) {
        temptadv<-temp[temp$tadvc >= zeroseps[sepnum] & temp$tadvc < zeroseps[sepnum+1],]
        zeroquants[((sepnum-1)*length(pquants)+1):(sepnum*length(pquants))]<-quantile(temptadv$rt,probs=pquants,na.rm=TRUE,names=FALSE)
      }
      zeroquants[is.na(zeroquants)] <- 0
      
      #set lower and upper bounds at 0 and 100s if "unbound" is turned on
      if (unbound==1){
        for (sepnum in 1:ntadvseps) {
          zeroquants[(sepnum-1)*length(pquants)+1]<-0
          zeroquants[sepnum*length(pquants)]<-100000 #some really big number
        }
      }
      
      zerodataprobs<-matrix(0,nrow=1,ncol=ntadvseps*(length(pquants)-1))
      for (sepnum in 1:ntadvseps) {
        temptadv<-temp[temp$tadvc >= zeroseps[sepnum] & temp$tadvc < zeroseps[sepnum+1],]
        zerodatahist<-hist(temptadv$rt,breaks=zeroquants[((sepnum-1)*length(pquants)+1):(sepnum*length(pquants))],plot=FALSE)
        zerodataprobs[((sepnum-1)*length(zerodatahist$counts)+1):(sepnum*(length(zerodatahist$counts)))]<-zerodatahist$counts
      }
    }
    
    ###now handle the simulations
    simsprobs<-matrix(0,nrow=length(absvaldiffs[absvaldiffs!=0]),ncol=(ntadvseps*(length(pquants)+length(nquants)-2)))
    
    for (i in start:length(absvaldiffs)){
      # use the same seps from the data, not from the sims
      tempdat<-data[data$absvaldiff==absvaldiffs[i],]
      posdat<-tempdat[((tempdat$valdiff>0) & tempdat$choice==1) | ((tempdat$valdiff<0) & tempdat$choice==0),]
      negdat<-tempdat[((tempdat$valdiff>0) & tempdat$choice==0) | ((tempdat$valdiff<0) & tempdat$choice==1),]
      posseps <- c(-1000000,quantile(posdat$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      negseps <- c(-1000000,quantile(negdat$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      #seps <- c(-1000000,quantile(tempdat$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      temp<-sims[sims$absvaldiff==absvaldiffs[i],]
      possim<-temp[((temp$valdiff>0) & temp$choice==1) | ((temp$valdiff<0) & temp$choice==0),]
      negsim<-temp[((temp$valdiff>0) & temp$choice==0) | ((temp$valdiff<0) & temp$choice==1),]
      #if (nrow(possim) > 0 & nrow(negsim) > 0) {
      for (sepnum in 1:ntadvseps) {
        possimtadv <- possim[possim$tadvc >= posseps[sepnum] & possim$tadvc < posseps[sepnum+1],]
        negsimtadv <- negsim[negsim$tadvc >= negseps[sepnum] & negsim$tadvc < negseps[sepnum+1],]
        posrthist<-hist(possimtadv$rt,breaks=posquants[(i-start+1),((sepnum-1)*length(pquants)+1):(sepnum*length(pquants))],plot=FALSE)
        posrtprobs<-posrthist$counts
        negrthist<-hist(negsimtadv$rt,breaks=negquants[(i-start+1),((sepnum-1)*length(nquants)+1):(sepnum*length(nquants))],plot=FALSE)
        negrtprobs<-negrthist$counts
        simsprobs[(i-start+1),((sepnum-1)*(length(posrtprobs)+length(negrtprobs))+1):(sepnum*(length(posrtprobs)+length(negrtprobs)))]<-c(posrtprobs,negrtprobs)
      }
      # Make sure at least 1 in each bin
      tempCS<-ifelse((any(simsprobs[(i-start+1),] > 0) & 
                        sum(postemp[(i-start+1),])>(ntadvseps*(length(pquants)-1)) & 
                        sum(negtemp[(i-start+1),])>(ntadvseps*(length(nquants)-1))),
                     chisq.test(dataprobs[(i-start+1),],p=simsprobs[(i-start+1),],rescale.p=TRUE)$statistic,0)
      tempML<-ifelse((any(simsprobs[(i-start+1),] > 0) & 
                        sum(postemp[(i-start+1),])>(ntadvseps*(length(pquants)-1)) & 
                        sum(negtemp[(i-start+1),])>(ntadvseps*(length(nquants)-1))),
                     sum(dataprobs[(i-start+1),]*log(simsprobs[(i-start+1),]/(sum(simsprobs[(i-start+1),])))),0)
      if (is.na(tempCS)) (tempCS <- 0)
      if (is.na(tempML)) (tempML <- 0)
      CS<-CS+tempCS
      ML<-ML+tempML
    }
    
    #case where there is zero value difference
    if (start==2){
      tempdat<-data[data$absvaldiff==0,]
      zeroseps <- c(-1000000,quantile(tempdat$tadvc,probs = tadvseps,na.rm = TRUE,names = FALSE),1000000)
      zerosimsprobs<-matrix(0,nrow=1,ncol=(ntadvseps*(length(pquants)-1)))
      temp<-sims[sims$absvaldiff==0,]	
      for (sepnum in 1:ntadvseps) {
        temptadv <- temp[temp$tadvc >= zeroseps[sepnum] & temp$tadvc < zeroseps[sepnum+1],]
        zerorthist<-hist(temptadv$rt,breaks=zeroquants[((sepnum-1)*length(pquants)+1):(sepnum*length(pquants))],plot=FALSE)
        zerosimsprobs[((sepnum-1)*(length(zerorthist$counts))+1):(sepnum*(length(zerorthist$counts)))]<-zerorthist$counts
      }
      # Make sure at least one in each bin!
      tempCS<-ifelse(any(zerosimsprobs > 0) & 
                       sum(zerodataprobs)>(ntadvseps*(length(pquants)-1)), 
                     chisq.test(zerodataprobs,p=zerosimsprobs,rescale.p=TRUE)$statistic,0) 
      tempML<-ifelse(any(zerosimsprobs > 0) & 
                       sum(zerodataprobs)>(ntadvseps*(length(pquants)-1)), 
                     sum(zerodataprobs*log(zerosimsprobs/(sum(zerosimsprobs)))),0)
      if (is.na(tempCS)) (tempCS <- 0)
      if (is.na(tempML)) (tempML <- 0)
      CS<-CS+tempCS
      ML<-ML+tempML
      
    }	
    gof<-c(CS,ML)
  } else {
    gof <- c(Inf,-Inf)
  }
  return(gof)
}