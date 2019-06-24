# 10.17: replace if statement in row 38 with an else statement, 
# because occasionally there's a dataset with both Lastfix and lastfix in it...?
# 1.17.17: NO. REVERTED BACK TO PREVIOUS DUE TO CALLING FN on sims dataset: fixed the *1000 issue in the fit function by dividing by 1000 in the fixsumsFN
# (but not the fixsumsFNmydata because my data is in seconds; other data is in ms)
# V3 has adaptation for 4 ROIs
fixsumsFN <- function(sims) {
  sims$leftfixtime<-(sims$roi%%2)*sims$fixdur
  sims$rightfixtime<-((sims$roi+1)%%2)*sims$fixdur
  
  #do a quick cumulative sum over those fixation durations over the whole dataset
  sims$tleft<-cumsum(sims$leftfixtime)
  sims$tright<-cumsum(sims$rightfixtime)
  
  #now only look at the last fixation of every trial
  temp<-sims[sims$revfixnum==1,]
  
  #subtract the sum at the end of previous trial from the sum at the end of the current trial - that gives you just the sum from the current trial
  templ<-c(0,temp$tleft[1:(nrow(temp)-1)])
  tempr<-c(0,temp$tright[1:(nrow(temp)-1)])
  temp$tleft<-temp$tleft-templ
  temp$tright<-temp$tright-tempr
  
  #calculate the final total left time advantage variable
  temp$leftadv<-temp$tleft-temp$tright
  head(temp)
  
  fixsums <- temp 
  return(fixsums)
}

fixsumsFNmydata <- function(sims) {
  sims$leftfixtime<-(sims$ROI%%2)*sims$Time
  sims$rightfixtime<-((sims$ROI+1)%%2)*sims$Time
  
  #do a quick cumulative sum over those fixation durations over the whole dataset
  sims$tleft<-cumsum(sims$leftfixtime)
  sims$tright<-cumsum(sims$rightfixtime)
  
  #now only look at the last fixation of every trial
  if (length(sims$Lastfix)>0) {
    temp<-sims[sims$Lastfix==1,]
  } else {
    temp<-sims[sims$lastfix==1,]
  }
  
  #subtract the sum at the end of previous trial from the sum at the end of the current trial - that gives you just the sume from the current trial
  templ<-c(0,temp$tleft[1:(nrow(temp)-1)])
  tempr<-c(0,temp$tright[1:(nrow(temp)-1)])
  temp$tleft<-temp$tleft-templ
  temp$tright<-temp$tright-tempr
  
  #calculate the final total left time advantage variable
  temp$leftadv<-temp$tleft-temp$tright
  head(temp)
  
  fixsums <- temp 
  fixsums$leftfixtime <- NULL
  fixsums$rightfixtime <- NULL
  fixsums$Lastfix <- NULL
  fixsums$Firstfix <- NULL
  fixsums$Midfix <- NULL
  colnames(fixsums)[(ncol(fixsums)-2):ncol(fixsums)] <- c("Fix1","Fix2","TAdvL")
  fixsums$tleft <- fixsums$Fix1
  fixsums$tright <- fixsums$Fix2
  return(fixsums)
}