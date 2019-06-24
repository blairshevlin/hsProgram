#this code expects a data file that has already been run through the "data_prep_for_addm" script and saved as "data_mat"; nsim is number of simulations, d is drift parameter, sigma is within trial variability in drift, theta is attention parameter, ndt is non-decision time in ms, dsigma is across trial variability in drift. It produces a simulated dataset called "sims".
sim_addm_fun_v2<-function(nsim,d,sigma,theta,ndt,dsigma,data_mat){
	# random number generator library
  	library(RcppZiggurat)
   
	ntrial<-nrow(data_mat) #number of unique trials in the dataset (get this from the data)
	maxfix<-(ncol(data_mat)-4)/2 #maximum number of fixations in a given trial
	
	#preallocate the simulation vectors
	rt<-rep(0,nsim*maxfix*ntrial) #reaction time vector
	choice<-rep(-1,nsim*maxfix*ntrial) #choice vector
	fixnum<-rep(-1,nsim*maxfix*ntrial) #fixation number vector
	revfixnum<-rep(-1,nsim*maxfix*ntrial) #reverse fixation number vector
	fixdur<-rep(-1,nsim*maxfix*ntrial) #fixation duration vector
	nfix<-rep(-1,nsim*maxfix*ntrial) #total number of fixations vector
	roi<-rep(-1,nsim*maxfix*ntrial) #roi vector
	leftval<-rep(-1,nsim*maxfix*ntrial) #left value vector
	rightval<-rep(-1,nsim*maxfix*ntrial) #right value vector
	
	roi_mat<-data_mat[,5:(4+maxfix)]
	fixlength_mat<-data_mat[,(5+maxfix):ncol(data_mat)]
	value_mat<-(roi_mat==1)*data_mat[,3]+(roi_mat==2)*data_mat[,4]
	lvalue_mat<-data_mat[,3]*matrix(1,nrow=nrow(data_mat),ncol=maxfix)
	rvalue_mat<-data_mat[,4]*matrix(1,nrow=nrow(data_mat),ncol=maxfix)
	drift_mat<-d*((roi_mat==1)*(lvalue_mat-theta*rvalue_mat)+(roi_mat==2)*(theta*lvalue_mat-rvalue_mat))
	
	counter<-1
	#loop over unique trials
	for (k in 1:ntrial){
		nroi<-roi_mat[k,]
		nfixlength<-fixlength_mat[k,]
		ndrift<-drift_mat[k,]
		ntemp = cumsum(nfixlength)  
		index = rep(0,max(ntemp))
		index[ntemp[1:length(ntemp)-1]+1] = 1
		index[1] = 1
		index = cumsum(index)
		bigdrift = ndrift[index]
		
		#presample vector of noise to sample from
		maxt<-max(ntemp) #maximum possible RT (in ms)
		noises<-zrnorm(maxt*nsim)*sigma #pre-assemble vector of noise to sample from 
		
		#presample across-trial noise
		dnoises<-zrnorm(nsim)*dsigma
		
		j<-0
		
		#loop over number of simulations per trial
		for (i in 1:nsim){
			trt<-Inf
			j<-j+1
			noise<-noises[((j-1)*maxt+1):((j)*maxt)] #sample noise realizations
			evidence<-bigdrift+noise+dnoises[i] #add drift to noise
			RDV<-cumsum(evidence)#sum up the evidence
			absRDV<-abs(RDV)
			trt<-min(which(absRDV>=1))#find the first time that the evidence exceeds a magnitude of 1
			ifelse(trt==Inf,trt<-maxt,trt<-trt)#set maxt as the finishing time if time runs out
			
			#save the trial data
			fix.num<-min(which(cumsum(nfixlength)>=trt))
			nfix[counter:(counter+fix.num-1)]<-fix.num
			choice[counter:(counter+fix.num-1)]<-ceiling(RDV[trt]/1000000) #was the RDV + or -?
			rt[counter:(counter+fix.num-1)]<-trt+ndt #add non-decision time to the RT
			fixnum[counter:(counter+fix.num-1)]<-(1:fix.num) 
			revfixnum[counter:(counter+fix.num-1)]<-(fix.num:1)
			roi[counter:(counter+fix.num-1)]<-nroi[1:fix.num]
			leftval[counter:(counter+fix.num-1)]<-lvalue_mat[k,1]
			rightval[counter:(counter+fix.num-1)]<-rvalue_mat[k,1]
			#fix the last fixation time
			tfix.length = nfixlength[1:fix.num]
			erro = sum(tfix.length)-trt
			tfix.length[fix.num] = tfix.length[fix.num]-erro
			fixdur[counter:(counter+fix.num-1)]<-tfix.length
			
			counter<-counter+fix.num
		}
	}
	sims<-data.frame(choice,rt,fixnum,nfix,revfixnum,roi,leftval,rightval,fixdur)	
	sims<-sims[sims$choice>=0,]	
	sims<-sims[sims$rt<(maxt+ndt),] #remove trials that did not finish
									rm(noises,absRDV,bigdrift,choice,counter,data_mat,drift_mat,erro,evidence,fix.num,fixdur,fixlength_mat,fixnum,i,index,j,k,leftval,lvalue_mat,maxfix,maxt,ndrift,nfix,nfixlength,noise,nroi,nsim,ntemp,ntrial,RDV,revfixnum,rightval,roi,roi_mat,rt,rvalue_mat,tfix.length,trt,value_mat)

	return(sims)
}

