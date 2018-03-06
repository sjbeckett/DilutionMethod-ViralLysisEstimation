# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.

# MIT License


ANALYSIS <- function(MODEL, parameters, model_initial_state, host_vector, dilute_vector, dilution_levels, saveon, PROJECTINFO) {
	#MODEL - dynamical system model to use with library deSolve
	#parameters - parameterisation of above dynamical system
	#model_initial_state - initial conditions in dyn. sys.  e.g. c(100,6000,4000)
	#host_vector -- the state variables that are counted as belonging to the plankton -i.e. what is counted as changing in abundance in dilution experiment e.g. c(1,2)
	#dilute_vector -- the state variables to dilute i.e. those that are assumed to not be able to pass through the filter into the filtrate. e.g. c(1,2,3)
	#dilution_levels -- what levels of dilution are each bottle at?
	#PROJECTINFO - name to save data to

	#length shortcuts
	LState = length(model_initial_state)
	LHvec = length(host_vector)
	LDvec = length(dilute_vector)

	# ERROR CHECKING BEGIN
	if(LHvec<1 || LHvec>LState) {
		stop("host_vector is too short or too long")
	}
	if(min(host_vector)<1 || max(host_vector)> LState){
		stop("host_vector contains bad indexes")
	}
	if(LDvec<1 || LDvec>LState) {
		stop("dilution_vector is too short or too long")
	}
	if(min(dilute_vector)<1 || max(dilute_vector)> LState){
		stop("dilution_vector contains bad indexes")
	}
	if(min(dilution_levels)<0){
		stop("dilution level is negative amount of WSW!")
	}
	if(max(dilution_levels)>1){
		warning("This is a concentration experiment. WSW is > 1")
	}
	if(saveon!=0 && saveon!=1){
		stop("saveon is not assigned as 1 (save) or 0 (no save)")
	}
	# ERROR CHECKING END

	library(deSolve)
	print(PROJECTINFO)

	#Simulation parameters
	TimeLength = 1.5*24 # (in hours)
	RecordingTime = (1/6)  # every 10 minutes (in hours)
	SimTimes = seq(0,TimeLength,RecordingTime)  # all the time points that are to be recorded (in hours)

	#Dilution levels
	FRACTIONS = dilution_levels # PROPORTIONS OF ORIGINAL SAMPLE TO USE  i.e. the proportions of WSW in each bottle.

	INITIALSTATE = t(t(FRACTIONS))%*%t((model_initial_state))  # create matrix of steady states. Each row is a different FRACTION, columns are the the different steady states
	#This is only applicable if diluting everything.....


	#If don't want to dilute some things then need to correct for this:
	DONTDILUTE = setdiff(1:LState,dilute_vector) # find the state variables do not want to dilute
	if(length(DONTDILUTE)>0) {  #if something needs editing so no dilution performed...
		for(www in 1:length(DONTDILUTE)) {
			INITIALSTATE[,DONTDILUTE[www]] = model_initial_state[DONTDILUTE[www]]
		}
	}
	
	#NOW INITIALSTATE contains initial conditions for each bottle; conditioned on whether they can pass through filter or not.


	#VectorLengthShortcuts
	LTim = length(SimTimes)
	LF = length(FRACTIONS)
	

	#COLLECTION
	RECORDALL = array( 0, dim=c(LTim,LF,LState) )
	
	#Simulation


	for(bb in 1:LF) {
		#setTxtProgressBar(pb,bb)
		state_use = INITIALSTATE[bb,]
		#print(state_use)
		EXP = ode(y = state_use, times = SimTimes, func = MODEL, parms = parameters, method='ode45') # run the experiment
		
		#place timeseries to an array
		for(cc in 1:LState) {
			RECORDALL[,bb,cc] = EXP[,(cc+1)]
		}
			
	}


	state_use = INITIALSTATE[LF,]



	##  Find apparent growth rates for each host (and find the average for all hosts) for each dilution level
	print("finding app. growth rates")
	pb = txtProgressBar(min = 1, max = LTim, style = 3)
	Appgrowth = array( 0, dim=c(LTim,LF,LHvec) )  #want app growth of all hosts
	AverageAppgrowth = array(0,dim=c(LTim,LF))


	HOST_TOT_t0 = rowSums(as.matrix(INITIALSTATE[,host_vector]))  # sum of host populations at t=0 for each dilution level
	hostAbnow = c() #initialise empty vector
	
	for(aa in 2:LTim) { #can't measure at t=0 for each time
		setTxtProgressBar(pb,aa)
		inverseTime = (1/SimTimes[aa])
		for(bb in 1:LF) { #for each fraction
			for(cc in 1:LHvec) { # for each host
				dd=host_vector[cc] # host index in record all
				hostAbnow[cc] = RECORDALL[aa,bb,dd]
				Appgrowth[aa,bb,cc] = inverseTime*log(hostAbnow[cc]/RECORDALL[1,bb,dd])
			}

				AverageAppgrowth[aa,bb] = inverseTime*log(sum(hostAbnow)/HOST_TOT_t0[bb])
		}
	}

	close(pb)
	


	##  Linear regressions to determine growth and decay
	print("performing full linear regression")
	pb = txtProgressBar(min = 1, max = LTim*2, style = 3)

	Mu0_INDIV = array(0,dim=c(LTim,LHvec))
	DECAY_INDIV = array(0,dim=c(LTim,LHvec))
	AdjR_INDIV = array(0,dim=c(LTim,LHvec))


	Mu0_ALL = array(0,dim=c(LTim))
	DECAY_ALL = array(0,dim=c(LTim))
	AdjR_ALL = array(0,dim=c(LTim))


	for(aa in 2:LTim) {
		setTxtProgressBar(pb,aa)
		for(cc in 1:LHvec) {
			#for each host individually				
			REGMODEL_INDIV = lm(Appgrowth[aa,,cc] ~ FRACTIONS)
			Mu0_INDIV[aa,cc] = REGMODEL_INDIV[[1]][[1]]
			DECAY_INDIV[aa,cc] = REGMODEL_INDIV[[1]][[2]]
			STATS = summary(REGMODEL_INDIV)
			AdjR_INDIV[aa,cc]  = STATS$adj.r.squared
		}
			
			#ALL
			REGMODEL_ALL = lm(AverageAppgrowth[aa,] ~ FRACTIONS)
			Mu0_ALL[aa] = REGMODEL_ALL[[1]][[1]]
			DECAY_ALL[aa] = REGMODEL_ALL[[1]][[2]]
			STATS = summary(REGMODEL_ALL)
			AdjR_ALL[aa]  = STATS$adj.r.squared
	}

	

	#repeat for "two point" dilution method
	TP = c(1,LF)  # first and last dilution

	Mu0_INDIV_TP = array(0,dim=c(LTim,LHvec))
	DECAY_INDIV_TP = array(0,dim=c(LTim,LHvec))
	AdjR_INDIV_TP = array(0,dim=c(LTim,LHvec))


	Mu0_ALL_TP = array(0,dim=c(LTim))
	DECAY_ALL_TP = array(0,dim=c(LTim))
	AdjR_ALL_TP = array(0,dim=c(LTim))


	for(aa in 2:LTim) {
		setTxtProgressBar(pb,LTim+aa)
		for(cc in 1:LHvec) {
			#for each host individually				
			REGMODEL_INDIV = lm(Appgrowth[aa,TP,cc] ~ FRACTIONS[TP])
			Mu0_INDIV_TP[aa,cc] = Appgrowth[aa,1,cc] # use the first dilution as intrinsic growth rate
			DECAY_INDIV_TP[aa,cc] = REGMODEL_INDIV[[1]][[2]]
			STATS = summary(REGMODEL_INDIV)
			AdjR_INDIV_TP[aa,cc]  = STATS$adj.r.squared
		}
			
			#ALL
			REGMODEL_ALL = lm(AverageAppgrowth[aa,TP] ~ FRACTIONS[TP])
			Mu0_ALL_TP[aa] = AverageAppgrowth[aa,1] # use first dilution as intrinsic growth rate
			DECAY_ALL_TP[aa] = REGMODEL_ALL[[1]][[2]]
			STATS = summary(REGMODEL_ALL)
			AdjR_ALL_TP[aa]  = STATS$adj.r.squared
	}

	close(pb)





	TOSAVE = list(info = PROJECTINFO, alltimeseries = RECORDALL, Timings = SimTimes, Ind_appG = Appgrowth, All_appG = AverageAppgrowth, Ind_mu0 = Mu0_INDIV, Ind_mu0_tp = Mu0_INDIV_TP, Ind_DY = DECAY_INDIV, Ind_DY_tp = DECAY_INDIV_TP, All_mu0 = Mu0_ALL, All_mu_tp = Mu0_ALL_TP, All_DY=DECAY_ALL, All_DY_tp = DECAY_ALL_TP, modelused= MODEL, parametersused = as.list(parameters), InitialState = model_initial_state, hostVector = host_vector, adjRsqALL = AdjR_ALL, adjRsqTP = AdjR_ALL_TP)
	
	if(saveon == 1) {
	save(TOSAVE, file=PROJECTINFO$filename)
	}

	return(TOSAVE)
}





