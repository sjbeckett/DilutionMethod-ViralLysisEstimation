# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#diversity_concept.R


# Parameterisation

#Experimental
Dil_levels = seq(0.1,1,0.1)
TIME1 = 2
TIME2 = 24

#Biological
r1 = 1/24 #plankton growth rate [ per hour ]
r2 = 2/24 #plankton growth rate [ per hour ]

K = 2.2*10^7 # plankton per ml

phi1 = 1*10^-8 # adsoprtion rate  [ ml/(virion . hr) ]
phi2 = 1*10^-7 # adsoprtion rate  [ ml/(virion . hr) ]

Beta1 = 50 # Burst size, [ new virions per cell ]
Beta2 = 50 # Burst size, [ new virions per cell ]

omega1 = 0.02 # inactivation of virus infectivity [ per hour ] 
omega2 = 0.02 # inactivation of virus infectivity [ per hour ]

a1 = 2*10^-6   # grazing/clearance rate on phytoplankton 1   [ ml/(grazer . hr) ]
a2 = a1   # grazing/clearance rate on phytoplankton 2  [ ml(grazer . hr) ]
G0 = 1000  # grazer density [ grazers per ml ]


#Want to explore strain space by changing relative growth rates and abundances (set by virus traits - see steady state conditions) between the two strains

SQUARESIZE = 50 # number of levels to try in each direction

#Fix P_2 trait and abundance
r2 = 0.2/24  # per day
PHI2 = 1*10^-10  # mililitres per virus. hour
phi2=PHI2

#Choose P_1 traits and abundance
r1s = seq(0.5,50,length.out=SQUARESIZE)*r2 #up to 10 per day
PHI1s = exp(seq(log(1), log(1000), length.out = SQUARESIZE))*PHI2 

#choose particular r1 and phi1
R_IND = 29 # in range of 1 - SQUARESIZE
P_IND = 21 # in range of 1 - SQUARESIZE

r1 = r1s[R_IND]
phi1 = PHI1s[P_IND]

#find steady states
		N1star = omega1/((Beta1-1)*phi1)
		N2star = omega2/((Beta2-1)*phi2)

		V1star = (1/phi1)*(r1*(1 - (N1star+N2star)/K) - a1*G0)
		V2star = (1/phi2)*(r2*(1 - (N1star+N2star)/K) - a2*G0)

#Run simulation

model_params = c(growth1= r1, growth2 = r2, carryingCap = K, AdsorptionRate1 = phi1 , AdsorptionRate2 = phi2, burst1 = Beta1, burst2 = Beta2, inactivation1 = omega1, inactivation2 = omega2 , GrazingRate1 = a1, GrazingRate2 = a2 )
		
PROJECTINFO = list(filename="model2.RData",modelused="n1v1_n2v2_g")

model_state = list(c( N1 = N1star, N2 = N2star, V1= V1star, V2 = V2star, G = G0))[[1]]

	GG1_b_cdil = ANALYSIS(NV2G,model_params,model_state,c(1,2),c(1,2,5),Dil_levels,0,PROJECTINFO) #classic dilution -- dilute P and G
	GG1_b_mdil = ANALYSIS(NV2G,model_params,model_state,c(1,2),c(1,2,3,4,5),Dil_levels,0,PROJECTINFO) #modified dilution  -- dilute V,P and G
	GG1_b_vdil = ANALYSIS(NV2G,model_params,model_state,c(1,2),c(3,4),Dil_levels,0,PROJECTINFO) #virus dilution -- dilute V

#Make conceptual plot(s)

OBJ = GG1_b_mdil
PRE_EXP_TIME = 12/(1/6) # timesteps before dilution occurs (24 hours is 24/(1/6) as 10 minute steps)
DIL_USED = 3 # which dilution level to show
H1 = rep(N1star,PRE_EXP_TIME)
H2 = rep(N2star,PRE_EXP_TIME)

H1 = c(H1, OBJ$alltimeseries[,DIL_USED,1]) # time, dilution, pop
H2 = c(H2, OBJ$alltimeseries[,DIL_USED,2])

BOTH = H1+H2
TIMES = c(rev(-1:-PRE_EXP_TIME*(OBJ$Timings[2]-OBJ$Timings[1])),OBJ$Timings)



#plot dilution plot, showing steady states and recovery

par(mfrow=c(1,2), mar=c(4,4,0.5,0.5), oma=c(0,0,0,0))

plot(0,0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",col="white")

par(mar=c(4,4.8,2,1))
par(oma=c(0,0,0,0))
LWD=5

YLIM=c(0,1.2*max(BOTH))
plot(TIMES,BOTH,type='l',bty='n',ylim=YLIM,xaxt='n',xlab="Time (h)", ylab=expression("Cells ml" ^-1),cex.axis=1.2,cex.lab=1.3,lwd=LWD)
lines(TIMES,H1,col='green',lwd=LWD)
lines(TIMES,H2,col='red',lwd=LWD)
axis(side=1,at=seq(0,120,24),cex.axis=1.2)


TXT=1.6

#text to show what is what
timeIND = 37  #12 hours before dil
time=TIMES[timeIND]
text(time-2,H1[timeIND]+200000,'P1',col='green',cex=TXT)
text(time,H2[timeIND]-200000,'P2',col='red',cex=TXT)
text(time-2,BOTH[timeIND]+200000,'P1+P2',cex=TXT)


## plot lysis rates

TIME= 24
INDEX = which( (abs(GG1_b_cdil$Timings-TIME)) == min((abs(GG1_b_cdil$Timings-TIME))) )
HOURTODAY = 24 #change to 1 if want to report as hourly rates rather than daily rates.

P1vir_act = HOURTODAY*(phi1*V1star)
P1vir_est_m = HOURTODAY*((-GG1_b_mdil$Ind_DY[INDEX,1]) - (-GG1_b_cdil$Ind_DY[INDEX,1]))
P1vir_est_v = HOURTODAY*(-GG1_b_vdil$Ind_DY[INDEX,1])
P2vir_act = HOURTODAY*(phi2*V2star)
P2vir_est_m = HOURTODAY*((-GG1_b_mdil$Ind_DY[INDEX,2]) - (-GG1_b_cdil$Ind_DY[INDEX,2]))
P2vir_est_v = HOURTODAY*(-GG1_b_vdil$Ind_DY[INDEX,2])

Comm_Vir_act = HOURTODAY*((phi1*V1star*N1star + phi2*V2star*N2star)/(N1star+N2star))
Comm_Vir_est_m = HOURTODAY*((-GG1_b_mdil$All_DY[INDEX]) - (-GG1_b_cdil$All_DY[INDEX]))
Comm_Vir_est_v = HOURTODAY*(-GG1_b_vdil$All_DY[INDEX])


data = cbind(c(Comm_Vir_act,Comm_Vir_est_m,Comm_Vir_est_v),c(P1vir_act,P1vir_est_m,P1vir_est_v),c(P2vir_act,P2vir_est_m,P2vir_est_v))
colnames(data)=c("Community","P1","P2")
rownames(data)=c("Model input","Modified dilution","Viral dilution")

barplot(data, col=colors()[c(23,12,616)] , border="white", font.axis=2, beside=T, legend=rownames(data), font.lab=2, args.legend=list(bty="n",ncol=1,x="topright"),ylab=paste("Viral lysis rates (per day)"),ylim=24*c(-0.005,5/24))	











dev.copy2eps(file=paste("../../figures/diversity_concept",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/diversity_concept",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/diversity_concept",".jpg",sep=""))
dev.off()
dev.off()












