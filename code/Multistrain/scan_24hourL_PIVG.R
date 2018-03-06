# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#scan_24hourL_PIVG.R


## Using the Phytoplankton-Infected-Viral-Grazer dynamic model - scan over all virus:grazer impact at several niche competition levels (as dictated in triangle of "doom" shown in schematics.R)

## Setting Life History Traits

#Phytoplankton
r1 = 1/24 # phytoplankton growth rate [per hour]
K = 2.2*10^7 # carrying capacity of phytoplankton [ cells per ml ]

#Viral
phi1 = 1*10^-8 # adsoprtion rate  [ ml/(virion . hr) ]
Beta1 = 50 # burst size, new [virions per cell]
omega1 = 0.02 # rate of inactivation of virus infectivity [per hour]    Suttle&Chen 1992 give rate in dark of 0.009-0.028 per hour; in light of 0.4 - 0.8 per hour
L1 = 24 # latent period (time between infection and lysis) [hours]
eta1 = 1/L1 #lysis rate for infected cells [per hour]

#Grazing
a1 = 2*10^-6   # grazing/clearance rate on phytoplankton   [ ml/(grazer . hr) ]
G0 = 1000  # grazer density [grazers per ml]


## Set Bottom-Up (Niche competition) levels to track

Xs = c(0.05,0.5,0.95) # low, medium, high niche comp.
LX = length(Xs)

##Set Top-Down (Lysis:Grazing) levels to track

Ys = seq(0.01,0.99,length.out=20) # lysis:grazing
LY = length(Ys)

## Dilution Simulation Parameters

#Choose dilution levels
Dil_levels = seq(0.1,1,0.1)

#Choose incubation times to measure at
TIME = 24  #incubation time of 24 h
TIMEs = 2  #incubation time of 2 h

#Storage objects
ACTUAL_LR_all = matrix(0,LX,LY)
modified_lr_all = matrix(0,LX,LY)
virdil_lr_all = matrix(0,LX,LY)
ACTUAL_LR_alls = matrix(0,LX,LY)
modified_lr_alls = matrix(0,LX,LY)
virdil_lr_alls = matrix(0,LX,LY)



for(aa in 1:LX){ # for each niche competition level
	for(bb in 1:LY){ # for each lysis:grazing level
		
		x = Xs[aa]
		y = Ys[bb]

		# condition on attack rate (a1) and carrying capacity(K)
		a1 = (eta1*(1-x)+r1*(1-x)) /(G0 *( x+ eta1*(1+y*r1/eta1)/(r1*(1-y))   +  (1-x)*(1+y*r1/eta1)/(1-y))) 
		Bhat = Beta1*eta1/(eta1+a1*G0) # to make consistent with PVG dynamics (Weitz, 2015)
		Pstar = omega1/(phi1*((Bhat*eta1)/(eta1 + a1*G0) - 1))
		K = (r1*Pstar*(1-x) + eta1*Pstar)/(x*(eta1+a1*G0))


		Nstar = omega1/( phi1*((Bhat*eta1)/(eta1+a1*G0)-1) )
		Vstar = (r1*(1-Nstar/K)-a1*G0)/( phi1*(1+r1*Nstar/(K*(eta1+a1*G0))) )
		Istar= phi1*Nstar*Vstar/(eta1+a1*G0)

		model_params = c(growth1= r1, carryingCap = K, AdsorptionRate1 = phi1 ,  burst1 = Bhat, inactivation1 = omega1, GrazingRate1=a1, LatentPeriod=L1)	
		PROJECTINFO = list(filename="model.RData",modelused="gnv")

		model_state = list(c( N1 = Nstar, I1 = Istar, V1= Vstar, G = G0))[[1]]
		GG1_b_cdil = ANALYSIS(GNIV,model_params,model_state,c(1,2),c(1,2,4),Dil_levels,0,PROJECTINFO) #classic dilution -- V not diluted, P is diluted, I is diluted, G is diluted
		GG1_b_mdil = ANALYSIS(GNIV,model_params,model_state,c(1,2),c(1,2,3,4),Dil_levels,0,PROJECTINFO) # modified dilution method -- everything diluted
		GG1_b_vdil = ANALYSIS(GNIV,model_params,model_state,c(1,2),c(3),Dil_levels,0,PROJECTINFO)  #virus dilution -- only viruses are diluted
	
		#@24h
		INDEX = which( (abs(GG1_b_cdil$Timings-TIME)) == min((abs(GG1_b_cdil$Timings-TIME))) ) # find time closest to incubation time
	
		ACTUAL_LR_all[aa,bb] = eta1*Istar/(Nstar+Istar) #lysis rate

		modified_lr_all[aa,bb] = -(GG1_b_mdil$All_DY[INDEX] - GG1_b_cdil$All_DY[INDEX])
	
		virdil_lr_all[aa,bb] = -GG1_b_vdil$All_DY[INDEX]

		#@2h
		INDEXs = which( (abs(GG1_b_cdil$Timings-TIMEs)) == min((abs(GG1_b_cdil$Timings-TIMEs))) )

		ACTUAL_LR_alls[aa,bb] = eta1*Istar/(Nstar+Istar) #lysis rate

		modified_lr_alls[aa,bb] = -(GG1_b_mdil$All_DY[INDEXs] - GG1_b_cdil$All_DY[INDEXs])

		virdil_lr_alls[aa,bb] = -GG1_b_vdil$All_DY[INDEXs]

	}
}


## Plotting


#@24 h
THICK=5
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), oma=c(0,0,0,0))
IND=1
plot(Ys*100,24*ACTUAL_LR_all[IND,],xlab="%Top-down mortality by viruses",ylim=c(0,1),ylab="Estimated Rate (per day)",lty=2,lwd=5,font.axis=2,font.lab=2)
legend("topleft",bty="n",c("Model lysis rate input","Modified dilution estimate","Viral dilution estimate"),pch=c(1,NA,NA),lty=c(NA,1,1),lwd=THICK, col=colors()[c(24,12,616)])
lines(Ys*100,24*modified_lr_all[IND,],col=colors()[c(12)],lwd=THICK)
lines(Ys*100,24*virdil_lr_all[IND,],col=colors()[c(616)],lwd=THICK)
plot(Ys*100,virdil_lr_all[IND,]/ACTUAL_LR_all[IND,],col=colors()[c(616)],lwd=THICK, ylab="Bias",xlab="%Top-down mortality by viruses",ylim=c(0,1),type="l",font.axis=2,font.lab=2)
lines(Ys*100,modified_lr_all[IND,]/ACTUAL_LR_all[IND,],col=colors()[c(12)],lwd=THICK)
abline(h=1, lwd=5, lty=2)

IND=2
plot(Ys*100,24*ACTUAL_LR_all[IND,],xlab="%Top-down mortality by viruses",ylim=c(0,1),ylab="Estimated Rate (per day)",lty=2,lwd=5,font.axis=2,font.lab=2)
lines(Ys*100,24*modified_lr_all[IND,],col=colors()[c(12)],lwd=THICK)
lines(Ys*100,24*virdil_lr_all[IND,],col=colors()[c(616)],lwd=THICK)
plot(Ys*100,virdil_lr_all[IND,]/ACTUAL_LR_all[IND,],col=colors()[c(616)],lwd=THICK, ylab="Bias",xlab="%Top-down mortality by viruses",ylim=c(0,1),type="l",font.axis=2,font.lab=2)
lines(Ys*100,modified_lr_all[IND,]/ACTUAL_LR_all[IND,],col=colors()[c(12)],lwd=THICK)
abline(h=1, lwd=5, lty=2)


IND=3
plot(Ys*100,24*ACTUAL_LR_all[IND,],xlab="%Top-down mortality by viruses",ylim=c(0,1),ylab="Estimated Rate (per day)",lty=2,lwd=5,font.axis=2,font.lab=2)
lines(Ys*100,24*modified_lr_all[IND,],col=colors()[c(12)],lwd=THICK)
lines(Ys*100,24*virdil_lr_all[IND,],col=colors()[c(616)],lwd=THICK)
plot(Ys*100,virdil_lr_all[IND,]/ACTUAL_LR_all[IND,],col=colors()[c(616)],lwd=THICK, ylab="Bias",xlab="%Top-down mortality by viruses",ylim=c(0,1),type="l",font.axis=2,font.lab=2)
lines(Ys*100,modified_lr_all[IND,]/ACTUAL_LR_all[IND,],col=colors()[c(12)],lwd=THICK)
abline(h=1, lwd=5, lty=2)

dev.copy2eps(file=paste("../../figures/scanPIVG_24hour_24h",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/scanPIVG_24hour_24h",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/scanPIVG_24hour_24h",".jpg",sep=""))
dev.off()
dev.off()


#@ 2nd time

THICK=5
par(mfrow=c(3,2), mar=c(4,4,0.5,0.5), oma=c(0,0,0,0))
IND=1
plot(Ys*100,24*ACTUAL_LR_alls[IND,],xlab="%Top-down mortality by viruses",ylim=c(0,1),ylab="Estimated Rate (per day)",lty=2,lwd=5,font.axis=2,font.lab=2)
legend("topleft",bty="n",c("Model lysis rate input","Modified dilution estimate","Viral dilution estimate"),pch=c(1,NA,NA),lty=c(NA,1,1),lwd=THICK, col=colors()[c(24,12,616)])
lines(Ys*100,24*modified_lr_alls[IND,],col=colors()[c(12)],lwd=THICK)
lines(Ys*100,24*virdil_lr_alls[IND,],col=colors()[c(616)],lwd=THICK)
plot(Ys*100,virdil_lr_alls[IND,]/ACTUAL_LR_alls[IND,],col=colors()[c(616)],lwd=THICK, ylab="Bias",xlab="%Top-down mortality by viruses",ylim=c(0,1),type="l",font.axis=2,font.lab=2)
lines(Ys*100,modified_lr_alls[IND,]/ACTUAL_LR_alls[IND,],col=colors()[c(12)],lwd=THICK)
abline(h=1, lwd=5, lty=2)

IND=2
plot(Ys*100,24*ACTUAL_LR_alls[IND,],xlab="%Top-down mortality by viruses",ylim=c(0,1),ylab="Estimated Rate (per day)",lty=2,lwd=5,font.axis=2,font.lab=2)
lines(Ys*100,24*modified_lr_alls[IND,],col=colors()[c(12)],lwd=THICK)
lines(Ys*100,24*virdil_lr_alls[IND,],col=colors()[c(616)],lwd=THICK)
plot(Ys*100,virdil_lr_alls[IND,]/ACTUAL_LR_alls[IND,],col=colors()[c(616)],lwd=THICK, ylab="Bias",xlab="%Top-down mortality by viruses",ylim=c(0,1),type="l",font.axis=2,font.lab=2)
lines(Ys*100,modified_lr_alls[IND,]/ACTUAL_LR_alls[IND,],col=colors()[c(12)],lwd=THICK)
abline(h=1, lwd=5, lty=2)


IND=3
plot(Ys*100,24*ACTUAL_LR_alls[IND,],xlab="%Top-down mortality by viruses",ylim=c(0,1),ylab="Estimated Rate (per day)",lty=2,lwd=5,font.axis=2,font.lab=2)
lines(Ys*100,24*modified_lr_alls[IND,],col=colors()[c(12)],lwd=THICK)
lines(Ys*100,24*virdil_lr_alls[IND,],col=colors()[c(616)],lwd=THICK)
plot(Ys*100,virdil_lr_alls[IND,]/ACTUAL_LR_alls[IND,],col=colors()[c(616)],lwd=THICK, ylab="Bias",xlab="%Top-down mortality by viruses",ylim=c(0,1),type="l",font.axis=2,font.lab=2)
lines(Ys*100,modified_lr_alls[IND,]/ACTUAL_LR_alls[IND,],col=colors()[c(12)],lwd=THICK)
abline(h=1, lwd=5, lty=2)


dev.copy2eps(file=paste("../../figures/scanPIVG_24hour_2h",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/scanPIVG_24hour_2h",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/scanPIVG_24hour_2h",".jpg",sep=""))
dev.off()
dev.off()




