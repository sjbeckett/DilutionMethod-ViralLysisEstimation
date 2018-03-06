# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

# LHS_PVG.R

# LHS scan
library(lhs)

#Parameter ranges



#LHS design

samples = 500

R_r = c(0.1/24,2/24)
R_K = 10^c(6,8)

R_B = c(10,100)
R_phi = 10^c(-10,-7)
R_omega = c(0.2,5)/24

R_a = 10^c(-7,-4)

G0 = 1000

SAMPLING = randomLHS(n=samples,k=6)
#gives matrix of n rows and k columns

#Storage

#Store results
	ACTUAL_LR_all = matrix(0,samples,1)
	modified_lr_all = matrix(0,samples,1)
	virdil_lr_all = matrix(0,samples,1)
	ACTUAL_LR_alls = matrix(0,samples,1)
	modified_lr_alls = matrix(0,samples,1)
	virdil_lr_alls = matrix(0,samples,1)


#Dilution levels
Dil_levels = seq(0.1,1,0.1)

#Choose incubation times to measure at
TIME = 24  #incubation time of 24 h
TIMEs = 2  #incubation time of 2 h
bb=1


for(aa in 1:samples){ #for each sample
	#Find parameters rom LHS
	r1 = min(R_r)+ diff(R_r)*SAMPLING[aa,1]
	K = min(R_K)+ diff(R_K)*SAMPLING[aa,2]

	Beta1 = min(R_B)+ diff(R_B)*SAMPLING[aa,3]
	phi1 = min(R_phi)+ diff(R_phi)*SAMPLING[aa,4]
	omega1 = min(R_omega)+ diff(R_omega)*SAMPLING[aa,5]

	a1 = min(R_a)+ diff(R_a)*SAMPLING[aa,6]

	model_params = c(growth1= r1, carryingCap = K, AdsorptionRate1 = phi1 ,  burst1 = Beta1, inactivation1 = omega1, GrazingRate1=a1)	

	#Set the steady states for this sim.
	Nstar = omega1/((Beta1-1)*phi1)
	Vstar = (1/phi1) * (r1*( 1- (omega1/((Beta1-1)*phi1*K)) ) - a1*G0)
	if(Nstar>0 & Vstar>0){

		#simulations
		PROJECTINFO = list(filename="model.RData",modelused="gnv")

		model_state = list(c( N1 = Nstar, V1= Vstar, G = G0))[[1]]
		GG1_b_cdil = ANALYSIS(GNV,model_params,model_state,c(1),c(1,3),Dil_levels,0,PROJECTINFO) #classic dilution --  V not diluted, P is diluted, G is diluted
		GG1_b_mdil = ANALYSIS(GNV,model_params,model_state,c(1),c(1,2,3),Dil_levels,0,PROJECTINFO) # modified dilution method -- everything diluted
		GG1_b_vdil = ANALYSIS(GNV,model_params,model_state,c(1),c(2),Dil_levels,0,PROJECTINFO)  #virus dilution -- only viruses are diluted

		INDEX = which( (abs(GG1_b_cdil$Timings-TIME)) == min((abs(GG1_b_cdil$Timings-TIME))) ) # find time closest to incubation time
	
		ACTUAL_LR_all[aa,bb] = Vstar*phi1 #lysis rate

		modified_lr_all[aa,bb] = -(GG1_b_mdil$All_DY[INDEX] - GG1_b_cdil$All_DY[INDEX])

		virdil_lr_all[aa,bb] = -GG1_b_vdil$All_DY[INDEX]

		INDEXs = which( (abs(GG1_b_cdil$Timings-TIMEs)) == min((abs(GG1_b_cdil$Timings-TIMEs))) )

		ACTUAL_LR_alls[aa,bb] = Vstar*phi1 #lysis rate

		modified_lr_alls[aa,bb] = -(GG1_b_mdil$All_DY[INDEXs] - GG1_b_cdil$All_DY[INDEXs])

		virdil_lr_alls[aa,bb] = -GG1_b_vdil$All_DY[INDEXs]
	}

}

#calculate bias

BIAS_M24 = modified_lr_all/ACTUAL_LR_all
BIAS_V24 = virdil_lr_all/ACTUAL_LR_all

BIAS_M2 = modified_lr_alls/ACTUAL_LR_alls
BIAS_V2 = virdil_lr_alls/ACTUAL_LR_alls

#24 hours

dens1x = density(BIAS_M24,na.rm=TRUE)
dens2x = density(BIAS_V24,na.rm=TRUE)

plot(dens1x,ylim = range(dens1x$y,dens2x$y),xlim=range(dens1x$x,dens2x$x),col=colors()[c(12)],lwd=5,main="",xlab=paste("Bias in estimation of viral lysis rate, N=",dens1x$n))
lines(dens2x,col=colors()[c(616)],lwd=5)
legend("topleft",c("Modified dilution estimate","Viral dilution estimate"),lty=c(1,1),lwd=c(5,5),col=colors()[c(12,616)],bty="n")
print(dens1x$n)


dev.copy2eps(file=paste("../../figures/densityPVG_24h",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/densityPVG_24h",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/densityPVG_24h",".jpg",sep=""))
dev.off()
dev.off()

#2 hours

dens1 = density(BIAS_M2,na.rm=TRUE)
dens2 = density(BIAS_V2,na.rm=TRUE)

plot(dens1,ylim = range(dens1$y,dens2$y),xlim=range(dens1$x,dens2$x),col=colors()[c(12)],lwd=5,main="",xlab=paste("Bias in estimation of viral lysis rate, N=",dens1$n))
lines(dens2,col=colors()[c(616)],lwd=5)
legend("topleft",c("Modified dilution estimate","Viral dilution estimate"),lty=c(1,1),lwd=c(5,5),col=colors()[c(12,616)],bty="n")
print(dens1$n)


dev.copy2eps(file=paste("../../figures/densityPVG_2h",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/densityPVG_2h",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/densityPVG_2h",".jpg",sep=""))
dev.off()
dev.off()

#Both on one plot
dev.new(width=8.625000,height= 4.716754)
par(mfrow=c(1,2), mar=c(4,4,1,0.5), oma=c(0,0,0,0))
plot(dens1,ylim = range(dens1$y,dens2$y),xlim=range(dens1$x,dens2$x),col=colors()[c(12)],lwd=5,main="2h incubation",xlab=paste("Bias in estimation of viral lysis rate, N=",dens1$n))
lines(dens2,col=colors()[c(616)],lwd=5)
legend("topleft",c("Modified dilution estimate","Viral dilution estimate"),lty=c(1,1),lwd=c(5,5),col=colors()[c(12,616)],bty="n")
print(dens1$n)

plot(dens1x,ylim = range(dens1x$y,dens2x$y),xlim=range(dens1x$x,dens2x$x),col=colors()[c(12)],lwd=5,main="24h incubation",xlab=paste("Bias in estimation of viral lysis rate, N=",dens1x$n))
lines(dens2x,col=colors()[c(616)],lwd=5)
legend("topleft",c("Modified dilution estimate","Viral dilution estimate"),lty=c(1,1),lwd=c(5,5),col=colors()[c(12,616)],bty="n")
print(dens1$n)

dev.copy2eps(file=paste("../../figures/densityPVG_both",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/densityPVG_both",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/densityPVG_both",".jpg",sep=""))
dev.off()
dev.off()





