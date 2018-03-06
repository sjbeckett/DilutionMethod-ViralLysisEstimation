# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#multistrain_PPVVG.R


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
r1s = seq(0.5,50,length.out=SQUARESIZE)*r2 #up to 50*r2 per day
PHI1s = exp(seq(log(1), log(1000), length.out = SQUARESIZE))*PHI2 

N1star = omega1/((Beta1-1)*PHI1s)
N2star = omega2/((Beta2-1)*phi2)

V1star = (1/PHI1s)*(r1s*(1 - (N1star+N2star)/K) - a1*G0)
V2star = (1/phi2)*(r2*(1 - (N1star+N2star)/K) - a2*G0)


#Storage objects
#Bulk community lysis rates
ACTUAL_LR_all_24h = matrix(0,SQUARESIZE,SQUARESIZE)
modified_lr_all_24h = matrix(0,SQUARESIZE,SQUARESIZE)
virdil_lr_all_24h = matrix(0,SQUARESIZE,SQUARESIZE)

ACTUAL_LR_all_2h = matrix(0,SQUARESIZE,SQUARESIZE)
modified_lr_all_2h = matrix(0,SQUARESIZE,SQUARESIZE)
virdil_lr_all_2h = matrix(0,SQUARESIZE,SQUARESIZE)

#Type 1 lysis rates
ACTUAL_LR_P1_24h = matrix(0,SQUARESIZE,SQUARESIZE)
modified_lr_P1_24h = matrix(0,SQUARESIZE,SQUARESIZE)
virdil_lr_P1_24h = matrix(0,SQUARESIZE,SQUARESIZE)

ACTUAL_LR_P1_2h = matrix(0,SQUARESIZE,SQUARESIZE)
modified_lr_P1_2h = matrix(0,SQUARESIZE,SQUARESIZE)
virdil_lr_P1_2h = matrix(0,SQUARESIZE,SQUARESIZE)

#Type 2 lysis rates
ACTUAL_LR_P2_24h = matrix(0,SQUARESIZE,SQUARESIZE)
modified_lr_P2_24h = matrix(0,SQUARESIZE,SQUARESIZE)
virdil_lr_P2_24h = matrix(0,SQUARESIZE,SQUARESIZE)

ACTUAL_LR_P2_2h = matrix(0,SQUARESIZE,SQUARESIZE)
modified_lr_P2_2h = matrix(0,SQUARESIZE,SQUARESIZE)
virdil_lr_P2_2h = matrix(0,SQUARESIZE,SQUARESIZE)

#ABUND
P1STAR = matrix(0,SQUARESIZE,SQUARESIZE)
P2STAR = matrix(0,SQUARESIZE,SQUARESIZE)
V1STAR = matrix(0,SQUARESIZE,SQUARESIZE)
V2STAR = matrix(0,SQUARESIZE,SQUARESIZE)



for(aa in 1:length(r1s)) { # for all growth rates considered
	for(bb in 1:length(PHI1s)) { # for all adsorption rates considered

		#assign P_1 and V_1 traits
		r1 = r1s[aa]
		phi1 = PHI1s[bb]

		#find steady states
		N1star = omega1/((Beta1-1)*phi1)
		N2star = omega2/((Beta2-1)*phi2)

		V1star = (1/phi1)*(r1*(1 - (N1star+N2star)/K) - a1*G0)
		V2star = (1/phi2)*(r2*(1 - (N1star+N2star)/K) - a2*G0)

		P1STAR[aa,bb] = N1star
		P2STAR[aa,bb] = N2star
		V1STAR[aa,bb] = V1star
		V2STAR[aa,bb] = V2star
		
		model_params = c(growth1= r1, growth2 = r2, carryingCap = K, AdsorptionRate1 = phi1 , AdsorptionRate2 = phi2, burst1 = Beta1, burst2 = Beta2, inactivation1 = omega1, inactivation2 = omega2 , GrazingRate1 = a1, GrazingRate2 = a2 )
		
		PROJECTINFO = list(filename="model2.RData",modelused="n1v1_n2v2_g")

		model_state = list(c( N1 = N1star, N2 = N2star, V1= V1star, V2 = V2star, G = G0))[[1]]

		GG1_b_cdil = ANALYSIS(NV2G,model_params,model_state,c(1,2),c(1,2,5),Dil_levels,0,PROJECTINFO) #classic dilution -- dilute P and G
		GG1_b_mdil = ANALYSIS(NV2G,model_params,model_state,c(1,2),c(1,2,3,4,5),Dil_levels,0,PROJECTINFO) #modified dilution  -- dilute V,P and G
		GG1_b_vdil = ANALYSIS(NV2G,model_params,model_state,c(1,2),c(3,4),Dil_levels,0,PROJECTINFO) #virus dilution -- dilute V


		#@2h
		INDEXs = which( (abs(GG1_b_cdil$Timings-TIME1)) == min((abs(GG1_b_cdil$Timings-TIME1))) )

		ACTUAL_LR_all_2h[aa,bb] = ((phi1*V1star)*N1star + (phi2*V2star)*N2star)/( N1star + N2star)
		ACTUAL_LR_P1_2h[aa,bb] = phi1*V1star
		ACTUAL_LR_P2_2h[aa,bb] = phi2*V2star

		modified_lr_all_2h[aa,bb] = -(GG1_b_mdil$All_DY[INDEXs] - GG1_b_cdil$All_DY[INDEXs])
		modified_lr_P1_2h[aa,bb] = -(GG1_b_mdil$Ind_DY[INDEXs,1] - GG1_b_cdil$Ind_DY[INDEXs,1])
		modified_lr_P2_2h[aa,bb] = -(GG1_b_mdil$Ind_DY[INDEXs,2] - GG1_b_cdil$Ind_DY[INDEXs,2])

		virdil_lr_all_2h[aa,bb] = -GG1_b_vdil$All_DY[INDEXs]
		virdil_lr_P1_2h[aa,bb] =  -GG1_b_vdil$Ind_DY[INDEXs,1]
		virdil_lr_P2_2h[aa,bb] =  -GG1_b_vdil$Ind_DY[INDEXs,2]

		#@24h
		INDEX = which( (abs(GG1_b_cdil$Timings-TIME2)) == min((abs(GG1_b_cdil$Timings-TIME2))) ) # find time closest to incubation time
	
		ACTUAL_LR_all_24h[aa,bb] = ((phi1*V1star)*N1star + (phi2*V2star)*N2star)/( N1star + N2star)
		ACTUAL_LR_P1_24h[aa,bb] = phi1*V1star
		ACTUAL_LR_P2_24h[aa,bb] = phi2*V2star

		modified_lr_all_24h[aa,bb] = -(GG1_b_mdil$All_DY[INDEX] - GG1_b_cdil$All_DY[INDEX])
		modified_lr_P1_24h[aa,bb] = -(GG1_b_mdil$Ind_DY[INDEX,1] - GG1_b_cdil$Ind_DY[INDEX,1])
		modified_lr_P2_24h[aa,bb] = -(GG1_b_mdil$Ind_DY[INDEX,2] - GG1_b_cdil$Ind_DY[INDEX,2])

		virdil_lr_all_24h[aa,bb] = -GG1_b_vdil$All_DY[INDEX]
		virdil_lr_P1_24h[aa,bb] =  -GG1_b_vdil$Ind_DY[INDEX,1]
		virdil_lr_P2_24h[aa,bb] =  -GG1_b_vdil$Ind_DY[INDEX,2]
		print(c(aa,bb))
	}
}



BIAS_M_all_2h = modified_lr_all_2h/ACTUAL_LR_all_2h
BIAS_M_all_24h = modified_lr_all_24h/ACTUAL_LR_all_24h
BIAS_M_P1_2h =  modified_lr_P1_2h/ACTUAL_LR_P1_2h
BIAS_M_P1_24h = modified_lr_P1_24h/ACTUAL_LR_P1_24h
BIAS_M_P2_2h = modified_lr_P2_2h/ACTUAL_LR_all_2h
BIAS_M_P2_24h = modified_lr_P2_24h/ACTUAL_LR_all_24h

BIAS_V_all_2h = virdil_lr_all_2h/ACTUAL_LR_all_2h
BIAS_V_all_24h = virdil_lr_all_24h/ACTUAL_LR_all_24h
BIAS_V_P1_2h = virdil_lr_P1_2h/ACTUAL_LR_P1_2h
BIAS_V_P1_24h = virdil_lr_P1_24h/ACTUAL_LR_P1_24h
BIAS_V_P2_2h = virdil_lr_P2_2h/ACTUAL_LR_P2_2h
BIAS_V_P2_24h = virdil_lr_P2_24h/ACTUAL_LR_P2_24h




save(list=ls(all.names=TRUE),file="SPACE.RData")


##PLOTTING




ZLIMcomm = range(BIAS_M_all_2h,BIAS_V_all_2h,1)
ZLIMind1 =  range(BIAS_M_P1_2h,BIAS_V_P1_2h,1)
ZLIMind2 =  range(BIAS_M_P2_2h,BIAS_V_P2_2h,1)
ZLIM24comm=range(BIAS_M_all_24h,BIAS_V_all_24h,1)
ZLIM24ind1=range(BIAS_M_P1_24h,BIAS_V_P1_24h,1)
ZLIM24ind2=range(BIAS_M_P2_24h,BIAS_V_P2_24h,1)

library(viridis)
CRP1 = viridis(50)

library(fields)

LBCX = 1 #size of contour labels
LVL = c(0,1)#c(0,0.5,0.8,1,1.2,2) # contour lines to draw
MET = "simple" #method for plotting contour labels
ConCol= c("red","white")


#2h

par(mfrow=c(3,2), mar=c(4,4.5,1,1), oma=c(0,0,0,0))

image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_M_all_2h,zlim=ZLIMcomm,col=CRP1,xlab="",ylab="",log="x",main="Modified dilution")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_M_all_2h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()

image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_V_all_2h,zlim=ZLIMcomm,col=CRP1,xlab="",ylab="",log="x",main="Virus dilution")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_V_all_2h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()

image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_M_P1_2h,zlim=ZLIMind1,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_M_P1_2h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()
image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_V_P1_2h,zlim=ZLIMind1,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_V_P1_2h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()
image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_M_P2_2h,zlim=ZLIMind2,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_M_P2_2h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()
image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_V_P2_2h,zlim=ZLIMind2,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_V_P2_2h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()



dev.copy2eps(file=paste("../../figures/PPVVG_2h",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/PPVVG_2h",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/PPVVG_2h",".jpg",sep=""))
dev.off()
dev.off()

#24

par(mfrow=c(3,2), mar=c(4,4.5,1,1), oma=c(0,0,0,0))

image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_M_all_24h,zlim=ZLIM24comm,col=CRP1,xlab="",ylab="",log="x",main="Modified dilution")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_M_all_24h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()

image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_V_all_24h,zlim=ZLIM24comm,col=CRP1,xlab="",ylab="",log="x",main="Virus dilution")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_V_all_24h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()

image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_M_P1_24h,zlim=ZLIM24ind1,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_M_P1_24h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()
image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_V_P1_24h,zlim=ZLIM24ind1,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_V_P1_24h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()
image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_M_P2_24h,zlim=ZLIM24ind2,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_M_P2_24h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()
image.plot(x=(P2STAR/P1STAR)[1,],y=r1s/r2,BIAS_V_P2_24h,zlim=ZLIM24ind2,col=CRP1,xlab="",ylab="",log="x")
title(xlab=expression("Relative abundance P"[2]*":P"[1]),line=3,font=2)
title(ylab=expression("Relative growth rate P"[1]*":P"[2]),line=3,font.lab=2)
contour(x=(P2STAR/P1STAR)[1,],y=r1s/r2,z=BIAS_V_P2_24h, add=TRUE, levels=LVL,labcex=LBCX ,method=MET,col=ConCol)
box()


dev.copy2eps(file=paste("../../figures/PPVVG_24h",".eps",sep=""))
dev.copy2pdf(file=paste("../../figures/PPVVG_24h",".pdf",sep=""))
dev.copy(jpeg,file=paste("../../figures/PPVVG_24h",".jpg",sep=""))
dev.off()
dev.off()








