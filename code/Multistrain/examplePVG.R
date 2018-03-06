# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#examplePVG.R



##  Example simulation of bottle dynamics during diluted incubations for a Phytoplankton, Viral and Grazer model. Phytoplankton are susceptible to being grazed and to viral infection.


## Setting Life History Traits

#Phytoplankton
r1 = 1/24 # phytoplankton growth rate [per hour]
K = 2.2*10^7 # carrying capacity of phytoplankton [ cells per ml ]

#Viral
phi = 1*10^-8 # adsoprtion rate  [ ml/(virion . hr) ]
Beta = 50 # burst size, new [virions per cell]
omega = 0.02 # rate of inactivation of virus infectivity [per hour]

#Grazing
a1 = 2*10^-6   # grazing/clearance rate on phytoplankton   [ ml/(grazer . hr) ]
G0 = 1000  # grazer density [grazers per ml]


## Steady State Calculations
# These are used to set the initial conditions for the dilution experiments

Gstar = G0
Nstar = omega/((Beta-1)*phi)
Vstar = (1/phi) * (r1*( 1- (omega/((Beta-1)*phi*K)) ) - a1*G0)


## Dilution Simulation Parameters

Dil_levels = seq(0.1,1,0.1)  # dilution levels (proportion of WSW F)

PROJECTINFO = list(filename="model.RData",modelused="gnv")
model_params = c(growth1= r1, carryingCap = K, AdsorptionRate1 = phi ,  burst1 = Beta, inactivation1 = omega, GrazingRate1=a1)
model_state = list(c( N1 = Nstar, V1= Vstar, G = G0))[[1]]

## Running the Simulation

#Classic
GG1_b_cdil = ANALYSIS(GNV,model_params,model_state,c(1),c(1,3),Dil_levels,0,PROJECTINFO) #classic dilution series --  V not diluted, P is diluted, G is diluted

#Modified
GG1_b_mdil = ANALYSIS(GNV,model_params,model_state,c(1),c(1,2,3),Dil_levels,0,PROJECTINFO) # modified dilution series -- everything diluted

#Viral
GG1_b_vdil = ANALYSIS(GNV,model_params,model_state,c(1),c(2),Dil_levels,0,PROJECTINFO)  #virus dilution series -- only viruses are diluted


## Plotting Outputs

#Function to find input and estimated rates from dilution experiments
FIND<-function(GG1_b_cdil,GG1_b_mdil,GG1_b_vdil,r1,a1,G0,K,Nstar,Vstar,phi1,TIME){
#TIME -Incubation time in hours
INDEX = which( (abs(GG1_b_cdil$Timings-TIME)) == min((abs(GG1_b_cdil$Timings-TIME))) ) # find time closest to incubation time

ACTUAL_mu = r1  #growth
ACTUAL_GR = a1*G0 #grazing rate
ACTUAL_LR = Vstar*phi1 #lysis rate

classic_mu = GG1_b_cdil$All_mu0[INDEX]
classic_gr = -GG1_b_cdil$All_DY[INDEX]
classic_lr = NA

modified_mu = GG1_b_mdil$All_mu0[INDEX]
modified_gr = classic_gr #classic grazing and modified grazing are the same
modified_lr = -(GG1_b_mdil$All_DY[INDEX] - GG1_b_cdil$All_DY[INDEX])
#modified_lr = GG1_b_mdil$All_mu0[INDEX]-GG1_b_cdil$All_mu0[INDEX]

virdil_mu = NA
virdil_gr = NA
virdil_lr = -GG1_b_vdil$All_DY[INDEX]

data = cbind(c(ACTUAL_mu,classic_mu,modified_mu,virdil_mu),c(ACTUAL_GR,classic_gr,modified_gr,virdil_gr),c(ACTUAL_LR,classic_lr,modified_lr,virdil_lr))
colnames(data)=c("Growth","Grazing","Lysis")
rownames(data)=c("Model input","Classical dilution","Modified dilution","Viral dilution")
return(data)
}


TIME=24 # Want to know rates after 24h incubation
data = FIND(GG1_b_cdil,GG1_b_mdil,GG1_b_vdil,r1,a1,G0,K,Nstar,Vstar,phi,TIME)
barplot(24*data, col=colors()[c(23,89,12,616)] , border="white", font.axis=2, beside=T, legend=rownames(data), font.lab=2, args.legend=list(bty="n",ncol=2,x="top"),ylab=paste("Rates (per day)"),main="24h incubation",ylim=24*c(-0.005,0.065))
	

dev.copy2eps(file="../../figures/examplePVG.eps")
dev.copy2pdf(file="../../figures/examplePVG.pdf")
dev.off()

### Apparent growth plot
F=seq(0.1,1,0.1)
G = seq(0,1,0.1)
LWD=4
INDEX = which( (abs(GG1_b_cdil$Timings-TIME)) == min((abs(GG1_b_cdil$Timings-TIME))) )


plot(F,24*GG1_b_mdil$Ind_appG[INDEX,,],xlab="Proportion of WSW (F)",ylab=expression("Apparent growth rate (d"^-1 *")"),col = colors()[12],ylim=c(0,1.05),xlim=c(0,1),lwd=LWD)
lines(G,24*GG1_b_mdil$Ind_mu0[INDEX]+G*24*GG1_b_mdil$Ind_DY[INDEX],col = colors()[12],lwd=LWD)
points(F,24*GG1_b_cdil$Ind_appG[INDEX,,],col=colors()[89],lwd=LWD)
lines(G,24*GG1_b_cdil$Ind_mu0[INDEX]+G*24*GG1_b_cdil$Ind_DY[INDEX],col = colors()[89],lwd=LWD)
points(F,24*GG1_b_vdil$Ind_appG[INDEX,,],col=colors()[616],lwd=LWD)
lines(G,24*GG1_b_vdil$Ind_mu0[INDEX]+G*24*GG1_b_vdil$Ind_DY[INDEX],col = colors()[616],lwd=LWD)


legend("topright",c("Classical dilution","Modified dilution","Viral dilution"),col=colors()[c(89,12,616)],bty="n",lty=1,lwd=LWD)


dev.copy2eps(file="../../figures/examplePVG_AppG.eps")
dev.copy2pdf(file="../../figures/examplePVG_AppG.pdf")
dev.off()

### Population plots

CC = colors()[89]
MM = colors()[12]
VV = colors()[616]
LWD=4



F=1
plot(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],xlab=expression("Cells ml" ^-1),ylab=expression("Viruses ml"^-1),xlim=c(0,93000),ylim=c(0,4250000),col=MM,type="l",lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=2
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=3
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=4
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=5
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=6
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=7
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=8
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=9
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)
arrows(GG1_b_mdil$alltimeseries[INDEX-1,F,1],GG1_b_mdil$alltimeseries[INDEX-1,F,2],GG1_b_mdil$alltimeseries[INDEX,F,1],GG1_b_mdil$alltimeseries[INDEX,F,2],length=0.1,col=MM,lwd=LWD)
F=10
lines(GG1_b_mdil$alltimeseries[1:INDEX,F,1],GG1_b_mdil$alltimeseries[1:INDEX,F,2],col=MM,lwd=LWD)



F=1
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=2
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=3
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=4
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=5
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=6
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=7
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=8
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=9
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)
arrows(GG1_b_cdil$alltimeseries[INDEX-1,F,1],GG1_b_cdil$alltimeseries[INDEX-1,F,2],GG1_b_cdil$alltimeseries[INDEX,F,1],GG1_b_cdil$alltimeseries[INDEX,F,2],length=0.1,col=CC,lwd=LWD)
F=10
lines(GG1_b_cdil$alltimeseries[1:INDEX,F,1],GG1_b_cdil$alltimeseries[1:INDEX,F,2],col=CC,lwd=LWD)


F=1
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=2
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=3
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=4
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=5
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=6
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=7
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=8
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=9
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)
arrows(GG1_b_vdil$alltimeseries[INDEX-1,F,1],GG1_b_vdil$alltimeseries[INDEX-1,F,2],GG1_b_vdil$alltimeseries[INDEX,F,1],GG1_b_vdil$alltimeseries[INDEX,F,2],length=0.1,col=VV,lwd=LWD)
F=10
lines(GG1_b_vdil$alltimeseries[1:INDEX,F,1],GG1_b_vdil$alltimeseries[1:INDEX,F,2],col=VV,lwd=LWD)

text(20000,4.1*10^6,"Classic dilution",col=CC)
text(18000,2.1*10^6,"Modified dilution",col=MM,srt=60)
text(55000,2.5*10^5,"Viral dilution",col=VV)
CX = 0.7
text(GG1_b_cdil$alltimeseries[1,1,1],GG1_b_cdil$alltimeseries[1,1,2],"F=0.1",pos=3,col=CC,cex=CX)
text(GG1_b_mdil$alltimeseries[1,1,1],GG1_b_mdil$alltimeseries[1,1,2],"F=0.1",pos=2,col=MM,cex=CX)
text(GG1_b_vdil$alltimeseries[1,1,1],GG1_b_vdil$alltimeseries[1,1,2],"F=0.1",pos=2,col=VV,cex=CX)

dev.copy2eps(file="../../figures/examplePVG_popDyn.eps")
dev.copy2pdf(file="../../figures/examplePVG_popDyn.pdf")
dev.off()



