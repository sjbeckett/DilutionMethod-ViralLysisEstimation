# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#examplePIVG.R



##  Example simulation of bottle dynamics during diluted incubations for a Phytoplankton,Infected phytoplankton, Viral and Grazer model. Phytoplankton are susceptible to being grazed and to viral infection.


## Setting Life History Traits

#Phytoplankton
r1 = 1/24 # phytoplankton growth rate [per hour]
K = 2.2*10^7 # carrying capacity of phytoplankton [ cells per ml ]

#Viral
phi1 = 1*10^-8 # adsoprtion rate  [ ml/(virion . hr) ]
Beta1 = 50 # burst size, new [virions per cell]
omega1 = 0.02 # rate of inactivation of virus infectivity [per hour]
L1 = 2 #latent period (time between infection and lysis) [hr]
eta1 = 1/L1 #inverse latent period -  rate of lysis of the infected cell class

#Grazing
a1 = 2*10^-6   # grazing/clearance rate on phytoplankton   [ ml/(grazer . hr) ]
G0 = 1000  # grazer density [grazers per ml]


## Steady State Calculations
# These are used to set the initial conditions for the dilution experiments

	Gstar = G0
	Nstar = omega1/( phi1*((Beta1*eta1)/(eta1+a1*G0)-1) )
	Vstar = (r1*(1-Nstar/K)-a1*G0)/( phi1*(1+r1*Nstar/(K*(eta1+a1*G0))) )
	Istar = phi1*Nstar*Vstar/(eta1+a1*G0)

## Dilution Simulation Parameters

Dil_levels = seq(0.1,1,0.1)  # dilution levels (proportion of WSW F)

PROJECTINFO = list(filename="model.RData",modelused="gniv")
model_params = c(growth1= r1, carryingCap = K, AdsorptionRate1 = phi1 ,  burst1 = Beta1, inactivation1 = omega1, GrazingRate1=a1, LatentPeriod=L1)
model_state = list(c( N1 = Nstar, I1 = Istar, V1= Vstar, G = G0))[[1]]

## Running the Simulation

#Classic
GG1_b_cdil = ANALYSIS(GNIV,model_params,model_state,c(1,2),c(1,2,4),Dil_levels,0,PROJECTINFO) #classic dilution series --  V not diluted, P is diluted, G is diluted

#Modified
GG1_b_mdil = ANALYSIS(GNIV,model_params,model_state,c(1,2),c(1,2,3,4),Dil_levels,0,PROJECTINFO) # modified dilution series -- everything diluted

#Viral
GG1_b_vdil = ANALYSIS(GNIV,model_params,model_state,c(1,2),c(3),Dil_levels,0,PROJECTINFO)  #virus dilution series -- only viruses are diluted


## Plotting Outputs

#Function to find input and estimated rates from dilution experiments

FIND2<-function(GG1_b_cdil,GG1_b_mdil,GG1_b_vdil,r1,a1,G0,K,Nstar,Istar,Vstar,phi1,eta1,TIME){
#TIME -Incubation time in hours
INDEX = which( (abs(GG1_b_cdil$Timings-TIME)) == min((abs(GG1_b_cdil$Timings-TIME))) ) # find time closest to incubation time

ACTUAL_mu = r1*Nstar/(Nstar+Istar)  #growth rate (all cells)
ACTUAL_GR = a1*G0 #grazing rate (all cells)
ACTUAL_LR = eta1*Istar/(Nstar+Istar) #lysis rate (all cells)

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
data = FIND2(GG1_b_cdil,GG1_b_mdil,GG1_b_vdil,r1,a1,G0,K,Nstar,Istar,Vstar,phi,eta1,TIME)
barplot(24*data, col=colors()[c(23,89,12,616)] , border="white", font.axis=2, beside=T, legend=rownames(data), font.lab=2, args.legend=list(bty="n",ncol=2,x="top"),ylab=paste("Rates (per day)"),main="24h incubation",ylim=24*c(-0.005,0.065))
	

dev.copy2eps(file="../../figures/examplePIVG.eps")
dev.copy2pdf(file="../../figures/examplePIVG.pdf")
dev.off()









