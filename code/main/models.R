# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#Dynamical model for anlaysis

NV2G <- function(t,state,params) { #logistic growth model with grazing 
		with(as.list(c(state, params)), {
		dN1dt <-  growth1*N1*(1-((N1+N2)/carryingCap)) - AdsorptionRate1*N1*V1 - GrazingRate1*G*N1
		dN2dt <-  growth2*N2*(1-((N1+N2)/carryingCap)) - AdsorptionRate2*N2*V2 - GrazingRate2*G*N2
		dV1dt <-  burst1*AdsorptionRate1*N1*V1-AdsorptionRate1*N1*V1 - inactivation1*V1
		dV2dt <-  burst2*AdsorptionRate2*N2*V2-AdsorptionRate2*N2*V2 - inactivation2*V2
		dGdt  <-  0
		list(c(dN1dt,dN2dt,dV1dt,dV2dt,dGdt))
		})
}


GNV <- function(t,state,params) { #grazer,phytoplankton,virus
		with(as.list(c(state, params)), {
		dN1dt <-  growth1*N1*(1-((N1)/carryingCap)) - AdsorptionRate1*N1*V1 - GrazingRate1*G*N1
		dV1dt <-  burst1*AdsorptionRate1*N1*V1-AdsorptionRate1*N1*V1 - inactivation1*V1
		dGdt  <-  0
		list(c(dN1dt,dV1dt,dGdt))
		})
}

#grazing on infected and susceptible

GNIV <- function(t,state,params) { #grazer,phytoplankton,infected,virus
		with(as.list(c(state, params)), {
		dN1dt <-  growth1*N1*(1-((N1+I1)/carryingCap)) - AdsorptionRate1*N1*V1 - GrazingRate1*G*N1
		dI1dt <-  AdsorptionRate1*N1*V1 - I1/LatentPeriod  -  GrazingRate1*G*I1
		dV1dt <-  burst1*(I1/LatentPeriod)-AdsorptionRate1*N1*V1 - inactivation1*V1
		dGdt  <-  0
		list(c(dN1dt,dI1dt,dV1dt,dGdt))
		})
}
