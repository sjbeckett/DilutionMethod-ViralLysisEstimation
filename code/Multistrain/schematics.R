# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#schematics.R

## contains code to create and save schematic images about:
## 1) finding viral lysis rate as the difference in slopes between the apparent growth curves of the classical and modified dilution series.
## 2) How steady state mortality is partitioned between grazing, viral lysis and niche competition.

#additional dependencies
library(fields)

##  1) Plot idealised apparent growth curves for the classical and modified dilution series to illustrate how they are used in practice

XS=seq(0,1,0.1)
CLAS = 0.7 - 0.6*XS  #INT - SLOPE
MOD = 0.9 - 0.8*XS 

plot(XS[-1],MOD[-1],xlim=c(0,1),ylim=c(0,1),xlab="Proportion of WSW (F)",ylab="Apparent growth rate (per day)",xaxs="i",yaxs="i",col="white")
points(XS[-1],CLAS[-1],pch=19,col="white")
lines(XS,CLAS,lty=2,lwd=3)
lines(XS,MOD,lwd=3)
legend("topright",c("Modified dilution","Classic dilution"),lty=c(1,2),lwd=c(3,3),bty="n")

dev.copy2eps(file="../../figures/appgrocurve.eps")
dev.copy2pdf(file="../../figures/appgrocurve.pdf")
dev.off()


##  2) Plot partitioning of steady state mortality

MAX=100 #mortality must sum to 100%
x = seq(0,MAX,length.out=800) #let x represent lysis
y = x # let y represent grazing
Len=length(x)

z=matrix(0,Len,Len)

for(aa in 1:Len){
	for(bb in 1:Len){

		z[aa,bb] = MAX - (x[aa] + y[bb])# let z represent niche competition

	}
}


z[z<0] = NA

CRP1 = colorRampPalette(c('lightblue','dark blue'), space="Lab")
CRP1=CRP1(50)


image.plot(x,y,z,xlab="% viral lysis",ylab="% grazing",legend.shrink=0.8,col=CRP1,legend.width=1.5,legend.args=list(text="% niche competition",side=4,font=1,line=2.2,cex=1.8),font=1,font.lab=1,cex.lab=1.8,bty="n")

NicheEXAMP = c(0.05,0.5,0.95)*MAX  #Example levels of niche competition to search across.


#Plotting
TXTSZ = 1.2
PS = 4
ANG=-45
C1="yellow"
lines(x, MAX-x-NicheEXAMP[1],col=C1)
text(PS,4,"95%", srt=ANG,col=C1,cex=TXTSZ)
lines(x, MAX-x-NicheEXAMP[2],col=C1)
text(PS,49,"50%", srt=ANG,col=C1,cex=TXTSZ)
lines(x, MAX-x-NicheEXAMP[3],col=C1)
text(PS,94," 5%", srt=ANG,col=C1,cex=TXTSZ)
dev.copy2pdf(file="../../figures/schematic.pdf")
dev.copy2eps(file="../../figures/schematic.eps")
dev.off()


