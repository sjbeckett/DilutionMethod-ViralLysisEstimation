# Part of the code used in:
# Beckett and Weitz. Code for: The effect of strain level diversity on robust inference of virus-induced mortality.
# 
# MIT License

#workflow.R


#CODE workflow

#1. load libraries required to run simulations and analysis. This is essential to run steps 2-8.
source("../main/LOADREQUIREMENTS.R")

#2. SCHEMATICS -- plot app. growth curves, and triangle of partitioned mortality (shown in figures 1 and 2)
source("schematics.R")

#3. EXAMPLE -- example of dilution based estimation
source("examplePVG.R")

#4. SCANNING THE TRIANGLE (P,V,G) -- show 2 hour and 24 hour incubation estimations of viral lysis
source("scanPVG.R")

#5. EXAMPLE -- example of dilution based estimation using an infected state
source("examplePIVG.R")

#6. SCANNING THE TRIANGLE (P,V,I,G) 15 minute latent period -- show 2 hour and 24 hour incubation estimations of viral lysis
source("scan_15minL_PIVG.R")

#7. SCANNING THE TRIANGLE (P,V,I,G) 4 hour latent period -- show 2 hour and 24 hour incubation estimations of viral lysis
source("scan_4hourL_PIVG.R")

#8. SCANNING THE TRIANGLE (P,V,I,G) 24 hour latent period -- show 2 hour and 24 hour incubation estimations of viral lysis
source("scan_4hourL_PIVG.R")

#9.Densities of PVG lysis rate bias
source("LHS_PVG.R")

#10.Densities of PVIG lysis rate bias when latent period is 15 minutes
source("LHS_PIVG_15.R")

#11.Densities of PVIG lysis rate bias when latent period is 4h
source("LHS_PIVG_4h.R")

#12.Densities of PVIG lysis rate bias when latent period is 24h
source("LHS_PIVG_24h.R")

#13.Diversity - concept and viral lysis rates when one phytoplankton type is relatively  low abundance, and fast growing.
source("diversity_concept.R")

#14.Diversity growth vs. abundance
source("multistrain_PPVVG.R")
