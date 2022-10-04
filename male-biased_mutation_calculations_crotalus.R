############################################################################
# Male-biased mutation rate calculations in Crotalus
############################################################################

### Goal: examine output from pixy - pi-within and dxy-between rattlesnake
### species. Use pi and dxy to calculate da (net nucleotide differences)
### and compare distributions between autosomes and the Z chromosome.
### Use autosome and Z-linked values as input into formulas for calculating
### the ratio of male-to-female mutation rate (alpha m) and the ratio of Z-to-
### autosome mutation rate, given alpha m (uZ/uA).

### Set working directory and load dependencies-----------------------------------

setwd('./pixy_results/')

library(scales)

### Clear environment--------------------------------------------------------------

rm(list = ls())

### Read in data-------------------------------------------------------------------

# Concatenated
cpi <- read.table('pilot_analysis_v2.male.intergenic.all_pi.txt',header=T)
cdxy <- read.table('pilot_analysis_v2.male.intergenic.all_dxy.txt',header=T)

# Pi
cvpi <- cpi[which(cpi$pop=='CV'),]
copi <- cpi[which(cpi$pop=='CO'),]
capi <- cpi[which(cpi$pop=='CA'),]
crpi <- cpi[which(cpi$pop=='CR'),]

# dxy
cvco.dxy <- cdxy[which(cdxy$pop1=='CV' & cdxy$pop2=='CO'),]
cvca.dxy <- cdxy[which(cdxy$pop1=='CA' & cdxy$pop2=='CV'),]
cvcr.dxy <- cdxy[which(cdxy$pop1=='CR' & cdxy$pop2=='CV'),]
coca.dxy <- cdxy[which(cdxy$pop1=='CA' & cdxy$pop2=='CO'),]
cocr.dxy <- cdxy[which(cdxy$pop1=='CR' & cdxy$pop2=='CO'),]

### Parse autosomes and Z chromosome----------------------------------------------

# Parse pi
cvpi.auto <- cvpi[which(cvpi$chromosome!='scaffold-Z'),]
copi.auto <- copi[which(copi$chromosome!='scaffold-Z'),]
capi.auto <- capi[which(capi$chromosome!='scaffold-Z'),]
crpi.auto <- crpi[which(crpi$chromosome!='scaffold-Z'),]

cvpi.z <- cvpi[which(cvpi$chromosome=='scaffold-Z'),]
cvpi.z <- cvpi.z[which(cvpi.z$window_pos_1<106800000),]
copi.z <- copi[which(copi$chromosome=='scaffold-Z'),]
copi.z <- copi.z[which(copi.z$window_pos_1<106800000),]
capi.z <- capi[which(capi$chromosome=='scaffold-Z'),]
capi.z <- capi.z[which(capi.z$window_pos_1<106800000),]
crpi.z <- crpi[which(crpi$chromosome=='scaffold-Z'),]
crpi.z <- crpi.z[which(crpi.z$window_pos_1<106800000),]

# Parse dxy
cvca.dxy.auto <- cvca.dxy[which(cvca.dxy$chromosome!='scaffold-Z'),]
cvcr.dxy.auto <- cvcr.dxy[which(cvcr.dxy$chromosome!='scaffold-Z'),]
coca.dxy.auto <- coca.dxy[which(coca.dxy$chromosome!='scaffold-Z'),]
cocr.dxy.auto <- cocr.dxy[which(cocr.dxy$chromosome!='scaffold-Z'),]

cvca.dxy.z <- cvca.dxy[which(cvca.dxy$chromosome=='scaffold-Z'),]
cvca.dxy.z <- cvca.dxy.z[which(cvca.dxy.z$window_pos_1<106800000),]
cvcr.dxy.z <- cvcr.dxy[which(cvcr.dxy$chromosome=='scaffold-Z'),]
cvcr.dxy.z <- cvcr.dxy.z[which(cvcr.dxy.z$window_pos_1<106800000),]
coca.dxy.z <- coca.dxy[which(coca.dxy$chromosome=='scaffold-Z'),]
coca.dxy.z <- coca.dxy.z[which(coca.dxy.z$window_pos_1<106800000),]
cocr.dxy.z <- cocr.dxy[which(cocr.dxy$chromosome=='scaffold-Z'),]
cocr.dxy.z <- cocr.dxy.z[which(cocr.dxy.z$window_pos_1<106800000),]

### Calculate mean pi & Z/A ratio--------------------------------------------------

mean(cvpi.auto$avg_pi)
mean(copi.auto$avg_pi)
mean(capi.auto$avg_pi)
mean(crpi.auto$avg_pi)

mean(cvpi.z$avg_pi)
mean(copi.z$avg_pi)
mean(capi.z$avg_pi)
mean(crpi.z$avg_pi)

### Calculate mean dxy between diamondback and Western rattlesnake clades----------

mean(cvca.dxy.auto$avg_dxy)
mean(cvcr.dxy.auto$avg_dxy)
mean(coca.dxy.auto$avg_dxy)
mean(cocr.dxy.auto$avg_dxy)

mean(cvca.dxy.z$avg_dxy)
mean(cvcr.dxy.z$avg_dxy)
mean(coca.dxy.z$avg_dxy)
mean(cocr.dxy.z$avg_dxy)

### Calculate da (net nucleotide differences)---------------------------------------

# These calculations are based on the standard da formula.

# Formula for da:
# da = dxy - (piX + piY)/2
# da ~ 2uT = d (i.e., pairwise sequence divergence)

# da can be used to estimate am ('alpha m'; the ratio of male to female mutation rate based on Z and A divergence)
# am can then be used to estimate uZ/uA (Z to autosome mutation rate):
# uZ/uA = (2/3am + 1/3) / (1/2am + 1/2)

# First, calculate da between C. viridis and C. atrox, and between C. oreganus and C. atrox
# Then, calculate da between C. viridis and C. ruber, and between C. oreganus and C. ruber

# Calculate da values and plot distributions

par(mfrow=c(2,2))

cvca.da.auto <- cvca.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + capi.auto$avg_pi)/2)
cvca.da.z <- cvca.dxy.z$avg_dxy - ((cvpi.z$avg_pi + capi.z$avg_pi)/2)
boxplot(cvca.da.auto,cvca.da.z,outline=F)

coca.da.auto <- coca.dxy.auto$avg_dxy - ((copi.auto$avg_pi + capi.auto$avg_pi)/2)
coca.da.z <- coca.dxy.z$avg_dxy - ((copi.z$avg_pi + capi.z$avg_pi)/2)
boxplot(coca.da.auto,coca.da.z,outline=F)

cvcr.da.auto <- cvcr.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + crpi.auto$avg_pi)/2)
cvcr.da.z <- cvcr.dxy.z$avg_dxy - ((cvpi.z$avg_pi + crpi.z$avg_pi)/2)
boxplot(cvcr.da.auto,cvcr.da.z,outline=F)

cocr.da.auto <- cocr.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + crpi.auto$avg_pi)/2)
cocr.da.z <- cocr.dxy.z$avg_dxy - ((cvpi.z$avg_pi + crpi.z$avg_pi)/2)
boxplot(cocr.da.auto,cocr.da.z,outline=F)

# Calculate alpha m for specific comparisons:
# alphaM <- (3Z/A -2) / (4 - 3Z/A)

am.cvca <- ((3*mean(cvca.da.z))/mean(cvca.da.auto) - 2) / (4 - ((3*mean(cvca.da.z))/mean(cvca.da.auto)))
am.cvcr <- ((3*mean(cvcr.da.z))/mean(cvcr.da.auto) - 2) / (4 - ((3*mean(cvcr.da.z))/mean(cvcr.da.auto)))
am.coca <- ((3*mean(coca.da.z))/mean(coca.da.auto) - 2) / (4 - ((3*mean(coca.da.z))/mean(coca.da.auto)))
am.cocr <- ((3*mean(cocr.da.z))/mean(cocr.da.auto) - 2) / (4 - ((3*mean(cocr.da.z))/mean(cocr.da.auto)))

cvca.pi.mean <- (mean(cvpi$avg_pi)+mean(capi$avg_pi))/2
cvcr.pi.mean <- (mean(cvpi$avg_pi)+mean(crpi$avg_pi))/2
coca.pi.mean <- (mean(copi$avg_pi)+mean(capi$avg_pi))/2
cocr.pi.mean <- (mean(copi$avg_pi)+mean(crpi$avg_pi))/2

# Then take the average of autosomal and Z-linked means between comparisons

mean(cvca.da.auto)
mean(cvca.da.z)
mean(coca.da.auto)
mean(coca.da.z)

mean(cvcr.da.auto)
mean(cvcr.da.z)
mean(cocr.da.auto)
mean(cocr.da.z)

da.auto <- (mean(cvca.da.auto) + mean(coca.da.auto) + mean(cvcr.da.auto) + mean(cocr.da.auto))/4
da.z <- (mean(cvca.da.z) + mean(coca.da.z) + mean(cvcr.da.z) + mean(cocr.da.z))/4

# Calculate am
# alphaM <- (3Z/A -2) / (4 - 3Z/A)

alphaM <- ((3*da.z)/da.auto - 2) / (4 - ((3*da.z)/da.auto))
alphaM
# 2.08

# Calculate uZ/uA
# uZ/uA = (2/3am + 1/3) / (1/2am + 1/2)

uZuA <- ((0.66*alphaM)+0.33) / ((0.5*alphaM)+0.5)
uZuA
# 1.106 (fairly similar to birds!)


### Calculate da with 4-species pi average------------------------------------------

# These calculations are based an adjusted da formula to approximate ancestral
# variation more closely.

# Formula for da (modified):
# da = dxy - (piCV + piCA + piCO + piCR)/4
# da ~ 2uT = d (i.e., pairwise sequence divergence)

# da can be used to estimate am ('alpha m'; the ratio of male to female mutation rate based on Z and A divergence)
# am can then be used to estimate uZ/uA (Z to autosome mutation rate):
# uZ/uA = (2/3am + 1/3) / (1/2am + 1/2)

# First, calculate da between C. viridis and C. atrox, and between C. oreganus and C. atrox
# Then, calculate da between C. viridis and C. ruber, and between C. oreganus and C. ruber

# Calculate da values and plot distributions:

par(mfrow=c(2,2))

cvca.da.auto <- cvca.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + capi.auto$avg_pi + copi.auto$avg_pi + crpi.auto$avg_pi)/4)
cvca.da.z <- cvca.dxy.z$avg_dxy - ((cvpi.z$avg_pi + capi.z$avg_pi + copi.z$avg_pi + crpi.z$avg_pi)/4)
boxplot(cvca.da.auto,cvca.da.z,outline=F)

coca.da.auto <- coca.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + capi.auto$avg_pi + copi.auto$avg_pi + crpi.auto$avg_pi)/4)
coca.da.z <- coca.dxy.z$avg_dxy - ((cvpi.z$avg_pi + capi.z$avg_pi + copi.z$avg_pi + crpi.z$avg_pi)/4)
boxplot(coca.da.auto,coca.da.z,outline=F)

cvcr.da.auto <- cvcr.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + capi.auto$avg_pi + copi.auto$avg_pi + crpi.auto$avg_pi)/4)
cvcr.da.z <- cvcr.dxy.z$avg_dxy - ((cvpi.z$avg_pi + capi.z$avg_pi + copi.z$avg_pi + crpi.z$avg_pi)/4)
boxplot(cvcr.da.auto,cvcr.da.z,outline=F)

cocr.da.auto <- cocr.dxy.auto$avg_dxy - ((cvpi.auto$avg_pi + capi.auto$avg_pi + copi.auto$avg_pi + crpi.auto$avg_pi)/4)
cocr.da.z <- cocr.dxy.z$avg_dxy - ((cvpi.z$avg_pi + capi.z$avg_pi + copi.z$avg_pi + crpi.z$avg_pi)/4)
boxplot(cocr.da.auto,cocr.da.z,outline=F)

# Calculate alpha m for specific comparisons:
# alphaM <- (3Z/A -2) / (4 - 3Z/A)

am.cvca <- ((3*mean(cvca.da.z))/mean(cvca.da.auto) - 2) / (4 - ((3*mean(cvca.da.z))/mean(cvca.da.auto)))
am.cvcr <- ((3*mean(cvcr.da.z))/mean(cvcr.da.auto) - 2) / (4 - ((3*mean(cvcr.da.z))/mean(cvcr.da.auto)))
am.coca <- ((3*mean(coca.da.z))/mean(coca.da.auto) - 2) / (4 - ((3*mean(coca.da.z))/mean(coca.da.auto)))
am.cocr <- ((3*mean(cocr.da.z))/mean(cocr.da.auto) - 2) / (4 - ((3*mean(cocr.da.z))/mean(cocr.da.auto)))

cvca.pi.mean <- (mean(cvpi$avg_pi)+mean(capi$avg_pi))/2
cvcr.pi.mean <- (mean(cvpi$avg_pi)+mean(crpi$avg_pi))/2
coca.pi.mean <- (mean(copi$avg_pi)+mean(capi$avg_pi))/2
cocr.pi.mean <- (mean(copi$avg_pi)+mean(crpi$avg_pi))/2

# Then take the average of autosomal and Z-linked means between comparisons

mean(cvca.da.auto)
mean(cvca.da.z)
mean(coca.da.auto)
mean(coca.da.z)

mean(cvcr.da.auto)
mean(cvcr.da.z)
mean(cocr.da.auto)
mean(cocr.da.z)

da.auto <- (mean(cvca.da.auto) + mean(coca.da.auto) + mean(cvcr.da.auto) + mean(cocr.da.auto))/4
da.z <- (mean(cvca.da.z) + mean(coca.da.z) + mean(cvcr.da.z) + mean(cocr.da.z))/4

# Calculate am
# alphaM <- (3Z/A -2) / (4 - 3Z/A)

alphaM <- ((3*da.z)/da.auto - 2) / (4 - ((3*da.z)/da.auto))
alphaM
# 2.026

# Calculate uZ/uA
# uZ/uA = (2/3am + 1/3) / (1/2am + 1/2)

uZuA <- ((0.66*alphaM)+0.33) / ((0.5*alphaM)+0.5)
uZuA
# 1.10 (Very little effect of mean pi in da calculations on overall am and uZ/uA calculations)

