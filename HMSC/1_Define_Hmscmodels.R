setwd() 

library(Hmsc)

#Load Y, X and S matrices
load(file="allData_final.R")

#Exclude 15 samples with > 34 days trapping from analysis (mostly offseason samples).
#119 samples are retained.
exclude = which(X$TrapDays>34)
otuY = otuY[-exclude,]
haploY = haploY[-exclude,]
X = droplevels(X[-exclude,])
S = droplevels(S[-exclude,])


# Check for absent (0) or ubiquitous taxa (1) and exclude those observed < 5 times.

#First for OTUs:
rareotu = which(colSums(1*(otuY>0))<5)
length(rareotu)
otuY = otuY[,-rareotu]

#Then for haplotypes:
rarehaplo = which(colSums(1*(haploY>0))<5)
length(rarehaplo)
haploY = haploY[,-rarehaplo]

X$Time.Period = as.factor(X$Time.Period)
X$FullSample = as.factor(X$FullSample)
summary(X)

#A subset of five variables from X are included in the Hmsc models:
X_for_model = X[, c("Time.Period", "ForestWood", "bio1", "TrapDays", "FullSample")]
XFormula = ~Time.Period + ForestWood + bio1 + TrapDays + FullSample

# Define the midpoint between the sampling start-date and end-date.
# This is used to assess temporal autocorrelation in the data (via a temporal random effect).
S$MidPoint = X$RunDay + X$TrapDays/2

#View sampling locations
plot(S$Longitude, S$Latitude)

#define one set of coordinates per site
uniquesites = unique(data.frame(S$siteID, S$Longitude, S$Latitude))
names(uniquesites) = c("siteID", "Longitude", "Latitude")
dim(uniquesites)
# 37 sites in total.

#These site locations are coded as a spatial random effect.
#Sampling date midpoint as a temporal random effect.
# Unique sample ID is also included as a random effect.

studyDesign = data.frame(site = as.factor(S$siteID), time = as.factor(S$MidPoint), sample = as.factor(S$Sample))

xy = uniquesites[, -1]
row.names(xy) = uniquesites[,1]
rL.site = HmscRandomLevel(sData = xy, longlat = TRUE)

Time = levels(studyDesign$time)
time = data.frame(time=as.numeric(Time))
rownames(time) = Time
rL.time = HmscRandomLevel(sData = time)

Sample = studyDesign$sample
rL.sample = HmscRandomLevel(units = levels(Sample))


#Define a hurdle Hmsc model of taxa presence-absence and log(abundance) for OTUs and Haplotypes in turn (four model structures defined in total).

otuYpa = 1*(otuY>0)
otuYabu = otuY
otuYabu[otuYabu==0] = NA
otuYabu=log(otuYabu)

haploYpa = 1*(haploY>0)
haploYabu = haploY
haploYabu[haploYabu==0] = NA
haploYabu=log(haploYabu)

m1 = Hmsc(Y=otuYpa, XData = X_for_model,  XFormula = XFormula,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time, "sample" = rL.sample)})

m2 = Hmsc(Y=otuYabu, YScale = TRUE,
          XData = X_for_model,  XFormula = XFormula,
          distr="normal",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time, "sample" = rL.sample)})

m3 = Hmsc(Y=haploYpa, XData = X_for_model,  XFormula = XFormula,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time, "sample" = rL.sample)})

m4 = Hmsc(Y=haploYabu, YScale = TRUE,
          XData = X_for_model,  XFormula = XFormula,
          distr="normal",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time, "sample" = rL.sample)})


models = list(m1, m2, m3, m4)
modelnames = c("OTUpres_abs", "OTUabu_cop", "HAPpres_abs", "HAPabu_cop")
save(models,modelnames,file = "unfitted_models.Rdata")


