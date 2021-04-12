# setwd("hartop_emily")
localDir = "."
ModelDir = file.path(localDir, "models")
DataDir = file.path(localDir, "data")
library(Hmsc)

#Load study design (S), X covariate and mOTU (Y) data.
load(file=file.path(DataDir,"allData3.R")) #S,X,Y

# Check for absent (0) or ubiquitous mOTUs (1).
range(colMeans(Y>0))
min(colSums(Y>0))
# =1.

# Exclude mOTUs observed in < 10 samples. Leaves 162 mOTUS.
raretaxa = which(colSums(1*(Y>0))<10)
length(raretaxa)
Y = Y[,-raretaxa]

hist(colMeans(Y>0),main="prevalence")

hist(as.matrix(log(Y[Y>0])),main="log abundance conditional on presence")

# Many taxa are rare, need a zero-inflated model. 
# Set up a hurdle model: taxa presence-absence and log(abundance conditional on presence).

summary(X)

X$FullSample = as.factor(X$FullSample)

#Re-named the four sampling time periods, so that they are sequentially ordered.
X$Time.Period[which(X$Time.Period=="Late Spring")] = "1_Late_Spring"
X$Time.Period[which(X$Time.Period=="Late Summer")] = "3_Late_Summer"
X$Time.Period[which(X$Time.Period=="Midsummer")] = "2_Midsummer"
X$Time.Period[which(X$Time.Period=="Offseason")] = "4_Offseason"
X$Time.Period = as.factor(X$Time.Period)

# From among the XData, I will include bio10 (Mean temperature of the warmest quarter), %Forest+woodland, %Agriculture,
# %Built, and the categorical variable "Time.Period" as initial X covariates (fixed effects).


#Test which potential X covariates (climatic and habitat-related) correlate best with
# a 2-dimensional NMDS of the mOTU abudance data.
library(vegan)
empty = which(rowSums(Y)==0)
nmds = metaMDS(vegdist(Y[-empty,]), k = 2, trymax=100)
cortable = round(cor(data.frame(nmds$points, X[-empty, 14:40])), digits = 3)[, 1:3]

head(cortable[order(cortable[, 1]), ]);tail(cortable[order(cortable[, 1]), ])
# First ordination axis correlates most strongly with the gradient in summer temperatures
# (Mean Temperature of Warmest Quarter) across Sweden: bio10 ( R = -0.8).

head(cortable[order(cortable[, 2]), ]);tail(cortable[order(cortable[, 2]), ])
# Ordination axis 2 is not clearly related to the climatic or land cover data.

# More to sampling date/season:
plot(X$RunDay[-empty], nmds$points[,2])
cor.test(X$RunDay[-empty], nmds$points[,2])
boxplot(nmds$points[,2]~X$Time.Period[-empty])
# Offseason samples (those for which sampling began in the late fall) are distinct in composition 
#from spring to summer samples.

# I will also include the number of trapping days to assess the impact of differences in field sampling effort
# and the binary variable "Full.seq.sample" to indicate whether samples were fully sequenced or not.

XFormula = ~Time.Period + ForestWood + Agriculture + bio10 + TrapDays + FullSample

head(S)
S$MidPoint = X$RunDay + X$TrapDays/2
# Calculate midpoint (in days) between the sampling start-date and end-date for each trapping sample to
# assess temporal autocorrelation in the data (via a temporal random effect).

plot(S$Longitude, S$Latitude, asp = 1)

uniquesites = unique(data.frame(S$siteID, S$Longitude, S$Latitude))
names(uniquesites) = c("siteID", "Longitude", "Latitude")
dim(uniquesites)
# = 37 coordinate pairs are unique.

# These site locations are coded as a spatial random effect.

studyDesign = data.frame(site = as.factor(S$siteID), time = as.factor(S$MidPoint))

xy = uniquesites[, -1]
row.names(xy) = uniquesites[,1]
rL.site = HmscRandomLevel(sData = xy, longlat = TRUE)

Time = levels(studyDesign$time)
time = data.frame(time=as.numeric(Time))
rownames(time) = Time
rL.time = HmscRandomLevel(sData = time)

Ypa = 1*(Y>0)
Yabu = Y
Yabu[Yabu==0] = NA
Yabu=log(Yabu)

m1 = Hmsc(Y=Ypa, XData = X,  XFormula = XFormula,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time)})

m2 = Hmsc(Y=Yabu, YScale = TRUE,
          XData = X,  XFormula = XFormula,
          distr="normal",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time)})

##############
#Import haplotype Y data for models 5-8:

haploSXY = read.csv2("data/SXYhaplo_new.csv")
names(haploSXY)[1:14]
# after checking sample order is the same here as in S, deleted all but haplotype columns:
#unique(match(haploSXY$Sample, S$Sample)-1:134)
haploY = haploSXY[,-c(1:14)]
# As above, I have set the threshold for inclusion in the model at 10 observations (samples).

rarehaplo = which(colSums(haploY>0)<10)
length(rarehaplo)
haploY = haploY[,-rarehaplo]
#2504 haplotypes are rare by this definition. Excluding these leaves 193 haplotypes in the dataset.


haploYpa = 1*(haploY>0)
haploYabu = haploY
haploYabu[haploYabu==0] = NA
haploYabu=log(haploYabu)

##############

m3 = Hmsc(Y=haploYpa, XData = X,  XFormula = XFormula,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time)})

m4 = Hmsc(Y=haploYabu, YScale = TRUE,
          XData = X,  XFormula = XFormula,
          distr="normal",
          studyDesign=studyDesign,
          ranLevels={list("site" = rL.site, "time" = rL.time)})


models = list(m1, m2, m3, m4)
modelnames = c("OTUpres_abs", "OTUabu_cop", "HAPpres_abs", "HAPabu_cop")

save(models,modelnames,file = file.path(ModelDir, "unfitted_models"))
