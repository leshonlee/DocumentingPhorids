###############################################################################
# SCRIPT for: 
# DARK TAXON SCIENCE: SPECIES RICHNESS, SPATIOTEMPORAL DISTRIBUTION, 
# AND COMMUNITY COMPOSITION OF A HYPERDIVERSE GROUP (DIPTERA: PHORIDAE)
#
# Analyses of β-diversity patterns
#
###############################################################################

# Libraries ####
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(readr)
library(hms)
library(ggplot2)
library(pipeR)
library(rlist)
library(purrr)
library(ggpubr)
library(scales)
library(chron)
library(gsubfn)
library(stringi)
library(slider)
library(tidyquant)
library(tidyverse)
library(zoo)
library(GGally)
library(psych)
library(psychTools)
library(GPArotation)
library(corrplot)
library(betapart)
library(vegan)
library(geosphere)
library(janitor)
library(Matrix)
library(ggplot2)
library(ggmap)
library(sp)
library(sf)
library(sp)
library(rnaturalearth)
library(ggrepel)
library(viridis)
library(vegan)

# Importing Species Matrices ####
path_to_data <- "//xxxx/DarkTaxa/Data/Original/"
setwd(path_to_data)
speciesmat <- read.delim("./output_commatrix1point7") 
haplomat <- read.delim("./output_commatrixhaplo") 
metadata <- read.csv("./metadata - copia.csv", sep = ";") 

path_save_s <- "//xxxx/DarkTaxa/Data/BetaMatrices/Spatial/"
path_save_t <- "//xxxx/DarkTaxa/Data/BetaMatrices/Temporal/"
path_save_st <- "//xxxx/DarkTaxa/Data/BetaMatrices/SpatioTemporal/"


# Preparing them for Betadiversity analyses: ####
speciesmat2 <- speciesmat %>%
  replace(., . > 0, 1) %>% # Replace values > 0 with 1
  mutate(SampleID = paste0("sp.", row_number())) %>% 
  select(SampleID, everything()) %>% t(.) %>% as.data.frame(.) %>%
  rownames_to_column(.) %>% row_to_names(row_number = 1) %>%
  mutate(Site = paste0("Site.",substr(SampleID,2,3))) %>%
  select(Site, SampleID, everything())

haplomat2 <- haplomat %>%
  replace(., . > 0, 1) %>% # Replace values > 0 with 1
  mutate(SampleID = paste0("hap.", row_number())) %>% 
  select(SampleID, everything()) %>% t(.) %>% as.data.frame(.) %>%
  rownames_to_column(.) %>% row_to_names(row_number = 1) %>%
  mutate(Site = paste0("Site.",substr(SampleID,2,3))) %>%
  select(Site, SampleID, everything())

# SPECIES SPATIAL MATRICES = adding all info of each site. ####
speciesmatrix_spatial<-speciesmat2 %>% group_by(Site) %>%
  summarize(across(starts_with("sp"), ~ ifelse(any(. == 1), 1, 0)))

haplomatrix_spatial<-haplomat2 %>% group_by(Site) %>%
  summarize(across(starts_with("hap"), ~ ifelse(any(. == 1), 1, 0)))

# Preparing metadata dataframe to calculare average week for each time period ####
metadata2 <- metadata %>% rename(SampleID = Sample.ID) %>%
  mutate(
    Site = paste0("Site.", substr(SampleID, 2, 3)),
    Temp = case_when(
      Time.period == "Offseason" ~ 1,
      Time.period == "Late Spring" ~ 2,
      Time.period == "Midsummer" ~ 3,
      Time.period == "Late Summer" ~ 4,
      TRUE ~ NA),
    Start = as.Date(Start, format = "%Y-%m-%d"),
    End = as.Date(End, format = "%Y-%m-%d"))

# Filter for End dates in the year 2019
metadata3 <- metadata2 %>% filter(format(End, "%Y") == "2019")

# Find the latest value from 2018 for each Time category
replacement_values <- metadata2 %>%
  filter(format(End, "%Y") == "2018") %>%
  group_by(Time.period) %>%
  summarize(Latest2018 = max(End))

# Join the filtered data with the replacement values
result_data <- metadata3 %>%
  left_join(replacement_values, by = "Time.period")

# Replace End with Latest2018 if it belongs to 2019
result_data <- result_data %>%
  mutate(End = ifelse(format(End, "%Y") == "2019", Latest2018, End),
         End = as.Date(End, format = "%Y-%m-%d")) %>%
  select(-Latest2018)

# Combine the filtered and updated data with the original data
metadata_final <- bind_rows(metadata2, result_data) %>%
  filter(format(End, "%Y") == "2018")
metadata_final <- metadata_final %>% 
  mutate(
    Week1 = strftime(Start, format = "%V"),
    Week2 = strftime(End, format = "%V")
  ) %>% select(Site, Site.ID, SampleID, Time.period, Temp, Start, End, Week1, Week2, everything())

generate_week_sequence <- function(Week1, Week2) {
  return(seq(Week1, Week2))
}
# Calculate first the mean week average for each sampling
metadataf <- metadata_final %>% 
  rowwise() %>% 
  mutate(week_sequence = list(generate_week_sequence(Week1, Week2))) %>%
  mutate(WeekAvg = mean(week_sequence, na.rm = TRUE), .after = Week2)
         
# Calculate the average week number for each time period      
periodweekavg<-metadataf %>%
  group_by(Time.period) %>%
  summarize(PeriodWeekAvg = round(mean(WeekAvg, na.rm = TRUE))) 
# Joining the average week number to the main metadata
metadataf <- metadataf %>%
  left_join(periodweekavg, by = "Time.period") %>%
  select(Site, Site.ID, SampleID, Time.period, Temp, Start, End, Week1, Week2,WeekAvg, PeriodWeekAvg, everything())



# SPECIES TEMPORAL MATRICES = adding all info of sampling time ####
# Run if I want to take away site 19.
# speciesmat2<-speciesmat2 %>% 
#  dplyr::filter(Site != "Site.19")
speciesmatrix_temp<-speciesmat2 %>% left_join(metadataf %>%  
  select(SampleID, Time.period, Temp, Start, End, PeriodWeekAvg), by = "SampleID") %>%
  rename(Time = Time.period) %>%
  select(Site, SampleID, Time, Temp, Start, End, PeriodWeekAvg, everything()) %>%
  group_by(PeriodWeekAvg) %>%
  summarize(across(starts_with("sp"), ~ ifelse(any(. == 1), 1, 0)))

# Run if I want to take away site 19.
# haplomat2<-haplomat2 %>% 
#   dplyr::filter(Site != "Site.19")
haplomatrix_temp<-haplomat2 %>%
  left_join(metadataf %>%
  select(SampleID, Time.period, Temp, Start, End, PeriodWeekAvg), by = "SampleID") %>%
  rename(Time = Time.period) %>%
  select(Site, SampleID, Time, Temp, Start, End, PeriodWeekAvg, everything()) %>%
  group_by(PeriodWeekAvg) %>%
  summarize(across(starts_with("hap"), ~ ifelse(any(. == 1), 1, 0)))

# SPECIES SPATIO-TEMPORAL MATRICES = adding all info of sampling time of each site ####
# Run if I want to take away site 19.
# speciesmat2<-speciesmat2 %>% 
#  dplyr::filter(Site != "Site.19")
speciesmatrix_spatemp <- speciesmat2 %>% 
  left_join(metadataf %>% 
  select(SampleID, Time.period, Temp, Start, End, PeriodWeekAvg), by = "SampleID") %>%
  rename(Time = Time.period) %>%
  select(Site, SampleID, Time, Temp, Start, End, PeriodWeekAvg, everything()) %>%
  unite(SampleID, Site, PeriodWeekAvg, sep = "-") %>%
  group_by(SampleID) %>%
  summarize(across(starts_with("sp"), ~ ifelse(any(. == 1), 1, 0)))

# Run if I want to take away site 19.
# haplomat2<-haplomat2 %>% 
#   dplyr::filter(Site != "Site.19")
haplomatrix_spatemp<-haplomat2 %>%
  left_join(metadataf %>%
  select(SampleID, Time.period, Temp, Start, End, PeriodWeekAvg), by = "SampleID") %>%
  rename(Time = Time.period) %>%
  select(Site, SampleID, Time, Temp, Start, End, PeriodWeekAvg, everything()) %>%
  unite(SampleID, Site, PeriodWeekAvg, sep = "-") %>%
  group_by(SampleID) %>%
  summarize(across(starts_with("hap"), ~ ifelse(any(. == 1), 1, 0)))





# SPATIAL DISTANCE MATRICES: ####
SpatialDistances <- metadataf %>%
  group_by(Site) %>%
  slice(1) %>%
  select(Site, Longitude, Latitude) %>%
  column_to_rownames("Site") %>%
  distm(., fun = distHaversine)

colnames(SpatialDistances) <- unique(metadataf$Site)
rownames(SpatialDistances) <- unique(metadataf$Site)

# Data frame with the pair plots and raw spatial distance
SpatialDistances_df <-
  as.data.frame(as.table(SpatialDistances)) %>%
  rename(SpatialDist = Freq)

# Data frame with the pair plots and standarized spatial distance
SpatialDistances_dfStd <-
  as.data.frame(as.table(SpatialDistances)) %>%
  rename(SpatialDist_std = Freq) %>%
  mutate_at(c("SpatialDist_std"), ~ (scale(.) %>% as.vector))

# Distance matrix with standarized distances
SpatialDistances_std <-
  as.matrix(with(SpatialDistances_dfStd, sparseMatrix(Var1, Var2, x = SpatialDist_std)))
colnames(SpatialDistances_std) <- rownames(SpatialDistances)
rownames(SpatialDistances_std) <- rownames(SpatialDistances)



# TEMPORAL DISTANCE MATRICES: ####
TemporalDistances <- metadataf %>%
  group_by(PeriodWeekAvg) %>%
  select(PeriodWeekAvg) %>%
  distinct(PeriodWeekAvg) %>% arrange(PeriodWeekAvg)

TempDistMat<- as.matrix(dist(TemporalDistances, upper = T))
colnames(TempDistMat) <- unique(TemporalDistances$PeriodWeekAvg)
rownames(TempDistMat) <- unique(TemporalDistances$PeriodWeekAvg)




# LATITUDINAL DISTANCE MATRICES: ####
LatitudinalDistances <- metadata %>%
  mutate(Site = paste0("Site.", substr(Sample.ID, 2,3))) %>%
  group_by(Site) %>%
  slice(1) %>%
  select(Site, Latitude)

LatDistMat<- as.matrix(dist(LatitudinalDistances, upper = T))
colnames(LatDistMat) <- unique(LatitudinalDistances$Site)
rownames(LatDistMat) <- unique(LatitudinalDistances$Site)





# Save all dataframes and matrices created: ####
# Species Matrices
write.table(speciesmatrix_spatial,paste0(path_save_s, "spmat_spatial.csv"), col.names = TRUE,
            row.names = FALSE, sep = ";", dec = ",")
write.table(haplomatrix_spatial,paste0(path_save_s, "haplomat_spatial.csv"), col.names = TRUE,
            row.names = FALSE, sep = ";", dec = ",")
write.table(speciesmatrix_temp,paste0(path_save_t, "spmat_temporal.csv"), col.names = TRUE,
            row.names = FALSE, sep = ";", dec = ",")
write.table(haplomatrix_temp,paste0(path_save_t, "haplomat_temporal.csv"), col.names = TRUE,
            row.names = FALSE, sep = ";", dec = ",")

write.table(speciesmatrix_spatemp,paste0(path_save_st, "spmat_spatiotemporal.csv"), col.names = TRUE,
            row.names = FALSE, sep = ";", dec = ",")
write.table(haplomatrix_spatemp,paste0(path_save_st, "haplomat_spatiotemporal.csv"), col.names = TRUE,
            row.names = FALSE, sep = ";", dec = ",")

# Spatial and Temporal Matrices 
write.table(SpatialDistances_std,paste0(path_save_s, "spatialdistances_std.csv"), col.names = TRUE,
            row.names = TRUE, sep = ";", dec = ".")
write.table(TempDistMat,paste0(path_save_t, "temporaldistances.csv"), col.names = TRUE,
            row.names = TRUE, sep = ";", dec = ".")

# Calculating beta diversity and importing distance matrices ####
# SPECIES AND HAPLOTYPE - SPATIAL BETADIVERSITY
setwd(path_save_s)
SpBetaSpatial <-
  read.csv("./spmat_spatial.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
HapBetaSpatial <-
  read.csv("./haplomat_spatial.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")

mean(SpBetaSpatioTemporal$beta.jtu)
mean(SpBetaSpatioTemporal$beta.jne)

mean(HapBetaSpatioTemporal$beta.jtu)
mean(HapBetaSpatioTemporal$beta.jne)

# SPECIES AND HAPLOTYPE - TEMPORAL BETADIVERSITY
setwd(path_save_t)
SpBetaTemporal <-
  read.csv("./spmat_temporal.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
HapBetaTemporal<-
  read.csv("./haplomat_temporal.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")

# SPECIES AND HAPLOTYPE - SPATIO-TEMPORAL BETADIVERSITY
setwd(path_save_st)
SpBetaSpatioTemporal <-
  read.csv("./spmat_spatiotemporal.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
HapBetaSpatioTemporal<-
  read.csv("./haplomat_spatiotemporal.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")

# SPATIAL AND TEMPORAL DISTANCE MATRICES
setwd(path_save_s)
SpatialDistMatrix <-
  read.csv("./spatialdistances_std.csv", sep = ";")
setwd(path_save_t)
TemporalDistMatrix <-
  read.csv("./temporaldistances.csv", sep = ";")
colnames(TemporalDistMatrix) <- rownames(TemporalDistMatrix)

#### ---------------------------------------------------------------------- ####
#### ------------------------- CREATING DATASETS -------------------------- ####
#### ---------------------------------------------------------------------- ####

# Dataframes of Spatial Analyses for analyses ####
# Species
#splatdist <- as.data.frame(as.table(as.matrix(LatDistMat))) %>% 
#  rename(LatitudinalDistance = Freq)

spdist <- 
  as.data.frame(as.table(as.matrix(SpatialDistMatrix))) %>% 
  rename(SpatialDistance = Freq)

sptur <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatial$beta.jtu))) %>% 
  rename(Turnover = Freq)
spnes <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatial$beta.jne))) %>% 
  rename(Nestedness = Freq)
sptot <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatial$beta.jac))) %>% 
  rename(Btot = Freq)

list_df_sp = list(spdist,sptur,spnes,sptot)

Species.df <- list_df_sp %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(
    PlotPair = paste0(Var1, "-", Var2),
    DataType = "Species", 
    Var1 = as.character(Var1),
    Var2 = as.character(Var2)
  ) %>% arrange(Var1) %>%
  # Here I delete pairs that are duplicated and pairs between same plot
  filter(Var1 > Var2)

Species.df <- list_df_sp %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(
    PlotPair = paste0(Var1, "-", Var2),
    DataType = "Species") %>%
# Here I delete pairs that are duplicated and pairs between same plot
distinct(.,
         Btot,
         SpatialDistance, .keep_all = TRUE) %>%
  slice(-1)
# Haplotypes
#haplatdist <- as.data.frame(as.table(as.matrix(LatDistMat))) %>% 
#  rename(LatitudinalDistance = Freq)

hapdist <- 
  as.data.frame(as.table(as.matrix(SpatialDistMatrix))) %>% 
  rename(SpatialDistance = Freq)

haptur <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatial$beta.jtu))) %>% 
  rename(Turnover = Freq)
hapnes <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatial$beta.jne))) %>% 
  rename(Nestedness = Freq)
haptot <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatial$beta.jac))) %>% 
  rename(Btot = Freq)

list_df_hap = list(hapdist,haptur,hapnes,haptot)
Haplotype.df <- list_df_hap %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(
    PlotPair = paste0(Var1, "-", Var2),
    DataType = "Haplotypes", 
    Var1 = as.character(Var1),
    Var2 = as.character(Var2)
  ) %>% arrange(Var2) %>%
  # Here I delete pairs that are duplicated and pairs between same plot
  filter(Var1 > Var2)

# Total
SpatialSpHap <- rbind(Species.df,Haplotype.df)

#### ---------------------------------------------------------------------- ####

# Dataframes of Temporal Analyses for analyses ####
# Species
spdist <- 
  as.data.frame(as.table(as.matrix(TemporalDistMatrix))) %>% 
  rename(TemporalDistance = Freq)
sptur <- 
  as.data.frame(as.table(as.matrix(SpBetaTemporal$beta.jtu))) %>% 
  rename(Turnover = Freq)
spnes <- 
  as.data.frame(as.table(as.matrix(SpBetaTemporal$beta.jne))) %>% 
  rename(Nestedness = Freq)
sptot <- 
  as.data.frame(as.table(as.matrix(SpBetaTemporal$beta.jac))) %>% 
  rename(Btot = Freq)

list_df_sp = list(spdist,sptur,spnes,sptot)
Species.df <- list_df_sp %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(
    PlotPair = paste0(Var1, "-", Var2),
    DataType = "Species"
  ) %>%
  # Here I delete pairs that are duplicated and pairs between same plot
  distinct(.,
           Btot,
           TemporalDistance, .keep_all = TRUE) %>%
  slice(-1)

# Haplotypes
hapdist <- 
  as.data.frame(as.table(as.matrix(TemporalDistMatrix))) %>% 
  rename(TemporalDistance = Freq)
haptur <- 
  as.data.frame(as.table(as.matrix(HapBetaTemporal$beta.jtu))) %>% 
  rename(Turnover = Freq)
hapnes <- 
  as.data.frame(as.table(as.matrix(HapBetaTemporal$beta.jne))) %>% 
  rename(Nestedness = Freq)
haptot <- 
  as.data.frame(as.table(as.matrix(HapBetaTemporal$beta.jac))) %>% 
  rename(Btot = Freq)

list_df_hap = list(hapdist,haptur,hapnes,haptot)
Haplotype.df <- list_df_hap %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(
    PlotPair = paste0(Var1, "-", Var2),
    DataType = "Haplotypes"
  ) %>%
  # Here I delete pairs that are duplicated and pairs between same plot
  distinct(.,
           Btot,
           TemporalDistance, .keep_all = TRUE) %>%
  slice(-1)

# Total
TemporalSpHap <- rbind(Species.df,Haplotype.df)

#### ---------------------------------------------------------------------- ####

# Dataframes of Spatio-Temporal Analyses per Site ####
# Species
sptur <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatioTemporal$beta.jtu))) %>% 
  rename(Turnover = Freq) %>% 
  filter(substring(Var1, 1, 7) == substring(Var2, 1, 7)) %>%
  mutate(PlotCode = substring(Var1, 1, 7), .before = Var1,
         Week1 = substring(Var1, 9, 10),
         Week2 = substring(Var2, 9, 10),
         Temporal_Dist = abs(as.numeric(Week1) - as.numeric(Week2)),
         WeekPair = paste0(Week1, "-", Week2)) %>%
  select(1:3,5:7,4,8)

spnes <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatioTemporal$beta.jne))) %>% 
  rename(Nestedness = Freq) %>% 
  filter(substring(Var1, 1, 7) == substring(Var2, 1, 7)) %>%
  mutate(PlotCode = substring(Var1, 1, 7), .before = Var1)

sptot <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatioTemporal$beta.jac))) %>% 
  rename(Btot = Freq) %>% 
  filter(substring(Var1, 1, 7) == substring(Var2, 1, 7)) %>%
  mutate(PlotCode = substring(Var1, 1, 7), .before = Var1)

list_df = list(sptur,spnes,sptot)

Species.df <-
  list_df %>% reduce(inner_join, by = c("PlotCode", "Var1", "Var2")) %>%
  mutate(DataType = "Species") %>%
  mutate(WeekPairSorted = ifelse(
    substring(WeekPair, 1, 2) < substring(WeekPair, 4, 5),
    WeekPair,
    str_replace_all(paste(
      substring(WeekPair, 4, 5), "-", substring(WeekPair, 1, 2)), " ", ""))) %>%
  group_by(PlotCode, WeekPairSorted) %>%
  slice(-1) %>% ungroup() %>%  mutate_all(~ ifelse(is.nan(.), 0, .))

# Haplotypes
haptur <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatioTemporal$beta.jtu))) %>% 
  rename(Turnover = Freq) %>% 
  filter(substring(Var1, 1, 7) == substring(Var2, 1, 7)) %>%
  mutate(PlotCode = substring(Var1, 1, 7), .before = Var1,
         Week1 = substring(Var1, 9, 10),
         Week2 = substring(Var2, 9, 10),
         Temporal_Dist = abs(as.numeric(Week1) - as.numeric(Week2)),
         WeekPair = paste0(Week1, "-", Week2)) %>%
  select(1:3,5:7,4,8)

hapnes <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatioTemporal$beta.jne))) %>% 
  rename(Nestedness = Freq) %>% 
  filter(substring(Var1, 1, 7) == substring(Var2, 1, 7)) %>%
  mutate(PlotCode = substring(Var1, 1, 7), .before = Var1)

haptot <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatioTemporal$beta.jac))) %>% 
  rename(Btot = Freq) %>% 
  filter(substring(Var1, 1, 7) == substring(Var2, 1, 7)) %>%
  mutate(PlotCode = substring(Var1, 1, 7), .before = Var1)

list_df = list(haptur,hapnes,haptot)

Haplotype.df <-
  list_df %>% reduce(inner_join, by = c("PlotCode", "Var1", "Var2")) %>%
  mutate(DataType = "Haplotypes") %>%
  mutate(WeekPairSorted = ifelse(
    substring(WeekPair, 1, 2) < substring(WeekPair, 4, 5),
    WeekPair,
    str_replace_all(paste(
      substring(WeekPair, 4, 5), "-", substring(WeekPair, 1, 2)), " ", ""))) %>%
  group_by(PlotCode, WeekPairSorted) %>%
  slice(-1) %>% ungroup() %>%  mutate_all(~ ifelse(is.nan(.), 0, .))

# Total and adding latitude
Latitude <-
  metadata %>% group_by(Site.ID) %>% slice(1) %>% select(Site.ID, Latitude)

SpatioTemporalSpHap <- rbind(Species.df, Haplotype.df) %>%
  mutate(Site.ID = as.numeric(gsub("[Site.]", "", PlotCode)), .before = PlotCode) %>%
  left_join(Latitude, by = "Site.ID")

#### ---------------------------------------------------------------------- ####

# Dataframes of Spatio-Temporal Analyses for all Site Combinations ####
# Species
sptur <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatioTemporal$beta.jtu))) %>% 
  rename(Turnover = Freq) %>% 
  filter(substring(Var1, 9, 10) == substring(Var2, 9, 10)) %>%
  mutate(Week1 = substring(Var1, 9, 10),
         Week2 = substring(Var2, 9, 10),
         Temporal_Dist = abs(as.numeric(Week1) - as.numeric(Week2)),
         WeekPair = paste0(Week1, "-", Week2),
         Site1 = substring(Var1, 1, 7),
         Site2 = substring(Var2, 1, 7),
         SitePair = paste0(Site1, "-", Site2)) %>%
  select(8,9,10,4,5,7,1,2,3,6)

spnes <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatioTemporal$beta.jne))) %>% 
  rename(Nestedness = Freq)
sptot <- 
  as.data.frame(as.table(as.matrix(SpBetaSpatioTemporal$beta.jac))) %>% 
  rename(Btot = Freq)

list_df = list(sptur,spnes,sptot)

Species.df <-
  list_df %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(DataType = "Species") %>%
  mutate(WeekPairSorted = ifelse(
    substring(WeekPair, 1, 2) < substring(WeekPair, 4, 5),
    WeekPair,
    str_replace_all(paste(
      substring(WeekPair, 4, 5), "-", substring(WeekPair, 1, 2)), " ", ""))) %>%
  mutate_all( ~ ifelse(is.nan(.), 0, .))

# Adding spatial distance
spdist <- as.data.frame(as.table(as.matrix(SpatialDistMatrix))) %>%
  rename(Spatial_Dist = Freq, Site1 = Var1, Site2 = Var2)

Species.df <- Species.df %>%
  left_join(spdist, Spatial_dist, by = c("Site1", "Site2") )

# Haplotypes
haptur <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatioTemporal$beta.jtu))) %>% 
  rename(Turnover = Freq) %>% 
  filter(substring(Var1, 9, 10) == substring(Var2, 9, 10)) %>%
  mutate(Week1 = substring(Var1, 9, 10),
         Week2 = substring(Var2, 9, 10),
         Temporal_Dist = abs(as.numeric(Week1) - as.numeric(Week2)),
         WeekPair = paste0(Week1, "-", Week2),
         Site1 = substring(Var1, 1, 7),
         Site2 = substring(Var2, 1, 7),
         SitePair = paste0(Site1, "-", Site2)) %>%
  select(8,9,10,4,5,7,1,2,3,6)

hapnes <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatioTemporal$beta.jne))) %>% 
  rename(Nestedness = Freq)
haptot <- 
  as.data.frame(as.table(as.matrix(HapBetaSpatioTemporal$beta.jac))) %>% 
  rename(Btot = Freq)

list_df = list(haptur,hapnes,haptot)

Haplotype.df <-
  list_df %>% reduce(inner_join, by = c("Var1", "Var2")) %>%
  mutate(DataType = "Haplotypes") %>%
  mutate(WeekPairSorted = ifelse(
    substring(WeekPair, 1, 2) < substring(WeekPair, 4, 5),
    WeekPair,
    str_replace_all(paste(
      substring(WeekPair, 4, 5), "-", substring(WeekPair, 1, 2)), " ", ""))) %>%
  mutate_all( ~ ifelse(is.nan(.), 0, .))

# Adding spatial distance
hapdist <- as.data.frame(as.table(as.matrix(SpatialDistMatrix))) %>%
  rename(Spatial_Dist = Freq, Site1 = Var1, Site2 = Var2)

Haplotype.df <- Haplotype.df %>%
  left_join(hapdist, Spatial_dist, by = c("Site1", "Site2"))

# Total and adding latitude
splatdist <- as.data.frame(as.table(as.matrix(LatDistMat))) %>% 
  rename(LatitudinalDistance = Freq) %>% mutate(SitePair = paste0(Var1, "-", Var2)) %>%
  select(4,3)

SpatioTemporalSpHap <- rbind(Species.df,Haplotype.df) %>%
  left_join(splatdist, LatitudinalDistance, by = "SitePair") %>%
  select(-10)

SpatioTemporalSpHap %>% group_by(DataType, WeekPair) %>%
  summarize(t = mean(Turnover),
            n = mean(Nestedness) ,
            tt = mean(Btot))
mean(SpatioTemporalSpHap$Turnover)
#### ---------------------------------------------------------------------- ####

#### ---------------------------------------------------------------------- ####
#### ------------------------- PLOTTING DATASETS -------------------------- ####
#### ---------------------------------------------------------------------- ####

#### PLOTTING REGRESSION LINES - SPATIAL ####
turnover.plot <-
  SpatialSpHap  %>% 
  ggplot(aes(x=SpatialDistance, y=Turnover, color=DataType)) + 
  xlab("") + 
  ylab("Turnover") +
  geom_point() +
  geom_smooth(method = "glm") +
  scale_fill_manual(values = c("#FFC20A","#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A","#0C7BDC"))  +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_grid(~DataType, scales="free", space="free_x")

nestedness.plot <-
  SpatialSpHap  %>% 
  ggplot(aes(x=SpatialDistance, y=Nestedness, color=DataType)) + 
  xlab("") + 
  ylab("Nestedness") +
  geom_point() +
  geom_smooth(method = "glm") +
  scale_fill_manual(values = c("#FFC20A","#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A","#0C7BDC"))  +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_grid(~DataType, scales="free", space="free_x")

totalbeta.plot <-
  SpatialSpHap  %>% 
  ggplot(aes(x=SpatialDistance, y=Btot, color=DataType)) + 
  xlab("Std Euclidean Distance") + 
  ylab("Total Betadiversity") +
  geom_point() +
  geom_smooth(method = "glm") +
  scale_fill_manual(values = c("#FFC20A","#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A","#0C7BDC"))  +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_grid(~DataType, scales="free", space="free_x") 
  geom_label_repel(aes(label = ifelse(Btot > 0.96, as.character(PlotPair), "")), 
                   box.padding = 0.5, 
                   segment.size = 0.2, 
                   stat = "unique", 
                   size = 3, 
                   max.overlaps = 1000)

ggarrange(turnover.plot,nestedness.plot,totalbeta.plot, ncol = 1) 


#### PLOTTING REGRESSION LINES - TEMPORAL ####
turnover.plot <-
  TemporalSpHap  %>% 
  ggplot(aes(x=TemporalDistance, y=Turnover, color=DataType)) + 
  xlab("") + 
  ylab("Turnover") +
  geom_point() +
  geom_smooth(method = "glm") +
  scale_fill_manual(values = c("#FFC20A","#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A","#0C7BDC"))  +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_grid(~DataType, scales="free", space="free_x")

nestedness.plot <-
  TemporalSpHap  %>% 
  ggplot(aes(x=TemporalDistance, y=Nestedness, color=DataType)) + 
  xlab("") + 
  ylab("Nestedness") +
  geom_point() +
  geom_smooth(method = "glm") +
  scale_fill_manual(values = c("#FFC20A","#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A","#0C7BDC"))  +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_grid(~DataType, scales="free", space="free_x")

totalbeta.plot <-
  TemporalSpHap  %>% 
  ggplot(aes(x=TemporalDistance, y=Btot, color=DataType)) + 
  xlab("Temporal Distance (Weeks)") + 
  ylab("Total Betadiversity") +
  geom_point() +
  geom_smooth(method = "glm") +
  scale_fill_manual(values = c("#FFC20A","#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A","#0C7BDC"))  +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_grid(~DataType, scales="free", space="free_x") 
geom_label_repel(aes(label = ifelse(Btot > 0.96, as.character(PlotPair), "")), 
                 box.padding = 0.5, 
                 segment.size = 0.2, 
                 stat = "unique", 
                 size = 3, 
                 max.overlaps = 1000)

ggarrange(turnover.plot,nestedness.plot,totalbeta.plot, ncol = 1) 


#### PLOTTING REGRESSION LINES - SPATIO-TEMPORAL: Plots compared by Time Period ####
slopes_tur <- SpatioTemporalSpHap %>%
  group_by(PlotCode) %>%
  summarize(SlopeTur = coef(lm(Turnover ~ Temporal_Dist))["Temporal_Dist"])

slopes_nes <- SpatioTemporalSpHap %>%
  group_by(PlotCode) %>% 
  summarize(SlopeNes = coef(lm(Nestedness ~ Temporal_Dist))["Temporal_Dist"])

slopes_tot <- SpatioTemporalSpHap %>%
  group_by(PlotCode) %>% 
  summarize(SlopeTot = coef(lm(Btot ~ Temporal_Dist))["Temporal_Dist"])

df <- left_join(SpatioTemporalSpHap, slopes_tur, by = "PlotCode") %>%
  left_join(., slopes_nes, by = "PlotCode") %>% left_join(., slopes_tot, by = "PlotCode") %>%
  filter(!Site.ID == "19")

# Taking away point 19: 
SpatioTemporalSpHap <- SpatioTemporalSpHap %>% filter(!Site.ID == "19") 

# Plots of betadiversity colored by latitude (blue = south, red = north).
# Temporal distance is in week.
colors <- colorRampPalette(c("darkblue","white","darkred"))(36)
colors <- colorRampPalette("Greens")(37)
colors <- viridis(option = "C", 37)
turnover.plot <-
  SpatioTemporalSpHap  %>% 
  ggplot(aes(x=Temporal_Dist, y=Turnover, color=as.factor(Latitude))) + 
  xlab("Temporal Distance (Weeks)") + 
  ylab("Turnover") +
  geom_point(size = 4, alpha = .60) +
  geom_smooth(method = "glm", alpha = .10, se = T) +
  scale_color_manual(values = colors)  +
  ylim(c(0,1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")
turnover.plot + facet_wrap(DataType~ifelse(df$SlopeTur > 0, "Positive Slope", "Negative Slope"))

nestedness.plot <-
  SpatioTemporalSpHap  %>% 
  ggplot(aes(x=Temporal_Dist, y=Nestedness, color=as.factor(Latitude))) + 
  xlab("Temporal Distance (Weeks)") + 
  ylab("Nestedness") +
  geom_point(size = 4, alpha = .60) +
  geom_smooth(method = "glm", alpha = .10, se = T) +
  #scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)  +
  ylim(c(0,1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")
nestedness.plot + facet_wrap(DataType~ifelse(df$SlopeNes > 0, "Positive Slope", "Negative Slope"))

totalbeta.plot <-
  SpatioTemporalSpHap  %>% 
  ggplot(aes(x=Temporal_Dist, y=Btot, color=as.factor(Latitude))) + 
  xlab("Temporal Distance (Weeks)") + 
  ylab("Total Betadiversity") +
  geom_point(size = 4, alpha = .60) +
  geom_smooth(method = "glm", alpha = .10, se = F) +
  scale_color_manual(values = colors)  +
  ylim(c(0.5,1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") 
totalbeta.plot + facet_wrap(DataType~ifelse(df$SlopeTot > 0, "Positive Slope", "Negative Slope"))


#### PLOTTING REGRESSION LINES - SPATIO-TEMPORAL: Plots by TimePeriod ####
# Taking away point 19, changing TimePeriod names, removing duplicates
SpatioTemporalSpHap <-
  SpatioTemporalSpHap %>% mutate(
    TimePeriod = case_when(
      Week1 == 21 ~ "Late Spring",
      Week1 == 25 ~ "Midsummer",
      Week1 == 31 ~ "Late Summer",
      Week1 == 44 ~ "Offseason",
      TRUE ~ NA
    )
  ) %>% filter(!Site1 == "Site.19") %>%
  distinct(., Turnover, Spatial_Dist, TimePeriod, DataType, .keep_all = TRUE) %>%
  slice(-1)

SpatioTemporalSpHap_filtered <- subset(SpatioTemporalSpHap, Site1 != Site2)

# Plots of betadiversity colored by Time Period
# Std Spatial distance. Plots are only compared between same time periods
colors <- colorRampPalette(c("darkblue","white","darkred"))(37)
colors <- colorRampPalette("Greens")(37)
colors <- viridis(option = "C", 4)

turnover.plot <-
  SpatioTemporalSpHap_filtered  %>% 
  mutate(TimePeriod = fct_relevel(TimePeriod, "Late Spring", "Midsummer", "Late Summer", "Offseason")) %>%
  ggplot(aes(x=Spatial_Dist, y=Turnover, color = as.factor(TimePeriod))) + 
  xlab("Std Spatial Distance") + 
  ylab("Turnover") +
  geom_point(size = 4, alpha = .30) +
  geom_smooth(method = "glm", alpha = .50, se = T) +
  scale_color_manual(values = colors)  +
  ylim(c(0,1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_wrap(DataType~TimePeriod, ncol = 4)

nestedness.plot <-
  SpatioTemporalSpHap_filtered  %>% 
  mutate(TimePeriod = fct_relevel(TimePeriod, "Late Spring", "Midsummer", "Late Summer", "Offseason")) %>% 
  ggplot(aes(x=Spatial_Dist, y=Nestedness, color=as.factor(TimePeriod))) + 
  xlab("Std Spatial Distance") + 
  ylab("Nestedness") +
  geom_point(size = 4, alpha = .30) +
  geom_smooth(method = "glm", alpha = .10, se = T) +
  scale_color_manual(values = colors)  +
  ylim(c(0,1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")  +
  facet_wrap(DataType~TimePeriod, ncol = 4)

totalbeta.plot <-
  SpatioTemporalSpHap_filtered  %>% 
  mutate(TimePeriod = fct_relevel(TimePeriod, "Late Spring", "Midsummer", "Late Summer", "Offseason")) %>% 
  ggplot(aes(x=Spatial_Dist, y=Btot, color=as.factor(TimePeriod))) + 
  xlab("Std Spatial Distance") + 
  ylab("Total Betadiversity") +
  geom_point(size = 4, alpha = .30) +
  geom_smooth(method = "glm", alpha = .10, se = F) +
  scale_color_manual(values = colors)  +
  ylim(c(0,1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none") +
  facet_wrap(DataType~TimePeriod, ncol = 4)

#### ---------------------------------------------------------------------- ####

#### BOXPLOTS FOR EACH PLOT - SPATIO-TEMPORAL: Plots by TimePeriod ####
# Turnover
SpatioTemporalSpHap %>%
  ggplot(., aes(
    PlotCode,
    y = Turnover,
    color = as.factor(Latitude)
  )) +
  #geom_violin(alpha = .4) +
  #geom_point() +
  geom_boxplot() +
  coord_flip() +
  facet_wrap( ~ DataType) +
  scale_color_manual(values = colors)  +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")
# Nestedness
SpatioTemporalSpHap %>%
  ggplot(., aes(
    PlotCode,
    y = Nestedness,
    color = as.factor(Latitude)
  )) +
  #geom_violin(alpha = .4) +
  #geom_point() +
  geom_boxplot() +
  coord_flip() +
  facet_wrap( ~ DataType) +
  scale_color_manual(values = colors)  +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")
# Betadiversity total
SpatioTemporalSpHap %>%
  ggplot(., aes(
    PlotCode,
    y = Btot,
    color = as.factor(Latitude)
  )) +
  #geom_violin(alpha = .4) +
  #geom_point() +
  geom_boxplot() +
  coord_flip() +
  facet_wrap( ~ DataType) +
  scale_color_manual(values = colors)  +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")

colors <- viridis(option = "C", 37)
#### BOXPLOTS FOR EACH PLOT - SPATIO-TEMPORAL: Plots by TimePeriod ####
# Turnover
SpatioTemporalSpHap_filtered %>%
  ggplot(., aes(
    TimePeriod,
    y = Turnover,
    color = as.factor(TimePeriod)
  )) +
  #geom_violin(alpha = .4) +
  #geom_point() +
  geom_boxplot() +
  coord_flip() +
  facet_wrap( ~ DataType) +
  scale_color_manual(values = colors)  +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")
# Nestedness
SpatioTemporalSpHap %>%
  ggplot(., aes(
    TimePeriod,
    y = Nestedness,
    color = as.factor(TimePeriod)
  )) +
  #geom_violin(alpha = .4) +
  #geom_point() +
  geom_boxplot() +
  coord_flip() +
  facet_wrap( ~ DataType) +
  scale_color_manual(values = colors)  +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")
# Betadiversity total
SpatioTemporalSpHap %>%
  ggplot(., aes(
    TimePeriod,
    y = Btot,
    color = as.factor(TimePeriod)
  )) +
  #geom_violin(alpha = .4) +
  #geom_point() +
  geom_boxplot() +
  coord_flip() +
  facet_wrap( ~ DataType) +
  scale_color_manual(values = colors)  +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light(base_size = 16, base_family = "serif") +
  theme(axis.text = element_text(size = 12, colour = "black"),
        strip.text = element_text(colour = 'black'),
        legend.position = "none")

colors <- viridis(option = "C", 4)


#### ---------------------------------------------------------------------- ####

#### ---------------------------------------------------------------------- ####
#### --------------------------   MANTEL TESTS   -------------------------- ####
#### ---------------------------------------------------------------------- ####

# Spatial Analyses ####
mantel(SpBetaSpatial$beta.jtu, SpatialDistMatrix)
mantel(SpBetaSpatial$beta.jne, SpatialDistMatrix)
mantel(SpBetaSpatial$beta.jac, SpatialDistMatrix)
mantel(HapBetaSpatial$beta.jtu, SpatialDistMatrix)
mantel(HapBetaSpatial$beta.jne, SpatialDistMatrix)
mantel(HapBetaSpatial$beta.jac, SpatialDistMatrix)

# Temporal Analyses ####
mantel(SpBetaTemporal$beta.jtu, TemporalDistMatrix)
mantel(SpBetaTemporal$beta.jne, TemporalDistMatrix)
mantel(SpBetaTemporal$beta.jac, TemporalDistMatrix)
mantel(HapBetaTemporal$beta.jtu, TemporalDistMatrix)
mantel(HapBetaTemporal$beta.jne, TemporalDistMatrix)
mantel(HapBetaTemporal$beta.jac, TemporalDistMatrix)

# Spatio-temporal Analyses ####
SpatioTemporalSpHap <-
  SpatioTemporalSpHap %>% mutate(
    TimePeriod = case_when(
      Week1 == 21 ~ "Late Spring",
      Week1 == 25 ~ "Midsummer",
      Week1 == 31 ~ "Late Summer",
      Week1 == 44 ~ "Offseason",
      TRUE ~ NA
    )
  ) %>% filter(!Site1 == "Site.19")



## Species Late Spring ####
spdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptur_latespring <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

spnes_latespring <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptot_latespring <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(spdist, sptur_latespring)
mantel(spdist, spnes_latespring)
mantel(spdist, sptot_latespring)

## Species Midsummer ####
spdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptur_midsummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

spnes_midsummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptot_midsummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(spdist, sptur_midsummer)
mantel(spdist, spnes_midsummer)
mantel(spdist, sptot_midsummer)
## Species Late Summer ####
spdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptur_latesummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

spnes_latesummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptot_latesummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(spdist, sptur_latesummer)
mantel(spdist, spnes_latesummer)
mantel(spdist, sptot_latesummer)



## Species Offseason ####
spdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptur_offseason <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

spnes_offseason <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

sptot_offseason <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Species" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(spdist, sptur_offseason)
mantel(spdist, spnes_offseason)
mantel(spdist, sptot_offseason)



## Haplotypes Late Spring ####
hapdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptur_latespring <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

hapnes_latespring <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptot_latespring <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Spring") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(hapdist, haptur_latespring)
mantel(hapdist, hapnes_latespring)
mantel(hapdist, haptot_latespring)

## Haplotypes Midsummer ####
hapdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptur_midsummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

hapnes_midsummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptot_midsummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Midsummer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(hapdist, haptur_midsummer)
mantel(hapdist, hapnes_midsummer)
mantel(hapdist, haptot_midsummer)
## Haplotypes Late Summer ####
hapdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptur_latesummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

hapnes_latesummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptot_latesummer <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Late Summer") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(hapdist, haptur_latesummer)
mantel(hapdist, hapnes_latesummer)
mantel(hapdist, haptot_latesummer)



## Haplotypes Offseason ####
hapdist <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Spatial_Dist) %>%
  pivot_wider(names_from = Site2, values_from = Spatial_Dist) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptur_offseason <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Turnover) %>%
  pivot_wider(names_from = Site2, values_from = Turnover) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

hapnes_offseason <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Nestedness) %>%
  pivot_wider(names_from = Site2, values_from = Nestedness) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

haptot_offseason <-
  SpatioTemporalSpHap %>% dplyr::filter(DataType == "Haplotypes" &
                                          TimePeriod == "Offseason") %>%
  mutate(Site1 = as.character(Site1),
         Site2 = as.character(Site2)) %>% #arrange(Site2) %>%
  select(Site1, Site2, Btot) %>%
  pivot_wider(names_from = Site2, values_from = Btot) %>%
  column_to_rownames("Site1")  %>%
  as.dist()

mantel(hapdist, haptur_offseason)
mantel(hapdist, hapnes_offseason)
mantel(hapdist, haptot_offseason)





#### ---------------------------------------------------------------------- ####



#### ---------------------------------------------------------------------- ####
#### ------------------------- MODELING DATASETS -------------------------- ####
#### ---------------------------------------------------------------------- ####


# CHECKING NORMALITY ####
# Kilpisjarvi + Varanger All data
names(SpatialSpHap)
vm<-c(4:7)
par(mfrow = c(2,3))
for (i in vm){
  hist(SpatialSpHap[, i],
       xlab = colnames(SpatialSpHap[i]),
       main = colnames(SpatialSpHap[i]),
       ylab = "Freq")
}

names(SpatialSpHap[, c(4:7)])
betakv <- SpatialSpHap[, c(4:7)] %>% keep(is.numeric) %>%
  gather() %>%
  mutate(name = fct_relevel(key, names(SpatialSpHap[, c(4:7)]))) %>%
  ggplot(aes(value)) +
  facet_wrap(~ name, nrow = 2, scales = "free") +
  labs(title = "BetaDiversity Kilpisjarvi and Varanger",
       y = "density",
       x = "values") +
  geom_density(color = "blue",
               fill = "red",
               alpha = 0.5) +
  geom_point(aes(x = value, y = 0), alpha = 0.5)

# Based on the distribution of Turnover, Nestedness and Beta total:
#  - I transform the variables. 
#  - First I try with the formula: logit(Y) = log(Y/(1-Y))
#  - The distribution of response variables seems to be good now! 
#  - However, there are 0s in the data so I need to transform them.
#     >Number of 0s in Turnover and Nestedness and transformation to avoid -Inf
#         colSums(VK_total[,4:6]==0) # Number of 0s
#     
#     > Turnover:
#         sort(unique(VK_total$Turnover))[1:5] 
#         log(0.1176471/(1-0.1176471)) # -2.014903
#         Chosen transformation closest to minimum value that is not 0: 0.1
#         log(0.1/(1-0.1)) # -2.197225
#
#     >Nestedness:
#         sort(unique(VK_total$Nestedness))[1:5]  
#         log(0.0002749009/(1-0.0002749009)) # -8.198825
#         Chosen transformation closest to minimum value that is not 0: 0.0001
#         log(0.0001/(1-0.0001)) # -9.21024

colnames(SpatialSpHap)
colSums(SpatialSpHap[,5:7]==0)

sort(unique(VK_total_trans$Turnover_logit))[1:10] 
log(0.0001239157 /(1-0.0001239157 )) # -2.014903 -8.995785
#Chosen transformation closest to minimum value that is not 0: 0.1
log(0.0001/(1-0.0001)) # -2.197225


VK_total_trans <- SpatialSpHap %>%
  mutate(across(Turnover, ~ replace(., . ==  1 , 0.9))) %>%
  mutate(across(Nestedness, ~ replace(., . ==  1 , 0.9))) %>%
  mutate(across(Btot, ~ replace(., . ==  1 , 0.9))) %>%
  mutate(across(Nestedness, ~ replace(., . ==  0 , 0.0001))) %>%
  mutate(Turnover_logit = log(Turnover/(1-Turnover)),
         Nestedness_logit = log(Nestedness/(1-Nestedness)),
         Btotal_logit = log(Btot/(1-Btot)))

colnames(VK_total_trans)
betakv <- VK_total_trans[, c(5:7, 10:12)] %>% keep(is.numeric) %>%
  gather() %>%
  mutate(name = fct_relevel(key, names(VK_total_trans[, c(5:7, 10:12)]))) %>%
  ggplot(aes(value)) +
  facet_wrap(~ name, nrow = 2, scales = "free") +
  labs(title = "BetaDiversity Kilpisjarvi and Varanger",
       y = "density",
       x = "values") +
  geom_density(color = "blue",
               fill = "red",
               alpha = 0.5) +
  geom_point(aes(x = value, y = 0), alpha = 0.5)

# MODELS: ####
#   Turnover ~ SamplingSite * Taxon * Elevational Diff + (1|PlotPair) 
#   Nestedness ~ SamplingSite * Taxon * Elevational Diff + (1|PlotPair) 
#   Total Betadiversity ~ SamplingSite * Taxon * Elevational Diff + (1|PlotPair) 
#
unique(VK_total_trans$Turnover_logit)
Jtur_Alt_m <- glmmTMB(
  Turnover_logit ~ DataType * SpatialDistance + (1|PlotPair),
  data = VK_total_trans,
  family = gaussian()
)
summary(Jtur_Alt_m)

Jnes_Alt_m <- glmmTMB(
  Nestedness_logit ~ DataType * SpatialDistance + (1|PlotPair) ,
  data = VK_total_trans,
  family = poisson()
)
summary(Jnes_Alt_m)

Jtot_Alt_m <- glmmTMB(
  Btotal_logit ~ DataType * SpatialDistance + (1|PlotPair),
  data = VK_total_trans,
  family = poisson()
)
summary(Jtot_Alt_m)
Anova(Jtot_Alt_m, type = 3)
check_model(Jtot_Alt_m)
# MODELS RESULTS EXPORT ####
library(sjPlot)
library(sjmisc)
library(sjlabelled)
# All models joint results
tab_model(Jtur_Alt_m, Jnes_Alt_m, Jtot_Alt_m,
          dv.labels = c(
            "Turnover~Elevation Difference",
            "Nestedness~Elevation Difference",
            "TotalBeta~Elevation Difference"
          ), file = "ElevationalDifferenceModels_all.doc")
# Individual models results
tab_model(Jnes_Alt_m,
          dv.labels = c("Nestedness~Elevational Difference"),
          file = "ElevationalDifferenceModels_Nestedness.doc")
tab_model(Jtur_Alt_m,
          dv.labels = c("Turnover~Elevational Difference"),
          file = "ElevationalDifferenceModels_Turnover.doc")
tab_model(Jtot_Alt_m,
          dv.labels = c("TotalBeta~Elevational Difference"),
          file = "ElevationalDifferenceModels_BetaDiversity.doc")


#### ---------------------------------------------------------------------- ####


# Map of the sampling sites in Sweden: ####

# Getting the Swedish map:
sweden_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin == "Sweden")

# Create a ggplot with the spatial data
# Sampling plot area map
MapSites <- ggplot() +
  geom_sf(data = sweden_map,
          fill = "white",
          color = "black") +
  geom_point(
    data = metadata,
    aes(x = Longitude, y = Latitude),
    color = "red",
    size = 3
  ) + geom_text_repel(
    data = metadata,
    aes(x = Longitude, y = Latitude, label = Site.ID),
    nudge_x = 0.1,
    nudge_y = 0.2,
    stat = "unique",
    max.overlaps = 1000
  )

#### ---------------------------------------------------------------------- ####

###
###
###
#
# SPATIAL SPECIES AND HAPLOTIPE BETADIVERSITY
setwd(path_save_s)
SpBetaSpatial <-
  read.csv("./spmat_spatial_19.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
HapBetaSpatial <-
  read.csv("./haplomat_spatial_19.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
# TEMPORAL SPECIES AND HAPLOTIPE BETADIVERSITY
setwd(path_save_t)
SpBetaTemporal <-
  read.csv("./spmat_temporal_19.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
HapBetaTemporal<-
  read.csv("./haplomat_temporal_19.csv", row.names = 1, sep = ";") %>%
  beta.pair(., index.family = "jaccard")
# SPATIAL AND TEMPORAL DISTANCE MATRICES
setwd(path_save_s)
SpatialDistMatrix <-
  read.csv("./spatialdistances_std_19.csv", sep = ";")
setwd(path_save_t)
TemporalDistMatrix <-
  read.csv("./temporaldistances.csv", sep = ";")
colnames(TemporalDistMatrix) <- rownames(TemporalDistMatrix)

vm<-c("beta.jtu","beta.jne","beta.jac")
par(mfrow = c(3,4))
for (i in vm){
  plot(as.dist(SpatialDistMatrix), as.dist(SpBetaSpatial[[i]]),ylab = paste("Species Beta Spatial",i))
  plot(as.dist(TemporalDistMatrix), as.dist(SpBetaTemporal[[i]]),ylab = paste("Species Beta Temporal",i))
  plot(as.dist(SpatialDistMatrix), as.dist(HapBetaSpatial[[i]]),ylab = paste("Haplo Beta Spatial",i))
  plot(as.dist(TemporalDistMatrix), as.dist(HapBetaTemporal[[i]]),ylab = paste("Haplo Beta Temporal",i))
}
