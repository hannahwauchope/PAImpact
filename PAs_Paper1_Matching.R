##Before/After Protection Timing Analysis
#### Find the high quality sites
options(repos = c(CRAN = "http://cran.rstudio.com"))

cluster <- FALSE

if(cluster==TRUE){
  install.packages("rjags")
  library(rjags)
  library(reshape2, lib.loc="/home/hsw34/my-R-libs/")
  library(data.table, lib.loc="/home/hsw34/my-R-libs/")
  library(pbmcapply, lib.loc="/home/hsw34/my-R-libs/")
  library(tidyverse, lib.loc="/home/hsw34/my-R-libs/")
  library(plyr, lib.loc="/home/hsw34/my-R-libs/")
  library(stringr, lib.loc="/home/hsw34/my-R-libs/")
  library(pbapply, lib.loc="/home/hsw34/my-R-libs/")
  library(StatMatch, lib.loc="/home/hsw34/my-R-libs/")
  
  ResultsFP <- "/rds/user/hsw34/hpc-work/results/c3/"
  DataFP <- "/rds/user/hsw34/hpc-work/data/"
  ncores <- 32 
} else {
  library(ggmap)
  library(data.table)
  library(geosphere)
  library(maptools)
  library(maps)
  library(mapproj)
  library(reshape2)
  library(rgdal)
  library(rgeos)
  library(sp)
  library(raster)
  library(PBSmapping)
  library(foreign)
  library(pbmcapply)
  library(mgcv)
  library(ncdf4)
  library(chron)
  library(lattice)
  library(LaF)
  library(MASS)
  library(pbapply)
  library(tidyverse)
  library(usdm)
  library(MatchIt)
  library(e1071)
  library(fields)
  library(plyr)
  library(resample)
  library(mice)
  library(StatMatch)
  library(lme4)
  
  ResultsFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter3/Analysis/"
  DataFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Data/"
  ncores <- 6 
}

#Test changing something

load(file=paste0(DataFP, "WaterbirdData_Tatsuya/FullDataSet_Edits/Hannah_Consolidation/Spec6_Cov.RData"))
GetMeanCovs <- function(dataset){
  Spec6_Shrink <- dataset
  Spec6_Shrink[,c("Season", "Species", "Order", "Family", "Genus", "Count", "MigCode")] <- NULL
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  Range <- function(x) {
    max(x) - min(x)
  }
  
  variables <- VariablesForMatchingByYear
  variables <- variables[!variables %in% c("SW", "Travel", "Slope")]
  staticvariables <- c("ISO3", "Country", "GeoRegion", "GeoSubRegion", "SW", "LOTW", "Travel", "Slope")
  SitesCovCast <- dcast(Spec6_Shrink, SiteCode ~ ., Mode, value.var="Anthrome")
  names(SitesCovCast) <- c("SiteCode", "Anthrome")
  SitesCovCast <- cbind(SitesCovCast, do.call(cbind, pblapply(variables, function(x){
    SitesCovCast <- dcast(Spec6_Shrink, SiteCode ~ ., mean, value.var=x, na.rm = TRUE)
    names(SitesCovCast) <- c("SiteCode", x)
    return(as.data.frame(SitesCovCast)[2])
  })))
  
  SitesCovStatic <- unique(Spec6_Shrink[,c("SiteCode", staticvariables)])
  SitesCovCast <- merge(SitesCovCast, SitesCovStatic, by="SiteCode")
  SitesCovCast$AnthRound <- as.factor(round_any(SitesCovCast$Anthrome, 10))
  
  return(SitesCovCast)
}
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpecSeason'), function(x) length(unique(Dataset[,c(x)])))
}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #grey, orange, light blue, green, yellow, dark blue, dark orange, pink

#### Extract sites within PAs, merge with all data ####
PP <- readOGR("/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Data/ProtectedPlanet/WDPA_June2017-shapefile-polygons.shp")
Countries <- readOGR("/Users/hannahwauchope/Documents/ArcGIS/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp")

##Make a dataframe of sites and years
SitesAndYears <- unique(Spec6_Cov[,c("SiteCode", "Longitude", "Latitude","Year")])
SiteByYear <- dcast(SitesAndYears, SiteCode + Longitude + Latitude ~ ., length, value.var="Year")
SiteByYearMin <- dcast(SitesAndYears, SiteCode ~ ., min, value.var="Year")
SiteByYearMax <- dcast(SitesAndYears, SiteCode ~ ., max, value.var="Year")
SitesMinMax <- merge(SiteByYearMin, SiteByYearMax, by = "SiteCode")
names(SitesMinMax) <- c("SiteCode", "MinYear", "MaxYear")
names(SiteByYear) <- c("SiteCode", "Longitude","Latitude","NumYears")
SiteByYear <- merge(SiteByYear, SitesMinMax, by="SiteCode")

sitepoints <- cbind(SiteByYear$Longitude, SiteByYear$Latitude)
sitepoints <- SpatialPointsDataFrame(sitepoints, SiteByYear, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

##Overlay with PA shapefile
PPOverSites <- over(sitepoints, PP) #Do we lose duplicate PA entries? Maybe? How to fix?
PPSites <- cbind(as.data.frame(SiteByYear), PPOverSites)
#PPSites <- na.omit(PPSites)

NOTPPSites <- unique(PPSites$SiteCode)
NOTPPSites <- Spec6[!Spec6$SiteCode %in% NOTPPSites,]
#write.csv(NOTPPSites, "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter3/NOT_PPSites.csv", row.names=FALSE)

##Clean
#Remove UNESCO biosphere reserves & proposed sites
PPSites <- subset(PPSites, DESIG!="UNESCO-MAB Biosphere Reserve")
PPSites <- subset(PPSites, DESIG_ENG!="UNESCO-MAB Biosphere Reserve")
PPSites <- subset(PPSites, STATUS!="Proposed")

#Round everything to an appropriate number
round_df <- function(x, digits) {
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
PPSites <- round_df(PPSites,4)

PPSites <- subset(PPSites, STATUS_YR!=0) #Remove sites without a status year

#Ok and this cleans so that for each site we just have one PA going on (they sometimes have multiple entries and stuff). 
#We take minimum year if status year varies, and if there are multiple entries of minimum year just the first occurence of a particular year, 
#but with al the desigs entered in and a binary value for multiple min year or not. 
DuplicatedSites <- PPSites[duplicated(PPSites$SiteCode),]
NotDuplicatedSites <- PPSites[!duplicated(PPSites$SiteCode),]
NotDuplicatedSites <- NotDuplicatedSites[!NotDuplicatedSites$SiteCode %in% DuplicatedSites$SiteCode,]

DuplicatedSitesData <- PPSites[PPSites$SiteCode %in% DuplicatedSites$SiteCode,]
DuplicatedSitesData$STATUS_YR <- as.numeric(as.character(DuplicatedSitesData$STATUS_YR))
DuplicatedSitesData$SiteStatYr <- paste0(DuplicatedSitesData$SiteCode, ".", DuplicatedSitesData$STATUS_YR)
DuplicatedSitesCast <- dcast(DuplicatedSitesData, SiteCode~., min, value.var="STATUS_YR") #In cases where duplicated PAs have different desingation years, take the earliest year
names(DuplicatedSitesCast) <- c("SiteCode", "STATUS_YR")
DuplicatedSitesCast$SiteStatYr <- paste0(DuplicatedSitesCast$SiteCode, ".", DuplicatedSitesCast$STATUS_YR)

DuplicatedSitesYear <- DuplicatedSitesData[DuplicatedSitesData$SiteStatYr %in% DuplicatedSitesCast$SiteStatYr,] #So this is now a dataframe of the duplicates, but cut down just the minimum year for each duplicate

#But sometimes we have multiple PAs designated in the same year (or possibly the same just with different statuses)
DuplicatedYear <- subset(dcast(DuplicatedSitesYear, SiteStatYr~., length, value.var="SiteCode"), .>1)
DuplicatedYear$SiteCode <- str_split_fixed(DuplicatedYear$SiteStatYr, "[.]",2)[,1]
DuplicatedYear <- DuplicatedSitesYear[DuplicatedSitesYear$SiteCode %in% DuplicatedYear$SiteCode,] #Get all entries for those sites

DuplicatedYearCast <- aggregate(DESIG~SiteCode, DuplicatedYear, paste) #Gather together all the DESIGs into one column, separated by commas
DuplicatedYearCast$Duplicates <- 2
DuplicatedYearCast2 <- aggregate(DESIG_ENG~SiteCode, DuplicatedYear, paste) #Ditto DESIG_ENG
DuplicatedYearCast <- merge(DuplicatedYearCast, DuplicatedYearCast2, by="SiteCode")
DuplicatedYearCast$DESIG <- as.character(DuplicatedYearCast$DESIG)
DuplicatedYearCast$DESIG_ENG <- as.character(DuplicatedYearCast$DESIG_ENG)

DuplicatedYear <- DuplicatedSitesYear[duplicated(DuplicatedSitesYear$SiteCode),] #Get the duplicate sites JUST TAKES FIRST ENTRY FOR EACH SO WE HAVE ALL THE DATA FOR DESIG ETC BUT NOT FOR EVERYTHING ELSE
DuplicatedYear <- DuplicatedYear[!duplicated(DuplicatedYear$SiteCode),] #Get the duplicate sites JUST TAKES FIRST ENTRY FOR EACH SO WE HAVE ALL THE DATA FOR DESIG ETC BUT NOT FOR EVERYTHING ELSE

DuplicatedYear[,c("DESIG", "DESIG_ENG")] <- NULL
DuplicatedYear <- merge(DuplicatedYear, DuplicatedYearCast, by="SiteCode")

NOTDuplicatedYear <- DuplicatedSitesYear[!DuplicatedSitesYear$SiteCode %in% DuplicatedYear$SiteCode,] 
NOTDuplicatedYear$Duplicates <- 1
NOTDuplicatedYear <- data.frame(lapply(NOTDuplicatedYear, function(x){
  if(is.factor(x)==TRUE){
    column <- as.character(x)
    return(column)
  } else {
    return(x)
  }
}), stringsAsFactors=FALSE)

NOTDuplicatedYear <- rbind(NOTDuplicatedYear, DuplicatedYear)

NotDuplicatedSites <- data.frame(lapply(NotDuplicatedSites, function(x){
  if(is.factor(x)==TRUE){
    column <- as.character(x)
    return(column)
  } else {
    return(x)
  }
}), stringsAsFactors=FALSE)
NotDuplicatedSites$Duplicates <- 0

All_Sites <- rbind(NotDuplicatedSites, NOTDuplicatedYear) #And combine non duplicate with cleaned duplicates

# 0 = no duplicates, 1 = duplicates but there was one min year, 2 = duplicates by there was more than one min year

## Add to bird data ##
#Add to bird data
PA_Counts <- merge(PPSites, Spec6, by=c("SiteCode", "Longitude", "Latitude"))
#PA_Counts <- merge(PA_Bird_Counts, specstats, by.x="Species", by.y="BLsciname")
PA_Counts[,c("MinYear", "MaxYear")] <- NULL
PA_Counts$SiteSpec <- paste0(PA_Counts$Species, ".", PA_Counts$SiteCode)

PA_Counts_MinYearSS <- dcast(PA_Counts, SiteSpec~., min, value.var="Year")
names(PA_Counts_MinYearSS) <- c("SiteSpec", "MinYearSS")
PA_Counts_MaxYearSS <- dcast(PA_Counts, SiteSpec~., max, value.var="Year")
names(PA_Counts_MaxYearSS) <- c("SiteSpec", "MaxYearSS")
PA_Counts <- merge(PA_Counts, PA_Counts_MinYearSS, by="SiteSpec")
PA_Counts <- merge(PA_Counts, PA_Counts_MaxYearSS, by="SiteSpec")
PA_Counts$FallsWithin <- ifelse(PA_Counts$STATUS_YR > PA_Counts$MinYearSS & PA_Counts$STATUS_YR < PA_Counts$MaxYearSS, "TRUE", "FALSE")

Birds_fallsw <- subset(PA_Counts, FallsWithin=="TRUE")

#Get actual measured years
Birds_fallsw$BeforeAfterStat <- ifelse(Birds_fallsw$Year<=Birds_fallsw$STATUS_YR, "BeforeStat", "AfterStat")
Birds_fallsw_countedyears <- unique(Birds_fallsw[,c("SiteSpec", "STATUS_YR", "Longitude", "Latitude", "MinYearSS", "MaxYearSS", "Year", "BeforeAfterStat")])
Birds_fallsw_countedyearscast <- dcast(Birds_fallsw_countedyears, SiteSpec~BeforeAfterStat, length, value.var="Year")
Birds_fallsw <- merge(Birds_fallsw, Birds_fallsw_countedyearscast, by="SiteSpec")

write.csv(Birds_fallsw, "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter3/PPSites2.csv", row.names=FALSE)


#### Prepare Site Matching Data ####
PPSitesOriginal <- fread(paste0(ResultsFP, "PP_Sites2.csv"))
NOTPPSitesOriginal <- fread(paste0(ResultsFP, "NOT_PP_Sites.csv"))

PPSitesOriginal$STATUS_YR <- as.numeric(PPSitesOriginal$STATUS_YR)
names(PPSitesOriginal)[names(PPSitesOriginal)=="NAME"] <- "PAName"
PPSitesOriginal$BeforeStatYears <- (PPSitesOriginal$STATUS_YR-PPSitesOriginal$MinYearSS)+1
PPSitesOriginal$AfterStatYears <- (PPSitesOriginal$MaxYearSS-PPSitesOriginal$STATUS_YR)
PPSitesOriginal$SiteSpecSeason <- paste0(PPSitesOriginal$SiteSpec, ".", PPSitesOriginal$Season)
PPSitesOriginal$Treatment <- 1

NOTPPSitesOriginal$SiteSpecSeason <- paste0(NOTPPSitesOriginal$Species, ".", NOTPPSitesOriginal$SiteCode, ".", NOTPPSitesOriginal$Season)
NOTPPSites_MinYearSS <- dcast(NOTPPSitesOriginal, SiteSpecSeason~., min, value.var="Year")
names(NOTPPSites_MinYearSS) <- c("SiteSpecSeason", "MinYearSS")
NOTPPSites_MaxYearSS <- dcast(NOTPPSitesOriginal, SiteSpecSeason~., max, value.var="Year")
names(NOTPPSites_MaxYearSS) <- c("SiteSpecSeason", "MaxYearSS")
NOTPPSitesOriginal <- merge(NOTPPSitesOriginal, NOTPPSites_MinYearSS, by="SiteSpecSeason")
NOTPPSitesOriginal <- merge(NOTPPSitesOriginal, NOTPPSites_MaxYearSS, by="SiteSpecSeason")
NOTPPSitesOriginal$Treatment <- 0

AllSites <- as.data.frame(rbind(PPSitesOriginal, NOTPPSitesOriginal, fill=TRUE))
TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#Just dec to feb
AllSites <- subset(AllSites, Season=="DectoFeb")
TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#Remove all zero counts
RemoveZeros <- dcast(AllSites, SiteSpecSeason~., sum, value.var="Count")
RemoveZeros <- subset(RemoveZeros, .==0)
AllSites <- AllSites[!AllSites$SiteSpecSeason %in% RemoveZeros$SiteSpecSeason,]

TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#CleanCBC
AllSitesCBC <- subset(AllSites, Dataset=="CBC")
AllSitesCBC$SpecPop <- paste(AllSitesCBC$Species, AllSitesCBC$Population, sep="_") #Combine species name and population into one variable

CountvsHoursLogHours <- rbindlist(pblapply(unique(AllSitesCBC$SpecPop), function(i){ #Cycle through each SpecPop
  #print(i)
  birdie <- subset(AllSitesCBC, SpecPop==i) #Subset dataframe to specpop
  birdiecast <- dcast(birdie, SiteSpecSeason~., sum, value.var="Count")
  birdie <- birdie[,c("SiteSpecSeason", "Count", "Hours")]
  glmmy <- tryCatch(glm.nb(Count~log(Hours), link=log, data=birdie), error=function(e){
    print(e)
    return(NULL)
    }) #Run a GLM of counts vs. log of hours
  if(is.null(glmmy)){return(NULL)}
  significant <- as.data.frame(ifelse(summary(glmmy)$coeff[-1,4]< 0.05, coef(glmmy, silent=TRUE)[[2]], NA)) #Return the coefficient IF the p value is less than 0.05, otherwise NA
  significant$SpecPop <- i
  names(significant) <- c("Slope", "SpecPop")
  return(significant)
}))

CountvsHoursLogHours <- CountvsHoursLogHours[!is.na(CountvsHoursLogHours$Slope),] #Remove the NAs (i.e. insignificant relations)
CountvsHoursLogHours <- subset(CountvsHoursLogHours, Slope>0) #Remove the negative slopes

AllSitesCBCCleaned <- unique(AllSitesCBC[!AllSitesCBC$SpecPop %in% CountvsHoursLogHours$SpecPop,][,c("SiteSpecSeason")])
AllSites <- AllSites[!AllSites$SiteSpecSeason %in% AllSitesCBCCleaned,]
TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#Remove sites with < 10 years
Buffer <- 5
NYearBuffer <- 3
AllSitesUnique <- unique(AllSites[,c("SiteCode", "Year")])
AllSitesUniqueMax <- dcast(AllSitesUnique, SiteCode~., max, value.var="Year")
AllSitesUniqueMin <- dcast(AllSitesUnique, SiteCode~., min, value.var="Year")
AllSitesUniqueMax$Diff <- (AllSitesUniqueMax$.-AllSitesUniqueMin$.)+1
AllSitesUniqueMax <- subset(AllSitesUniqueMax, Diff>=(Buffer*2))
AllSites <- AllSites[AllSites$SiteCode %in% AllSitesUniqueMax$SiteCode,]

AllSitesTREAT <- subset(AllSites[AllSites$Treatment==1,], AfterStatYears>=Buffer & BeforeStatYears>=Buffer & BeforeStat>=NYearBuffer & AfterStat>=NYearBuffer)
AllSites <- rbind(subset(AllSites, Treatment==0), AllSitesTREAT)

TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#Standardise so Not PP sites only contain species present in PP sites
PPSites <- subset(AllSites, Treatment==1)
NOTPPSites <- subset(AllSites, Treatment==0)
NOTPPSites <- NOTPPSites[NOTPPSites$Species %in% PPSites$Species,]
AllSites <- rbind(PPSites, NOTPPSites)
TheNumbers(subset(AllSites, Treatment==0))

#
AllSitesData <- unique(AllSites[,c("SiteSpecSeason", "Treatment", "MinYearSS", "MaxYearSS", "AfterStat", "BeforeStat", "BeforeStatYears", "AfterStatYears", "WDPAID", "PAName", "DESIG_TYPE", "IUCN_CAT", "GIS_AREA", "NO_TAKE", "STATUS", "STATUS_YR", "METADATAID")])

Spec6_AllSites <- merge(Spec6_Cov, AllSitesData, by="SiteSpecSeason")
Spec6_AllSites$AnthRound <- as.factor(round_any(Spec6_AllSites$Anthrome, 10))
Spec6_AllSites$SiteSpecSeasonYear <- paste0(Spec6_AllSites$SiteSpecSeason, ".", Spec6_AllSites$Year)
GovMean <- dcast(Spec6_AllSites, SiteCode ~ ., mean, value.var="Gov", na.rm = TRUE) #Get a mean governance value for each country
  names(GovMean) <- c("SiteCode", "GovMean")
  Spec6_AllSites <- merge(Spec6_AllSites, GovMean, by="SiteCode", all=T)

##Find collinear variables
ContinuousVariables <- Spec6_AllSites[,c("SiteSpecSeasonYear", "Treatment", "PrecipAnnual", "PrecipSeason", "MeanTempAnnual", "MinTempAnnual", "MaxTempAnnual", "MeanTempSeason", "MinTempSeason", "MaxTempSeason",
                                    "Nitr", "Phos", "grazing", "ir_norice", "ir_rice", "pasture", "rangeland", "rf_norice", "rf_rice", "popd", "uopp", "popc", "GovMean", 
                                    "SW", "Travel", "Slope")] #Get just continuous variables

ContinuousVariables <- ContinuousVariables[complete.cases(ContinuousVariables),]
Spec6_AllSites <- Spec6_AllSites[Spec6_AllSites$SiteSpecSeasonYear %in% ContinuousVariables$SiteSpecSeasonYear,]

#A few species got knocked out in that cull, so remove them from the control sites
SpeciesCheck <- dcast(Spec6_AllSites, Species~Treatment, length, value.var="SiteCode")
names(SpeciesCheck) <- c("Species", "Cont", "Treat")
SpeciesToRemove <- subset(SpeciesCheck, Treat==0)
Spec6_AllSites <- Spec6_AllSites[!Spec6_AllSites$Species %in% SpeciesToRemove$Species,]

TheNumbers(subset(Spec6_AllSites, Treatment==0)) #          8925            262         182180 (SiteCode, Species, SiteSpecSeason)
TheNumbers(subset(Spec6_AllSites, Treatment==1)) #            547            262           9693  

ContinuousVariables <- Spec6_AllSites[,c("SiteSpecSeasonYear", "Treatment", "PrecipAnnual", "PrecipSeason", "MeanTempAnnual", "MinTempAnnual", "MaxTempAnnual", "MeanTempSeason", "MinTempSeason", "MaxTempSeason",
                                         "Nitr", "Phos", "grazing", "ir_norice", "ir_rice", "pasture", "rangeland", "rf_norice", "rf_rice", "popd", "uopp", "popc", "GovMean", 
                                         "SW", "Travel", "Slope")] #Get just continuous variables

#Find collinearity
CorrPred <- ContinuousVariables
CorrPred[,c("SiteSpecSeasonYear", "Treatment")] <- NULL
SumCheck <- colSums(CorrPred)
SumCheck <- SumCheck[SumCheck>0]
SumCheck
CorrPred2 <- CorrPred[,which(names(CorrPred) %in% names(SumCheck))]
vif(CorrPred2)
CorrPred2[,c("grazing", "MinTempSeason")] <- NULL
vif(CorrPred2)
CorrPred2[,c("popd", "MeanTempSeason")] <- NULL
vif(CorrPred2)
CorrPred2[,c("MeanTempAnnual")] <- NULL
vif(CorrPred2)
CorrPred2[,c("MaxTempSeason")] <- NULL
vif(CorrPred2)

CorrelationPearson <- cor(CorrPred2, method="pearson")
apply(CorrelationPearson, 2, function(x) x[x>0.70])
CorrPred2$Phos <- NULL
CorrelationPearson <- cor(CorrPred2, method="pearson")
apply(CorrelationPearson, 2, function(x) x[x>0.70])
VariablesForMatchingByYear <- colnames(CorrelationPearson)

#Finish prepping
ContinuousVariablesToRemove <- ContinuousVariables[, !colnames(ContinuousVariables) %in% c("SiteSpecSeasonYear", "Treatment", VariablesForMatchingByYear)]
Spec6_AllSites <- Spec6_AllSites[,!colnames(Spec6_AllSites) %in% colnames(ContinuousVariablesToRemove)]

write.csv(Spec6_AllSites, paste0(ResultsFP, "Spec6_AllSites.csv"), row.names = FALSE)
save(VariablesForMatchingByYear, file=paste0(ResultsFP, "VariablesForMatchingByYear.RData"))

#### Conduct Matching ####
Spec6_AllSites <- as.data.frame(fread(paste0(ResultsFP, "Spec6_AllSites.csv")))
PPSites <- subset(Spec6_AllSites, Treatment==1)
NOTPPSites <- subset(Spec6_AllSites, Treatment==0)
load(file=paste0(ResultsFP, "VariablesForMatchingByYear.RData"))
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpecSeason'), function(x) length(unique(Dataset[,c(x)])))
}

### Calculate direction of trends
PPSites$BA <- ifelse(PPSites$Year>PPSites$STATUS_YR, 1, 0)
PPSitesBefore <- subset(PPSites, BA==0)
BeforeSlopes <- pbmclapply(unique(PPSitesBefore$SiteSpecSeason), function(SS){
  BeforeDat <- subset(PPSitesBefore, SiteSpecSeason==SS)[,c("Count", "Year", "Hours", "Dataset")]
  if(length(unique(BeforeDat$Count))==1 & unique(BeforeDat$Count)==0){
    return(as.data.frame(cbind(NA, NA, SS)))
  }
  if(unique(BeforeDat$Dataset=="CBC")){
    glmmy <- tryCatch(glm.nb(Count~Year + offset(log(x$Hours)), link=log, data=BeforeDat), error=function(e){NA})
  } else {
    glmmy <- tryCatch(glm.nb(Count~Year, link=log, data=BeforeDat), error=function(e){NA})
  }
  slope <- tryCatch(coef(glmmy, silent=TRUE)[[2]], error=function(e){NA})  
  significant <- tryCatch(summary(glmmy)$coeff[-1,4], error=function(e){NA})
  return(as.data.frame(cbind(slope, significant, SS)))
}, mc.cores=4)

BeforeSlopes2 <- rbindlist(BeforeSlopes[unlist(sapply(BeforeSlopes, function(x) ifelse(ncol(x)==2, FALSE, TRUE)))])
BeforeSlopes2$significant <- as.character(BeforeSlopes2$significant)
BeforeSlopes2$slope <- as.character(BeforeSlopes2$slope)
NumYears <- dcast(PPSitesBefore, SiteSpecSeason~., length, value.var="Year")
BeforeSlopes2 <- merge(BeforeSlopes2, NumYears, by.x="SS", by.y="SiteSpecSeason")
BeforeSlopes2$Class <- ifelse(is.na(BeforeSlopes2$significant), 0, #If the significance is NA make it zero
                              ifelse(BeforeSlopes2$.>6, ifelse(BeforeSlopes2$slope>0,1,-1), #If the significance is not NA and there's more than 6 years, make it the slope
                                     ifelse(BeforeSlopes2$significant>0.05, 0, #If the significance is not NA, it's 6 or less years and the significance is >0.05 make it zero
                                            ifelse(BeforeSlopes2$slope>0, 1, -1)))) #Otherwise make it the slope

BeforeSlopes2 <- BeforeSlopes2[,c("SS", "Class")]

PPSites <- merge(PPSites, BeforeSlopes2, by.x="SiteSpecSeason", by.y="SS")
TheNumbers(PPSites)
PPSites <- subset(PPSites, STATUS_YR>1970) #Not enough data to get a good mahalanobis estimate. ONLY ONE SITE
TheNumbers(PPSites)

#Create MD matrices for each designation year
AllYears <- unique(PPSites$STATUS_YR)
MDYearList <- pbmclapply(AllYears, function(YEAR){
  PPSitesCovs <- GetMeanCovs(subset(PPSites, Year<= YEAR & STATUS_YR==YEAR))
  NOTPPSitesCovs <- GetMeanCovs(subset(NOTPPSites, Year<= YEAR))
  MD <- as.data.frame(mahalanobis.dist(data.x= NOTPPSitesCovs[,colnames(NOTPPSitesCovs) %in% VariablesForMatchingByYear], data.y=PPSitesCovs[,colnames(PPSitesCovs) %in% VariablesForMatchingByYear]))
  rownames(MD) <- NOTPPSitesCovs$SiteCode
  colnames(MD) <- PPSitesCovs$SiteCode
  MDScale <- scale(MD,center=rep(0.5, ncol(MD)), scale=rep(10, ncol(MD)))
  return(MDScale)
}, mc.cores=ncores)

names(MDYearList) <- AllYears
save(MDYearList, file=paste0(ResultsFP, "MDYearList.RData"))
save(PPSites, file=paste0(ResultsFP, "PPSites.RData"))
save(NOTPPSites, file=paste0(ResultsFP, "NOTPPSites.RData"))

##Begin Matching
load(file=paste0(ResultsFP, "VariablesForMatchingByYear.RData"))
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpecSeason'), function(x) length(unique(Dataset[,c(x)])))
}
load(file=paste0(ResultsFP, "MDYearList.RData"))
load(file=paste0(ResultsFP, "PPSites.RData"))
load(file=paste0(ResultsFP, "NOTPPSites.RData"))

#Create empty MD matrix of all sites
PPSiteNames <- unique(PPSites$SiteCode)
NOTPPSiteNames <- unique(NOTPPSites$SiteCode)

MD <- matrix(nrow = length(NOTPPSiteNames), ncol = length(PPSiteNames))
rownames(MD) <- NOTPPSiteNames
colnames(MD) <- PPSiteNames

Spec <- "Actitis hypoleucos"

Matching <- rbindlist(pblapply(unique(PPSites$Species), function(Spec){
  if(file.exists(paste0(ResultsFP, "MatchFinal/Matched", Spec,".csv"))){
    return(NULL)
  }
  
  print(Spec)
  write.csv(NULL, paste0(ResultsFP, "MatchFinal/Matched", Spec,".csv"))
  
  #Reduce treatment to only sites with buffer years before and after
  PPSitesSpec <- subset(PPSites, Species==Spec & Treatment==1 & AfterStatYears>=5 & BeforeStatYears>=5)
  #PPSiteCovs <- unique(PPSitesSpec[,c("SiteCode", "Treatment", "Year", VariablesForMatchingByYear)])
  
  if(nrow(PPSitesSpec)==0){return(NULL)}
  NOTPPSitesSpec <- subset(NOTPPSites, Species==Spec & Treatment==0)
  SpecAllData <- data.table::rbindlist(list(PPSitesSpec, NOTPPSitesSpec), use.names=TRUE, fill=TRUE)
  SpecAllData$Match <- "Not Matched"
  if(length(unique(PPSitesSpec$SiteCode))<2){
    return(NULL)
  }
  MatchMatrix <- MD[,colnames(MD) %in% PPSitesSpec$SiteCode]
  MatchMatrix <- MatchMatrix[rownames(MatchMatrix) %in% NOTPPSitesSpec$SiteCode,]
  MatchMatrix <- as.data.frame(MatchMatrix)
  
  #Now filter out exact matched sites (by giving the difference in propensity scores/MD a value of 1000)
  for(i in colnames(MatchMatrix)){
    print(i)
    StatYr <- unique(subset(PPSitesSpec, SiteCode==i)$STATUS_YR)
    ProtCovs <- GetMeanCovs(subset(PPSitesSpec, SiteCode==i & Year<= StatYr))
    if(nrow(subset(NOTPPSitesSpec, Year<= StatYr))==0){next}
    UnProtCovs <- GetMeanCovs(subset(NOTPPSitesSpec, Year<= StatYr))
    
    Anth <- as.character(ProtCovs$AnthRound)
    Region <- as.character(ProtCovs$GeoRegion)
    SlopeClass <- unique(subset(PPSitesSpec, SiteCode==i)$Class)
    NotProtSub <- subset(UnProtCovs, AnthRound==Anth & GeoRegion==Region)
    NotProtSub <- NOTPPSitesSpec[NOTPPSitesSpec$SiteCode %in% NotProtSub$SiteCode,]
    NotProtSub$BeforeStatYears <- (StatYr-NotProtSub$MinYearSS)+1 
    NotProtSub$AfterStatYears <- (NotProtSub$MaxYearSS-StatYr)
    
    #This block gets us to cases with enough data before and after (at least 3 years sampled, at least 5 years covered), now we need to run GLMs
    NotProtSub <- subset(NotProtSub, BeforeStatYears>=5 & AfterStatYears>=5)
    if(nrow(NotProtSub)==0){
      MatchMatrix[,i] <- 10
      next
    }
    NotProtSub$BA <- ifelse(NotProtSub$Year<=StatYr, 0, 1)
    NotProtSubYears <- dcast(NotProtSub, SiteCode~BA, length, value.var="Year")
    names(NotProtSubYears) <- c("SiteCode", "Before", "After")
    NotProtSubYears <- subset(NotProtSubYears, Before>=3 & After>=3)
    NotProtSub <- NotProtSub[NotProtSub$SiteCode %in% NotProtSubYears$SiteCode,]
    
    NotProtSub <- subset(NotProtSub, BA==0)
    if(nrow(NotProtSub)==0){
      MatchMatrix[,i] <- 10
      next
    }
    BeforeSlopes <- pbmclapply(unique(NotProtSub$SiteSpecSeason), function(SS){
      BeforeDat <- subset(NotProtSub, SiteSpecSeason==SS)[,c("Count", "Year", "Hours", "Dataset")]
      if(length(unique(BeforeDat$Count))==1 & unique(BeforeDat$Count)==0){
        return(as.data.frame(cbind(NA, NA, SS)))
      }
      if(unique(BeforeDat$Dataset=="CBC")){
        glmmy <- tryCatch(glm.nb(Count~Year + offset(log(x$Hours)), link=log, data=BeforeDat), error=function(e){NA})
      } else {
        glmmy <- tryCatch(glm.nb(Count~Year, link=log, data=BeforeDat), error=function(e){NA})
      }
      slope <- tryCatch(coef(glmmy, silent=TRUE)[[2]], error=function(e){NA})  
      significant <- tryCatch(summary(glmmy)$coeff[-1,4], error=function(e){NA})
      return(as.data.frame(cbind(slope, significant, SS)))
    }, mc.cores=4)
    BeforeSlopes2 <- as.data.frame(rbindlist(BeforeSlopes[unlist(sapply(BeforeSlopes, function(x) ifelse(ncol(x)==2, FALSE, TRUE)))]))
    names(BeforeSlopes2) <- c("slope", "significant", "SS")
    BeforeSlopes2$significant <- as.character(BeforeSlopes2$significant)
    BeforeSlopes2$slope <- as.character(BeforeSlopes2$slope)
    NumYears <- dcast(NotProtSub, SiteSpecSeason~., length, value.var="Year")
    BeforeSlopes2 <- merge(BeforeSlopes2, NumYears, by.x="SS", by.y="SiteSpecSeason")
    BeforeSlopes2$Class <- ifelse(is.na(BeforeSlopes2$significant), 0, #If the significance is NA make it zero
                                  ifelse(BeforeSlopes2$.>6, ifelse(BeforeSlopes2$slope>0,1,-1), #If the significance is not NA and there's more than 6 years, make it the slope
                                         ifelse(BeforeSlopes2$significant>0.05, 0, #If the significance is not NA, it's 6 or less years and the significance is >0.05 make it zero
                                                ifelse(BeforeSlopes2$slope>0, 1, -1)))) #Otherwise make it the slope
    BeforeSlopes3 <- subset(BeforeSlopes2, Class==SlopeClass)
    print(c(SlopeClass, unique(BeforeSlopes3$Class)))
    if(nrow(BeforeSlopes3)!=0){
      FilteredSites <- str_split_fixed(BeforeSlopes3$SS, "[.]",3)[,2]
    } else {
      FilteredSites <- NULL
    }
    UnProtCovs <- UnProtCovs[UnProtCovs$SiteCode %in% FilteredSites,]
    
    ### Ok and now get the distances and add to matrix
    DistanceValues <- MDYearList[[paste0(StatYr)]]
    if(ncol(DistanceValues)==1){
      DistanceValues <- DistanceValues[rownames(DistanceValues) %in% UnProtCovs$SiteCode,]
      if(length(DistanceValues)==1){
        k <- paste0(UnProtCovs$SiteCode)
        MatchMatrix[k, i] <- DistanceValues
        MatchMatrix[!rownames(MatchMatrix) %in% UnProtCovs$SiteCode,i] <- 10
      } else {
        MatchMatrix[c(names(DistanceValues)), i] <- DistanceValues
        MatchMatrix[!rownames(MatchMatrix) %in% UnProtCovs$SiteCode,i] <- 10
      }
    } else {
      DistanceValues <- DistanceValues[rownames(DistanceValues) %in% UnProtCovs$SiteCode,]
      if(is.null(nrow(DistanceValues))){
        k <- paste0(UnProtCovs$SiteCode)
        MatchMatrix[k, i] <- DistanceValues[[i]]
        MatchMatrix[!rownames(MatchMatrix) %in% UnProtCovs$SiteCode,i] <- 10
      } else {
        MatchMatrix[c(names(DistanceValues[,i])), i] <- DistanceValues[,i]
        MatchMatrix[!rownames(MatchMatrix) %in% UnProtCovs$SiteCode,i] <- 10
      }
    }
  }
  
  MatchMatrix$Idiot <- "Hah, Idiot" #Because r is stupid and if there's only one column it won't return the rowname of the min value
  #Assess competition
  #Turn it into a loop
  Greedy <- pbmclapply(1:1000, function(iteration){
    MatchMatrix3 <- MatchMatrix
    MatchedSites <- data.frame("Treatment" = colnames(MatchMatrix3), "Control" = NA, "MD" = NA)
    Sample <- sample.int(ncol(MatchMatrix3)-1)
    for(x in Sample){
      paste(x)
      if(min(MatchMatrix3[,x])==10){
        next()
      } else {
        Treatment <- rownames(MatchMatrix3[MatchMatrix3[,x] == min(MatchMatrix3[,x]),])
        if(length(Treatment)>1){Treatment <- Treatment[1]}
        MatchedSites[x,2] <- Treatment
        MatchedSites[x,3] <- MatchMatrix3[MatchMatrix3[,x] == min(MatchMatrix3[,x]),x][1]
        MatchMatrix3[Treatment, ] <- 10
      }
    }
    return(MatchedSites)
  }, mc.cores=ncores)
  
  GlobalDistance <- sapply(Greedy, function(x){
    Greedyx <- x[complete.cases(x),]
    return(sum(Greedyx$MD))
  })
  Matched <- Greedy[[match(min(GlobalDistance),GlobalDistance)]]
  Matched <- subset(Matched, Treatment!="Idiot")
  Matched <- Matched[!is.na(Matched$Control),]
  if(nrow(Matched)==0){
    return(NULL)
  }
  Matched$MatchID <- 1:nrow(Matched)
  Matched$Treatment <- as.character(Matched$Treatment)
  Matched$Control <- as.character(Matched$Control)
  Matched <- melt(Matched, id.vars=c("MatchID", "MD"))
  names(Matched) <- c("MatchID", "MD", "Treatment", "SiteCode")
  Matched$Treatment <- NULL
  Matched <- merge(Matched, SpecAllData, by="SiteCode", all.x=TRUE)
  Matched$Match <- "Matched"
  
  SpecAllData$MatchID <- NA
  SpecAllData$MD <- NA
  SpecAllData <- as.data.frame(apply(SpecAllData, 2, as.character))
  Matched <- as.data.frame(apply(Matched, 2, as.character))
  BindTogether <- list(Matched, SpecAllData)
  SpecSCCAll <- do.call(rbind, lapply(BindTogether, function(x) x[match(names(BindTogether[[1]]), names(x))]))
  SpecSCCAll$Species <- Spec
  write.csv(SpecSCCAll, paste0(ResultsFP, "MatchFinal/Matched", Spec,".csv"), row.names = FALSE)
  return(SpecSCCAll)
}))

#### Assess matching ####
MatchingFinal <- as.data.frame(subset(rbindlist(lapply(list.files(path=paste0(ResultsFP, "MatchFinal/"), full.names=TRUE), fread)), Match=="Matched"))
MatchingFinal <- MatchingFinal[complete.cases(MatchingFinal$SiteCode),]
MatchingFinal$MatchID <- as.numeric(MatchingFinal$MatchID)
MatchingFinal$SpecMatch <- paste0(MatchingFinal$Species, ".", MatchingFinal$MatchID)
MatchingCovariates <- rbindlist(pbmclapply(unique(MatchingFinal$SpecMatch), function(SpecMatchID){
  try <- subset(MatchingFinal, SpecMatch==SpecMatchID)
  StatYr <- unique(try$STATUS_YR)[complete.cases(unique(try$STATUS_YR))]
  ProtCovs <- GetMeanCovs(subset(try, Treatment==1 & Year<= StatYr))
  ProtCovs$Treatment <- 1
  UnProtCovs <- GetMeanCovs(subset(try, Treatment==0 & Year<= StatYr))
  UnProtCovs$Treatment <- 0
  AllCovs <- rbind(ProtCovs, UnProtCovs)
  AllCovs$SpecMatch <- SpecMatchID
  return(AllCovs)
}, mc.cores=4)) #Get the covariate values for predesignation years (now that we know which PA each unprotected site is matched to)
MatchingCovariates$Species <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,1]
MatchingCovariates$MatchID <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,2]

Distances <- unique(MatchingFinal[,c("SpecMatch", "MD")])
MatchingCovariates <- merge(MatchingCovariates, Distances, by="SpecMatch", all=T)

AbsDist <- function(dataset){
  Distance <- rbindlist(lapply(unique(dataset$Species), function(Spec){
    ok <- as.data.frame(subset(dataset, Species==Spec))
    ok <- ok[complete.cases(ok$SiteCode),]
    ok2 <- ok[,colnames(ok) %in% c(VariablesForMatchingByYear, "Treatment")]
    ok2$MatchID <- ok$MatchID
    ok2 <- melt(ok2, id.var=c("MatchID", "Treatment"))
    hey <- dcast(ok2, MatchID + Treatment~variable, value.var="value")
    hey2 <- rbindlist(lapply(unique(hey$MatchID), function(ID){
      hey2 <- subset(hey, MatchID==ID)
      if(nrow(hey2)<2){
        return(NULL)
      }
      hey2 <- as.data.frame(t(as.data.frame(apply(hey2, 2, function(x) diff(as.numeric(x))))))
      hey2$MatchID <- ID
      row.names(hey2) <- NULL
      return(hey2)
    }))
    hey3 <- melt(hey2, id.vars="MatchID")
    hey4 <- dcast(hey3, variable~., mean, value.var="value")
    hey4$Species <- Spec
    names(hey4) <- c("Variable", "MeanDist", "Species")
    return(hey4)
  }))
  return(Distance)
}
StDiffMean <- function(dataset){
  #SpecSites <- subset(dataset, Species==Spec)
  SpecSites <- dataset[complete.cases(dataset$MD),]
  if(nrow(SpecSites)<=2){
    return(NULL)
  }
  SiteCast <- do.call(cbind,lapply(c(1,0), function(x){ #Get the mean values for each variable in each treatment
    TheVariables <- as.data.frame(subset(SpecSites, Treatment==x))
    TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear]
    TheVariables$Treatment <- NULL
    TheVariables <- apply(TheVariables, 2, as.numeric)
    return(as.data.frame(colMeans(TheVariables, na.rm = TRUE)))
  }))
  names(SiteCast) <- c("Treat", "Cont")
  SiteCast$MeanDiff <- abs(SiteCast$Treat - SiteCast$Cont)
  SiteCastVar <- do.call(cbind,lapply(c(1,0), function(x){ #Get the sd for each variable in each treatment
    TheVariables <- as.data.frame(subset(SpecSites, Treatment==x))
    TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear]
    TheVariables$Treatment <- NULL
    TheVariables <- apply(TheVariables, 2, as.numeric)
    return(as.data.frame(colVars(TheVariables, na.rm = TRUE)))
  }))
  names(SiteCastVar) <- c("Treat", "Cont")
  SiteCast$VarComp <- sqrt((SiteCastVar$Treat + SiteCastVar$Cont)/2)
  SiteCast$d <- SiteCast$MeanDiff/SiteCast$VarComp
  SiteCast[is.na(SiteCast)] <- 0
  SiteCast$Cov <- row.names(SiteCast)
  SiteCast <- SiteCast[,c(5:6)]
  #SiteCast$Species <- Spec
  #SiteCast <- dcast(SiteCast, Cov~Match, value.var="d")
  #return(mean(SiteCast$Match))
  return(SiteCast)
}
CompareMatching <- function(dataset, DiffThresh, MatchType){
  DeleteOldFiles <- c(list.files(paste0(ResultsFP, "Summaries/", MatchType), full.names=TRUE), list.files(paste0(ResultsFP, "MatchData/", MatchType), full.names=TRUE))
  sapply(DeleteOldFiles, unlink)
  Comp <- pbmclapply(unique(dataset$Species), function(i){
    print(i)
    MatchingSub <- subset(dataset, Species==i)
    MatchingSub <- MatchingSub[!is.na(MatchingSub$MD),]
    if(nrow(MatchingSub)<=4){
      return(NULL)
    }
    Iteration <- 1
    TheDistance <- StDiffMean(MatchingSub)
    TheDistance$Iteration <- Iteration
    Threshold <- subset(TheDistance, d<DiffThresh)
    if(nrow(Threshold)!=0){
      ThresholdCast <- dcast(Threshold, Iteration~Cov, value.var="d")
      ThresholdCast <- ThresholdCast[complete.cases(ThresholdCast),]
      if(ncol(ThresholdCast)!=length(VariablesForMatchingByYear)){
        ThresholdCast <- data.frame()
      }
    } else {
      ThresholdCast <- data.frame()
    }
    
    while(nrow(ThresholdCast)==0 & nrow(MatchingSub)>4){
      #print(Iteration)
      MAX <- max(MatchingSub$MD)
      MatchingSub <- subset(MatchingSub, MD!=MAX)
      TheDistance <- StDiffMean(MatchingSub)
      TheDistance$Iteration <- Iteration+1
      Iteration <- Iteration+1
      Threshold <- subset(TheDistance, d<DiffThresh)
      if(nrow(Threshold)==0){
        ThresholdCast <- data.frame()
        next
      }
      ThresholdCast <- dcast(Threshold, Iteration~Cov, value.var="d")
      ThresholdCast <- ThresholdCast[complete.cases(ThresholdCast),]
      if(ncol(ThresholdCast)!=length(VariablesForMatchingByYear)){
        ThresholdCast <- data.frame()
      }
    }
    if(nrow(ThresholdCast)>0){
      ThresholdCast <- as.data.frame(t(ThresholdCast))
      names(ThresholdCast) <- "StDiff"
      ThresholdCast$Variable <- row.names(ThresholdCast)
      AbsoluteDist <- as.data.frame(AbsDist(MatchingSub))
      Summary <- merge(ThresholdCast, AbsoluteDist, by="Variable")
      write.csv(Summary, paste0(ResultsFP, "Summaries/Final/Summary_", i, ".csv"), row.names=FALSE)
      write.csv(MatchingSub, paste0(ResultsFP, "MatchData/Final/MatchData_", i, ".csv"), row.names=FALSE)
      
      return(list(Summary, MatchingSub))
    } else {
      return(NULL)
    }
  }, mc.cores=ncores)
}
CompMatch <- CompareMatching(MatchingCovariates, 0.25, "Final") #Run through and iteratively remove matches until we have the matched set for each species

SummariesFinal <- rbindlist(lapply(list.files(path=paste0(ResultsFP, "Summaries/Final/"), full.names=TRUE), fread)) #Read in that summary material
names(SummariesFinal) <- c("Variable", "StDiff", "MeanDist", "Species")
MatchDataFinal <- rbindlist(lapply(list.files(path=paste0(ResultsFP, "MatchData/Final/"), full.names=TRUE), fread))
MatchingFinalCleaned <- MatchingFinal[MatchingFinal$SpecMatch %in% MatchDataFinal$SpecMatch,]

TheNumbers(subset(MatchingFinalCleaned, Treatment==1))
TheNumbers(subset(MatchingFinalCleaned, Treatment==0))
TheNumbers(PPSites)
TheNumbers(NOTPPSites)

#Extract species stats for comparison
TaxaData <- function(Dataset){
  AllSitesStats <- unique(Dataset[,c("Species", "SiteCode", "Treatment")])
  AllSitesStats <- merge(AllSitesStats, SpecStats, by="Species")
  AllSitesStatsFamSite <- unique(AllSitesStats[,c("Order", "Family", "SiteCode", "Treatment")])
  FamiliesSite <- dcast(AllSitesStatsFamSite, Order + Family~Treatment, length, value.var="SiteCode")
  AllSitesStatsFamSpec <- unique(AllSitesStats[,c("Order", "Family","Genus", "Species", "Treatment")])
  FamiliesSpec <- dcast(AllSitesStatsFamSpec, Order + Family~., length, value.var="Species")
  AllSitesStatsFamGen <- unique(AllSitesStats[,c("Order", "Family", "Genus", "Treatment")])
  FamiliesGen <- dcast(AllSitesStatsFamGen, Order + Family~., length, value.var="Genus")
  FamilyData <- merge(FamiliesSite, FamiliesGen, all=TRUE)
  names(FamilyData) <- c("Order", "Family", "UnprotectedSites", "ProtectedSites", "Genera")
  FamilyData <- merge(FamilyData, FamiliesSpec, all=TRUE)
  names(FamilyData) <- c("Order", "Family", "UnprotectedSites", "ProtectedSites", "Genera", "Species")
  NOrder <- length(unique(FamilyData$Order))
  NFamily <- length(unique(FamilyData$Family))
  NGenus <- length(unique(AllSitesStatsFamSpec$Genus))
  NSpecies <- length(unique(AllSitesStatsFamSpec$Species))
  return(list(FamilyData, NOrder,NFamily,NGenus, NSpecies))
}
AllSites2 <- rbind(PPSites[,c("SiteCode", "Species", "Treatment")], NOTPPSites[,c("SiteCode", "Species", "Treatment")])
AllTaxa <- TaxaData(AllSites2)
AllTaxa
MatchTaxa <- TaxaData(MatchDataFinal)
MatchTaxa

BothTaxa <- merge(AllTaxa[[1]], MatchTaxa[[1]], by=c("Order", "Family"), all=TRUE)
BothTaxa$UnprotectedSites <- ifelse(is.na(BothTaxa$UnprotectedSites.y), BothTaxa$UnprotectedSites.x, paste0(BothTaxa$UnprotectedSites.x, " (", BothTaxa$UnprotectedSites.y,")"))
BothTaxa$ProtectedSites <- ifelse(is.na(BothTaxa$ProtectedSites.y), BothTaxa$ProtectedSites.x, paste0(BothTaxa$ProtectedSites.x, " (", BothTaxa$ProtectedSites.y,")"))
BothTaxa$Genera <- ifelse(is.na(BothTaxa$Genera.y), BothTaxa$Genera.x, paste0(BothTaxa$Genera.x, " (", BothTaxa$Genera.y,")"))
BothTaxa$Species <- ifelse(is.na(BothTaxa$Species.y), BothTaxa$Species.x, paste0(BothTaxa$Species.x, " (", BothTaxa$Species.y,")"))
BothTaxa <- BothTaxa[,c("Order", "Family", "Genera", "Species","ProtectedSites", "UnprotectedSites")]
BothTaxa$Order <- as.character(BothTaxa$Order)
BothTaxa <- BothTaxa[order(BothTaxa$Order),]

AllOrdAllFam <- as.data.frame(matrix(nrow=1, ncol=6, c("All Orders", "All Families", paste0(AllTaxa[[4]], " (", MatchTaxa[[4]], ")"), paste0(AllTaxa[[5]], " (", MatchTaxa[[5]], ")"), 
                                                       paste0(TheNumbers(PPSites)[[1]], " (", TheNumbers(subset(MatchingFinalCleaned, Treatment==1))[[1]], ")"), 
                                                       paste0(TheNumbers(NOTPPSites)[[1]], " (", TheNumbers(subset(MatchingFinalCleaned, Treatment==0))[[1]], ")"))))
names(AllOrdAllFam) <- names(BothTaxa)
BothTaxa <- rbind(BothTaxa, AllOrdAllFam)
write.csv(BothTaxa, "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter3/Figures/TaxaCompTable.txt", row.names=FALSE)