############################################################################################################################################################################
### This code was written by Hannah Wauchope to analyse data for the paper "Robust study design shows variable impact of protected areas on waterbirds"
### This is script 4 of 4 in the workflow
### Last edited 27th August, 2021
### Please direct queries to hannah.wauchope@gmail.com
###
### This script takes the data, and matches, obtained from Script 3 and does the following for each of the 20 Latin Hypercube analyses plus the focal analysis):
### Assesses the quality of matches for BACI and CI populations, reduces the data to only adequate matches
### Calculates the impact of protected areas on populations according to BA, CI and BACI frameworks, and classifies these into the various impact categories defined in Figure 3 and Extended Data Figures 1 and 2
### For BACI, runs models to see how various predictors correlate with protected area impact (see Extended Data Table 2 for details of predictors, sources of data, and explanations of data transformations/cleaning)
###
### The script then produces the various figures and tables found throughout the paper and supplementary material, and summarises any statistics reported in the paper
###
### NOTES: Throughout these scripts I use the field "SiteSpec" to mean Population, i.e. a particular species at a particular site. This is the unit that analyses are performed on for BA 
### NOTES: In a similar vein, "SpecMatch" is used to refer to two paired SiteSpecs, a protected and unprotected SiteSpec. In this way, a SpecMatch is the unit of analysis for CI and BACI (i.e. counts of particular species at one protected site and its matched unprotected site)
### NOTES: In the paper, we refer to three types of change that a protected area can elicit on a population: immediate, trend and average. In earlier versions of papers, we called these terms level, slope and mean (respectively). In this code I mostly use the latter set of names.
### NOTES: For reasons known early to past me, I labelled the cumulative link models that correlate predictors with protected area impacts (i.e. for Fig 4 etc) as "Category Models". Unsure why, but they are referred to as such throughout this script.
############################################################################################################################################################################

#### Initialise ####
cluster <- FALSE
if(cluster==TRUE){
  library(data.table, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(rgdal, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(sp, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(pbapply, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(dplyr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(stringr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(tidyr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(plyr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(pbmcapply, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(MASS, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(usdm, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(rlist, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(resample, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(foreign, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(lme4, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(scales, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(ordinal, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(carData, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/") 
  library(car, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/") 
  library(DHARMa, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/") 
  library(lmtest, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/") 
  library(nlme, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/") 
  library(lme4, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/") 
  library(stringr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(purrr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(dplyr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  ResultsFP <- "/gpfs/ts0/projects/Research_Project-T115802/PAs/"
  DataFP <- "/gpfs/ts0/projects/Research_Project-T115802/Data/"
  ncores <- 14
} else {
  library(data.table)
  library(rgdal)
  library(sp)
  library(wdpar)
  library(pbapply)
  library(tidyverse)
  library(MASS)
  library(plyr)
  library(pbmcapply)
  library(usdm)
  library(StatMatch)
  library(rlist)
  library(resample)
  library(foreign)
  library(lme4)
  library(scales)
  library(rredlist)
  library(ordinal)
  library(ggbeeswarm)
  library(fields)
  library(car) 
  library(DHARMa) 
  library(lmtest) 
  library(nlme) 
  library(lme4) 
  library(glmmTMB) #not on cluster
  library(ggalt) #not on cluster
  library(ggalluvial) #not on cluster
  library(ggeffects) #not on cluster
  library(cowplot) #Not on cluster
  library(RColorBrewer) #Not on cluster
  library(taxize) #Not on cluster
  library(gridExtra) #Not on cluster
  
  ncores <- 4
  
  ###File paths
  DataFP <- "/Users/hannahwauchope/DropBox/Work/Projects/PhD/Data/"
  ResultsFP <- "/Users/hannahwauchope/DropBox/Work/Projects/PhD/Chapter2_4_PAs/Analysis/"
  FiguresFP <-"/Users/hannahwauchope/DropBox/Work/Projects/PhD/Chapter2_4_PAs/Figures/"
  
}
options(scipen=999)

#### Define functions ####
#Define map projections
MollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
WGSCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#Defint font size to be used in plots
plotfontsize <- 7

#Function to round a number to a specified degree of freedom
round_df <- function(x, digits) {
  x <- as.data.frame(x)
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[,numeric_columns], digits)
  x
}

#Function to calculate mode of a categorical variable
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#Function to assess models for temporal autocorrelataion (TAC) - doesn't work on models with random effects. Requires mode object (Mod), the data input into the model (ModDat), whether the model is a BACI model or no (defaults to no), and also the TotalYearsBuffer
TACCheck <- function(Mod, ModDat, BABACI="BA", TotalYearsBuffer){
  ModDat$Year2 <- ModDat$Year
  
  if(BABACI!="BA"){
    ModDat[ModDat$CI==1,]$Year2 <- ModDat[ModDat$CI==1,]$Year2+(TotalYearsBuffer*2) #This is just a workaround that doesn't affect results. In order for DHARMa to calculate Durbin Watson properly it needs a vector of the order of years in the dataset, so this creates them with no year duplicates, by making dummy values for CI equals 1. e.g. the years of the model might have initially been 2002, 2003, 2004 for both CI =0 and CI =1, now they would be 2002, 2003, 2004 (for CI=0), and 2012, 2013, 2014 (for CI=1). 
  }
  
  #DHARMa
  simulationOutput <- simulateResiduals(fittedModel = Mod, n=1000)
  DW_DHARMa <- tryCatch(testTemporalAutocorrelation(simulationOutput, time=ModDat$Year2, plot=FALSE, alternative = "greater")$p, error=function(e){1})
  
  #Car
  DW_car <- tryCatch(durbinWatsonTest(Mod, alternative = "positive")$p, error=function(e){1})
  
  #lmtest
  DW_lmtest <- tryCatch(dwtest(Mod, alternative = "greater", iterations=100)$p, error=function(e){1})
  
  AllDW <- c(DW_DHARMa, DW_car, DW_lmtest)
  AllDW[is.na(AllDW)] <- 1
  
  return(AllDW)
}

#This function is used in the matching process to calculate 1 value for each matching covariate from the span of values over a number of years
#For instance, for a particular site, the matching process may require an estimate of mean temperature and the mode land use type over a 10 year period etc etc. This calculates it.
#"dataset" is the given data, structured with a column for the Site Codes, a column for year, and columns for each covariate. "VariablesForMatchingByYear" is a list of the covariates that means/modes are required for.
#The function returns a dataset with a column for each site code and a column for the mean/mode of each required covariate.
GetMeanCovs <- function(dataset, VariablesForMatchingByYear){
  CountsSub <- as.data.table(dataset) 
  
  #Function for calculating mode
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  #Function for calculating range
  Range <- function(x) {
    max(x) - min(x)
  }
  
  variables <- VariablesForMatchingByYear
  variables <- variables[!variables %in% c("SWOccurrence", "Travel", "Slope", "GovMean")] #Remove these variables as they do not change year by year
  SitesCovCast <- dcast(CountsSub, SiteCode ~ ., Mode, value.var="Anthrome", drop=FALSE) #Calculate the mode for anthrome
  names(SitesCovCast) <- c("SiteCode", "Anthrome")
  
  SitesCovCast <- cbind(SitesCovCast, do.call(cbind, pblapply(variables, function(x){ #Calculate the mean for all required variables 
    SitesCovCast <- dcast(CountsSub, SiteCode ~ ., mean, value.var=x, drop=FALSE)
    names(SitesCovCast) <- c("SiteCode", x)
    return(as.data.frame(SitesCovCast)[2])
  })))
  SitesCovStatic <- CountsSub[match(unique(CountsSub$SiteCode), CountsSub$SiteCode),][,c("SiteCode", "ISO3", "Country", "GeoRegion", "GeoSubRegion", "SWOccurrence", "LOTWID", "Travel", "Slope", "GovMean")] #Get the values for all covariates that don't change through time
  SitesCovCast <- merge(SitesCovCast, SitesCovStatic, by="SiteCode") #Bring data together
  SitesCovCast$AnthRound <- as.factor(round_any(SitesCovCast$Anthrome, 10)) #Get the higher order value for anthrome (e.g. 21, 22, 23 and 24 are all anthrome values for various types of village, this rounds to 20 which just means "village")
  return(SitesCovCast)
}

#This function is used to calculate the slope of one (or many) population time series using negative binomial GLMs. The name is a slight misnomer as I no long only use it to calculate slopes in the "before" period, but also the "after" period.
#The function recieves "before data" which is the main datastructure, requiring at minimum a column for SiteCode, Species, SiteSpec, Year, Count, Hours (the effort term for CBC data), CI (1 for protected sites, 0 for unprotected sites)
#It returns a column for each population ("sitespec"), a column called "BeforeSlope" which is 1 for significantly positive trends, 0 for insignificant trends, -1 for significantly negative trends, and "AllZero" if all zeroes. 
#ZeroThresh is one of the values defined in the latin hypercube sampling, a value between 0.6 and 0.9. If the proportion of zeroes in the time series if above the threshold, it is classified as "All Zeroes" (see Wauchope et al., 2021, Trends in ecology and evolution for an explanation of the rationale behind this)
#Parallise is logical, if true the function is parallelised. 
#actualvalues is a logical, if true then in addition to the "BeforeSlope" column, columns are returned that give the actual slope estimates and significances from the model.
#PThresh is included for a supplementary analysis that relaxes the p<0.05 cut off for significance, to see if we can better detected PA impact with slight trends
#TotalYearsBuffer required for the nested function that checks temporal autocorrelation.
CalculateBeforeSlopes <- function(BeforeData, ZeroThresh, parallelise, actualvalues, PThresh, TotalYearsBuffer){
  if(parallelise==TRUE){
    ncoresslopes <- ncores/2
  } else {
    ncoresslopes <- 1
  }
  
  if(nrow(BeforeData[is.na(BeforeData$Hours),])!=0){
    BeforeData[is.na(BeforeData$Hours),]$Hours <- 1 #For all the IWC sites, make "Hours" 1 (which means for modelling it just disappears, which we want as for IWC data effort is standardised)
  }  
  
  AllCounts <- AllZeroes(BeforeData, ZeroThresh) #Use the function "AllZeroes" (defined below) to identify the "AllZero" time series
  NonZero <- subset(AllCounts, AllZero=="Counts")
  AllZero <- subset(AllCounts, AllZero=="AllZero")
  names(AllZero) <- c("SiteSpec", "CI", "BeforeSlope")
  if(nrow(NonZero)==0){
    return(AllCounts)
  }
  NonZero$SiteSpec <- paste0(NonZero$SiteSpec, "_", NonZero$CI) #Combine CI and SiteSpec to get an individual name for each timeseries (e.g. the before period for a population, cf the after period, as these are calculated separately)
  BeforeData$SiteSpec <- paste0(BeforeData$SiteSpec, "_", BeforeData$CI)
  SS <- "ZA00455.Spatula smithii_1" #(here for checking)
  #For all the populations that aren't "All Zero" use GLMs to calculate the slope of each population (in cases where the model fails, return NA)
  BeforeSlopes <- pbmclapply(unique(NonZero$SiteSpec), function(SS){
    #print(SS)
    BeforeDat <- subset(BeforeData, SiteSpec==SS)[,c("Count", "Year", "Hours", "Dataset")]
    BeforeDat <- BeforeDat[order(BeforeDat$Year),]
    BeforeDat$Year <- BeforeDat$Year - min(BeforeDat$Year) + 1 #Scale years in case of need to use random function
    if(length(unique(BeforeDat$Count))==1){if(unique(BeforeDat$Count)==0){ #AllZero function misses time periods that are 100% zeroes, isolate any of these
      return(as.data.frame(cbind(slope="AllZero", significant="AllZero", SiteSpec=SS)))}
    }
    glmmy <- suppressWarnings(tryCatch(glm.nb(Count~Year + offset(log(BeforeDat$Hours)), link=log, data=BeforeDat), error=function(e){NA}))
    slope <- tryCatch(coef(glmmy, silent=TRUE)[[2]], error=function(e){NA})  #Extract the slope
    significant <- tryCatch(summary(glmmy)$coeff[-1,4], error=function(e){NA}) #Extract the significance
    random <- 0
    
    #Check for temporal autocorrelation. If it exists, run model again with a random factor on year and recalculate slope and significance. (first checking if the first model failed, if so, skip)
    if(length(glmmy)!=1){
      if(min(TACCheck(glmmy, BeforeDat, BABACI = "BA", TotalYearsBuffer)) < 0.05){
        glmmy <- suppressWarnings(tryCatch(glmer.nb(Count~Year + (1|Year) + offset(log(BeforeDat$Hours)), data=BeforeDat), error=function(e){NA}))
        slope <- tryCatch(summary(glmmy)$coeff[2,1], error=function(e){NA})  #Extract the slope
        significant <- tryCatch(summary(glmmy)$coeff[2,4], error=function(e){NA}) #Extract the significance
        random <- 1
      }
    }
    
    return(as.data.frame(cbind(slope, significant, SiteSpec=SS, random)))
  }, mc.cores=ncoresslopes) #
  
  #Bring all the data together
  BeforeSlopes <- as.data.frame(rbindlist(BeforeSlopes[unlist(sapply(BeforeSlopes, function(x) ifelse(ncol(x)==2, FALSE, TRUE)))]))
  ThresholdCheck <- ifelse(TotalYearsBuffer<7, 0.75, 0.95) #A lot of populations fail in this model if there isn't much data, so we lower the threshold in the following check for those cases
  if(length(unique(BeforeSlopes[complete.cases(BeforeSlopes),]$SiteSpec))/length(unique(NonZero$SiteSpec))<ThresholdCheck){stop("the beforeslopes function has gone wrong and you're losing a lot of site spec, possbly to do with glmer.nb")}
  if(nrow(subset(BeforeSlopes, random=="1"))==0){stop("there are no populations with temporal autocorrelation coming through, something's gone wrong with TACCheck or glmer")}
  BeforeSlopes$significant <- as.numeric(as.character(BeforeSlopes$significant))
  BeforeSlopes$slope <- as.numeric(as.character(BeforeSlopes$slope))
  NumYears <- dcast(as.data.table(BeforeData), SiteSpec~., length, value.var="Year") #Calculate the number of years of data that are  in each time series.
  BeforeSlopes <- merge(BeforeSlopes, NumYears, by="SiteSpec")
  
  #Get the values for "Before Slope" (1 for significantly positive, 0 for insignificant, -1 for significantly negative, "AllZero" for "All Zeroes)
  BeforeSlopes$BeforeSlope <- ifelse(is.na(BeforeSlopes$significant), 0, #If the significance is NA make it zero
                                     ifelse(BeforeSlopes$slope=="AllZero", "AllZero", #If there's AllZero, make it all zero
                                            ifelse(BeforeSlopes$significant<PThresh, ifelse(BeforeSlopes$slope>0,1,-1), 0)))
  
  names(BeforeSlopes)[names(BeforeSlopes)=="slope"] <- "SlopeEstimate"
  names(BeforeSlopes)[names(BeforeSlopes) == "significant"] <- "SlopeSignificance"
  BeforeSlopes$CI <- str_split_fixed(BeforeSlopes$SiteSpec, "[_]", 2)[,2]
  BeforeSlopes$SiteSpec <- str_split_fixed(BeforeSlopes$SiteSpec, "[_]", 2)[,1]
  
  if(actualvalues==TRUE){
    BeforeSlopes <- rbindlist(list(BeforeSlopes[,c("SiteSpec", "CI", "BeforeSlope", "SlopeEstimate", "SlopeSignificance")], AllZero), fill=TRUE)
  } else {
    BeforeSlopes <- rbind(BeforeSlopes[,c("SiteSpec", "CI", "BeforeSlope")], AllZero)
  }
  if(nrow(BeforeSlopes)!=nrow(BeforeSlopes[complete.cases(BeforeSlopes[,c("SiteSpec", "CI", "BeforeSlope")]),])){stop("there are NAs in before slopes")}
  return(BeforeSlopes)
}

#This function is used to identify time series that are "All Zero". Is zero threshold is one of the values defined in the latin hypercube sampling, a value between 0.6 and 0.9. If the proportion of zeroes in the time series if above the threshold, it is classified as "All Zeroes" (see Wauchope et al., 2021, Trends in ecology and evolution for an explanation of the rationale behind this)
#The function recieves time series ("Count Data) requiring at minimum a column for SiteCode, Species, SiteSpec, Year, Count, Hours (the effort term for CBC data) and CI (1 for protected sites, 0 for unprotected sites)
#It returns a dataframe with a column for "SiteSpec" and a column called "AllZero" which is either "AllZero" (for AllZero cases) or "Counts" (for non all zero cases)
#ZeroThresh is a value between 0.6 and 0.9, defined by latin hypercube sampling
AllZeroes <- function(CountData, ZeroThresh){
  #If the minimum value in the time series is not 0, jump straight to defining them as counts
  if(min(CountData$Count)==0){
    ZeroCounts <- dcast(as.data.table(subset(CountData, Count==0)), SiteSpec + CI~., length, value.var="Count") #Count the number of zero counts for each population/CI combination
    names(ZeroCounts) <- c("SiteSpec", "CI", "ZeroCounts")
    AllCounts <- dcast(as.data.table(CountData), SiteSpec + CI ~., length, value.var="Count") #Count the number of counts (Zero or no) for eaching population/CI combination
    names(AllCounts) <- c("SiteSpec", "CI", "AllCounts")
    
    AllCounts <- as.data.frame(merge(AllCounts, ZeroCounts, by=c("SiteSpec", "CI"), all=T)) #Merge the two together
    if(nrow(AllCounts)!=nrow(AllCounts[complete.cases(AllCounts),])){
      AllCounts[is.na(AllCounts$ZeroCounts),]$ZeroCounts <- 0
    }
    AllCounts$AllZero <- ifelse(AllCounts$ZeroCounts/AllCounts$AllCounts<=ZeroThresh, "Counts", "AllZero") #Calculate the proportion of zeroes for each population, if above the threshold "AllZero", if equal to or below, "Counts"
    AllCounts <- AllCounts[,c("SiteSpec", "CI", "AllZero")]
    
  } else { #If there's no zeroes, they're all counts
    AllCounts <- unique(CountData[,c("SiteSpec", "CI")])
    AllCounts$AllZero <- "Counts"
    AllCounts$SiteSpec <- as.character(AllCounts$SiteSpec)
  }
  return(AllCounts)
}

#These three functions clean and prepare data (matched data in the case of CI and BACI) for analysis of PA effectiveness. 
#Each recieves BA, BACI, or CI data in the form of the count data (including the columns for sitecode, species, sitespec, year, count, hours, etc, plus covariates. And in the case of BACI and CI a match ID)
#They also recieve the "ZeroThresh" - one of the values defined in the latin hypercube sampling, a value between 0.6 and 0.9. If the proportion of zeroes in the time series if above the threshold, it is classified as "All Zeroes" (see Wauchope et al., 2021, Trends in ecology and evolution for an explanation of the rationale behind this)
#PThresh is included for a supplmentary analysis that relaxes the p<0.05 cut off for significance, to see if we can better detected PA impact with slight trends
CleanBA <- function(BAData, ZeroThresh, PThresh, TotalYearsBuffer){
  BA <- BAData
  
  AfterSlopes <- CalculateBeforeSlopes(subset(BA, BA==1), ZeroThresh, parallelise=TRUE, actualvalues = FALSE, PThresh, TotalYearsBuffer) #This calculates the trend of the BA data in the after period, and returns it as 1 (positive trend), 0 (insignificant trend) or -1 (negative trend) (The before slope has already been calculated in the matching phase)
  names(AfterSlopes) <- c("SiteSpec", "CI", "AfterSlope")
  BA <- merge(BA, AfterSlopes, by=c("SiteSpec", "CI"), all=T)
  if(nrow(BA)!=nrow(BAData)){stop("rows lost")}
  
  BA$AfterSlope2 <- ifelse(BA$AfterSlope=="AllZero", "AllZero", "Counts") #Create a field that just marks the after slope as either all zero or counts
  BA$BeforeSlope2 <- ifelse(BA$BeforeSlope=="AllZero", "AllZero", "Counts") #Create a field that just marks the before slope as either all zero or counts
  
  BA <- subset(BA, AfterSlope2!="AllZero"|BeforeSlope2!="AllZero") #Remove cases where before the Before and After periods are AllZero (i.e. basically the species is functionally not present at the site either before or after)
  
  BA[is.na(BA$Hours),]$Hours <- 1 #For all the IWC sites, make "Hours" 1 (which means for modelling it just disappears, which we want as for IWC data effort is standardised)
  
  BA$ModelCat <- ifelse(BA$BeforeSlope2=="AllZero" | BA$AfterSlope=="AllZero", "Categorise", "Model") #Split the data into either "categorise" or "model". Categorise data has all zeroes, and so is just categorised as either an immigration or extinction. Model data is then used later to calculate the slope change from before to after. 
  return(BA) #Return the data
}
CleanBACI <- function(BACIData, ZeroThresh, PThresh, TotalYearsBuffer){
  BACI <- BACIData
  
  #Get slopes/all zeroes before and after
  BeforeSlopes <- unique(BACI[,c("SpecMatch", "BeforeSlope")])[complete.cases(unique(BACI[,c("SpecMatch", "BeforeSlope")])),]
  BACI$BeforeSlope <- NULL
  BACI <- merge(BACI, BeforeSlopes, all=T)
  
  AfterSlopes <- CalculateBeforeSlopes(subset(BACI, BA==1), ZeroThresh, parallelise=TRUE, actualvalues = FALSE, PThresh, TotalYearsBuffer)
  names(AfterSlopes) <- c("SiteSpec", "CI", "AfterSlope")
  BACI <- merge(BACI, AfterSlopes, all=T)
  if(nrow(BACI)!=nrow(BACIData)){stop("rows lost")}

  #Remove cases with all zeroes across the board
  BACI$BeforeSlope2 <- ifelse(BACI$BeforeSlope=="AllZero", "AllZero", "Counts")
  BACI$AfterSlope2 <- ifelse(BACI$AfterSlope=="AllZero", "AllZero", "Counts")
  
  AllZeroBACI <- melt(as.data.table(unique(BACI[,c("SpecMatch", "SiteSpec", "BeforeSlope2", "AfterSlope2")])), id.vars=c("SpecMatch", "SiteSpec"))
  AllZeroBACI <- dcast(AllZeroBACI, SpecMatch~value, length, value.var="variable")
  AllZeroBACI$ModelCat <- ifelse(AllZeroBACI$AllZero==4, "Remove", ifelse(AllZeroBACI$AllZero>0, "Categorise", "Model")) #Split the data into either "categorise" or "model". Categorise data has all zeroes, and so is just categorised as either an immigration or extinction. Model data is then used later to calculate the slope change from before to after. 
  if(nrow(AllZeroBACI)/length(unique(BACI$SpecMatch))!=1){stop("SpecMatchLost")}
  BACI <- merge(BACI, AllZeroBACI[,c("SpecMatch", "ModelCat")])
  BACI <- subset(BACI, ModelCat!="Remove")
  return(BACI)
}
CleanCI <- function(CIData, ZeroThresh, PThresh, TotalYearsBuffer){
  CI <- CIData
  AfterSlopes <- AllZeroes(CI, ZeroThresh) #Check there are no AllZeroes in CI (we can't deal with them)
  if(length(unique(AfterSlopes$AllZero))!=1){stop("There are still all zeroes in CI")}
  if(unique(AfterSlopes$AllZero)!="Counts"){stop("There are still all zeroes in CI")}
  
  CISlopes <- CalculateBeforeSlopes(CI, ZeroThresh, parallelise=TRUE, actualvalues = FALSE, PThresh, TotalYearsBuffer) #Calculate the slope of the CI populations
  names(CISlopes) <- c("SiteSpec", "CI", "ProtSlope")
  CI <- merge(CI, CISlopes, by=c("SiteSpec", "CI"), all=T)
  if(nrow(CI)!=nrow(CIData)){stop("rowslost!")}
  
  CI$ModelCat <- "Model" #All CI data is model, there is no categorise (cf BA, BACI)
  return(CI)
}

#These three functions categorise every population in the BACI, BA and CI datasets in the categories defined in Figre 3, and Extended Data Figures 1 and 2 (respectively)
#They each receive the data output from the "Clean" functions above, plus a filepath where output should be saved
#PThresh is included for a supplmentary analysis that relaxes the p<0.05 cut off for significance, to see if we can better detected PA impact with slight trends
#TotalYearsBuffer required for the nested function that checks temporal autocorrelation.

BACIPopulationCategorise <- function(BACI, MatchedFP, PThresh, TotalYearsBuffer){
  ### First run models that calculate the average and trend change from before to after between control and intervention sites.
  BACIModel <- subset(BACI, ModelCat=="Model") #Take the "BACI" data that is for modelling (cf the all zero cases that are to be categorised)
  #First we run a null model - this compares models with and without the BA term and see which fit better. If the model *without* the BA term fits better, then the data doesn't show any change from before to after protection and so we categorise it as no impact.
  #We also use this chance to double check that the parallel assumption is still fulfilled when we run a full model. If the interaction term between CI and Year is significant (indicating non-parallel before trends), the population is discarded.
  #The function runs through by population pair (i.e. each matched protected and unprotected population)
  NullModel <- rbindlist(pbmclapply(unique(BACIModel$SpecMatch), function(x){
    print(x)
    BACIPop <- subset(BACIModel, SpecMatch==x) #Subset to the relevant "specmatch" i.e. protected and unprotected population pair
    BACIPop$Year <- BACIPop$Year-BACIPop$STATUS_YR #Centre the year around 0 (see Wauchope et al., 2021, Trends in Ecology and Evolution for an explanation for why)
    BACIPop <- BACIPop[order(BACIPop$CI, BACIPop$Year),]
    random <- 0
    
    NullModel <- tryCatch(glm.nb(Count~Year + CI + CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL}) #Run the null model without the BA term
    if(is.null(NullModel)){
      return(data.frame(SpecMatch=x, Model="Failed", Random = random))
    }
    #Check for temporal autocorrelation. If it exists, run model again with a random factor on year and recalculate slope and significance. (first checking if the first model failed, if so, skip)
    if(min(TACCheck(NullModel, BACIPop, BABACI = "BACI", TotalYearsBuffer)) < 0.05){
      NullModel <- tryCatch(glmer.nb(Count~Year + CI + CI*Year + (1|CI:Year) + offset(log(Hours)), data=BACIPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL}) #Run the null model without the BA term
      random <- 1
    }
    if(is.null(NullModel)){
      return(data.frame(SpecMatch=x, Model="Failed", Random = random))
    }
    
    if(random == 0){
      FullModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL}) #Run the full model with BA term
    } else {
      FullModel <- tryCatch(glmer.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + (1|CI:Year) + offset(log(Hours)), data=BACIPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL}) #Run the full model with BA term
    }
    
    if(is.null(FullModel)|tryCatch(nrow(summary(FullModel)$coeff), error=function(e){2})<8){
      return(data.frame(SpecMatch=x, Model="Failed", Random = random))
    }
    if(summary(FullModel)$coeff[7,4]<0.05){ #Mark cases that are not parallel as such
      return(data.frame(SpecMatch=x, Model="NotParallel", Random = random))
    }
    if(AIC(FullModel)<AIC(NullModel)){ #Compare the AIC of the null and full model. If the Null model fits better, mark as "Null"
      return(data.frame(SpecMatch=x, Model="Model", Random = random))
    } else{
      return(data.frame(SpecMatch=x, Model="Null", Random = random))
    }
  }, mc.cores=ncores)) #
  if(nrow(subset(NullModel, Random==1))==0){stop("there are no random specmatch in the BACI null models which is odd, implies model is breaking")}
  BACIModel2 <- subset(NullModel, Model=="Model") #Get those cases where the Null model did NOT fit better (we'll now run the full model on these)
  NullModel <- subset(NullModel, Model!="Model") #Get the other cases
  NullModel <- merge(NullModel, unique(subset(BACI, CI==1)[,c("SpecMatch", "AfterSlope")]), all.x=T) #Bring in the "After slopes" of the null cases, to categorise:
  NullModel[NullModel$Model=="Null"]$Model <- ifelse(NullModel[NullModel$Model=="Null"]$AfterSlope=="-1", "No Impact (population declining)", ifelse(NullModel[NullModel$Model=="Null"]$AfterSlope=="1", "No Impact (population increasing)", "No Impact (no population trend)"))
  
  #Now run the full model on the population pairs where the Null model *didn't* fit better, extract model output (this is not very efficient, should just extract from the above, but it is what it is)
  BACIByPopulation <- pbmclapply(unique(BACIModel2$SpecMatch), function(x){ 
    BACIPop <- subset(BACIModel, SpecMatch==x)
    BACIPop$Year <- BACIPop$Year-BACIPop$STATUS_YR
    BACIPop <- BACIPop[order(BACIPop$CI, BACIPop$Year),]
    
    SlopeModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL})
    random <- 0
    if(is.null(SlopeModel)){
      return(NULL)
    }
    if(min(TACCheck(SlopeModel, BACIPop, BABACI = "BACI", TotalYearsBuffer)) < 0.05){
      SlopeModel <- tryCatch(glmer.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + (1|CI:Year) + offset(log(Hours)), data=BACIPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL})
      random <- 1
    }
    if(is.null(SlopeModel)){
      return(NULL)
    }
    ModelOutput <- as.data.frame(summary(SlopeModel)$coeff)
    ModelOutput$Coef <- row.names(ModelOutput)
    ModelOutput$SpecMatch <- x
    ModelOutput$Class <- unique(BACIPop$BeforeSlope)
    ModelOutput$Random <- random
    return(ModelOutput)
  }, mc.cores=ncores)
  
  if(length(BACIByPopulation[!sapply(BACIByPopulation, is.null)])/length(BACIByPopulation)<0.95){stop("There are still NUll models coming through")}
  BACIByPopulation <- BACIByPopulation[!sapply(BACIByPopulation, is.null)]
  BACIByPopulation <- rbindlist(BACIByPopulation)
  if(nrow(subset(BACIByPopulation, Random==1))==0){stop("there are no random specmatch in the BACI models which is odd, implies model is breaking")}
  BACIByPopulation$Random <- NULL
  names(BACIByPopulation) <- c("Estimate", "Error", "z", "P", "Coef", "SpecMatch", "Class")
  if(nrow(subset(BACIByPopulation, Coef=="Year:CI" & P<0.05))/nrow(BACIByPopulation)>0.05){stop("There are still non parallel slopes coming through")}
  BACILevelSlope <- subset(BACIByPopulation, Coef=="BA:CI" | Coef=="Year:BA:CI") #Pull out the coefficients that tell us the level and slope change (now called immediate and trend change in the paper)
  dir.create(file.path(MatchedFP, "ModelOutput/"), showWarnings = FALSE)
  write.csv(BACILevelSlope, file=paste0(MatchedFP, "ModelOutput/BACILevelSlope.csv"), row.names=FALSE) #Write out the actual model estimates
  BACILevelSlope <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACILevelSlope.csv")) #Read back in (sometimes I need to start here to fix code)
  
  BACILevelSlope$Category <- ifelse(BACILevelSlope$P>PThresh, "Insig", ifelse(BACILevelSlope$Estimate>0, "Pos", "Neg")) #Class the estimates into insig, positive or negative
  BACILevelSlope <- dcast(as.data.table(BACILevelSlope), SpecMatch~Coef, value.var="Category") #Cast so that we have a column for level change and a column for slope change
  BACILevelSlope <- merge(BACILevelSlope, unique(subset(BACI, CI==1)[,c("SpecMatch", "AfterSlope")]), all.x=T) #Merge with the after slopes
  names(BACILevelSlope)[names(BACILevelSlope)=="BA:CI"] <- "Level"
  names(BACILevelSlope)[names(BACILevelSlope)=="Year:BA:CI"] <- "Slope"
  BACILevelSlope$Outcome <- ifelse(BACILevelSlope$Level =="Pos" | BACILevelSlope$Slope =="Pos", "Positive Impact", ifelse(BACILevelSlope$Level=="Insig" & BACILevelSlope$Slope=="Insig", ifelse(BACILevelSlope$AfterSlope=="-1", "No Impact (population declining)", ifelse(BACILevelSlope$AfterSlope=="1", "No Impact (population increasing)", "No Impact (no population trend)")), "Negative Impact")) #Categorise
  
  #Bring together all the modelled data
  BACIModelledOutcomes <- rbind(BACILevelSlope[,c("SpecMatch", "Outcome")], NullModel[,c("SpecMatch", "Model")], use.names=FALSE)
  if(nrow(BACIModelledOutcomes)/length(unique(BACIModel$SpecMatch))<0.95){stop("Some specmatch have been lost!")}
  if("Null" %in% unique(BACIModelledOutcomes$Outcome)){stop("there are NULLs")}
  BACIModelledOutcomes$OutcomeSummary <- BACIModelledOutcomes$Outcome
  BACIModelledOutcomes[BACIModelledOutcomes$OutcomeSummary=="Failed" | BACIModelledOutcomes$OutcomeSummary=="NotParallel" ,]$OutcomeSummary <- "Excluded"
  BACIModelledOutcomes[grepl(pattern = "No Impact", BACIModelledOutcomes$OutcomeSummary),]$OutcomeSummary <- "No Impact"
  
  ### Now categorise all the All Zeroes
  BACICat <- subset(BACI, ModelCat=="Categorise") #Get the "categorise" populations, i.e. those that don't need modelling, just categorisation
  
  #Do an annoying bit of data switching so that we have a dataframe with a column SpecMatch and then ones for BeforeControl (BC), BeforeIntervention (BI), AfterControl (AC) and AfterIntervention (AI), and where the cells say whether each of these has counts or All Zeroes
  BACICat <- unique(BACICat[,c("SpecMatch", "CI", "BA", "AfterSlope2", "BeforeSlope2")])
  BACICatBefore <- subset(BACICat, BA==0)[,c("SpecMatch", "CI", "BA", "BeforeSlope2")]
  names(BACICatBefore) <- c("SpecMatch", "CI", "BA","Slope")
  BACICatAfter <- subset(BACICat, BA==1)[,c("SpecMatch", "CI", "BA", "AfterSlope2")]
  names(BACICatAfter) <- c("SpecMatch", "CI", "BA","Slope")
  BACICat <- rbind(BACICatBefore, BACICatAfter)
  BACICat <- dcast(as.data.table(BACICat), SpecMatch~BA+CI, value.var="Slope")
  names(BACICat) <- c("SpecMatch", "BC", "BI", "AC", "AI")
  
  #Categorise the outcome based on where the all zeroes are
  BACICat$Outcome <- ifelse(BACICat$BC=="AllZero" & BACICat$BI=="AllZero" & BACICat$AC=="Counts" & BACICat$AI=="Counts", "ImmigratedBoth",
                            ifelse(BACICat$BC=="AllZero" & BACICat$BI=="AllZero" & BACICat$AC=="AllZero" & BACICat$AI=="Counts", "ImmigratedIn",
                                   ifelse(BACICat$BC=="AllZero" & BACICat$BI=="AllZero" & BACICat$AC=="Counts" & BACICat$AI=="AllZero", "ImmigratedOut",
                                          ifelse(BACICat$BC=="Counts" & BACICat$BI=="Counts" & BACICat$AC=="AllZero" & BACICat$AI=="AllZero", "ExtinctBoth",
                                                 ifelse(BACICat$BC=="Counts" & BACICat$BI=="Counts" & BACICat$AC=="AllZero" & BACICat$AI=="Counts","ExtinctOut", "ExtinctIn")))))
  
  #Give the outcome summaries (e.g. immigrated both is a case of no impact)
  BACICat$OutcomeSummary <- ifelse(BACICat$Outcome=="ExtinctBoth", "No Impact",
                                   ifelse(BACICat$Outcome=="ImmigratedBoth", "No Impact", 
                                          ifelse(BACICat$Outcome=="ExtinctOut" | BACICat$Outcome=="ImmigratedIn", "Positive Impact", "Negative Impact")))
  
  ### Combine, summarise and check we haven't lost anything
  BACISummarise <- rbind(BACIModelledOutcomes[,c("SpecMatch", "Outcome", "OutcomeSummary")], BACICat[,c("SpecMatch", "Outcome", "OutcomeSummary")])
  
  if(length(unique(BACIModel[!BACIModel$SpecMatch %in% BACISummarise$SpecMatch,]$SpecMatch))>0){
    LostSpecMatch <- data.frame("SpecMatch"=unique(BACIModel[!BACIModel$SpecMatch %in% BACISummarise$SpecMatch,]$SpecMatch), "Outcome"="Failed", "OutcomeSummary"="Excluded")
    BACISummarise <- rbind(BACISummarise, LostSpecMatch, fill=TRUE)
  }
  
  if(nrow(BACISummarise)!=length(unique(BACISummarise$SpecMatch))){stop("There are duplicates")}
  if(nrow(BACISummarise)/length(unique(BACI$SpecMatch))!=1){stop("We've lost some specmatch")}
  
  BACISummarise$Species <- str_split_fixed(BACISummarise$SpecMatch, "[.]", 2)[,1]
  BACISummarise <- merge(BACISummarise, unique(subset(BACI, CI==1)[,c("SpecMatch", "SiteCode")]))
  
  write.csv(BACISummarise, file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"), row.names=FALSE)
  return("Done")
}
BAPopulationCategorise <- function(BA, MatchedFP, PThresh, TotalYearsBuffer){
  ### First run models that calculate the mean and slope change from before to after
  BAModel <- subset(BA, ModelCat=="Model")
  BAModel$Year <- BAModel$Year-BAModel$STATUS_YR #Centre the year around 0 (see Wauchope et al., 2021, Trends in Ecology and Evolution for an explanation for why)
  
  #First we run a null model - this compares models with and without the BA term and see which fit better. If the model *without* the BA term fits better, then the data doesn't show any change from before to after protection and so we categorise it as no impact.
  #The function runs through by population
  NullModel <- rbindlist(pbmclapply(unique(BAModel$SiteSpec), function(x){
    #print(x)
    BAPop <- subset(BAModel, SiteSpec==x) #Subset to relevant population
    BAPop <- BAPop[order(BAPop$Year),]

    NullModel <- tryCatch(glm.nb(Count~Year + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL}) #Run null model
    random <- 0
    if(is.null(NullModel)){
      return(data.frame(SiteSpec=x, Model="Failed", Random=random))
    }
    
    #Check for temporal autocorrelation. If it exists, run model again with a random factor on year and recalculate slope and significance. (first checking if the first model failed, if so, skip)
    if(min(TACCheck(NullModel, BAPop, BABACI = "BA", TotalYearsBuffer)) < 0.05){
      NullModel <- tryCatch(glmer.nb(Count~Year + (1|Year) + offset(log(Hours)), data=BAPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL}) #Run null model
      random <- 1
    }
    
    if(is.null(NullModel)){
      return(data.frame(SiteSpec=x, Model="Failed", Random=random))
    }
    
    if(random == 0){
      FullModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL}) #Run full model
    } else {
      FullModel <- tryCatch(glmer.nb(Count~Year + BA + BA*Year + (1|Year) + offset(log(Hours)), data=BAPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL}) #Run full model
    }
 
    if(is.null(FullModel)){
      return(data.frame(SiteSpec=x, Model="Failed", Random=random))
    }
    if(AIC(FullModel)<AIC(NullModel)){ #Check which fits better
      return(data.frame(SiteSpec=x, Model="Model", Random=random))
    } else{
      return(data.frame(SiteSpec=x, Model="Null", Random=random))
    }
  }, mc.cores=ncores))
  if(nrow(subset(NullModel, Random==1))==0){stop("there are no random specmatch in the BA null models which is odd, implies model is breaking")}
  BAModel2 <- subset(NullModel, Model=="Model") #Get those cases where the full model fits better (to be modelled again in a second)
  NullModel <- subset(NullModel, Model!="Model") #Extract the Null models
  NullModel <- merge(NullModel, unique(BA[,c("SiteSpec", "AfterSlope")]), all.x=T) #Get the after slopes of these to then categorise:
  NullModel[NullModel$Model=="Null"]$Model <- ifelse(NullModel[NullModel$Model=="Null"]$AfterSlope=="-1", "No Impact (population declining)", ifelse(NullModel[NullModel$Model=="Null"]$AfterSlope=="1", "No Impact (population increasing)", "No Impact (no population trend)"))
  
  #Now run the full model on the populations where the Null model *didn't* fit better, extract model output
  BAByPopulation <- pblapply(unique(BAModel2$SiteSpec), function(x){
    BAPop <- subset(BAModel, SiteSpec==x)
    BAPop <- BAPop[order(BAPop$Year),]

    SlopeModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL})
    random <- 0
    if(is.null(SlopeModel)){
      return(NULL)
    }
    
    if(min(TACCheck(SlopeModel, BAPop, BABACI = "BA", TotalYearsBuffer)) < 0.05){
      SlopeModel <- tryCatch(glmer.nb(Count~Year + BA + BA*Year + (1|Year) + offset(log(Hours)), data=BAPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL}) #Run full model
      random <- 1
    }
    
    if(is.null(SlopeModel)){
      return(NULL)
    }
    
    ModelOutput <- as.data.frame(summary(SlopeModel)$coeff)
    ModelOutput$Coef <- row.names(ModelOutput)
    ModelOutput$SiteSpec <- x
    ModelOutput$Class <- unique(BAPop$BeforeSlope)
    ModelOutput$Random <- random
    return(ModelOutput)
  })

  if(length(BAByPopulation[!sapply(BAByPopulation, is.null)])/length(BAByPopulation)<0.95){stop("There are still NUll models coming through")}
  BAByPopulation <- rbindlist(BAByPopulation)
  if(nrow(subset(BAByPopulation, Random==1))==0){stop("there are no random sitespec in the BA models which is odd, implies model is breaking")}
  BAByPopulation$Random <- NULL
  names(BAByPopulation) <- c("Estimate", "Error", "z", "P", "Coef", "SiteSpec", "Class")
  BALevelSlope <- subset(BAByPopulation, Coef=="BA" | Coef=="Year:BA")
  
  dir.create(file.path(MatchedFP, "ModelOutput/"), showWarnings = FALSE)
  write.csv(BALevelSlope, file=paste0(MatchedFP, "ModelOutput/BALevelSlope.csv"), row.names=FALSE) #Write out data
  BALevelSlope <- read.csv(paste0(MatchedFP, "ModelOutput/BALevelSlope.csv"))
  
  BALevelSlope$Category <- ifelse(BALevelSlope$P>PThresh, "Insig", ifelse(BALevelSlope$Estimate>0, "Pos", "Neg")) #Categorise the estimates as pos, neg or insig
  BALevelSlope <- dcast(as.data.table(BALevelSlope), SiteSpec~Coef, value.var="Category") #Cast so there's a column for level and a column for slope
  names(BALevelSlope)[names(BALevelSlope)=="BA"] <- "Level"
  names(BALevelSlope)[names(BALevelSlope)=="Year:BA"] <- "Slope"
  
  BALevelSlope <- merge(BALevelSlope, unique(BA[,c("SiteSpec", "AfterSlope")])) #After in after slope, categorise:
  BALevelSlope$Outcome <- ifelse(BALevelSlope$Level =="Pos" | BALevelSlope$Slope =="Pos", "Positive Impact", ifelse(BALevelSlope$Level=="Insig" & BALevelSlope$Slope=="Insig", ifelse(BALevelSlope$AfterSlope=="-1", "No Impact (population declining)", ifelse(BALevelSlope$AfterSlope=="1", "No Impact (population increasing)", "No Impact (no population trend)")), "Negative Impact"))
  
  BAModelledOutcomes <- rbind(BALevelSlope[,c("SiteSpec", "Outcome")], NullModel[,c("SiteSpec", "Model")], use.names=FALSE)
  if(nrow(BAModelledOutcomes)/length(unique(BAModel$SiteSpec))<0.95){stop("Some specmatch have been lost!")}
  BAModelledOutcomes$OutcomeSummary <- BAModelledOutcomes$Outcome #Sumamrise the outcomes
  BAModelledOutcomes[BAModelledOutcomes$OutcomeSummary=="Failed",]$OutcomeSummary <- "Excluded"
  BAModelledOutcomes[grepl(pattern = "No Impact", BAModelledOutcomes$OutcomeSummary),]$OutcomeSummary <- "No Impact"
  
  ### All zeroes
  BACat <- subset(BA, ModelCat=="Categorise") #Get the cases of all zero to categorise (unless there aren't any)
  if(nrow(BACat)!=0){
    #Wangle the data so there's a column for SiteSpec, a column for before and a column for after, where the cells say whether it is counts or all zeroes
    BACat <- unique(BACat[,c("SiteSpec", "BA", "AfterSlope2", "BeforeSlope2")])
    BACatBefore <- subset(BACat, BA==0)[,c("SiteSpec", "BA", "BeforeSlope2")]
    names(BACatBefore) <- c("SiteSpec", "BA","Slope")
    BACatAfter <- subset(BACat, BA==1)[,c("SiteSpec", "BA", "AfterSlope2")]
    names(BACatAfter) <- c("SiteSpec", "BA","Slope")
    BACat <- rbind(BACatBefore, BACatAfter)
    
    BACat <- dcast(as.data.table(BACat), SiteSpec~BA, value.var="Slope")
    names(BACat) <- c("SiteSpec", "Before", "After")
    if(nrow(subset(BACat, Before=="AllZero" & After=="AllZero"))){stop("All zeroes before and after coming through")}
    BACat$Outcome <- ifelse(BACat$Before=="AllZero" & BACat$After=="Counts", "Immigrated", "Extinct") #Classify as immigration or extinction
    
    BACat$OutcomeSummary <- ifelse(BACat$Outcome=="Immigrated", "Positive Impact", "Negative Impact")
    BASummarise <- rbind(BAModelledOutcomes[,c("SiteSpec", "Outcome", "OutcomeSummary")], BACat[,c("SiteSpec", "Outcome", "OutcomeSummary")])
  } else {
    BASummarise <- BAModelledOutcomes[,c("SiteSpec", "Outcome", "OutcomeSummary")]
  }
  
  ### Combine, summarise and check we haven't lost anything
  if(length(unique(BAModel[!BAModel$SiteSpec %in% BASummarise$SiteSpec,]$SiteSpec))>0){
    LostSpecMatch <- data.frame("SiteSpec"=unique(BAModel[!BAModel$SiteSpec %in% BASummarise$SiteSpec,]$SiteSpec), "Outcome"="Failed", "OutcomeSummary"="Excluded")
    BASummarise <- rbind(BASummarise, LostSpecMatch, fill=TRUE)
  }
  
  if(nrow(BASummarise)!=length(unique(BASummarise$SiteSpec))){stop("There are duplicates")}
  if(nrow(BASummarise)/length(unique(BA$SiteSpec))!=1){stop("We've lost some specmatch")}
  
  BASummarise$Species <- str_split_fixed(BASummarise$SiteSpec, "[.]", 2)[,2]
  BASummarise$SiteCode <- str_split_fixed(BASummarise$SiteSpec, "[.]", 2)[,1]
  write.csv(BASummarise, file=paste0(MatchedFP, "ModelOutput/BACategories.csv"), row.names=FALSE)
  return("Done")
}
CIPopulationCategorise <- function(CI, CIMatchedFP, PThresh, TotalYearsBuffer){
  ### Run models that calculate the difference in slope between control and intervention populations
  CIModel <- subset(CI, ModelCat=="Model")
  CIModel$Year <- CI$Year-CI$STATUS_YR #Centre the year around 0 (see Wauchope et al., 2021, Trends in Ecology and Evolution for an explanation for why)
  
  CIByPopulation <- pbmclapply(unique(CIModel$SpecMatch), function(x){
    #print(x)
    CIPop <- subset(CIModel, SpecMatch==x)
    CIPop <- CIPop[order(CIPop$CI, CIPop$Year),]
    
    SlopeModel <- tryCatch(glm.nb(Count~Year + CI + CI*Year + offset(log(Hours)), link=log, data=CIPop), error=function(e){NULL})
    random <- 0
    if(is.null(SlopeModel)){
      return(NULL)
    }
    
    if(min(TACCheck(SlopeModel, CIPop, BABACI = "CI", TotalYearsBuffer)) < 0.05){
      SlopeModel <- tryCatch(glmer.nb(Count~Year + CI + CI*Year + (1|CI:Year) + offset(log(Hours)), data=CIPop, control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000))), error=function(e){NULL})
      random <- 1
    }
    
    if(is.null(SlopeModel)){
      return(NULL)
    }
    
    ModelOutput <- as.data.frame(summary(SlopeModel)$coeff)
    ModelOutput$Coef <- row.names(ModelOutput)
    ModelOutput$SpecMatch <- x
    ModelOutput$Random <- random
    return(ModelOutput)
  }, mc.cores=ncores) #
  
  #Get the coefficients etc
  if(length(CIByPopulation[!sapply(CIByPopulation, is.null)])/length(CIByPopulation)<0.95){stop("There are lots of NUll models coming through (CI)")}
  CIByPopulation <- rbindlist(CIByPopulation)
  if(nrow(subset(CIByPopulation, Random==1))==0){stop("there are no random specmatch in the CI models which is odd, implies model is breaking")}
  CIByPopulation$Random <- NULL
  names(CIByPopulation) <- c("Estimate", "Error", "z", "P", "Coef", "SpecMatch")
  CISlope <- subset(CIByPopulation, Coef=="Year" | Coef=="Year:CI")
  
  dir.create(file.path(CIMatchedFP, "ModelOutput/"), showWarnings = FALSE)
  write.csv(CISlope, file=paste0(CIMatchedFP, "ModelOutput/CISlope.csv"), row.names=FALSE)
  CISlope <- read.csv(paste0(CIMatchedFP, "ModelOutput/CISlope.csv"))
  CISlope$Category <- ifelse(CISlope$P>PThresh, "Insig", ifelse(CISlope$Estimate>0, "Pos", "Neg"))
  CISlope <- dcast(as.data.table(CISlope), SpecMatch~Coef, value.var="Category")
  names(CISlope)[names(CISlope)=="Year"] <- "UnprotSlope"
  names(CISlope)[names(CISlope)=="Year:CI"] <- "SlopeChange"
  
  #Classify the CI data
  CISlope <- merge(CISlope, unique(subset(CIModel, CI==1)[,c("SpecMatch","ProtSlope")]), by="SpecMatch", all=T)
  CISlope$Outcome <- ifelse(CISlope$Slope =="Pos", "Positive Impact", ifelse(CISlope$Slope=="Insig", ifelse(CISlope$ProtSlope=="-1", "No Impact (population declining)", ifelse(CISlope$ProtSlope=="1", "No Impact (population increasing)", "No Impact (no population trend)")), "Negative Impact"))
  
  #Check we haven't lost any specmatch (other than those where the model failed)
  CISlope[is.na(CISlope$Outcome),]$Outcome <- "Failed"
  if(nrow(CISlope)/length(unique(CIModel$SpecMatch))<0.95){stop("Some specmatch have been lost!")}
  if(length(unique(CIModel[!CIModel$SpecMatch %in% CISlope$SpecMatch,]$SpecMatch))>0){
    LostSpecMatch <- data.frame("SpecMatch" = unique(CIModel[!CIModel$SpecMatch %in% CISlope$SpecMatch,]$SpecMatch), "Outcome"="Failed", "OutcomeSummary"="Excluded")
    CISlope <- rbind(CISlope, LostSpecMatch, fill=TRUE)
  }
  CISlope$OutcomeSummary <- CISlope$Outcome
  CISlope[CISlope$Outcome=="Failed",]$OutcomeSummary <- "Excluded"
  CISlope[grepl(pattern = "No Impact", CISlope$OutcomeSummary),]$OutcomeSummary <- "No Impact"
  
  ### Combine and summarise
  CISummarise <- CISlope[,c("SpecMatch", "Outcome", "OutcomeSummary")]
  if(nrow(CISummarise)!=length(unique(CISummarise$SpecMatch))){stop("There are duplicates")}
  if(nrow(CISummarise)/length(unique(CI$SpecMatch))!=1){stop("We've lost some specmatch")}
  
  CISummarise$Species <- str_split_fixed(CISummarise$SpecMatch, "[.]", 2)[,1]
  CISummarise <- merge(CISummarise, unique(subset(CI, CI==1)[,c("SpecMatch", "SiteCode", "SiteSpec")]), by="SpecMatch")
  write.csv(CISummarise, file=paste0(CIMatchedFP, "ModelOutput/CICategories.csv"), row.names=FALSE)
  return("Done")
}

#A function to calculate the area of each of PA for assessing predictors of protected area effectiveness
#Function require "referencedat" which is a dataframe with a column for sitecode and a column for CI (at minimum). The function will return the area for all Sites in the SiteCode column
GetPAArea <- function(ReferenceDat){
  ReferenceDat <- subset(ReferenceDat, CI==1) #Subset to just protected sites
  ProtectedSites <- fread(paste0(ResultsFP, "ProtectedSites_Uncleaned.csv")) #Read in the initial file from the Matching Script when we first collected all the PA data
  ProtectedSites <- ProtectedSites[ProtectedSites$SiteCode %in% ReferenceDat$SiteCode] #Reduce protected sites to just the sites in reference dat
  ProtectedSites <- merge(ProtectedSites, unique(ReferenceDat[,c("STATUS_YR", "SiteCode")]), by="SiteCode") #Add the used sitecode to the data (just to check there aren't any errors where we've used not the earliest designated protected area)
  names(ProtectedSites)[names(ProtectedSites)=="STATUS_YR.y"] <- "UsedStatYr"
  names(ProtectedSites)[names(ProtectedSites)=="STATUS_YR.x"] <- "STATUS_YR"
  
  #Subset to cases where any other PAs were designated in the survey period
  if(nrow(subset(ProtectedSites, STATUS_YR<UsedStatYr))>0){stop("There are cases of earlier status years - note this error will remain till you rerun matching")} #Check for errors
  
  ProtectedSites <- subset(ProtectedSites, STATUS_YR<(UsedStatYr+TotalYearsBuffer))
  
  GetPercentageofManyProts <- FALSE #This is to check how many sites occur in areas with overlapping protected areas
  if(GetPercentageofManyProts){
    ManyProt <- dcast(ProtectedSites, SiteCode~., length)
    nrow(subset(ManyProt, .>1))/nrow(ManyProt)
  }
  
  OriginalArea <- unique(subset(ProtectedSites, STATUS_YR==UsedStatYr)[,c("SiteCode", "GIS_AREA")]) #Get the area of the first designated PA
  OriginalArea <- dcast(OriginalArea, SiteCode~., max, value.var="GIS_AREA")
  names(OriginalArea) <- c("SiteCode", "OriginalArea")
  ProtectedSites <- merge(ProtectedSites, OriginalArea)
  ProtectedSites$LargerThanOriginal <- ifelse(ProtectedSites$GIS_AREA>ProtectedSites$OriginalArea, 1, 0)
  
  #Remove later designated PAs that are smaller than the original (as by being smaller they don't matter)
  ProtectedSites <- subset(ProtectedSites, STATUS_YR==UsedStatYr | GIS_AREA!=OriginalArea & LargerThanOriginal==1)
  
  #Take the mean of any remaining sites
  ProtectedSitesMean <- dcast(ProtectedSites, SiteCode~., mean, value.var="GIS_AREA")
  names(ProtectedSitesMean) <- c("SiteCode", "MeanArea")
  ProtectedSitesMean <- merge(ProtectedSitesMean, unique(ProtectedSites[,c("SiteCode", "OriginalArea")]))
  ProtectedSitesMean$Diff <- abs(ProtectedSitesMean$MeanArea - ProtectedSitesMean$OriginalArea)/ProtectedSitesMean$OriginalArea
  ProtectedSitesMean <- subset(ProtectedSitesMean, Diff<0.1)
  
  PAArea <- ProtectedSitesMean[,c("SiteCode", "MeanArea")]
  names(PAArea) <- c("SiteCode", "PAArea")
  return(PAArea)
}

#A function to gather all predictor data to assess protected area effectiveness. It requires the InputFP (to save output), the output from the GetPAArea function (called "PAArea") and RefData, which is a dataframe with a column for sitecode and a column for CI (at minimum)
#It will return a dataframe with predictor estimates for each protected site
GetPredictors <- function(InputFP, PAArea, RefData){
  load(file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData")) #Read in the protected counts from the beginning of analysis, as this still has all the covariates
  SpeciesSiteCovs <- unique(ProtectedCountData[,c("SiteCode", "Species", "MigStatus", "Family", "Order", "GovMean", "STATUS_YR", "ISO3", "Country", "Latitude", "Longitude")]) #Reduce down to only the predictors we care about that are already in the data (the rest we'll get now)
  RefData <- subset(RefData, CI==1) #Reduce refdat to protected sites
  SpeciesSiteCovs <- SpeciesSiteCovs[SpeciesSiteCovs$SiteCode %in% RefData$SiteCode,] #Reduce to the sites in refdat
  
  #Get IUCN redlist data
  RedListData <- read.csv(paste0(ResultsFP, "RedListSpeciesData.csv")) #This file is created later in this code file
  RedListData$RedList <- factor(RedListData$RedList)
  RedListData$RedList <- recode_factor(RedListData$RedList, "LC"="Not Threatened", "NT"="Threatened", "VU"="Threatened", "EN"="Threatened", "CR"="Threatened") #Rename
  if(length(unique(RedListData$Species))/length(unique(SpeciesSiteCovs$Species))<0.95){stop("SpeciesLost")}
  SpeciesSiteCovs <- merge(SpeciesSiteCovs, RedListData) #Merge with the covariate data
  
  #Get body size data (used to be generation length data! Genlength is now body size, but we're gonna leave the label as "GenLength" to save headaches in switching all the code terms)
  # GenLength <- as.data.frame(fread(paste0(DataFP, "GenerationLength/genlengthbirdconsbiol2020table4.csv")))
  GenLength <- as.data.frame(fread(paste0(DataFP, "Bodysize/BirdFuncDat_UpdatedSpecies.csv"))) #This file is created later in this code file
  GenLength <- GenLength[GenLength$Species %in% SpeciesSiteCovs$Species,]
  if(length(unique(GenLength$Species))/length(unique(SpeciesSiteCovs$Species))<0.95){stop("SpeciesLost")}
  GenLength <- GenLength[,c("Species", "Diet-5Cat", "BodyMass-Value")]
  names(GenLength) <- c("Species", "Diet", "GenLength")
  SpeciesSiteCovs <- merge(SpeciesSiteCovs, GenLength, by="Species", all.x=T) #Merge with covariate data (kept in diet for some reason, this isn't used)
  
  #PA Area
  SpeciesSiteCovs <- merge(SpeciesSiteCovs, PAArea, all=T) #Merge in the PA area file
  
  #Re organise Anthrome data
  RefData$AnthRound <- as.numeric(as.character(RefData$AnthRound))
  AnthMode <- dcast(as.data.table(RefData), SiteCode~., Mode, value.var="AnthRound")
  names(AnthMode) <- c("SiteCode", "AnthMode")
  if(length(unique(AnthMode$SiteCode))/length(unique(SpeciesSiteCovs$SiteCode))<0.95){stop("SpeciesLost")}
  SpeciesSiteCovs <- merge(SpeciesSiteCovs, AnthMode)
  
  #Get sites that are Ramsar Sites or EU Special Protected Area - Birds Directive sites
  ProtectedSites <- fread(paste0(ResultsFP, "ProtectedSites_Uncleaned.csv"))
  ManagedSitesTerms <- c("ramsar", "Birds Directive")
  
  ManagedSites <- ProtectedSites[grepl(paste(ManagedSitesTerms, collapse="|"), ProtectedSites$DESIG_ENG, ignore.case=TRUE)]
  CheckUnmanaged <- ProtectedSites[!grepl(paste(ManagedSitesTerms, collapse="|"), ProtectedSites$DESIG_ENG, ignore.case=TRUE)]

  SpeciesSiteCovs$Managed <- "Other"
  SpeciesSiteCovs[SpeciesSiteCovs$SiteCode %in% ManagedSites$SiteCode,]$Managed <- "Managed"

  #Clean up
  SpeciesSiteCovs$AnthMode <- as.character(SpeciesSiteCovs$AnthMode)
  SpeciesSiteCovs$PAAreaLog <- log(SpeciesSiteCovs$PAArea)
  SpeciesSiteCovs$RedList <- factor(SpeciesSiteCovs$RedList, levels=c("Not Threatened", "Threatened"))
  SpeciesSiteCovs$MigStatus <- factor(SpeciesSiteCovs$MigStatus)
  SpeciesSiteCovs$MigStatus <- recode_factor(SpeciesSiteCovs$MigStatus, "Resident"="Resident", "Breeding"="Migrant", "Non-breeding"="Migrant", "Passage"="Migrant")
  SpeciesSiteCovs$Family <- factor(SpeciesSiteCovs$Family)
  SpeciesSiteCovs$AnthMode <- factor(SpeciesSiteCovs$AnthMode, levels=c('10', '20', '30', '40', '50', '60'))
  SpeciesSiteCovs$AnthMode <- recode_factor(SpeciesSiteCovs$AnthMode, '10'='Urban', '20'='Village', '30'='Croplands', '40'='Rangeland', '50'='Semi-natural', '60'='Wild')
  SpeciesSiteCovs$Managed <- factor(SpeciesSiteCovs$Managed, levels=c("Other", "Managed"))  
  SpeciesSiteCovs$ISO3 <- factor(SpeciesSiteCovs$ISO3)
  return(SpeciesSiteCovs)
}

#This function prepares data for running the cumulative link predictor models (which I labeled "CategoryModels" early on for some reason and never changed it)
#The function requires the data output from the "PopulationCategorise" functions, the data output from "GetPredictors" a Filepath of where to save output, and a logical for whether to save output or not
PrepCategoryModelData <- function(CategoryData, Predictors, SaveFP, SaveOutput=TRUE){
  CatModelData <- merge(CategoryData, Predictors, by=c("SiteCode", "Species")) #Add predictors to the category data
  CatModelData$OutcomeSummary <- as.character(CatModelData$OutcomeSummary)
  CatModelData$SiteCode <- as.factor(CatModelData$SiteCode)
  CatModelData <- subset(CatModelData, OutcomeSummary!="Excluded") #Remove any populations that were excluded
  CatModelData$OutcomeSummary <- factor(CatModelData$OutcomeSummary, levels = c("Negative Impact", "No Impact", "Positive Impact"))
  CatModelData$STATUS_YR <- CatModelData$STATUS_YR-min(CatModelData$STATUS_YR) #Centre years around 0
  CatModelData$GenLength <- log(CatModelData$GenLength) #Log gen length (now body size)
  dir.create(file.path(paste0(SaveFP, "Drivers/")), showWarnings = FALSE)
  if(SaveOutput==TRUE){
    write.csv(CatModelData, paste0(SaveFP, "Drivers/CategoryModelsData.csv"), row.names = FALSE)
  }
  return(CatModelData)
}

#This function creates the formula for running the category models. It requires the Model Data output from PrepCategoryModelData, AND a logical for PAArea for whether to include this in the model or not (we run models with and without this as a predictor because of the number of sites that don't have area data)
CategoryModelFunction <- function(ModelData, PAArea, SitePairDist, Redlist = FALSE, Summary2 = FALSE){
  ModelData$Note <- "" #Create a blank note column in the data, to log anything that happens in this function
  ModelData$Order <- factor(ModelData$Order)
  ##To account for spatial autocorrelation we will now assign sites to grid cells
  if(cluster==FALSE){
    Countries <- readOGR("/Users/hannahwauchope/Dropbox/Work/Data/GISLayers/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp")
  } else {
    Countries <- readOGR(paste0(DataFP, "WorldSHP/TM_WORLD_BORDERS-0.3-SSudan.shp"))
  }  
  CountriesMoll <- spTransform(Countries, MollCRS)
  
  WorldRas <- raster(ext=extent(CountriesMoll), resolution=SitePairDist*1000*2, crs=MollCRS) #Create a blank world raster with grid cells the size of SOMETHING
  WorldRas$Grid <- 1:length(WorldRas)
  
  sitecoordinates <- unique(ModelData[,c("SiteCode", "Latitude", "Longitude")])
  sitepoints <- cbind(sitecoordinates$Longitude, sitecoordinates$Latitude)
  sitepoints <- SpatialPointsDataFrame(sitepoints, sitecoordinates, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  sitepoints <- spTransform(sitepoints, MollCRS)
  sitecoordinates$Grid <- extract(WorldRas, sitepoints)
  
  ModelData <- merge(ModelData, sitecoordinates[,c("SiteCode", "Grid")], by="SiteCode")
  ModelData$Grid <- as.factor(ModelData$Grid)
  
  if(Redlist == TRUE){
    #Add in EU27 IUCN status and subset to relevant countries for the EU member states analysis
    EURL <- read.csv(paste0(DataFP, "EURedList/EuropeanRedList.csv"), header=TRUE)
    EuMembers <- read.csv(paste0(DataFP, "EURedList/EUMemberStates.csv"), header=FALSE)
    ModelData <- ModelData[ModelData$Country %in% EuMembers$V1,]
    ModelData <- merge(ModelData, EURL, by="Species")
    ModelData <- subset(ModelData, RedList_eu27!="NE")
    ModelData$RedList_eu27 <- factor(ModelData$RedList_eu27, levels=c("LC", "NT", "VU", "EN"))
    ModelData$RedList_eu27 <- recode_factor(ModelData$RedList_eu27, "LC" = "Least Concern", "NT" = "Threatened", "VU" = "Threatened", "EN" = "Threatened")
    ModelData$Grid <- factor(ModelData$Grid, levels=unique(ModelData$Grid))
    varstocheck <- c("MigStatus", "Managed", "Order", "AnthMode", "RedList_eu27")
  } else {
    varstocheck <- c("MigStatus", "Managed", "Order", "AnthMode")
  }
  ModelData <- droplevels(ModelData)
  for(var in varstocheck){ #For the categorical variables, remove any levels that have less than 5 entries, and if all have less than four then remove the variable entirely (And make a note of this in the note column)
    if(min(subset(as.data.frame(table(ModelData[,c(var)])), Freq>0)$Freq)<5){
      ModelData$Note <- paste0(ModelData$Note, paste0(subset(as.data.frame(table(ModelData[,c(var)])), Freq<5)$Var1, collapse=","), " from ", var, " removed_")
      ModelData <- ModelData[!ModelData[,c(var)] %in% subset(as.data.frame(table(ModelData[,c(var)])), Freq<5)$Var1,]
    }
    if(nrow(subset(as.data.frame(table(ModelData[,c(var)])), Freq>0))==1){
      ModelData$Note <- paste0(ModelData$Note, ", ", var, "removed entirely")
      ModelData[,c(var)] <- NULL
    }
  }
  
  if(PAArea==FALSE){ #If PAArea is false remove this from the data
    ModelData$PAAreaLog <- NULL
  }
  
  #Create the formula from what variables remain in the ModelData dataset (based on what was removed above)
  Response <- ifelse(Summary2 == TRUE, "OutcomeSummary2", "OutcomeSummary")
  
  if(Redlist == TRUE){
    Formula <- as.formula(paste(Response, paste(c(paste0(c(names(ModelData)[names(ModelData) %in% c("MigStatus", "GenLength", "GovMean", "AnthMode", "Managed", "Order", "PAAreaLog", "RedList_eu27")])), 
                                                          "(1|ISO3:SiteCode)", "(1|Grid:SiteCode)", "(1|Species)"), 
                                                        collapse = " + "), sep = " ~ "))
    FormulaNoGrid <- as.formula(paste(Response, paste(c(paste0(c(names(ModelData)[names(ModelData) %in% c("MigStatus", "GenLength", "GovMean", "AnthMode", "Managed", "Order", "PAAreaLog", "RedList_eu27")])), 
                                                                "(1|ISO3:SiteCode)", "(1|Species)"), 
                                                              collapse = " + "), sep = " ~ "))
  } else {
    Formula <- as.formula(paste(Response, paste(c(paste0(c(names(ModelData)[names(ModelData) %in% c("MigStatus", "GenLength", "GovMean", "AnthMode", "Managed", "Order", "PAAreaLog")])), 
                                                          "(1|ISO3:SiteCode)", "(1|Grid:SiteCode)", "(1|Species)"), 
                                                        collapse = " + "), sep = " ~ "))
    FormulaNoGrid <- as.formula(paste(Response, paste(c(paste0(c(names(ModelData)[names(ModelData) %in% c("MigStatus", "GenLength", "GovMean", "AnthMode", "Managed", "Order", "PAAreaLog")])), 
                                                          "(1|ISO3:SiteCode)", "(1|Species)"), 
                                                        collapse = " + "), sep = " ~ "))
  }
  return(list(ModelData, Formula, FormulaNoGrid)) #Return the data, formula
}

#Now, finally, a function to run the actual category models, that is the cumulative link models that correlate protected area effectiveness on each population with site and species covariates.
#It requires the output from "CategoryModelFunction", a name for the model, a logical for whether it should be saved, and if this is true, a Save Filepath
CategoryModel <- function(CatList, ModelName, SaveMe = FALSE, SaveFP=NULL){
  CatModel <- clmm(CatList[[2]], data=CatList[[1]]) #Run the CLMM
  ModelRun <- "WithGrid"
  if(is.na(summary(CatModel)[[1]][1,2])){
    CatModel <- clmm(CatList[[3]], data=CatList[[1]]) #Run the CLMM
    ModelRun <- "NoGrid"
  }
  if(SaveMe==TRUE){
    save(CatModel, file=paste0(SaveFP, "BACIAll.RData"))
  }
  ModelSummary <- as.data.frame(summary(CatModel)[[1]]) #Pull out model coefficients and save
  ModelSummary$Coef <- row.names(ModelSummary)
  ModelSummary$Model <- ModelName
  names(ModelSummary) <- c("Estimate", "Error", "Z", "P", "Coef", "Model")
  ModelSummary$Notes <- unique(CatList[[1]]$Note)
  ModelSummary$ModelRun <- ModelRun
  return(ModelSummary)
}

#This is a hideous function but it does the job. It is the function to create the Bar Graph as seen in Figure 3 and Extended Data Figures 1 and 2. The function scales the width of the bars to the number of sites each species occurs in, or the number of species in a site. 
BarPlotWidths <- function(WidthData, Log=TRUE, Cov, FactorLevels, xlabel, AngledxLabs=FALSE, ColourList, SuccessValues, MainResults, aspectrat, includelegend=TRUE, removexaxislabel=FALSE, AxisPosition="bottom", keeplegend=TRUE){
  WidthData <- as.data.frame(WidthData)
  WidthData$Cov <- WidthData[,c(paste0(Cov))]
  if(MainResults==TRUE | MainResults=="MainWithSpecies"){
    WidthDataCast <- as.data.frame(dcast(as.data.table(WidthData), Cov~variable, value.var="value"))
    WidthDataCast$SuccessSum <- rowSums(WidthDataCast[names(WidthDataCast) %in% SuccessValues])/rowSums(WidthDataCast[,c(2:ncol(WidthDataCast))])
    WidthDataCast <- WidthDataCast[order(WidthDataCast$SuccessSum, decreasing = TRUE),] #, decreasing = TRUE
  }
  
  BarWidths <- as.data.frame(dcast(as.data.table(WidthData), Cov~., sum, value.var="value"))
  names(BarWidths) <- c("Cov", "Width")
  WidthData <- merge(WidthData, BarWidths, by="Cov")
  WidthData$PercentVal <- WidthData$value/WidthData$Width
  if(Log==TRUE){
    WidthData$LogWidth <- log(WidthData$Width+1)
  } else {
    WidthData$LogWidth <- WidthData$Width
  }
  if(MainResults==TRUE | MainResults=="MainWithSpecies"){
    WidthData$Cov <- factor(WidthData$Cov, levels=WidthDataCast$Cov)
  }
  WidthData <- WidthData[order(WidthData$Cov),]
  WidthDataUnique <- unique(WidthData[,c("Cov", "LogWidth")])
  pos <- 0.5 * (cumsum(WidthDataUnique$LogWidth) + cumsum(c(0, WidthDataUnique$LogWidth[-length(WidthDataUnique$LogWidth)])))
  WidthDataUnique$Position <- pos
  WidthData <- merge(WidthData, WidthDataUnique, by=c("Cov", "LogWidth"))
  #WidthData$variable <- factor(WidthData$variable, levels=c(1:max(as.numeric(as.character(WidthData$variable)))))
  WidthData$variable <- factor(WidthData$variable, levels=FactorLevels)
  
  if(MainResults==TRUE){
    xtext <- element_blank()
  } else if(MainResults==FALSE) {
    xtext <- element_text(size=plotfontsize, colour="black", angle=ifelse(AngledxLabs==FALSE, 0, 45))
  } else if(MainResults=="MainWithSpecies"){
    xtext <- element_text(size=6, colour="black", angle=90)
  }
  
  xlabel <- Cov
  
  if(removexaxislabel==FALSE){
    xtitle <- element_text(size=plotfontsize, colour="black")
  } else {
    xtitle <- element_blank()
  }
  
  Plot <- ggplot(WidthData, aes(x=Position, y=PercentVal, fill=variable, colour=variable))+
    geom_bar(stat="identity", position="stack", width=WidthData$LogWidth, size=0, colour="transparent")+ 
    scale_fill_manual(values = ColourList, name=("Category"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
    scale_colour_manual(values = ColourList, name=("Category"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
    scale_x_continuous(expand = c(0, 0), breaks=WidthDataUnique$Position, labels=WidthDataUnique$Cov, position = AxisPosition) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("Proportion")+
    xlab(xlabel)+
    theme(aspect.ratio=aspectrat, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=12),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, size=12),
          legend.title=element_text(size=12),
          axis.text.x = xtext,
          axis.ticks = element_line(colour="black"),
          axis.title.x = xtitle,
          legend.text = element_text(size=12, colour="black"),
          legend.justification = "top",
          legend.position = ifelse(includelegend==TRUE,"right", "none"))
  
  return(list(Plot, WidthData))
}

#Define all the colours (Don't use all these anymore, but they're handy)
BlueCols <- colorRampPalette(c("#060437ff", "#3e6aabff", "#e0f3f8ff"))
RedCols <- colorRampPalette(c("#fffdbfff", "#d73027ff", "#680e18ff"))

cols <- colorRampPalette(c("lightslategrey", "aliceblue"))
ColourListPlot <- c(BlueCols(5)[1:3],cols(3)[[2]], "aliceblue", "lightgoldenrod1", "goldenrod1", RedCols(5)[3:5])
PosCol <- rgb(red=rowMeans(as.data.frame(col2rgb(BlueCols(5)[1:3])))[[1]], green=rowMeans(as.data.frame(col2rgb(BlueCols(5)[1:3])))[[2]], blue=rowMeans(as.data.frame(col2rgb(BlueCols(5)[1:3])))[[3]], maxColorValue = 255)
NegCol <- rgb(red=rowMeans(as.data.frame(col2rgb(RedCols(5)[3:5])))[[1]], green=rowMeans(as.data.frame(col2rgb(RedCols(5)[3:5])))[[2]], blue=rowMeans(as.data.frame(col2rgb(RedCols(5)[3:5])))[[3]], maxColorValue = 255)
NeuPosCol <- rgb(red=rowMeans(as.data.frame(col2rgb(c(cols(3)[[2]], "aliceblue"))))[[1]], green=rowMeans(as.data.frame(col2rgb(c(cols(3)[[2]], "aliceblue"))))[[2]], blue=rowMeans(as.data.frame(col2rgb(c(cols(3)[[2]], "aliceblue"))))[[3]], maxColorValue = 255)
NeuNegCol <- rgb(red=rowMeans(as.data.frame(col2rgb(c("lightgoldenrod1", "goldenrod1"))))[[1]], green=rowMeans(as.data.frame(col2rgb(c("lightgoldenrod1", "goldenrod1"))))[[2]], blue=rowMeans(as.data.frame(col2rgb(c("lightgoldenrod1", "goldenrod1"))))[[3]], maxColorValue = 255)
NeutCol <- "#ffd553ff"
OutcomeSummaryCols <- c("Negative Impact"= NegCol, "No Impact (population declining or extinct)" = ColourListPlot[[7]], "No Impact (population stable, increasing or has immigrated)" = ColourListPlot[[4]], "Positive Impact" = PosCol, "Excluded" = "grey80")

GetMeanColour <- function(ColVec){rgb(red=rowMeans(as.data.frame(col2rgb(ColVec)))[[1]], 
                                      green=rowMeans(as.data.frame(col2rgb(ColVec)))[[2]], 
                                      blue=rowMeans(as.data.frame(col2rgb(ColVec)))[[3]], maxColorValue = 255)}

ProtCol <- "#009F2B"
UnprotCol <- "#002400"
DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_21/")
Scenarios <- read.csv(paste0(ResultsFP, "LatinSquareSamples.csv"))
Scenarios$Scenario <- 1:nrow(Scenarios)
Scen <- 21

#### Get some final predictor data ####
### Get Redlist data
GetRedList <- "Done"

if(GetRedList=="NotDone"){
  load(paste0(DatabaseFP, "CountDataCleaned.RData"))
  AllSpecies <- unique(CountData$Species)
  rm(CountData)
  IUCN_REDLIST_KEY <- "412cdf89c6a8241d53e22f14d47c409c9a4baa87783ca3f3fc85dbc029228e2b"
  RedListData <- rbindlist(pbmclapply(AllSpecies, function(x){
    print(x)
    rl_search(name = x, key = IUCN_REDLIST_KEY)$result
  }, mc.cores=ncores))
  RedListData <- RedListData[,c("scientific_name", "category", "eoo_km2")]
  names(RedListData) <- c("Species", "RedList", "EOO")
  write.csv(RedListData, paste0(DatabaseFP, "RedListSpeciesData.csv"), row.names=FALSE)
}

### Get Body size data
GetBodySize <- "Done"
if(GetBodySize != "Done"){
  library(taxize)
  load(file=paste0(ResultsFP, "CountData.RData"))
  GenLength <- as.data.frame(fread(paste0(DataFP, "Bodysize/BirdFuncDat.csv")))
  SpeciesListAll <- as.data.frame(unique(subset(CountData, Database=="Waterbird")$Species))
  names(SpeciesListAll) <- "Species"
  
  SpeciesList <- as.data.frame(unique(subset(CountData, Database=="Waterbird")$Species)[!unique(subset(CountData, Database=="Waterbird")$Species) %in% GenLength$Scientific])
  names(SpeciesList) <- "Species"
  TreeofLifeID <- pbmclapply(SpeciesList$Species, function (x) tryCatch(as.data.frame(get_tolid_(as.character(x))), error=function(e){"NoMatch"}), mc.cores=8)
  
  names(TreeofLifeID) <- SpeciesList$Species
  
  #This function is from "WaterbirdDataCollateCleanRawData.R", I've just it to get any species synonym names so we can match the waterbird data to the body size data
  TreeofLifeCleaning <- function(TreeofLifeID, SpeciesList){
    #Pull out those with only one match, get all synonyms
    TreeofLifeIDOneMatch <- TreeofLifeID[lapply(TreeofLifeID, nrow)==1]
    TreeofLifeIDOneMatch <- rbindlist(lapply(1:length(TreeofLifeIDOneMatch), function(x){
      TOL <- TreeofLifeIDOneMatch[[x]]
      TOL$CountsSpecies <- names(TreeofLifeIDOneMatch)[x]
      names(TOL) <- c("UniqueName", "MatchedName", "IsApproxMatch", "IsSynonym", "NomenclatureCode", "Score", "Flags", "SuppressedFromSynth", "OTTID", "Rank", "Source", "CountsSpecies")
      return(TOL)
    }))
    
    TreeofLifeIDOneMatch <- rbindlist(lapply(unique(TreeofLifeIDOneMatch$CountsSpecies), function(x){
      ok <- subset(TreeofLifeIDOneMatch, CountsSpecies == x)
      Names <- as.data.frame(t(as.data.frame(c(unique(c(ok$CountsSpecies, ok$UniqueName, ok$MatchedName))))))
      row.names(Names) <- NULL
      if(ncol(Names)==1){
        names(Names) <- "Species"
      } else {
        names(Names) <- c("Species", paste0("Synonym_", 1:(ncol(Names)-1)))
      }
      return(Names)
    }), fill=TRUE)
    
    # Pull out those with many matches
    TreeofLifeIDManyMatch <- TreeofLifeID[lapply(TreeofLifeID, nrow)>1]
    TreeofLifeIDManyMatch <- rbindlist(lapply(1:length(TreeofLifeIDManyMatch), function(x){
      TOL <- TreeofLifeIDManyMatch[[x]]
      TOL$CountsSpecies <- names(TreeofLifeIDManyMatch)[x]
      names(TOL) <- c("UniqueName", "MatchedName", "IsApproxMatch", "IsSynonym", "NomenclatureCode", "Score", "Flags", "SuppressedFromSynth", "OTTID", "Rank", "Source", "CountsSpecies")
      return(TOL)
    }))
    
    TreeofLifeIDSynonyms <- subset(TreeofLifeIDManyMatch, IsApproxMatch=="FALSE")
    TreeofLifeIDSynonyms <- TreeofLifeIDSynonyms[ifelse(TreeofLifeIDSynonyms$UniqueName==TreeofLifeIDSynonyms$CountsSpecies,FALSE,TRUE),c(1,2)]
    names(TreeofLifeIDSynonyms) <- c("Synonym_1", "Species")
    TreeofLifeIDSynonyms <- rbindlist(lapply(unique(TreeofLifeIDSynonyms$Species), function(x){
      if(nrow(subset(TreeofLifeIDSynonyms, Species==x))==1){
        return(subset(TreeofLifeIDSynonyms, Species==x)[,c(2,1)])} else{
          Synonyms <- dcast(subset(TreeofLifeIDSynonyms, Species==x), Species~Synonym_1, value.var="Synonym_1")
          names(Synonyms) <- c("Species", paste0("Synonym_", 1:(ncol(Synonyms)-1)))
          return(Synonyms)
        }
    }), fill=TRUE)
    
    #Leave the approx matches to be manually checked
    
    ### Bring all the synonyms together
    SpeciesSynonyms <- rbind(TreeofLifeIDOneMatch, TreeofLifeIDSynonyms, fill=TRUE)
    
    SpeciesNoMatch <- as.data.frame(SpeciesList[!SpeciesList$Species %in% SpeciesSynonyms$Species,1])
    names(SpeciesNoMatch) <- "Species"
    SpeciesSynonyms <- rbind(SpeciesSynonyms, SpeciesNoMatch, fill=TRUE)
    
    if(nrow(SpeciesSynonyms)!=nrow(SpeciesList)){stop("Somethings gone wrong with synonym matching!")}
    
    #Find all cases where original name, updated name, TOLID_ID, or a synonym matches a shapefile
    SpeciesSynonyms$FinalName <- apply(SpeciesSynonyms, 1, function (x){
      if(x["Species"] %in% GenLength$Scientific){
        x["Species"]
      } else if(x["Synonym_1"] %in% GenLength$Scientific){
        x["Synonym_1"]
      } else if(x["Synonym_2"] %in% GenLength$Scientific){
        x["Synonym_2"]
      } else {"NoMatch"}
    })
    
    return(SpeciesSynonyms)
  }
  
  #Get any synonyms, cross check all names with Birds of the World, return a list of either the correct name (with original name + matched name, which will either be the same or a synonym), or just the incorrect name with "nomatch"
  SpeciesSynonyms <- TreeofLifeCleaning(TreeofLifeID, SpeciesList)
  SpeciesListAll <- merge(SpeciesListAll, SpeciesSynonyms[,c(1,4)], by="Species", all=T)
  SpeciesListAll[is.na(SpeciesListAll$FinalName),]$FinalName <- SpeciesListAll[is.na(SpeciesListAll$FinalName),]$Species
  SpeciesList2 <- subset(SpeciesListAll, FinalName=="NoMatch")
  #write.csv(SpeciesList2, "/Users/hannahwauchope/Desktop/BodySizeCheck.csv", row.names=FALSE)
  
  UpdatedSpecNames <- read.csv("/Users/hannahwauchope/Desktop/BodySizeCheck.csv")
  SpeciesListAll2 <- rbind(subset(SpeciesListAll, FinalName!="NoMatch"), UpdatedSpecNames)
  
  SpeciesListAllCheck <- as.data.frame(SpeciesListAll2[!unique(SpeciesListAll2$FinalName) %in% GenLength$Scientific,])
  GenLength <- merge(GenLength, SpeciesListAll2, by.x="Scientific", by.y="FinalName")
  names(GenLength)[names(GenLength)=="Scientific"] <- "BodySizeSpecies"
  write.csv(GenLength, paste0(DataFP, "Bodysize/BirdFuncDat_UpdatedSpecies.csv"))
}

#### Assess Matching, establish PA impact on each population/population pair, and the see how effectiveness correlates with predictors for BACI ####
Scenarios <- read.csv(paste0(ResultsFP, "LatinSquareSamples.csv")) #Read in the LHC parameters
Scenarios$Scenario <- 1:nrow(Scenarios)

#Loop through each of the LHC scenarios to assess matching, calculate effectiveness
for(Scen in c(1:nrow(Scenarios))){
  Sys.sleep(sample(1:50,1))
  if(!file.exists(paste0(ResultsFP, "Organise_", Scen, ".csv"))){
    write.csv(NULL, paste0(ResultsFP, "Organise_", Scen, ".csv"))
    #### Define Scenario based on LHC ####
    Database <- "Waterbird"
    CBCBuffer <- "CBCBuffer" #NoCBCbuffer #CBCBuffer
    PABuffer <- Scenarios[Scen,]$PABuffer
    ZeroThresh <- Scenarios[Scen,]$ZeroThresh
    TotalYearsBuffer <- Scenarios[Scen,]$TotalYearsBuffer #5, 10, 20
    MeasuredYearsBuffer <- Scenarios[Scen,]$MeasuredYearsBuffer #3, 7, 15
    ImputeFlag <- "ZeroesImputed"
    StDiffThresh <- Scenarios[Scen,]$StDiffThresh
    SitePairDist <- Scenarios[Scen,]$SitePairDist
    PropPairsThresh <- Scenarios[Scen,]$PropPairsThresh
    PThresh <- Scenarios[Scen,]$PThresh
    dir.create(file.path(paste0(ResultsFP, "LHC/")), showWarnings = FALSE)
    DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
    dir.create(file.path(DatabaseFP), showWarnings = FALSE)
    
    Sys.sleep(sample(1:20,1))
    
    #### BA ####
    print(paste0("BA", Scen))
    InputFP <- paste0(DatabaseFP, "BA/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    
    # Assess effectiveness
    print("Clean and Categorise BA")
    load(file=paste0(DatabaseFP, "BACI/ProtectedCountsReadyForMatching.RData")) #I know this says the BACI filepath, but this contains all the before after data BEFORE it was reduced for BACI matching
    BA <- CleanBA(ProtectedCountData, ZeroThresh, PThresh, TotalYearsBuffer) #Run the CleanBA function, and save the output (need it for later as well as here)
    save(BA, file=paste0(InputFP, "BA.RData"))
    load(file=paste0(InputFP, "BA.RData"))
    
    #Categorise populations (according to categories shown in Extended Data Figure 1)
    BAPopulationCategorise(BA, InputFP, PThresh, TotalYearsBuffer) #(this saves "BACategories.csv")

    #### CI ####
    print(paste0("CI", Scen))
    InputFP <- paste0(DatabaseFP, "CI/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    ### Assess matching CI
    #This section takes the matches we've run in the previous script and assesses how well pairs are matched, discarding any that aren't matched well enough. This is Step F of Extended Data Figure 6. 
    
    print("Assess Matching")
    load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData")) #Load the variables used for matching in this particular LHC run
    MatchingFinal <- as.data.frame(rbindlist(lapply(list.files(path=paste0(InputFP, "SpeciesMatch/"), full.names = TRUE), fread), use.names=TRUE)) #Read in all the matched data output from the Matching script
    MatchingFinal$SpecMatch <- paste0(MatchingFinal$Species, ".", MatchingFinal$MatchID) #Create a unique identified for each matched pair, by combining species and match id. e.g. SpecMatch = "Actitis hypoleucos.1" is the SpecMatch for that species at a particular pair of one protected and one unprotected site.
    
    #Clean the dataset 
    StatYrs <- unique(MatchingFinal[,c("SpecMatch", "STATUS_YR")])[complete.cases(unique(MatchingFinal[,c("SpecMatch", "STATUS_YR")])),] #Now that we've paired protected and unprotected sites, we can add the designation year to the unprotected sites to make sure we just take the after years data. (e.g. if the protected site was designated in 1999, we'll make the designation year 1999 for the unprotected site as well)
    MatchingFinal$STATUS_YR <- NULL
    MatchingFinal <- merge(MatchingFinal, StatYrs, all=T)
    MatchingFinal[is.na(MatchingFinal$Hours),]$Hours <- 1 #Make hours 1 for IWC data, this means it just doesn't count in the models (because IWC data is standardised for effort)
    
    #Reduce unprotected sites to right number of years after (and check that protected were already correct) REmember that because it's CI, we're only using data on the years AFTER designation
    RowCheck <- nrow(subset(MatchingFinal, CI==1))
    MatchingFinal <- subset(MatchingFinal, Year<= (STATUS_YR+TotalYearsBuffer) & Year>STATUS_YR)
    
    ##Checks
    #Protected sites have right years
    if(nrow(subset(MatchingFinal, CI==1))!=RowCheck){stop("There are some protected sites that don't have the right number of years, which means something is wrong with year subs pre matching")}
    
    #Unprotected have right years
    CIMeasuredYears <- dcast(as.data.table(MatchingFinal), SiteCode + CI~., length, value.var="Count")
    if(nrow(CIMeasuredYears)!=nrow(CIMeasuredYears[CIMeasuredYears$.>=MeasuredYearsBuffer])){stop("There are some CI sites that don't have the right number of years, which means something is wrong with year subs in matching")}
    
    Check <- dcast(as.data.table(MatchingFinal), SpecMatch + CI ~. , length, value.var="Year")
    if(max(Check$.)>TotalYearsBuffer | min(Check$.)<MeasuredYearsBuffer){stop("there is something amiss with the number of counts")}
    
    #Site codes
    if(nrow(MatchingFinal[complete.cases(MatchingFinal$SiteCode),])!=nrow(MatchingFinal)){stop("There are NAs in site codes")}
    
    ##Extract covariates for each site, whether protected or unprotected. We need to do this to calculate the Standardised Difference in Means - we get one value for each matching covariate for the years AFTER designation, because it's CI
    if(!file.exists(file=paste0(InputFP, "MatchingCovs.RData"))){
      MatchingCovariates <- rbindlist(pbmclapply(unique(MatchingFinal$SpecMatch), function(SpecMatchID){
        SpecSub <- subset(MatchingFinal, SpecMatch==SpecMatchID) #Subset to the relevant specmatch
        StatYr <- unique(SpecSub$STATUS_YR)[complete.cases(unique(SpecSub$STATUS_YR))] #Get the stat year
        ProtCovs <- GetMeanCovs(subset(SpecSub, CI==1), VariablesForMatchingByYear) #Get covariates for protected site
        ProtCovs$CI <- 1
        UnProtCovs <- GetMeanCovs(subset(SpecSub, CI==0), VariablesForMatchingByYear) #Get covariates for unprotected site
        UnProtCovs$CI <- 0
        AllCovs <- rbind(ProtCovs, UnProtCovs) #Bring together
        AllCovs$SpecMatch <- SpecMatchID
        return(AllCovs)
      }, mc.cores=ncores))
      save(MatchingCovariates, file=paste0(InputFP, "MatchingCovs.RData"))
    }
    
    load(file=paste0(InputFP, "MatchingCovs.RData")) #Load the covariate data (created just above)
    MatchingCovariates$Species <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,1] #Get species names back
    MatchingCovariates$MatchID <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,2] #Get match ID back
    names(MatchingCovariates)[names(MatchingCovariates) %in% "Treatment"] <- "CI"
    
    Distances <- unique(MatchingFinal[,c("SpecMatch", "MD")]) #Get the mahalanobis distances for each SpecMatch from the MatchingFinal dataset
    MatchingCovariates <- merge(MatchingCovariates, Distances, by="SpecMatch", all=T) #Add these to the covariate data
    if(nrow(MatchingCovariates)!=nrow(MatchingCovariates[complete.cases(MatchingCovariates$MD),])){stop("there are NAs)")}
    
    #Here's three functions to calculate the standardised difference in means and assess matching
    
    #This function calculates the absolute distance between protected and unprotected site covariates for a species. It's used in the CompareMatching Function (it's not longer used for matching, but I used it while comparing propensity score and MD matching, as discussed in the supp material)
    AbsDist <- function(dataset){
      Distance <- rbindlist(lapply(unique(dataset$Species), function(Spec){
        CovariateData <- as.data.frame(subset(dataset, Species==Spec))
        CovariateData <- CovariateData[,colnames(CovariateData) %in% c(VariablesForMatchingByYear, "CI", "MatchID")]
        AbsDiff <- rbindlist(lapply(unique(CovariateData$MatchID), function(ID){
          IDSub <- subset(CovariateData, MatchID==ID)
          if(nrow(IDSub)<2){
            return(NULL)
          }
          IDSub <- as.data.frame(t(as.data.frame(apply(IDSub, 2, function(x) abs(diff(as.numeric(x)))))))
          IDSub$MatchID <- ID
          row.names(IDSub) <- NULL
          IDSub$CI <- NULL
          return(IDSub)
        }))
        AbsDiff <- melt(AbsDiff, id.vars="MatchID")
        AbsDiff <- dcast(AbsDiff, variable~., mean, value.var="value")
        AbsDiff$Species <- Spec
        names(AbsDiff) <- c("Variable", "MeanDist", "Species")
        return(AbsDiff)
      }))
      return(Distance)
    }
    
    #This function calculates the standardised difference in means (it is used below). It takes a dataset with a column for Species, SpecMatch, CI and the covariates in VariablesForMatchingByYear. This should only take data for one species at a time (this is done in the CompareMathching function)
    StDiffMean <- function(dataset){
      SpecSites <- dataset
      if(nrow(SpecSites)<=2){ #If there's 2 or less sites (i.e. one protected, one protected) return nothing as that's not enough data. 
        return(NULL)
      }
      SiteCast <- do.call(cbind,lapply(c(1,0), function(x){ #Get the mean values for each variable in each treatment
        TheVariables <- as.data.frame(subset(SpecSites, CI==x))
        TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear]
        TheVariables <- apply(TheVariables, 2, as.numeric)
        return(as.data.frame(colMeans(TheVariables, na.rm = TRUE)))
      }))
      names(SiteCast) <- c("Treat", "Cont")
      SiteCast$MeanDiff <- abs(SiteCast$Treat - SiteCast$Cont) #Get the difference in means between treatment and control sites (numerator for the SDiM function)
      SiteCastVar <- do.call(cbind,lapply(c(1,0), function(x){ #Get the sd for each variable in each treatment
        TheVariables <- as.data.frame(subset(SpecSites, CI==x))
        TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear]
        TheVariables <- apply(TheVariables, 2, as.numeric)
        return(as.data.frame(colVars(TheVariables, na.rm = TRUE)))
      }))
      names(SiteCastVar) <- c("Treat", "Cont")
      SiteCast$VarComp <- sqrt((SiteCastVar$Treat + SiteCastVar$Cont)/2) #Get denominator for SDiM function
      SiteCast$d <- SiteCast$MeanDiff/SiteCast$VarComp #Calculate the SDiM
      SiteCast[is.na(SiteCast)] <- 0
      SiteCast$Cov <- row.names(SiteCast)
      SiteCast <- SiteCast[,c(5:6)]
      return(SiteCast)
    }
    
    #This function works through each species, and for that species iteratively removes the worst matched sites until the Standardised Difference in Means is below the threshold value for all covariates (specified by LHC, "DiffThresh"), or there are less than two matches left. 
    CompareMatching <- function(dataset, DiffThresh){
      dir.create(file.path(paste0(InputFP, "MatchSummaries/")), showWarnings = FALSE) #Create a folder to save the summary of which sites are remaining in the matched dataset for each species (once inadequate matches are removed)
      dir.create(file.path(paste0(InputFP, "MatchData/")), showWarnings = FALSE) #Create a folder to save the full final dataset for each species, containing all retained matched sites and the counts over the years etc. 
      DeleteOldFiles <- c(list.files(paste0(InputFP, "MatchSummaries/"), full.names=TRUE), list.files(paste0(InputFP, "MatchData/"), full.names=TRUE)) #Clear old files if this has run before
      sapply(DeleteOldFiles, unlink)  #Clear old files if this has run before
      
      #Now we run through by species, assess the SDiM, if it's not below the threshold for all covariates, remove the worst matched pair, and keep going until we get below the threshold or run out of sites. 
      Comp <- pbmclapply(unique(dataset$Species), function(i){ 
        print(i)
        MatchingSub <- subset(dataset, Species==i) #Subset data to relevant species
        NumPairs <- length(unique(MatchingSub$MatchID)) #Count how many matched pairs we have
        if(nrow(MatchingSub)<=4){ #If there are less than two matched pairs (i.e. four sites) the species is out
          return(NULL)
        }
        Iteration <- 1 #Mark the iteration we're on
        TheDistance <- StDiffMean(MatchingSub) #Calculate the SDiM for each covariate
        TheDistance$Iteration <- Iteration #Add the iteration
        Threshold <- subset(TheDistance, d<DiffThresh) #Subset TheDistance to only covariates below the SDiM threshold
        if(nrow(Threshold)!=0){ #If we still have covariates left, count how many we have, and see if it's *all* of them, or just some of them. 
          ThresholdCast <- dcast(as.data.table(Threshold), Iteration~Cov, value.var="d") #Create a Cast frame which is basically just one row, and a column for each covariate (With the SDiM as the values)
          if(ncol(ThresholdCast)!=length(VariablesForMatchingByYear)){ #If not all, covariates are below the threshold, make a blank entry for ThresholdCast
            ThresholdCast <- data.frame()
          }
        } else {
          ThresholdCast <- data.frame()
        }
        
        while(nrow(ThresholdCast)==0 & nrow(MatchingSub)>4){ #While ThresholdCast is blank, and there are still at least 4 sites, remove the matched pair with the largest mahalanobis distance, and calculate SDiM again, keep looping till we achieve this
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
          ThresholdCast <- dcast(as.data.table(Threshold), Iteration~Cov, value.var="d")
          ThresholdCast <- ThresholdCast[complete.cases(ThresholdCast),]
          if(ncol(ThresholdCast)!=length(VariablesForMatchingByYear)){
            ThresholdCast <- data.frame()
          }
        }
        if(nrow(ThresholdCast)>0){ #Once the while loop has finished, either ThresholdCast will still be blank (meaning that we ended up with less than four sites, and in that case return null) or it will have 1 row
          ThresholdCast <- as.data.frame(t(ThresholdCast)) #If it's done, transpose threshold cast so there's a column for "Covariate" and a column for "distance"
          names(ThresholdCast) <- "StDiff"
          ThresholdCast$Variable <- row.names(ThresholdCast)
          AbsoluteDist <- as.data.frame(AbsDist(MatchingSub)) #Calculate the absolute distance (no longer really used)
          Summary <- merge(ThresholdCast, AbsoluteDist, by="Variable")
          Summary$PropPairs <- length(unique(MatchingSub$MatchID))/NumPairs
          write.csv(Summary, paste0(InputFP, "MatchSummaries/Summary_", i, ".csv"), row.names=FALSE) #Write out "Summary" which just tells us the standardised difference in means values/absolute distance for the species
          write.csv(MatchingSub, paste0(InputFP, "MatchData/MatchData_", i, ".csv"), row.names=FALSE) #Write out "MatchingSub" which tells us which sites remain as matched pairs for each species
          return(i)
        } else {
          return(NULL)
        }
      }, mc.cores=ncores)
    }
    
    CompMatch <- CompareMatching(MatchingCovariates, StDiffThresh) #Run through and iteratively remove matches until we have the matched set for each species
    
    SummariesFinal <- rbindlist(lapply(list.files(path=paste0(InputFP, "MatchSummaries/"), full.names=TRUE), fread)) #Read in that summary material
    nrow(subset(SummariesFinal, PropPairs>PropPairsThresh))/nrow(SummariesFinal) #See for how many species we have at least x% of protected populations with a match (one of the values defined in the LHC)
    SummariesFinal <- subset(SummariesFinal, PropPairs>PropPairsThresh) #Subset to only those species
    MatchDataFinal <- rbindlist(lapply(list.files(path=paste0(InputFP, "MatchData/"), full.names=TRUE), fread)) #Read in the actual data on matched sites/species
    MatchingFinalCleaned <- MatchingFinal[MatchingFinal$SpecMatch %in% MatchDataFinal$SpecMatch,] #Subset our full dataset (of counts per year etc) down to only the sites/species that remain in the final matched datasets
    MatchingFinalCleaned <- MatchingFinalCleaned[MatchingFinalCleaned$Species %in% unique(SummariesFinal$Species),] #Subset to just the species where we more than the threshold level of matched populations. 

    save(MatchingFinalCleaned, file=paste0(InputFP, "MatchingFinalCleaned.RData")) #Save
    
    if(length(unique(MatchingFinalCleaned$Species))<3){ #If there are less than 3 species we stop here because we can't run the predictor models
      write.csv(NULL, paste0(ResultsFP, "Organise_", Scen, "_ONLYONESPECIES_CI.csv"))
    } else {
      ### Clean Matched Data, calculate protected area impact
      print("clean and model data")
      load(file=paste0(InputFP, "MatchingFinalCleaned.RData")) #Read in matched data 
      CI <- CleanCI(MatchingFinalCleaned, ZeroThresh, PThresh, TotalYearsBuffer) #Clean the data
      save(CI, file=paste0(InputFP, "CI.RData")) #Save output (used later)
      load(file=paste0(InputFP, "CI.RData"))
      CIPopulationCategorise(CI, InputFP, PThresh, TotalYearsBuffer) #Categorise CI populations into the categories shown in Extended Data Figure 2
    }
    #### BACI ####
    print(paste0("BACI", Scen))
    InputFP <- paste0(DatabaseFP, "BACI/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    BAInputFP <- paste0(DatabaseFP, "BA/")
    CIInputFP <- paste0(DatabaseFP, "CI/")
    
    #### BACI Assess Matching ####
    #See comments on CI for this section, as they are basically identical. The only difference is that in BACI we assess matching on covariates in the BEFORE years rather than the AFTER years. 
    
    print("Assess Matching")
    load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))
    MatchingFinal <- as.data.frame(rbindlist(lapply(list.files(path=paste0(InputFP, "SpeciesMatch/"), full.names = TRUE), fread), use.names=TRUE))
    MatchingFinal$SpecMatch <- paste0(MatchingFinal$Species, ".", MatchingFinal$MatchID)
    
    #Clean the dataset 
    StatYrs <- unique(MatchingFinal[,c("SpecMatch", "STATUS_YR")])[complete.cases(unique(MatchingFinal[,c("SpecMatch", "STATUS_YR")])),]
    MatchingFinal$STATUS_YR <- NULL
    MatchingFinal <- merge(MatchingFinal, StatYrs, all=T)
    MatchingFinal$BA <- ifelse(MatchingFinal$Year>MatchingFinal$STATUS_YR, 1, 0)
    MatchingFinal[is.na(MatchingFinal$Hours),]$Hours <- 1
    
    #Reduce unprotected sites to right number of years before and after (and check that protected were already correct)
    RowCheck <- nrow(subset(MatchingFinal, CI==1))
    MatchingFinal <- subset(MatchingFinal, Year<= (STATUS_YR+TotalYearsBuffer) & Year>(STATUS_YR-TotalYearsBuffer))
    
    ##Checks
    #Protected sites have right years
    if(nrow(subset(MatchingFinal, CI==1))!=RowCheck){stop("There are some protected sites that don't have the right number of years, which means something is wrong with year subs pre matching")}
    
    #Unprotected have right years
    BACIMeasuredYears <- dcast(as.data.table(MatchingFinal), SiteCode + CI~BA, length, value.var="Count")
    if(nrow(BACIMeasuredYears)!=nrow(BACIMeasuredYears[BACIMeasuredYears$'0'>=MeasuredYearsBuffer & BACIMeasuredYears$'1'>=MeasuredYearsBuffer,])){stop("There are some BACI sites that don't have the right number of years, which means something is wrong with year subs in matching")}
    
    Check <- dcast(as.data.table(MatchingFinal), SpecMatch + CI + BA ~. , length, value.var="Year")
    if(max(Check$.)>TotalYearsBuffer | min(Check$.)<MeasuredYearsBuffer){stop("there is something amiss with the number of counts")}
    
    #Site codes
    if(nrow(MatchingFinal[complete.cases(MatchingFinal$SiteCode),])!=nrow(MatchingFinal)){stop("There are NAs in site codes")}
    
    ##Extract covariates for each Site
    if(!file.exists(file=paste0(InputFP, "MatchingCovs.RData"))){
      MatchingCovariates <- rbindlist(pbmclapply(unique(MatchingFinal$SpecMatch), function(SpecMatchID){
        SpecSub <- subset(MatchingFinal, SpecMatch==SpecMatchID)
        StatYr <- unique(SpecSub$STATUS_YR)[complete.cases(unique(SpecSub$STATUS_YR))]
        ProtCovs <- GetMeanCovs(subset(SpecSub, CI==1 & Year<= StatYr), VariablesForMatchingByYear)
        ProtCovs$CI <- 1
        UnProtCovs <- GetMeanCovs(subset(SpecSub, CI==0 & Year<= StatYr), VariablesForMatchingByYear)
        UnProtCovs$CI <- 0
        AllCovs <- rbind(ProtCovs, UnProtCovs)
        AllCovs$SpecMatch <- SpecMatchID
        return(AllCovs)
      }, mc.cores=ncores)) #Get the covariate values for predesignation years (now that we know which PA each unprotected site is matched to)
      save(MatchingCovariates, file=paste0(InputFP, "MatchingCovs.RData"))
    }
    
    load(file=paste0(InputFP, "MatchingCovs.RData"))
    MatchingCovariates$Species <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,1]
    MatchingCovariates$MatchID <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,2]
    names(MatchingCovariates)[names(MatchingCovariates) %in% "Treatment"] <- "CI"
    
    Distances <- unique(MatchingFinal[,c("SpecMatch", "MD")])
    MatchingCovariates <- merge(MatchingCovariates, Distances, by="SpecMatch", all=T)
    if(nrow(MatchingCovariates)!=nrow(MatchingCovariates[complete.cases(MatchingCovariates$MD),])){stop("there are NAs)")}
    
    AbsDist <- function(dataset){
      Distance <- rbindlist(lapply(unique(dataset$Species), function(Spec){
        CovariateData <- as.data.frame(subset(dataset, Species==Spec))
        CovariateData <- CovariateData[,colnames(CovariateData) %in% c(VariablesForMatchingByYear, "CI", "MatchID")]
        AbsDiff <- rbindlist(lapply(unique(CovariateData$MatchID), function(ID){
          IDSub <- subset(CovariateData, MatchID==ID)
          if(nrow(IDSub)<2){
            return(NULL)
          }
          IDSub <- as.data.frame(t(as.data.frame(apply(IDSub, 2, function(x) abs(diff(as.numeric(x)))))))
          IDSub$MatchID <- ID
          row.names(IDSub) <- NULL
          IDSub$CI <- NULL
          return(IDSub)
        }))
        AbsDiff <- melt(AbsDiff, id.vars="MatchID")
        AbsDiff <- dcast(AbsDiff, variable~., mean, value.var="value")
        AbsDiff$Species <- Spec
        names(AbsDiff) <- c("Variable", "MeanDist", "Species")
        return(AbsDiff)
      }))
      return(Distance)
    }
    StDiffMean <- function(dataset){
      SpecSites <- dataset
      if(nrow(SpecSites)<=2){
        return(NULL)
      }
      SiteCast <- do.call(cbind,lapply(c(1,0), function(x){ #Get the mean values for each variable in each treatment
        TheVariables <- as.data.frame(subset(SpecSites, CI==x))
        TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear]
        TheVariables <- apply(TheVariables, 2, as.numeric)
        return(as.data.frame(colMeans(TheVariables, na.rm = TRUE)))
      }))
      names(SiteCast) <- c("Treat", "Cont")
      SiteCast$MeanDiff <- abs(SiteCast$Treat - SiteCast$Cont)
      SiteCastVar <- do.call(cbind,lapply(c(1,0), function(x){ #Get the sd for each variable in each treatment
        TheVariables <- as.data.frame(subset(SpecSites, CI==x))
        TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear]
        TheVariables <- apply(TheVariables, 2, as.numeric)
        return(as.data.frame(colVars(TheVariables, na.rm = TRUE)))
      }))
      names(SiteCastVar) <- c("Treat", "Cont")
      SiteCast$VarComp <- sqrt((SiteCastVar$Treat + SiteCastVar$Cont)/2)
      SiteCast$d <- SiteCast$MeanDiff/SiteCast$VarComp
      SiteCast[is.na(SiteCast)] <- 0
      SiteCast$Cov <- row.names(SiteCast)
      SiteCast <- SiteCast[,c(5:6)]
      return(SiteCast)
    }
    CompareMatching <- function(dataset, DiffThresh){
      dir.create(file.path(paste0(InputFP, "MatchSummaries/")), showWarnings = FALSE)
      dir.create(file.path(paste0(InputFP, "MatchData/")), showWarnings = FALSE)
      DeleteOldFiles <- c(list.files(paste0(InputFP, "MatchSummaries/"), full.names=TRUE), list.files(paste0(InputFP, "MatchData/"), full.names=TRUE))
      sapply(DeleteOldFiles, unlink)
      
      Comp <- pbmclapply(unique(dataset$Species), function(i){
        print(i)
        MatchingSub <- subset(dataset, Species==i)
        NumPairs <- length(unique(MatchingSub$MatchID))
        if(nrow(MatchingSub)<=4){
          return(NULL)
        }
        Iteration <- 1
        TheDistance <- StDiffMean(MatchingSub)
        TheDistance$Iteration <- Iteration
        Threshold <- subset(TheDistance, d<DiffThresh)
        if(nrow(Threshold)!=0){
          ThresholdCast <- dcast(as.data.table(Threshold), Iteration~Cov, value.var="d")
          if(ncol(ThresholdCast)!=length(VariablesForMatchingByYear)){
            ThresholdCast <- data.frame()
          }
        } else {
          ThresholdCast <- data.frame()
        }
        
        while(nrow(ThresholdCast)==0 & nrow(MatchingSub)>4){
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
          ThresholdCast <- dcast(as.data.table(Threshold), Iteration~Cov, value.var="d")
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
          Summary$PropPairs <- length(unique(MatchingSub$MatchID))/NumPairs
          write.csv(Summary, paste0(InputFP, "MatchSummaries/Summary_", i, ".csv"), row.names=FALSE)
          write.csv(MatchingSub, paste0(InputFP, "MatchData/MatchData_", i, ".csv"), row.names=FALSE)
          return(i)
        } else {
          return(NULL)
        }
      }, mc.cores=ncores)
    }
    
    CompMatch <- CompareMatching(MatchingCovariates, StDiffThresh) #Run through and iteratively remove matches until we have the matched set for each species
    
    SummariesFinal <- rbindlist(lapply(list.files(path=paste0(InputFP, "MatchSummaries/"), full.names=TRUE), fread)) #Read in that summary material
    nrow(subset(SummariesFinal, PropPairs>PropPairsThresh))/nrow(SummariesFinal)
    SummariesFinal <- subset(SummariesFinal, PropPairs>PropPairsThresh)
    MatchDataFinal <- rbindlist(lapply(list.files(path=paste0(InputFP, "MatchData/"), full.names=TRUE), fread))
    MatchingFinalCleaned <- MatchingFinal[MatchingFinal$SpecMatch %in% MatchDataFinal$SpecMatch,]
    MatchingFinalCleaned <- MatchingFinalCleaned[MatchingFinalCleaned$Species %in% unique(SummariesFinal$Species),]
    save(MatchingFinalCleaned, file=paste0(InputFP, "MatchingFinalCleaned.RData"))
    if(length(unique(MatchingFinalCleaned$Species))<3){
      write.csv(NULL, paste0(ResultsFP, "Organise_", Scen, "_ONLYONESPECIES_BACI.csv"))
      next
    }
    #### BACI Clean Matched Data, get the data for covariates that predict effectiveness ####
    print("Clean and Model")
    load(file=paste0(InputFP, "MatchingFinalCleaned.RData")) #Load matched data
    if(nrow(MatchingFinalCleaned)==0){next}
    BACI <- CleanBACI(MatchingFinalCleaned, ZeroThresh, PThresh, TotalYearsBuffer) #Clean data
    save(BACI, file=paste0(InputFP, "BACI.RData")) #Save (used later)
    load(file=paste0(InputFP, "BACI.RData"))
    
    BACIPopulationCategorise(BACI, InputFP, PThresh, TotalYearsBuffer) #Categorise effecitveness outcomes as shown in Figure 3. 
    
    load(file=paste0(InputFP, "BACI.RData"))
    BACIArea <- GetPAArea(BACI) #Get area of protected areas
    BACIPredictors <- GetPredictors(InputFP, BACIArea, BACI) #Get other predictor data
    save(BACIPredictors, file=paste0(InputFP, "BACIPredictors.RData")) #Save
    
    #### BACI Category Models ####
    #BACI Full Models
    BACICategories <- read.csv(file=paste0(InputFP, "ModelOutput/BACICategories.csv")) #Load data on the category of protected area impact on each SpecMatch
    load(file=paste0(InputFP, "BACIPredictors.RData")) #Load predictor data
    BACIModelData <- PrepCategoryModelData(BACICategories, BACIPredictors, InputFP, SaveOutput = FALSE) #Prepare data for modelling
    BACIModelData$SiteSpec <- paste0(BACIModelData$SiteCode, ".", BACIModelData$Species)
    BACIModelData$Model <- "All"
    if(length(BACIModelData$Species)<10){ #If there are less than ten species in this LHC run the model doesn't run (too many covariates), so skip
      ModelAll <- NULL
      ModelNoPA <- NULL
      ModelRedlist <- NULL
    } else {
      ModelAll <- CategoryModel(CategoryModelFunction(BACIModelData[!is.na(BACIModelData$PAAreaLog),], PAArea=TRUE, SitePairDist), "All", SaveMe=TRUE, SaveFP = paste0(InputFP, "Drivers/")) #Run cumulative link models understanding relationship between predictors and protected area effectiveness
      ModelNoPA <- CategoryModel(CategoryModelFunction(BACIModelData, PAArea=FALSE, SitePairDist), "NoPA") #Run models again without protected area (as we lose a lot of Matched Pairs due to a lack of data)
      ModelRedlist <- CategoryModel(CategoryModelFunction(BACIModelData, PAArea=FALSE, SitePairDist, Redlist=TRUE), "RedList") #Run models using just eu27 data to assess variation in red list status (we don't include PA area in this to get the maximum number of species/sites, given we're interested in red list status in this model)
    }
    
    AllModels <- rbindlist(list(ModelNoPA, ModelAll, ModelRedlist)) #Bring together model output, save
    write.csv(AllModels, paste0(InputFP, "Drivers/CategoryModelsOutput.csv"), row.names = FALSE)
    
    AllModelData <- rbindlist(list(BACIModelData), fill=TRUE) #Bring together the data that goes into the models, save
    write.csv(AllModelData, paste0(InputFP, "Drivers/CategoryModelsData.csv"), row.names = FALSE)
  }
    #### Mark as Done ####
    write.csv(NULL, paste0(ResultsFP, "Organise_", Scen, "_DONE.csv"))
} #nrow(Scenarios)

#### Site Location Maps (Figure 1, Figure S1) ####
#First, make a map of site locations for every LHC scenario (this is not included in the paper, except for the Focal Analysis, Scenario 21, the map is shown in Figure S1)
for(Scen in c(1:21)){
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
  load(paste0(DatabaseFP, "BACI/BACI.RData")) #Load BACI data
  BACISites <- unique(BACI[,c("SiteCode", "CI", "Latitude", "Longitude")])
  BACISites$CI <- factor(BACISites$CI)
  BACISites$CI <- recode_factor(BACISites$CI, "1"="Protected", "0"="Unprotected")
  
  mapWorld <- borders("world", colour="gray80", fill="gray80")
  mapEurope <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey80", fill="grey80")
  ProtCols1 <- c("Protected"="#4aa257ff","Unprotected"="#664e94ff")

  World  <- ggplot(BACISites) + mapWorld +
    geom_point(aes(x=Longitude, y=Latitude, colour=CI), size=1, shape=1, stroke=0.8)+
    scale_colour_manual(values=ProtCols1)+ #, guide=FALSE
    coord_proj("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")+
    theme(
      aspect.ratio=0.53, 
      panel.grid = element_blank(), 
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0),
      legend.justification = "top",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background=element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.key=element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.text = element_text(size=7),
      legend.title = element_blank(),
      legend.position = c(0.05, 1))#c(0.935, 1)0.6
  #World
  Europe  <- ggplot(BACISites) + mapEurope +
    geom_point(aes(x=Longitude, y=Latitude, colour=CI), size=1, shape=1, stroke=0.8)+
    scale_colour_manual(values=ProtCols1, name="Site Type", guide=FALSE)+ #
    coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(size = 1, fill = NA, colour="black"),
      legend.justification = "top",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background=element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"))
  #Europe
  dir.create(paste0(FiguresFP, "SiteMaps/"))
  pdf(paste0(FiguresFP, "SiteMaps/BACISiteMap_Scen",Scen,".pdf"), (183/25.4), (95/25.4), bg="transparent") #Save as an image res=300, For some reason pdf ony recieves measurements in inches so ahve to convert
  grid.newpage()
  v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
  v2 <- viewport(width = 0.51, height = 0.47, x = -0.12, y = 0.07, just=c("left", "bottom")) #plot area for the inset map
  print(World,vp=v1) 
  print(Europe,vp=v2, mar=c(0,0,0,0))
  dev.off()
}

#### Run the maps again, but this time with separate maps for protected and unprotected
for(Scen in c(1:21)){
  for(ProtStat in c("Protected", "Unprotected")){
    DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
    load(paste0(DatabaseFP, "BACI/BACI.RData")) #Load BACI data
    BACISites <- unique(BACI[,c("SiteCode", "CI", "Latitude", "Longitude")])
    BACISites$CI <- factor(BACISites$CI)
    BACISites$CI <- recode_factor(BACISites$CI, "1"="Protected", "0"="Unprotected")
    BACISites <- subset(BACISites, CI==ProtStat)
    mapWorld <- borders("world", colour="gray80", fill="gray80")
    mapEurope <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey80", fill="grey80")
    ProtCols1 <- c("Protected"="#4aa257ff","Unprotected"="#664e94ff")
    
    if(ProtStat=="Protected"){
      LegendPosition <- c(0.05, 1)
    } else {
      LegendPosition <- "none"
    }
    
    World  <- ggplot(BACISites) + mapWorld +
      geom_point(aes(x=Longitude, y=Latitude, colour=CI), size=1, shape=1, stroke=0.8)+
      scale_colour_manual(values=ProtCols1)+ #, guide=FALSE
      coord_proj("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")+
      theme(
        aspect.ratio=0.53, 
        panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        legend.justification = "top",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.background=element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        legend.key=element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = LegendPosition)#c(0.935, 1)0.6
    World
    Europe  <- ggplot(BACISites) + mapEurope +
      geom_point(aes(x=Longitude, y=Latitude, colour=CI), size=1, shape=1, stroke=0.8)+
      scale_colour_manual(values=ProtCols1, name="Site Type", guide="none")+ #
      coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
      theme(
        panel.grid = element_blank(), 
        panel.border = element_rect(size = 1, fill = NA, colour="black"),
        legend.justification = "top",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.background=element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"))
    #Europe
    dir.create(paste0(FiguresFP, "SiteMaps/"))
    pdf(paste0(FiguresFP, "SiteMaps/BACISiteMap_Scen",Scen,"_", ProtStat,".pdf"), (183/25.4), (95/25.4), bg="transparent") #Save as an image res=300, For some reason pdf ony recieves measurements in inches so ahve to convert
    grid.newpage()
    v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
    v2 <- viewport(width = 0.51, height = 0.47, x = -0.12, y = 0.07, just=c("left", "bottom")) #plot area for the inset map
    print(World,vp=v1) 
    print(Europe,vp=v2, mar=c(0,0,0,0))
    dev.off()
  }
}

#Now, we want to make a map that colours sites by the number of LHC scenarios they occur in. 
#First get the site locations for each LHC run
AllScenSites <- rbindlist(pblapply(c(1:21), function(Scen){
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
  load(paste0(DatabaseFP, "BACI/BACI.RData"))
  BACISites <- unique(BACI[,c("SiteCode", "CI", "Latitude", "Longitude")])
  BACISites$CI <- factor(BACISites$CI)
  BACISites$CI <- recode_factor(BACISites$CI, "1"="Protected", "0"="Unprotected")
  BACISites$Scenario <- Scen
  return(BACISites)
}))

AllScenSitesCast <- dcast(AllScenSites, SiteCode + CI + Latitude + Longitude ~ . ,length) #Get the number of LHC scenarios for each site
AllScenSitesCast <- AllScenSitesCast[order(AllScenSitesCast$CI),]
AllScenSitesCast$. <- AllScenSitesCast$./21 #Divide the . column (i.e. the number of LHC scenarios) by 21, to get a proportion
names(AllScenSitesCast)[[5]] <- "Scenarios" #Rename
AllScenSitesCast$ScenBin <- cut(AllScenSitesCast$Scenarios, 42) #Cut the proportion column into 42 bins
AllScenSitesCast$ScenBinCI <- paste0(AllScenSitesCast$CI, AllScenSitesCast$ScenBin) #Make a column for CI plus the bins
AllScenOrder <- order(AllScenSitesCast$CI, AllScenSitesCast$ScenBin, decreasing=c(TRUE, FALSE), method="radix") #Get the right order so it plots right
AllScenSitesCast <- AllScenSitesCast[AllScenOrder,] #Reorder

#To make the legend we need to manually provide the colours values at various alpha values
highcolUP <- "#4f3384ff"      
highcolP <- "#056d16ff"      

#Get alphas for unprotected sites
AllAlphasUP <- c(sapply(seq(0.4, 1, length.out=20), function(amin){
  lowcol.hexUP <- as.hexmode(round(col2rgb(highcolUP) * amin + 255 * (1 - amin)))
  lowcolUP <- paste0("#",   sep = "",
                     paste(format(lowcol.hexUP, width = 2), collapse = ""))
  return(lowcolUP)
}), highcolUP)

#Get alphas for protected sites
AllAlphasP <- c(sapply(seq(0.4, 1, length.out=20), function(amin){
  lowcol.hexP <- as.hexmode(round(col2rgb(highcolP) * amin + 255 * (1 - amin)))
  highcolP <- paste0("#",   sep = "",
                     paste(format(lowcol.hexP, width = 2), collapse = ""))
  return(highcolP)
}), highcolP)

mapWorld <- borders("world", colour="gray80", fill="gray80")
mapEurope <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey80", fill="grey80")

#Plot world map
World  <- ggplot(AllScenSitesCast) + mapWorld +
  geom_point(aes(x=Longitude, y=Latitude, colour=ScenBinCI), size=1, shape=1, stroke=0.8)+
  scale_colour_manual(values=c(AllAlphasP, AllAlphasUP))+ #, guide=FALSE
  coord_proj("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")+
  theme(
    aspect.ratio=0.53, 
    panel.grid = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    legend.justification = "top",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.background=element_rect(fill = "white"),
    panel.background = element_rect(fill = "transparent"),
    legend.key=element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.text = element_text(size=25),
    legend.title = element_blank(),
    legend.position = "none")#c(0.935, 1)0.6 #c(0.07, 1)
World

#Plot Europe inset
Europe  <- ggplot(AllScenSitesCast) + mapEurope +
  geom_point(aes(x=Longitude, y=Latitude, colour=ScenBinCI), size=1, shape=1, stroke=0.8)+
  scale_colour_manual(values=c(AllAlphasP, AllAlphasUP))+ #, guide=FALSE
  coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(size = 1, fill = NA, colour="black"),
    legend.justification = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.background=element_rect(fill = "transparent", colour=NA),
    panel.background = element_rect(fill = "white"))
Europe

#Make legend
GLegend <- ggplot(unique(AllScenSitesCast[,c("CI", "ScenBin", "ScenBinCI")]), aes(CI,ScenBin))+
  geom_tile(aes(fill=ScenBinCI))+
  ylab("Number of analyses")+
  scale_fill_manual(values=c(AllAlphasP, AllAlphasUP))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), labels=c(1, rep("", 19), 21))+
  coord_flip()+
  theme(legend.position="none", aspect.rat=0.4, 
        axis.title=element_text(size=7),
        plot.margin=margin(t=10,b=10,l=10),
        axis.text.x = element_text(size=5, colour="black"),
        axis.text.y = element_text(size=5, colour="black"),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())
GLegend

#Save
pdf(paste0(FiguresFP, "SiteMaps/AllScen.pdf"), (183/25.4), (95/25.4), bg="transparent") #Save as an image res=300, For some reason pdf ony recieves measurements in inches so ahve to convert
grid.newpage()
v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2 <- viewport(width = 0.51, height = 0.47, x = -0.12, y = 0.07, just=c("left", "bottom")) #plot area for the inset map
v3 <- viewport(width = 0.69, height = 0.23, x = 0.53, y = 0.23, just=c("left", "top"))
print(World,vp=v1) 
print(Europe,vp=v2, mar=c(0,0,0,0))
print(GLegend,vp=v3, mar=c(0,0,0,0))
dev.off()

#### Sankey Plots (Figure 2, Extended Data Figure 1) ####
print("Sankey")

### Focal Analysis Sankey Plot (Figure 2)
#This function creates one Sankey panel for comparing either BA or CI to BACI
#It requires "Comp" which is either "BA" or "CI", and logicals for whether to include excluded populations, a legend, and a y axis label
SankeyPlot <- function(Comp, InclNotAssessed=FALSE, InclLegend=TRUE, YaxisLabel=TRUE, Yaxistext=TRUE){
  Scen <- 21
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
  
  BACISummarise <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv")) #Read in BACI data
  BACISummarise$SiteSpec <- paste0(BACISummarise$SiteCode, ".", BACISummarise$Species)
  
  CompSummarise <- read.csv(file=paste0(DatabaseFP, Comp, "/ModelOutput/", Comp, "Categories.csv")) #Read in BA or CI data
  if(InclNotAssessed==FALSE){ #Subset BA/CI and BACI to the same sitespec
    CompSummarise <- CompSummarise[CompSummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, CompSummarise$SiteSpec),]
    BACISummarise <- BACISummarise[BACISummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, CompSummarise$SiteSpec),]
  }
  
  #Reduce datasets to relevant columns
  BACISummariseShort <- BACISummarise[,c("SiteSpec", "OutcomeSummary")]
  names(BACISummariseShort) <- c("SiteSpec", "BACI")
  
  CompSummarise <- CompSummarise[,c("SiteSpec", "OutcomeSummary")]
  names(CompSummarise) <- c("SiteSpec", "Comp")
  
  SummariseMerge <- merge(BACISummariseShort, CompSummarise, by = "SiteSpec", all=T) #Merge datasets together
  if(InclNotAssessed==TRUE){
    SummariseMerge[is.na(SummariseMerge)] <- "Excluded"
    SummariseMerge$BACI <- as.character(SummariseMerge$BACI)
    SummariseMerge$Comp <- as.character(SummariseMerge$Comp)
    SummariseMerge[SummariseMerge$BACI=="Excluded",]$BACI <- "Not Assessed"
    SummariseMerge[SummariseMerge$Comp=="Excluded",]$Comp <- "Not Assessed"
  } else {
    SummariseMerge <- subset(SummariseMerge, BACI!="Excluded" & Comp!="Excluded")
  }
  
  SummariseMergeCast <- dcast(as.data.table(SummariseMerge), Comp + BACI ~., length, value.var="SiteSpec") #Get the number of populations in each combination of BA/CI and BACI
  names(SummariseMergeCast) <- c("Comp", "BACI", "Freq")
  if(InclNotAssessed==TRUE){
    SummariseMergeCast <- subset(SummariseMergeCast, BACI!="Not Assessed" | Comp!="Not Assessed")
  }
  SummariseMergeCast$FreqProp <- SummariseMergeCast$Freq/sum(SummariseMergeCast$Freq) #Get the number of populations in each combination as a proportion
  
  SummariseMergeCast$Comp <- factor(SummariseMergeCast$Comp, levels=c("Not Assessed", "Negative Impact", "No Impact",  "Positive Impact")) #Make levels a factor
  SummariseMergeCast$BACI <- factor(SummariseMergeCast$BACI, levels=c("Not Assessed", "Negative Impact", "No Impact",  "Positive Impact")) #Make levels a factor
   
  OutcomeSummaryCols <- c("Negative Impact"= NegCol, "No Impact" = NeutCol, "Positive Impact" = PosCol) #Define colours
  
  #Plot Specs
  LegPos <- ifelse(InclLegend=="TRUE", "top", "none")
  if(YaxisLabel==TRUE){
    YaxLab <- element_text(size=plotfontsize) #family="Helvetica Neue", 
  } else {
    YaxLab <- element_blank()
  }
  if(Yaxistext==TRUE){
    Ytext <- element_text(size=plotfontsize-2) #family="Helvetica Neue", 
  } else {
    Ytext <- element_blank()
  }

  SummariseMergeCast2 <- as.data.table(SummariseMergeCast)
  SummariseMergeCast2$ID <- 1:nrow(SummariseMergeCast2)
  SummariseMergeCast2 <- melt(SummariseMergeCast2, id.vars=c("ID", "Freq", "FreqProp"))
  SummariseMergeCast2$variable <- recode_factor(SummariseMergeCast2$variable, "Comp" = Comp, "BACI" = "BACI")
  #Run a geom_alluvium plot
  Sankey <- ggplot(SummariseMergeCast2, aes(x=variable, y=FreqProp, stratum=value, alluvium=ID, fill=value))+
    geom_alluvium(reverse=FALSE, alpha=0.9)+
    geom_stratum(width = 0.125, reverse = FALSE, colour="black", size=0.4)+
    scale_fill_manual(values=OutcomeSummaryCols, name="Outcome")+
    ylab("Proportion of populations")+
    #xlab(paste0(Comp, " to BACI"))+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0.0025, 0.0025)) +
    guides(fill = guide_legend(reverse=TRUE))+
    theme(panel.background = element_blank(),
          aspect.ratio = 1.3,
          panel.grid = element_blank(),
          legend.position = LegPos,
          legend.text = element_text(size=plotfontsize), #family="Helvetica Neue", 
          legend.title = element_blank(),
          axis.text.x = element_text(size=plotfontsize-1),
          axis.ticks.x = element_blank(),
          axis.title.y = YaxLab,
          axis.text.y = Ytext,
          axis.title.x = element_blank(), #family="Helvetica Neue", 
          panel.border = element_rect(size = 1, fill = NA),
          legend.justification = c(0.41, 0),
          legend.key.size = unit(0.4, 'cm'),
          plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  return(list(Sankey, SummariseMerge))
}

#Use BA to just get the legend exported
BASankeyLegend <- SankeyPlot("BA", InclLegend=TRUE) 

#Get BA Sankey
BASankey <- SankeyPlot("BA", InclLegend=FALSE, YaxisLabel=TRUE) 
nrow(BASankey[[2]]) #Number of populations
#Get CI sankey
CISankey <- SankeyPlot("CI", InclLegend=FALSE, YaxisLabel=FALSE, Yaxistext = FALSE)
nrow(CISankey[[2]]) #Number of populations
#Save plot with BA and CI, as per Figure 2. 
grobs <- ggplotGrob(BASankeyLegend[[1]])$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
pgrid <- plot_grid(BASankey[[1]], CISankey[[1]], ncol = 2, align="hv", labels = c("a", "b"), label_size = 8, label_x = 0.03, label_y = 1.013) #rel_widths = c(0.8, 0.7), 
p <- plot_grid(pgrid, legend, ncol = 1, rel_heights = c(0.7, .06))
ggsave(paste0(FiguresFP, "Fig2_Sankey.pdf"), p, device="pdf", width = 89, height = 60, units = "mm") #Save as a pdf) #, dpi=1000
ggsave(paste0(FiguresFP, "Fig2_Sankey.png"), p, device="png", width = 89, height = 60, units = "mm", dpi=300) #Save as a pdf) #, dpi=1000

### Sankey Sensitivity (ED Figure 4)
#Get the data on the proportion of changes from Positive to No impact etc (all combinations) for each LHC
SankeyAllRuns <- rbindlist(lapply(c(1:21), function(SankScen){
  print(SankScen)
  #This is the same function as above, just minus the plot (i.e. we export the data used in the plot)
  SankeyStatsFunc <- function(Comp, InclNotAssessed=FALSE, InclLegend=TRUE, YaxisLabel=TRUE){
    BACISummarise <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
    BACISummarise$SiteSpec <- paste0(BACISummarise$SiteCode, ".", BACISummarise$Species)
    
    CompSummarise <- read.csv(file=paste0(DatabaseFP, Comp, "/ModelOutput/", Comp, "Categories.csv"))
    if(InclNotAssessed==FALSE){
      CompSummarise <- CompSummarise[CompSummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, CompSummarise$SiteSpec),]
      BACISummarise <- BACISummarise[BACISummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, CompSummarise$SiteSpec),]
    }
    
    if(nrow(CompSummarise)==0){
      return(NULL)
    }
    
    #Make sankey plots
    BACISummariseShort <- BACISummarise[,c("SiteSpec", "OutcomeSummary")]
    names(BACISummariseShort) <- c("SiteSpec", "BACI")
    
    CompSummarise <- CompSummarise[,c("SiteSpec", "OutcomeSummary")]
    names(CompSummarise) <- c("SiteSpec", "Comp")
    
    SummariseMerge <- merge(BACISummariseShort, CompSummarise, by = "SiteSpec", all=T)
    if(InclNotAssessed==TRUE){
      SummariseMerge[is.na(SummariseMerge)] <- "Excluded"
      SummariseMerge$BACI <- as.character(SummariseMerge$BACI)
      SummariseMerge$Comp <- as.character(SummariseMerge$Comp)
      SummariseMerge[SummariseMerge$BACI=="Excluded",]$BACI <- "Not Assessed"
      SummariseMerge[SummariseMerge$Comp=="Excluded",]$Comp <- "Not Assessed"
    } else {
      SummariseMerge <- subset(SummariseMerge, BACI!="Excluded" & Comp!="Excluded")
    }
    return(SummariseMerge)
  }
  
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", SankScen, "/")
  
  #Get the numbers for BA and CI
  BASankey <- SankeyStatsFunc("BA", InclLegend=FALSE) 
  CISankey <- SankeyStatsFunc("CI", InclLegend=FALSE, YaxisLabel=FALSE)
  
  if(is.null(BASankey) | is.null(CISankey)){
    return(NULL)
  }
  #Reorganise data
  BASankeyStats <- BASankey
  BASankeyStats <- dcast(as.data.table(BASankeyStats), Comp + BACI ~., length, value.var="SiteSpec")
  names(BASankeyStats) <- c("Comp", "BACI", "Freq")
  ChangedBA <- subset(BASankeyStats, Comp!=BACI)
  PosToOtherBA <- subset(BASankeyStats, Comp=="Positive Impact" & BACI!="Positive Impact")
  NegToOtherBA <- subset(BASankeyStats, Comp=="Negative Impact" & BACI!="Negative Impact")
  NeuToOtherBA <- subset(BASankeyStats, Comp=="No Impact" & BACI!="No Impact")
  NeuToPosBA <- subset(BASankeyStats, Comp=="No Impact" & BACI=="Positive Impact")

  CISankeyStats <- CISankey
  CISankeyStats <- dcast(as.data.table(CISankeyStats), Comp + BACI ~., length, value.var="SiteSpec")
  names(CISankeyStats) <- c("Comp", "BACI", "Freq")
  ChangedCI <- subset(CISankeyStats, Comp!=BACI)
  PosToOtherCI <- subset(CISankeyStats, Comp=="Positive Impact" & BACI!="Positive Impact")
  NegToOtherCI <- subset(CISankeyStats, Comp=="Negative Impact" & BACI!="Negative Impact")
  NeuToOtherCI <- subset(CISankeyStats, Comp=="No Impact" & BACI!="No Impact")
  NeuToPosCI <- subset(CISankeyStats, Comp=="No Impact" & BACI=="Positive Impact")


  #Return a dataframe on the number of populations in each combination
  return(data.frame("Scenario" = rep(SankScen, 2), "BACI" = c("BA", "CI"),
                    "PropChanged"= c(sum(ChangedBA$Freq)/sum(BASankeyStats$Freq), sum(ChangedCI$Freq)/sum(CISankeyStats$Freq)),
                    "PosToOther" = c(sum(PosToOtherBA$Freq)/sum(subset(BASankeyStats, Comp=="Positive Impact")$Freq),
                                     sum(PosToOtherCI$Freq)/sum(subset(CISankeyStats, Comp=="Positive Impact")$Freq)),
                    "NegToOther" = c(sum(NegToOtherBA$Freq)/sum(subset(BASankeyStats, Comp=="Negative Impact")$Freq),
                                     sum(NegToOtherCI$Freq)/sum(subset(CISankeyStats, Comp=="Negative Impact")$Freq)),
                    "NeuToOther" = c(sum(NeuToOtherBA$Freq)/sum(subset(BASankeyStats, Comp=="No Impact")$Freq),
                                     sum(NeuToOtherCI$Freq)/sum(subset(CISankeyStats, Comp=="No Impact")$Freq)),
                    "NeutoPos" = c(sum(NeuToPosBA$Freq)/sum(subset(BASankeyStats, Comp=="No Impact")$Freq),
                                   sum(NeuToPosCI$Freq)/sum(subset(CISankeyStats, Comp=="No Impact")$Freq)),
                    "NPop" = c(nrow(BASankey), nrow(CISankey))))
}))

#Reformat dataframe for plotting
SankeyAllRuns <- melt(SankeyAllRuns, id.vars=c("Scenario", "BACI")) 
SankeyAllRuns$variable <- factor(SankeyAllRuns$variable)
SankeyAllRuns <- subset(SankeyAllRuns, variable!="NeutoPos")
SankeyAllRuns$variable <- recode_factor(SankeyAllRuns$variable, "PropChanged"="All populations", "PosToOther" = "Positive impact populations", "NeuToOther" = "No impact populations", "NegToOther" = "Negative impact populations", "NPop" = "Number of populations")

#Make the various plot facets as seen in ED Figure 4

#Total Change
All <- ggplot(subset(SankeyAllRuns, variable=="All populations"))+
  geom_boxplot(aes(x=BACI, y=value), outlier.shape=NA, fill="grey10", colour="grey30", alpha=0.2)+#fill=BACI, 
  geom_jitter(aes(x=BACI, y=value), colour="grey30", alpha=0.5, width=0.25)+
  geom_point(data=subset(SankeyAllRuns, Scenario==21 & variable=="All populations"), aes(x=BACI, y=value), size=3, colour="grey30")+
  xlab("All populations")+
  ylab("Proportion of populations")+
  scale_y_continuous(limits=c(0, 1), breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  theme(panel.background = element_blank(),
        aspect.ratio = 1.3,
        panel.grid = element_blank(),
        legend.text = element_text(size=plotfontsize), #family="Helvetica Neue", 
        legend.title = element_blank(),
        axis.title.y = element_text(size=plotfontsize, colour="black"), #family="Helvetica Neue",
        axis.title.x = element_text(size=plotfontsize, colour="black"),
        axis.text.x = element_text(size=6, colour="black"),
        axis.text.y = element_text(size=6, colour="black"),
        panel.border = element_rect(size = 1, fill = NA),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

#Populations changing from BACI = Positive
Pos <- ggplot(subset(SankeyAllRuns, variable=="Positive impact populations"))+
  geom_boxplot(aes(x=BACI, y=value), outlier.shape=NA, fill=PosCol, colour=PosCol, alpha=0.2)+#fill=BACI, 
  geom_jitter(aes(x=BACI, y=value), colour=PosCol, alpha=0.5, width=0.25)+
  geom_point(data=subset(SankeyAllRuns, Scenario==21 & variable=="Positive impact populations"), aes(x=BACI, y=value), size=3, colour=PosCol)+
  scale_y_continuous(limits=c(0, 1.01), breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  xlab("Positive impact populations")+
  ylab("Proportion")+
  theme(panel.background = element_blank(),
        aspect.ratio = 1.3,
        panel.grid = element_blank(),
        legend.text = element_text(size=plotfontsize), #family="Helvetica Neue", 
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=plotfontsize, colour="black"),
        axis.text.x = element_text(size=6, colour="black"),        
        axis.text.y = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))

#Populations changing from BACI = No impact
Neu <- ggplot(subset(SankeyAllRuns, variable=="No impact populations"))+
  geom_boxplot(aes(x=BACI, y=value), outlier.shape=NA, fill=NeutCol, colour=NeutCol, alpha=0.2)+#fill=BACI, 
  geom_jitter(aes(x=BACI, y=value), colour=NeutCol, alpha=0.5, width=0.25)+
  geom_point(data=subset(SankeyAllRuns, Scenario==21 & variable=="No impact populations"), aes(x=BACI, y=value), size=3, colour=NeutCol)+
  scale_y_continuous(limits=c(0, 1.01), breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  xlab("No impact populations")+
  theme(panel.background = element_blank(),
        aspect.ratio = 1.3,
        panel.grid = element_blank(),
        legend.text = element_text(size=plotfontsize), #family="Helvetica Neue", 
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=plotfontsize, colour="black"),
        axis.text.x = element_text(size=6, colour="black"),
        axis.text.y = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))

#Populations changing from BACI = Negative Impact
Neg <- ggplot(subset(SankeyAllRuns, variable=="Negative impact populations"))+
  geom_boxplot(aes(x=BACI, y=value), outlier.shape=NA, fill=NegCol, colour=NegCol, alpha=0.2)+#fill=BACI, 
  geom_jitter(aes(x=BACI, y=value), colour=NegCol, alpha=0.5, width=0.25)+
  geom_point(data=subset(SankeyAllRuns, Scenario==21 & variable=="Negative impact populations"), aes(x=BACI, y=value), size=3, colour=NegCol)+
  scale_y_continuous(limits=c(0, 1.01), breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
  xlab("Negative impact populations")+
  theme(panel.background = element_blank(),
        aspect.ratio = 1.3,
        panel.grid = element_blank(),
        legend.text = element_text(size=plotfontsize), #family="Helvetica Neue", 
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=plotfontsize, colour="black"),
        axis.text.x = element_text(size=6, colour="black"),
        axis.text.y = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))

pgrid <- plot_grid(All, Pos, Neu, Neg, ncol=4, align="hv", labels = c("a", "b", "c", "d"), label_size = 8, label_x=0.025, label_y=1)
ggsave(paste0(FiguresFP, "EDFig4_SankeySensitivity.pdf"), pgrid, 
       device="pdf", width = 183, height = 60, units = "mm", dpi=1000) #Save as a pdf)

#### Stacked Bar Plots (Figure 3, Extended Data Figures 3, 4) ####
print("Stacked Bar")
Scen <- 21
DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_21/")

### BACI ###
BACISummarise <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv")) #Read in data on BACI impacts per matched species/site matched pair
#BACISummariseSitesPerSpec <- dcast(as.data.table(unique(BACISummarise[,c("Species", "SiteCode")])), Species~., length, value.var="SiteCode") #Find the number of sites that each species occurs in (used to determine width of bars in plot) (this is now done in the function, leaving here just in case I've missed something hah)
ColourListPlot <- c("#0f0b3fff", "#485c92ff", "#759bc7ff", "#808d9cff", "#b9c1caff", "#f9faf7ff", "#ffeb88ff", "#ffcb48ff", "#dc5851ff", "#a71111ff", "#550f17ff") #Get the right colours (runs from dark blue to adrk red)
FactorLevels <- c("ImmigratedIn", "ExtinctOut", "Positive Impact", "ImmigratedBoth", "No Impact (population increasing)", "No Impact (no population trend)", "No Impact (population declining)", "ExtinctBoth","Negative Impact", "ImmigratedOut", "ExtinctIn") #Define what the categories are for BACI

#By Species
BACISummariseForPlot <- melt(dcast(as.data.table(subset(BACISummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SpecMatch"), id.vars="Species") #Count the number of sites in each outcome per species
BACIResultsPlot <- BarPlotWidths(BACISummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Positive Impact", "ImmigratedIn", "ExtinctOut"), MainResults=TRUE, aspectrat=0.2) #Run the barplotwidths function
BACIResultsPlot[[1]]
ggsave(paste0(FiguresFP, "Inkscape/BACIBars_BySpecies.png"), BACIResultsPlot[[1]], device="png", width = 350, height = 80, units = "mm", dpi=1000) #Save as a pdf) #Save

#Range of proportion positively impacted by species
RangeEst <- BACIResultsPlot[[2]]
RangeEst$Positive <- ifelse(RangeEst$variable=="ImmigratedIn" | RangeEst$variable=="ExtinctOut" | RangeEst$variable=="Positive Impact", "Yes", "No")
RangeEst <- dcast(as.data.table(RangeEst), Species ~ Positive, sum, value.var="value")
RangeEst$PropPos <- RangeEst$Yes/rowSums(RangeEst[,c("Yes", "No")])
c(max(RangeEst$PropPos), min(RangeEst$PropPos), mean(RangeEst$PropPos), sd(RangeEst$PropPos))

#By Site
BACISummariseForPlot <- melt(dcast(as.data.table(subset(BACISummarise, OutcomeSummary!="Excluded")), SiteCode~Outcome, length, value.var="SpecMatch"), id.vars="SiteCode") #Count the number of species in each outcome per site
BACIResultsPlot <- BarPlotWidths(BACISummariseForPlot, Log=FALSE, Cov="SiteCode", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Positive Impact", "ImmigratedIn", "ExtinctOut"), MainResults=TRUE, aspectrat=0.2)
BACIResultsPlot[[1]]
ggsave(paste0(FiguresFP, "Inkscape/BACIBars_BySite.png"), BACIResultsPlot[[1]], device="png", width = 350, height = 80, units = "mm", dpi=1000) #Save as a pdf)

RangeEst <- BACIResultsPlot[[2]]
RangeEst$Positive <- ifelse(RangeEst$variable=="ImmigratedIn" | RangeEst$variable=="ExtinctOut" | RangeEst$variable=="Positive Impact", "Yes", "No")
RangeEst <- dcast(as.data.table(RangeEst), SiteCode ~ Positive, sum, value.var="value")
RangeEst$PropPos <- RangeEst$Yes/rowSums(RangeEst[,c("Yes", "No")])
c(max(RangeEst$PropPos), min(RangeEst$PropPos), mean(RangeEst$PropPos), sd(RangeEst$PropPos))

### BA ###
BASummarise <- read.csv(file=paste0(DatabaseFP, "BA/ModelOutput/BACategories.csv")) #Read in BA data, subset to the same SiteSpec as in BACI (to enable accurate comparison)
BACISummarise <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
BACISummarise$SiteSpec <- paste0(BACISummarise$SiteCode, ".", BACISummarise$Species)
BASummarise <- BASummarise[BASummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, BASummarise$SiteSpec),]
BASummarise <- subset(BASummarise, OutcomeSummary!="Excluded")

FactorLevels <- c("Immigrated", "Positive Impact", "No Impact (population increasing)", "No Impact (no population trend)", "No Impact (population declining)", "Negative Impact", "Extinct") #Define the categories of BA
BAColourListPlot <- c("#1a275eff","#759bc7ff", GetMeanColour(c("#808d9cff", "#b9c1caff")), "#f9faf7ff", GetMeanColour(c("#ffeb88ff", "#ffcb48ff")), "#dc5851ff", "#a71111ff") #DEfine colours

#Get stats for paper
nrow(BASummarise) #npop 
length(unique(BASummarise$Species)) #nspec
length(unique(BASummarise$SiteCode)) #nsite

#By Species
BASummariseForPlot <- melt(dcast(as.data.table(subset(BASummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SiteSpec"), id.vars="Species")
BAResultsPlot <- BarPlotWidths(BASummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=BAColourListPlot, SuccessValues=c("Immigrated", "Positive Impact"), MainResults=TRUE, aspectrat=0.2)
BAResultsPlot[[1]]
ggsave(paste0(FiguresFP, "Inkscape/BABars_BySpecies.png"), BAResultsPlot[[1]], device="png", width = 350, height = 80, units = "mm", dpi=1000) #Save as a pdf)

#By Site
BASummariseForPlot <- melt(dcast(as.data.table(subset(BASummarise, OutcomeSummary!="Excluded")), SiteCode~Outcome, length, value.var="SiteSpec"), id.vars="SiteCode")
BAResultsPlot <- BarPlotWidths(BASummariseForPlot, Log=FALSE, Cov="SiteCode", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=BAColourListPlot, SuccessValues=c("Immigrated", "Positive Impact"), MainResults=TRUE, aspectrat=0.2)
BAResultsPlot[[1]]
ggsave(paste0(FiguresFP, "Inkscape/BABars_BySite.png"), BAResultsPlot[[1]], device="png", width = 350, height = 80, units = "mm", dpi=1000) #Save as a pdf)

### CI ###
CISummarise <- read.csv(file=paste0(DatabaseFP, "CI/ModelOutput/CICategories.csv")) #Read in CI data, subset to the same SiteSpec as in BACI (to enable accurate comparison)
BACISummarise <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
BACISummarise$SiteSpec <- paste0(BACISummarise$SiteCode, ".", BACISummarise$Species)
CISummarise <- CISummarise[CISummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, CISummarise$SiteSpec),]
CISummarise <- subset(CISummarise, OutcomeSummary!="Excluded")

FactorLevels <- c("Positive Impact", "No Impact (population increasing)", "No Impact (no population trend)", "No Impact (population declining)", "Negative Impact")
CIColourListPlot <- c("#3e6a9bff", GetMeanColour(c("#808d9cff", "#b9c1caff")), "#f9faf7ff", GetMeanColour(c("#ffeb88ff", "#ffcb48ff")), RedCols(5)[3])

#Get stats for paper
nrow(CISummarise) #npop
length(unique(CISummarise$Species)) #nspec
length(unique(CISummarise$SiteCode)) #nsite

#By Species
CISummariseForPlot <- melt(dcast(as.data.table(subset(CISummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SpecMatch"), id.vars="Species")
CIResultsPlot <- BarPlotWidths(CISummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=CIColourListPlot, SuccessValues=c("Positive Impact"), MainResults=TRUE, aspectrat=0.2)
CIResultsPlot[[1]]
ggsave(paste0(FiguresFP, "Inkscape/CIBars_BySpecies.png"), CIResultsPlot[[1]], device="png", width = 350, height = 80, units = "mm", dpi=1000) #Save as a pdf)

#By Site
CISummariseForPlot <- melt(dcast(as.data.table(subset(CISummarise, OutcomeSummary!="Excluded")), SiteCode~Outcome, length, value.var="SpecMatch"), id.vars="SiteCode")
CIResultsPlot <- BarPlotWidths(CISummariseForPlot, Log=FALSE, Cov="SiteCode", FactorLevels = FactorLevels, xlabel="SiteCode", AngledxLabs=FALSE, ColourList=CIColourListPlot, SuccessValues=c("Positive Impact"), MainResults=TRUE, aspectrat=0.2)
CIResultsPlot[[1]]
ggsave(paste0(FiguresFP, "Inkscape/CIBars_BySite.png"), CIResultsPlot[[1]], device="png", width = 350, height = 80, units = "mm", dpi=1000) #Save as a pdf)

#### Plot Predictor Models (Figure 4, Extended Data Figure 5, Figure S4) ####
### Get model output and organise data
#Read in model output from predictor models
AllModels <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "LHC/"), pattern = "CategoryModelsOutput.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  if(ncol(Dat)==1){
    return(NULL)
  }
  Dat$Scenario <- str_split_fixed(str_split_fixed(x, "[/]", 14)[,12], "[_]", 2)[,2]
  return(Dat)
}))
AllModels <- AllModels[!AllModels$Scenario %in% c("22", "23", "24", "25")]
AllModels <- subset(AllModels, Model != "RedList")
#Read in data to get some metadata (don't use this much anymore but here it is)
AllData <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "LHC/"), pattern = "CategoryModelsData.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  return(data.frame(Scenario=str_split_fixed(str_split_fixed(x, "[/]", 14)[,12], "[_]", 2)[,2], 
                    NumRowAll=nrow(Dat[complete.cases(Dat$PAAreaLog),]), NumRowNoPA = nrow(Dat),
                    NumSpecNoPA=length(unique(Dat$Species)), NumSiteNoPA=length(unique(Dat$SiteCode)),
                    NumSpec=length(unique(Dat[complete.cases(Dat$PAAreaLog),]$Species)), NumSite=length(unique(Dat[complete.cases(Dat$PAAreaLog),]$SiteCode))))
}), fill=TRUE)
AllData <- AllData[!AllData$Scenario %in% c("22", "23", "24", "25")]

AllModels <- merge(AllModels, AllData, by=c("Scenario")) #Bring together
AllModels <- subset(AllModels, Coef!= "Negative Impact|No Impact" & Coef != "No Impact|Positive Impact") #Remove the model output for the categories (we're not interested in this)
AllModels[AllModels$Scenario==21]$Scenario <- "Focal Analysis"

### Stacked Bars of BACI estimates for all LHCs
AllModelsMS <- subset(AllModels, Model=="All")
AllModelsMS <- AllModelsMS[complete.cases(AllModelsMS$Error),]
CoefNames <- read.csv(paste0(ResultsFP, "CoefNames.csv")) #This is just a csv that has full names for all the coefficients (e.g. "Protected area size" rather than "PA Area"), plus a marker for the order they should be in on the plot
AllModelsMS <- merge(AllModelsMS, CoefNames, by="Coef") #Add to model data
AllModelsMS <- AllModelsMS[order(AllModelsMS$OrderMarker),]

#This is a bunch of fairly inefficient cleaning of all the predictor covariates, how to order the factors etc, creating binary columns to mark whether estimates are significant or not, and positive or negative
AllModelsMS$CoefName <- factor(AllModelsMS$CoefName, levels=unique(AllModelsMS$CoefName))
AllModelsMS$CoefGroup <- factor(AllModelsMS$CoefGroup, levels=unique(AllModelsMS$CoefGroup))
AllModelsMS$Scenario <- factor(AllModelsMS$Scenario, levels=c(c(20:1), "Focal Analysis"))
AllModelsMS$PosNeg <- ifelse(AllModelsMS$Estimate>0, "Pos", "Neg")
AllModelsMS$Sig <- ifelse(AllModelsMS$P<0.05, "Sig", "Insig")
AllModelsMS$CoefName <- recode_factor(AllModelsMS$CoefName, "PA Area"="Protected area size", "Ramsar or SPA Site"="Managed for waterbirds")
AllModelsMS$ContVar <- "Categorical"
AllModelsMS[AllModelsMS$CoefName=="Body Size" | AllModelsMS$CoefName=="Protected area size" |AllModelsMS$CoefName=="Governance"]$ContVar <- "Continuous"

AllModelsMS$PosNegSig <- paste0(AllModelsMS$PosNeg, "_", AllModelsMS$Sig) #Create one value that expressed whether the value is positive or negative and significance or insignificant
AllModelsMS2 <- dcast(AllModelsMS, CoefName + CoefGroup ~ PosNegSig, length) #Count the number of LHC scenarios that have each kind of response (Pos/Neg & Significant/Insig) by covariate
AllModelsMS2 <- melt(AllModelsMS2, id.vars=c("CoefName", "CoefGroup"))
AllModelsMS2$variable <- factor(AllModelsMS2$variable, levels=c("Neg_Sig", "Neg_Insig", "Pos_Insig", "Pos_Sig")) #Define as factor
AllModelsMS2$variable <- recode_factor(AllModelsMS2$variable, "Neg_Sig" = "- (p < 0.05)", "Neg_Insig"="- (p > 0.05)","Pos_Insig" = "+ (p > 0.05)", "Pos_Sig"="+ (p < 0.05)") #Rename Factor levels for legend
OutcomeSummaryCols2 <- c("#9F1F1F", "goldenrod1", "#B3C0CC", "#223771")
names(OutcomeSummaryCols2) <- c("- (p < 0.05)", "- (p > 0.05)", "+ (p > 0.05)", "+ (p < 0.05)") #Name colours
AllModelsMS2$CoefGroup <- factor(AllModelsMS2$CoefGroup, levels=c("Area", "Managed", "Governance", "Body Size", "Migrant Status","Anthrome", "Order"))

#Plot
Plot <- ggplot(AllModelsMS2, aes(x=CoefName, y=value, fill=variable, colour=variable))+
  geom_bar(stat="identity", position="stack", size=0, colour="transparent", width=1)+
  facet_grid(CoefGroup~., scales="free", space='free')+
  scale_fill_manual(values = OutcomeSummaryCols2, name=("Estimate"))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  ylab("Number of analyses")+
  coord_flip()+
  theme(#aspect.ratio=1, 
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    text = element_text(size=6),
    panel.border = element_rect(size = 0.8, fill = NA, colour="grey60"),
    #panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(hjust = 0, size=6),
    strip.text.y = element_blank(),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black"),
    axis.title.x = element_text(size=6, colour="black"),
    axis.title.y = element_blank(),
    #axis.line = element_line(),
    legend.text = element_text(size=6, colour="black"),
    legend.title = element_text(size=6, colour="black"),
    legend.justification = "top",
    legend.key.size = unit(0.4, 'cm'))
Plot
ggsave(paste0(FiguresFP, "Inkscape/CoefficientsStacked.pdf"), Plot, device="pdf", width = 89, height = 80, units = "mm") #Save as a pdf , dpi=1000

### All LHC coefficient estimates with actual estimates and standard errors (ED Figure 5, Figure S3)

#This is a function to create the plot
AllCoefEstimates <- function(AllModels, NoPA=FALSE){
  #Reduce to NoPA model output if desired
  if(NoPA==FALSE){
    AllModels2 <- subset(AllModels, Model=="All")
  } else {
    AllModels2 <- subset(AllModels, Model=="NoPA")
    AllModels2$NumRowAll <- AllModels2$NumRowNoPA
    AllModels2$NumSpec <- AllModels2$NumSpecNoPA
    AllModels2$NumSite <- AllModels2$NumSiteNoPA
  }
  AllModels2 <- AllModels2[complete.cases(AllModels2$Error),] #Remove any errors

  #Organise data
  CoefNames <- read.csv(paste0(ResultsFP, "CoefNames.csv")) #This is just a csv that has full names for all the coefficients (e.g. "Protected area size" rather than "PA Area"), plus a marker for the order they should be in on the plot
  AllModels2 <- merge(AllModels2, CoefNames, by="Coef")
  AllModels2 <- AllModels2[order(AllModels2$OrderMarker),]
  AllModels2$CoefName <- factor(AllModels2$CoefName, levels=unique(AllModels2$CoefName))
  AllModels2$CoefGroup <- factor(AllModels2$CoefGroup, levels=unique(AllModels2$CoefGroup))
  AllModels2$Scenario <- factor(AllModels2$Scenario, levels=c(c(20:1), "Focal Analysis"))
  AllModels2$PosNeg <- ifelse(AllModels2$Estimate>0, "Pos", "Neg")
  
  #Plot
  ggplot(AllModels2, aes(x=Scenario, y=exp(Estimate), colour=CoefGroup, fill=CoefGroup))+
    geom_point(shape=1, alpha=0.7)+
    geom_errorbar(aes(ymin=exp(Estimate-(1.96*Error)), ymax=exp(Estimate+(1.96*Error))), width=.12, alpha=0.7)+
    geom_point(data=subset(AllModels2, P<0.05), aes(x=Scenario, y=exp(Estimate)))+
    geom_errorbar(data = subset(AllModels2, P<0.05), aes(ymin=exp(Estimate-(1.96*Error)), ymax=exp(Estimate+(1.96*Error))), width=.12)+
    facet_wrap(~CoefName, scales="free_x", nrow=3)+ #
    geom_hline(yintercept = 1.00, linetype="dashed")+
    scale_y_continuous(breaks = trans_breaks(identity, identity, n = 3))+
    scale_fill_manual(values=brewer.pal(length(levels(AllModels2$CoefGroup)), name="Dark2"), name="Covariate Group")+
    scale_colour_manual(values=brewer.pal(length(levels(AllModels2$CoefGroup)), name="Dark2"), name="Covariate Group")+
    xlab("Analysis")+
    ylab("Odds Ratio")+
    coord_flip()+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    theme(aspect.ratio=1.25, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize-2),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, size=plotfontsize-2),
          legend.title=element_text(size=plotfontsize),
          legend.key = element_blank(),
          axis.ticks = element_line(colour="black"),
          axis.text.x = element_text(size=plotfontsize-2, colour="black"),
          axis.text.y = element_text(size=plotfontsize-2, colour="black"),
          axis.title.y = element_text(size=plotfontsize, colour="black"),
          axis.title.x = element_text(size=plotfontsize, colour="black"),
          legend.text = element_text(size=plotfontsize-1, colour="black"),
          legend.justification = "left",
          legend.direction = "horizontal",
          legend.position = c(-0.10,-0.13),
          legend.spacing.x = unit(0.1, 'cm'),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.5, "cm"),
          plot.margin=unit(c(0.2,0.2,1.5,0.2),"cm"))
}

#With PA as a predictor
BACICoefEstimates <- AllCoefEstimates(AllModels)
ggsave(paste0(FiguresFP, "EDFig5_PredictorsAllEstimates.pdf"), BACICoefEstimates, device="pdf", width = 183, height = 147, units = "mm") #Save as a pdf

#Without PA as a predictor
BACICoefEstimatesNoPA <- AllCoefEstimates(AllModels, NoPA=TRUE)
ggsave(paste0(FiguresFP, "FigS4_PredictorsAllEstimateNoPA.pdf"), BACICoefEstimatesNoPA, device="pdf", width = 183, height = 147, units = "mm") #Save as a pdf

#### Extended Data Figure 2 ####
AllBACICategories <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "LHC/"), pattern = "Categories.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  Dat$BABACI <- str_split_fixed(x, "[/]", 14)[,13]
  Dat$Scenario <- str_split_fixed(str_split_fixed(x, "[/]", 14)[,12], "[_]", 2)[,2]
  return(Dat)
}), fill=TRUE)

AllBACICategories <- AllBACICategories[!AllBACICategories$Scenario %in% c("22", "23", "24", "25")]

#Reduce BA and CI to only sitespec occurring in BACI
AllBACICategories[AllBACICategories$BABACI=="BACI"]$SiteSpec <- paste0(AllBACICategories[AllBACICategories$BABACI=="BACI"]$SiteCode, ".", AllBACICategories[AllBACICategories$BABACI=="BACI"]$Species)
AllBACICategories <- rbindlist(lapply(unique(AllBACICategories$Scenario), function(Scenny){
  ScenDat <- subset(AllBACICategories, Scenario==Scenny)
  ScenBACI <- subset(ScenDat, BABACI=="BACI")
  ScenBA <- subset(ScenDat, BABACI=="BA")
  ScenCI <- subset(ScenDat, BABACI=="CI")
  ScenBA <- ScenBA[ScenBA$SiteSpec %in% intersect(ScenBACI$SiteSpec, ScenBA$SiteSpec),]
  ScenCI <- ScenCI[ScenCI$SiteSpec %in% intersect(ScenBACI$SiteSpec, ScenCI$SiteSpec),]
  return(rbind(ScenBACI, ScenBA, ScenCI))
}))

AllBACICategoriesProp <- dcast(AllBACICategories, Scenario + BABACI ~ OutcomeSummary, length, value.var="SiteCode")
AllBACICategoriesProp <- cbind(AllBACICategoriesProp[,c(1:2)], prop.table(as.matrix(AllBACICategoriesProp[,-c(1:2)]), margin = 1))
AllBACICategoriesProp <- melt(as.data.table(AllBACICategoriesProp), id.vars=c("Scenario", "BABACI"))
PlotCols <- c("Negative Impact"= NegCol, "No Impact" = "#ffbf27ff", "Positive Impact" = PosCol, "Excluded" = "grey60")

AllBACICategoriesProp$variable <- factor(AllBACICategoriesProp$variable, levels=rev(levels(AllBACICategoriesProp$variable)))
AllBACICategoriesProp$BABACI <- factor(AllBACICategoriesProp$BABACI, levels=c("BACI", "BA", "CI"))
#AllBACICategoriesProp$BABACI <- recode_factor(AllBACICategoriesProp$BABACI, "BACI" = "bold('a ') (BACI)", "BA" = "bold('b ') (BA)", "CI" = "bold('c ') (CI)")

OutcomePlot <- function(Grouping, YaxisLabel=FALSE){
  if(YaxisLabel==TRUE){
    YaxLab <- element_text(size=plotfontsize) #family="Helvetica Neue", 
    Ytext <- element_text(size=plotfontsize-1)
  } else {
    YaxLab <- element_blank()
    Ytext <- element_blank()
  }
  ggplot()+
    geom_boxplot(data=subset(AllBACICategoriesProp, BABACI==Grouping), aes(x=variable, y=as.numeric(value), colour=variable, fill=variable), size=0.5, outlier.shape=NA, alpha=0.2)+
    geom_jitter(data=subset(AllBACICategoriesProp, Scenario!=21 & BABACI==Grouping), aes(x=variable, y=as.numeric(value), colour=variable), size=1, alpha=0.5)+
    geom_point(data=subset(AllBACICategoriesProp, Scenario==21 & BABACI==Grouping), aes(x=variable, y=as.numeric(value), colour=variable, group=variable), size=3)+
    #facet_wrap(~BABACI, ncol=3, labeller=label_parsed)+
    scale_colour_manual(values=PlotCols, name="Outcome")+
    scale_fill_manual(values=PlotCols, name="Outcome")+
    ylab("Proportion of populations")+
    xlab(Grouping)+
    scale_y_continuous(limits=c(0, 0.805), breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
    theme(aspect.ratio=1.1, 
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, size=8),
          legend.title=element_text(size=plotfontsize),
          axis.ticks = element_line(colour="black"),
          axis.text.x = element_text(size=6),
          axis.title.x=element_text(size=plotfontsize),
          axis.title.y = YaxLab,
          axis.text.y = Ytext,
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top",
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
}
BACIPlot <- OutcomePlot("BACI", YaxisLabel = TRUE)
BAPlot <- OutcomePlot("BA")
CIPlot <- OutcomePlot("CI")
AlPlot <- plot_grid(BACIPlot, BAPlot, CIPlot, ncol=3, labels=c("a", "b", "c"), label_size = 8, align="hv", label_x = 0.02, label_y = 1.01)
ggsave(paste0(FiguresFP, "EDFig2_AllOutcomes.pdf"), AlPlot, device="pdf", width = 183, height = 65, units = "mm") #Save as a pdf

#### Stats reported in paper ####
### Get the proportion of populations that fall into the various impact categories for BACI
BACISummarise <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))

BACISummarise$OutcomeSummary <- as.character(BACISummarise$OutcomeSummary)
BACISummarise$OutcomeSummary2 <- as.character(BACISummarise$OutcomeSummary)

BACISummarise$OutcomeSummary2 <- ifelse(BACISummarise$Outcome=="No Impact (population increasing)" | BACISummarise$Outcome == "ImmigratedBoth", "No Impact (Pos)", BACISummarise$OutcomeSummary2)
BACISummarise$OutcomeSummary2 <- ifelse(BACISummarise$Outcome=="No Impact (population declining)" | BACISummarise$Outcome == "ExtinctBoth", "No Impact (Neg)", BACISummarise$OutcomeSummary2)
BACISummarise$OutcomeSummary2 <- ifelse(BACISummarise$Outcome=="No Impact (no population trend)", "No Impact (no population trend)", BACISummarise$OutcomeSummary2)

prop.table(table(BACISummarise$OutcomeSummary2)) #Proportion across all categories
prop.table(table(BACISummarise$OutcomeSummary)) #Proportion across summary categories (Positive, No Impact, Negative)
prop.table(table(subset(BACISummarise, OutcomeSummary=="No Impact")$OutcomeSummary2)) #Propotion just within No impact

### N for Latin HyperCubes
#Get number of sites/species/populations in each latin hypercube sample
AllN <- rbindlist(pblapply(c(1:21), function(Scen){
  load(paste0(ResultsFP, "LHC/Scenario_", Scen, "/BACI/BACI.RData"))
  return(data.frame(Scenario=Scen, NSites=length(unique(subset(BACI, CI==1)$SiteCode)), NSpec=length(unique(BACI$Species)), NPop=length(unique(BACI$SpecMatch))))
}))

AllN <- AllN[c(21, 1:20),] #Put focal analysis at the top
write.csv(AllN, paste0(FiguresFP, "AllN.csv"), row.names=FALSE) #Save

#Get number of sites/species/populations across all analyses
AllDat <- rbindlist(pblapply(c(1:21), function(Scen){
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
  load(paste0(DatabaseFP, "BACI/BACI.RData"))
  BACI$Scenario <- Scen
  return(BACI)
}), fill=TRUE)

data.frame(NProtectedSites=length(unique(subset(AllDat, CI==1)$SiteCode)), NUnprotectedSites=length(unique(subset(AllDat, CI==0)$SiteCode)), NSpec=length(unique(AllDat$Species)), NPop=length(unique(AllDat$SpecMatch)), NCountries=length(unique(AllDat$ISO3)), NContinent = length(unique(AllDat$GeoRegion)))

### Proportion of sites in Managed vs Unmanaged PAs
load(file=paste0(DatabaseFP, "BACI/BACIPredictors.RData"))
load(file=paste0(DatabaseFP, "BA/BAPredictors.RData"))
load(file=paste0(DatabaseFP, "CI/CIPredictors.RData"))

prop.table(table(BACIPredictors$Managed))
prop.table(table(BAPredictors$Managed))
prop.table(table(CIPredictors$Managed))

#How many "negative impact" populations are in unmanaged sites? 
ManageProps <- rbindlist(pblapply(c(1:21), function(Scenny){
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scenny, "/")
  
  BACICategories <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
  load(file=paste0(DatabaseFP, "BACI/BACIPredictors.RData"))
  
  BACIModelData <- PrepCategoryModelData(BACICategories, BACIPredictors, paste0(DatabaseFP, "BACI/"), SaveOutput = FALSE)
  BACIModelData$SiteSpec <- paste0(BACIModelData$SiteCode, ".", BACIModelData$Species)
  
  BACIModelDataNeg <- subset(BACIModelData, OutcomeSummary=="Negative Impact" | Outcome=="No Impact (population declining)" | Outcome=="ExtinctBoth")
  Props <- as.data.table(prop.table(table(BACIModelDataNeg$Managed)))
  Props$Scen <- Scenny
  Props$PropAllMan <- nrow(subset(BACIModelDataNeg, Managed=="Managed"))/nrow(BACIModelData)
  Props$PropAllUn <- nrow(subset(BACIModelDataNeg, Managed=="Other"))/nrow(BACIModelData)
  return(Props)
}))

ManageProps <- dcast(ManageProps, Scen + PropAllUn + PropAllMan ~ V1, value.var="N")
ManageProps$PropUnmanaged <- ManageProps$Other/(ManageProps$Managed + ManageProps$Other)
hist(ManageProps$PropUnmanaged)

### Number of species and sites that we started out with pre any data subsetting
load(file=paste0(DataFP, "WaterbirdData_2020/WaterbirdCounts_AllCovariates.RData"))
length(unique(WaterbirdCounts$SiteCode))
length(unique(WaterbirdCounts$Species))
max(CountData$Year)
min(CountData$Year)

#### Taxonomic breakdown of datasets (Table S1) ####
#Get the count data for the three groups
Scen <- 21
DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")

load(paste0(DatabaseFP, "BACI/BACI.RData"))
load(paste0(DatabaseFP, "BA/BA.RData"))
load(paste0(DatabaseFP, "CI/CI.RData"))

BA <- BA[BA$SiteSpec %in% intersect(BA$SiteSpec, BACI$SiteSpec),]
CI <- CI[CI$SiteSpec %in% intersect(CI$SiteSpec, BACI$SiteSpec),]

#Function to extract taxonomic data
TaxaData <- function(Dataset, Group){
  AllSitesStats <- as.data.table(unique(Dataset[,c("Species","Genus", "Family", "Order", "SiteCode", "CI")]))
  
  #Sites per Order
  AllSitesStatsSite <- unique(AllSitesStats[,c("Order", "SiteCode", "CI")])
  Sites <- dcast(AllSitesStatsSite, Order ~CI, length, value.var="SiteCode")
  names(Sites)[names(Sites)==1] <- paste0(Group, " Protected Sites")
  names(Sites)[names(Sites)==0] <- paste0(Group, " Unprotected Sites")
  
  #Species per Order
  AllSitesStatsSpec <- unique(AllSitesStats[,c("Order", "Genus", "Species", "CI")])
  Species <- dcast(AllSitesStatsSpec, Order ~., length, value.var="Species")
  names(Species) <- c("Order", paste0(Group, " Species"))
  
  #Genera per Order
  AllSitesStatsGen <- unique(AllSitesStats[,c("Order", "Genus", "CI")])
  Genera <- dcast(AllSitesStatsGen, Order ~., length, value.var="Genus")
  names(Genera) <- c("Order", paste0(Group, " Genera"))
  
  #Families per Order
  AllSitesStatsFam <- unique(AllSitesStats[,c("Order", "Family", "CI")])
  Families <- dcast(AllSitesStatsFam, Order ~., length, value.var="Family")
  names(Families) <- c("Order", paste0(Group, " Family"))
  
  OrderData <- list(Sites, Species, Genera, Families) %>% reduce(left_join, by = c("Order"))
  
  NOrder <- length(unique(AllSitesStats$Order))
  NFamily <- length(unique(AllSitesStats$Family))
  NGenus <- length(unique(AllSitesStats$Genus))
  NSpecies <- length(unique(AllSitesStats$Species))
  TaxaList <- list(OrderData, NOrder,NFamily,NGenus, NSpecies)
  names(TaxaList) <- c("OrderData", "NumOrders", "NumFamilies", "NumGenera", "NumSpecies")
  return(TaxaList)
}

#Extract for each group
BATaxa <- TaxaData(BA, "BA")
CITaxa <- TaxaData(CI, "CI")
BACITaxa <- TaxaData(BACI, "BACI")

AllTaxa <- list(BATaxa[[1]], CITaxa[[1]], BACITaxa[[1]]) %>% reduce(left_join, by = c("Order")) #Bring together
AllTaxa[is.na(AllTaxa)] <- 0
AllTaxa$'BA Unprotected Sites' <- NA
AllTaxa <- AllTaxa[,c(1,5,10,15,4,9,14,3,8,13,2,7,12,16,6,11)]

write.csv(AllTaxa, paste0(FiguresFP, "TaxaMatchStats.csv"), row.names=FALSE) #Save

#### Protected area types (Table S2) ####
#Count how many sites are in protected areas of each type, as per Table S2
load(file=paste0(DatabaseFP, "BACI/BACI.RData"))

ProtectedSites <- fread(paste0(ResultsFP, "ProtectedSites_Uncleaned.csv")) #Read in protected area data
ManagedSitesTerms <- c("ramsar", "Birds Directive") #Define the terms we use to classify managed sites (same term as in the "GetPredictors" function)

SitesToRemove <- unique(subset(ProtectedSites, DESIG=="UNESCO-MAB Biosphere Reserve" | DESIG_ENG=="UNESCO-MAB Biosphere Reserve" | STATUS=="Proposed" | STATUS=="Not Reported" | STATUS_YR==0)$SiteCode) #Remove these sites (done earlier for main analysis)
ProtectedSites <- ProtectedSites[!ProtectedSites$SiteCode %in% SitesToRemove,]
ProtectedSites <- ProtectedSites[ProtectedSites$SiteCode %in% subset(BACI, CI==1)$SiteCode,]

ManagedSites <- ProtectedSites[grepl(paste(ManagedSitesTerms, collapse="|"), ProtectedSites$DESIG_ENG, ignore.case=TRUE)] #Extract the managed sites
UnmanagedSites <- ProtectedSites[!grepl(paste(ManagedSitesTerms, collapse="|"), ProtectedSites$DESIG_ENG, ignore.case=TRUE)] #And then the unmanaged sites
UnmanagedSites <- UnmanagedSites[!UnmanagedSites$SiteCode %in% ManagedSites$SiteCode,]

#Reorganise data
ManagedSites <- unique(ManagedSites[,c("DESIG_ENG", "SiteCode")])
UnmanagedSites <- unique(UnmanagedSites[,c("DESIG_ENG", "SiteCode")])
ManagedSites[ManagedSites$DESIG_ENG=="RAMSAR Sites"]$DESIG_ENG <- "Ramsar Site, Wetland of International Importance"

#Count the number of sites in protected areas of each type
ManagedSitesCast <- dcast(ManagedSites, DESIG_ENG~., length)
ManagedSitesCast <- ManagedSitesCast[order(ManagedSitesCast$., decreasing=TRUE),]

UnmanagedSitesCast <- dcast(UnmanagedSites, DESIG_ENG~., length)
UnmanagedSitesCast <- UnmanagedSitesCast[order(UnmanagedSitesCast$., decreasing=TRUE),]

#Bring together, save
AllManage <- rbind(ManagedSitesCast, UnmanagedSitesCast)
write.csv(AllManage, paste0(FiguresFP, "NSitesByManagement.csv"), row.names=FALSE)

#### Regression to mean (Figure S3) ####
#Calculate the extend to which the change in trend is predicted by the before trend for both BA and BACI (in the populations where we can model change from before to after, cf populations with All Zeroes)
#Then plot as per Figure S2
#To do this we need the model output from the change in trend models, and we also need to calculate what the slope is in the before period
Scen <- 21
DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_21/")

#Load data
BACILevelSlope <- as.data.table(read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACILevelSlope.csv"))) #Read in the change in slope models for BACI
load(file=paste0(DatabaseFP, "BACI/BACI.RData")) #Read in the actual count data for BACI
BALevelSlope <- as.data.table(read.csv(file=paste0(DatabaseFP, "BA/ModelOutput/BALevelSlope.csv"))) #Read in the change in slope models for BA
load(file=paste0(DatabaseFP, "BA/BA.RData")) #Read in the actual count data for BA

#BA
BA <- BA[BA$SiteSpec %in% intersect(BACI$SiteSpec, BA$SiteSpec),]
BALevelSlope <- BALevelSlope[BALevelSlope$SiteSpec %in% BA$SiteSpec,]
BABeforeSlope <- CalculateBeforeSlopes(subset(BA, BA==0), Scenarios[Scen,]$ZeroThresh, parallelise = TRUE, actualvalues = TRUE, Scenarios[Scen,]$PThresh, Scenarios[Scen,]$TotalYearsBuffer) #Calcuate what the trend is in the before period for BA
BABeforeSlope$BeforeSlope <- NULL
names(BABeforeSlope) <- c("SiteSpec", "CI", "BeforeSlope", "BeforeP") 
BALevelSlope <- merge(BALevelSlope, BABeforeSlope, by=c("SiteSpec")) #Add before slope data to slope model data
BALevelSlope <- subset(BALevelSlope, Coef=="Year:BA")
BALevelSlope <- BALevelSlope[,c("SiteSpec", "Estimate", "P", "BeforeSlope", "BeforeP")]
BALevelSlope$BABACI <- "BA"

#BACI get data
BACINonZero <- subset(subset(BACI, BA==0 & CI==1), ModelCat=="Model")
BACIBeforeSlope <- CalculateBeforeSlopes(BACINonZero, Scenarios[Scen,]$ZeroThresh, parallelise = TRUE, actualvalues = TRUE, Scenarios[Scen,]$PThresh, Scenarios[Scen,]$TotalYearsBuffer) #Calcuate what the trend is in the before period for BACI
BACIBeforeSlope$BeforeSlope <- NULL
names(BACIBeforeSlope) <- c("SiteSpec", "CI", "BeforeSlope", "BeforeP")
BACI$CI <- as.character(BACI$CI)
BACINonZero$CI <- as.character(BACINonZero$CI)
BACIBeforeSlope <- merge(BACIBeforeSlope, unique(BACINonZero[,c("SiteSpec", "SpecMatch", "CI")]), by=c("CI", "SiteSpec"))
BACILevelSlope <- merge(BACILevelSlope, BACIBeforeSlope, by="SpecMatch") #Add before slope data to slope model data
BACILevelSlope <- subset(BACILevelSlope, Coef=="Year:BA:CI")
BACILevelSlope <- BACILevelSlope[,c("SiteSpec", "Estimate", "P", "BeforeSlope", "BeforeP")]
BACILevelSlope$BABACI <- "BACI"

#Bring BA and BACI together (just makes it easy for plotting/calculations)
LevelSlope <- rbind(BALevelSlope, BACILevelSlope)

#Classify each population based on whether the before slope is significant or not, and whether the change in slope is significant or not
LevelSlope$Significance <- ifelse(LevelSlope$P<0.05, ifelse(LevelSlope$BeforeP<0.05, "Both Significant", "Change Significant"), ifelse(LevelSlope$BeforeP<0.05, "BeforeSlope Significant", "Neither Significant"))
LevelSlope$Significance <- factor(LevelSlope$Significance, levels=c("Both Significant", "BeforeSlope Significant", "Change Significant", "Neither Significant"))
LevelSlope <- LevelSlope[order(LevelSlope$Significance, decreasing=TRUE),]
LevelSlope <- LevelSlope[complete.cases(LevelSlope),]

#Subset to only cases where the before slope and the change is significant (we're not interested in the other cases if there hasn't been significant change)
BothSig <- subset(LevelSlope, Significance=="Change Significant" | Significance=="Both Significant")

summary(lm(Estimate~BeforeSlope, data=subset(LevelSlope, BABACI=="BA"))) #Get the r2 for BA
summary(lm(Estimate~BeforeSlope, data=subset(LevelSlope, BABACI=="BACI"))) #Get the r2 for BACI

Rsquares <- data.frame(BABACI=c("BA", "BACI"), 
                       rsq=c(summary(lm(Estimate~BeforeSlope, data=subset(BothSig, BABACI=="BA")))[[9]], summary(lm(Estimate~BeforeSlope, data=subset(BothSig, BABACI=="BACI")))[[9]]),
                       BeforeSlope = c(3, 3), Estimate = c(2.2, 3.3)) #Get the rsquares and define where they should occur on the plot
Rsquares$rsq <- round_any(Rsquares$rsq, 0.01)
Rsquares$rsq <- paste("R^2 == ", Rsquares$rsq)

BothSig$Significance <- recode_factor(BothSig$Significance, "Both Significant" = "Significant before slope", "Change Significant" = "Insignificant before slope") #Code factors

#Plot
RegressToMean <- ggplot(data=BothSig, aes(x=BeforeSlope, y=Estimate))+ #, fill=Result
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_point(shape=16, aes(colour=Significance))+
  geom_smooth(method = "lm", colour=NegCol, size=0.5, fill=NegCol)+
  scale_colour_manual(values=c("black", "grey"), name="")+
  facet_wrap(~BABACI, scales="free")+
  xlab("Slope in years before designation")+
  ylab("Change in slope after designation")+ 
  geom_text(data = Rsquares,aes(label=rsq), parse=TRUE, size=2.5)+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text = element_text(size=plotfontsize, colour="black"),
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "left",
        legend.key=element_blank(),
        legend.position = "bottom")
RegressToMean
ggsave(paste0(FiguresFP, "FigS3_RegressToMean.pdf"), RegressToMean, device="pdf", width = 183, height = 110, units = "mm") #Save as a pdf)

#Colour BA by total slope estimate
BAColour <- FALSE
if(BAColour==TRUE){
  BATotalSlope <- CalculateBeforeSlopes(BA, Scenarios[Scen,]$ZeroThresh, parallelise = TRUE, actualvalues = TRUE, PThresh, TotalYearsBuffer)
  BATotalSlope <- BATotalSlope[,c("SiteSpec", "SlopeEstimate", "SlopeSignificance")]
  names(BATotalSlope) <- c("SiteSpec", "TotalSlope", "TotalSignificance")
  BALevelSlope <- merge(BALevelSlope, BATotalSlope, by="SiteSpec")
  BALevelSlope$TotalSignal <- ifelse(BALevelSlope$TotalSignificance>0.05, "Insig", ifelse(BALevelSlope$TotalSlope<0, "Neg", "Pos"))
  
  RegressToMeanBAColour <- ggplot(data=BALevelSlope, aes(x=BeforeSlope, y=Estimate))+ #, fill=Result
    geom_vline(xintercept = 0, linetype="dashed")+
    geom_hline(yintercept = 0, linetype="dashed")+
    geom_point(alpha=0.3, aes(color=TotalSignal))+
    geom_smooth(method = "lm", colour=NegCol, size=0.5, fill=NegCol)+
    facet_wrap(~BABACI, scales="free")+
    xlab("Slope in years before designation")+
    ylab("Change in slope after designation")+ 
    #xlim(-20,20)+
    #ylim(-20,20)+
    #geom_text(data = Rsquares,aes(label=rsq), parse=TRUE)+
    theme(aspect.ratio=1, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, size=plotfontsize),
          legend.title=element_text(size=plotfontsize),
          axis.text = element_text(size=plotfontsize, colour="black"),
          axis.ticks = element_line(colour="black"),
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top")
  RegressToMeanBAColour
}

#### Abundance estimates of analysis species vs. all waterbirds (Figure S6) ####
#Get species in each LHC run
AnalysisTax <- rbindlist(lapply(c(1:21), function(Scenny){
  InputFP <- paste0(ResultsFP, "LHC/Scenario_", Scenny, "/BACI/")
  BACICategories <- read.csv(file=paste0(InputFP, "ModelOutput/BACICategories.csv"))
  return(data.frame("Species" = unique(BACICategories$Species), "Group" = Scenny))
}))

#Get abundance data, merge
AbundanceData <- read.csv(paste0(DataFP, "BirdAbundance/BirdAbundance.csv"))
AnalysisTax <- merge(AnalysisTax, AbundanceData[,c("Scientific.name", "Abundance.estimate")], by.x="Species", by.y="Scientific.name", all.x=TRUE)
AnalysisTax$Group <- as.character(AnalysisTax$Group)
AnalysisTax[AnalysisTax$Group==21]$Group <- "Focal Analysis"

#Get abundance data on all waterbirdspecies, add to dataset
WaterbirdFamilies <- read.csv(paste0(DataFP, "WaterbirdData_2020/WaterbirdFamilies.csv"))
AbundanceData$Family <- str_split_fixed(AbundanceData$Family, "[ (]", 2)[,1]
AbundanceDataWater <- AbundanceData[AbundanceData$Family %in% WaterbirdFamilies$WaterbirdFamilies,]
AbundanceDataWater$Group <- "All Waterbirds"
AnalysisTax <- rbind(AnalysisTax, AbundanceDataWater[,c("Scientific.name", "Group", "Abundance.estimate")], use.names=FALSE)

#Create a dataframe to specify the width of bars (so that focal analysis and all waterbirds are a bit separate from the LHC runs)
WidthData <- data.frame("Width" = c(rep(1, 20),3, 3), "Position" = c(1:21, 23), "Names"=c(1:22))
WidthData$Names2 <- c(c(1:20), "Focal Analysis", "All Waterbirds")
WidthData$Position2 <- 0.5 * (cumsum(WidthData$Width) + cumsum(c(0, WidthData$Width[-length(WidthData$Width)])))

AnalysisTax <- merge(AnalysisTax, WidthData[,c("Names2", "Position2", "Width")], by.x="Group", by.y="Names2")
AnalysisTax$Names2 <- factor(AnalysisTax$Group, levels=c(c(1:20), "Focal Analysis", "All Waterbirds"))
AnalysisTax <- AnalysisTax[order(AnalysisTax$Names2),]
AnalysisTax$Names3 <- str_wrap(AnalysisTax$Names2, width=12)
AnalysisTax$Names3 <- factor(AnalysisTax$Names3, levels=unique(AnalysisTax$Names3))

#Plot
AbundancePlot <-ggplot(AnalysisTax)+
  geom_boxplot(aes(x=Position2, y=Abundance.estimate, group=Names3), width=0.8)+ # x=Names3, weight=weight,  , varwidth = TRUE
  scale_x_continuous(expand = c(0.01, 0.01), breaks=WidthData$Position2, labels=WidthData$Names2)+
  scale_y_log10(expand = c(0,0), breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), 
                labels=c("1", "10", "100", "1000", "10,000", "100,000", "1,000,000", "10,000,000", "100,000,000", "1,000,000,000"))+
  ylab("Abundance Estimate")+
  theme(aspect.ratio=0.5, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.ticks = element_line(colour="black"),
        axis.title.x = element_blank(),
        axis.text = element_text(size=plotfontsize),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")

ggsave(paste0(FiguresFP, "FigS6_AbundanceComparisons.pdf"), AbundancePlot, device="pdf", width = 183, height = 85, units = "mm") #Save as a pdf)

#How many species without abundance estimates
AnalysisTax$NAFlag <- ifelse(is.na(AnalysisTax$Abundance.estimate), 1, 0)
PropNas <- dcast(AnalysisTax, Group ~ NAFlag, length)
PropNas$Prop <- PropNas$'1'/sum(PropNas[,c(2,3)])*100

#### Analysis of differences in IUCN status for European birds (Figure S7) ####

#We run the models to asses EU red list when we run all cumulative link models. This code now pulls out the results and plots. 
AllModels <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "LHC/"), pattern = "CategoryModelsOutput.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  if(ncol(Dat)==1){
    return(NULL)
  }
  Dat$Scenario <- str_split_fixed(str_split_fixed(x, "[/]", 14)[,12], "[_]", 2)[,2]
  return(Dat)
}))
AllModels <- AllModels[!AllModels$Scenario %in% c("22", "23", "24", "25")]
AllModels <- subset(AllModels, Model == "RedList")
AllModels[AllModels$Scenario==21]$Scenario <- "Focal Analysis"
AllModels$Scenario <- factor(AllModels$Scenario, levels=c(c(20:1), "Focal Analysis"))
AllModels$PosNeg <- ifelse(AllModels$Estimate>0, "Pos", "Neg")
AllModels <- subset(AllModels, Coef=="RedList_eu27Threatened")
AllModels <- subset(AllModels, Scenario != "20")
RedList <- ggplot(AllModels, aes(x=Scenario, y=Estimate))+
  geom_point(shape=1, alpha=0.7)+
  geom_errorbar(aes(ymin=Estimate-(1.96*Error), ymax=Estimate+(1.96*Error)), width=.12, alpha=0.7)+
  geom_point(data=subset(AllModels, P<0.05), aes(x=Scenario, y=Estimate))+
  geom_errorbar(data = subset(AllModels, P<0.05), aes(ymin=Estimate-(1.96*Error), ymax=Estimate+(1.96*Error)), width=.12)+
  geom_hline(yintercept = 0.00, linetype="dashed")+
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 3), limits = c(-1.8, 1.8))+
  coord_flip()+
  theme(aspect.ratio=1.55, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=10),
        legend.title=element_text(size=plotfontsize),
        axis.ticks = element_line(colour="black"),
        axis.text.x = element_text(size=6, colour="black"),
        axis.title.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=6, colour="black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
RedList
ggsave(paste0(FiguresFP, "FigS7_RedListSensitivities.pdf"), RedList, device="pdf", width = 65, height = 80, units = "mm") #Save as a pdf
ggsave(paste0(FiguresFP, "FigS7_RedListSensitivities.png"), RedList, device="png", width = 65, height = 80, units = "mm", dpi=1000) #Save as a pdf

#Find out how many threatened/non threatened species are in Scenario 20 (given its estimates are insane)
Scen20Data <- read.csv(list.files(paste0(ResultsFP, "LHC/Scenario_20/"), pattern = "CategoryModelsData.csv$", recursive = TRUE, full.names = TRUE))
EURL <- read.csv(paste0(DataFP, "EURedList/EuropeanRedList.csv"), header=TRUE)
EuMembers <- read.csv(paste0(DataFP, "EURedList/EUMemberStates.csv"), header=FALSE)
Scen20Data <- Scen20Data[Scen20Data$Country %in% EuMembers$V1,]
Scen20Data <- merge(Scen20Data, EURL, by="Species")
Scen20Data <- subset(Scen20Data, RedList_eu27!="NE")
Scen20Data$RedList_eu27 <- factor(Scen20Data$RedList_eu27, levels=c("LC", "NT", "VU", "EN"))
Scen20Data$RedList_eu27 <- recode_factor(Scen20Data$RedList_eu27, "LC" = "Least Concern", "NT" = "Threatened", "VU" = "Threatened", "EN" = "Threatened")
table(unique(Scen20Data[,c("Species", "RedList_eu27")])$RedList_eu27)

#### Sensitivity for length of time sampled and strictness for identifying a significant change in trend (Figure S8, 9)####
### Proportion of populations in positive, no and negative impact categories for BA, CI and BACI

#These two additional analyses are run with all the main analyses as scenarios 22, 23, 24 and 25. 
#Now we read in the data and plot, following the same format as in other code sections

AllBACICategories <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "LHC/"), pattern = "Categories.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  Dat$BABACI <- str_split_fixed(x, "[/]", 14)[,13]
  Dat$Scenario <- str_split_fixed(str_split_fixed(x, "[/]", 14)[,12], "[_]", 2)[,2]
  return(Dat)
}), fill=TRUE)

#Reduce BA and CI to only sitespec occurring in BACI
AllBACICategories[AllBACICategories$BABACI=="BACI"]$SiteSpec <- paste0(AllBACICategories[AllBACICategories$BABACI=="BACI"]$SiteCode, ".", AllBACICategories[AllBACICategories$BABACI=="BACI"]$Species)
AllBACICategories <- rbindlist(lapply(unique(AllBACICategories$Scenario), function(Scenny){
  ScenDat <- subset(AllBACICategories, Scenario==Scenny)
  ScenBACI <- subset(ScenDat, BABACI=="BACI")
  ScenBA <- subset(ScenDat, BABACI=="BA")
  ScenCI <- subset(ScenDat, BABACI=="CI")
  ScenBA <- ScenBA[ScenBA$SiteSpec %in% intersect(ScenBACI$SiteSpec, ScenBA$SiteSpec),]
  ScenCI <- ScenCI[ScenCI$SiteSpec %in% intersect(ScenBACI$SiteSpec, ScenCI$SiteSpec),]
  return(rbind(ScenBACI, ScenBA, ScenCI))
}))

AllBACICategoriesProp <- dcast(AllBACICategories, Scenario + BABACI ~ OutcomeSummary, length, value.var="SiteCode")
AllBACICategoriesProp <- cbind(AllBACICategoriesProp[,c(1:2)], prop.table(as.matrix(AllBACICategoriesProp[,-c(1:2)]), margin = 1))
AllBACICategoriesProp <- melt(as.data.table(AllBACICategoriesProp), id.vars=c("Scenario", "BABACI"))
OutcomeSummaryCols <- c("Negative Impact"= NegCol, "No Impact" = "#ffbf27ff", "Positive Impact" = PosCol, "Excluded" = "grey80")
OutcomeSummaryCols2 <- c("Negative Impact"= "blue", "No Impact" = "blue", "Positive Impact" = "blue", "Excluded" = "blue")

AllBACICategoriesProp$variable <- factor(AllBACICategoriesProp$variable, levels=rev(levels(AllBACICategoriesProp$variable)))
AllBACICategoriesProp$BABACI <- factor(AllBACICategoriesProp$BABACI, levels=c("BACI", "BA", "CI"))
AllBACICategoriesProp$Scenario <- as.character(AllBACICategoriesProp$Scenario)

AllBACICategoriesProp$ScenarioName <- "Nothing"
AllBACICategoriesProp[AllBACICategoriesProp$Scenario == "22"]$ScenarioName <- "11 to 15 years sampled before and after protection"
AllBACICategoriesProp[AllBACICategoriesProp$Scenario == "23"]$ScenarioName <- "16 to 20 years sampled before and after protection"
AllBACICategoriesProp[AllBACICategoriesProp$Scenario == "24"]$ScenarioName <- "p < 0.1 to detect impact of protection"
AllBACICategoriesProp[AllBACICategoriesProp$Scenario == "25"]$ScenarioName <- "p < 0.2 to detect impact of protection"
AllBACICategoriesProp$ScenarioName <- factor(AllBACICategoriesProp$ScenarioName, levels=unique(AllBACICategoriesProp$ScenarioName))

OutcomePlotSensitivities <- function(Grouping, YaxisLabel=FALSE){
  if(YaxisLabel==TRUE){
    YaxLab <- element_text(size=plotfontsize) #family="Helvetica Neue", 
    Ytext <- element_text(size=plotfontsize-1)
  } else {
    YaxLab <- element_blank()
    Ytext <- element_blank()
  }
  ggplot()+
    geom_boxplot(data=subset(AllBACICategoriesProp[!AllBACICategoriesProp$Scenario %in% c("22", "23", "24", "25")], BABACI==Grouping), aes(x=variable, y=as.numeric(value)), size=0.5, outlier.shape=NA, colour="grey80")+
    geom_beeswarm(data=subset(AllBACICategoriesProp[AllBACICategoriesProp$Scenario %in% c("22", "23", "24", "25")], BABACI==Grouping), aes(x=variable, y=as.numeric(value), colour=variable, fill=variable, shape=ScenarioName), dodge.width=1, size=3, stroke=1)+
    scale_colour_manual(values=OutcomeSummaryCols, name="Sensitivity Test", guide="none")+
    scale_fill_manual(values=OutcomeSummaryCols, name="Sensitivity Test", guide="none")+
    scale_shape_manual(values=c(2, 17, 0, 15), name="Sensitivity Test")+
    ylab("Proportion of populations")+
    xlab(Grouping)+
    scale_y_continuous(limits=c(0, 0.805), breaks=c(0,0.2,0.4, 0.6, 0.8, 1))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
    theme(aspect.ratio=1.1, 
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, size=8),
          legend.title=element_text(size=plotfontsize),
          axis.ticks = element_line(colour="black"),
          axis.text.x = element_text(size=6),
          axis.title.x=element_text(size=plotfontsize),
          axis.title.y = YaxLab,
          axis.text.y = Ytext,
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top",
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))
}
BACIPlot <- OutcomePlotSensitivities("BACI", YaxisLabel = TRUE)
BAPlot <- OutcomePlotSensitivities("BA")
CIPlot <- OutcomePlotSensitivities("CI")
MakeLegend <- ggplot()+
  geom_beeswarm(data=AllBACICategoriesProp[AllBACICategoriesProp$Scenario %in% c("22", "23", "24", "25")], aes(x=variable, y=as.numeric(value)*100, shape=ScenarioName), dodge.width=1, size=3, stroke=1)+
  scale_shape_manual(values=c(2, 17, 0, 15), name="Sensitivity Test")+
  guides(shape = guide_legend(ncol = 2))+
  theme(legend.text = element_text(size=6), #family="Helvetica Neue", 
        legend.key=element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.32,0.6))
AllPlot <- plot_grid(BACIPlot, BAPlot, CIPlot, ncol=3, labels=c("a", "b", "c"), label_size = 8, align="hv", label_x = 0.02, label_y = 1.01)
AllPlot2 <- plot_grid(AllPlot, get_legend(MakeLegend), nrow=2, rel_heights = c(0.8, 0.2))
ggsave(paste0(FiguresFP, "FigS8_SensitivityBABACI.pdf"), AllPlot2, device="pdf", width = 183, height = 90, units = "mm") #Save as a pdf

### Number of species and sites in sensitivity analyses

### N for Latin HyperCubes
#Get number of sites/species/populations in each latin hypercube sample
AllN <- rbindlist(pblapply(c(22:25), function(Scen){
  load(paste0(ResultsFP, "LHC/Scenario_", Scen, "/BACI/BACI.RData"))
  return(data.frame(Scenario=Scen, NSites=length(unique(subset(BACI, CI==1)$SiteCode)), NSpec=length(unique(BACI$Species)), NPop=length(unique(BACI$SpecMatch))))
}))

AllN$Scenario <- as.character(AllN$Scenario)

AllN[AllN$Scenario == "22"]$Scenario <- "11 to 15 years sampled before and after protection"
AllN[AllN$Scenario == "23"]$Scenario <- "16 to 20 years sampled before and after protection"
AllN[AllN$Scenario == "24"]$Scenario <- "p < 0.1 to detect impact of protection"
AllN[AllN$Scenario == "25"]$Scenario <- "p < 0.2 to detect impact of protection"

write.csv(AllN, paste0(FiguresFP, "AllN_Sensitivity.csv"), row.names=FALSE) #Save

### Site locations for sensitivity analysis
for(Scens in list(c(22, 23), c(24, 25))){
  AllScenSites <- rbindlist(pblapply(Scens, function(Scen){
    DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
    load(paste0(DatabaseFP, "BACI/BACI.RData"))
    BACISites <- unique(BACI[,c("SiteCode", "CI", "Latitude", "Longitude")])
    BACISites$CI <- factor(BACISites$CI)
    BACISites$CI <- recode_factor(BACISites$CI, "1"="Protected", "0"="Unprotected")
    BACISites$Scenario <- Scen
    return(BACISites)
  }))
  
  AllScenSitesCast <- dcast(AllScenSites, SiteCode + CI + Latitude + Longitude ~ . ,length) #Get the number of LHC scenarios for each site
  AllScenSitesCast <- AllScenSitesCast[order(AllScenSitesCast$CI),]
  AllScenSitesCast$. <- AllScenSitesCast$./21 #Divide the . column (i.e. the number of LHC scenarios) by 21, to get a proportion
  names(AllScenSitesCast)[[5]] <- "Scenarios" #Rename
  AllScenSitesCast$ScenBin <- cut(AllScenSitesCast$Scenarios, 42) #Cut the proportion column into 42 bins
  AllScenSitesCast$ScenBinCI <- paste0(AllScenSitesCast$CI, AllScenSitesCast$ScenBin) #Make a column for CI plus the bins
  AllScenOrder <- order(AllScenSitesCast$CI, AllScenSitesCast$ScenBin, decreasing=c(TRUE, FALSE), method="radix") #Get the right order so it plots right
  AllScenSitesCast <- AllScenSitesCast[AllScenOrder,] #Reorder
  
  #To make the legend we need to manually provide the colours values at various alpha values
  highcolUP <- "#4f3384ff"      
  highcolP <- "#056d16ff"      
  
  #Get alphas for unprotected sites
  AllAlphasUP <- c(sapply(seq(0.4, 1, length.out=1), function(amin){
    lowcol.hexUP <- as.hexmode(round(col2rgb(highcolUP) * amin + 255 * (1 - amin)))
    lowcolUP <- paste0("#",   sep = "",
                       paste(format(lowcol.hexUP, width = 2), collapse = ""))
    return(lowcolUP)
  }), highcolUP)
  
  #Get alphas for protected sites
  AllAlphasP <- c(sapply(seq(0.4, 1, length.out=1), function(amin){
    lowcol.hexP <- as.hexmode(round(col2rgb(highcolP) * amin + 255 * (1 - amin)))
    highcolP <- paste0("#",   sep = "",
                       paste(format(lowcol.hexP, width = 2), collapse = ""))
    return(highcolP)
  }), highcolP)
  
  mapWorld <- borders("world", colour="gray80", fill="gray80")
  mapEurope <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey80", fill="grey80")
  
  #Plot world map
  World  <- ggplot(AllScenSitesCast) + mapWorld +
    geom_point(aes(x=Longitude, y=Latitude, colour=ScenBinCI), size=1, shape=1, stroke=0.8)+
    scale_colour_manual(values=c(AllAlphasP, AllAlphasUP))+ #, guide=FALSE
    coord_proj("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")+
    theme(
      aspect.ratio=0.53, 
      panel.grid = element_blank(), 
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0),
      legend.justification = "top",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background=element_rect(fill = "white"),
      panel.background = element_rect(fill = "transparent"),
      legend.key=element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.text = element_text(size=25),
      legend.title = element_blank(),
      legend.position = "none")#c(0.935, 1)0.6 #c(0.07, 1)
  World
  
  #Plot Europe inset
  Europe  <- ggplot(AllScenSitesCast) + mapEurope +
    geom_point(aes(x=Longitude, y=Latitude, colour=ScenBinCI), size=1, shape=1, stroke=0.8)+
    scale_colour_manual(values=c(AllAlphasP, AllAlphasUP))+ #, guide=FALSE
    coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(size = 1, fill = NA, colour="black"),
      legend.justification = "none",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.background=element_rect(fill = "transparent", colour=NA),
      panel.background = element_rect(fill = "white"))
  Europe
  #Make legend
  GLegend <- ggplot(unique(AllScenSitesCast[,c("CI", "ScenBin", "ScenBinCI")]), aes(CI,ScenBin))+
    geom_tile(aes(fill=ScenBinCI))+
    ylab("Number of analyses")+
    scale_fill_manual(values=c(AllAlphasP, AllAlphasUP))+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0), labels=c(1, 2))+
    coord_flip()+
    theme(legend.position="none", aspect.rat=1, 
          axis.title=element_text(size=7),
          plot.margin=margin(t=10,b=10,l=10),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.ticks = element_blank(),
          axis.title.y=element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank())
  GLegend
  
  #Save
  pdf(paste0(FiguresFP, "SiteMaps/Sensitivities_", paste0(Scens, collapse="_"), ".pdf"), (183/25.4), (95/25.4), bg="transparent") #Save as an image res=300, For some reason pdf ony recieves measurements in inches so ahve to convert
  grid.newpage()
  v1 <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
  v2 <- viewport(width = 0.51, height = 0.47, x = -0.12, y = 0.07, just=c("left", "bottom")) #plot area for the inset map
  v3 <- viewport(width = 0.69, height = 0.23, x = 0.53, y = 0.23, just=c("left", "top"))
  print(World,vp=v1) 
  print(Europe,vp=v2, mar=c(0,0,0,0))
  print(GLegend,vp=v3, mar=c(0,0,0,0))
  dev.off()
}  

#### CLMM Proportional Odds Ratio - Check Assumptions (Figure S10, 11) ####
Database <- "Waterbird"

#Write a little function to estimate log odds ratios (with thanks to https://stats.idre.ucla.edu/r/dae/ordinal-logistic-regression/)
sf <- function(y) {
  c('Y3' = qlogis(mean(y >= 3)), 'Y2' = qlogis(mean(y >= 2)))
}

#Read in the data we input into the cumulative link models
BACIModelDataAll <- rbindlist(pbmclapply(c(1:21), function(Scenny){
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scenny, "/")
  BACICategories <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
  load(file=paste0(DatabaseFP, "BACI/BACIPredictors.RData"))
  
  BACIModelData <- PrepCategoryModelData(BACICategories, BACIPredictors, paste0(DatabaseFP, "BACI/"), SaveOutput = FALSE)
  BACIModelData$SiteSpec <- paste0(BACIModelData$SiteCode, ".", BACIModelData$Species)
  BACIModelData$Model <- "All"
  BACIModelData$Scenario <- Scenny
  return(BACIModelData)
}, mc.cores=ncores))
BACIModelDataAll <- subset(BACIModelDataAll, MigStatus != "Uncertain")

#Calculate the proportional odds ratios for each covariate
PropOddsAll <- rbindlist(pblapply(1:21, function(Scenny){
  CoefNames <- read.csv(paste0(ResultsFP, "CoefNames_PO.csv")) #This is just a csv that has full names for all the coefficients (e.g. "Protected area size" rather than "PA Area"), plus a marker for the order they should be in on the plot
  
  #Subset data to the latin hypercube scenario
  BACIModelDataAll2 <- as.data.frame(subset(BACIModelDataAll, Scenario==Scenny)) #Subset model 
  
  #Loop through each covariate and calculate
  PropOddsScen <- rbindlist(lapply(unique(CoefNames$Covar), function(Covarry){
    names(BACIModelDataAll2)[names(BACIModelDataAll2) == (Covarry)] <- "TargetVar"
    BACIModelDataAll2 <- BACIModelDataAll2[complete.cases(BACIModelDataAll2$TargetVar),]
    ContLevels <- NULL
    
    #For continuous covariates create 5 levels spanning 5 quintiles
    if(Covarry %in% c("PAAreaLog", "GenLength", "GovMean")){
      CutBreaks <- unique(quantile(BACIModelDataAll2$TargetVar))
      CutBreaks[[1]] <- CutBreaks[[1]]-0.1
      BACIModelDataAll2$TargetVar <- cut(BACIModelDataAll2$TargetVar, breaks=CutBreaks)
      ContLevels <- levels(BACIModelDataAll2$TargetVar)
      levels(BACIModelDataAll2$TargetVar) <- c(1:length(unique(BACIModelDataAll2$TargetVar)))
    }
    
    #Split the data into each level of the covariate (e.g. waterbird managed, mixed management)
    SplitDat <- split(BACIModelDataAll2, BACIModelDataAll2$TargetVar)
    
    #Use our log odds function to calculate log odds
    LogDat <- rbindlist(lapply(c(1:length(SplitDat)), function (x){
      datty <- SplitDat[[x]]
      ORs <- exp(sf(as.numeric(datty$OutcomeSummary)))
      data.frame(LevelDiff = ORs[2]/ORs[1], Level=names(SplitDat)[[x]])
    }))
    LogDat$Coef <- Covarry
    
    #For categorical covariates calculate the difference between target level and reference level
    if(Covarry=="AnthMode"){
      LogDat$Range <- abs(LogDat$LevelDiff - LogDat[Level=="Urban"]$LevelDiff)
    } else if(Covarry=="Order"){
      LogDat$Range <- abs(LogDat$LevelDiff - LogDat[Level=="Anseriformes"]$LevelDiff)
    } else if(Covarry=="Managed"){
      LogDat$Range <- abs(LogDat$LevelDiff - LogDat[Level=="Other"]$LevelDiff)
    } else if(Covarry=="MigStatus"){
      LogDat$Range <- abs(LogDat$LevelDiff - LogDat[Level=="Resident"]$LevelDiff)
    } else {
      LogDat$Range <- max(LogDat$LevelDiff)-min(LogDat$LevelDiff) #For continuous covariates calculate the range
    }
    return(LogDat)
  }))
  PropOddsScen$Scenario <- Scenny
  return(PropOddsScen)
}))


#Inefficient code to prepare data for plotting
AllModels <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "LHC/"), pattern = "CategoryModelsOutput.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  if(ncol(Dat)==1){
    return(NULL)
  }
  Dat$Scenario <- str_split_fixed(str_split_fixed(x, "[/]", 14)[,12], "[_]", 2)[,2]
  return(Dat)
}))
AllModels <- AllModels[!AllModels$Scenario %in% c("22", "23", "24", "25")]
AllModels <- subset(AllModels, Model != "RedList")
PropOddsAll$Coef2 <- PropOddsAll$Coef
PropOddsAll[!PropOddsAll$Coef2 %in% c("PAAreaLog", "GenLength", "GovMean"), "Coef2"] <- paste0(PropOddsAll[!PropOddsAll$Coef2 %in% c("PAAreaLog", "GenLength", "GovMean"),]$Coef2, PropOddsAll[!PropOddsAll$Coef2 %in% c("PAAreaLog", "GenLength", "GovMean")]$Level)
PropOddsAll$Scenario <- as.character(PropOddsAll$Scenario)
PropOddsAll <- subset(PropOddsAll, LevelDiff!="NaN")
PropOddsAll <- merge(PropOddsAll, subset(AllModels, Model=="All"), by.x=c("Coef2", "Scenario"), by.y=c("Coef", "Scenario"))
PropOddsAll$Scenario <- factor(PropOddsAll$Scenario, levels=c(20:1, 21))
levels(PropOddsAll$Scenario) <- c(20:1, "Focal Analysis")

PropOddsAll <- subset(PropOddsAll, Range!="Inf")
CoefNames <- read.csv(paste0(ResultsFP, "CoefNames.csv")) #This is just a csv that has full names for all the coefficients (e.g. "Protected area size" rather than "PA Area"), plus a marker for the order they should be in on the plot
PropOddsAll <- merge(PropOddsAll, CoefNames, by.x="Coef2", by.y="Coef")
PropOddsAll <- PropOddsAll[order(PropOddsAll$OrderMarker),]
PropOddsAll$CoefName <- factor(PropOddsAll$CoefName, levels=unique(PropOddsAll$CoefName))
PropOddsAll$CoefGroup <- factor(PropOddsAll$CoefGroup, levels=unique(PropOddsAll$CoefGroup))

#Plot
PropOddsRange <- ggplot(PropOddsAll)+
  geom_point(aes(x=Scenario, y=Range, colour=CoefGroup))+
  facet_wrap(~Coef2)+
  facet_wrap(~CoefName, nrow=3)+ #
  #geom_hline(yintercept = 0.25, linetype="dashed")+
  scale_y_continuous(breaks = trans_breaks(identity, identity, n = 3))+
  scale_fill_manual(values=brewer.pal(length(levels(PropOddsAll$CoefGroup)), name="Dark2"), name="Covariate Group")+
  scale_colour_manual(values=brewer.pal(length(levels(PropOddsAll$CoefGroup)), name="Dark2"), name="Covariate Group")+
  xlab("Analysis")+
  ylab("Range in Log Odds Ratios")+
  coord_flip()+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+
  theme(aspect.ratio=1.25, 
        panel.background = element_blank(),
        panel.grid.major.y = element_blank() ,
        panel.grid.major.x = element_line(color = "grey90") ,
        panel.grid.minor.x = element_line(color = "grey90") ,
        panel.grid.minor.y = element_blank() ,        
        text = element_text(size=plotfontsize-2),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize-2),
        legend.title=element_text(size=plotfontsize),
        legend.key = element_blank(),
        axis.ticks = element_line(colour="black"),
        axis.text.x = element_text(size=plotfontsize-2, colour="black"),
        axis.text.y = element_text(size=plotfontsize-2, colour="black"),
        axis.title.y = element_text(size=plotfontsize, colour="black"),
        axis.title.x = element_text(size=plotfontsize, colour="black"),
        legend.text = element_text(size=plotfontsize-1, colour="black"),
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.position = c(-0.10,-0.14),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(0.5, "cm"),
        plot.margin=unit(c(0.2,0.2,1.5,0.2),"cm"))
PropOddsRange

ggsave(paste0(FiguresFP, "FigS10_PropOddsRange.pdf"), PropOddsRange, device="pdf", width = 183, height = 137, units = "mm") #Save as a pdf

PropOddsAll$Thresh <- ifelse(PropOddsAll$Range > 0.5, 1, 0)
PropAboveThresh <- dcast(PropOddsAll, Coef2 ~ Thresh, length)
PropAboveThresh$Prop <- PropAboveThresh$`1`/rowSums(PropAboveThresh[,c(2,3)])

### Test for whether our conclusions in the main paper are affected by less accurate models
#Regress P value/coefficient estimate by range in proportional odds (see Supporting Information 13)
AssumptionsAssoc <- rbindlist(lapply(unique(PropOddsAll$Coef2), function(x){
  LMDat <- subset(PropOddsAll, Coef2==x)
  LMDat$PosNeg <- ifelse(LMDat$Estimate>0, 1, 0)
  Est <- as.data.frame(t(coef(summary(lm(LMDat$Estimate ~ LMDat$Range)))[2,c(1,4)]))
  Est$Corr <- "Estimate"
  
  EstBin <- as.data.frame(t(coef(summary(glm(LMDat$PosNeg ~ LMDat$Range, family="binomial")))[2,c(1,4)]))
  EstBin$Corr <- "EstBinom"
  
  P <- as.data.frame(t(coef(summary(lm(LMDat$P ~ LMDat$Range)))[2,c(1,4)]))
  P$Corr <- "P"
  
  LMOut <- rbindlist(list(Est, P, EstBin), use.names=FALSE)
  LMOut$Coef2 <- x
  names(LMOut) <- c("Estimate", "P", "Corr", "Coef2")
  return(LMOut)
}))

#For models with 50% of points above 0.5 difference in prop odds ratio OR with a significant binomial relationship between posneg and Range (meaning the assumption affects whether we make a positive or negative estimate), run separate binomial models
union(subset(AssumptionsAssoc, P<0.05 & Corr=="EstBinom")$Coef2, subset(PropAboveThresh, Prop>0.5)$Coef2)

#This code is basically the same as the cumulative link models, just split so taht there's one model for Negative to No impact, and one for no to positive impact
BinomLogisticModels <- pbmclapply(c(1:21), function(Scenny){
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scenny, "/")
  BACICategories <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
  load(file=paste0(DatabaseFP, "BACI/BACIPredictors.RData"))
  
  BACIModelData <- PrepCategoryModelData(BACICategories, BACIPredictors, paste0(DatabaseFP, "BACI/"), SaveOutput = FALSE)
  BACIModelData$SiteSpec <- paste0(BACIModelData$SiteCode, ".", BACIModelData$Species)
  BACIModelData$Model <- "All"
  
  OutcomeData <- CategoryModelFunction(BACIModelData[!is.na(BACIModelData$PAAreaLog),], PAArea=TRUE, Scenarios[Scenny,]$SitePairDist)
  
  NegtoNeutDat <- OutcomeData[[1]]
  NegtoNeutDat <- subset(NegtoNeutDat, OutcomeSummary != "Positive Impact")
  NegtoNeut <- glmer(OutcomeSummary ~ MigStatus + Order + GenLength + AnthMode + (1 | ISO3:SiteCode) + (1 | Grid:SiteCode) + 
                       (1 | Species), data=NegtoNeutDat, family = binomial(link = "logit"))
  NegtoNeutDat <- as.data.frame(summary(NegtoNeut)$coefficients)
  NegtoNeutDat$Coef <- row.names(NegtoNeutDat)
  NegtoNeutDat$Model <- "NegtoNeut"
  
  NeuttoPosDat <- OutcomeData[[1]]
  NeuttoPosDat <- subset(NeuttoPosDat, OutcomeSummary != "Negative Impact")
  NeuttoPos <- glmer(OutcomeSummary ~ MigStatus + Order + GenLength + AnthMode + (1 | ISO3:SiteCode) + (1 | Grid:SiteCode) + 
                       (1 | Species), data=NeuttoPosDat, family = binomial(link = "logit"))
  NeuttoPosDat <- as.data.frame(summary(NeuttoPos)$coefficients)
  NeuttoPosDat$Coef <- row.names(NeuttoPosDat)
  NeuttoPosDat$Model <- "NeuttoPos"

  AllModelDat <- rbind(NegtoNeutDat, NeuttoPosDat) #, NegtoPosDat
  AllModelDat$Scenario <- Scenny
  return(AllModelDat)
}, mc.cores=5)
BinomLogisticModels2 <- rbindlist(BinomLogisticModels)
BinomLogisticModels2 <- subset(BinomLogisticModels2, Coef!="(Intercept)")
names(BinomLogisticModels2) <- c("Estimate", "Error", "Z", "P", "Coef", "Model", "Scenario")
BinomLogisticModels2$Scenario <- as.character(BinomLogisticModels2$Scenario)
BinomLogisticModels2[BinomLogisticModels2$Scenario==21]$Scenario <- "Focal Analysis"

#This is a bunch of fairly (very) inefficient cleaning of all the predictor covariates, how to order the factors etc, creating binary columns to mark whether estimates are significant or not, and positive or negative
CoefNames <- read.csv(paste0(ResultsFP, "CoefNames.csv")) #This is just a csv that has full names for all the coefficients (e.g. "Protected area size" rather than "PA Area"), plus a marker for the order they should be in on the plot
BinomLogisticModels2 <- merge(BinomLogisticModels2, CoefNames, by="Coef") #Add to model data

BinomLogisticModels2$CoefName <- factor(BinomLogisticModels2$CoefName, levels=unique(BinomLogisticModels2$CoefName))
BinomLogisticModels2$CoefGroup <- factor(BinomLogisticModels2$CoefGroup, levels=unique(BinomLogisticModels2$CoefGroup))
BinomLogisticModels2$Scenario <- factor(BinomLogisticModels2$Scenario, levels=c(c(20:1), "Focal Analysis"))
BinomLogisticModels2$PosNeg <- ifelse(BinomLogisticModels2$Estimate>0, "Pos", "Neg")
BinomLogisticModels2$Sig <- ifelse(BinomLogisticModels2$P<0.05, "Sig", "Insig")
BinomLogisticModels2$CoefName <- recode_factor(BinomLogisticModels2$CoefName, "PA Area"="Protected area size", "Ramsar or SPA Site"="Managed for waterbirds")
BinomLogisticModels2$ContVar <- "Categorical"
BinomLogisticModels2[BinomLogisticModels2$CoefName=="Body Size" | BinomLogisticModels2$CoefName=="Protected area size" |BinomLogisticModels2$CoefName=="Governance"]$ContVar <- "Continuous"

BinomLogisticModels2$PosNegSig <- paste0(BinomLogisticModels2$PosNeg, "_", BinomLogisticModels2$Sig) #Create one value that expressed whether the value is positive or negative and significance or insignificant
BinomLogisticModels2 <- dcast(BinomLogisticModels2, CoefName + CoefGroup + Model ~ PosNegSig, length) #Count the number of LHC scenarios that have each kind of response (Pos/Neg & Significant/Insig) by covariate
BinomLogisticModels2 <- melt(BinomLogisticModels2, id.vars=c("CoefName", "CoefGroup", "Model"))
BinomLogisticModels2$variable <- factor(BinomLogisticModels2$variable, levels=c("Neg_Sig", "Neg_Insig", "Pos_Insig", "Pos_Sig")) #Define as factor
BinomLogisticModels2$variable <- recode_factor(BinomLogisticModels2$variable, "Neg_Sig" = "- (p < 0.05)", "Neg_Insig"="- (p > 0.05)","Pos_Insig" = "+ (p > 0.05)", "Pos_Sig"="+ (p < 0.05)") #Rename Factor levels for legend
OutcomeSummaryCols2 <- OutcomeSummaryCols[c(1:4)] #Get colours
names(OutcomeSummaryCols2) <- c("- (p < 0.05)", "- (p > 0.05)", "+ (p > 0.05)", "+ (p < 0.05)") #Name colours
BinomLogisticModels2$CoefGroup <- factor(BinomLogisticModels2$CoefGroup, levels=c("Area", "Managed", "Governance", "Body Size", "Migrant Status","Anthrome", "Order"))

BinomLogisticModels2$Model <- factor(BinomLogisticModels2$Model)
BinomLogisticModels2$Model <- recode_factor(BinomLogisticModels2$Model, "NegtoNeut" = "No vs Negative Impact", "NegttoPos" = "Positive vs Negative Impact", "NeuttoPos" = "Positive vs No Impact")

#Plot (This plot is the same as Figure 4 but now split by the three models - Negative vs. No Impact, Negative vs Positive Impact and No Impact vs Positive Impact. )
Plot <- ggplot(BinomLogisticModels2, aes(x=CoefName, y=value, fill=variable, colour=variable))+
  geom_bar(stat="identity", position="stack", size=0, colour="transparent", width=1)+
  facet_grid(CoefGroup~Model, scales="free", space='free')+
  scale_fill_manual(values = OutcomeSummaryCols2, name=("Estimate"))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  ylab("Number of analyses")+
  coord_flip()+
  theme(#aspect.ratio=1, 
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    text = element_text(size=plotfontsize),
    panel.border = element_rect(size = 0.8, fill = NA, colour="grey60"),
    #panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(hjust = 0, size=6),
    strip.text.y = element_blank(),
    axis.ticks = element_line(colour="black"),
    axis.text.x = element_text(size=6, colour="black"),
    axis.text.y = element_text(size=6, colour="black"),
    axis.title.x = element_text(size=7, colour="black"),
    axis.title.y = element_blank(),
    #axis.line = element_line(),
    legend.text = element_text(size=6, colour="black"),
    legend.title = element_text(size=7, colour="black"),
    legend.justification = "top",
    legend.key.size = unit(0.4, 'cm'))
Plot
ggsave(paste0(FiguresFP, "Inkscape/BinomialLogisticModels.pdf"), Plot, device="pdf", width = 183, height = 73, units = "mm", dpi=1000) #Save as a pdf

#Estimates of PA Size and Management are significantly negatively associate with POR range. 
#Therefore, we are most confident about the models with small POR differences. Let's get estimates of effect size and weight them by POR difference. 

#Get quintiles for PA area
PASizeQuantiles <- paste0(quantile(BACIModelDataAll[complete.cases(BACIModelDataAll$PAAreaLog),]$PAAreaLog, c(0.05, 0.25, 0.5, 0.75, 0.95)), collapse = ",")

MSDone <- "Done"
if(MSDone != "Done"){
  
  #Re rerun the model for each scenario, and then use ggpredict to get effect sizes (these are output as the proportion likelihood of a negative/no/positive impact occurring, dependign on the predictor)
  ManagementSizeEffectSizes <- pbmclapply(c(1:21), function(Scenny){
    DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scenny, "/")
    BACICategories <- read.csv(file=paste0(DatabaseFP, "BACI/ModelOutput/BACICategories.csv"))
    load(file=paste0(DatabaseFP, "BACI/BACIPredictors.RData"))
    
    BACIModelData <- PrepCategoryModelData(BACICategories, BACIPredictors, paste0(DatabaseFP, "BACI/"), SaveOutput = FALSE)
    BACIModelData$SiteSpec <- paste0(BACIModelData$SiteCode, ".", BACIModelData$Species)
    BACIModelData$Model <- "All"
    
    OutcomeData <- CategoryModelFunction(BACIModelData[!is.na(BACIModelData$PAAreaLog),], PAArea=TRUE, Scenarios[Scenny,]$SitePairDist)
    CatModel <- clmm(OutcomeData[[2]], data=OutcomeData[[1]]) #Run the CLMM
    if(is.na(summary(CatModel)[[1]][1,2])){
      CatModel <- clmm(OutcomeData[[3]], data=OutcomeData[[1]]) #Run the CLMM
      ModelRun <- "NoGrid"
    }
    
    PredDat <- as.data.frame(ggpredict(CatModel, c("Managed", paste0("PAAreaLog[", PASizeQuantiles, "]")), type="re"))
    PredDat$Scenario <- Scenny
    return(PredDat)
  }, mc.cores=ncores)
  ManagementSizeEffectSizes <- rbindlist(ManagementSizeEffectSizes)
  save(ManagementSizeEffectSizes, file=paste0(ResultsFP, "LHC/", "ManageSizeEffectSizes.RData"))
}
load(file=paste0(ResultsFP, "LHC/", "ManageSizeEffectSizes.RData"))

#Rename data, calculate difference in likelihood of positive outcome (Level 3) between small mixed management and large managed (Across each scenario)
names(ManagementSizeEffectSizes) <- c("Managed", "Pred", "Error", "ConfLow", "ConfHigh", "OutcomeSummary", "PAAreaLog", "Scenario")
EffectSize <- dcast(ManagementSizeEffectSizes, OutcomeSummary + Scenario ~ Managed + as.numeric(PAAreaLog), value.var="Pred")
EffectSizePos <- subset(EffectSize, OutcomeSummary=="3")
EffectSizePos$SmallUnBigMan <- EffectSizePos$Managed_5 - EffectSizePos$Other_1
EffectSizePos$Scenario <- as.character(EffectSizePos$Scenario)
EffectSizePos <- merge(EffectSizePos, unique(subset(PropOddsAll[,c("Scenario", "Range", "Coef")], Coef=="Managed" | Coef=="PAAreaLog")), by="Scenario")

EffectSizePosLevel <- dcast(EffectSizePos[,c("Scenario", "SmallUnBigMan", "Range", "Coef")], Scenario + SmallUnBigMan ~ Coef, value.var="Range")
EffectSizePosLevel$RangeWeight <- 1/(EffectSizePosLevel$Managed * EffectSizePosLevel$PAAreaLog)
EffectSizePosLevel$RangeWeight2 <- 1-((EffectSizePosLevel$RangeWeight-(min(EffectSizePosLevel$RangeWeight)))/(max(EffectSizePosLevel$RangeWeight) - min(EffectSizePosLevel$RangeWeight)))

#Calculate max, mini and weight mean
max(EffectSizePosLevel$SmallUnBigMan)
min(EffectSizePosLevel$SmallUnBigMan)
mean(EffectSizePosLevel$SmallUnBigMan)
weighted.mean(EffectSizePosLevel$SmallUnBigMan, w=EffectSizePosLevel$RangeWeight2)

