############################################################################################################################################################################
### This code was written by Hannah Wauchope to analyse data for the paper "Protected areas have a mixed impact on waterbirds, but management helps"
### This is script 3 of 4 in the workflow
### Last edited 17th October, 2021
### Please direct queries to hannah.wauchope@gmail.com
###
### This script subsets the data to the requirements for analysis, identifies the protected areas that sites fall within and divides data into 20 Latin Hypercube Sample analyses, plus a focal analysis
### Finally, it conducts matching according BACI and CI protocols. The quality of the matches are assessed in Script 4
###
### NOTES: Throughout these scripts I use the field "SiteSpec" to mean Population, i.e. a particular species at a particular site. This is generally the unit that analyses are performed on
############################################################################################################################################################################

#### Initialise ####
cluster <- FALSE
if(cluster==TRUE){
  library(data.table, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(rgdal, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(sp, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(pbapply, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(MASS, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(plyr, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(pbmcapply, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(usdm, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(rlist, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(resample, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(scales, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(foreign, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(DOS, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
  library(fields, lib.loc="/gpfs/ts0/projects/Research_Project-T115802/RPackages/")
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
  DataFP <- "/gpfs/ts0/projects/Research_Project-T115802/data/"
  ncores <- 16 
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
  library(scales)
  library(maptools)
  library(DOS)
  library(rgeos) #not on cluster
  library(fields) #not on cluster
  library(car) 
  library(DHARMa) 
  library(lmtest) 
  library(nlme) 
  library(lme4) 

  ncores <- 6
  
  ###File paths
  DataFP <- "/Users/hannahwauchope/Dropbox/Work/Projects/PhD/Data/"
  ResultsFP <- "/Users/hannahwauchope/Dropbox/Work/Projects/PhD/Chapter2_4_PAs/Analysis/"
}

#### Define functions ####
#Define map projections
MollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
WGSCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#Defint font size to be used in plots
plotfontsize <- 12

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

#Function for me to just check the number of Sites/Species/Populations ("SiteSpec") as I go. 
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpec'), function(x) length(unique(Dataset[,c(x)])))
}

#Function to assess models for temporal autocorrelataion (TAC) - doesn't work on models with random effects. Requires mode object (Mod), the data input into the model (ModDat), whether the model is a BACI model or no (defaults to no), and also the TotalYearsBuffer
TACCheck <- function(Mod, ModDat, BABACI="BA", TotalYearsBuffer){
  ModDat$Year2 <- ModDat$Year
  
  if(BABACI=="BACI"){
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

#This is a function used in the matching process to remove collinear variables. It recieves a dataset with a column for SiteCode, a column for Year, and a column for each of the various matching covariates. It uses Variance Inflation Factor (VIF; values above 4) and Pearson's Correlation Coefficient (values above 0.7) to identify collinear variables
#The function returns a list of two elements, the first is a dataframe reduced to only non-collinear variables, and the second is a vector listing those variable names
#Scale is logical, if true all covariates are rescaled to within 1 and 100.
#Variables is a vector specify the variable names to be considered by the function
CollVIF <- function(dataset, variables, scale){
  CorrPred <- dataset[colnames(dataset) %in% variables]
  CorrPred <- CorrPred[complete.cases(CorrPred),]
  if(nrow(CorrPred)==0){return(NULL)}
  
  #First, check for any variables that are zero across the board and remove
  SumCheck <- colSums(CorrPred, na.rm = TRUE)
  SumCheck <- SumCheck[SumCheck>0]
  CorrPred <- CorrPred[,(names(CorrPred) %in% names(SumCheck))]
  
  #Next, check for any variables that are the same across the board and remove
  CheckforSameValue <- apply(CorrPred, 2, function(x) length(unique(x)))
  CheckforSameValue <- CheckforSameValue[CheckforSameValue==1]
  CorrPred <- CorrPred[,!(names(CorrPred) %in% names(CheckforSameValue))]
  
  #Next, check for collinearity using VIF. Iteratively remove variables where VIF = Inf, NaN, or is above 4 (removing the variables with the highest VIFs first)
  TheVIF <- as.data.frame(usdm::vif(CorrPred))
  while(nrow(subset(TheVIF, VIF=="Inf"))>0){
    TheVIF <- subset(TheVIF, VIF=="Inf")[1,]
    CorrPred <- CorrPred[,!(names(CorrPred) %in% TheVIF$Variables)]
    TheVIF <- as.data.frame(usdm::vif(CorrPred))
  }
  while(nrow(subset(TheVIF, VIF=="NaN"))>0){
    TheVIF <- subset(TheVIF, VIF=="NaN")[1,]
    CorrPred <- CorrPred[,!(names(CorrPred) %in% TheVIF$Variables)]
    TheVIF <- as.data.frame(usdm::vif(CorrPred))
  }
  while(max(TheVIF$VIF)>4){
    TheVIF <- subset(TheVIF, VIF==max(TheVIF$VIF))
    CorrPred <- CorrPred[,!(names(CorrPred) %in% TheVIF$Variables)]
    TheVIF <- as.data.frame(usdm::vif(CorrPred))
  }
  
  #Next, check for collinearity using Pearson's correlation coefficient, remove one of any pairs where this is above 0.7
  CorrelationPearson <- cor(CorrPred, method="pearson")
  Pearson <- apply(CorrelationPearson, 2, function(x) x[x>0.70])
  while(max(sapply(Pearson, length))>1){
    PearsonCheck <- sapply(Pearson, function(x) if(length(x)>1){return(NULL)} else {return(unique(x))})
    Pearson <- Pearson[sapply(PearsonCheck, is.null)]
    CorrPred <- CorrPred[,!(names(CorrPred) %in% names(Pearson)[[1]])]
    CorrelationPearson <- cor(CorrPred, method="pearson")
    Pearson <- apply(CorrelationPearson, 2, function(x) x[x>0.70])
  }
  
  if(scale==TRUE){
    CorrPredScale <- dataset[colnames(dataset) %in% colnames(CorrPred)]
    CorrPred <- as.data.frame(apply(CorrPredScale, 2, function (x) rescale(x, to = c(0, 100))))
  }
  
  #Remove any variables that are mostly the same value (at least 80% of values are equal, occurs in a few unique cases)
  SameValue <- apply(CorrPred, 2, function(x) length(x[x==Mode(x)])/length(x))
  SameValue <- SameValue[SameValue<0.8]
  VariablesForModel <- names(SameValue)
  CorrPred <- CorrPred[names(CorrPred) %in% VariablesForModel]
  
  dataset2 <- cbind(dataset[!colnames(dataset) %in% variables], CorrPred)
  return(list(dataset2, VariablesForModel))
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

#This is a function to ascertain whether two time series are parallel in slope or not in the before period (to meet an assumption of matching, see Wauchope et al., 2021, Trends in Ecology and Evolution)
#The function compares Comp1, ONE (only 1!) population time series where CI = 1, to Comp2, one or more population time series where CI = 0, and calculate whether the slopes are parallel between Comp1 and each of the population time series in Comp2
#It returns ONLY those SiteSpec from Comp2 that are parallel to Comp1 (if there are no parallel SiteSpec in Comp2, returns NULL)
#This is done by iteratively running glms between comp1 and each time series in comp2 and seeing if there is a significant difference between the two.
#TotalYearsBuffer required for the nested function that checks temporal autocorrelation.
#ZeroThresh is one of the values defined in the latin hypercube sampling, a value between 0.6 and 0.9. If the proportion of zeroes in the time series if above the threshold, it is classified as "All Zeroes" (see Wauchope et al., 2021, Trends in ecology and evolution for an explanation of the rationale behind this)

ParallelSlopes <- function(Comp1, Comp2, ParallelP, ZeroThresh, parallelise, TotalYearsBuffer){
  if(parallelise==TRUE){
    ncoresslopes <- ncores
  } else {
    ncoresslopes <- 1
  }
  
  Counts1 <- subset(Comp1, BA==0) #This function is only ever relevant to the before period, so subset to then
  Counts2 <- subset(Comp2, BA==0)
  Counts1Zeroes <- AllZeroes(Counts1, ZeroThresh) #First identify any cases that are "AllZero
  Counts2Zeroes <- AllZeroes(Counts2, ZeroThresh)
  if(Counts1Zeroes$AllZero=="AllZero"){ #If Comp1 is a case of AllZeroes, subset to only cases of Comp2 that are also AllZero, return these 
    Counts2Zeroes <- subset(Counts2Zeroes, AllZero=="AllZero")
    return(Counts2Zeroes)
  }
  
  Counts2Zeroes <- subset(Counts2Zeroes, AllZero=="Counts")
  if(nrow(Counts2Zeroes)==0){
    return(Counts2Zeroes)
  }
  
  #Iteratively run through the population time series in Comp2 and run GLMs with each with Comp1. If CI is significant in the model, then the slopes are significantly different.
  #Returns a dataframe with each SiteSpec in Comp2, and the estimate and significance of the CI term
  Parallel <- rbindlist(pbmclapply(unique(Counts2Zeroes$SiteSpec), function(SS){
    CountsSub <- subset(Counts2, SiteSpec==SS)
    Counts <- rbindlist(list(Counts1, CountsSub), fill=TRUE)[,c("Count", "Year", "CI", "BA", "Hours", "Dataset")]
    Counts[is.na(Counts$Hours),]$Hours <- 1
    Counts <- Counts[order(Counts$CI, Counts$Year),]
    
    Counts$Year <- Counts$Year - min(Counts$Year) + 1 #Scale years in case of need to use random function

    CountsComp <- suppressWarnings(tryCatch(glm.nb(Count~ Year + CI + Year*CI + offset(log(Hours)), link=log, data=Counts), error=function(e){NA}))
    slope <- tryCatch(coef(CountsComp, silent=TRUE)[[4]], error=function(e){NA})  
    significant <- tryCatch(summary(CountsComp)$coeff[4,4], error=function(e){NA})
    random <- 0
    #Check for temporal autocorrelation. If it exists, run model again with a random factor on year and recalculate slope and significance. (first checking if the first model failed, if so, skip)
    if(length(CountsComp)!=1){
      if(min(TACCheck(CountsComp, Counts, BABACI = "BACI", TotalYearsBuffer)) < 0.05){
        CountsComp <- suppressWarnings(tryCatch(glmer.nb(Count~ Year + CI + Year*CI + (1|CI:Year) + offset(log(Counts$Hours)), data=Counts), error=function(e){NA}))
        slope <- tryCatch(summary(CountsComp)$coeff[4,1], error=function(e){NA})  #Extract the slope
        significant <- tryCatch(summary(CountsComp)$coeff[4,4], error=function(e){NA}) #Extract the significance
        random <- 1
      }
    }

    return(data.frame("Sig"=significant, "Slope"=slope, "SiteSpec"=SS, "random" = random))
  }, mc.cores=ncoresslopes))
  Parallel <- subset(Parallel[complete.cases(Parallel),], Sig>ParallelP) #Subset to only site spec from Comp2 that are not significantly different from Comp1 (ParallelP always equals 0.05, defined below)
  return(Parallel)
}

#### Load and Prepare Count Data ####
load(file=paste0(DataFP, "WaterbirdData_2020/WaterbirdCounts_AllCovariates.RData"))
WaterbirdCounts$Database <- "Waterbird"
WaterbirdCounts$Class <- "Aves"
WaterbirdCounts$SiteSpec <- paste0(WaterbirdCounts$SiteCode, ".", WaterbirdCounts$Species)

write.csv(unique(WaterbirdCounts[,c("SiteCode", "Latitude", "Longitude")]), paste0(DataFP, "WaterbirdData_2020/WaterbirdSiteCoordinates.csv"), row.names=FALSE)

CountData <- subset(WaterbirdCounts, Season=="DectoFeb")
CountData$SiteSpec <- paste0(CountData$SiteCode, ".", CountData$Species)

#Remove SiteSpec with less than 6 years, just to speed everything up (we never use populations that have less than this, because of our minimum requirement of at least 3 years before and after protection, = 6 years all up)
YearsPerSiteSpec <- dcast(CountData, SiteSpec~., length, value.var="Year")
YearsPerSiteSpec <- subset(YearsPerSiteSpec, . >= 6)
CountData <- CountData[CountData$SiteSpec %in% YearsPerSiteSpec$SiteSpec,]

save(CountData, file=paste0(ResultsFP, "CountData.RData"))

#### Extract Protected Area Data ####
#Get site locations and convert to spatial points dataframe
load(file=paste0(ResultsFP, "CountData.RData"))

Sites <- unique(CountData[,c("SiteCode", "Longitude", "Latitude")])
SitesPoints <- cbind(Sites$Longitude, Sites$Latitude)
SitesPoints <- SpatialPointsDataFrame(SitesPoints, Sites, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
rm(CountData) #This object is massive so remove to speed things up
gc()

##Snap sites to Protected Planet Base Layer. This is the land cover layer used by protected planet - some coordinates in our site data lie on the coast just off the land cover layer and so might be marked as unprotected despite being on the edge of a coastal PA. To resolve this, all sites are snapped with the protected planet base layer, with any snapped more than 10km removed.
#Snapped to the base layer is the right thing to do, rather than snapping to the actual protected area shapefile - this would mean that many sites would be snapped that do NOT occur in protected areas!
Snap <- "Done" #This takes a LONG time so is commented out so it doesn't run unless I want it to
if(Snap!="Done"){
  Base <- readOGR(paste0(DataFP,"ProtectedPlanet/BaseLayer/gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp")) #Read in the base layer used by protected planet  NOAA. Global Self-consistent, Hierarchical, High-resolution Geography Database (GSHHG). (2017).)
  BaseOverlap <- over(SitesPoints, Base) #Overlay the site points with this
  BaseOverlap$SiteCode <- SitesPoints$SiteCode
  SitesPointsNA <- SitesPoints[SitesPoints$SiteCode %in% BaseOverlap[is.na(BaseOverlap$id),]$SiteCode,] #Identify any sites that do not fall within the base layer
  SitesPointsNAMoll <- spTransform(SitesPointsNA, MollCRS) #Transform to Mollweide projection (equal area, to calcuate the 10km cut off)
  
  BaseMoll <- spTransform(Base, MollCRS) #Transform the base layer to mollweide
  BaseMollLines <- as(BaseMoll, "SpatialLinesDataFrame") #Convert to spatial lines
  
  Snap <- snapPointsToLines(SitesPointsNAMoll, BaseMollLines) #Snap the points to the lines
  save(Snap, file=paste0(DataFP, "ProtectedPlanet/BaseLayer/Snap/SnappedPoints.RData")) #Save this file
  
  SnapWGS <- spTransform(Snap, WGSCRS) #Convert back to WGS projection
  writeOGR(SnapWGS, paste0(DataFP, "ProtectedPlanet/BaseLayer/Snap/"), "SnapCheck", driver="ESRI Shapefile", overwrite=TRUE) #Save
}

#Next, calculate exactly how far each site is from each protected area (this is necessary as we remove sites that are close to, but not in, protected areas to account for spill over effects)
DistancetoPP <- "Done"
if(DistancetoPP!="Done"){
  #First, remove any snapped sites that were snapped more than 10km. Removing snaps of >10km keeps 5296/5389 sites.
  SnapWGS <- readOGR(paste0(DataFP, "ProtectedPlanet/BaseLayer/Snap/SnapCheck.shp")) #Read in the snapcheck file from above
  SnappedPoints <- SnapWGS[SnapWGS$SiteCod %in% Sites$SiteCode,]
  SnappedPointsFilter <- as.data.frame(subset(SnappedPoints, snp_dst<10000))[,c("SiteCod", "Latitud", "Longitd", "coords.x1", "coords.x2")]
  names(SnappedPointsFilter) <- c("SiteCode", "Latitude", "Longitude", "PPLongitude", "PPLatitude")
  
  SitesPointsNotSnap <- as.data.frame(SitesPoints[!SitesPoints$SiteCode %in% SnappedPoints$SiteCod,])
  SitesPointsNotSnap$PPLongitude <- SitesPointsNotSnap$Longitude
  SitesPointsNotSnap$PPLatitude <- SitesPointsNotSnap$Latitude
  SitesPointsNotSnap[,c("coords.x1", "coords.x2")] <- NULL
  
  SitesPoints <- rbind(SnappedPointsFilter, SitesPointsNotSnap) #Make a new points object of the non-snapped points (i.e. points that already overlap with the land layer of protected planet) and snapped points that were snapped <10km
  SitesPoints$PPLatitude <- as.numeric(SitesPoints$PPLatitude)
  SitesPoints$PPLongitude <- as.numeric(SitesPoints$PPLongitude)
  SitesPointsCoords <- cbind(SitesPoints$PPLongitude, SitesPoints$PPLatitude)
  SitesPoints <- SpatialPointsDataFrame(SitesPointsCoords, SitesPoints, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  writeOGR(SitesPoints, paste0(DataFP, "ProtectedPlanet/BaseLayer/Snap/"), "SitesPoints", driver="ESRI Shapefile", overwrite=TRUE) #Save this
  
  #The actual calculation I did in QGIS because it took much too long in R. This uses the NNjoin function which calculates distance between points and polygons My protocol was thus:
  
  #In QGIS, import SitesPoints and ProtectedPlanet shapefile. Convert both to mollweide (Saving in "Distance to PP" folder), restart QGIS, conduct nnjoin (make sure the whole project and both shapefiles are in mollweide, sometimes importing them makes everything go a bit janky), save as PPSitesJoined.csv
  #Protected planet shapefile is the full shapefile download from protectedplanet.net. 
  #THIS TAKES 5-6 DAYS
}

#Read in distance shapefile
PPDistance <- fread(paste0(DataFP, "ProtectedPlanet/BaseLayer/DistancetoPP/PPSitesJoined.csv"))

PPDistance <- PPDistance[,c("SiteCod", "Latitud", "Longitd", "PPLngtd", "PPLattd", "distance")]
names(PPDistance) <- c("SiteCode", "Latitude", "Longitude", "PPLongitude", "PPLatitude", "distance")

#Get just the protected sites (plus a bit of a buffer, to be safe, this isn't actually necessary and gets dealt with later) and overlay with protected planet to get ALL protected areas they overlap with, plus their data (size, status, name etc). (Currently we just have a binary "distance to PA is 0m" or not)
PPDistanceProt <- subset(PPDistance, distance<1000)

PPDistanceProtCoords <- cbind(PPDistanceProt$PPLongitude, PPDistanceProt$PPLatitude)
PPDistanceProtPoints <- SpatialPointsDataFrame(PPDistanceProtCoords, PPDistanceProt, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Read in Protected Planet (PP)
PP <- readOGR(paste0(DataFP, "ProtectedPlanet/WDPA_2020/WDPA_Jan2020-shapefile/WDPA_Jan2020-shapefile-polygons.shp"))
PPData <- read.dbf(paste0(DataFP, "ProtectedPlanet/WDPA_2020/WDPA_Jan2020-shapefile/WDPA_Jan2020-shapefile-polygons.dbf"))

#Clean PP data according to wdpar
#Ensure there are only polygons (no points) in the file
PPSf <- st_as_sf(PP)
PolygonOrPoint <- vapply(sf::st_geometry(PPSf), inherits, logical(1), c("POINT", "MULTIPOINT"))
if(unique(PolygonOrPoint)!=FALSE){stop("there are points in the WDPA polygon")}

#Overlay Sites with Protected planet. Return as a nested list: each site is a list element, and within that is a list of all the protected areas the site intersects (This is because often PAs overlap)
ProtectedSites <- over(PPDistanceProtPoints, PP, returnList = TRUE) 

#Add the site names to the ProtectedSites list
PPDistanceSites <- as.data.frame(PPDistanceProtPoints)[,c("SiteCode", "Longitude", "Latitude", "distance")]
ProtectedSites <- pblapply(1:length(ProtectedSites), function(x){
  PS <- ProtectedSites[[x]]
  if(nrow(PS)==0){
    return(NULL)
  }
  Site <- PPDistanceSites[x,]
  PS <- cbind(PS, Site)
  return(PS)
})

ProtectedSites <- rbindlist(ProtectedSites) #Turn the list into a dataframe. In doing so any sites that do not intersect with any PA are removed
write.csv(ProtectedSites, paste0(ResultsFP, "ProtectedSites_Uncleaned.csv"), row.names=FALSE)

#Get unprotected sites (by process of elimination)
UnprotectedSites <- PPDistance[!PPDistance$SiteCode %in% unique(ProtectedSites$SiteCode),] 
write.csv(UnprotectedSites, paste0(ResultsFP, "UnprotectedSites_Uncleaned.csv"), row.names=FALSE)

#### Change GeoRegions ####
### This should have been done in the covariate code, but got missed out so just updating here. Just extracting the country and "GeoRegion" (roughly continent) that each site occurs in
load(file=paste0(ResultsFP, "CountData.RData"))

Countries <- readOGR("/Users/hannahwauchope/DropBox/Work/Data/GISLayers/CountryContinentWorldSHPS/World/TMBordersSSudanSthAmericaRegion/TMBordersSSudanSthAmericaRegion.shp")
Sites <- unique(CountData[,c("SiteCode", "Latitude", "Longitude")])

SitePoints <- cbind(Sites$Longitude, Sites$Latitude)
SitePoints <- SpatialPointsDataFrame(SitePoints, Sites, proj4string = WGSCRS)

CountriesOverSites <- over(SitePoints, Countries)
SiteCountries <- cbind(as.data.frame(Sites), CountriesOverSites)

#Take the points that fall just outside of country polygons
NASites <- SiteCountries[is.na(SiteCountries$ISO3),]
NAPoints <- cbind(NASites$Longitude, NASites$Latitude)
NAPoints <- SpatialPointsDataFrame(NAPoints, NASites, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

##  Project data into an equal area projection
CountriesMoll <- spTransform(Countries, MollCRS)
NAPointsMoll <- spTransform(NAPoints, MollCRS)

## For each point, find name (and details) of nearest country
NASnap <- rbindlist(pblapply(1:nrow(NAPointsMoll), function(i) as.data.frame(CountriesMoll[which.min(gDistance(NAPointsMoll[i,], CountriesMoll, byid=TRUE)),])))
NAPoints <- as.data.frame(NAPointsMoll)
NAPoints <- cbind(NAPoints[,c(1:3)], NASnap)

SiteCountries <- rbind(SiteCountries[!is.na(SiteCountries$ISO3),], NAPoints)[,c("SiteCode", "ISO3", "NAME", "REGION", "SUBREGION")]
names(SiteCountries) <- c("SiteCode", "ISO3", "Country", "GeoRegion", "GeoSubRegion")

CountData$GeoRegion <- NULL
CountData <- merge(CountData, SiteCountries[,c("SiteCode", "GeoRegion")], by="SiteCode", all=T)

save(CountData, file=paste0(ResultsFP, "CountData.RData"))

### Great, proceed!

#### Create Latin Hypercube Samples ####
#This function creates the latin hypercube samples for the various parameters that are varied in matching and analysis (see main paper)
#Because it randomly chooses new values each time, this function is "Done" and won't be run again, and the output has been saved.
#Scenarios 22 and 23 are to check for whether PA impact can only be detected in very long time series (keeping all other parameters as in focal analysis)
#Scenarios 24 and 25 are to check for whether a p<0.05 cut off is too strict to detect PA impact (i.e. maybe the change is very slight), keeping all other parameters as in the focal analysis
CreateLHC <- "Done"
if(CreateLHC != "Done"){
  library(lhs)
  Scale <- function(Column, MinVal, MaxVal){(((Column-min(Column))/(max(Column)-min(Column))) * (MaxVal-MinVal)) + MinVal}
  Square <- randomLHS(20,7)
  Square[,1] <- round(Scale(Square[,1], 5, 15)) #Total years buffer
  Square[,2] <- Scale(Square[,2], 0.5, 5) #1 to 5km PA buffer
  Square[,3] <- Scale(Square[,3], 0.1, 0.25) #0.1 to 0.25 St Diff Threshold
  Square[,4] <- Scale(Square[,4], 0.5, 0.9) #0.5 to 0.9 PropPairs Threshold
  Square[,5] <- Scale(Square[,5], 0.6, 0.8) #0.6 to 0.8 Zero Threshold
  Square[,6] <- Scale(Square[,6], 0.6, 1) #0.6 to 1 Measured years buffer
  Square[,7] <- Scale(Square[,7], 100, 2500) #100km to 2500km Distance between matched sites
  
  Square <- as.data.frame(Square)
  names(Square) <- c("TotalYearsBuffer", "PABuffer", "StDiffThresh", "PropPairsThresh", "ZeroThresh", "MeasuredYearsBuffer", "SitePairDist")
  Square$PThresh <- 0.05
  Square[21,] <- c(10, 1, 0.25, 0.7, 0.7, 0.7, 500, 0.05)
  Square[22,] <- c(15, 1, 0.25, 0.7, 0.7, 0.7, 500, 0.05)
  Square[23,] <- c(20, 1, 0.25, 0.7, 0.7, 0.7, 500, 0.05)
  Square[24,] <- c(10, 1, 0.25, 0.7, 0.7, 0.7, 500, 0.1)
  Square[25,] <- c(10, 1, 0.25, 0.7, 0.7, 0.7, 500, 0.2)
  
  Square$MeasuredYearsBuffer <- round(Square$TotalYearsBuffer*Square$MeasuredYearsBuffer)
  Square$Scenario <- 1:nrow(Square)
  write.csv(Square, paste0(ResultsFP, "LatinSquareSamples.csv"))
}

#### Conduct Matching ####
#We now have a dataset of which sites occur within protected areas (ProtectedSites_Uncleaned.csv), one of which sites don't (UnprotectedSites_Uncleaned.csv), and a dataframe of counts of waterbirds through time at each site ("CountData.RData"). We also have 22 Latin Hypercube Sample Scenarios to analyse. 
#This section of the script I run on the cluster (along with the "Initialise and Define Functions" sections) as it takes a VERY long time, around 200 cluster hours with 16 cores. 

Scenarios <- read.csv(paste0(ResultsFP, "LatinSquareSamples.csv")) #Read in latin hypercube sample data
Scenarios$Scenario <- 1:nrow(Scenarios) #Number scenarios

#Now, run matching, looping through each scenario
for(Scen in c(1:(nrow(Scenarios)))){
  #### Define Scenario and parameters, mostly from the LHC (i.e. the "Scenarios" object) ####
  PABuffer <- Scenarios[Scen,]$PABuffer
  ZeroThresh <- Scenarios[Scen,]$ZeroThresh
  TotalYearsBuffer <- Scenarios[Scen,]$TotalYearsBuffer 
  MeasuredYearsBuffer <- Scenarios[Scen,]$MeasuredYearsBuffer 
  SitePairDist <- Scenarios[Scen,]$SitePairDist
  ParallelP <- 0.05
  dir.create(file.path(paste0(ResultsFP, "LHC/")), showWarnings = FALSE)
  DatabaseFP <- paste0(ResultsFP, "LHC/Scenario_", Scen, "/")
  dir.create(file.path(DatabaseFP), showWarnings = FALSE)
  
  Sys.sleep(sample(1:20,1)) #This is here so that parallel runs don't write over each other
  
  #I use these markers to tell the cores in parallelisation to skip a scenario if another core is already working on it
  if(file.exists(paste0(ResultsFP, "Scen", Scen, "Marker.csv"))){
    next
  }
  write.csv(NULL, paste0(ResultsFP, "Scen", Scen, "Marker.csv"))
  
  #### Prepare Data ####
  #This opens with an if clause to skip this step if it's already been done ("VariablesForMatchingByYear.RData" is a final output of this step so serves as the object to check on)
  if(!file.exists(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))){
    print(paste0(Scen, " prepare data"))
    #### Clean Protected Sites ####
    ProtectedSites <- fread(paste0(ResultsFP, "ProtectedSites_Uncleaned.csv"))
    UnprotectedSites <- fread(paste0(ResultsFP, "UnprotectedSites_Uncleaned.csv"))
    load(file=paste0(ResultsFP, "CountData.RData"))
    
    CountData <- subset(CountData, Database=="Waterbird") #Cut count data down to just waterbird data
    CountData <- subset(CountData, ISO3!="RUS") #Remove Russian sites as we do not have permission to use them

    ProtectedSites <- ProtectedSites[ProtectedSites$SiteCode %in% CountData$SiteCode,] #Cut protected sites to the same subsets
    UnprotectedSites <- UnprotectedSites[UnprotectedSites$SiteCode %in% CountData$SiteCode,] #Cut unprotected sites to the same subsets
    
    ##Clean Unprotected Sites
    UnprotectedSites <- merge(UnprotectedSites, unique(CountData[,c("SiteCode", "Dataset")]), by="SiteCode") #Pull back in which dataset the unprotected counts come from
    UnprotectedSites <- subset(UnprotectedSites, Dataset=="IWC" | Dataset=="CBC" & distance>24140/2) #For CBC sites, make sure they're at least the radius of the search circle away from a PA (see the second paragraph of "Protected (and Unprotected) Area Data" in methods of paper)
    UnprotectedSites <- subset(UnprotectedSites, distance>PABuffer*1000) #Remove unprotected sites that are closer than the PABuffer (*1000 because the buffer parameter is given in KM, but the value is measured in metres)
    UnprotectedCountData <- CountData[CountData$SiteCode %in% UnprotectedSites$SiteCode,] #Add in the counts for the relevant sites
    save(UnprotectedCountData, file=paste0(DatabaseFP, "UnprotectedCountData.RData"))
    rm(UnprotectedCountData, UnprotectedSites)
    
    ##Clean protected sites
    #Find the loss of UNESCO biosphere, proposed sites, or sites without a status or status year:
    ProtectedSitesUnique <- unique(ProtectedSites[,c("SiteCode", "WDPAID", "DESIG", "STATUS", "STATUS_YR")])
    length(unique(subset(ProtectedSitesUnique, DESIG=="UNESCO-MAB Biosphere Reserve"))$SiteCode) #563
    length(unique(subset(ProtectedSitesUnique, STATUS=="Proposed"))$SiteCode) #144
    length(unique(subset(ProtectedSitesUnique, STATUS=="Not Reported"))$SiteCode) #7
    length(unique(subset(ProtectedSitesUnique, STATUS_YR==0))$SiteCode) #246
    
    #First, remove those with no STATUS_YR, plus UNESCO biosphere reserved, proposed sites and sites without a reported status
    SitesToRemove <- unique(subset(ProtectedSites, DESIG=="UNESCO-MAB Biosphere Reserve" | DESIG_ENG=="UNESCO-MAB Biosphere Reserve" | STATUS=="Proposed" | STATUS=="Not Reported" | STATUS_YR==0)$SiteCode)
    ProtectedSites <- ProtectedSites[!ProtectedSites$SiteCode %in% SitesToRemove,]
    
    #Now, for the rest we can extract "unprotected" subsets of the data - for all years before the designation year of the PA. This will give us more options for matching later. (See methods) 
    ProtectedSitesMinYear <- dcast(ProtectedSites, SiteCode~., min, value.var="STATUS_YR") #Get the min desig year
    if(nrow(ProtectedSitesMinYear[is.na(ProtectedSitesMinYear$.),])>0){stop("Stop, something's introduced NAs to Protected Site Min Status Year")}
    names(ProtectedSitesMinYear) <- c("SiteCode", "MinDesigYear")
    UnprotectedInitially <- CountData[CountData$SiteCode %in% ProtectedSitesMinYear$SiteCode,]
    UnprotectedInitially <- subset(UnprotectedInitially, Dataset!="CBC") #Remove CBC records from this because we can't be sure they weren't in the vicinity of an protected area

    UnprotectedInitially <- merge(UnprotectedInitially, ProtectedSitesMinYear, by=c("SiteCode"), all.x=TRUE)
    if(nrow(UnprotectedInitially[is.na(UnprotectedInitially$MinDesigYear),])!=0){stop("You've lost some sites!")}
    UnprotectedInitially$DeleteRow <- ifelse(UnprotectedInitially$Year<UnprotectedInitially$MinDesigYear,0,1)
    UnprotectedInitially <- subset(UnprotectedInitially, DeleteRow==0)
    UnprotectedInitially$DeleteRow <- NULL
    save(UnprotectedInitially, file=paste0(DatabaseFP, "UnprotectedInitially.RData"))
    
    #Round numeric columns to an appropriate number
    ProtectedSites <- round_df(ProtectedSites,4)
    
    #Now, in cases where we have multiple overlapping PAs for one site, reduce to PA with earliest status year, but record other incidences of area, IUCN category and PA ID
    STATUS_YR <- dcast(as.data.table(unique(ProtectedSites[,c("SiteCode", "STATUS_YR")])), SiteCode~., min, value.var="STATUS_YR")
    names(STATUS_YR) <- c("SiteCode", "STATUS_YR")
    WDPAID <- dcast(as.data.table(unique(ProtectedSites[,c("SiteCode", "WDPAID")])), SiteCode~., fun.aggregate=function(x) paste(x, collapse = ", "), value.var="WDPAID")
    names(WDPAID) <- c("SiteCode", "WDPAID")
    IUCN_CAT <- dcast(as.data.table(unique(ProtectedSites[,c("SiteCode", "IUCN_CAT")])), SiteCode~., fun.aggregate=function(x) paste(x, collapse = ", "), value.var="IUCN_CAT")
    names(IUCN_CAT) <- c("SiteCode", "IUCN_CAT")
    GIS_AREA <- dcast(as.data.table(unique(ProtectedSites[,c("SiteCode", "GIS_AREA")])), SiteCode~., fun.aggregate=function(x) paste(x, collapse = ", "), value.var="GIS_AREA")
    names(GIS_AREA) <- c("SiteCode", "GIS_AREA")
    DESIG <- dcast(as.data.table(unique(ProtectedSites[,c("SiteCode", "DESIG")])), SiteCode~., fun.aggregate=function(x) paste(x, collapse = ", "), value.var="DESIG")
    names(DESIG) <- c("SiteCode", "DESIG")
    
    if(length(unique(nrow(STATUS_YR), nrow(WDPAID), nrow(IUCN_CAT), nrow(GIS_AREA), nrow(DESIG)))!=1){stop("Your reducing to one PA per site has gone wrong")}
    
    ProtectedSites <- list(STATUS_YR, WDPAID, IUCN_CAT, GIS_AREA, DESIG) %>% reduce(left_join, by = "SiteCode")
    
    ### Add back into count data
    ProtectedCountData <- merge(CountData, ProtectedSites, by="SiteCode")
    rm(CountData)
    ProtectedCountData$SiteSpec <- paste0(ProtectedCountData$SiteCode, ".", ProtectedCountData$Species)
    
    #Subset to cases where the PA was designated within the survey period (so that we can assess trends both before and after)
    ProtectedCountDataMinYear <- dcast(ProtectedCountData, SiteSpec~., min, value.var="Year")
    names(ProtectedCountDataMinYear) <- c("SiteSpec", "MinYear")
    ProtectedCountDataMaxYear <- dcast(ProtectedCountData, SiteSpec~., max, value.var="Year")
    names(ProtectedCountDataMaxYear) <- c("SiteSpec", "MaxYear")
    ProtectedCountData <- list(ProtectedCountData, ProtectedCountDataMinYear, ProtectedCountDataMaxYear) %>% reduce(left_join, by = "SiteSpec")
    ProtectedCountData$FallsWithin <- ifelse(ProtectedCountData$STATUS_YR > ProtectedCountData$MinYear & ProtectedCountData$STATUS_YR < ProtectedCountData$MaxYear, "TRUE", "FALSE")
    ProtectedCountData <- subset(ProtectedCountData, FallsWithin=="TRUE")
    
    save(ProtectedCountData, file=paste0(DatabaseFP, "ProtectedCountData.RData"))
    
    #### Clean CBC species and rearrange data ####
    #Here we remove any CBC species that do not show a relationship between number of counts and effort (see methods). And also do a bit of general cleaning
    
    #Load protected and unprotected data, assign CI (1 for protected data, 0 for unprotected data)
    load(file=paste0(DatabaseFP, "ProtectedCountData.RData"))
    ProtectedCountData$CI <- 1
    load(file=paste0(DatabaseFP, "UnprotectedCountData.RData"))
    UnprotectedCountData$CI <- 0
    CountData <- rbindlist(list(ProtectedCountData, UnprotectedCountData), fill=TRUE)
    rm(ProtectedCountData, UnprotectedCountData)
    
    #Remove CBC species without a relationship between effort and number of individuals

    #CleanCBC
    CountDataCBC <- subset(CountData, Dataset=="CBC")
    
    CountvsLogHours <- rbindlist(pbmclapply(unique(CountDataCBC$Species), function(i){ #Cycle through each SpecPop
      #print(i)
      Species <- subset(CountDataCBC, Species==i)[,c("SiteSpec", "Count", "Hours")] #Subset dataframe to specpop
      SpeciesvsError <- tryCatch(glm.nb(Count~log(Hours), link=log, data=Species), error=function(e){
        print(e)
        return(NULL)
      }) #Run a GLM of counts vs. log of hours
      if(is.null(SpeciesvsError)){return(NULL)}
      SpeciesvsErrorSignificant <- as.data.frame(ifelse(summary(SpeciesvsError)$coeff[-1,4]< 0.05, coef(SpeciesvsError, silent=TRUE)[[2]], NA)) #Return the coefficient IF the p value is less than 0.05, otherwise NA
      if(nrow(SpeciesvsErrorSignificant)==0){return(NULL)}
      SpeciesvsErrorSignificant$Species <- i
      names(SpeciesvsErrorSignificant) <- c("Slope", "Species")
      return(SpeciesvsErrorSignificant)
    }, mc.cores=ncores))
    
    CountvsLogHours <- CountvsLogHours[!is.na(CountvsLogHours$Slope),] #Remove the NAs (i.e. insignificant relations)
    CountvsLogHours <- subset(CountvsLogHours, Slope>0) #Remove the negative slopes
    #144 species left down from 310. Ouch. 
    
    CountDataCBC <- unique(CountDataCBC[!CountDataCBC$Species %in% CountvsLogHours$Species,][,c("SiteSpec")]) #Identify those to remove
    CountData <- CountData[!CountData$SiteSpec %in% CountDataCBC,] #Remove from main dataset


    #Standardise so that the same species are in protected and unprotected sites 
    Protected <- subset(CountData, CI==1)
    Unprotected <- subset(CountData, CI==0)
    rm(CountData)
    Unprotected <- Unprotected[Unprotected$Species %in% Protected$Species,]
    Protected <- Protected[Protected$Species %in% Unprotected$Species,]
    CountData <- rbind(Protected, Unprotected)
    save(CountData, file=paste0(DatabaseFP, "CountDataCleaned.RData"))
    
    #### Clean covariates and remove collinear variables ####
    load(file=paste0(DatabaseFP, "CountDataCleaned.RData"))
    
    #Make a subset for testing coding (this takes a while with the full dataset)
    Test <- FALSE
    if(Test){
      RandomSpecies <- unique(CountData$Species)[sample(1:length(unique(CountData$Species)), 100)]
      CountData <- CountData[CountData$Species %in% RandomSpecies,]
      CountData <- subset(CountData, Database=="Waterbird") 
    }
    
    CountData$AnthRound <- as.factor(round_any(CountData$Anthrome, 10)) #Get the higher order value for anthrome (e.g. 21, 22, 23 and 24 are all anthrome values for various types of village, this rounds to 20 which just means "village")
    CountData$SiteSpecYear <- paste0(CountData$SiteSpec, ".", CountData$Year) #Create a unique ID for each count entry (count of a particular species at a particular site in a particular year)
    GovMean <- dcast(CountData, SiteCode ~ ., mean, value.var="Gov", na.rm = TRUE) #Get a mean governance value for each site based on the years it was surveyed
    names(GovMean) <- c("SiteCode", "GovMean")
    #Add the mean governance values to the main dataset (check nothing is lost)
    Check <- nrow(CountData)
    CountData <- merge(CountData, GovMean, by="SiteCode")
    if(Check!=nrow(CountData)){stop("You've lost rows!")}
    CountData <- as.data.frame(CountData)
    rm(GovMean)
    
    #We're now going to reudce the predcitor covariates to only those that don't covary using the ColVIF function that I define at the start of the script. This function doesn't cope if any covariates have NA values, so first we'll pull out any entries where any covariate is an NA, and add a flag these entries as such. (Nb, this is a very small proportion of the full dataset)
    #First we pull out of the dataset all the variables that are continuous (And therefore that we need to assess for covariation)
    ContinuousVariableNames <- c("SiteSpecYear", "MeanAnnualTemp", "MinAnnualTemp", "MaxAnnualTemp", "TotalAnnualPrecip", "Nitr", "Phos", "ConvRangeland", "Cropland",
                                 "Grazing", "IrrigatedNoRice", "IrrigatedRice", "Pastute", "Rangeland", "RainfedNoRice", "RainfedRice", "TotalIrrigatedLand", "TotalRainfedLand",
                                 "TotalRice", "PopulationCount", "PopulationDensity", "RuralPopulationCount", "TotalBuiltUpArea", "UrbanPopulationCount", "Travel", "Slope", "Pasture", "TotalBuiltupArea")
    ContinuousVariableNames <- c(ContinuousVariableNames, c("Gov", "MeanSeasonalTemp",
                                                              "MinSeasonalTemp", "MaxSeasonalTemp", "TotalSeasonalPrecip", "SWOccurrence", "SWExtent"))
    
    ContinuousVariables <- CountData[,names(CountData) %in% ContinuousVariableNames] #Get just continuous variables
    #Make "ContinuousVariables" only include cases where all the continuous variables are *not* NA
    ContinuousVariables <- ContinuousVariables[complete.cases(ContinuousVariables),]
    ContinuousVariables <- ContinuousVariables[match(unique(ContinuousVariables$SiteSpecYear), ContinuousVariables$SiteSpecYear),]
    ContinuousVariables$NAFlag <- 0
    if(nrow(ContinuousVariables)!=length(unique(ContinuousVariables$SiteSpecYear))){stop("There's duplicates!")}
    
    #Flag any entries in the main dataset that do have NA values in the continuous variables
    CountData$NAFlag <- ifelse(CountData$SiteSpecYear %in% ContinuousVariables$SiteSpecYear, 1,0)
    
    #Find collinearity in the non NA entries
    CorrPred <- ContinuousVariables
    CorrPred[,c("SiteSpecYear", "NAFlag")] <- NULL
    
    RemoveCollinearity <- CollVIF(CorrPred, ContinuousVariableNames, scale=FALSE) #Run the function CollVIFto cut down to uncorrelated variables
    
    VariablesForMatchingByYear <- RemoveCollinearity[[2]] #Get the names of the reduced set of variables
    
    #Finish prepping
    ContinuousVariablesToRemove <- ContinuousVariableNames[!ContinuousVariableNames %in% c("SiteSpecYear", VariablesForMatchingByYear)]
    CountData <- CountData[,!colnames(CountData) %in% ContinuousVariablesToRemove] #Remove the variables that were covarying
    
    #Separate into protected and unprotected sets
    ProtectedCountData <- subset(CountData, CI==1)
    UnprotectedCountData <- subset(CountData, CI==0)
    rm(CountData)
    
    #Combine unprotected data with unprotected initially (we pulled these out earlier - this is data where a site becomes protected in say, 2000, but that means we can use the counts in years before that as unprotected counts)
    load(file=paste0(DatabaseFP, "UnprotectedInitially.RData"))
    UnprotectedInitially$MinDesigYear <- NULL
    UnprotectedInitially$CI <- 0
    UnprotectedInitially$AnthRound <- as.factor(round_any(UnprotectedInitially$Anthrome, 10))
    GovMean <- dcast(UnprotectedInitially, SiteCode ~ ., mean, value.var="Gov", na.rm = TRUE) #Get a mean governance value for each Site
    names(GovMean) <- c("SiteCode", "GovMean")
    UnprotectedInitially <- as.data.frame(merge(UnprotectedInitially, GovMean, by="SiteCode"))
    UnprotectedInitially <- UnprotectedInitially[,!colnames(UnprotectedInitially) %in% ContinuousVariablesToRemove]
    UnprotectedInitially$UnprotInitialFlag <- 1 #Add a flag just so we know which cases these are
    
    UnprotectedCountData <- rbindlist(list(UnprotectedCountData, UnprotectedInitially), fill=TRUE)
    rm(UnprotectedInitially)
    
    #Add BA to protected, that is, add a column that is 0 for years before the protected area is designated, and 1 after
    ProtectedCountData$BA <- ifelse(ProtectedCountData$Year>ProtectedCountData$STATUS_YR, 1, 0)
    
    #Save the protected data, the unprotected data, and the list of variable names to be used for matching
    save(ProtectedCountData, file=paste0(DatabaseFP, "ProtectedCountDataCollinearRemoved.RData"))
    save(UnprotectedCountData, file=paste0(DatabaseFP, "UnprotectedCountDataCollinearRemoved.RData"))
    VariablesForMatchingByYear <- append(VariablesForMatchingByYear[!VariablesForMatchingByYear %in% "Gov"], "GovMean")
    save(VariablesForMatchingByYear, file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))
  }
  #### BACI Matching ####
  #First, a function to do final data preparation for matching. This is the section where we start cutting the data to the required parameters for the particular LHC run. Such as cutting the data to the right number of surveyed years before and after protection.
  #It's also where we calculate the mahalanobis distance matrix to use in matching
  PrepData <- function(ZeroThresh, TotalYearsBuffer, MeasuredYearsBuffer, ParallelP, Scen, Database){
    InputFP <- paste0(DatabaseFP, "BACI/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    
    if(file.exists(paste0(InputFP, "MDYearList.RData"))&file.exists(paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))&file.exists(paste0(InputFP, "UnprotectedCountsReadyForMatching.RData"))){
      return(NULL)
    }
    
    #Load the data
    load(file=paste0(DatabaseFP, "ProtectedCountDataCollinearRemoved.RData"))
    load(file=paste0(DatabaseFP, "UnprotectedCountDataCollinearRemoved.RData"))
    load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))
    UnprotectedCountData$CI <- 0

    MaxStatusYr <- max(ProtectedCountData$Year)-MeasuredYearsBuffer
    
    #Get only cases where we have covariate data in the before years (i.e. pre Max Status Year); we lost a tiny amount of data here where the covariate data I used didn't extend to 2017, but we require the covariate data for matching so have to cut these cases
    UnprotectedCountData <- as.data.frame(UnprotectedCountData)
    UnprotectedCountData$SiteSpecYear <- paste0(UnprotectedCountData$SiteCode, ".", UnprotectedCountData$Species, ".", UnprotectedCountData$Year)
    UnprotectedCountsCovs <- subset(UnprotectedCountData[,colnames(UnprotectedCountData) %in% c("SiteSpecYear","Year", VariablesForMatchingByYear)], Year<=MaxStatusYr)
    UnprotectedCountsRemove <- UnprotectedCountsCovs[!complete.cases(UnprotectedCountsCovs),]
    nrow(UnprotectedCountsRemove)/nrow(UnprotectedCountsCovs) #This is to check the proportion of data lost, is always tiny. 
    UnprotectedCountData <- UnprotectedCountData[!UnprotectedCountData$SiteSpecYear %in% UnprotectedCountsRemove$SiteSpecYear,]
    
    ProtectedCountsCovs <- subset(ProtectedCountData[,colnames(ProtectedCountData) %in% c("SiteSpecYear","Year", VariablesForMatchingByYear)], Year<=MaxStatusYr)
    ProtectedCountsRemove <- ProtectedCountsCovs[!complete.cases(ProtectedCountsCovs),]
    nrow(ProtectedCountsRemove)/nrow(ProtectedCountsCovs)
    ProtectedCountData <- ProtectedCountData[!ProtectedCountData$SiteSpecYear %in% ProtectedCountsRemove$SiteSpecYear,]
    
    #Subset protected data to cases with the appropriate number of years before and after based on that particular LHC run
    ProtectedCountData <- subset(ProtectedCountData, Year<= (STATUS_YR+TotalYearsBuffer) & Year>(STATUS_YR-TotalYearsBuffer))
    ProtectedCountsMeasured <- dcast(as.data.table(ProtectedCountData), SiteSpec~BA, length, value.var="Count")
    ProtectedCountsMeasured <- ProtectedCountsMeasured[ProtectedCountsMeasured$'0'>=MeasuredYearsBuffer & ProtectedCountsMeasured$'1'>=MeasuredYearsBuffer,]
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$SiteSpec %in% ProtectedCountsMeasured$SiteSpec,]
    
    #Calculate before slopes of protected sites. This is used in the BA dataset in the results code file
    ProtectedCountsBefore <- subset(ProtectedCountData, BA==0)
    BeforeSlopes <- CalculateBeforeSlopes(ProtectedCountsBefore, ZeroThresh, parallelise=TRUE, actualvalues = FALSE, PThresh = 0.05, TotalYearsBuffer)
    if(nrow(BeforeSlopes)!=length(unique(ProtectedCountsBefore$SiteSpec))){stop("SiteSpec lost!")}
    ProtectedCountData <- merge(ProtectedCountData, BeforeSlopes, by=c("SiteSpec", "CI"))
    
    #Create mahalanobis distance matrices for each year of protected area designations (see Extended Data Figure 6b).

    #Get same species in protected and unprotected data
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$Species %in% UnprotectedCountData$Species,]
    UnprotectedCountData <- UnprotectedCountData[UnprotectedCountData$Species %in% ProtectedCountData$Species,]
    
    #Get the years we have at least 2 PAs designated (otherwise the mahalanobis distance function breaks)
    SitesPerDesigYear <- dcast(as.data.table(unique(ProtectedCountData[,c("SiteCode", "STATUS_YR")])), STATUS_YR~., length, value.var="SiteCode")
    SitesPerDesigYear <- subset(SitesPerDesigYear, .>2) #Subset to years with at least 2 sites
    AllYears <- SitesPerDesigYear$STATUS_YR
    
    print("begin mahal calculations")
    
    #Run through all the years where PAs have been designated, and create a mahalanobis distance matrix for each one
    MDYearList <- pblapply(AllYears, function(YEAR){
      print(YEAR)
      ProtectedSitesCovs <- as.data.frame(GetMeanCovs(subset(ProtectedCountData, Year<= YEAR & STATUS_YR==YEAR & Year>(YEAR-TotalYearsBuffer)), VariablesForMatchingByYear)) #On the protected data, use the "Get Mean Covs" function to get one value for each covariate for the x years before the designation year of this loop of the apply (where x is defined by this LHC run) [we don't care about the years after designation cos we only match on before years]
      if(nrow(ProtectedSitesCovs)!=length(unique(ProtectedSitesCovs$SiteCode))){stop("There are duplicates in ProtectedCounts!")}
      
      ProtectedSitesCovs <- ProtectedSitesCovs[,colnames(ProtectedSitesCovs) %in% c("SiteCode", VariablesForMatchingByYear)] #Reduce the protected dataset to just the covariates we're using for matching
      if(nrow(ProtectedSitesCovs)!=nrow(ProtectedSitesCovs[complete.cases(ProtectedSitesCovs),])){stop("There are NAs in Protected!")}
      
      UnprotectedCountsSub <- as.data.table(subset(UnprotectedCountData, Year<= YEAR  & Year>(YEAR-TotalYearsBuffer))) #Reduce the unprotected data to just x years before the designation year
      UnprotectedCountsSub <- UnprotectedCountsSub[UnprotectedCountsSub$SiteSpec %in% subset(dcast(UnprotectedCountsSub, SiteSpec~., length, value.var="Count"), . >= MeasuredYearsBuffer)$SiteSpec,] #Check that there are enough years in this time period (e.g. maybe we're cutting the data to the 10 years before and after designation, but for one population only 3 years in that 10 year period were measured, and we require at least 7)
      UnprotectedSitesCovs <- as.data.frame(GetMeanCovs(UnprotectedCountsSub, VariablesForMatchingByYear))
      if(nrow(UnprotectedSitesCovs)!=length(unique(UnprotectedSitesCovs$SiteCode))){stop("There are duplicates in UnprotectedCounts!")}
      
      UnprotectedSitesCovs <- UnprotectedSitesCovs[,colnames(UnprotectedSitesCovs) %in% c("SiteCode", VariablesForMatchingByYear)] #Reduce the unprotected dataset to just the covariates we're using for matching
      if(nrow(UnprotectedSitesCovs)!=nrow(UnprotectedSitesCovs[complete.cases(UnprotectedSitesCovs),])){stop("There are NAs in Unprotected!")}
      
      #In extremely unusual cases, one of the covariates is zero across the board for the protected data (for instance maybe there's only three PAs designated in this particular year, and all are in places where the "grazing" land use type is zero)
      #These will break the mahalanobis distance function (And are also useless to matching if all are zero) so we identify and remove those covariates
      SumCheck <- colSums(ProtectedSitesCovs[,c(2:ncol(ProtectedSitesCovs))])
      SumCheck <- SumCheck[SumCheck>0]
      ProtectedSitesCovs <- ProtectedSitesCovs[,names(ProtectedSitesCovs) %in% c("SiteCode", names(SumCheck))]
      UnprotectedSitesCovs <- UnprotectedSitesCovs[,names(UnprotectedSitesCovs) %in% c("SiteCode", names(SumCheck))]
      
      AllData <- rbind(ProtectedSitesCovs,UnprotectedSitesCovs) #Bring together protected and unprotected data
      AllData$SiteCode <- NULL
      MD <- as.data.frame(t(mahal(c(rep(1, nrow(ProtectedSitesCovs)), rep(0, nrow(UnprotectedSitesCovs))), AllData))) #Run the mahalanobis distance function (from package DOS)
      if(identical(MD[is.na(unique(MD))], numeric(0))==FALSE){stop("There are NAs in Mahalanobis")}
      
      #MD is a matrix where the rows are the unprotected sites and the columns are the protected sites. Each cell is the mahalanobis distance between the sites. 
      rownames(MD) <- UnprotectedSitesCovs$SiteCode
      colnames(MD) <- ProtectedSitesCovs$SiteCode
      return(MD)
    }) #
    if(min(sapply(MDYearList, ncol))<2){stop("One MD year list only has one protected site")}
    print("save mahal")
    names(MDYearList) <- AllYears
    save(MDYearList, file=paste0(InputFP, "MDYearList.RData")) #Save the list of MD matrices (one for each designation year)

    save(ProtectedCountData, file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData")) #Save the data
    save(UnprotectedCountData, file=paste0(InputFP, "UnprotectedCountsReadyForMatching.RData")) #Save the data
  }
  
  #Next, a function to run matching
  RunMatching <- function(ZeroThresh, TotalYearsBuffer, MeasuredYearsBuffer, ParallelP, SitePairDist, Scen){
    #Create file paths, load in all the data
    InputFP <- paste0(DatabaseFP, "BACI/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    
    load(file=paste0(InputFP, "MDYearList.RData"))
    load(file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))
    load(file=paste0(InputFP, "UnprotectedCountsReadyForMatching.RData"))
    load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))
    
    #Filter protected sites to only years that are included in MDYearList
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$STATUS_YR %in% names(MDYearList),]
    
    #Create empty MD matrix of all sites (we'll populate it in the matching function; this is essentially creating the matrix from Extended Data Figure 6c, but empty as of yet)
    ProtectedSiteNames <- unique(ProtectedCountData$SiteCode)
    UnprotectedSiteNames <- unique(UnprotectedCountData$SiteCode)
    MD <- matrix(nrow = length(UnprotectedSiteNames), ncol = length(ProtectedSiteNames))
    rownames(MD) <- UnprotectedSiteNames
    colnames(MD) <- ProtectedSiteNames
    
    #Create a blank version of the output data (this comes in handy later)
    MatchNames <- c("MatchID", "MD", "Match", unique(c(names(ProtectedCountData), names(UnprotectedCountData))))
    
    dir.create(file.path(paste0(InputFP, "SpeciesMatch/")), showWarnings = FALSE)
    
    Spec <- "Anas platyrhynchos" #Just here for testing
    #Now, run through matching species by species
    Matching <- pblapply(unique(ProtectedCountData$Species), function(Spec){
      if(file.exists(paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"))){
        return(NULL)
      }
      
      print(Spec)
      write.csv(NULL, paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"))
      
      #Reduce data to the relevant species
      ProtectedCountsSpecies <- subset(ProtectedCountData, Species==Spec)
      UnprotectedCountsSpecies <- subset(UnprotectedCountData, Species==Spec)
      
      #If the species occurs at only one protected site, discard it (write out the blank "MatchNames" dataframe)
      if(length(unique(ProtectedCountsSpecies$SiteCode))<2){
        write.csv(setNames(data.frame(matrix(ncol = length(MatchNames), nrow = 0)), MatchNames), paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"), row.names = FALSE)
        return(NULL)   
      }
      
      #Bring together protected and unprotected datasets
      CountsSpecies <- data.table::rbindlist(list(as.data.frame(ProtectedCountsSpecies), as.data.frame(UnprotectedCountsSpecies)), use.names=TRUE, fill=TRUE)
      CountsSpecies$Match <- "Not Matched"
      
      #Reduce the empty mahalanobis distance matrix to just the sites where this species occurs
      MatchMatrix <- MD[,colnames(MD) %in% unique(ProtectedCountsSpecies$SiteCode)]
      MatchMatrix <- MatchMatrix[rownames(MatchMatrix) %in% unique(UnprotectedCountsSpecies$SiteCode),]
      
      #If there's only one unprotected site the dataframe goes weird so we just add a dummy column to keep it in the right format (I'm sure there's a smarter way to do this..)
      if(class(MatchMatrix)[[1]]=="logical"){
        MatchMatrix <- as.data.frame(rbind(MatchMatrix, rep("Dummy", length(MatchMatrix))))
        rownames(MatchMatrix) <- c(unique(UnprotectedCountsSpecies$SiteCode), "Dummy")
      }
      MatchMatrix <- as.data.frame(MatchMatrix)

      #Now this really massive loop is to remove all the cases where exact matches remove potential unprotected/protected site matches (This is the step at Extended Data Figure 6d)
      #We do this by looping through each protected site in turn, and knock out any unprotected sites that are not exact matches by marking them as "NoMatch" in the empty MD matrix
      #This loop also then populations the MD matrix with the actual mahalanobis distances (Calculated at the end of the PrepData function just above)
      MatchMatrix2 <- pbmclapply(colnames(MatchMatrix), function(i){
        print(i)
        
        #Prep the data
        MatchColumn <- as.data.frame(MatchMatrix[,i]) #Get the column of the protected site plus all the potential unprotected sites
        MatchColumn$UnprotectedSites <- rownames(MatchMatrix) #Add a column for each of the unprotected sites
        names(MatchColumn) <- c(i, "UnprotectedSites")
        MatchColumn <- subset(MatchColumn, UnprotectedSites!="Dummy") #Remove the dummy column if it's there
        
        StatYr <- unique(subset(ProtectedCountsSpecies, SiteCode==i)$STATUS_YR) #Identify the status year of the protected site, i
        
        ProtectedCountsSpeciesSite <- subset(ProtectedCountsSpecies, SiteCode==i) #Now, get the data for protected site i (from "ProtectedCountsSpecies" - we subset this earlier)
        UnprotectedCountsSpeciesSite <- subset(UnprotectedCountsSpecies, Year<= (StatYr+TotalYearsBuffer) & Year>(StatYr-TotalYearsBuffer)) #Then get the data for all the unprotected sites (from "UnrotectedCountsSpecies" - we subset this earlier)
        if(nrow(UnprotectedCountsSpeciesSite)==0){ #If the species doesn't occur at any of those unprotected sites then the protected site doesn't have any matches, so mark all the unprotected sites as "NoMatch"
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        UnprotectedCountsSpeciesSite$BA <- ifelse(UnprotectedCountsSpeciesSite$Year>StatYr, 1, 0) #Add a BA column to the unprotected sites (we can do this now because we're working off only 1 protected area, which we know the designation year of. Make the unprotected entries 0 in the years before the designation year and 1 after)
        UnprotectedMeasuredYears <- dcast(as.data.table(UnprotectedCountsSpeciesSite), SiteCode~BA, length, value.var="Count") #Count the number of entries in the before and after years
        UnprotectedMeasuredYears <- UnprotectedMeasuredYears[UnprotectedMeasuredYears$'0'>=MeasuredYearsBuffer & UnprotectedMeasuredYears$'1'>=MeasuredYearsBuffer,] #Check there are enough measured years before and after
        UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% UnprotectedMeasuredYears$SiteCode,] #Remove any unprotected sites that don't have enough measured years
        
        if(nrow(UnprotectedCountsSpeciesSite)==0){ #If we're out of unprotected sites, mark the protected site as having no matches
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        #Conduct exact matching
        ProtCovs <- GetMeanCovs(subset(ProtectedCountsSpeciesSite, BA==0), VariablesForMatchingByYear) #Get one value for each covariate in the before years of the protected data
        UnProtCovs <- GetMeanCovs(subset(UnprotectedCountsSpeciesSite, BA==0), VariablesForMatchingByYear) #Ditto unprotected
        
        Anth <- as.character(ProtCovs$AnthRound) #Protected anthrome
        Region <- as.character(ProtCovs$GeoRegion) #Protected region
        MigStat <- unique(as.character(ProtectedCountsSpeciesSite$MigStatus)) #Migratory status of species at protected site
        UnProtCovsSub <- subset(UnProtCovs, AnthRound==Anth & GeoRegion==Region) #Subset unprotected sites to correct anthrome and region
        UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% UnProtCovsSub$SiteCode,] #Subset to correct anthrome and region
        UnprotectedCountsSpeciesSite <- subset(UnprotectedCountsSpeciesSite, MigStatus==MigStat) #Subset unprotected sites to correct migratory status
        if(nrow(UnprotectedCountsSpeciesSite)==0){ #If we're out of unprotected sites, mark the protected site as having no matches
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        if(length(unique(UnprotectedCountsSpeciesSite$SiteCode))!=length(unique(UnprotectedCountsSpeciesSite$SiteSpec))){stop("Something's gone wrong, unprotected sitecode ! = sitespec")}
        
        #Get parallel slopes. This is to satisfy the parallel trends assumption (see paper). We do this using a function I defined in "Functions"
        ParallelSites <- ParallelSlopes(ProtectedCountsSpeciesSite, UnprotectedCountsSpeciesSite, ParallelP, ZeroThresh, parallelise=FALSE, TotalYearsBuffer)
        
        #If there are no unprotected sites that fulfil the parallel trends assumption, mark the protected site as having no matches
        if(nrow(ParallelSites)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        ParallelSites$SiteCode <- str_split_fixed(ParallelSites$SiteSpec, "[.]",3)[,1]
        
        #Calculate the physical distance between the sites, so that we can exclude any that are more than x kms away (where x is defined by the relevant LHC run)
        PhysDist <- as.vector(rdist.earth(as.matrix(unique(ProtectedCountsSpeciesSite[,c("Longitude", "Latitude")])), as.matrix(unique(UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% ParallelSites$SiteCode, c("Longitude", "Latitude")])), miles=FALSE))
        PhysDistSites <- unique(UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% ParallelSites$SiteCode, c("Longitude", "Latitude", "SiteCode")])
        PhysDistSites$PhysDist <- PhysDist
        ParallelSites <- merge(ParallelSites, PhysDistSites[,c("SiteCode", "PhysDist")], by="SiteCode")
        ParallelSites <- subset(ParallelSites, PhysDist<= SitePairDist) #Sub to only sites that are less than x distance away
        
        if(nrow(ParallelSites)==0){ #If we're out of unprotected sites, mark the protected site as having no matches
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        FilteredSites <- ParallelSites$SiteCode #Get the names of all the unprotected sites left
        
        ### Get the mahalanobis distances and add to the MD matrix
        DistanceValues <- as.data.frame(MDYearList[[paste0(StatYr)]]) #Get the mahalanobis distance matrix for the relevant designation year, name it "Distance Values"
        DistanceValues$UnprotectedSites <- rownames(DistanceValues)
        
        #Subset distance values to the relevant sites
        DistanceValues <- DistanceValues[rownames(DistanceValues) %in% FilteredSites,c(i, "UnprotectedSites")]
        MatchColumn[,i] <- "NoMatch"
        
        MatchColumn <- rbind(MatchColumn[!MatchColumn$UnprotectedSites %in% FilteredSites,], DistanceValues)
        if(nrow(MatchColumn)!=nrow(MatchMatrix)){stop("differeing row numbers in match column and match matrix")}
        MatchColumn$UnprotectedSites <- factor(MatchColumn$UnprotectedSites, levels=rownames(MatchMatrix))
        MatchColumn <- MatchColumn[order(MatchColumn$UnprotectedSites),]
        
        if(length(MatchColumn[,i])!=length(MatchColumn[,i][complete.cases(MatchColumn[,i])])){stop("Still NAs in matchcolumn")}
        return(MatchColumn) #Return the column for protected site i, with the values for all the unprtoected sites filled in (either "NoMatch" of the mahalanobis distance)
      }, mc.cores=ncores) # 
      MatchMatrix3 <- join_all(MatchMatrix2, by='UnprotectedSites', type='full') #Bring all the columns returned from MatchMatrix 2 together into one big matrix (this completes step d of  Extended Data Figure 6)
      
      #Do a few checks to make sure everything worked right
      rownames(MatchMatrix3) <- MatchMatrix3$UnprotectedSites
      if(ncol(MatchMatrix3)!=ncol(MatchMatrix)+1){stop("Columns have been lost!")}
      if(nrow(MatchMatrix3[complete.cases(MatchMatrix3),])!=nrow(MatchMatrix3)){stop("There are NAs in MatchMatrix!")}
      MatchMatrix3 <- MatchMatrix3[,c("UnprotectedSites", names(MatchMatrix3)[!names(MatchMatrix3) %in% "UnprotectedSites"])]
      
      #Run the greedy algorithm to get the optimal match for each site (Extended Data Figure 6 step e)
      print("start greedy")
      Greedy <- pbmclapply(1:1000, function(iteration){ #Run 1000 iterations
        MatchMatrix4 <- MatchMatrix3
        MatchedSites <- data.frame("Treatment" = colnames(MatchMatrix4), "Control" = NA, "MD" = NA)
        #The first column of matchmatrix4 is the names of the unprotected sites, so now we make a vector called Sample that is the values from 2 to the number of columns (corresponding to all the protected sites in MatchMatrix), randomised
        #This is not an elegant way to achieve that but it works
        if(ncol(MatchMatrix4)==2){ 
          Sample <- 2
        } else {
          Sample <- sample(c(2:ncol(MatchMatrix4)))
        }
        #Run through each value in sample (i.e. each protected site in the random order of sample) and assign an unprotected site for each. Once that unprotected site is assigned, knock it out, and go to the next value in sample
        for(x in Sample){
          #paste(x)
          if(min(MatchMatrix4[,x])=="NoMatch"){ #If there's only NoMatch options for that protected site, skip it
            next()
          } else { #Otherwise identify the minimum value, and then mark that unprotected site as "NoMatch" for the rest of the protected sites
            Treatment <- rownames(MatchMatrix4[MatchMatrix4[,x] == min(MatchMatrix4[,x]),])
            if(length(Treatment)>1){Treatment <- Treatment[[1]]}
            MatchedSites[x,2] <- Treatment
            MatchedSites[x,3] <- min(MatchMatrix4[,x])
            MatchMatrix4[Treatment, ] <- "NoMatch"
          }
        }
        MatchedSites$MD <- as.numeric(MatchedSites$MD)
        return(MatchedSites)
      }, mc.cores=(ncores-4))
      print("calculate global distance")
      
      #For each of the 1000 iterations of the greedy algorithm, calculate the global (i.e. total) mahalanobis distance (end of ED Fig 6 step e)
      GlobalDistance <- sapply(Greedy, function(x){
        if(is.null(tryCatch(x[complete.cases(x),],  error=function(e){NULL}))){print(c(head(x), class(x)))}
        Greedyx <- x[complete.cases(x),]
        return(sum(Greedyx$MD))
      })
      Matched <- Greedy[[match(min(GlobalDistance),GlobalDistance)]] #Choose the greedy algorithm with the smallest global distance
      Matched <- Matched[!is.na(Matched$Control),]
      
      #If there's no matched protected sites for this species, write out an empty dataframe for this species
      if(nrow(Matched)==0){
        write.csv(setNames(data.frame(matrix(ncol = length(MatchNames), nrow = 0)), MatchNames), paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"), row.names = FALSE)
        return(NULL) 
      }
      
      #Do some final cleaning to create a dataframe 
      Matched$MatchID <- 1:nrow(Matched)
      Matched$Treatment <- as.character(Matched$Treatment)
      Matched$Control <- as.character(Matched$Control)
      Matched <- melt(as.data.table(Matched), id.vars=c("MatchID", "MD"), variable.name = "CI", value.name="SiteCode")
      Matched$CI <- ifelse(Matched$CI=="Treatment",1,0)
      
      #Add the counts in. Now we have all the counts for the protected sites and unprotected sites, and we've linked the matched pairs together via the "MatchID" value (e.g. Match ID would be 1 for a protected and unprotected site pair)
      MatchedCounts <- merge(Matched, CountsSpecies, by=c("SiteCode", "CI"))
      if(length(unique(MatchedCounts$SiteCode))!=length(unique(Matched$SiteCode))){stop("Merging matched with counts has lost sites")}
      MatchedCounts$Match <- "Matched"
      write.csv(MatchedCounts, paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"), row.names = FALSE) #Write out all the counts
      return(MatchedCounts)
    })
  }
  
  #Create the folder/file path for matching
  InputFP <- paste0(DatabaseFP, "BACI/")
  dir.create(file.path(InputFP), showWarnings = FALSE)
  print(paste0(Scen, "BeginBACI"))
  
  #Run the prep data function
  PrepData(ZeroThresh=ZeroThresh, TotalYearsBuffer=TotalYearsBuffer, MeasuredYearsBuffer=MeasuredYearsBuffer, ParallelP=ParallelP, Scen=Scen, Database=Database)
  
  #Run the matching function
  RunMatching(ZeroThresh=ZeroThresh, TotalYearsBuffer=TotalYearsBuffer, MeasuredYearsBuffer=MeasuredYearsBuffer, ParallelP=ParallelP, SitePairDist=SitePairDist, Scen=Scen)
  
  #### CI Matching ####
  #These functions are near identical to the BACI matching functions. I have only commented sections that have changed
  
  PrepDataCI <- function(ZeroThresh, TotalYearsBuffer, MeasuredYearsBuffer){
    InputFP <- paste0(DatabaseFP, "CI/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    
    if(file.exists(paste0(InputFP, "MDYearList.RData")) & file.exists(paste0(InputFP, "ProtectedCountsReadyForMatching.RData")) & file.exists(paste0(InputFP, "UnprotectedCountsReadyForMatching.RData"))){
      return(NULL)
    }
    
    load(file=paste0(DatabaseFP, "ProtectedCountDataCollinearRemoved.RData"))
    load(file=paste0(DatabaseFP, "UnprotectedCountDataCollinearRemoved.RData"))
    load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))
    
    MaxStatusYr <- max(ProtectedCountData$Year)-MeasuredYearsBuffer
    
    #Get only cases where we have complete covariate data in the before years (i.e. pre Max Status Year)
    UnprotectedCountData <- as.data.frame(UnprotectedCountData)
    UnprotectedCountData$SiteSpecYear <- paste0(UnprotectedCountData$SiteCode, ".", UnprotectedCountData$Species, ".", UnprotectedCountData$Year)
    UnprotectedCountsCovs <- UnprotectedCountData[,colnames(UnprotectedCountData) %in% c("SiteSpecYear","Year", VariablesForMatchingByYear)]
    UnprotectedCountsRemove <- UnprotectedCountsCovs[!complete.cases(UnprotectedCountsCovs),]
    nrow(UnprotectedCountsRemove)/nrow(UnprotectedCountsCovs)
    UnprotectedCountData <- UnprotectedCountData[!UnprotectedCountData$SiteSpecYear %in% UnprotectedCountsRemove$SiteSpecYear,]
    
    ProtectedCountsCovs <- ProtectedCountData[,colnames(ProtectedCountData) %in% c("SiteSpecYear","Year", VariablesForMatchingByYear)]
    ProtectedCountsRemove <- ProtectedCountsCovs[!complete.cases(ProtectedCountsCovs),]
    nrow(ProtectedCountsRemove)/nrow(ProtectedCountsCovs)
    ProtectedCountData <- ProtectedCountData[!ProtectedCountData$SiteSpecYear %in% ProtectedCountsRemove$SiteSpecYear,]
    
    #Subset protected data to cases with the appropriate number of years before and after
    ProtectedCountData <- subset(ProtectedCountData, Year<= (STATUS_YR+TotalYearsBuffer) & Year>(STATUS_YR)) #Now we only take the years AFTER designation (As in CI we're ignoring the before period)
    ProtectedCountsMeasured <- dcast(as.data.table(ProtectedCountData), SiteSpec~., length, value.var="Count")
    ProtectedCountsMeasured <- ProtectedCountsMeasured[ProtectedCountsMeasured$.>=MeasuredYearsBuffer,]
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$SiteSpec %in% ProtectedCountsMeasured$SiteSpec,]
    
    #Remove cases of all zeroes
    ProtectedCountDataZeroes <- AllZeroes(ProtectedCountData, ZeroThresh)
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$SiteSpec %in% subset(ProtectedCountDataZeroes, AllZero=="Counts")$SiteSpec,]
    
    #Get same species in protected and unprotected data
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$Species %in% UnprotectedCountData$Species,]
    UnprotectedCountData <- UnprotectedCountData[UnprotectedCountData$Species %in% ProtectedCountData$Species,]
    
    #Get the years we have at least 2 PAs designated
    SitesPerDesigYear <- dcast(as.data.table(unique(ProtectedCountData[,c("SiteCode", "STATUS_YR")])), STATUS_YR~., length, value.var="SiteCode")
    SitesPerDesigYear <- subset(SitesPerDesigYear, .>2) #Subset to years with at least 2 sites
    AllYears <- SitesPerDesigYear$STATUS_YR
    
    MDYearList <- pblapply(AllYears, function(YEAR){
      print(YEAR)
      
      ProtectedSitesCovs <- as.data.frame(GetMeanCovs(subset(ProtectedCountData, STATUS_YR==YEAR), VariablesForMatchingByYear))
      if(nrow(ProtectedSitesCovs)!=length(unique(ProtectedSitesCovs$SiteCode))){stop("There are duplicates in ProtectedCounts!")}
      
      ProtectedSitesCovs <- ProtectedSitesCovs[,colnames(ProtectedSitesCovs) %in% c("SiteCode", VariablesForMatchingByYear)]
      ProtectedSitesCovs <- ProtectedSitesCovs[complete.cases(ProtectedSitesCovs),]
      if(nrow(ProtectedSitesCovs)!=nrow(ProtectedSitesCovs[complete.cases(ProtectedSitesCovs),])){stop("There are NAs in Protected!")}
      if(nrow(ProtectedSitesCovs)<2){return(NULL)}
      
      UnprotectedCountsSub <- as.data.table(subset(UnprotectedCountData, Year>YEAR  & Year <= (YEAR+TotalYearsBuffer))) #Now subset unprotected data to just AFTER years, because it's CI (but to a max of x years after designation, according to the TotalYearsBuffer)
      UnprotectedCountsSub <- UnprotectedCountsSub[UnprotectedCountsSub$SiteSpec %in% subset(dcast(UnprotectedCountsSub, SiteSpec~., length, value.var="Count"), . >= MeasuredYearsBuffer)$SiteSpec,]
      if(nrow(UnprotectedCountsSub)==0){return(NULL)}
      UnprotectedSitesCovs <- as.data.frame(GetMeanCovs(UnprotectedCountsSub, VariablesForMatchingByYear))
      if(nrow(UnprotectedSitesCovs)!=length(unique(UnprotectedSitesCovs$SiteCode))){stop("There are duplicates in UnprotectedCounts!")}
      
      UnprotectedSitesCovs <- UnprotectedSitesCovs[,colnames(UnprotectedSitesCovs) %in% c("SiteCode", VariablesForMatchingByYear)]
      UnprotectedSitesCovs <- UnprotectedSitesCovs[complete.cases(UnprotectedSitesCovs),]
      if(nrow(UnprotectedSitesCovs)!=nrow(UnprotectedSitesCovs[complete.cases(UnprotectedSitesCovs),])){stop("There are NAs in Unprotected!")}
      if(nrow(UnprotectedSitesCovs)<2){return(NULL)}
      #Add in a remove zeroes function
      SumCheck <- colSums(ProtectedSitesCovs[,c(2:ncol(ProtectedSitesCovs))])
      SumCheck <- SumCheck[SumCheck>0]
      ProtectedSitesCovs <- ProtectedSitesCovs[,names(ProtectedSitesCovs) %in% c("SiteCode", names(SumCheck))]
      UnprotectedSitesCovs <- UnprotectedSitesCovs[,names(UnprotectedSitesCovs) %in% c("SiteCode", names(SumCheck))]
      
      AllData <- rbind(ProtectedSitesCovs,UnprotectedSitesCovs)
      AllData$SiteCode <- NULL
      MD <- as.data.frame(t(mahal(c(rep(1, nrow(ProtectedSitesCovs)), rep(0, nrow(UnprotectedSitesCovs))), AllData)))
      if(identical(MD[is.na(unique(MD))], numeric(0))==FALSE){stop("There are NAs in Mahalanobis")}
      rownames(MD) <- UnprotectedSitesCovs$SiteCode
      colnames(MD) <- ProtectedSitesCovs$SiteCode
      return(MD)
    }) #
    names(MDYearList) <- AllYears
    MDYearList <- MDYearList[!sapply(MDYearList, is.null)]
    if(min(sapply(MDYearList, ncol))<2){stop("One MD year list only has one protected site")}
    save(MDYearList, file=paste0(InputFP, "MDYearList.RData"))
    save(ProtectedCountData, file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))
    save(UnprotectedCountData, file=paste0(InputFP, "UnprotectedCountsReadyForMatching.RData"))
  }
  CIMatching <- function(ZeroThresh, TotalYearsBuffer, MeasuredYearsBuffer){
    InputFP <- paste0(DatabaseFP, "CI/")
    dir.create(file.path(InputFP), showWarnings = FALSE)
    
    load(file=paste0(InputFP, "MDYearList.RData"))
    load(file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))
    load(file=paste0(InputFP, "UnprotectedCountsReadyForMatching.RData"))
    load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))
    
    #Filter protected sites to only years that are included in MDYearList
    ProtectedCountData <- ProtectedCountData[ProtectedCountData$STATUS_YR %in% names(MDYearList),]
    
    #Create empty MD matrix of all sites
    ProtectedSiteNames <- unique(ProtectedCountData$SiteCode)
    UnprotectedSiteNames <- unique(UnprotectedCountData$SiteCode)
    
    MD <- matrix(nrow = length(UnprotectedSiteNames), ncol = length(ProtectedSiteNames))
    rownames(MD) <- UnprotectedSiteNames
    colnames(MD) <- ProtectedSiteNames
    
    #Create a dummy variable:
    MatchNames <- c("MatchID", "MD", "Match", unique(c(names(ProtectedCountData), names(UnprotectedCountData))))
    
    dir.create(file.path(paste0(InputFP, "SpeciesMatch/")), showWarnings = FALSE)
    
    Spec <- "Podiceps auritus"
    Spec <- "Anas acuta"
    
    Matching <- pblapply(unique(ProtectedCountData$Species), function(Spec){
      if(file.exists(paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"))){
        return(NULL)
      }
      
      print(Spec)
      write.csv(NULL, paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"))
      
      #Reduce treatment to only sites with buffer years before and after
      ProtectedCountsSpecies <- subset(ProtectedCountData, Species==Spec & BA==1) #Subset the protected data to only where BA == 1 (i.e. the After Years)
      UnprotectedCountsSpecies <- subset(UnprotectedCountData, Species==Spec)
      
      ProtectedCountsSpeciesZeroes <- AllZeroes(ProtectedCountsSpecies, ZeroThresh)
      ProtectedCountsSpecies <- ProtectedCountsSpecies[ProtectedCountsSpecies$SiteSpec %in% subset(ProtectedCountsSpeciesZeroes, AllZero=="Counts")$SiteSpec,]
      
      if(length(unique(ProtectedCountsSpecies$SiteCode))<2){
        write.csv(setNames(data.frame(matrix(ncol = length(MatchNames), nrow = 0)), MatchNames), paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"), row.names = FALSE)
        return(NULL)   
      }
      
      CountsSpecies <- data.table::rbindlist(list(ProtectedCountsSpecies, UnprotectedCountsSpecies), use.names=TRUE, fill=TRUE)
      CountsSpecies$Match <- "Not Matched"
      
      MatchMatrix <- MD[,colnames(MD) %in% unique(ProtectedCountsSpecies$SiteCode)]
      MatchMatrix <- MatchMatrix[rownames(MatchMatrix) %in% unique(UnprotectedCountsSpecies$SiteCode),]
      if(class(MatchMatrix)[[1]]=="logical"){
        MatchMatrix <- as.data.frame(rbind(MatchMatrix, rep("Dummy", length(MatchMatrix))))
        rownames(MatchMatrix) <- c(unique(UnprotectedCountsSpecies$SiteCode), "Dummy")
        MatchMatrix <- MatchMatrix[1,]
      } else {
        MatchMatrix <- as.data.frame(MatchMatrix)
      }
      
      #Now filter out exact matched sites (by giving the difference in propensity scores/MD a value of 10)
      MatchMatrix2 <- pbmclapply(colnames(MatchMatrix), function(i){
        print(i)
        MatchColumn <- as.data.frame(MatchMatrix[,i])
        MatchColumn$UnprotectedSites <- rownames(MatchMatrix)
        names(MatchColumn) <- c(i, "UnprotectedSites")
        MatchColumn <- subset(MatchColumn, UnprotectedSites!="Dummy")
        
        StatYr <- unique(subset(ProtectedCountsSpecies, SiteCode==i)$STATUS_YR)
        
        ProtectedCountsSpeciesSite <- subset(ProtectedCountsSpecies, SiteCode==i)
        
        UnprotectedCountsSpeciesSite <- subset(UnprotectedCountsSpecies, Year<= (StatYr+TotalYearsBuffer) & Year>(StatYr)) #Again, only the after years
        if(nrow(UnprotectedCountsSpeciesSite)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        UnprotectedMeasuredYears <- dcast(as.data.table(UnprotectedCountsSpeciesSite), SiteCode~., length, value.var="Count")
        UnprotectedMeasuredYears <- UnprotectedMeasuredYears[UnprotectedMeasuredYears$.>=MeasuredYearsBuffer,]
        UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% UnprotectedMeasuredYears$SiteCode,]
        
        if(nrow(UnprotectedCountsSpeciesSite)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        #Remove cases of all zeroes (we cannot compare these for CI, see methods)
        UnprotectedCountsSpeciesSiteZeroes <- AllZeroes(UnprotectedCountsSpeciesSite, ZeroThresh)
        UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteSpec %in% subset(UnprotectedCountsSpeciesSiteZeroes, AllZero=="Counts")$SiteSpec,]
        
        if(nrow(UnprotectedCountsSpeciesSite)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        #Conduct exact matching
        ProtCovs <- as.data.frame(GetMeanCovs(subset(ProtectedCountsSpeciesSite), VariablesForMatchingByYear))
        ProtCovsCont <- ProtCovs[,colnames(ProtCovs) %in% c("SiteCode", VariablesForMatchingByYear)]
        ProtCovsCont <- ProtCovsCont[complete.cases(ProtCovsCont),]
        if(nrow(ProtCovsCont)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        UnProtCovs <- GetMeanCovs(UnprotectedCountsSpeciesSite, VariablesForMatchingByYear)
        
        Anth <- as.character(ProtCovs$AnthRound) #Protected anthrome
        Region <- as.character(ProtCovs$GeoRegion) #Protected region
        MigStat <- unique(as.character(ProtectedCountsSpeciesSite$MigStatus))
        UnProtCovsSub <- as.data.frame(subset(UnProtCovs, AnthRound==Anth & GeoRegion==Region)) #Subset to correct antrhome and region
        UnProtCovsSub <- UnProtCovsSub[,colnames(UnProtCovsSub) %in% c("SiteCode", VariablesForMatchingByYear)]
        UnProtCovsSub <- UnProtCovsSub[complete.cases(UnProtCovsSub),]
        
        UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% UnProtCovsSub$SiteCode,] #Subset to correct antrhome and region
        UnprotectedCountsSpeciesSite <- subset(UnprotectedCountsSpeciesSite, MigStatus==MigStat)
        if(nrow(UnprotectedCountsSpeciesSite)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        if(length(unique(UnprotectedCountsSpeciesSite$SiteCode))!=length(unique(UnprotectedCountsSpeciesSite$SiteSpec))){stop("Something's gone wrong, unprotected sitecode ! = sitespec")}
        
        #Calculate the physical distance between the sites
        FilteredSites <- unique(UnprotectedCountsSpeciesSite[, c("SiteCode", "Longitude", "Latitude")])
        FilteredSites$PhysDist <- as.vector(rdist.earth(as.matrix(unique(ProtectedCountsSpeciesSite[,c("Longitude", "Latitude")])), as.matrix(unique(FilteredSites[, c("Longitude", "Latitude")])), miles=FALSE))
        FilteredSites <- subset(FilteredSites, PhysDist<= SitePairDist)
        
        if(nrow(FilteredSites)==0){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        
        FilteredSites <- unique(FilteredSites$SiteCode)
        
        #Note there's now no section to identify parallel trends because we can't do this with CI
        
        ### Get the mahalanobis distances and add to matrix
        if(is.null(MDYearList[[paste0(StatYr)]])){
          MatchColumn[,i] <- "NoMatch"
          return(MatchColumn)
        }
        DistanceValues <- as.data.frame(MDYearList[[paste0(StatYr)]])
        DistanceValues$UnprotectedSites <- rownames(DistanceValues)
        
        #Subset distance values to the relevant sites
        DistanceValues <- DistanceValues[rownames(DistanceValues) %in% FilteredSites,c(i, "UnprotectedSites")]
        MatchColumn[,i] <- "NoMatch"
        
        MatchColumn <- rbind(MatchColumn[!MatchColumn$UnprotectedSites %in% FilteredSites,], DistanceValues)
        if(nrow(MatchColumn)!=nrow(MatchMatrix)){stop("differing row numbers in match column and match matrix")}
        MatchColumn$UnprotectedSites <- factor(MatchColumn$UnprotectedSites, levels=rownames(MatchMatrix))
        MatchColumn <- MatchColumn[order(MatchColumn$UnprotectedSites),]
        
        if(length(MatchColumn[,i])!=length(MatchColumn[,i][complete.cases(MatchColumn[,i])])){stop("Still NAs in matchcolumn")}
        return(MatchColumn)
      }, mc.cores=ncores) # 
      MatchMatrix3 <- join_all(MatchMatrix2, by='UnprotectedSites', type='full')
      rownames(MatchMatrix3) <- MatchMatrix3$UnprotectedSites
      if(ncol(MatchMatrix3)!=ncol(MatchMatrix)+1){stop("Columns have been lost!")}
      if(nrow(MatchMatrix3[complete.cases(MatchMatrix3),])!=nrow(MatchMatrix3)){stop("There are NAs in MatchMatrix!")}
      MatchMatrix3 <- MatchMatrix3[,c("UnprotectedSites", names(MatchMatrix3)[!names(MatchMatrix3) %in% "UnprotectedSites"])]
      
      #Find the optimal match for each protected site
      Greedy <- pbmclapply(1:1000, function(iteration){
        MatchMatrix4 <- MatchMatrix3
        MatchedSites <- data.frame("Treatment" = colnames(MatchMatrix4), "Control" = NA, "MD" = NA)
        Sample <- sample(c(2:ncol(MatchMatrix4)))
        for(x in Sample){
          paste(x)
          if(min(MatchMatrix4[,x])=="NoMatch"){
            next()
          } else {
            Treatment <- rownames(MatchMatrix4[MatchMatrix4[,x] == min(MatchMatrix4[,x]),])
            if(length(Treatment)>1){Treatment <- Treatment[1]}
            MatchedSites[x,2] <- Treatment
            MatchedSites[x,3] <- min(MatchMatrix4[,x])
            MatchMatrix4[Treatment, ] <- "NoMatch"
          }
        }
        MatchedSites$MD <- as.numeric(MatchedSites$MD)
        return(MatchedSites)
      }, mc.cores=(ncores-4))
      
      GlobalDistance <- sapply(Greedy, function(x){
        Greedyx <- x[complete.cases(x),]
        return(sum(Greedyx$MD))
      })
      Matched <- Greedy[[match(min(GlobalDistance),GlobalDistance)]]
      Matched <- Matched[!is.na(Matched$Control),]
      if(nrow(Matched)==0){
        write.csv(setNames(data.frame(matrix(ncol = length(MatchNames), nrow = 0)), MatchNames), paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"), row.names = FALSE)
        return(NULL) 
      }
      Matched$MatchID <- 1:nrow(Matched)
      Matched$Treatment <- as.character(Matched$Treatment)
      Matched$Control <- as.character(Matched$Control)
      Matched <- melt(as.data.table(Matched), id.vars=c("MatchID", "MD"), variable.name = "CI", value.name="SiteCode")
      Matched$CI <- ifelse(Matched$CI=="Treatment",1,0)
      MatchedCounts <- merge(Matched, CountsSpecies, by=c("SiteCode", "CI"))
      if(length(unique(MatchedCounts$SiteCode))!=length(unique(Matched$SiteCode))){stop("Merging matched with counts has lost sites")}
      MatchedCounts$Match <- "Matched"
      write.csv(MatchedCounts, paste0(InputFP, "SpeciesMatch/Matched", Spec,".csv"), row.names = FALSE)
      return(MatchedCounts)
    })
  }
  
  print(Scen)
  InputFP <- paste0(DatabaseFP, "CI/")
  dir.create(file.path(InputFP), showWarnings = FALSE)
  
  #Run the functions
  PrepDataCI(ZeroThresh=ZeroThresh, TotalYearsBuffer=TotalYearsBuffer, MeasuredYearsBuffer=MeasuredYearsBuffer)
  CIMatching(ZeroThresh=ZeroThresh, TotalYearsBuffer=TotalYearsBuffer, MeasuredYearsBuffer=MeasuredYearsBuffer)
  
  #### Mark as done ####
  write.csv(NULL, paste0(ResultsFP, "Scen", Scen, "MarkerFINISHED.csv"))
}

#### The final stage of matching is conducted in the next script, once the cluster run has finished ####