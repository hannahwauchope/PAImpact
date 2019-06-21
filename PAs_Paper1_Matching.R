#BACI PA analysis
#8456
options(repos = c(CRAN = "http://cran.rstudio.com"))

cluster <- FALSE

if(cluster==TRUE){
  library(reshape2, lib.loc="/home/hsw34/my-R-libs/")
  library(data.table, lib.loc="/home/hsw34/my-R-libs/")
  library(pbmcapply, lib.loc="/home/hsw34/my-R-libs/")
  library(tidyverse, lib.loc="/home/hsw34/my-R-libs/")
  library(plyr, lib.loc="/home/hsw34/my-R-libs/")
  library(stringr, lib.loc="/home/hsw34/my-R-libs/")
  library(pbapply, lib.loc="/home/hsw34/my-R-libs/")
  library(StatMatch, lib.loc="/home/hsw34/my-R-libs/")
  library(rgdal, lib.loc="/home/hsw34/my-R-libs/")
  library(rgeos, lib.loc="/home/hsw34/my-R-libs/")
  library(sp, lib.loc="/home/hsw34/my-R-libs/")
  library(raster, lib.loc="/home/hsw34/my-R-libs/")
  library(MASS, lib.loc="/home/hsw34/my-R-libs/")

  ResultsFP <- "/rds/user/hsw34/hpc-work/results/c3/"
  DataFP <- "/rds/user/hsw34/hpc-work/data/"
  ncores <- 32 
} else {
  library(ggmap)
  library(stringr)
  library(pbapply)
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
  variables <- variables[!variables %in% c("SW", "Travel", "Slope", "GovMean")]
  Spec6_ShrinkVariables <- Spec6_Shrink[colnames(Spec6_Shrink) %in% variables]
  if(nrow(Spec6_ShrinkVariables)==1){
    Spec6_ShrinkVariables <- as.data.frame(t(as.data.frame(apply(Spec6_ShrinkVariables, 2, function(x) as.numeric(as.character(x))))))
    row.names(Spec6_ShrinkVariables) <- NULL
  } else {
    Spec6_ShrinkVariables <- as.data.frame(apply(Spec6_ShrinkVariables, 2, function(x) as.numeric(as.character(x))))
  }
  
  Spec6_ShrinkVariables$SiteCode <- Spec6_Shrink$SiteCode
  staticvariables <- c("ISO3", "Country", "GeoRegion", "GeoSubRegion", "SW", "LOTW", "Travel", "Slope", "GovMean")
  SitesCovCast <- dcast(Spec6_Shrink, SiteCode ~ ., Mode, value.var="Anthrome")
  names(SitesCovCast) <- c("SiteCode", "Anthrome")
  SitesCovCast <- cbind(SitesCovCast, do.call(cbind, pblapply(variables, function(x){
    SitesCovCast <- dcast(Spec6_ShrinkVariables, SiteCode ~ ., mean, value.var=x, na.rm = TRUE)
    names(SitesCovCast) <- c("SiteCode", x)
    return(as.data.frame(SitesCovCast)[2])
  })))
  SitesCovStatic <- unique(Spec6_Shrink[,c("SiteCode", staticvariables)])
  SitesCovCast <- merge(SitesCovCast, SitesCovStatic, by="SiteCode")
  SitesCovCast$AnthRound <- as.factor(round_any(as.numeric(as.character(SitesCovCast$Anthrome)), 10))
  
  return(SitesCovCast)
}
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpecSeason'), function(x) length(unique(Dataset[,c(x)])))
}
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #grey, orange, light blue, green, yellow, dark blue, dark orange, pink

SiteCoordinates <- unique(Spec6_Cov[,c("SiteCode", "Latitude", "Longitude")])

AssessMap <- function(PASites){
  if(!"Treatment" %in% colnames(PASites)){
    PASites$Treatment <- 1
  }
  if(!"Latitude" %in% colnames(PASites)){
    PASites <- merge(PASites, SiteCoordinates, by="SiteCode")
  }
  PASites <- as.data.frame(PASites)
  PASites <- unique(PASites[,c("SiteCode", "Latitude", "Longitude", "Treatment")])
  PASites$Treatment <- factor(as.character(PASites$Treatment), levels=c(1,0))
  PASites$Treatment <- recode_factor(PASites$Treatment, '1' = "Protected", '0' ="Unprotected")
  PASites <- PASites[order(PASites$Treatment, decreasing=TRUE),]
  
  #Making an inset map watch the fuck out
  mapWorld <- borders(database="world", colour="grey70", fill="grey70")
  mapWorld2 <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey70", fill="grey70")
  
  ggplot() + mapWorld2 +
  #geom_polygon(data = mapWorld4, aes(x=mapWorld4$long, y=mapWorld4$lat, map_id=mapWorld4$group)) +
  #geom_polygon(aes(x=CountryDeets$long, y=CountryDeets$lat))+
  geom_point(aes(x=PASites$Longitude, y=PASites$Latitude, colour=PASites$Treatment), size=1.3, shape=19, alpha=0.9)+ #, colour="dodgerblue4" #, colour="#009E73"
  coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
  scale_colour_manual(values=c("#0B4EA2", "olivedrab3"), guide=FALSE)+
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(size = 1, fill = NA, colour="black"),
    legend.justification = "top",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.background=element_rect(fill = "transparent"),
    panel.background = element_rect(fill = "transparent"))
}

#### Extract sites within PAs, merge with all data ####
PP <- readOGR(paste0(DataFP, "/ProtectedPlanet/WDPA_June2017-shapefile-polygons.shp")) #Read in protected area file (full global PA shapefile, downloaded from WDPA)
Countries <- readOGR("/Users/hannahwauchope/Documents/ArcGIS/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp") #Read through a shapefile of the world's countries

##Make a dataframe of sites and years
SitesAndYears <- unique(Spec6_Cov[,c("SiteCode", "Longitude", "Latitude","Year")]) #Get a dataset of just the Sites and the year they were designation
SiteByYear <- dcast(SitesAndYears, SiteCode + Longitude + Latitude ~ ., length, value.var="Year") #Cast to find the number of years surveyed at each site
SiteByYearMin <- dcast(SitesAndYears, SiteCode ~ ., min, value.var="Year") #Find the minimum survey year at each site
SiteByYearMax <- dcast(SitesAndYears, SiteCode ~ ., max, value.var="Year") #Find the maximum survey year at each site
SitesMinMax <- merge(SiteByYearMin, SiteByYearMax, by = "SiteCode") #Bring the data together
names(SitesMinMax) <- c("SiteCode", "MinYear", "MaxYear") #Rename
names(SiteByYear) <- c("SiteCode", "Longitude","Latitude","NumYears") 
SiteByYear <- merge(SiteByYear, SitesMinMax, by="SiteCode") #Bring the data together

sitepoints <- cbind(SiteByYear$Longitude, SiteByYear$Latitude) #Get just the latitude and longitude
sitepoints <- SpatialPointsDataFrame(sitepoints, SiteByYear, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) #Convert to a spatial points dataframe for overlaying

##Overlay with PA shapefile
#Let's split it up for cluster
ThemsTheBreaks <- seq(0,nrow(sitepoints), 100) #Get breaks so we can split the site points into 100 groups
ThemsTheBreaks <- c(ThemsTheBreaks[1:(length(ThemsTheBreaks)-1)], (nrow(sitepoints))) #Add on the length of sitespoints as the last break
sitepointssplit <- lapply(1:(length(ThemsTheBreaks)-1), function(x){ #Split Sitepoints into a list of smaller chunks (using ThemsTheBreaks)
  sitepointssub <- sitepoints[(ThemsTheBreaks[x]+1):(ThemsTheBreaks[x+1]),] 
  return(sitepointssub)
})
PPOverlay <- pbmclapply(1:length(sitepointssplit), function(i){ #Now, loop through that list to overlay protected areas
  if(file.exists(paste0(ResultsFP, "PPOverlays/PPOverlay", "_", i, ".csv"))){ #Check if the cluster has already done it
    return(NULL)
  }
  if(file.exists(paste0(ResultsFP, "PPOverlays/PPOverlay", "_", i, "_CHECK.csv"))){ #Check if the cluster is working on it
    return(NULL)
  }
  write.csv(NULL, paste0(ResultsFP, "PPOverlays/PPOverlay", "_", i, "_CHECK.csv"), row.names = FALSE) #If not, mark that it's working on it now
  PPOverSites <- over(sitepointssplit[[i]], PP, returnList = TRUE) #Overlay the points with protected planet. We must do this one by one to get all PAs that each point intersects with. It returns a list with data for each PA
  if(length(PPOverSites)==0){ #If the point doesn't intersect, then just write out a null file
    write.csv(NULL, paste0(ResultsFP, "PPOverlays/PPOverlay", "_", i, ".csv"), row.names = FALSE)
  }
  sitepointssplitdf <- as.data.frame(sitepointssplit[[i]]) #Make a dataframe of the points
  GetSiteData <- rbindlist(lapply(1:length(PPOverSites), function(x){ #Extact the PA data
    if(nrow(PPOverSites[[x]])==0){
      return(NULL)
    }
    OverSites <- PPOverSites[[x]] 
    OverSites$SiteCode <- as.character(sitepointssplitdf[x,c("SiteCode")]) #Add the site code
    OverSites$Duplicated <- ifelse(nrow(OverSites)>1, 1, 0) #If there is more than one row (i.e. more than one PA) then mark that it's a duplicated entry with more than 1 PA
    return(OverSites)
  }))
  write.csv(GetSiteData, paste0(ResultsFP, "PPOverlays/PPOverlay", "_", i, ".csv"), row.names = FALSE) #Write out the data
}, mc.cores=ncores)

## Calculate distance from each site to a PA to make a buffer for unprotected sites
sitepointsmoll <- spTransform(sitepoints, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")) #Convert the points to an equal area projection
#PPmoll <- spTransform(PP, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")) #Convert the PAs to an equal area projection
#save(PPmoll, file=paste0(DataFP, "/ProtectedPlanet/PPmoll.RData")) #Save cos this takes forever
load(file=paste0(DataFP, "/ProtectedPlanet/PPmoll.RData")) #Load it back in
AllDistances <- pbmclapply(1:nrow(sitepointsmoll), function(x){ #Loop through each site
  if(x<5001){ #Make seperate folders because otherwise there's too many. This is definitely not the most efficient way to do this but ah well. 
    RunFile <- 5000
  } else if(x<10001){
    RunFile <- 10000
  } else if(x<15001){
    RunFile <- 15000
  } else if(x<20001){
    RunFile <- 20000
  } else {
    RunFile <- 25000
  }
  if(file.exists(paste0(ResultsFP, "PointDistances/",RunFile, "/", x, "_Distance.csv"))){ #If the cluster is already working on it, skip
    return(NULL)
  }
  write.csv(NULL, paste0(ResultsFP, "PointDistances/",RunFile, "/", x, "_Distance.csv")) #Otherwise mark as such
  Site <- sitepointsmoll[x,] #Define the point
  AllDistances <- gDistance(Site,PPmoll) #Find distances - just returns distance to nearest WDPA polygon
  MinDist <- data.frame(AllDistances, unique(Site$SiteCode)) #Make a dataframe with distance and site name
  names(MinDist) <- c("PointDist","SiteCode") #Name
  write.csv(MinDist, paste0(ResultsFP, "PointDistances/",RunFile, "/", x, "_Distance.csv"), row.names=FALSE) #Save
}, mc.cores=ncores)

AllDistances <- rbindlist(lapply(c(5000,10000,15000, 20000,25000), function(RunFile){ #Loop through the folders, read in the file, bind together
  print(RunFile)
  Dist <- rbindlist(lapply(list.files(path = paste0(ResultsFP, "PointDistances/",RunFile, "/"), full.names=TRUE), read.csv))
  return(Dist)
}))

NOTPPSitesOriginal <- Spec6_Cov[Spec6_Cov$SiteCode %in% subset(AllDistances, PointDist>1000)$SiteCode,] #Subset our main dataset to sites over 1km from PAs
save(NOTPPSitesOriginal, file=paste0(ResultsFP, "NOTPPSitesOriginal.RData"))

##Clean Protected Sites
PPSites <- as.data.frame(rbindlist(lapply(list.files(path=paste0(ResultsFP, "PPOverlays/"), pattern=paste0("*.csv$"), full.names=TRUE), fread))) #Read in the sites that overlay PAs (note, have to use this rather than just the distances, as this includes data on all the stacked PAs the point overlaps with)
#And by this method there are 5829 protected sites!!! YAY!!!!! AND same by the distance method! YAY!

#Remove UNESCO biosphere reserves & proposed sites
PPSites <- subset(PPSites, DESIG!="UNESCO-MAB Biosphere Reserve")
PPSites <- subset(PPSites, DESIG_ENG!="UNESCO-MAB Biosphere Reserve")
PPSites <- subset(PPSites, STATUS!="Proposed")

#Round everything to an appropriate number
round_df <- function(x, digits) {
  numeric_columns <- sapply(x, class) == 'numeric'
  x[,numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
PPSites <- round_df(PPSites,4)

PPSites <- subset(PPSites, STATUS_YR!=0) #Remove sites without a status year

#Ok and this cleans so that for each site we just have one PA going on (they sometimes have multiple entries and stuff). 
#We take minimum year if status year varies, and if there are multiple entries of minimum year just the first occurence of a particular year, 
#but with all the desigs entered in and a binary value for multiple min year or not. 
DuplicatedSites <- subset(PPSites, Duplicated>0) #We marked them as duplicated when extracting
DuplicatedSites$Duplicated <- NULL #Remove the duplicated column
NotDuplicatedSites <- subset(PPSites, Duplicated==0) #Get the non-duplicates sites
NotDuplicatedSites$Duplicated <- NULL

DuplicatedSites$STATUS_YR <- as.numeric(as.character(DuplicatedSites$STATUS_YR)) #Convert status_yr to numeric
DuplicatedSites$SiteStatYr <- paste0(DuplicatedSites$SiteCode, ".", DuplicatedSites$STATUS_YR) #Make a field that marks both site and status year
DuplicatedSitesCast <- dcast(DuplicatedSites, SiteCode~., min, value.var="STATUS_YR") #In cases where duplicated PAs have different desingation years, take the earliest year
names(DuplicatedSitesCast) <- c("SiteCode", "STATUS_YR") #Rename
DuplicatedSitesCast$SiteStatYr <- paste0(DuplicatedSitesCast$SiteCode, ".", DuplicatedSitesCast$STATUS_YR) #Make a field that marks both site and status year

DuplicatedSitesYear <- DuplicatedSites[DuplicatedSites$SiteStatYr %in% DuplicatedSitesCast$SiteStatYr,] #So this is now a dataframe of the duplicates, but cut down just the minimum year for each duplicate

#But sometimes we have multiple PAs designated in the same year (or possibly the same PA just with different statuses)
DuplicatedYear <- subset(dcast(DuplicatedSitesYear, SiteCode + SiteStatYr~., length, value.var="SiteCode"), .>1)
DuplicatedYear <- DuplicatedSitesYear[DuplicatedSitesYear$SiteCode %in% DuplicatedYear$SiteCode,] #Get all entries for those sites

DuplicatedYearCast <- aggregate(DESIG~SiteCode, DuplicatedYear, paste) #Gather together all the DESIGs into one column, separated by commas
DuplicatedYearCast$Duplicated <- 2
DuplicatedYearCast2 <- aggregate(DESIG_ENG~SiteCode, DuplicatedYear, paste) #Ditto DESIG_ENG
DuplicatedYearCast <- merge(DuplicatedYearCast, DuplicatedYearCast2, by="SiteCode")
DuplicatedYearCast3 <- aggregate(GIS_AREA~SiteCode, DuplicatedYear, paste) #Ditto GIS_Area
DuplicatedYearCast <- merge(DuplicatedYearCast, DuplicatedYearCast3, by="SiteCode")
DuplicatedYearCast$DESIG <- as.character(DuplicatedYearCast$DESIG)
DuplicatedYearCast$DESIG_ENG <- as.character(DuplicatedYearCast$DESIG_ENG)
DuplicatedYearCast$GIS_AREA <- as.character(DuplicatedYearCast$GIS_AREA)

DuplicatedYear <- DuplicatedSitesYear[duplicated(DuplicatedSitesYear$SiteCode),] #Get the duplicate sites JUST TAKES FIRST ENTRY FOR EACH SO WE HAVE ALL THE DATA FOR DESIG, DESIG_ENG and GIS_AREA BUT NOT FOR EVERYTHING ELSE
DuplicatedYear <- DuplicatedYear[!duplicated(DuplicatedYear$SiteCode),] #Get the duplicate sites JUST TAKES FIRST ENTRY FOR EACH SO WE HAVE ALL THE DATA FOR DESIG ETC BUT NOT FOR EVERYTHING ELSE

DuplicatedYear[,c("DESIG", "DESIG_ENG", "GIS_AREA")] <- NULL
DuplicatedYear <- merge(DuplicatedYear, DuplicatedYearCast, by="SiteCode") #Add in our merged metadata

NOTDuplicatedYear <- DuplicatedSitesYear[!DuplicatedSitesYear$SiteCode %in% DuplicatedYear$SiteCode,] #Get those sites that were in the duplicate list, but had different status years
NOTDuplicatedYear$Duplicated <- 1 #mark them as duplicate = 1 (i.e. duplicates, but only one designated first)
NOTDuplicatedYear <- data.frame(lapply(NOTDuplicatedYear, function(x){ #Convert factors to characters
  if(is.factor(x)==TRUE){
    column <- as.character(x)
    return(column)
  } else {
    return(x)
  }
}), stringsAsFactors=FALSE)

NOTDuplicatedYear <- rbind(NOTDuplicatedYear, DuplicatedYear) #Combine the duplicates

NotDuplicatedSites <- data.frame(lapply(NotDuplicatedSites, function(x){ #For those with no duplicates at all, convert factors to characters
  if(is.factor(x)==TRUE){
    column <- as.character(x)
    return(column)
  } else {
    return(x)
  }
}), stringsAsFactors=FALSE)
NotDuplicatedSites$Duplicated <- 0
NOTDuplicatedYear$SiteStatYr <- NULL

AllSites <- rbind(NotDuplicatedSites, NOTDuplicatedYear) #And combine non duplicate with cleaned duplicates
# 0 = no duplicates, 1 = duplicates but there was one min year, 2 = duplicates by there was more than one min year

## Add to bird data ##
#Add to bird data
PA_Counts <- merge(AllSites, Spec6_Cov, by=c("SiteCode")) #Merge with the count data
PA_Counts[,c("MinYear", "MaxYear")] <- NULL #Remove the min year and max year (so these were no longer necessary to calculate, but leaving in as they do no harm)
PA_Counts$SiteSpec <- paste0(PA_Counts$Species, ".", PA_Counts$SiteCode) #Make a "SiteSpec" combined value, so we can track actual counts for species, not just sites

PA_Counts_MinYearSS <- dcast(PA_Counts, SiteSpec~., min, value.var="Year") #Cast to get minimum year for each sitespec
names(PA_Counts_MinYearSS) <- c("SiteSpec", "MinYearSS") #Rename
PA_Counts_MaxYearSS <- dcast(PA_Counts, SiteSpec~., max, value.var="Year") #Cast to get maximum year
names(PA_Counts_MaxYearSS) <- c("SiteSpec", "MaxYearSS")
PA_Counts <- merge(PA_Counts, PA_Counts_MinYearSS, by="SiteSpec") #Bring together
PA_Counts <- merge(PA_Counts, PA_Counts_MaxYearSS, by="SiteSpec")
PA_Counts$MaxYearSS <- as.numeric(PA_Counts$MaxYearSS) #Convert to numeric
PA_Counts$MinYearSS <- as.numeric(PA_Counts$MinYearSS)
PA_Counts$STATUS_YR <- as.numeric(as.character(PA_Counts$STATUS_YR))

PA_Counts$FallsWithin <- ifelse(PA_Counts$STATUS_YR > PA_Counts$MinYearSS & PA_Counts$STATUS_YR < PA_Counts$MaxYearSS, "TRUE", "FALSE") #See whether the PA designation year falls within the min and max year

Birds_fallsw <- subset(PA_Counts, FallsWithin=="TRUE") #Reduce to those that do

#Get actual measured years
Birds_fallsw$BeforeAfterStat <- ifelse(Birds_fallsw$Year<=Birds_fallsw$STATUS_YR, "BeforeStat", "AfterStat") #Get the number of counted years before and after, not just  the absolute min and max
Birds_fallsw_countedyears <- unique(Birds_fallsw[,c("SiteSpec", "STATUS_YR", "Longitude", "Latitude", "MinYearSS", "MaxYearSS", "Year", "BeforeAfterStat")]) #Reduce
Birds_fallsw_countedyearscast <- dcast(Birds_fallsw_countedyears, SiteSpec~BeforeAfterStat, length, value.var="Year") #Cast
Birds_fallsw <- merge(Birds_fallsw, Birds_fallsw_countedyearscast, by="SiteSpec") #Add the values in
PPSitesOriginal <- Birds_fallsw
save(PPSitesOriginal, file=paste0(ResultsFP, "PPSitesOriginal.RData"))

#To Check Points
write.csv(unique(PPSitesOriginal[,c("SiteCode", "Latitude", "Longitude")]), "/Users/hannahwauchope/Desktop/PPSITES.csv", row.names=FALSE) #Write these out to check points if we need to
write.csv(unique(NOTPPSitesOriginal[,c("SiteCode", "Latitude", "Longitude")]), "/Users/hannahwauchope/Desktop/NOTPPSITES.csv", row.names=FALSE)

#### Prepare Site Matching Data ####
load(file=paste0(ResultsFP, "PPSitesOriginal.RData"))
load(file=paste0(ResultsFP, "NOTPPSitesOriginal.RData"))

PPSitesOriginal$STATUS_YR <- as.numeric(PPSitesOriginal$STATUS_YR) #Make the status year numeric
names(PPSitesOriginal)[names(PPSitesOriginal)=="NAME"] <- "PAName" #Change NAME to PANAME to avoid confusion later
PPSitesOriginal$BeforeStatYears <- (PPSitesOriginal$STATUS_YR-PPSitesOriginal$MinYearSS)+1 #Get the number of years before designation. We INCLUDE the year of designation with the before.
PPSitesOriginal$AfterStatYears <- (PPSitesOriginal$MaxYearSS-PPSitesOriginal$STATUS_YR) #Get the number of years after designation
PPSitesOriginal$SiteSpecSeason <- paste0(PPSitesOriginal$SiteSpec, ".", PPSitesOriginal$Season) #add a site spec season
PPSitesOriginal$Treatment <- 1 #Mark that these are treated sites

NOTPPSitesOriginal$SiteSpecSeason <- paste0(NOTPPSitesOriginal$Species, ".", NOTPPSitesOriginal$SiteCode, ".", NOTPPSitesOriginal$Season) #Add a site spec season to unprotected sites
NOTPPSites_MinYearSS <- dcast(NOTPPSitesOriginal, SiteSpecSeason~., min, value.var="Year") #Find min year of unprotected sitespec
names(NOTPPSites_MinYearSS) <- c("SiteSpecSeason", "MinYearSS")
NOTPPSites_MaxYearSS <- dcast(NOTPPSitesOriginal, SiteSpecSeason~., max, value.var="Year") #Find max year
names(NOTPPSites_MaxYearSS) <- c("SiteSpecSeason", "MaxYearSS")
NOTPPSitesOriginal <- merge(NOTPPSitesOriginal, NOTPPSites_MinYearSS, by="SiteSpecSeason") #Bring together
NOTPPSitesOriginal <- merge(NOTPPSitesOriginal, NOTPPSites_MaxYearSS, by="SiteSpecSeason")
NOTPPSitesOriginal$Treatment <- 0 #Mark that these are untreated sites

AllSites <- as.data.frame(data.table::rbindlist(list(PPSitesOriginal, NOTPPSitesOriginal), fill=TRUE)) #Bring data together
TheNumbers(subset(AllSites, Treatment==1)) #Count number of sites, species, and populations in each
TheNumbers(subset(AllSites, Treatment==0))

#Just dec to feb, to make sure the covariates are appropriate at all times (we lose very little)
AllSites <- subset(AllSites, Season=="DectoFeb")
TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#Remove all zero counts
RemoveZeros <- dcast(AllSites, SiteSpecSeason~., sum, value.var="Count")
RemoveZeros <- subset(RemoveZeros, .==0)
AllSites <- AllSites[!AllSites$SiteSpecSeason %in% RemoveZeros$SiteSpecSeason,]

TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#CleanCBC - See methods, we can only retain CBC species that show a predictable relationship between count and effort
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

AllSitesCBCCleaned <- unique(AllSitesCBC[!AllSitesCBC$SpecPop %in% CountvsHoursLogHours$SpecPop,][,c("SiteSpecSeason")]) #Remove failed CBC species
AllSites <- AllSites[!AllSites$SiteSpecSeason %in% AllSitesCBCCleaned,]
TheNumbers(subset(AllSites, Treatment==1)) #Check how many we lost (not too many)
TheNumbers(subset(AllSites, Treatment==0))

#Remove site species combinations with < 10 years
Buffer <- 5 #We require 5 years before and after designation
NYearBuffer <- 3 #And at least 3 years measured before and after
AllSitesUnique <- unique(AllSites[,c("SiteSpecSeason", "Year")]) #Get the site species and the unique years each has surveys for
AllSitesUniqueMax <- dcast(AllSitesUnique, SiteSpecSeason~., max, value.var="Year") #Find max
AllSitesUniqueMin <- dcast(AllSitesUnique, SiteSpecSeason~., min, value.var="Year") #Find min
AllSitesUniqueMax$Diff <- (AllSitesUniqueMax$.-AllSitesUniqueMin$.)+1 #Find distance spanned
AllSitesUniqueMax <- subset(AllSitesUniqueMax, Diff>=(Buffer*2)) #Keep only those with at least ten years spanned (just to cut it down, we'll do further filters later)
AllSites <- AllSites[AllSites$SiteSpecSeason %in% AllSitesUniqueMax$SiteSpecSeason,] #Cut our dataset to these site species

PPSites <- subset(AllSites[AllSites$Treatment==1,], AfterStatYears>=Buffer & BeforeStatYears>=Buffer & BeforeStat>=NYearBuffer & AfterStat>=NYearBuffer) #For treatment sites, we can cut to exactly number of years before and after designation
NOTPPSites <- subset(AllSites, Treatment==0) #Add to unprotected sites

TheNumbers(PPSites)
TheNumbers(NOTPPSites)

#Standardise so Not PP sites only contain species present in PP sites and vice versa
NOTPPSites <- NOTPPSites[NOTPPSites$Species %in% PPSites$Species,]
PPSites <- PPSites[PPSites$Species %in% NOTPPSites$Species,]

AllSites <- rbind(PPSites, NOTPPSites)
TheNumbers(subset(AllSites, Treatment==1))
TheNumbers(subset(AllSites, Treatment==0))

#Get the data (not including counts)
AllSites$AnthRound <- as.factor(round_any(AllSites$Anthrome, 10)) #Round the anthrome
AllSites$SiteSpecSeasonYear <- paste0(AllSites$SiteSpecSeason, ".", AllSites$Year)
GovMean <- dcast(AllSites, SiteCode ~ ., mean, value.var="Gov", na.rm = TRUE) #Get a mean governance value for each country
names(GovMean) <- c("SiteCode", "GovMean")
AllSites <- merge(AllSites, GovMean, by="SiteCode", all=T)

##Find collinear variables
ContinuousVariables <- AllSites[,c("SiteSpecSeasonYear", "Treatment", "PrecipAnnual", "PrecipSeason", "MeanTempAnnual", "MinTempAnnual", "MaxTempAnnual", "MeanTempSeason", "MinTempSeason", "MaxTempSeason",
                                    "Nitr", "Phos", "grazing", "ir_norice", "ir_rice", "pasture", "rangeland", "rf_norice", "rf_rice", "popd", "uopp", "popc", "GovMean", 
                                    "SW", "Travel", "Slope")] #Get just continuous variables

ContinuousVariables <- ContinuousVariables[complete.cases(ContinuousVariables),] #Remove any cases with NAs for variables (not many)
AllSites <- AllSites[AllSites$SiteSpecSeasonYear %in% ContinuousVariables$SiteSpecSeasonYear,]

#A few species got knocked out in that cull, so remove them from the control sites
SpeciesCheck <- dcast(AllSites, Species~Treatment, length, value.var="SiteCode")
names(SpeciesCheck) <- c("Species", "Cont", "Treat")
SpeciesToKeep <- subset(SpeciesCheck, Treat!=0)
SpeciesToKeep <- subset(SpeciesToKeep, Cont!=0)
AllSites <- AllSites[AllSites$Species %in% SpeciesToKeep$Species,]

TheNumbers(subset(AllSites, Treatment==0))
TheNumbers(subset(AllSites, Treatment==1))

ContinuousVariables <- AllSites[,c("SiteSpecSeasonYear", "Treatment", "PrecipAnnual", "PrecipSeason", "MeanTempAnnual", "MinTempAnnual", "MaxTempAnnual", "MeanTempSeason", "MinTempSeason", "MaxTempSeason",
                                         "Nitr", "Phos", "grazing", "ir_norice", "ir_rice", "pasture", "rangeland", "rf_norice", "rf_rice", "popd", "uopp", "popc", "GovMean", 
                                         "SW", "Travel", "Slope")] #Get just continuous variables again now that we've cleaned

#Find collinearity
CorrPred <- ContinuousVariables
CorrPred[,c("SiteSpecSeasonYear", "Treatment")] <- NULL
SumCheck <- colSums(CorrPred)
SumCheck <- SumCheck[SumCheck>0] #Check none of the variables are zero across the board
SumCheck
CorrPred2 <- CorrPred[,which(names(CorrPred) %in% names(SumCheck))]
vif(CorrPred2) #Check the vif, iteratively remove any variables with a vif greater than 4
CorrPred2[,c("grazing", "MinTempSeason")] <- NULL
vif(CorrPred2)
CorrPred2[,c("popd", "MeanTempSeason")] <- NULL
vif(CorrPred2)
CorrPred2[,c("MeanTempAnnual")] <- NULL
vif(CorrPred2)
CorrPred2[,c("MaxTempSeason")] <- NULL
vif(CorrPred2)

CorrelationPearson <- cor(CorrPred2, method="pearson") #Check correlations with pearson correlation coefficient, remove any >0.70
apply(CorrelationPearson, 2, function(x) x[x>0.70])
CorrPred2$Phos <- NULL
CorrelationPearson <- cor(CorrPred2, method="pearson")
apply(CorrelationPearson, 2, function(x) x[x>0.70])
VariablesForMatchingByYear <- colnames(CorrelationPearson)

#Finish prepping
ContinuousVariablesToRemove <- ContinuousVariables[, !colnames(ContinuousVariables) %in% c("SiteSpecSeasonYear", "Treatment", VariablesForMatchingByYear)] #Remove collinear variables
AllSites <- AllSites[,!colnames(AllSites) %in% colnames(ContinuousVariablesToRemove)]

save(AllSites, file = paste0(ResultsFP, "AllSites.RData")) #Write out data
save(VariablesForMatchingByYear, file=paste0(ResultsFP, "VariablesForMatchingByYear.RData"))

#Finally, calculate trends of protected populations
load(file = paste0(ResultsFP, "AllSites.RData")) #Read in data
PPSites <- subset(AllSites, Treatment==1) #Subset
NOTPPSites <- subset(AllSites, Treatment==0) #Subet
load(file=paste0(ResultsFP, "VariablesForMatchingByYear.RData"))

TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpecSeason'), function(x) length(unique(Dataset[,c(x)])))
} #Little function to calculate number of sites, species and populations in each dataset

### Calculate direction of trends
PPSites$BA <- ifelse(PPSites$Year>PPSites$STATUS_YR, 1, 0) #Make a column with 0s for years before designation and 1 for years after
PPSitesBefore <- subset(PPSites, BA==0) #Subset to just years before designation
BeforeSlopes <- pbmclapply(unique(PPSitesBefore$SiteSpecSeason), function(SS){ #Loop through the populations
  BeforeDat <- subset(PPSitesBefore, SiteSpecSeason==SS)[,c("Count", "Year", "Hours", "Dataset")] #Subset to this population
  if(unique(BeforeDat$Dataset=="CBC")){ #If the dataset is CBC, run a glm.nb including the hour offset
    glmmy <- tryCatch(glm.nb(Count~Year + offset(log(x$Hours)), link=log, data=BeforeDat), error=function(e){NA})
  } else { #Otherwise run a simple glm.nb
    glmmy <- tryCatch(glm.nb(Count~Year, link=log, data=BeforeDat), error=function(e){NA})
  }
  slope <- tryCatch(coef(glmmy, silent=TRUE)[[2]], error=function(e){NA}) #Extract the slope
  significant <- tryCatch(summary(glmmy)$coeff[-1,4], error=function(e){NA}) #And the significance
  return(as.data.frame(cbind(slope, significant, SS)))
}, mc.cores=4)

BeforeSlopes2 <- rbindlist(BeforeSlopes) #Bind together
BeforeSlopes2$significant <- as.character(BeforeSlopes2$significant)
BeforeSlopes2$slope <- as.character(BeforeSlopes2$slope)
NumYears <- dcast(PPSitesBefore, SiteSpecSeason~., length, value.var="Year") #Find number of years each population was surveyed for post designation
BeforeSlopes2 <- merge(BeforeSlopes2, NumYears, by.x="SS", by.y="SiteSpecSeason") #Add number of years to the slopes data

#Now classify as stable (0), positive (1) or negative (-1). See paper for jusitification. 
BeforeSlopes2$Class <- ifelse(is.na(BeforeSlopes2$significant), 0, #If the significance is NA make it zero
                              ifelse(BeforeSlopes2$.>6, ifelse(BeforeSlopes2$slope>0,1,-1), #If the significance is not NA and there's more than 6 years, make it the slope
                                     ifelse(BeforeSlopes2$significant>0.05, 0, #If the significance is not NA, it's 6 or less years and the significance is >0.05 make it zero
                                            ifelse(BeforeSlopes2$slope>0, 1, -1)))) #Otherwise make it the slope

BeforeSlopes2 <- BeforeSlopes2[,c("SS", "Class")] #Just keep the class

PPSites <- merge(PPSites, BeforeSlopes2, by.x="SiteSpecSeason", by.y="SS") #Add to PPSites Data
TheNumbers(PPSites)
TheNumbers(NOTPPSites)

PPSitesBA <- PPSites
save(PPSitesBA, file=paste0(ResultsFP, "PPSitesBA.RData")) #This is our full before/after dataset

#Check number of sites at each desig year, as Mahalanobis distance needs at least two protected sites of each designation year. We lose only a handful of sites from this
PPSitesDC <- unique(PPSites[,c("SiteCode", "STATUS_YR")]) #DC = Desig Check
PPSitesDCCast <- dcast(PPSitesDC, STATUS_YR~., length, value.var="SiteCode") #Cast to find number of protected sites of each Status_Year
PPSitesDCCast1 <- subset(PPSitesDCCast, .==1) #Get cases with only one site
PPSites <- PPSites[!PPSites$STATUS_YR %in% PPSitesDCCast1$STATUS_YR,] #Remove these
TheNumbers(PPSites)

#Create MD matrices for each designation year
AllYears <- unique(PPSites$STATUS_YR) #Get all the designation years
MDYearList <- pbmclapply(AllYears, function(YEAR){
  print(YEAR)
  PPSitesCovs <- GetMeanCovs(subset(PPSites, Year<= YEAR & STATUS_YR==YEAR)) #Get the means of covariates in years pre designation (protected sites)
  NOTPPSitesCovs <- GetMeanCovs(subset(NOTPPSites, Year<= YEAR)) #Get the means of covariates in years pre designation (unprotected sites)
  
  CheckForZeros <- rbind(PPSitesCovs, NOTPPSitesCovs) #Check for any covariates where there are zeros across the board (This breaks mahalanobis distance)
  CheckForZeros <- as.data.frame(t(as.data.frame(t(colSums(CheckForZeros[,colnames(CheckForZeros) %in% VariablesForMatchingByYear]))))) #Get sum of each variable
  CheckForZeros$Var <- row.names(CheckForZeros) #Add rownames
  CheckForZeros <- subset(CheckForZeros, V1==0) #Reduce to only cases where the sum is zero
  PPSitesCovs <- PPSitesCovs[!colnames(PPSitesCovs) %in% CheckForZeros$Var] #Remove these variables
  NOTPPSitesCovs <- NOTPPSitesCovs[!colnames(NOTPPSitesCovs) %in% CheckForZeros$Var] #Remove these variables
  
  MD <- as.data.frame(mahalanobis.dist(data.x= NOTPPSitesCovs[,colnames(NOTPPSitesCovs) %in% VariablesForMatchingByYear], data.y=PPSitesCovs[,colnames(PPSitesCovs) %in% VariablesForMatchingByYear])) #Calculate Mahalanobis distance
  rownames(MD) <- NOTPPSitesCovs$SiteCode
  colnames(MD) <- PPSitesCovs$SiteCode
  MDScale <- scale(MD,center=rep(0.5, ncol(MD)), scale=rep(10, ncol(MD))) #Scale values
  return(MDScale)
}, mc.cores=ncores)

names(MDYearList) <- AllYears #Give each matrix the name of a year

#Save
save(MDYearList, file=paste0(ResultsFP, "MDYearList.RData"))
save(PPSites, file=paste0(ResultsFP, "PPSites.RData"))
save(NOTPPSites, file=paste0(ResultsFP, "NOTPPSites.RData"))

#### Conduct Matching #### 
load(file=paste0(ResultsFP, "VariablesForMatchingByYear.RData")) #Load files
load(file=paste0(ResultsFP, "MDYearList.RData"))
load(file=paste0(ResultsFP, "PPSites.RData"))
load(file=paste0(ResultsFP, "NOTPPSites.RData"))
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpecSeason'), function(x) length(unique(Dataset[,c(x)])))
} #This is just a function to output and count # sites, species and site*species combinations

#Create empty MD matrix of all sites
PPSiteNames <- unique(PPSites$SiteCode) #Get the names of all the protected sites
NOTPPSiteNames <- unique(NOTPPSites$SiteCode) #Ditto unprotected sites

MD <- matrix(nrow = length(NOTPPSiteNames), ncol = length(PPSiteNames)) #Create an empty matrix, with protected sites as column names, and unprotected sites as rownames. This will be used to create a distance matrix for each species in the loop below
rownames(MD) <- NOTPPSiteNames
colnames(MD) <- PPSiteNames

Spec <- "Anser albifrons" #Used to test a particular species in the below loop

Matching <- pblapply(unique(PPSites$Species), function(Spec){ #Loop through all the protected species in the dataset
  if(file.exists(paste0(ResultsFP, "MatchFinal/Matched", Spec,"CHECK.csv"))){
    return(NULL)
  } #Check if another cluster is working on this species, if so skip
  write.csv(NULL, paste0(ResultsFP, "MatchFinal/Matched", Spec,"CHECK.csv")) #Mark that this cluster is working on the species
  
  PPSitesSpec <- subset(PPSites, Species==Spec) #Subset the PPSites dataset (with counts at all sites for all protected species) to the relevant species
  if(length(unique(PPSitesSpec$SiteCode))<2){ #If there are less than two protected sites, skip this species (And write out a null file to mark that this is finished)
    write.csv(NULL, paste0(ResultsFP, "MatchFinal/Matched", Spec,".csv"))
    return(NULL)
  }
  
  NOTPPSitesSpec <- subset(NOTPPSites, Species==Spec) #Subset the unprotected counts to the relevant species
  SpecAllData <- data.table::rbindlist(list(PPSitesSpec, NOTPPSitesSpec), use.names=TRUE, fill=TRUE) #Combine protected and unprotected data together, filling unoverlapping columns with NAs 

  MatchMatrix <- as.data.frame(MD[,colnames(MD) %in% PPSitesSpec$SiteCode]) #Subset the empty matrix (Created outside the loop) to only the relevant protected site columns
  MatchMatrix <- MatchMatrix[rownames(MatchMatrix) %in% NOTPPSitesSpec$SiteCode,] #Ditto unprotected site rows 
  
  #Now, we go through each protected site individually, first filtering out any sites that arent exact matches and then populating the MD matrix with the relevant distances from the MDYearList based on that protected site's designation year
  for(i in colnames(MatchMatrix)){ #I.e. for each protected population (because we've subset to both one site and one species)
    #print(i)
    StatYr <- unique(subset(PPSitesSpec, SiteCode==i)$STATUS_YR) #Get the designation year
    ProtCovs <- GetMeanCovs(subset(PPSitesSpec, SiteCode==i & Year<= StatYr)) #Get the mean covariate values for the protected population for all years pre-designation
    
    #We will exact match on these three covariates
    Anth <- as.character(ProtCovs$AnthRound) #Get the anthrome value for predesignation
    Region <- as.character(ProtCovs$GeoRegion) #Get the region value
    SlopeClass <- unique(subset(PPSitesSpec, SiteCode==i)$Class) #Get the "Class" value (i.e. whether the pre-designation trend is positive, negative or stable)
    
    #Then get the unprotected data
    if(nrow(subset(NOTPPSitesSpec, Year<= StatYr))==0){#If there are no unprotected sites with survey years before the designation year, skip.
      MatchMatrix[,i] <- 10
      next
    } 
    UnProtCovs <- GetMeanCovs(subset(NOTPPSitesSpec, Year<= StatYr)) #Get the mean covariate values for predesignation years for unprotected sites
    NotProtSub <- subset(UnProtCovs, AnthRound==Anth & GeoRegion==Region) #Subset these sites to only cases with the same AnthRound and GeoRegion as the protected site
    NotProtSub <- NOTPPSitesSpec[NOTPPSitesSpec$SiteCode %in% NotProtSub$SiteCode,] #Pull in the full count data for the filtered sites
    NotProtSub$BeforeStatYears <- (StatYr-NotProtSub$MinYearSS)+1 #Find how far the earliest survey year is from the designation year
    NotProtSub$AfterStatYears <- (NotProtSub$MaxYearSS-StatYr) #Ditto the latest survey year
    NotProtSub <- subset(NotProtSub, BeforeStatYears>=5 & AfterStatYears>=5) #Subset to cases where we cover at least 5 years before and after
    if(nrow(NotProtSub)==0){ #If no unprotected sites are left, skip
      MatchMatrix[,i] <- 10
      next
    }
    NotProtSub$BA <- ifelse(NotProtSub$Year<=StatYr, 0, 1) #Mark rows as before or after designation
    NotProtSubYears <- dcast(NotProtSub, SiteCode~BA, length, value.var="Year") #Cast to find number of years before and after designation
    names(NotProtSubYears) <- c("SiteCode", "Before", "After") #Rename
    NotProtSubYears <- subset(NotProtSubYears, Before>=3 & After>=3) #Subset to only cases with at least 3 surveyed years before and after designation
    NotProtSub <- NotProtSub[NotProtSub$SiteCode %in% NotProtSubYears$SiteCode,] #Pull in full count data from relevant sites
    
    NotProtSub <- subset(NotProtSub, BA==0) #Subset to just years before designation
    if(nrow(NotProtSub)==0){ #If no unprotected sites are left, skip
      MatchMatrix[,i] <- 10
      next
    }
    
    #Finally we need to exact match on class (i.e. slope direction pre-designation)
    BeforeSlopes <- pbmclapply(unique(NotProtSub$SiteSpecSeason), function(SS){ #Work through each unprotected population
      BeforeDat <- subset(NotProtSub, SiteSpecSeason==SS)[,c("Count", "Year", "Hours", "Dataset")] #Subset to just the relevant data
      if(length(unique(BeforeDat$Count))==1){ #If there are just zero counts in all years, return an NULL so that that site is removed as an option
        if(unique(BeforeDat$Count)==0){
          return(NULL)
        }
      }
      if(unique(BeforeDat$Dataset=="CBC")){ #If the dataset is CBC, run a glm.nb with an offset for hours, otherwise just run a simple glm.nb
        glmmy <- tryCatch(glm.nb(Count~Year + offset(log(x$Hours)), link=log, data=BeforeDat), error=function(e){NA})
      } else {
        glmmy <- tryCatch(glm.nb(Count~Year, link=log, data=BeforeDat), error=function(e){NA})
      }
      slope <- tryCatch(coef(glmmy, silent=TRUE)[[2]], error=function(e){NA}) #Record the slope
      significant <- tryCatch(summary(glmmy)$coeff[-1,4], error=function(e){NA}) #Record the significance level
      return(as.data.frame(cbind(slope, significant, SS))) #Return data
    }, mc.cores=4)
     
    BeforeSlopes2 <- BeforeSlopes[!sapply(BeforeSlopes, is.null)] #Remove Null cases
    if(length(BeforeSlopes2)==0){ #If no unprotected sites are left, skip
      MatchMatrix[,i] <- 10
      next
    }
    if(length(unique(sapply(BeforeSlopes2, ncol)))!=1){
      print(paste0("PANIC BEFORESLOPES2 HAS MORE THAN 1 NCOL ", Spec, i))
    }
    BeforeSlopes2 <- as.data.frame(rbindlist(BeforeSlopes2))
    names(BeforeSlopes2) <- c("slope", "significant", "SiteSpecSeason") #Rename
    BeforeSlopes2$significant <- as.character(BeforeSlopes2$significant)
    BeforeSlopes2$slope <- as.character(BeforeSlopes2$slope)
    NumYears <- dcast(NotProtSub, SiteSpecSeason~., length, value.var="Year") #Cast to find number of years predesignation
    BeforeSlopes2 <- merge(BeforeSlopes2, NumYears, by="SiteSpecSeason") #Add number of years to the slopes data
    BeforeSlopes2$Class <- ifelse(is.na(BeforeSlopes2$significant), 0, #If the significance is NA make it zero
                                  ifelse(BeforeSlopes2$.>6, ifelse(BeforeSlopes2$slope>0,1,-1), #If the significance is not NA and there's more than 6 years, make it the slope
                                         ifelse(BeforeSlopes2$significant>0.05, 0, #If the significance is not NA, it's 6 or less years and the significance is >0.05 make it zero
                                                ifelse(BeforeSlopes2$slope>0, 1, -1)))) #Otherwise make it the slope
    BeforeSlopes3 <- subset(BeforeSlopes2, Class==SlopeClass) #Reduce to only cases that match the slope class of our protected population
    
    if(nrow(BeforeSlopes3)!=0){ #If there is at least one unprotected site left
      print(c(SlopeClass, unique(BeforeSlopes3$Class))) #Print the classes (for code checking)
      FilteredSites <- str_split_fixed(BeforeSlopes3$SiteSpecSeason, "[.]",3)[,2] #Get the site code
    } else {
      FilteredSites <- NULL
    }
    UnProtCovs <- UnProtCovs[UnProtCovs$SiteCode %in% FilteredSites,] #Filter the unprotected covariate dataset down to just the remaining unprotected sites
    
    ### Ok and now get the distances and add to matrix
    DistanceValues <- MDYearList[[paste0(StatYr)]] #Select the right MD (built based on the right Status_Yr)
    DistanceValues <- DistanceValues[rownames(DistanceValues) %in% UnProtCovs$SiteCode,] #Reduce to only the relevant unprotected sites

    if(is.null(nrow(DistanceValues))){ #I.e. if there's only one unprotected site
      k <- paste0(UnProtCovs$SiteCode) #Get the name of that site
      MatchMatrix[k, i] <- DistanceValues[[i]] #For the particular cell of the relevant protected site (i) and the unprotected site (k), assign the distance value
      MatchMatrix[!rownames(MatchMatrix) %in% UnProtCovs$SiteCode,i] <- 10 #For the rest of the unprotected sites in the protected column, assign them a value of 10 (i.e. NULL)
    } else { #If there's more than one unprotected site
      MatchMatrix[c(names(DistanceValues[,i])), i] <- DistanceValues[,i] #Assign the distance values
      MatchMatrix[!rownames(MatchMatrix) %in% UnProtCovs$SiteCode,i] <- 10 #For the rest of the unprotected sites, assign them a value of 10 (i.e. NULL)
    }
    #End loop of protected populations
  }
  
  MatchMatrix$PlaceHolder <- "PlaceHolder" #Because r is stupid and if there's only one column it won't return the rowname of the min value

  Greedy <- pbmclapply(1:1000, function(iteration){ #Run a greedy algorithm to select the best match for each species
    MatchMatrix3 <- MatchMatrix #Create a Match Matrix just for this loop
    MatchedSites <- data.frame("Treatment" = colnames(MatchMatrix3), "Control" = NA, "MD" = NA) #Create a new dataframe which we will population with the matched pairs - a column for protected sites, a column for unprotected sites and a column for the match distance
    Sample <- sample.int(ncol(MatchMatrix3)-1) #Randomly sample values to the number of columns (i.e. protected sites) in match matrix. This is the order the greedy algorithm will work with for this particular iteration. Skip the last column (placeholde)
    for(x in Sample){ #Work through the protected columns in the order defined by "Sample"
      if(min(MatchMatrix3[,x])==10){ #If the minumum distance value in x's column is 10, it means it has no appropriate matches (they've all been blocked out). Skip. 
        next()
      } else {
        Treatment <- rownames(MatchMatrix3[MatchMatrix3[,x] == min(MatchMatrix3[,x]),]) #This is a bit of a weird way of doing it, but it gets the rowname (i.e. unprotected site name) of the minimum distance to the relevant protected column
        if(length(Treatment)>1){Treatment <- Treatment[1]} #If multiple unprotected sites have the same minimum distance, take the first one
        MatchedSites[x,2] <- Treatment #Put the relevant unprotected site in our matched dataframe
        MatchedSites[x,3] <- MatchMatrix3[MatchMatrix3[,x] == min(MatchMatrix3[,x]),x][1] #Add the distance
        MatchMatrix3[Treatment, ] <- 10 #Mark this site as "chosen" so future protected sites cant take it
      }
    }
    return(MatchedSites[complete.cases(MatchedSites),]) #Return all matched cases (i.e. ones without an NA)
  }, mc.cores=ncores)
  
  GlobalDistance <- sapply(Greedy, function(x){return(sum(x$MD))}) #Calculate the global distance (i.e. summed distance) for each greedy run 
  Matched <- Greedy[[match(min(GlobalDistance),GlobalDistance)]] #Take the Greedy run with the minimum global distance. This is our matched pairings
  if(nrow(Matched)==0){ #If there are less than two protected sites, skip this species (And write out a null file to mark that this is finished)
    write.csv(NULL, paste0(ResultsFP, "MatchFinal/Matched", Spec,".csv"))
    return(NULL)
  }
  Matched$MatchID <- 1:nrow(Matched) #Assign each protected/unprotected pair a match ID
  Matched$Treatment <- as.character(Matched$Treatment)
  Matched$Control <- as.character(Matched$Control)
  Matched <- melt(Matched, id.vars=c("MatchID", "MD")) #Melt down the matched dataset so all site codes (protected or unprotected) are in one column
  names(Matched) <- c("MatchID", "MD", "Treatment", "SiteCode") #Rename the columns
  Matched$Treatment <- NULL #Remove the treatment column (it's already defined in SpecAllData)
  Matched <- merge(Matched, SpecAllData, by="SiteCode") #Add the matched values (+ distance, and the MatchID which defines the pairs) to the full dataframe
  write.csv(Matched, paste0(ResultsFP, "MatchFinal/Matched", Spec,".csv"), row.names = FALSE) #Write out the data
  return(Matched)
})

#### FOR MEETING KILL ME LATER ####
PPSitesCountry <- unique(PPSites[,c("SiteCode","Country")])
PPSitesCountry <- dcast(PPSitesCountry, Country~., length, value.var="SiteCode")
PPSitesCountry <- PPSitesCountry[order(PPSitesCountry$.,decreasing = TRUE),]
names(PPSitesCountry) <- c("Country","NProtectedSites")

NOTPPSitesCountry <- unique(NOTPPSites[,c("SiteCode","Country")])
NOTPPSitesCountry <- dcast(NOTPPSitesCountry, Country~., length, value.var="SiteCode")
NOTPPSitesCountry <- NOTPPSitesCountry[order(NOTPPSitesCountry$.,decreasing = TRUE),]
names(NOTPPSitesCountry) <- c("Country","NUnprotectedSites")

Both <- merge(PPSitesCountry, NOTPPSitesCountry, by="Country")
Both <- Both[order(Both$NProtectedSites,decreasing = TRUE),]

#### Assess matching ####
load(file=paste0(ResultsFP, "PPSites.RData")) #Read in the original Protected data
load(file=paste0(ResultsFP, "NOTPPSites.RData")) #And the original unprotected data
load(file=paste0(ResultsFP, "VariablesForMatchingByYear.RData")) #And the variables we used for matching

#Check that all the species matching runs went through - first read in the "Check" files that are empty, but marked that a cluster run had started working on the species
MatchingFinal <- list.files(path=paste0(ResultsFP, "MatchFinal/"), pattern=paste0("*CHECK.csv$"), full.names=TRUE)
MatchingFinalNames <- str_split_fixed(MatchingFinal, "Matched",2)[,2] #Get the species names
MatchingFinalNames <- str_split_fixed(MatchingFinalNames, "CHECK.csv",2)[,1]

MatchingFinal2 <- list.files(path=paste0(ResultsFP, "MatchFinal/"), full.names=TRUE) #Then read in the actual matching files (well actually, just read in all files)
MatchingFinal2 <- MatchingFinal2[!MatchingFinal2 %in% MatchingFinal] #Then reduce to those that are the actual matching files
MatchingFinal2Names <- str_split_fixed(MatchingFinal2, "Matched",2)[,2] #Get the species names
MatchingFinal2Names <- str_split_fixed(MatchingFinal2Names, ".csv",2)[,1]

if(length(MatchingFinal2Names[!MatchingFinal2Names%in%MatchingFinalNames])>0){print("PANIC A SPECIES IS MISSING!")} #If we don't have the same length of check and actual species files, there's an error!

MatchingFinal <- pblapply(MatchingFinal2, function(x){ #If all is ok, read in the actual matching files
  return(fread(x, fill=TRUE))
})

MatchingFinal <- lapply(MatchingFinal, function(x){ if(ncol(x)==1){return(NULL)} else {return(x)}}) #Remove any empty files
MatchingFinal <- MatchingFinal[!sapply(MatchingFinal, is.null)]
MatchingFinal <- as.data.frame(rbindlist(MatchingFinal)) #Bind all species together into one object

#To assess the matching, we need to assess the standardisd difference in means for covariates between protected and unprotected groups. We therefore need to re-get the mean covariate values for all years predesignation
MatchingFinal$MatchID <- as.numeric(MatchingFinal$MatchID) #Make sure the MatchID is numeric
MatchingFinal$SpecMatch <- paste0(MatchingFinal$Species, ".", MatchingFinal$MatchID) #Make a unique Match ID for each species (these values should now only have one protected and unprotected site associated with them)
MatchingCovariates <- rbindlist(pbmclapply(unique(MatchingFinal$SpecMatch), function(SpecMatchID){ #Get the covariate values for predesignation years (now that we know which PA each unprotected site is matched to)
  print(SpecMatchID)
  try <- subset(MatchingFinal, SpecMatch==SpecMatchID) #Get the values for the protected and unprotected site of the relevant MatchID
  StatYr <- unique(try$STATUS_YR)[complete.cases(unique(try$STATUS_YR))] #Define the StatusYear
  ProtCovs <- GetMeanCovs(subset(try, Treatment==1 & Year<= StatYr)) #Get protected covariate means in years pre designation
  ProtCovs$Treatment <- 1 #Mark as protected
  UnProtCovs <- GetMeanCovs(subset(try, Treatment==0 & Year<= StatYr)) #Get unprotected covariate means in years pre designation
  UnProtCovs$Treatment <- 0 #Mark as unprotected
  AllCovs <- rbind(ProtCovs, UnProtCovs) #Bring covariates together
  AllCovs$SpecMatch <- SpecMatchID #Mark the match ID
  return(AllCovs)
}, mc.cores=ncores))

MatchingCovariates$Species <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,1] #Split back so we have the Species name
MatchingCovariates$MatchID <- str_split_fixed(MatchingCovariates$SpecMatch, "[.]",2)[,2] #And the match ID

Distances <- unique(MatchingFinal[,c("SpecMatch", "MD")]) #Get just the distances from the match final dataset
MatchingCovariates <- merge(MatchingCovariates, Distances, by="SpecMatch", all=T) #Add to the matching covariates dataset

#Now we can compare the standardised differences in means
StDiffMean <- function(dataset){
  SpecSites <- dataset
  if(nrow(SpecSites)<=2){
    return(NULL)
  }
  SiteCast <- do.call(cbind,lapply(c(1,0), function(x){ #Get the mean values for each variable in each treatment
    TheVariables <- as.data.frame(subset(SpecSites, Treatment==x)) #Subset to the treatment
    TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear] #Reduce to numeric columns
    TheVariables$Treatment <- NULL #Remove treatment
    TheVariables <- apply(TheVariables, 2, as.numeric) #Make sure they're all numeric
    return(as.data.frame(colMeans(TheVariables, na.rm = TRUE))) #Return the column means
  }))
  names(SiteCast) <- c("Treat", "Cont") #Add appropriate names
  SiteCast$MeanDiff <- abs(SiteCast$Treat - SiteCast$Cont) #Find the difference in means between treatment and control
  SiteCastVar <- do.call(cbind,lapply(c(1,0), function(x){ #Get the sd for each variable in each treatment
    TheVariables <- as.data.frame(subset(SpecSites, Treatment==x)) #Susbet to the treatment
    TheVariables <- TheVariables[,colnames(TheVariables) %in% VariablesForMatchingByYear] #Reduce to numeric columns
    TheVariables$Treatment <- NULL #Remove treatment
    TheVariables <- apply(TheVariables, 2, as.numeric) #Make sure they're all numeric
    return(as.data.frame(colVars(TheVariables, na.rm = TRUE))) #Return the column standard deviations
  }))
  names(SiteCastVar) <- c("Treat", "Cont")
  SiteCast$VarComp <- sqrt((SiteCastVar$Treat + SiteCastVar$Cont)/2) #Get standard deviation
  SiteCast$d <- SiteCast$MeanDiff/SiteCast$VarComp #Get standardised dsitance (mean/sd)
  SiteCast[is.na(SiteCast)] <- 0 #Change any NA's to zeros (occurs in cases where variable is 0 in all cases)
  SiteCast$Cov <- row.names(SiteCast) #Add variable names
  SiteCast <- SiteCast[,c(5:6)] #Just keep SDiM and variable name
  return(SiteCast)
}

CompareMatching <- function(dataset, DiffThresh, MatchType){
  DeleteOldFiles <- c(list.files(paste0(ResultsFP, "Summaries/", MatchType), full.names=TRUE), list.files(paste0(ResultsFP, "MatchData/", MatchType), full.names=TRUE))
  sapply(DeleteOldFiles, unlink)
  Comp <- pbmclapply(unique(dataset$Species), function(i){
    print(i)
    MatchingSub <- subset(dataset, Species==i) #Subset to the relevant species 
    if(nrow(MatchingSub)<=4){ #Kill the species with less than 2 matched pairs
      return(NULL)
    }

    TheDistance <- StDiffMean(MatchingSub) #Calculate the standardised difference in means for each covariate
    Threshold <- subset(TheDistance, d<DiffThresh) #Subset to all covariates that are less than the "difference threshold"
    
    while(nrow(Threshold)!=length(VariablesForMatchingByYear) & nrow(MatchingSub)>4){ #While not all variables have SDiM below the threshold
      #print(Iteration)
      MAX <- max(MatchingSub$MD) #Get the biggest MD distance
      MatchingSub <- subset(MatchingSub, MD!=MAX) #Remove the link that has that distance
      TheDistance <- StDiffMean(MatchingSub) #Recalculate the SDiM
      Threshold <- subset(TheDistance, d<DiffThresh)
    }
    if(nrow(Threshold)==length(VariablesForMatchingByYear) & nrow(MatchingSub)>=4){ #If we reached a point where all covariates were below the threshold
      names(Threshold) <- c("StDiff", "Variable") #rename  to StDiff
      Threshold$Species <- i
      write.csv(Threshold, paste0(ResultsFP, "Summaries/Final/Summary_", i, ".csv"), row.names=FALSE) #Save stdiffs in means
      write.csv(MatchingSub, paste0(ResultsFP, "MatchData/Final/MatchData_", i, ".csv"), row.names=FALSE) #Save actual datasets
      }
  }, mc.cores=ncores)
}
CompareMatching(MatchingCovariates, 0.25, "Final") #Run through and iteratively remove matches until we have the matched set for each species

SummariesFinal <- rbindlist(lapply(list.files(path=paste0(ResultsFP, "Summaries/Final/"), full.names=TRUE), fread)) #Read in that summary material
MatchCovariateMeans <- rbindlist(lapply(list.files(path=paste0(ResultsFP, "MatchData/Final/"), full.names=TRUE), fread))
MatchingFinalCleaned <- MatchingFinal[MatchingFinal$SpecMatch %in% MatchCovariateMeans$SpecMatch,]

#save(MatchingFinalCleaned, file=paste0(ResultsFP, "MatchData/MatchingFinalCleaned.RData"))
#load(file=paste0(ResultsFP, "MatchData/MatchingFinalCleaned.RData"))

TheNumbers(subset(MatchingFinalCleaned, Treatment==1))
TheNumbers(subset(MatchingFinalCleaned, Treatment==0))

load(file=paste0(ResultsFP, "PPSitesBA.RData")) #Read in the original Protected data
TheNumbers(PPSitesBA)
TheNumbers(NOTPPSites)

#Extract species stats for comparison
SpecStats <- read.csv(file=paste0(DataFP, "WaterbirdData_Tatsuya/FullDataSet_Edits/Hannah_Consolidation/SpeciesStats.csv"))
SpecStats$Order <- as.factor(SpecStats$Order)
SpecStats$Order <- recode_factor(SpecStats$Order, "CHARADRIIFORMES" = "Charadriiformes", "ANSERIFORMES" = "Anseriformes", "GRUIFORMES" = "Gruiformes", "CICONIIFORMES" = "Ciconiiformes",
                                 "SULIFORMES" = "Suliformes", "PELECANIFORMES" = "Pelicaniformes", "PROCELLARIIFORMES" = "Procellariiformes", "GAVIIFORMES" = "Gaviformes", 
                                 "PHOENICOPTERIFORMES" = "Phoenicopteriformes", "PODICIPEDIFORMES" = "Podicipediformes")
SpecStats$Genus <- str_split_fixed(SpecStats$Species, "[ ]", 2)[,1]

TaxaData <- function(Dataset){
  AllSitesStats <- unique(Dataset[,c("Species", "SiteCode", "Treatment")])
  AllSitesStats <- merge(AllSitesStats, SpecStats, by="Species")
  AllSitesStatsFamSite <- unique(AllSitesStats[,c("Order", "Family", "SiteCode", "Treatment")])
  FamiliesSite <- dcast(AllSitesStatsFamSite, Order + Family~Treatment, length, value.var="SiteCode")
  AllSitesStatsFamSpec <- unique(AllSitesStats[,c("Order", "Family","Species")])
  FamiliesSpec <- dcast(AllSitesStatsFamSpec, Order + Family~., length, value.var="Species")
  AllSitesStatsFamGen <- unique(AllSitesStats[,c("Order", "Family", "Genus")])
  FamiliesGen <- dcast(AllSitesStatsFamGen, Order + Family~., length, value.var="Genus")
  FamilyData <- merge(FamiliesSite, FamiliesGen, all=TRUE)
  names(FamilyData) <- c("Order", "Family", "UnprotectedSites", "ProtectedSites", "Genera")
  FamilyData <- merge(FamilyData, FamiliesSpec, all=TRUE)
  names(FamilyData) <- c("Order", "Family", "UnprotectedSites", "ProtectedSites", "Genera", "Species")
  NOrder <- length(unique(AllSitesStats$Order))
  NFamily <- length(unique(AllSitesStats$Family))
  NGenus <- length(unique(AllSitesStats$Genus))
  NSpecies <- length(unique(AllSitesStats$Species))
  return(list(FamilyData, NOrder,NFamily,NGenus, NSpecies))
}
AllSites2 <- rbind(PPSitesBA[,c("SiteCode", "Species", "Treatment")], NOTPPSites[,c("SiteCode", "Species", "Treatment")])
AllTaxa <- TaxaData(AllSites2)
AllTaxa
MatchTaxa <- TaxaData(MatchCovariateMeans)
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

#### Maps ####
library(ggplot2)
library(ggmap)
library(rgeos)
library(plyr)
library(lattice)
library(RColorBrewer)
library(pbapply)
library(maps)
library(ggalt)

MapVersion <- "Matched" #"Filtered" "Matched
for(MapVersion in c("Matched", "Filtered")){ #"NotFiltered"
  
  if(MapVersion=="NotFiltered"){
    NOTPPSitesNOTFILTERED <- fread(paste0(ResultsFP, "NOTPPSites.csv"))
    PPSitesNOTFILTERED <- as.data.frame(rbindlist(lapply(list.files(path=paste0(ResultsFP, "PPOverlays/"), pattern=paste0("*.csv$"), full.names=TRUE), fread)))
    PPSitesNOTFILTERED <- Spec6_Cov[Spec6_Cov$SiteCode %in% PPSitesNOTFILTERED$SiteCode,]
    PPSitesNOTFILTERED$Treatment <- 1
    NOTPPSitesNOTFILTERED$Treatment <- 0
    #write.csv(unique(PPSitesNOTFILTERED[,c("SiteCode", "Latitude", "Longitude", "Treatment")]), "/Users/hannahwauchope/Desktop/birdtest/PPNotFiltered.csv", row.names=FALSE)
    #write.csv(unique(NOTPPSitesNOTFILTERED[,c("SiteCode", "Latitude", "Longitude", "Treatment")]), "/Users/hannahwauchope/Desktop/birdtest/NOTPPNotFiltered.csv", row.names=FALSE)
    
    #SiteNamesPPNotFiltered <- SiteNames[SiteNames$SiteCode %in% PPSitesNOTFILTERED$SiteCode,]
    #SiteNamesNOTPPNotFiltered <- SiteNames[SiteNames$SiteCode %in% NOTPPSitesNOTFILTERED$SiteCode,]
    #hey <- SiteNamesNOTPPNotFiltered[grep(" ", SiteNamesNOTPPNotFiltered$SiteName),]
    #write.csv(hey, "/Users/hannahwauchope/Desktop/birdtest/NOTPPNotFiltered.csv", row.names=FALSE)
    
    
    MapSites <- rbind(PPSitesNOTFILTERED[,c("SiteCode", "Latitude", "Longitude", "Treatment")], NOTPPSitesNOTFILTERED[,c("SiteCode", "Latitude", "Longitude", "Treatment")])
  } 
  if(MapVersion=="Filtered"){
    load(file=paste0(ResultsFP, "PPSitesBA.RData"))
    load(file=paste0(ResultsFP, "NOTPPSites.RData"))
    #write.csv(unique(PPSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")]), "/Users/hannahwauchope/Desktop/birdtest/PPFiltered.csv", row.names=FALSE)
    #write.csv(unique(NOTPPSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")]), "/Users/hannahwauchope/Desktop/birdtest/NOTPPFiltered.csv", row.names=FALSE)
    
    MapSites <- rbind(PPSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")], NOTPPSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")])
  }
  if(MapVersion=="Matched"){
    load(file=paste0(ResultsFP, "PPSites.RData"))
    load(file=paste0(ResultsFP, "NOTPPSites.RData"))
    load(file=paste0(ResultsFP, "MatchData/MatchingFinalCleaned.RData"))
    MapSites <- rbind(PPSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")], NOTPPSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")])
    MapSites <- MapSites[MapSites$SiteCode %in% MatchingFinalCleaned$SiteCode,]
  }
  
  PASites <- unique(MapSites[,c("SiteCode", "Latitude", "Longitude", "Treatment")])
  PASites$Treatment <- factor(as.character(PASites$Treatment), levels=c(1,0))
  PASites$Treatment <- recode_factor(PASites$Treatment, '1' = "Protected", '0' ="Unprotected")
  PASites <- PASites[order(PASites$Treatment, decreasing=TRUE),]
  
  #Making an inset map watch the fuck out
  mapWorld <- borders(database="world", colour="grey70", fill="grey70")
  mapWorld2 <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey70", fill="grey70")
  
  mp  <- ggplot() + mapWorld +
    #geom_polygon(data = mapWorld4, aes(x=mapWorld4$long, y=mapWorld4$lat, map_id=mapWorld4$group)) +
    #geom_polygon(aes(x=CountryDeets$long, y=CountryDeets$lat))+
    geom_point(aes(x=PASites$Longitude, y=PASites$Latitude, colour=PASites$Treatment), size=1.8, shape=19, alpha=0.9)+ #, colour="dodgerblue4" #, colour="#009E73"
    coord_proj("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")+
    scale_colour_manual(values=c("#0B4EA2", "olivedrab3"), guide=FALSE)+
    theme(
      aspect.ratio=0.53, 
      panel.grid = element_blank(), 
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0),
      legend.justification = "top",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background=element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"))
  mp
  
  mp2  <- ggplot() + mapWorld2 +
    #geom_polygon(data = mapWorld4, aes(x=mapWorld4$long, y=mapWorld4$lat, map_id=mapWorld4$group)) +
    #geom_polygon(aes(x=CountryDeets$long, y=CountryDeets$lat))+
    geom_point(aes(x=PASites$Longitude, y=PASites$Latitude, colour=PASites$Treatment), size=1.3, shape=19, alpha=0.9)+ #, colour="dodgerblue4" #, colour="#009E73"
    coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
    scale_colour_manual(values=c("#0B4EA2", "olivedrab3"), guide=FALSE)+
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(size = 1, fill = NA, colour="black"),
      legend.justification = "top",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background=element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"))
  mp2
  
  png(paste0("/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter3/Figures/PreAnalysisPlan/", MapVersion, "_Map.png"), 19, 9, units="in", res=300, bg="transparent") #Save as an image
  grid.newpage()
  v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
  v2<-viewport(width = 0.44, height = 0.39, x = 0, y = 0.02, just=c("left", "bottom")) #plot area for the inset map
  print(mp,vp=v1) 
  print(mp2,vp=v2, mar=c(0,0,0,0))
  dev.off()
}

#### Graveyard ####
###ForTesting
PPOverSites <- over(subset(sitepoints, SiteCode=="15693"), PP, returnList = TRUE)
i <- 45 #PPOverlay 45 is missing!!!!
ok <- sitepointssplit[[i]]
SC <- "15693"
subset(ok, SiteCode=="15693")

PPOverSitesALL <- over(sitepoints, PP)
PPOverSitesALL2 <- PPOverSitesALL[!is.na(PPOverSitesALL$WDPAID),]
#BY THIS METHOD there are 5829 protected sites

#The three with NULL = 153, 159, 179
sitepointssub1 <- sitepoints[(ThemsTheBreaks[153]+1):(ThemsTheBreaks[153+1]),] 
sitepointssub2 <- sitepoints[(ThemsTheBreaks[159]+1):(ThemsTheBreaks[159+1]),] 
sitepointssub3 <- sitepoints[(ThemsTheBreaks[179]+1):(ThemsTheBreaks[179+1]),] 
sitepointssub1 <- rbind(sitepointssub1, sitepointssub2, sitepointssub3)
write.csv(sitepointssub1, "/Users/hannahwauchope/Desktop/THEDEADONES.csv", row.names=FALSE)
### End testing

#For 1 column MDYearList
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
}

#### METT Testing ####
METT <- read.csv("/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter3/METT/Tbl_Spatial.csv")
METTMerge <- unique(merge(PPSitesBA, METT, by.x="WDPAID", by.y="WDPA.ID")[,c("SiteCode", "ISO3", "Country", "WDPAID")])
length(unique(METTMerge2$SiteCode))
length(unique(METTMerge$WDPAID))

#Ok no METT data. 

#### Distance Testing ####
AllDistances2 <- subset(AllDistances, PointDist>0)

#Proportional cutoffs
#Protected site proportions
ProtProp <- dcast(PPSites, ISO3~., length, value.var="SiteCode")
ProtProp <- ProtProp[order(ProtProp$., decreasing= TRUE),]
ProtProp$Prop <- ProtProp$./sum(ProtProp$.)
ProtProp$. <- NULL
names(ProtProp) <- c("ISO3", "Protected")

#Distances proportions
nrow(subset(AllDistances2, PointDist<1000))/nrow(AllDistances2)
nrow(subset(AllDistances2, PointDist<5000))/nrow(AllDistances2)
nrow(subset(AllDistances2, PointDist<10000))/nrow(AllDistances2)
nrow(subset(AllDistances2, PointDist<15000))/nrow(AllDistances2)

CountryCast <- function(DATA, Buffer){
  ok <- dcast(DATA, ISO3~., length, value.var="SiteCode",fill=0)
  ok <- ok[order(ok$.,decreasing = TRUE),]
  ok$Proportion <- round(ok$./sum(ok$.),3)
  #ok <- subset(ok,Proportion>0.00)
  names(ok) <- c("ISO3","NSites","Prop")
  ok$Buffer <- Buffer
  return(ok)
}

SitesCountries <- unique(Spec6_Cov[,c("SiteCode" ,"ISO3")])
AllDistances3 <- merge(AllDistances2, SitesCountries, by="SiteCode")

NoBuff <- CountryCast(AllDistances3,"NoBuffer")
OneKm <- CountryCast(subset(AllDistances3, PointDist>1000),"1km")
FiveKm <- CountryCast(subset(AllDistances3, PointDist>5000),"5km")
TenKm <- CountryCast(subset(AllDistances3, PointDist>10000),"10km")

All <- rbind(NoBuff, OneKm, FiveKm,TenKm)
AllCast <- dcast(All, ISO3~Buffer, value.var="Prop", fill=0)
AllCast <- AllCast[order(AllCast$NoBuffer,decreasing=TRUE),]
AllCast <- AllCast[,c(1,ncol(AllCast),2:(ncol(AllCast)-1))]
rownames(AllCast) <-NULL

AllCast$Sum <- rowSums(AllCast[,c(2:ncol(AllCast))])
AllCast2 <- head(AllCast, 10)
AllCast3 <- subset(AllCast, Sum>0.1)
AllCast4 <- AllCast[AllCast$ISO3 %in% ProtProp$ISO3,]

AllCast5 <- AllCast4

AllCast5 <- merge(AllCast5,ProtProp, by="ISO3")
AllCast5$Sum <- NULL
AllCast5 <- merge(AllCast5,Countries[,c("ISO3", "NAME")], by="ISO3")
AllCast5$ISO3 <- NULL
AllCast5 <- AllCast5[order(AllCast5$Protected, decreasing=TRUE),]
AllCast5 <- head(AllCast5, 10)
AllCastMelt <- melt(AllCast5, id.vars = "NAME")
AllCastMelt$value <- AllCastMelt$value*100
AllCastMelt$variable <- factor(AllCastMelt$variable, levels=c("Protected","NoBuffer","1km","5km","10km"))

ggplot(data=AllCastMelt, aes(x=variable,y=value,group=NAME))+
  geom_point(aes(color=as.character(NAME)))+
  geom_line(aes(color=as.character(NAME)))+ #, linetype="dashed"
  #geom_text(aes(label = AllCastMelt$NAME, colour = AllCastMelt$NAME, x = 5, y = value), hjust = -.1)
  geom_dl(aes(label = NAME,color=NAME), method = list(dl.combine("last.points"), cex = 0.9)) +
  #scale_colour_discrete(guide = 'none')  +  
  xlab("Buffer")+
  ylab("Percentage of Sites")+
  theme(
    panel.background = element_blank(),
    plot.margin = unit(c(0,0,0,0), "lines"),
    panel.grid = element_blank(), 
    text = element_text(size=10, colour = "black"),
    axis.text.x = element_text(size=10, colour = "black"),
    axis.text.y = element_text(size=10, colour = "black"),
    panel.border = element_rect(size = 1, fill = NA),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, size=10, colour = "black"),
    aspect.ratio = 0.7)
