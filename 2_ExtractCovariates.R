############################################################################################################################################################################
### This code was written by Hannah Wauchope to prepare data for the paper "Protected areas have a mixed impact on waterbirds, but management helps"
### This is script 2 of 4 in the workflow
### Last edited 27th August, 2021
### Please direct queries to hannah.wauchope@gmail.com
###
### This script takes the data cleaned in Script 1 and extracts site and species specific covariates: migrant status, climatic variables, landuse and human population variables, governance values, distance to nearest city, slope and surface water
### See Extended Data Table 1 for details of the sources of these covariates, and explanations for any transformations etc
###
### NOTE: Throughout these scripts I use the field "SiteSpec" to mean Population, i.e. a particular species at a particular site. This is generally the unit that analyses are performed on
############################################################################################################################################################################

#### Initialise ####
cluster <- FALSE

if(cluster==FALSE){
  options(repos = c(CRAN = "http://cran.rstudio.com"))
  library(data.table)
  library(pbapply)
  library(pbmcapply)
  library(dplyr)
  library(rgdal)
  library(rgeos)
  library(ncdf4)
  library(chron)
  library(raster)
  library(rlist)
  library(stringr)
  library(tidyverse)
  library(abind)
  library(plyr)
  library(ClusterR)
  DataFP <- "/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Data/"
  
  ncores <- 4
} else {
  library(data.table, lib.loc="/home/hsw34/my-R-libs/")
  library(pbapply, lib.loc="/home/hsw34/my-R-libs/")
  library(pbmcapply, lib.loc="/home/hsw34/my-R-libs/")
  library(dplyr, lib.loc="/home/hsw34/my-R-libs/")
  library(rgdal, lib.loc="/home/hsw34/my-R-libs/")
  library(rgeos, lib.loc="/home/hsw34/my-R-libs/")
  library(ncdf4, lib.loc="/home/hsw34/my-R-libs/")
  library(chron, lib.loc="/home/hsw34/my-R-libs/")
  library(raster, lib.loc="/home/hsw34/my-R-libs/")
  library(rlist, lib.loc="/home/hsw34/my-R-libs/")
  library(stringr, lib.loc="/home/hsw34/my-R-libs/")
  library(tidyverse, lib.loc="/home/hsw34/my-R-libs/")
  DataFP <- "/rds/user/hsw34/hpc-work/data/"
  ResultsFP <- "/rds/user/hsw34/hpc-work/results/c3/"
  ncores <- 32
}
#library(rgdal, lib.loc="/home/hsw34/R/x86_64-redhat-linux-gnu-library/3.3/")

WGSCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
MollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

#### Snap function ####
#I use these two functions to snap sites to the nearest raster cell for the various covariates. Because R is slow, these functions involve exporting data and running part of the analysis in QGIS (all explained in functions)
SnapPoints1 <- function(Variable, TemplateFunction, WaterbirdCounts){
  #Create template raster of same extent/rest as the variable data, but with a unique value for each gridcell
  Template <- TemplateFunction
  Template[!is.na(Template)] <- 1:length(Template[!is.na(Template)])
  crs(Template) <- WGSCRS
  writeRaster(Template, file=paste0(DataFP, "WaterbirdData_2020/Covariates/Snap/", Variable, "/Template.tif"), format="GTiff", overwrite=TRUE)

  #Make a points file of sites
  Sites <- unique(WaterbirdCounts[,c("SiteCode", "Latitude", "Longitude")])
  SitesPoints <- cbind(Sites$Longitude, Sites$Latitude)
  SitesPoints <- SpatialPointsDataFrame(SitesPoints, Sites, proj4string = WGSCRS)
  
  #Find which points don't overlap and write to csv
  RasterExtract <- as.vector(raster::extract(Template, SitesPoints, method="simple"))
  Sites$Polygon <- RasterExtract
  SitesNA <- Sites[is.na(Sites$Polygon),]
  
  SitesNASHP <- cbind(SitesNA$Longitude, SitesNA$Latitude)
  SitesNASHP <- SpatialPointsDataFrame(SitesNASHP, SitesNA, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  SitesNASHPMoll <- spTransform(SitesNASHP, MollCRS)
  
  writeOGR(SitesNASHPMoll, paste0(DataFP, "WaterbirdData_2020/Covariates/Snap/", Variable, "/"), "NASitesMoll", driver="ESRI Shapefile", overwrite=TRUE)
  print("Then in QGIS: Import template tif, convert to mollweide (Raster-->Projections-->Warp Raster) and save as TemplateMoll.tif, restart QGIS, then convert to shapefile (Raster-->Conversion-->Polgonise), name the field UniqueValue . *** Import NASitesMoll. Find nearest polygon for each point and distance to the edge of the polygon (Vector-->NNJoin [make sure NNJoin plugin is installed]). *** Export joined points as csv, NASitesJoined.csv, (Right click on layer-->Export-->Save features as-->change to .csv")
  return(Sites)
}
SnapPoints2 <- function(Variable, WaterbirdCounts, SitesSnap){
  TemplatePoints <- as.data.frame(rasterToPoints(raster(paste0(DataFP, "WaterbirdData_2020/Covariates/Snap/", Variable, "/Template.tif"))))
  names(TemplatePoints) <- c("Longitude", "Latitude", "UniqueValue")
  
  NAPointsSnap <- read.csv(paste0(DataFP, "WaterbirdData_2020/Covariates/Snap/", Variable, "/NASitesJoined.csv"))
  
  NAPointsSnap <- merge(NAPointsSnap[,c("SiteCode", "join_UniqueValu", "distance")], TemplatePoints, by.x="join_UniqueValu", by.y="UniqueValue")
  NAPointsSnap$join_UniqueValu <- NULL
  
  #Now take sites less than 50km from the nearest grid cell
  NAPointsSnapCut <- subset(NAPointsSnap, distance<50000)
  NAPointsSnapCut$distance <- NULL
  SitesSnap$Polygon <- NULL
  Sites <- SitesSnap[!SitesSnap$SiteCode %in% NAPointsSnap$SiteCode]
  Sites <- rbind(Sites, NAPointsSnapCut)
  names(Sites) <- c("SiteCode", paste0("Lat", Variable), paste0("Lon", Variable))
  write.csv(NAPointsSnapCut, paste0(DataFP, "WaterbirdData_2020/Covariates/Snap/", Variable, "/SnapCheck.csv"), row.names=FALSE)
  WaterbirdCounts <- merge(WaterbirdCounts, Sites, by="SiteCode")
  print("Be sure to import the SnapCheck csv into QGIS to check everything looks right. Bear in mind that everything will look off because of how QGIS displays the raster, but if you overlay the original template raster (not TemplateMoll), it will look right")
  return(WaterbirdCounts)
}

#### Initial Clean + Migratory Status ####
BirdCounts <- fread(paste0(DataFP, "WaterbirdData_2020/BirdCounts.csv"))
BirdCounts[,c("V1", "SubnationalCode", "CountryCode", "IWCCountry")] <- NULL

#Have discovered CBC sometimes has multiple coordinates for a site. Let's clean this. 
Sites <- unique(BirdCounts[,c("SiteCode", "Latitude", "Longitude", "Dataset")])
SitesDupl <- dcast(unique(BirdCounts[,c("SiteCode", "Latitude", "Longitude", "Dataset")]), SiteCode + Dataset ~ ., length, value.var="SiteCode")
SitesDupl <- Sites[Sites$SiteCode %in% subset(SitesDupl, .==2)$SiteCode,]

SitesDuplDiff <- rbindlist(lapply(unique(SitesDupl$SiteCode), function(x){
  ok <- subset(SitesDupl, SiteCode==x)
  return(data.frame(SiteCode=x, LatDiff=ok[1, "Latitude"]-ok[2, "Latitude"], LonDiff=ok[1, "Longitude"]-ok[2, "Longitude"]))
}))

#Ok generally the coordinate differences are very small (<0.1), so let's just take the mean of the two. i.e. midpoint
SitesDuplMean <- rbindlist(lapply(unique(SitesDupl$SiteCode), function(x){
  ok <- subset(SitesDupl, SiteCode==x)
  return(data.frame(SiteCode=x, LatitudeUpdate=mean(ok$Latitude), LongitudeUpdate=mean(ok$Longitude)))
}))

SitesUpdate <- merge(Sites, SitesDuplMean, all=T)
SitesUpdate$LatitudeUpdate2 <- ifelse(is.na(SitesUpdate$LatitudeUpdate), SitesUpdate$Latitude, SitesUpdate$LatitudeUpdate)
SitesUpdate$LongitudeUpdate2 <- ifelse(is.na(SitesUpdate$LongitudeUpdate), SitesUpdate$Longitude, SitesUpdate$LongitudeUpdate)

Sites <- unique(SitesUpdate[,c("SiteCode", "LatitudeUpdate2", "LongitudeUpdate2")])
if(nrow(Sites)!=length(unique(Sites$SiteCode))){stop("There are still duplicate coordinates!")}
names(Sites) <- c("SiteCode", "Latitude", "Longitude")
Sites <- Sites[complete.cases(Sites),]

BirdCounts[,c("Latitude", "Longitude")] <- NULL
BirdCounts2 <- merge(BirdCounts, Sites, by="SiteCode")

WaterbirdFamilies <- read.csv(paste0(DataFP, "WaterbirdData_2020/WaterbirdFamilies.csv"))

WaterbirdCounts <- BirdCounts2[BirdCounts2$Family %in% WaterbirdFamilies$WaterbirdFamilies,]
rm(BirdCounts, BirdCounts2)

SpeciesPoints <- unique(WaterbirdCounts[,c("Species","SISRecID", "SiteCode", "Longitude", "Latitude")])
UniqueSpecies <- unique(SpeciesPoints$SISRecID)

SpeciesPoints$Check <- paste0(SpeciesPoints$Species, "_", SpeciesPoints$SISRecID, "_", SpeciesPoints$SiteCode)
SpeciesPointsDup <- SpeciesPoints[duplicated(SpeciesPoints$Check),]
if(nrow(SpeciesPointsDup)!=0){stop("there are duplicates!")}

#Read in distribution polygons, and extract the relevant ones. Nb, these are split into 3 shapefiles because the geodatabase couldn't be read in R
#Polygon extraction lines hashed out cos they take forever so only run once

#BirdPolygons1 <- readOGR(paste0(DataFP, "BirdlifePolygons/BOTW_Shape/BOTW_Shape_1_3000/BOTW_Shape_1_3000.shp"))
#BirdPolygons2 <- readOGR(paste0(DataFP, "BirdlifePolygons/BOTW_Shape/BOTW_Shape_3001_10000/BOTW_Shape_3001_10000.shp"))
#BirdPolygons3 <- readOGR(paste0(DataFP, "BirdlifePolygons/BOTW_Shape/BOTW_Shape_10001_17463/BOTW_Shape_10001_17463.shp"))

#Extract only the relevant polygons from the shapefiles
#SpeciesPolygons1 <- BirdPolygons1[BirdPolygons1$SISRecID %in% UniqueSpecies,]
#SpeciesPolygons2 <- BirdPolygons2[BirdPolygons2$SISRecID %in% UniqueSpecies,]
#SpeciesPolygons3 <- BirdPolygons3[BirdPolygons3$SISRecID %in% UniqueSpecies,]
#rm(BirdPolygons1, BirdPolygons2, BirdPolygons3)

#SpeciesPolygons <- list(SpeciesPolygons1, SpeciesPolygons2, SpeciesPolygons3)
#save(SpeciesPolygons, file=paste0(DataFP, "BirdlifePolygons/WaterbirdSpeciesPolygons.RData")) #Save this out immediately cos R threatens breaking (need a new laptop...)

#load(file=paste0(DataFP, "BirdlifePolygons/WaterbirdSpeciesPolygons.RData")) #Read back in and bind the three dataframes together, save out again
#SpeciesPolygons <- do.call(rbind, SpeciesPolygons)
#save(SpeciesPolygons, file=paste0(DataFP, "BirdlifePolygons/WaterbirdSpeciesPolygons.RData"))

load(file=paste0(DataFP, "BirdlifePolygons/WaterbirdSpeciesPolygons.RData"))

#This loop gets the nearest seasonal polygon for every point for every species.
SitesWithSeasonalFunc <- function(SpeciesPolygons){
  rbindlist(pbmclapply(unique(SpeciesPolygons$SISRecID), function(i){
    SpecPolygon <- subset(SpeciesPolygons, SISRecID==i)
    SpecPoints <- subset(SpeciesPoints, SISRecID==i)
    SitesPoints <- cbind(SpecPoints$Longitude, SpecPoints$Latitude)
    SitesPoints <- SpatialPointsDataFrame(SitesPoints, SpecPoints, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    closestpolygon <- suppressWarnings(gDistance(SpecPolygon, SitesPoints, byid=TRUE)) #No you DON'T want Hausdorff distance, you checked!
    colnames(closestpolygon) <- SpecPolygon$SEASONAL
    if(nrow(SpecPolygon)>1 & nrow(SpecPoints)>1){
      SpecPoints$Season <- colnames(closestpolygon[,apply(closestpolygon,1,function(x)return(array(which.min(x))))])
    } else if(nrow(SpecPoints)>1 & nrow(SpecPolygon)==1){
      SpecPoints$Season <- unique(SpecPolygon$SEASONAL)
    } else if(nrow(SpecPolygon)>1 & nrow(SpecPoints)==1){
      SpecPoints$Season <- colnames(closestpolygon)[apply(closestpolygon,1,which.min)]
    } else {
      SpecPoints$Season <- unique(SpecPolygon$SEASONAL)
    }
    return(SpecPoints)
  }, mc.cores=4))
}
SitesWithSeasonal <- SitesWithSeasonalFunc(SpeciesPolygons)

SitesWithSeasonal$MigStatus <- as.factor(SitesWithSeasonal$Season)
SitesWithSeasonal$MigStatus <- recode_factor(SitesWithSeasonal$MigStatus, "1"="Resident", "2"="Breeding", "3"="Non-breeding", "4"="Passage", "5"="Uncertain")
SitesWithSeasonal <- SitesWithSeasonal[,c("SISRecID", "SiteCode", "MigStatus")]

WaterbirdCounts <- merge(WaterbirdCounts, SitesWithSeasonal, by = c("SISRecID", "SiteCode"))

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatus.RData"))

#### Climate ####
#Data from https://catalogue.ceda.ac.uk/uuid/10d3e3640f004c578403419aac167d82

rastmonyear <- function(month, year, climtype){
  climate_output <- nc_open(list.files(path=paste0(DataFP, "/Climate_CRU/ClimateData"),pattern=paste0("*", climtype, ".dat.nc"), full.names=TRUE))
  
  #Get data
  lon <- ncvar_get(climate_output,"lon")
  lat <- ncvar_get(climate_output,"lat",verbose=F)
  time <- ncvar_get(climate_output,"time")
  tunits <- ncatt_get(climate_output,"time","units")
  fillvalue <- ncatt_get(climate_output, climtype,"_FillValue")
  clim_array <- ncvar_get(climate_output,climtype)
  clim_array[clim_array==fillvalue$value] <- NA #Change NAs to appropriate thing
  dlname <- ncatt_get(climate_output,climtype,"long_name")
  dunits <- ncatt_get(climate_output, climtype,"units")
  
  #Get metadata
  title <- ncatt_get(climate_output,0,"title")
  institution <- ncatt_get(climate_output,0,"institution")
  datasource <- ncatt_get(climate_output,0,"source")
  references <- ncatt_get(climate_output,0,"references")
  history <- ncatt_get(climate_output,0,"history")
  Conventions <- ncatt_get(climate_output,0,"Conventions")
  nc_close(climate_output)
  
  #Get the right times
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  dates <- chron(time,origin=c(tmonth, tday, tyear))
  dates <- rbindlist(lapply(1:length(dates), function(x) as.data.frame(t(strsplit(as.character(dates[x]), "/")[[1]][c(1,3)]))))
  names(dates) <- c("Month", "Year")
  dates$RasNum <- rownames(dates)
  
  #Set dates to start from 1920
  dates <- dates[(19*12+1):nrow(dates),]
  
  #Extract raster of designated month/year
  yearcode <- substr(as.numeric(as.character(year)), nchar(as.numeric(as.character(year)))-2+1, nchar(as.numeric(as.character(year))))
  m <- as.numeric(as.character(dates[dates$Month==month & dates$Year==yearcode,]$RasNum))
  clim_slice <- clim_array[,,m]
  #Ok so that gets alll our data out
  # matrix (nlon*nlon rows by 2 cols) of lons and lats
  lonlat <- as.matrix(expand.grid(lon,lat))
  # vector of values
  clim_vec <- as.vector(clim_slice)
  clim_df <- data.frame(cbind(lonlat,clim_vec))
  names(clim_df) <- c("lon","lat",paste(climtype,as.character(m), sep="_"))
  clim_ras <- rasterFromXYZ(clim_df)  #Convert first two columns as lon-lat and third as value                
  return(clim_ras)
} #All function variables need to be in quotations. This extracts a global raster of the correct month/year from the netcdf files

load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatus.RData"))
names(WaterbirdCounts)[names(WaterbirdCounts)=="Season.x"] <- "Season"
names(WaterbirdCounts)[names(WaterbirdCounts)=="Season.y"] <- "MigStatus"

#Snap coordinates and ditch any that don't fit
SitesSnap <- SnapPoints1("Clim", rastmonyear("01", paste(1999), "tmp"), WaterbirdCounts)
WaterbirdCounts <- SnapPoints2("Clim", WaterbirdCounts, SitesSnap) 

#Create a YearClim variable, so that yearly rainfall etc is taken from the past not the future.
WaterbirdCounts <-subset(WaterbirdCounts, Year>1939) #Start from 1940 when we have at least 100 sites
WaterbirdCounts$YearClim <- WaterbirdCounts$Year
WaterbirdCounts[WaterbirdCounts$Season=="DectoFeb",]$YearClim <- (WaterbirdCounts[WaterbirdCounts$Season=="DectoFeb",]$Year)-1

SitesYearSeason <- unique(WaterbirdCounts[,c("YearClim", "SiteCode", "LonClim", "LatClim", "Season")])
SitesYearSeason <- subset(SitesYearSeason, YearClim<2018)
SitesYearSeason <- SitesYearSeason[order(SitesYearSeason$YearClim, decreasing=TRUE),]
SitesYearSeason$YearClim <- as.factor(SitesYearSeason$YearClim)

TempExtractFun <- function(TempList, SitePoints){
  list.cbind(lapply(c("Mean", "Min", "Max"), function(summary){
    Seasons <- as.data.frame(lapply(names(TempList), function(season){
      TempRast <- overlay(TempList[[season]], fun=paste0(tolower(summary)))
      ExtractTemp <- raster::extract(TempRast, SitePoints, method="bilinear")
      return(ExtractTemp)
    }))
    names(Seasons) <- paste0(summary, names(TempList), "Temp")
    return(Seasons)
  }))
}
PreExtractFun <- function(PreList, SitePoints){
  PrecipitationExtract <- as.data.frame(lapply(names(PreList), function(season){
    PreRast <- overlay(PreList[[season]], fun=sum)
    ExtractPre <- raster::extract(PreRast, SitePoints, method="bilinear")
    return(ExtractPre)
  }))
  names(PrecipitationExtract) <- paste0("Total", names(PreList), "Precip")
  return(PrecipitationExtract)
}

ClimExtract <- pbmclapply(as.numeric(as.character(unique(SitesYearSeason$YearClim))), function(x){
  if(file.exists(paste0(DataFP, "/WaterbirdData_2020/Covariates/ClimExtract/Extract_", x, ".csv"))){return(NULL)}
  SitesDectoFeb <- subset(SitesYearSeason, Season=="DectoFeb" & YearClim==x)
  if(nrow(unique(SitesDectoFeb))!=nrow(SitesDectoFeb) | nrow(SitesDectoFeb) != nrow(SitesDectoFeb[complete.cases(SitesDectoFeb),])){stop("There's an issue with sites subset")}
  
  SitesPointsDectoFeb <- cbind(SitesDectoFeb$LonClim, SitesDectoFeb$LatClim)
  SitesPointsDectoFeb <- SpatialPointsDataFrame(SitesPointsDectoFeb, SitesDectoFeb, proj4string = WGSCRS)
  
  RastersTemp <- stack(lapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(month) rastmonyear(month, paste(x), "tmp")))
  RastersTempSeason <- stack(append(RastersTemp[[12]], lapply(c("01","02"), function(month) rastmonyear(month, paste(x+1), "tmp"))))
  TempList <- list(RastersTemp, RastersTempSeason)
  names(TempList) <- c("Annual", "Seasonal")
  
  RastersPre <- stack(lapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(month) rastmonyear(month, paste(x), "pre")))
  RastersPreSeason <- stack(append(RastersPre[[12]], lapply(c("01","02"), function(month) rastmonyear(month, paste(x+1), "pre"))))
  PreList <- list(RastersPre, RastersPreSeason)
  names(PreList) <- c("Annual", "Seasonal")

  TemperatureExtract <- TempExtractFun(TempList, SitesPointsDectoFeb)
  PrecipitationExtract <- PreExtractFun(PreList, SitesPointsDectoFeb)
  
  SitesExtracted <- as.data.frame(do.call(cbind, list(SitesDectoFeb, TemperatureExtract, PrecipitationExtract)))

  SitesJuntoAug <- subset(SitesYearSeason, Season=="JuntoAug" & YearClim==(x))
  
  if(nrow(SitesJuntoAug)!=0){
    if(nrow(unique(SitesJuntoAug))!=nrow(SitesJuntoAug) | nrow(SitesJuntoAug) != nrow(SitesJuntoAug[complete.cases(SitesJuntoAug),])){stop("There's an issue with jun to aug sites subset")}
    SitesPointsJuntoAug <- cbind(SitesJuntoAug$LonClim, SitesJuntoAug$LatClim)
    SitesPointsJuntoAug <- SpatialPointsDataFrame(SitesPointsJuntoAug, SitesJuntoAug, proj4string = WGSCRS)
    
    RastersTemp1 <- lapply(c("07","08","09","10","11","12"), function(month) rastmonyear(month, paste(x-1), "tmp"))
    RastersTemp2 <- lapply(c("01","02","03","04","05","06"), function(month) rastmonyear(month, paste(x), "tmp"))
    RastersTemp <- stack(append(RastersTemp1, RastersTemp2))
    
    RastersTempSeason <- stack(append(RastersTemp[[12]], lapply(c("07","08"), function(month) rastmonyear(month, paste(x), "tmp"))))
    TempList <- list(RastersTemp, RastersTempSeason)
    names(TempList) <- c("Annual", "Seasonal")
    
    RastersPre1 <- lapply(c("07","08","09","10","11","12"), function(month) rastmonyear(month, paste(x-1), "pre"))
    RastersPre2 <- lapply(c("01","02","03","04","05","06"), function(month) rastmonyear(month, paste(x), "pre"))
    RastersPre <- stack(append(RastersPre1, RastersPre2))
    
    RastersPreSeason <- stack(append(RastersPre[[12]], lapply(c("07","08"), function(month) rastmonyear(month, paste(x), "pre"))))
    PreList <- list(RastersPre, RastersPreSeason)
    names(PreList) <- c("Annual", "Seasonal")
    
    TemperatureExtract <- TempExtractFun(TempList, SitesPointsJuntoAug)
    PrecipitationExtract <- PreExtractFun(PreList, SitesPointsJuntoAug)
    
    SitesExtractedJuntoAug <- as.data.frame(do.call(cbind, list(SitesJuntoAug, TemperatureExtract, PrecipitationExtract)))
    
    SitesExtracted <- rbind(SitesExtracted, SitesExtractedJuntoAug)
  }
  write.csv(SitesExtracted, paste0(DataFP, "/WaterbirdData_2020/Covariates/ClimExtract/Extract_", x, ".csv"), row.names=FALSE)
  return(SitesExtracted)
}, mc.cores=2) #and this overlays points and extracts data

ClimExtract <- rbindlist(lapply(list.files(path=paste0(DataFP, "/WaterbirdData_2020/Covariates/ClimExtract/"), full.names=TRUE), function(x) fread(x)))
ClimExtract[,c("LonClim", "LatClim")] <- NULL

WaterbirdCounts <- merge(WaterbirdCounts, ClimExtract, by=c("YearClim", "SiteCode", "Season"))

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimate.RData"))

#### Fertiliser ####
load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimate.RData"))

SitesSnap <- SnapPoints1("Fert", raster(paste0(DataFP, "/Lu-Tian_2016_NP/Nfer_ASCII/nfery2000.asc")), WaterbirdCounts)
WaterbirdCounts <- SnapPoints2("Fert", WaterbirdCounts, SitesSnap) 

SitesYearSeason <- unique(WaterbirdCounts[,c("YearClim", "Season", "SiteCode", "LonFert", "LatFert")])
SitesYearSeason <- SitesYearSeason[order(SitesYearSeason$YearClim),]
SitesYearSeasonFert <- subset(SitesYearSeason, YearClim<2014)

NitrogenPhosphorous <- rbindlist(pbmclapply(unique(SitesYearSeasonFert$YearClim), function(x){
  print(x)
  Nitrogen <- raster(paste0(DataFP, "/Lu-Tian_2016_NP/Nfer_ASCII/nfery", (x), ".asc"))
  Phosphorous <- raster(paste0(DataFP, "Lu-Tian_2016_NP/Pfer_ASCII/pfery", (x), ".asc"))
  
  Sites <- subset(SitesYearSeasonFert, YearClim==x)
  if(nrow(Sites)!=unique(nrow(Sites))){stop("There's duplicates")}
  
  SitesDectoFeb <- subset(Sites, Season=="DectoFeb")

  SitesDectoFebPoints <- cbind(SitesDectoFeb$LonFert, SitesDectoFeb$LatFert)
  SitesDectoFebPoints <- SpatialPointsDataFrame(SitesDectoFebPoints, SitesDectoFeb, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  NitrogenExtract <- raster::extract(Nitrogen, SitesDectoFebPoints, method="bilinear")
  PhosphorousExtract <- raster::extract(Phosphorous, SitesDectoFebPoints, method="bilinear")
  
  SitesDectoFeb$Nitr <- NitrogenExtract
  SitesDectoFeb$Phos <- PhosphorousExtract
  
  SitesJuntoAug <- subset(Sites, Season=="JuntoAug")
  
  if(nrow(SitesJuntoAug)!=0){
    SitesJuntoAugPoints <- cbind(SitesJuntoAug$LonFert, SitesJuntoAug$LatFert)
    SitesJuntoAugPoints <- SpatialPointsDataFrame(SitesJuntoAugPoints, SitesJuntoAug, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
    NitrogenExtract <- raster::extract(Nitrogen, SitesJuntoAugPoints, method="bilinear")
    PhosphorousExtract <- raster::extract(Phosphorous, SitesJuntoAugPoints, method="bilinear")

    NitrogenYearBefore <- raster(paste0(DataFP, "/Lu-Tian_2016_NP/Nfer_ASCII/nfery", (x-1), ".asc"))
    PhosphorousYearBefore <- raster(paste0(DataFP, "Lu-Tian_2016_NP/Pfer_ASCII/pfery", (x-1), ".asc"))
    NitrogenYearBeforeExtract <- raster::extract(NitrogenYearBefore, SitesJuntoAugPoints, method="bilinear")
    PhosphorousYearBeforeExtract <- raster::extract(PhosphorousYearBefore, SitesJuntoAugPoints, method="bilinear")
    
    SitesJuntoAug$NitrYear <- NitrogenExtract
    SitesJuntoAug$PhosYear <- PhosphorousExtract
    
    SitesJuntoAug$NitrYearBefore <- NitrogenYearBeforeExtract
    SitesJuntoAug$PhosYearBefore <- PhosphorousYearBeforeExtract
    
    SitesJuntoAug$Nitr <- rowMeans(SitesJuntoAug[,c("NitrYear", "NitrYearBefore")])
    SitesJuntoAug$Phos <- rowMeans(SitesJuntoAug[,c("PhosYear", "PhosYearBefore")])
    SitesDectoFeb <- rbind(SitesDectoFeb, SitesJuntoAug[,c("YearClim", "Season", "SiteCode", "LonFert", "LatFert", "Nitr", "Phos")])
  }
  return(SitesDectoFeb)
}, mc.cores=4))

NitrogenPhosphorousNAs <- subset(SitesYearSeason, YearClim>=2014)
NitrogenPhosphorousNAs[,c("Nitr", "Phos")] <- NA
NitrogenPhosphorous <- rbind(NitrogenPhosphorous, NitrogenPhosphorousNAs)
NitrogenPhosphorous[,c("LonFert", "LatFert")] <- NULL

WaterbirdCounts <- merge(WaterbirdCounts, NitrogenPhosphorous, by=c("YearClim", "SiteCode", "Season"))

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliser.RData"))

#### HYDE Interpolation ####
#HYDE Data from: https://easy.dans.knaw.nl/ui/datasets/id/easy-dataset:74467

InterpolateValues <- function(HYDEType){
  if(HYDEType=="LandUse"){
    Variables <- sapply(list.files(paste0(DataFP, "HYDE/LandUse/1940AD_lu/")), function(x) str_split_fixed(x, "[1940AD.]",2)[,1])
    HYDETypeShort <- "lu"
  } else {
    Variables <- sapply(list.files(paste0(DataFP, "HYDE/Population/1940AD_pop/")), function(x) str_split_fixed(x, "[1940AD.]",2)[,1])
    HYDETypeShort <- "pop"
  }
  
  VarByYear <- expand.grid(Variables, c(seq(1940, 2000, 10)))
  pbmclapply(c(1:nrow(VarByYear)), function(row){
    var <- VarByYear[row,]$Var1
    Year <- VarByYear[row,]$Var2

    if(file.exists(paste0(DataFP, "HYDE/", HYDEType, "/InterpolatedValues/", var, Year-1, "AD.asc"))){
      return(NULL)
    }
    
    print(paste0("Begin ", var, " ", Year))
    
    #Import rasters from the decades you want to interpolate between
    range1 <- raster(paste0(DataFP, "HYDE/", HYDEType, "/", Year, "AD_", HYDETypeShort,"/", var, Year, "AD.asc"))
    range2 <- raster(paste0(DataFP, "HYDE/", HYDEType, "/", (Year-10), "AD_", HYDETypeShort,"/", var, (Year-10), "AD.asc"))
    
    range1points <- as.data.frame(rasterToPoints(range1))
    range2points <- as.data.frame(rasterToPoints(range2))
    
    Ranges <- list(range1points, range2points) %>% reduce(left_join, by=c("x", "y"))
    names(Ranges) <- c("Lon", "Lat", Year, (Year-10))
    
    #Interpolate each value
    Interpolate <- rbindlist(lapply(1:(nrow(Ranges)), function(x){
      RangePixel <- as.data.frame(t(Ranges[x,c(3:ncol(Ranges))]))
      RangePixel$Year <- as.numeric(rownames(RangePixel))
      names(RangePixel) <- c("Var", "Year")
      model <- lm(Var~Year, RangePixel)
      
      InterpolatedValues <- data.frame("Year"=c((Year-10):Year))
      InterpolatedValues$Values <- predict(model, InterpolatedValues)
      
      InterpolatedValues <- as.data.frame(t(InterpolatedValues))
      names(InterpolatedValues) <- InterpolatedValues[1,]
      InterpolatedValues <- InterpolatedValues[2,]
      InterpolatedValues$RasterID <- x
      return(InterpolatedValues)
    }))
    
    Coordinates <- Ranges[,c(1,2)]
    RasterData <- cbind(Coordinates, Interpolate)
    print(paste0("Write rasters ", var, " ", Year))
    
    #Write out the interpolated rasters
    for(x in c(4:(ncol(RasterData)-2))){
      RasterYear <- RasterData[,c(1,2,x)]
      RasterYear <- rasterFromXYZ(RasterYear)
      writeRaster(RasterYear, filename=paste0(DataFP, "HYDE/", HYDEType, "/InterpolatedValues/", var, names(RasterData)[x], "AD.asc"), format="ascii", overwrite=TRUE)
    }
  }, mc.cores=ncores)
}
InterpolateValues("LandUse")
InterpolateValues("Population")

#### HYDE ####

##Land Use
#croplandyear.asc     (total cropland area, in km2 per grid cell), after 1960 identical to FAO's category 'Arable land and permanent crops'.
#grazingyear.asc      (total land used for grazing, in km2 per grid cell), after 1960 identical to FAO's category 'Permanent Pasture'.
#pastureyear.asc      (total pasture area, in km2 per grid cell), defined as Grazing land with an aridity index > 0.5, assumed to be more intensively managed.
#rangelandyear.asc    (total pasture area, in km2 per grid cell), defined as Grazing land with an aridity index > 0.5, assumed to be less or not managed.
#ir_riceyear.asc    (total irrigated rice area, in km2 per grid cell).
#rf_riceyear.asc    (total rainfed area, in km2 per grid cell).
#ir_noriceyear.asc      (total irrigated other crops area (no rice), in km2 per grid cell).
#rf_noriceyear.asc      (total rainfed other crops area (no rice), in km2 per grid cell).
#tot_irriyear.asc     (total actual irrigated area, in km2 per grid cell).
#tot_rainfedyear.asc  (total rainfed area, in km2 per grid cell).
#tot_riceyear.asc     (total rice area, in km2 per grid cell).

##Population
#popc_year.asc (population counts, in inhabitants/gridcell)
#popd_year.asc (population density, in inhabitants/km2 per gridcell)
#rurc_year.asc (rural population counts, in inh/gridcell)
#urb_year.asc  (urban population counts, in inh/gridcell)
#uopp_year.asc (total built-up area, such as towns, cities, etc, in km2 per grid cell)
#We don't care about rural and urban population counts, as these are just a division of "popc" into rural and urban areas

#"grazing", "ir_norice", "ir_rice", "pasture", "rangeland", "rf_norice", "rf_rice", "popd", "uopp", "popc"

##Anthrome
#"11 Urban"
#"12 Dense settlements"

#"21 Village, Rice"
#"22 Village, Irrigated"
#"23 Village, Rainfed"
#"24 Village, Pastoral"

#"31 Croplands, residential irrigated"
#"32 Croplands, residential rainfed"
#"33 Croplands, populated"
#"34 Croplands, pastoral"

#"41 Rangeland, residential"
#"42 Rangeland, populated"
#"43 Rangeland, remote"

#"51 Semi-natural woodlands, residential"
#"52 Semi-natural woodlands, populated"
#"53 Semi-natural woodlands, remote"
#"54 Semi-natural treeless and barren lands"

#"61 Wild, remote - woodlands"
#"62 Wild, remote - treeless & barren"
#"63 Wild, remote - ice"

#"70 No definition"

#Unlike populations and landuse it doesn't make sense to linearly interpolate anthromes as they're categorical. Therefore, we'll take the value from the nearest decade. 
#So all years (below 2000, cos we have yearly anthrome data after that) are rounded to the nearest ten and the anthrome data taken from that 

load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliser.RData"))

SitesSnap <- SnapPoints1("HYDE", raster(paste0(DataFP, "/HYDE/Anthrome/2017AD_anthromes/anthromes2017AD.asc")), WaterbirdCounts)
WaterbirdCounts <- SnapPoints2("HYDE", WaterbirdCounts, SitesSnap) 

SitesYearSeason <- unique(WaterbirdCounts[,c("YearClim", "SiteCode", "Season", "LonHYDE", "LatHYDE")])
SitesYearSeason <- SitesYearSeason[order(SitesYearSeason$YearClim),]

HYDE <- rbindlist(pbmclapply(unique(SitesYearSeason$YearClim), function(x){
  print(x)
  SitesDectoFeb <- subset(SitesYearSeason, YearClim==x & Season=="DectoFeb")
  SitesPoints <- cbind(SitesDectoFeb$LonHYDE, SitesDectoFeb$LatHYDE)
  SitesPoints <- SpatialPointsDataFrame(SitesPoints, SitesDectoFeb, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  LURasters <- lapply(list.files(path = paste0(DataFP, "HYDE/LandUse/"),
                                           pattern = paste0(x, "AD.asc$"), recursive = TRUE, full.names = TRUE), function(x){
                                             ras <- raster(x)
                                             crs(ras) <- WGSCRS
                                             return(ras)})
  PopRasters <- lapply(list.files(path = paste0(DataFP, "HYDE/Population/"),
                                            pattern = paste0(x, "AD.asc$"), recursive = TRUE, full.names = TRUE), function(x){
                                              ras <- raster(x)
                                              crs(ras) <- WGSCRS
                                              return(ras)})
  HYDERasters <- c(LURasters, PopRasters)

  RasterExtract <- lapply(HYDERasters, function(y) raster::extract(y, SitesPoints, method="bilinear"))
  RasterExtractDF <- as.data.frame(do.call(cbind, RasterExtract))
  RasterNames <- sapply(HYDERasters, names)
  RasterNames <- str_split_fixed(RasterNames, paste(x), 2)[,1]
  names(RasterExtractDF) <- RasterNames
  SitesDectoFeb <- cbind(SitesDectoFeb, RasterExtractDF)
  
  SitesJuntoAug <- subset(SitesYearSeason, YearClim==x & Season=="JuntoAug")
  if(nrow(SitesJuntoAug)!=0){
    SitesPoints <- cbind(SitesJuntoAug$LonHYDE, SitesJuntoAug$LatHYDE)
    SitesPoints <- SpatialPointsDataFrame(SitesPoints, SitesJuntoAug, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    LURastersYearBefore <- lapply(list.files(path = paste0(DataFP, "HYDE/LandUse/"),
                                   pattern = paste0(x-1, "AD.asc$"), recursive = TRUE, full.names = TRUE), function(x){
                                     ras <- raster(x)
                                     crs(ras) <- WGSCRS
                                     return(ras)})
    PopRastersYearBefore <- lapply(list.files(path = paste0(DataFP, "HYDE/Population/"),
                                    pattern = paste0(x-1, "AD.asc$"), recursive = TRUE, full.names = TRUE), function(x){
                                      ras <- raster(x)
                                      crs(ras) <- WGSCRS
                                      return(ras)})
    HYDERastersYearBefore <- c(LURastersYearBefore, PopRastersYearBefore)
    
    RasterExtract <- lapply(HYDERasters, function(y) raster::extract(y, SitesPoints, method="bilinear"))
    RasterExtractDF <- as.data.frame(do.call(cbind, RasterExtract))
    RasterNames <- sapply(HYDERasters, names)
    RasterNames <- str_split_fixed(RasterNames, paste(x), 2)[,1]
    names(RasterExtractDF) <- RasterNames
    
    RasterExtractYearBefore <- lapply(HYDERastersYearBefore, function(y) raster::extract(y, SitesPoints, method="bilinear"))
    RasterExtractDFYearBefore <- as.data.frame(do.call(cbind, RasterExtractYearBefore))
    RasterNames <- sapply(HYDERastersYearBefore, names)
    RasterNames <- str_split_fixed(RasterNames, paste(x-1), 2)[,1]
    names(RasterExtractDFYearBefore) <- RasterNames
    
    RasterExtractDFMean <- abind(RasterExtractDF, RasterExtractDFYearBefore, along = 3)
    RasterExtractDFMean <- as.data.frame(rowMeans(RasterExtractDFMean, dims = 2))
    
    SitesJuntoAug <- cbind(SitesJuntoAug, RasterExtractDFMean)
    
    SitesDectoFeb <- rbind(SitesDectoFeb, SitesJuntoAug)
  }
  return(SitesDectoFeb)
}, mc.cores=ncores)) #and this overlays points and extracts data

SitesYearSeason$YearRound <- ifelse(SitesYearSeason$YearClim<2000, round_any(SitesYearSeason$YearClim, 10), SitesYearSeason$YearClim)
HYDEAnthromes <- rbindlist(pblapply(unique(SitesYearSeason$YearRound), function(x){
  SitesSub <- subset(SitesYearSeason, YearRound==x)
  SitesPoints <- cbind(SitesSub$LonHYDE, SitesSub$LatHYDE)
  SitesPoints <- SpatialPointsDataFrame(SitesPoints, SitesSub, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  AnthromeRaster <- raster(list.files(path=paste0(DataFP, "HYDE/Anthrome/"), pattern=paste0("*", x, "AD.asc"), full.names=TRUE, recursive = TRUE))
  
  RasterExtractDF <- as.data.frame(raster::extract(AnthromeRaster, SitesPoints, method="simple"))
  names(RasterExtractDF) <- "Anthrome"
  SitesSub <- cbind(SitesSub, RasterExtractDF)
  return(SitesSub)
})) #and this overlays points and extracts data

HYDE <- merge(HYDE, HYDEAnthromes)
if(nrow(HYDE)!=nrow(SitesYearSeason)){stop("Lost rows!")}

HYDE[,c("LonHYDE", "LatHYDE", "YearRound")] <- NULL
WaterbirdCounts <- merge(WaterbirdCounts, HYDE, by=c("YearClim", "SiteCode", "Season"))

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDE.RData"))

#### Countries and Governance ####
load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDE.RData"))
Countries <- readOGR("/Users/hannahwauchope/Documents/ArcGIS/CountryContinentWorldSHPS/World/TM_WORLD_BORDERS-0.3-SSudan.shp")

Sites <- unique(WaterbirdCounts[,c("SiteCode", "Latitude", "Longitude")])
SitePoints <- cbind(Sites$Longitude, Sites$Latitude)
SitePoints <- SpatialPointsDataFrame(SitePoints, Sites, proj4string = WGSCRS)

CountriesOverSites <- over(SitePoints, Countries)
WaterbirdCountries <- cbind(as.data.frame(Sites), CountriesOverSites)

#Take the points that fall just outside of country polygons
NASites <- WaterbirdCountries[is.na(WaterbirdCountries$ISO3),]
NAPoints <- cbind(NASites$Longitude, NASites$Latitude)
NAPoints <- SpatialPointsDataFrame(NAPoints, NASites, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

##  Project data into an equal area projection
CountriesMoll <- spTransform(Countries, MollCRS)
NAPointsMoll <- spTransform(NAPoints, MollCRS)

## For each point, find name (and details) of nearest country
NASnap <- rbindlist(pbmclapply(1:nrow(NAPointsMoll), function(i) as.data.frame(CountriesMoll[which.min(gDistance(NAPointsMoll[i,], CountriesMoll, byid=TRUE)),]), mc.cores=ncores))
NAPoints <- as.data.frame(NAPointsMoll)
NAPoints <- cbind(NAPoints[,c(1:3)], NASnap)

WaterbirdCountries <- rbind(WaterbirdCountries[!is.na(WaterbirdCountries$ISO3),], NAPoints)[,c("SiteCode", "ISO3", "NAME", "REGION", "SUBREGION")]
names(WaterbirdCountries) <- c("SiteCode", "ISO3", "Country", "GeoRegion", "GeoSubRegion")

#Now get governance data
Governance <- read.csv(file=paste0(DataFP, "Governance/WGI_csv/WGIData.csv"))
Governance <- Governance[grep("*.EST", Governance$Indicator.Code),]
Governance[,c("Country.Name", "Indicator.Name", "X")] <- NULL
Governance <- melt(as.data.table(Governance), id.vars = c("Country.Code", "Indicator.Code"))
Governance <- Governance[complete.cases(Governance),]
Governance <- dcast(Governance, Country.Code~variable, mean, value.var="value")
Governance <- melt(Governance, id.vars = "Country.Code")
Governance <- subset(Governance, value!="NaN")
Governance$variable <- str_split_fixed(Governance$variable, "X", 2)[,2]
names(Governance) <- c("ISO3", "Year", "Gov")

#Need to impute values for 1997, 1999 and 2001 - will take the mean

ImputedYears <- rbindlist(lapply(c(1997, 1999, 2001), function(x){
  Impute <- subset(Governance, Year==x-1 | Year==x+1)
  Impute <- dcast(Impute, ISO3~Year, value.var="Gov")
  Impute$NewYear <- sapply(1:nrow(Impute), function(y) mean(c(as.numeric(Impute[y, 2]), as.numeric(Impute[y,3]))))
  Impute <- Impute[,c(1,4)]
  Impute$Year <- x
  names(Impute) <- c("ISO3", "Gov", "Year")
  Impute <- Impute[,c(1,3,2)]
  return(Impute)
}))
Governance <- rbind(Governance, ImputedYears)

GovernanceJuntoAug <- rbindlist(lapply(unique(Governance$Year)[2:length(unique(Governance$Year))], function(x){
  x <- as.numeric(x)
  JuntAug <- subset(Governance, Year==x-1 | Year==x)
  JuntAug <- dcast(JuntAug, ISO3~Year, value.var="Gov")
  JuntAug$NewYear <- sapply(1:nrow(JuntAug), function(y) mean(c(as.numeric(JuntAug[y, 2]), as.numeric(JuntAug[y,3]))))
  JuntAug <- JuntAug[,c(1,4)]
  JuntAug$Year <- x
  names(JuntAug) <- c("ISO3", "Gov", "Year")
  JuntAug <- JuntAug[,c(1,3,2)]
  JuntAug$Season <- "JuntoAug"
  return(JuntAug)
}))
Governance$Season <- "DectoFeb"
Governance <- rbind(Governance, GovernanceJuntoAug)

#Combine governance data with site data
WaterbirdGovernance <- merge(WaterbirdCountries, Governance, all.x=TRUE)
WaterbirdLostCountries <- WaterbirdGovernance[is.na(WaterbirdGovernance$Gov),]
print(unique(WaterbirdLostCountries$Country)) #We lose four countries: Isle of Man, Aland Islands, Western Sahara and Guernsey

names(WaterbirdGovernance)[names(WaterbirdGovernance)=="Year"] <- "YearClim"
WaterbirdGovernance$YearClim <- as.numeric(as.character(WaterbirdGovernance$YearClim))
WaterbirdGovernance <- WaterbirdGovernance[,c("SiteCode", "YearClim", "Season", "Gov")]

WaterbirdCounts <- merge(WaterbirdCounts, WaterbirdGovernance, by=c("SiteCode", "Season", "YearClim"), all.x=TRUE)
NumRows <- nrow(WaterbirdCounts)
WaterbirdCounts <- merge(WaterbirdCounts, WaterbirdCountries, by=c("SiteCode"))
if(nrow(WaterbirdCounts)!=NumRows){stop("You've lost data!")}

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDECountryGov.RData"))

#Governance values centred by year, is this ok? Or are they - we haven't taken to rankings but the absolute values, so maybe this is fine?

#### Travel distance and slope ####
#Resample to the same resolution as HYDE
Template <- raster(paste0(DataFP, "/HYDE/Anthrome/2017AD_anthromes/anthromes2017AD.asc"))
Template[!is.na(Template)] <- 1:length(Template[!is.na(Template)])

Travel <- raster(paste0(DataFP, "WorldPop/1km_mosaics/mastergrid_traveltime50k_1km.tif"))
Travel <- raster::resample(Travel, Template, method="bilinear")
crs(Travel) <- WGSCRS
writeRaster(Travel, paste0(DataFP, "WorldPop/TravelTimeResampled.tif"))

Slope <- raster(paste0(DataFP, "WorldPop/1km_mosaics/mastergrid_slope_1km.tif"))
Slope <- raster::resample(Slope, Template, method="bilinear")
crs(Slope) <- WGSCRS
writeRaster(Slope, paste0(DataFP, "WorldPop/SlopeResampled.tif"))

#SnapPoints
load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDECountryGov.RData"))

SitesSnap <- SnapPoints1("WorldPop", Travel, WaterbirdCounts)
WaterbirdCounts <- SnapPoints2("WorldPop", WaterbirdCounts, SitesSnap) 

SitesYearSeason <- unique(WaterbirdCounts[,c("SiteCode", "LonWorldPop", "LatWorldPop")])
SitesPoints <- cbind(SitesYearSeason$LonWorldPop, SitesYearSeason$LatWorldPop)
SitesPoints <- SpatialPointsDataFrame(SitesPoints, SitesYearSeason, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

SitesYearSeason$Travel <- raster::extract(Travel, SitesPoints, method="bilinear")
SitesYearSeason$Slope <- raster::extract(Slope, SitesPoints, method="bilinear")

WaterbirdCounts <- merge(WaterbirdCounts, SitesYearSeason[,c("SiteCode", "Travel", "Slope")], by="SiteCode")

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDECountryGovWorldPop.RData"))

#### Surface Water #### 
#If you need to redownload SW data, do it in cluster. Navigate to parent folder of surface water, create a folder to hold the download then run:
#python "/rds/user/hsw34/hpc-work/data/Global_Surface_Water_Change/downloadWaterData.py" "Extent" "extent"
#Ensuring downloadWaterData.py is in the right spot. Be careful with spaces, they're important.

#Aggregate surface water up to 5'
#Remember in the slurm_submit file on cluster, change to 12GB of CPU "#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem": "
#Extract values
#Rules for rasters: 
#East runs 0-170, rounds down (e.g. 0 = 0-10)
#West runs 10-180, rounds down (e.g. 10 = -10 - 0)
#North runs 0-80, rounds up (so 0 is 10South to 0)
#South runs 10-50, rounds up (so 10 is -20 to -10)

load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDECountryGovWorldPop.RData"))

Sites <- unique(WaterbirdCounts[,c("SiteCode", "Latitude", "Longitude")])

#Add codes so we can read the appropriate raster in (cos there are 500)
Sites$LatCode <- round_any(Sites$Latitude, 10, f=ceiling)
Sites$LatCode <- paste0(abs(Sites$LatCode), ifelse(Sites$LatCode<0, "S", "N"))
Sites$LonCode <- round_any(Sites$Longitude, 10, f=floor)
Sites$LonCode <- paste0(abs(Sites$LonCode), ifelse(Sites$LonCode<0, "W", "E"))
Sites$SWTile <- paste0(Sites$LonCode, "_", Sites$LatCode)
SWTiles <- unique(Sites$SWTile)
save(SWTiles, file=paste0(DataFP, "Global_Surface_Water_Change/SWTiles.RData"))

ClusterRun <- "Done"
if(ClusterRun!="Done"){
  load(file=paste0(DataFP, "Global_Surface_Water_Change/SWTiles.RData"))
  rasterOptions(maxmemory = 1e+11)
  lapply(c("Occurrence", "Extent"), function(Category){
    pblapply(SWTiles, function(x){
      if(file.exists(paste0(DataFP, "Global_Surface_Water_Change/", Category, "Aggregate/", tolower(Category), "_", x, "_v1_1.tif"))){return(NULL)}
      print(paste0(Category, x))
      writeRaster(raster(matrix(c(0,0),1,1)), paste0(DataFP, "Global_Surface_Water_Change/", Category, "Aggregate/", tolower(Category), "_", x, "_v1_1.tif"), format="GTiff") #Write a blank raster
      Water <- raster(paste0(DataFP, "Global_Surface_Water_Change/", Category, "/", tolower(Category), "_", x, "_v1_1.tif"))
      WaterReclass <- raster::reclassify(Water, matrix(c(254.5, 255.5, 0), ncol=3))
      AggregateFun <- ifelse(Category=="Occurrence", "mean", "sum")
      WaterAggregate <- raster::aggregate(WaterReclass, fact=333, fun=AggregateFun) #brings resolution to 0.08325 (i.e. 4.995 arc minutes)
      writeRaster(WaterAggregate, paste0(DataFP, "Global_Surface_Water_Change/", Category, "Aggregate/", tolower(Category), "_", x, "_v1_1.tif"), format="GTiff", overwrite=TRUE)
      print(paste0(Category, x, "DONE"))
      return(x)
    })
    return(Category)
  })
}

#Transfer Aggregate folders to local folder

#Extract values
SurfaceWater <- rbindlist(pbmclapply(SWTiles, function(x){
  SitesTile <- subset(Sites, SWTile==x)
  SitesPoints <- cbind(SitesTile$Longitude, SitesTile$Latitude)
  SitesPoints <- SpatialPointsDataFrame(SitesPoints, SitesTile, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  Occurrence <- raster(paste0(DataFP, "Global_Surface_Water_Change/OccurrenceAggregate/occurrence_", x, "_v1_1.tif"))
  Extent <- raster(paste0(DataFP, "Global_Surface_Water_Change/ExtentAggregate/extent_", x, "_v1_1.tif"))
  
  SitesPoints$SWOccurrence <- raster::extract(Occurrence, SitesPoints)
  SitesPoints$SWExtent <- raster::extract(Extent, SitesPoints, method="bilinear")
  
  return(as.data.frame(SitesPoints)[,c("SiteCode", "SWOccurrence", "SWExtent")])
}, mc.cores=ncores))
if(nrow(SurfaceWater)!=nrow(Sites)){stop("There are duplicates!")}

WaterbirdCounts2 <- merge(WaterbirdCounts, SurfaceWater, by="SiteCode")
if(nrow(WaterbirdCounts2)!=nrow(WaterbirdCounts)){stop("We've lost sites in the merge")}
WaterbirdCounts <- WaterbirdCounts2
save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDECountryGovWorldPopSurfaceWater.RData"))

#### Clean and finalise ####
load(file=paste0(DataFP, "/WaterbirdData_2020/Covariates/WaterbirdCounts_MigStatusClimateFertiliserHYDECountryGovWorldPopSurfaceWater.RData"))

#Remove snapped coordinates
WaterbirdCounts[,c("LatClim", "LonClim", "LatHYDE", "LonHYDE", "LatFert", "LonFert", "LatWorldPop", "LonWorldPop")] <- NULL
WaterbirdCounts[,c("SISRecID", "YearClim", "BLTaxCode")] <- NULL

#Rename HYDE (Should've done this earlier, ugh)
names(WaterbirdCounts)[names(WaterbirdCounts)=="rf_norice"] <- "RainfedNoRice"
names(WaterbirdCounts)[names(WaterbirdCounts)=="ir_norice"] <- "IrrigatedNoRice"
names(WaterbirdCounts)[names(WaterbirdCounts)=="tot_irri"] <- "TotalIrrigatedLand"
names(WaterbirdCounts)[names(WaterbirdCounts)=="tot_rainfed"] <- "TotalRainfedLand"
names(WaterbirdCounts)[names(WaterbirdCounts)=="tot_rice"] <- "TotalRice"
names(WaterbirdCounts)[names(WaterbirdCounts)=="cropland"] <- "Cropland"
names(WaterbirdCounts)[names(WaterbirdCounts)=="rangeland"] <- "Rangeland"
names(WaterbirdCounts)[names(WaterbirdCounts)=="rf_rice"] <- "RainfedRice"
names(WaterbirdCounts)[names(WaterbirdCounts)=="ir_rice"] <- "IrrigatedRice"
names(WaterbirdCounts)[names(WaterbirdCounts)=="conv_rangeland"] <- "ConvRangeland"
names(WaterbirdCounts)[names(WaterbirdCounts)=="grazing"] <- "Grazing"
names(WaterbirdCounts)[names(WaterbirdCounts)=="pasture"] <- "Pasture"

names(WaterbirdCounts)[names(WaterbirdCounts)=="popc_"] <- "PopulationCount"
names(WaterbirdCounts)[names(WaterbirdCounts)=="popd_"] <- "PopulationDensity"
names(WaterbirdCounts)[names(WaterbirdCounts)=="rurc_"] <- "RuralPopulationCount"
names(WaterbirdCounts)[names(WaterbirdCounts)=="urbc_"] <- "UrbanPopulationCount"
names(WaterbirdCounts)[names(WaterbirdCounts)=="uopp_"] <- "TotalBuiltupArea"

save(WaterbirdCounts, file=paste0(DataFP, "/WaterbirdData_2020/WaterbirdCounts_AllCovariates.RData"))
