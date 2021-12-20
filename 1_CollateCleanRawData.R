############################################################################################################################################################################
### This code was written by Hannah Wauchope to clean data for the paper "Protected areas have a mixed impact on waterbirds, but management helps"
### This is script 1 of 4 in the workflow
### Last edited 27th August, 2021
### Please direct queries to hannah.wauchope@gmail.com
###
### This script takes raw data from: The International Waterbird Census (From Wetlands International) and the Christmas Bird Count (from Audubon)
### Data from wetlands international was acquired via http://iwc.wetlands.org/index.php/requestingdata . Once an account has been created, it's easiest to get in touch with the Wetlands International data managers to request the full dataset, rather than obtaining from individual countries. This involves explaining the purpose for the data usage, and this is then sent out to National Coordinators for approval. Russian National Coordinators declined for Russian data to be used in this paper.
### I received the dataset used in this publication from Wetlands international on the 5th March, 2020.
### Data from Audubon was acquired via http://netapp.audubon.org/cbcobservation/ . Here, individual species or site data can be downloaded, but in the case of this study it's easier to email Audubon (using the email on the website) to request the full dataset. This again involves explaining the purpose for the study and acquiring permissions. 
### I received the dataset used in this publication from Audubon on 9th August, 2019.

### This script cleans the data by taxonomically standardising species names by the currently used names in Birdlife International, sorts out cases of multiple counts, removes any sites with suspicious coordinats, and imputes zeroes
###
### We thank the coordinators, thousands of volunteer counters, and funders of the International Waterbird Census. 
### CBC Data is provided by National Audubon Society and through the generous efforts of Bird Studies Canada and countless volunteers across the western hemisphere to whom we are most grateful
############################################################################################################################################################################

#### Initialise ####
library(sp)
library(rgdal)
library(ggplot2)
library(raster)
library(data.table)
library(pbmcapply)
library(reshape2)
library(rgeos)
library(plyr)
library(ncdf4)
library(pbapply)
library(maps)
library(ggalt)
library(taxize)
library(stringr)
library(rredlist)
library(dplyr)

DataFP <- "/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Data/"

#Sort out taxonomy and shapefiles
Tax <- fread(paste0(DataFP, "Taxonomy/HBW_BLI_v4.csv"))
Tax <- Tax[,c("Order", "Family name", "CommonName", "ScientificName", "BirdLifeTaxonomicTreatment", "IUCN2019", "SISRecID")]
names(Tax) <- c("Order", "Family", "CommonName", "Species", "BLTaxCode", "IUCN", "SISRecID")
Tax$Genus <- str_split_fixed(Tax$Species, "[ ]",2)[,1] #Get genus
Tax$Order <- sapply(Tax$Order, function (x) paste(toupper(substring(tolower(x), 1,1)), substring(tolower(x), 2), sep="", collapse=" ")) #Change order to lower case(with capital first letter)

#BirdPolygons1 <- readOGR(paste0(DataFP, "BirdlifePolygons/BOTW_Shape/BOTW_Shape_1_3000/BOTW_Shape_1_3000.shp"))
#BirdPolygons2 <- readOGR(paste0(DataFP, "BirdlifePolygons/BOTW_Shape/BOTW_Shape_3001_10000/BOTW_Shape_3001_10000.shp"))
#BirdPolygons3 <- readOGR(paste0(DataFP, "BirdlifePolygons/BOTW_Shape/BOTW_Shape_10001_17463/BOTW_Shape_10001_17463.shp"))
BOTWNames <- fread(paste0(DataFP, "BirdlifePolygons/BOTW_Attributes.csv"))

#### Functions ####
TreeofLifeCleaning <- function(TreeofLifeID, SpeciesList){
  #Pull out those with only one match, get all synonyms
  TreeofLifeIDOneMatch <- TreeofLifeID[lapply(TreeofLifeID, nrow)==1]
  TreeofLifeIDOneMatch <- rbindlist(lapply(1:length(TreeofLifeIDOneMatch), function(x){
    TOL <- TreeofLifeIDOneMatch[[x]]
    TOL$CountsSpecies <- names(TreeofLifeIDOneMatch)[x]
    names(TOL) <- c("UniqueName", "MatchedName", "IsSynonym", "Score", "NomenclatureCode", "IsApproxMatch", "Flags", "OTTID", "Rank", "CountsSpecies")
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
    names(TOL) <- c("UniqueName", "MatchedName", "IsSynonym", "Score", "NomenclatureCode", "IsApproxMatch", "Flags", "OTTID", "Rank", "CountsSpecies")
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
    if(x["Species"] %in% BOTWNames$SCINAME){
      x["Species"]
    } else if(x["Synonym_1"] %in% BOTWNames$SCINAME){
      x["Synonym_1"]
    } else if(x["Synonym_2"] %in% BOTWNames$SCINAME){
      x["Synonym_2"]
    } else {"NoMatch"}
  })
  
  return(SpeciesSynonyms)
}
PostUpdateTaxonCleaning <- function(NameUpdates, SpeciesSynonyms){
  SpeciesSynonyms2 <- as.data.frame(rbind(subset(SpeciesSynonyms, FinalName!="NoMatch"), NameUpdates, fill=T))
  SpeciesSynonyms2 <- SpeciesSynonyms2[,c(1, 4:ncol(SpeciesSynonyms2))]
  
  #Now let's check again to see if any more are unmatched
  SpeciesSynonymsCheck <- subset(SpeciesSynonyms2, NameChange1Hybrid2!=2)
  SpeciesSynonymsCheck$FinalName <- apply(SpeciesSynonymsCheck, 1, function (x){
    if(x["FinalName"] %in% BOTWNames$SCINAME){
      x["FinalName"]
    } else {"NoMatch"}
  })
  
  if(nrow(subset(SpeciesSynonymsCheck, FinalName=="NoMatch"))!=0){stop("some species still aren't matched!")}
  
  #Sort out the hybrids
  SpeciesSynonyms2[is.na(SpeciesSynonyms2$NameChange1Hybrid2),]$NameChange1Hybrid2 <- 1
  SpeciesSynonymsHybrids <- subset(SpeciesSynonyms2, NameChange1Hybrid2==2) %>% mutate_all(as.character)
  SpeciesSynonymsHybrids <- melt(SpeciesSynonymsHybrids, id.vars=c("Species", "FinalName", "NameChange1Hybrid2"))
  SpeciesSynonymsHybrids <- subset(SpeciesSynonymsHybrids, value!="")
  names(SpeciesSynonymsHybrids) <- c("GivenSpecies", "Species", "HybridFlag", "HybridLevel", "HybridValue")
  SpeciesSynonymsHybrids$HybridFlag <- 1
  SpeciesSynonymsHybrids$Species <- NA
  
  SpeciesSynonymsNotHybrid <- subset(SpeciesSynonyms2, NameChange1Hybrid2==1)[,c(1:3)]
  names(SpeciesSynonymsNotHybrid) <- c("GivenSpecies", "Species", "HybridFlag")
  SpeciesSynonymsNotHybrid$HybridFlag <- 0
  
  SpeciesSynonymsFinal <- data.table::rbindlist(list(SpeciesSynonymsHybrids, SpeciesSynonymsNotHybrid), fill=TRUE)
  return(SpeciesSynonymsFinal)
}
Map <- function(Coordinates){
  mapWorld <- borders("world", colour="gray70", fill="gray70")
  
  Coordinates <- subset(CBCDistAll, Cat=="Lost")
  ggplot() + mapWorld +
    geom_point(aes(x=Coordinates$Longitude, y=Coordinates$Latitude), size=0.6, colour="#0B4EA2", shape=19, alpha=0.9)+
    theme(aspect.ratio=0.53, panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(hjust = 0), legend.justification = "top", axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(), panel.border = element_rect(size = 1, fill = NA, colour="black"),
          plot.background=element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent"))
}

#### IWC Initial Clean ####
IWCCounts <- fread(paste0(DataFP, "WaterbirdData_2020/IWC/IWC170220.csv"))

IWCCountsPoints <- unique(IWCCounts[,c("countrycode", "sitename", "site", "x", "y")])

IWCCountsPoints$Latitude <- IWCCountsPoints$y
IWCCountsPoints$Longitude <- IWCCountsPoints$x
IWCCountsPoints[,c("x", "y")] <- NULL
Map(IWCCountsPoints)

#There are a handful of sites where the lat and long entries are flipped! Flip back
IWCCountsPointsWrong <- subset(IWCCountsPoints, Latitude>72)
names(IWCCountsPointsWrong) <- c("countrycode", "sitename", "site", "Longitude", "Latitude")
IWCCountsPointsWrong <- IWCCountsPointsWrong[,c(1,2,3,5,4)]

Map(IWCCountsPointsWrong)

IWCCountsWrong <- subset(IWCCounts, y>72)
names(IWCCountsWrong) <- c("iwccountry", "countrycode", "site", "sitename", "y", "x", "geometrytype", "species", "reporting_name", "count", "visit_date")
IWCCountsWrong <- IWCCountsWrong[,c("iwccountry", "countrycode", "site", "sitename", "x","y", "geometrytype", "species", "reporting_name", "count", "visit_date")]

IWCCounts <- subset(IWCCounts, y<72)
IWCCounts <- rbind(IWCCounts, IWCCountsWrong)

#Update names
names(IWCCounts) <- c("IWCCountry", "CountryCode", "SiteCode", "SiteName", "Longitude", "Latitude", "Geometry", "SpeciesCode", "Species", "Count", "visit_date")

#GetDate
IWCCounts$Year <- sapply(as.character(IWCCounts$visit_date), function (x) strsplit(x,c("[-]+"))[[1]][1])
IWCCounts$Month <- sapply(as.character(IWCCounts$visit_date), function (x) strsplit(x,c("[-]+"))[[1]][2])
IWCCounts$Day <- sapply(as.character(IWCCounts$visit_date), function (x) strsplit(x,c("[-]+"))[[1]][3])
IWCCounts$visit_date <- NULL

#Convert Dates to Seasons
IWCCounts$Month <- as.numeric(as.character(IWCCounts$Month))
IWCCounts$Year <- as.numeric(as.character(IWCCounts$Year))
IWCCounts$Day <- as.numeric(as.character(IWCCounts$Day))

IWCCounts <- IWCCounts[IWCCounts$Month %in% c(12,1,2,6,7,8),] #Only take strict summer and winter months
IWCCounts[(IWCCounts$Month==12)]$Year <- IWCCounts[(IWCCounts$Month==12)]$Year+1 #Flip counts in December to the year after

IWCCounts$Season <- "DectoFeb"
IWCCounts[(IWCCounts$Month==6 | IWCCounts$Month==7 | IWCCounts$Month==8)]$Season <- "JuntoAug"
if(unique(subset(IWCCounts, Month==12 | Month==1 | Month==2)$Season)!="DectoFeb"){stop("Something's messed up with seasons")}
if(unique(subset(IWCCounts, Month==6 | Month==7 | Month==8)$Season)!="JuntoAug"){stop("Something's messed up with seasons")}

#### CBC Initial Clean ####
#Get counts
Years1_50 <- fread(paste0(DataFP, "WaterbirdData_2020/CBC/RawData/1-50-all-CBC_Circle_Species_Report_SQL_updated.csv"))
Years117_118 <- fread(paste0(DataFP, "WaterbirdData_2020/CBC/RawData/117-118-CBC_Circle_Species_Report_SQL_updated.csv"))
CBCCounts <- rbindlist(lapply(list.files(paste0(DataFP, "WaterbirdData_2020/CBC/RawData/"), pattern="all-CBC_Circle_Species_Report.txt$", full.names=TRUE), function(x){
  print(x)
  return(read.delim(x, fill=TRUE, header=TRUE, quote=""))
}))
names(CBCCounts) <- c("Abbrev", "Name", "Latitude", "Longitude", "Subnational_code", "Country_code", "Count_yr", "Cnt_dt", "COM_NAME", "SCI_NAME", "how_many", "TotalSpecies", "Editor_comment", "SORT_CBC")
CBCCounts$Abbrev <- gsub('"', '', as.character(CBCCounts$Abbrev))
CBCCounts$Subnational_code <- gsub('"', '', as.character(CBCCounts$Subnational_code))
CBCCounts$Country_code <- gsub('"', '', as.character(CBCCounts$Country_code))
CBCCounts$Cnt_dt <- gsub('"', '', as.character(CBCCounts$Cnt_dt))
CBCCounts$COM_NAME <- gsub('"', '', as.character(CBCCounts$COM_NAME))
CBCCounts$Name <- gsub('"', '', as.character(CBCCounts$Name))
CBCCounts$SCI_NAME <- gsub('"', '', as.character(CBCCounts$SCI_NAME))

CBCCounts <- rbind(CBCCounts, Years1_50)
CBCCounts <- rbind(CBCCounts, Years117_118)
CBCCounts$SORT_CBC <- as.numeric(gsub(',', '', as.character(CBCCounts$SORT_CBC)))

#### Taxonomic Cleaning IWC ####
## Rule for choosing correct species name: match to birdlife shapefiles
SpeciesList <- as.data.frame(unique(IWCCounts$Species))
names(SpeciesList) <- c("Species")

#The aim here is to get all possible synonyms for each species
#Get alternative species names using taxize
TreeofLifeID <- pbmclapply(SpeciesList$Species, function (x) tryCatch(as.data.frame(get_tolid_(as.character(x))), error=function(e){"NoMatch"}), mc.cores=8)
names(TreeofLifeID) <- SpeciesList$Species

#Get any synonyms, cross check all names with Birds of the World, return a list of either the correct name (with original name + matched name, which will either be the same or a synonym), or just the incorrect name with "nomatch"
SpeciesSynonyms <- TreeofLifeCleaning(TreeofLifeID, SpeciesList)

#Write out the "no matches" to manually clean
write.csv(subset(SpeciesSynonyms, FinalName=="NoMatch"), paste0(DataFP, "WaterbirdData_2020/IWC/Cleaning/NoMatchSpecies.csv"), row.names=FALSE)
NameUpdates <- read.csv(paste0(DataFP, "WaterbirdData_2020/IWC/Cleaning/NoMatchSpeciesResearchedIWC.csv"))

#Produce a dataframe with EITHER original name + matched name, or original name + higher level taxonomy if no species ID available. This will be used for zero imputation later
IWCTaxonCleaned <- PostUpdateTaxonCleaning(NameUpdates, SpeciesSynonyms) #Ignore warning

#Add in higher order taxonomy
GetSISID <- merge(IWCTaxonCleaned, unique(BOTWNames[,c("SCINAME", "SISRecID")]), by.x="Species", by.y="SCINAME", all.x=TRUE)
if(nrow(GetSISID) != nrow(IWCTaxonCleaned)){stop("Matching issue!")}
if(nrow(GetSISID[!is.na(GetSISID$SISRecID),]) != nrow(subset(IWCTaxonCleaned, HybridFlag==0))){stop("Matching issue!")}

GetSISID2 <- as.data.frame(merge(GetSISID, Tax, by=c("SISRecID", "Species"), all.x=TRUE))
if(nrow(GetSISID2) != nrow(IWCTaxonCleaned)){stop("Matching issue!")}
if(nrow(GetSISID2[!is.na(GetSISID2$SISRecID),]) != nrow(subset(IWCTaxonCleaned, HybridFlag==0))){stop("Matching issue!")}

GetSISID2$Genus <- NA
GetSISID2[GetSISID2$HybridFlag==0,]$Genus <- str_split_fixed(GetSISID2[GetSISID2$HybridFlag==0,]$Species, "[ ]", 2)[,1]

###Higher order matching check###
HybridGenus <- unique(subset(GetSISID2, HybridLevel=="Genus")$HybridValue)
MissingGenus <- HybridGenus[!HybridGenus %in% unique(Tax$Genus)]

HybridFamily <- unique(subset(GetSISID2, HybridLevel=="Family")$HybridValue)
HybridFamily <- HybridFamily[!HybridFamily %in% unique(Tax$Family)]

HybridOrder <- unique(subset(GetSISID2, HybridLevel=="Order")$HybridValue)
HybridOrder <- HybridOrder[!HybridOrder %in% unique(Tax$Order)]

IWCCounts2 <- merge(GetSISID2, IWCCounts, by.x="GivenSpecies", by.y="Species")
if(nrow(IWCCounts2)!=nrow(IWCCounts)){stop("Something's gone wrong!")}
write.csv(IWCCounts2, paste0(DataFP, "WaterbirdData_2020/IWC/IWCCounts_SpeciesCleaned.csv"), row.names = FALSE)

#### Taxonomic Cleaning CBC ####
CBCSpecies <- as.data.frame(unique(CBCCounts$SCI_NAME))
names(CBCSpecies) <- "Species"

CBCSpecies$FinalName <- apply(CBCSpecies, 1, function (x){
  if(x["Species"] %in% BOTWNames$SCINAME){
    x["Species"]
  } else {"NoMatch"}
})

CBCNoMatch <- subset(CBCSpecies, FinalName=="NoMatch")

#TreeofLifeID <- pbmclapply(CBCNoMatch$Species, function (x) tryCatch(as.data.frame(get_tolid_(as.character(x))), error=function(e){"NoMatch"}), mc.cores=8)
#names(TreeofLifeID) <- CBCNoMatch$Species
#save(TreeofLifeID, file=paste0(DataFP, "WaterbirdData_2020/CBC/TreeofLifeID.RData"))
load(file=paste0(DataFP, "WaterbirdData_2020/CBC/TreeofLifeID.RData"))

SpeciesSynonyms <- TreeofLifeCleaning(TreeofLifeID, CBCNoMatch)
write.csv(subset(SpeciesSynonyms, FinalName=="NoMatch"), paste0(DataFP, "WaterbirdData_2020/CBC/Cleaning/NoMatchSpecies.csv"), row.names=FALSE)

#Load up cleaned species
NameUpdates <- read.csv(paste0(DataFP, "WaterbirdData_2020/CBC/Cleaning/NoMatchSpeciesResearchedCBC.csv"))

CBCTaxonCleaned <- PostUpdateTaxonCleaning(NameUpdates, SpeciesSynonyms)

#Now add back into the CBC that was cleaned to begin with:
names(CBCSpecies) <- c("GivenSpecies", "Species")

CBCTaxonCleaned <- rbind(subset(CBCSpecies, Species!="NoMatch"), CBCTaxonCleaned, fill=TRUE)
if(nrow(CBCTaxonCleaned)!=nrow(CBCSpecies)){stop("Somethings gone wrong with synonym matching!")}
CBCTaxonCleaned[is.na(CBCTaxonCleaned$HybridFlag)]$HybridFlag <- 0

#Add in higher order taxonomy
GetSISID <- merge(CBCTaxonCleaned, unique(BOTWNames[,c("SCINAME", "SISRecID")]), by.x="Species", by.y="SCINAME", all.x=TRUE)
if(nrow(GetSISID) != nrow(CBCTaxonCleaned)){stop("Matching issue!")}
if(nrow(GetSISID[!is.na(GetSISID$SISRecID),]) != nrow(subset(CBCTaxonCleaned, HybridFlag==0))){stop("Matching issue!")}

GetSISID2 <- as.data.frame(merge(GetSISID, Tax, by=c("SISRecID", "Species"), all.x=TRUE))
if(nrow(GetSISID2) != nrow(CBCTaxonCleaned)){stop("Matching issue!")}
if(nrow(GetSISID2[!is.na(GetSISID2$SISRecID),]) != nrow(subset(CBCTaxonCleaned, HybridFlag==0))){stop("Matching issue!")}

GetSISID2$Genus <- NA
GetSISID2[GetSISID2$HybridFlag==0,]$Genus <- str_split_fixed(GetSISID2[GetSISID2$HybridFlag==0,]$Species, "[ ]", 2)[,1]

###Higher order matching check###
HybridGenus <- unique(subset(GetSISID2, HybridLevel=="Genus")$HybridValue)
HybridGenus <- HybridGenus[!HybridGenus %in% unique(Tax$Genus)]

HybridFamily <- unique(subset(GetSISID2, HybridLevel=="Family")$HybridValue)
HybridFamily <- HybridFamily[!HybridFamily %in% unique(Tax$Family)]

HybridOrder <- unique(subset(GetSISID2, HybridLevel=="Order")$HybridValue)
HybridOrder <- HybridOrder[!HybridOrder %in% unique(Tax$Order)]

#Add back to count data
CBCCounts2 <- merge(GetSISID2, CBCCounts, by.x="GivenSpecies", by.y="SCI_NAME", all=T)
if(nrow(CBCCounts2)!=nrow(CBCCounts)){stop("Something's gone wrong!")}
write.csv(CBCCounts2, paste0(DataFP, "WaterbirdData_2020/CBC/CBCCounts_SpeciesCleaned.csv"), row.names = FALSE)

#### CBC Effort Cleaning + Double Counts ####
#Read in count data
CBCCounts <- fread(paste0(DataFP, "WaterbirdData_2020/CBC/CBCCounts_SpeciesCleaned.csv"))

#Get Effort
Effort1_116 <- fread(paste0(DataFP, "WaterbirdData_2020/CBC/RawData/1-116-all-CBC_Effort_Report_2.csv"))
names(Effort1_116) <- c("Abbrev", "Name", "Count_yr", "Country_code", "Distance", "Hours")
Effort117_118 <- fread(paste0(DataFP, "WaterbirdData_2020/CBC/RawData/117-118-CBC_Effort_Report_SQL_updated-2.csv"))
names(Effort117_118) <- c("Abbrev", "Name", "Count_yr", "Country_code", "Method", "Distance", "Units", "Hours")
CBCEffort <- rbind(Effort1_116, Effort117_118, fill=TRUE)

#Distance and hours are correlated
CBCEffort$Distance <- as.numeric(CBCEffort$Distance)
CBCEffort$Hours <- as.numeric(CBCEffort$Hours)
CBCEffort <- CBCEffort[complete.cases(CBCEffort$Hours),]
CBCEffort <- CBCEffort[complete.cases(CBCEffort$Distance),]
CBCEffortThou <- CBCEffort[sample(1:nrow(CBCEffort), 1000, replace=FALSE),]
plot(CBCEffortThou$Hours, CBCEffortThou$Distance)

#Remove distance and zero count hours
CBCEffort[,c("Distance", "Units", "Method", "Name", "Country_code")] <- NULL
CBCEffort <- subset(CBCEffort, Hours!=0)

#Now combine hours
keys <- colnames(CBCEffort)[!colnames(CBCEffort) %in% "Hours"]
CBCEffort2 <- CBCEffort[,list(Hours=sum(Hours)),keys]

Check <- dcast(CBCEffort2, Abbrev + Count_yr~., length, value.var="Hours")
if(max(unique(Check$.))!=1){stop("There are still cases of double hours!")}

CBCCounts$SiteYear <- paste0(CBCCounts$Abbrev,"_", CBCCounts$Count_yr)
CBCEffort2$SiteYear <- paste0(CBCEffort2$Abbrev,"_", CBCEffort2$Count_yr)
CBCEffort2[,c("Abbrev", "Count_yr")] <- NULL

CBCCounts <- merge(CBCCounts, CBCEffort2, by="SiteYear")

#Now sort out dates
CBCCounts[,c("Year", "Month", "Day")] <- as.data.frame(str_split_fixed(CBCCounts$Cnt_dt, "[-]", 3))
CBCCounts$Month <- as.numeric(as.character(CBCCounts$Month))
CBCCounts[is.na(CBCCounts$Month),]$Month <- 0
CBCCounts <- subset(CBCCounts, Month!=3) #Remove march counts (only Dec, Jan, Feb)

CBCCounts$Year <- 1900+CBCCounts$Count_yr

CBCCounts$Day <- as.numeric(as.character(CBCCounts$Day))
CBCCounts[is.na(CBCCounts$Day),]$Day <- 0

CBCCounts[,c("Count_yr", "Cnt_dt")] <- NULL

#Now sort out double counts (sum by days counted, then mean by year)
CBCCounts$SiteSpecYear <- paste0(CBCCounts$Abbrev, "_", CBCCounts$Species, "_", CBCCounts$Year)
CBCCountsNoHybrids <- subset(CBCCounts, HybridFlag==0)

#Sum by date
CBCCountsNoHybrids[,c("COM_NAME", "GivenSpecies", "SORT_CBC", "Editor_comment", "TotalSpecies")] <- NULL
keys <- colnames(CBCCountsNoHybrids)[!colnames(CBCCountsNoHybrids) %in% "how_many"]
CBCCountsNoHybrids2 <- CBCCountsNoHybrids[,list(Count=sum(how_many)),keys]

#Mean by year
CBCCountsNoHybrids2[,c("Month", "Day")] <- NULL

keys <- colnames(CBCCountsNoHybrids2)[!colnames(CBCCountsNoHybrids2) %in% "Count"]
CBCCountsNoHybrids3 <- CBCCountsNoHybrids2[,list(Count=round(mean(Count))),keys]

#Check
DoubleCountCheck <- dcast(CBCCountsNoHybrids3, SiteSpecYear~., length, value.var="Count")
if(max(unique(DoubleCountCheck$.))!=1){stop("There are still double counts!")}

#Recombine Data
CBCCountsHybrids <- subset(CBCCounts, HybridFlag==1)
CBCCountsHybrids[,c("COM_NAME", "GivenSpecies", "SORT_CBC", "Editor_comment", "TotalSpecies", "Month", "Day", "how_many")] <- NULL
CBCCountsHybrids$Count <- NA

CBCCounts <- rbind(CBCCountsHybrids, CBCCountsNoHybrids3)
CBCCounts[,c("SiteSpecYear", "SiteYear")] <- NULL

#Remove multiple sites at the site coordinate:
CBCCounts$Latitude <- round(CBCCounts$Latitude, 4)
CBCCounts$Longitude <- round(CBCCounts$Longitude, 4)

CBCCounts$Coordinates <- paste0(CBCCounts$Latitude, "_", CBCCounts$Longitude)
CBCCountsCast <- dcast(unique(CBCCounts[,c("Coordinates", "Abbrev")]), Coordinates~., length, value.var="Abbrev")

CBCCounts <- CBCCounts[CBCCounts$Coordinates %in% subset(CBCCountsCast, .==1)$Coordinates,]
CBCCounts$Coordinates <- NULL

#Check on a map
Map(unique(CBCCounts[,c("Latitude", "Longitude", "Abbrev")]))

#Lovely
write.csv(CBCCounts, paste0(DataFP, "WaterbirdData_2020/CBC/CBCCounts_SpeciesEffortCleaned.csv"), row.names=FALSE)

#### Location Cleaning IWC ####
IWCCounts <- fread(paste0(DataFP, "WaterbirdData_2020/IWC/IWCCounts_SpeciesCleaned.csv"))

#Remove sites with no coordinates
IWCCounts <- subset(IWCCounts, Longitude!=0 & Latitude!=0)

#Remove error sites that Tatsuya found (in process of Nature paper: "Successful conservation of global waterbird populations depends on effective governance" 2017) (error by very bogus coordinates). WI have cleaned a lot of these now! So identify the few unchanged ones and remove those
ErrorSites <- read.csv("/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Data/WaterbirdData_Tatsuya/ErrorsRemoved/SitesToBeRemoved.csv")
ErrorSites <- ErrorSites[,c("SiteName", "SiteLon", "SiteLat")]
names(ErrorSites) <- c("SiteCode", "LongitudeError", "LatitudeError")

ErrorSites <- ErrorSites[ErrorSites$SiteCode %in% IWCCounts$SiteCode,]
Check <- unique(IWCCounts[IWCCounts$SiteCode %in% ErrorSites$SiteCode,c("SiteCode", "Latitude", "Longitude")])
ErrorSitesMerge <- merge(ErrorSites, Check, by="SiteCode")
ErrorSitesMerge <- subset(ErrorSitesMerge, LongitudeError==Longitude | LatitudeError == Latitude)

IWCCounts <- IWCCounts[!IWCCounts$SiteCode %in% ErrorSitesMerge$SiteCode]
IWCCounts$Latitude <- round(IWCCounts$Latitude, 4)
IWCCounts$Longitude <- round(IWCCounts$Longitude, 4)

###Now let's remove duplicate sites (i.e. multiple site codes with the same coordinate)
IWCSites <- unique(IWCCounts[,c("SiteCode", "SiteName", "Longitude", "Latitude")]) #Start with 46053 sites
IWCSites$Coordinates <- paste0(IWCSites$Latitude, ".", IWCSites$Longitude)
IWCSitesDuplCast <- dcast(IWCSites, Coordinates ~., length, value.var="SiteCode")
IWCSitesDuplCast <- subset(IWCSitesDuplCast, .>1)

IWCSitesDupl <- IWCSites[IWCSites$Coordinates %in% IWCSitesDuplCast$Coordinates,] #4162 sites
Map(IWCSites)
Map(IWCSitesDupl) 
IWCSites <- IWCSites[!IWCSites$SiteCode %in% IWCSitesDupl$SiteCode,]
Map(IWCSites) #Now 41,891 sites
IWCCounts <- IWCCounts[IWCCounts$SiteCode %in% IWCSites$SiteCode,]

#What about multiple counts per species per year
#First, remove hybrid counts
IWCCounts <- as.data.table(IWCCounts)

IWCCountsNoHybrids <- subset(IWCCounts, HybridFlag==0)
IWCCountsNoHybrids$SiteSpecSeasonYear <- paste0(IWCCountsNoHybrids$SiteCode, "_", IWCCountsNoHybrids$Species, "_", IWCCountsNoHybrids$Season, "_", IWCCountsNoHybrids$Year)

#Sum by day
IWCCountsNoHybrids[,c("SpeciesCode", "GivenSpecies")] <- NULL

keys <- colnames(IWCCountsNoHybrids)[!colnames(IWCCountsNoHybrids) %in% "Count"]
IWCCountsNoHybrids <- IWCCountsNoHybrids[,list(Count=sum(Count)),keys]

#Mean by year
IWCCountsNoHybrids[,c("Month", "Day")] <- NULL
keys <- colnames(IWCCountsNoHybrids)[!colnames(IWCCountsNoHybrids) %in% "Count"]
IWCCountsNoHybrids2 <- IWCCountsNoHybrids[,list(Count=round(mean(Count))),keys]

IWCCountsDuplCast <- dcast(IWCCountsNoHybrids2, SiteSpecSeasonYear~., length, value.var="Count")
if(max(IWCCountsDuplCast$.)>1){stop("There are still duplicate counts!")}
IWCCountsNoHybrids2$SiteSpecSeasonYear <- NULL

IWCCounts <- as.data.frame(subset(IWCCounts, HybridFlag==1))
IWCCounts[,c("SpeciesCode", "GivenSpecies", "Month", "Day")] <- NULL

IWCCounts <- IWCCounts[,c(names(IWCCounts)[!names(IWCCounts) %in% "Count"], "Count")]
IWCCounts <- rbind(IWCCounts, IWCCountsNoHybrids2)

write.csv(IWCCounts, paste0(DataFP, "WaterbirdData_2020/IWC/IWCCounts_LocationAndCountsCleaned.csv"), row.names=FALSE)

#### Combine Datasets, Impute Zeroes #### 
IWCCounts <- fread(paste0(DataFP, "WaterbirdData_2020/IWC/IWCCounts_LocationAndCountsCleaned.csv"))
CBCCounts <- fread(paste0(DataFP, "WaterbirdData_2020/CBC/CBCCounts_SpeciesEffortCleaned.csv"))

IWCCounts[IWCCounts$Count<0 & HybridFlag==0,]$HybridLevel <- "Species"
IWCCounts[IWCCounts$Count<0 & HybridFlag==0,]$HybridValue <- IWCCounts[IWCCounts$Count<0 & HybridFlag==0,]$Species
IWCCounts[IWCCounts$Count<0 & HybridFlag==0,]$Species <- NA
IWCCounts[IWCCounts$Count<0 & HybridFlag==0,]$HybridFlag <- 1

names(IWCCounts)
names(CBCCounts)

#Align names
names(CBCCounts)[names(CBCCounts) %in% "Abbrev"] <- "SiteCode"
names(CBCCounts)[names(CBCCounts) %in% "Country_code"] <- "CountryCode"
names(CBCCounts)[names(CBCCounts) %in% "Subnational_code"] <- "SubnationalCode"
names(CBCCounts)[names(CBCCounts) %in% "Name"] <- "SiteName"
CBCCounts$Season <- "DectoFeb"

names(IWCCounts)[!names(IWCCounts) %in% names(CBCCounts)]
names(CBCCounts)[!names(CBCCounts) %in% names(IWCCounts)]

IWCCounts$Dataset <- "IWC"
CBCCounts$Dataset <- "CBC"

BirdCounts <- rbind(CBCCounts, IWCCounts, fill=TRUE)

BirdCountsZeroImpute <- rbindlist(pbmclapply(unique(BirdCounts$SiteCode), function(x){
  SiteCounts <- subset(BirdCounts, SiteCode==x & HybridFlag==0)
  SiteCounts$ImputeFlag <- 0
  SiteYears <- unique(SiteCounts$Year)
  SiteSpecies <- unique(SiteCounts$Species)
  SiteHybrids <- subset(BirdCounts, SiteCode==x & HybridFlag==1)
  SiteHybridsSpecies <- unique(SiteHybrids[,c("HybridLevel", "HybridValue", "Year")])
  ZeroesSiteSpecies <- rbindlist(lapply(SiteSpecies, function(y){
    SiteSpeciesCounts <- subset(SiteCounts, HybridFlag==0 & Species==y)
    HybridMatch <- SiteHybridsSpecies[SiteHybridsSpecies$HybridValue %in% c(unique(SiteSpeciesCounts$Genus), unique(SiteSpeciesCounts$Family), unique(SiteSpeciesCounts$Order), unique(SiteSpeciesCounts$Class))]
    ZeroYears <- SiteYears[!SiteYears %in% unique(c(SiteSpeciesCounts$Year, HybridMatch$Year))]
    ZeroCounts <- do.call("rbind", replicate(length(ZeroYears), SiteSpeciesCounts[1,], simplify = FALSE))
    if(is.null(ZeroCounts)){
      return(SiteSpeciesCounts)
    }
    ZeroCounts$Count <- 0
    ZeroCounts$Year <- ZeroYears
    ZeroCounts$ImputeFlag <- 1
    SiteSpeciesCounts <- rbind(SiteSpeciesCounts, ZeroCounts)
    return(SiteSpeciesCounts)
  }))
  write.csv(ZeroesSiteSpecies, paste0(DataFP, "WaterbirdData_2020/ZeroImputedSites/Site", x, ".csv"), row.names=FALSE)
  return(ZeroesSiteSpecies)
}, mc.cores=6))

#Some sites are lost, but only those without any counted species (either only hybrids, or only negative obs)
BirdCountsZeroImpute <- rbindlist(pblapply(list.files(path=paste0(DataFP, "WaterbirdData_2020/ZeroImputedSites/"), pattern="*.csv", full.names = TRUE), function(x) fread(x)), fill=TRUE)
BirdCountsZeroImpute[,c("HybridFlag", "HybridLevel", "HybridValue")] <- NULL

write.csv(BirdCountsZeroImpute, paste0(DataFP, "WaterbirdData_2020/BirdCounts.csv"), row.names=FALSE)


