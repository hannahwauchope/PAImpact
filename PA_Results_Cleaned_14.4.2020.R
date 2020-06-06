#### Initialise ####
cluster <- FALSE
if(cluster==TRUE){
  library(data.table, lib.loc="/home/hsw34/my-R-libs/")
  library(rgdal, lib.loc="/home/hsw34/my-R-libs/")
  library(sp, lib.loc="/home/hsw34/my-R-libs/")
  library(pbapply, lib.loc="/home/hsw34/my-R-libs/")
  library(tidyverse, lib.loc="/home/hsw34/my-R-libs/")
  library(MASS, lib.loc="/home/hsw34/my-R-libs/")
  library(plyr, lib.loc="/home/hsw34/my-R-libs/")
  library(pbmcapply, lib.loc="/home/hsw34/my-R-libs/")
  library(MASS, lib.loc="/home/hsw34/my-R-libs/")
  library(usdm, lib.loc="/home/hsw34/my-R-libs/")
  library(StatMatch, lib.loc="/home/hsw34/my-R-libs/")
  library(rlist, lib.loc="/home/hsw34/my-R-libs/")
  library(resample, lib.loc="/home/hsw34/my-R-libs/")
  library(foreign, lib.loc="/home/hsw34/my-R-libs/")
  library(lme4, lib.loc="/home/hsw34/my-R-libs/")
  library(scales, lib.loc="/home/hsw34/my-R-libs/")
  
  ResultsFP <- "/rds/user/hsw34/hpc-work/results/c3/"
  DataFP <- "/rds/user/hsw34/hpc-work/data/"
  ncores <- 32 
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
  library(glmmTMB) #not on cluster
  library(ggalt) #not on cluster
  
  ncores <- 4
  
  ###File paths
  DataFP <- "/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Data/"
  ResultsFP <- "/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Chapter3/Analysis/"
  FiguresFP <-"/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Chapter3/Figures/Results/"
  
}

#99, 10, 7, zeroes imputed, 0.05, 0.25, 0.8. Redlist didn't run
#ERRORS: 7, 7, 5, 3, St Diff Thresh 0.1 = BA spec match lost
#77, 10, 76, zeroies imputed, parallel P 0.01, 01 = Mean Slope Models Failed

#### Define functions ####
#Round to degrees of freedom
MollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
WGSCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plotfontsize <- 12
round_df <- function(x, digits) {
  x <- as.data.frame(x)
  numeric_columns <- sapply(x, class) == 'numeric'
  x[numeric_columns] <-  round(x[,numeric_columns], digits)
  x
}
TheNumbers <- function(Dataset){
  Dataset <- as.data.frame(Dataset)
  sapply(c('SiteCode', 'Species', 'SiteSpec', 'Database'), function(x) length(unique(Dataset[,c(x)])))
}
GetMeanCovs <- function(dataset, VariablesForMatchingByYear){
  CountsSub <- as.data.table(dataset)
  #Spec6_Shrink[,c("Season", "Species", "Order", "Family", "Genus", "Count", "MigCode")] <- NULL
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  Range <- function(x) {
    max(x) - min(x)
  }
  
  variables <- VariablesForMatchingByYear
  variables <- variables[!variables %in% c("SWOccurrence", "Travel", "Slope", "GovMean")]
  SitesCovCast <- dcast(CountsSub, SiteCode ~ ., Mode, value.var="Anthrome", drop=FALSE)
  names(SitesCovCast) <- c("SiteCode", "Anthrome")
  SitesCovCast <- cbind(SitesCovCast, do.call(cbind, pblapply(variables, function(x){
    SitesCovCast <- dcast(CountsSub, SiteCode ~ ., mean, value.var=x, drop=FALSE)
    names(SitesCovCast) <- c("SiteCode", x)
    return(as.data.frame(SitesCovCast)[2])
  })))
  SitesCovStatic <- CountsSub[match(unique(CountsSub$SiteCode), CountsSub$SiteCode),][,c("SiteCode", "ISO3", "Country", "GeoRegion", "GeoSubRegion", "SWOccurrence", "LOTWID", "Travel", "Slope", "GovMean")]
  SitesCovCast <- merge(SitesCovCast, SitesCovStatic, by="SiteCode")
  SitesCovCast$AnthRound <- as.factor(round_any(SitesCovCast$Anthrome, 10))
  
  return(SitesCovCast)
}
CalculateBeforeSlopes <- function(BeforeData, ZeroThresh, parallelise, keepinsig){
  if(parallelise==TRUE){
    ncoresslopes <- ncores
  } else {
    ncoresslopes <- 1
  }
  
  AllCounts <- AllZeroes(BeforeData, ZeroThresh)
  NonZero <- subset(AllCounts, AllZero=="Counts")
  AllZero <- subset(AllCounts, AllZero=="AllZero")
  names(AllZero) <- c("SiteSpec", "BeforeSlope")
  if(nrow(NonZero)==0){
    return(AllCounts)
  }
  
  BeforeSlopes <- pbmclapply(unique(NonZero$SiteSpec), function(SS){
    BeforeDat <- subset(BeforeData, SiteSpec==SS)[,c("Count", "Year", "Hours", "Dataset")]
    if(length(unique(BeforeDat$Count))==1){if(unique(BeforeDat$Count)==0){
      return(as.data.frame(cbind(slope="AllZero", significant="AllZero", SiteSpec=SS)))}
    }
    if(unique(BeforeDat$Dataset=="CBC")){
      glmmy <- suppressWarnings(tryCatch(glm.nb(Count~Year + offset(log(BeforeDat$Hours)), link=log, data=BeforeDat), error=function(e){NA}))
    } else {
      glmmy <- suppressWarnings(tryCatch(glm.nb(Count~Year, link=log, data=BeforeDat), error=function(e){NA}))
    }
    slope <- tryCatch(coef(glmmy, silent=TRUE)[[2]], error=function(e){NA})  
    significant <- tryCatch(summary(glmmy)$coeff[-1,4], error=function(e){NA})
    return(as.data.frame(cbind(slope, significant, SiteSpec=SS)))
  }, mc.cores=ncoresslopes)
  
  BeforeSlopes <- as.data.frame(rbindlist(BeforeSlopes[unlist(sapply(BeforeSlopes, function(x) ifelse(ncol(x)==2, FALSE, TRUE)))]))
  BeforeSlopes$significant <- as.numeric(as.character(BeforeSlopes$significant))
  BeforeSlopes$slope <- as.numeric(as.character(BeforeSlopes$slope))
  NumYears <- dcast(as.data.table(BeforeData), SiteSpec~., length, value.var="Year")
  BeforeSlopes <- merge(BeforeSlopes, NumYears, by="SiteSpec")
  if(keepinsig==TRUE){
    BeforeSlopes$BeforeSlope <- ifelse(is.na(BeforeSlopes$significant), 0, #If the significance is NA make it zero
                                       ifelse(BeforeSlopes$slope=="AllZero", "AllZero", #If there's AllZero, make it all zero
                                              ifelse(BeforeSlopes$significant<0.05, ifelse(BeforeSlopes$slope>0,1,-1), 0)))
  } else {
    BeforeSlopes$BeforeSlope <- ifelse(is.na(BeforeSlopes$significant), 0, #If the significance is NA make it zero
                                       ifelse(BeforeSlopes$slope=="AllZero", "AllZero", #If there's AllZero, make it all zero
                                              ifelse(BeforeSlopes$.>6, ifelse(BeforeSlopes$slope>0,1,-1), #If the significance is not NA and there's more than 6 years, make it the slope
                                                     ifelse(BeforeSlopes$significant>0.05, 0, #If the significance is not NA, it's 6 or less years and the significance is >0.05 make it zero
                                                            ifelse(BeforeSlopes$slope>0, 1, -1))))) #Otherwise make it the slope    
  }
  
  BeforeSlopes <- rbind(BeforeSlopes[,c("SiteSpec", "BeforeSlope")], AllZero)
  if(nrow(BeforeSlopes)!=nrow(BeforeSlopes[complete.cases(BeforeSlopes),])){stop("there are NAs in before slopes")}
  return(BeforeSlopes)
}
AllZeroes <- function(CountData, ZeroThresh){
  if(min(CountData$Count)==0){
    ZeroCounts <- dcast(as.data.table(subset(CountData, Count==0)), SiteSpec~., length, value.var="Count")
    names(ZeroCounts) <- c("SiteSpec", "ZeroCounts")
    AllCounts <- dcast(as.data.table(CountData), SiteSpec ~., length, value.var="Count")
    names(AllCounts) <- c("SiteSpec", "AllCounts")
    
    AllCounts <- as.data.frame(merge(AllCounts, ZeroCounts, by="SiteSpec", all=T))
    if(nrow(AllCounts)!=nrow(AllCounts[complete.cases(AllCounts),])){
      AllCounts[is.na(AllCounts$ZeroCounts),]$ZeroCounts <- 0
    }
    AllCounts$AllZero <- ifelse(AllCounts$ZeroCounts/AllCounts$AllCounts<=ZeroThresh, "Counts", "AllZero")
    AllCounts <- AllCounts[,c("SiteSpec", "AllZero")]
    
  } else {
    AllCounts <- data.frame(SiteSpec=unique(CountData$SiteSpec), AllZero="Counts")
    AllCounts$SiteSpec <- as.character(AllCounts$SiteSpec)
  }
  return(AllCounts)
}
ParallelSlopes <- function(Comp1, Comp2, ParallelP, ZeroThresh, parallelise){
  if(parallelise==TRUE){
    ncoresslopes <- ncores
  } else {
    ncoresslopes <- 1
  }
  
  Counts1 <- subset(Comp1, BA==0)
  Counts2 <- subset(Comp2, BA==0)
  Counts1Zeroes <- AllZeroes(Counts1, ZeroThresh)
  Counts2Zeroes <- AllZeroes(Counts2, ZeroThresh)
  if(Counts1Zeroes$AllZero=="AllZero"){
    Counts2Zeroes <- subset(Counts2Zeroes, AllZero=="AllZero")
    return(Counts2Zeroes)
  }
  
  Counts2Zeroes <- subset(Counts2Zeroes, AllZero=="Counts")
  if(nrow(Counts2Zeroes)==0){
    return(Counts2Zeroes)
  }
  
  Parallel <- rbindlist(pbmclapply(unique(Counts2Zeroes$SiteSpec), function(SS){
    CountsSub <- subset(Counts2, SiteSpec==SS)
    Counts <- rbindlist(list(Counts1, CountsSub), fill=TRUE)[,c("Count", "Year", "CI", "BA", "Hours", "Dataset")]
    Counts[is.na(Counts$Hours),]$Hours <- 1
    CountsComp <- suppressWarnings(tryCatch(glm.nb(Count~ Year + CI + Year*CI + offset(log(Hours)), link=log, data=Counts), error=function(e){NA}))
    slope <- tryCatch(coef(CountsComp, silent=TRUE)[[4]], error=function(e){NA})  
    significant <- tryCatch(summary(CountsComp)$coeff[4,4], error=function(e){NA})
    return(data.frame("Sig"=significant, "Slope"=slope, "SiteSpec"=SS))
  }, mc.cores=ncoresslopes))
  Parallel <- subset(Parallel[complete.cases(Parallel),], Sig>ParallelP)
  return(Parallel)
}
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
  
  #Next, check for collinearity
  TheVIF <- as.data.frame(vif(CorrPred))
  while(nrow(subset(TheVIF, VIF=="Inf"))>0){
    TheVIF <- subset(TheVIF, VIF=="Inf")[1,]
    CorrPred <- CorrPred[,!(names(CorrPred) %in% TheVIF$Variables)]
    TheVIF <- as.data.frame(vif(CorrPred))
  }
  while(nrow(subset(TheVIF, VIF=="NaN"))>0){
    TheVIF <- subset(TheVIF, VIF=="NaN")[1,]
    CorrPred <- CorrPred[,!(names(CorrPred) %in% TheVIF$Variables)]
    TheVIF <- as.data.frame(vif(CorrPred))
  }
  while(max(TheVIF$VIF)>4){
    TheVIF <- subset(TheVIF, VIF==max(TheVIF$VIF))
    CorrPred <- CorrPred[,!(names(CorrPred) %in% TheVIF$Variables)]
    TheVIF <- as.data.frame(vif(CorrPred))
  }
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
  
  #Remove any variables that are mostly the same value
  SameValue <- apply(CorrPred, 2, function(x) length(x[x==Mode(x)])/length(x))
  SameValue <- SameValue[SameValue<0.8]
  VariablesForModel <- names(SameValue)
  CorrPred <- CorrPred[names(CorrPred) %in% VariablesForModel]
  
  dataset2 <- cbind(dataset[!colnames(dataset) %in% variables], CorrPred)
  return(list(dataset2, VariablesForModel))
}
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
BarPlotWidths <- function(WidthData, Log=TRUE, Cov, FactorLevels, xlabel, AngledxLabs=FALSE, ColourList, SuccessValues, MainResults, aspectrat){
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
  
  Plot <- ggplot(WidthData, aes(x=Position, y=PercentVal, fill=variable))+
    geom_bar(stat="identity", position="stack", width=WidthData$LogWidth)+ 
    scale_fill_manual(values = ColourList, name=("Category"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
    scale_x_continuous(expand = c(0, 0), breaks=WidthDataUnique$Position, labels=WidthDataUnique$Cov) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("Proportion")+
    xlab(xlabel)+
    theme(aspect.ratio=aspectrat, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, size=plotfontsize),
          legend.title=element_text(size=plotfontsize),
          axis.text.x = xtext,
          axis.ticks = element_line(colour="black"),
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top")
  return(list(Plot, WidthData))
}
Map <- function(Coordinates){
  mapWorld <- borders("world", colour="gray70", fill="gray70")
  Coordinates <- unique(Coordinates[,c("Latitude", "Longitude", "SiteCode", "CI"),])
  Coordinates$CI <- factor(Coordinates$CI)
  Coordinates$CI <- recode_factor(Coordinates$CI, '1' = "Protected", '0' ="Unprotected")
  Coordinates <- Coordinates[order(Coordinates$CI, decreasing=TRUE),]
  ggplot() + mapWorld +
    geom_point(aes(x=Coordinates$Longitude, y=Coordinates$Latitude,  colour=Coordinates$CI), size=0.6, shape=19, alpha=0.9)+
    scale_colour_manual(values=c("#0072B2", "#CC79A7"), guide=FALSE)+
    theme(aspect.ratio=0.53, panel.grid = element_blank(), strip.background = element_blank(),
          strip.text = element_text(hjust = 0), legend.justification = "top", axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_blank(), panel.border = element_rect(size = 1, fill = NA, colour="black"),
          plot.background=element_rect(fill = "transparent"), panel.background = element_rect(fill = "transparent"))
}
CleanBA <- function(BAData, ZeroThresh){
  BA <- BAData
  
  AfterZeroes <- AllZeroes(subset(BA, BA==1), ZeroThresh)
  names(AfterZeroes) <- c("SiteSpec", "AfterSlope")
  BA <- merge(BA, AfterZeroes, by="SiteSpec")
  BA$BeforeSlope2 <- ifelse(BA$BeforeSlope=="AllZero", "AllZero", "Counts")
  
  AllZero <- unique(BA[,c("SiteSpec", "BeforeSlope", "AfterSlope")])
  AllZero <- subset(AllZero, BeforeSlope=="AllZero" & AfterSlope=="AllZero")
  
  BA <- BA[!BA$SiteSpec %in% AllZero$SiteSpec,]
  BA[is.na(BA$Hours),]$Hours <- 1
  
  BA$ModelCat <- ifelse(BA$BeforeSlope2=="AllZero" | BA$AfterSlope=="AllZero", "Categorise", "Model")
  return(BA)
}
CleanBACI <- function(BACIData, ZeroThresh, RandomisationTest){
  BACI <- BACIData
  
  #Get slopes/all zeroes before and after
  BeforeSlopes <- unique(BACI[,c("SpecMatch", "BeforeSlope")])[complete.cases(unique(BACI[,c("SpecMatch", "BeforeSlope")])),]
  BACI$BeforeSlope <- NULL
  BACI <- merge(BACI, BeforeSlopes, all=T)
  
  AfterSlopes <- AllZeroes(subset(BACI, BA==1), ZeroThresh)
  names(AfterSlopes) <- c("SiteSpec", "AfterSlope")
  BACI <- merge(BACI, AfterSlopes, all=T)
  
  #Randomisation tests
  if(RandomisationTest=="Random"){
    BACI <- transform(BACI, Count = sample(Count) )
  }
  if(RandomisationTest=="NegBin"){
    Zeroes <- nrow(subset(BACI, Count==0))
    BACI$Count <- rnbinom(nrow(BACI), 50, 0.5)
    BACI[sample(nrow(BACI),Zeroes),]$Count <- 0
  }
  if(RandomisationTest=="NegBinZeroesConstant"){
    Zeroes <- nrow(subset(BACI, Count!=0))
    BACI[BACI$Count!=0,]$Count <- rnbinom(Zeroes, 50, 0.5)
  }
  
  #Remove cases with all zeroes across the board
  BACI$BeforeSlope2 <- ifelse(BACI$BeforeSlope=="AllZero", "AllZero", "Counts")
  AllZeroBACI <- melt(as.data.table(unique(BACI[,c("SpecMatch", "SiteSpec", "BeforeSlope2", "AfterSlope")])), id.vars=c("SpecMatch", "SiteSpec"))
  AllZeroBACI <- dcast(AllZeroBACI, SpecMatch~value, length, value.var="variable")
  AllZeroBACI$ModelCat <- ifelse(AllZeroBACI$AllZero==4, "Remove", ifelse(AllZeroBACI$AllZero>0, "Categorise", "Model"))
  if(nrow(AllZeroBACI)!=length(unique(BACI$SpecMatch))){stop("SpecMatchLost")}
  BACI <- merge(BACI, AllZeroBACI[,c("SpecMatch", "ModelCat")])
  BACI <- subset(BACI, ModelCat!="Remove")
  return(BACI)
} #Randomisation test can be "None, Random, NegBin or NegBinZeroesConstant
BACIPopulationCategorise <- function(BACI, MatchedFP){
  ### Mean Slope Models
  BACIModel <- subset(BACI, ModelCat=="Model")
  NullModel <- rbindlist(pbmclapply(unique(BACIModel$SpecMatch), function(x){
    print(x)
    BACIPop <- subset(BACIModel, SpecMatch==x)
    NullModel <- tryCatch(glm.nb(Count~Year + CI + CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL})
    if(is.null(NullModel)){
      return(data.frame(SpecMatch=x, Model="Failed"))
    }
    FullModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL})
    if(is.null(FullModel)){
      return(data.frame(SpecMatch=x, Model="Failed"))
    }
    if(summary(FullModel)[[11]][7,4]<0.05){
      return(data.frame(SpecMatch=x, Model="NotParallel"))
    }
    if(AIC(FullModel)<AIC(NullModel)){
      return(data.frame(SpecMatch=x, Model="Model"))
    } else{
      return(data.frame(SpecMatch=x, Model="Null"))
    }
  }, mc.cores=ncores))
  BACIModel2 <- subset(NullModel, Model=="Model")
  NullModel <- subset(NullModel, Model!="Model")
  NullModel <- merge(NullModel, unique(BACI[,c("SpecMatch", "BeforeSlope")]), all.x=T)
  NullModel[NullModel$Model=="Null"]$Model <- ifelse(NullModel[NullModel$Model=="Null"]$BeforeSlope=="-1", "NoChangeNeg", "NoChangePosStable")
  
  BACIByPopulation <- pbmclapply(unique(BACIModel2$SpecMatch), function(x){
    BACIPop <- subset(BACIModel, SpecMatch==x)
    BACIPop$Year <- BACIPop$Year-BACIPop$STATUS_YR
    
    SlopeModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL})
    if(is.null(SlopeModel)){
      return(NULL)
    }
    SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff)
    SlopeOutput$Coef <- row.names(SlopeOutput)
    SlopeOutput$Model <- "Slope"
    
    MeanModel <- tryCatch(glm.nb(Count~ BA + CI + BA*CI + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL})
    if(is.null(MeanModel)){
      return(NULL)
    }
    MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
    MeanOutput$Coef <- row.names(MeanOutput)
    MeanOutput$Model <- "Mean"
    
    ModelOutput <- rbind(SlopeOutput, MeanOutput)
    
    ModelOutput$SpecMatch <- x
    ModelOutput$Class <- unique(BACIPop$BeforeSlope)
    return(ModelOutput)
  }, mc.cores=ncores)
  
  if(length(BACIByPopulation[!sapply(BACIByPopulation, is.null)])/length(BACIByPopulation)!=1){stop("There are still NUll models coming through")}
  BACIByPopulation <- rbindlist(BACIByPopulation)
  names(BACIByPopulation) <- c("Estimate", "Error", "z", "P", "Coef", "Model", "SpecMatch", "Class")
  if(nrow(subset(BACIByPopulation, Coef=="Year:CI" & P<0.05))>0){stop("There are still non parallel slopes coming through")}
  BACIMeanSlope <- subset(BACIByPopulation, Model=="Mean" & Coef=="BA:CI" | Coef=="Year:BA:CI" & Model=="Slope")
  dir.create(file.path(MatchedFP, "ModelOutput/"), showWarnings = FALSE)
  write.csv(BACIMeanSlope, file=paste0(MatchedFP, "ModelOutput/BACIMeanSlope.csv"), row.names=FALSE)
  
  BACIMeanSlope$Category <- ifelse(BACIMeanSlope$P>0.05, "Insig", ifelse(BACIMeanSlope$Estimate>0, "Pos", "Neg"))
  BACIMeanSlope <- dcast(BACIMeanSlope, SpecMatch~Model, value.var="Category")
  BACIMeanSlope <- merge(BACIMeanSlope, unique(BACI[,c("SpecMatch", "BeforeSlope")]), all.x=T)
  BACIMeanSlope$Outcome <- ifelse(BACIMeanSlope$Mean =="Pos" | BACIMeanSlope$Slope =="Pos", "Success", ifelse(BACIMeanSlope$Mean=="Insig" & BACIMeanSlope$Slope=="Insig", ifelse(BACIMeanSlope$BeforeSlope=="-1", "NoChangeNeg", "NoChangePosStable"), "Failure"))
  
  BACIModelledOutcomes <- rbind(BACIMeanSlope[,c("SpecMatch", "Outcome")], NullModel[,c("SpecMatch", "Model")], use.names=FALSE)
  if(nrow(BACIModelledOutcomes)!=length(unique(BACIModel$SpecMatch))){stop("Some specmatch have been lost!")}
  if("Null" %in% unique(BACIModelledOutcomes$Outcome)){stop("there are NULLs")}
  BACIModelledOutcomes$OutcomeSummary <- ifelse(BACIModelledOutcomes$Outcome =="Failure", "Failure",
                                                ifelse(BACIModelledOutcomes$Outcome == "Success", "Success",
                                                       ifelse(BACIModelledOutcomes$Outcome =="NoChangeNeg" | BACIModelledOutcomes$Outcome =="NoChangePosStable", "Neutral", "Excluded")))
  
  ### All Zero Data
  BACICat <- subset(BACI, ModelCat=="Categorise")
  BACICat <- unique(BACICat[,c("SpecMatch", "CI", "BA", "AfterSlope", "BeforeSlope2")])
  BACICatBefore <- subset(BACICat, BA==0)[,c("SpecMatch", "CI", "BA", "BeforeSlope2")]
  names(BACICatBefore) <- c("SpecMatch", "CI", "BA","Slope")
  BACICatAfter <- subset(BACICat, BA==1)[,c("SpecMatch", "CI", "BA", "AfterSlope")]
  names(BACICatAfter) <- c("SpecMatch", "CI", "BA","Slope")
  BACICat <- rbind(BACICatBefore, BACICatAfter)
  
  BACICat <- dcast(as.data.table(BACICat), SpecMatch~BA+CI, value.var="Slope")
  names(BACICat) <- c("SpecMatch", "BC", "BI", "AC", "AI")
  BACICat$Outcome <- ifelse(BACICat$BC=="AllZero" & BACICat$BI=="AllZero" & BACICat$AC=="Counts" & BACICat$AI=="Counts", "ColonisedBoth",
                            ifelse(BACICat$BC=="AllZero" & BACICat$BI=="AllZero" & BACICat$AC=="AllZero" & BACICat$AI=="Counts", "ColonisedIn",
                                   ifelse(BACICat$BC=="AllZero" & BACICat$BI=="AllZero" & BACICat$AC=="Counts" & BACICat$AI=="AllZero", "ColonisedOut",
                                          ifelse(BACICat$BC=="Counts" & BACICat$BI=="Counts" & BACICat$AC=="AllZero" & BACICat$AI=="AllZero", "ExtinctBoth",
                                                 ifelse(BACICat$BC=="Counts" & BACICat$BI=="Counts" & BACICat$AC=="AllZero" & BACICat$AI=="Counts","ExtinctOut", "ExtinctIn")))))
  
  BACICat$OutcomeSummary <- ifelse(BACICat$Outcome=="ExtinctinBoth" | BACICat$Outcome=="ColonisedBoth", "Neutral", 
                                   ifelse(BACICat$Outcome=="ExtinctOut" | BACICat$Outcome=="ColonisedIn", "Success", "Failure"))
  
  ### Combine and summarise
  BACISummarise <- rbind(BACIModelledOutcomes[,c("SpecMatch", "Outcome", "OutcomeSummary")], BACICat[,c("SpecMatch", "Outcome", "OutcomeSummary")])
  if(nrow(BACISummarise)!=length(unique(BACISummarise$SpecMatch))){stop("There are duplicates")}
  if(nrow(BACISummarise)!=length(unique(BACI$SpecMatch))){stop("We've lost some specmatch")}
  
  BACISummarise$Species <- str_split_fixed(BACISummarise$SpecMatch, "[.]", 2)[,1]
  BACISummarise <- merge(BACISummarise, unique(subset(BACI, CI==1)[,c("SpecMatch", "SiteCode")]))
  
  write.csv(BACISummarise, file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"), row.names=FALSE)
  return("Done")
}
BAPopulationCategorise <- function(BA, MatchedFP){
  ### Mean Slope Models
  BAModel <- subset(BA, ModelCat=="Model")
  NullModel <- rbindlist(pbmclapply(unique(BAModel$SiteSpec), function(x){
    print(x)
    BAPop <- subset(BAModel, SiteSpec==x)
    NullModel <- tryCatch(glm.nb(Count~Year + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL})
    if(is.null(NullModel)){
      return(data.frame(SiteSpec=x, Model="Failed"))
    }
    FullModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL})
    if(is.null(FullModel)){
      return(data.frame(SiteSpec=x, Model="Failed"))
    }
    if(AIC(FullModel)<AIC(NullModel)){
      return(data.frame(SiteSpec=x, Model="Model"))
    } else{
      return(data.frame(SiteSpec=x, Model="Null"))
    }
  }, mc.cores=ncores))
  BAModel2 <- subset(NullModel, Model=="Model")
  NullModel <- subset(NullModel, Model!="Model")
  NullModel <- merge(NullModel, unique(BA[,c("SiteSpec", "BeforeSlope")]), all.x=T)
  NullModel[NullModel$Model=="Null"]$Model <- ifelse(NullModel[NullModel$Model=="Null"]$BeforeSlope=="-1", "NoChangeNeg", "NoChangePosStable")
  
  BAByPopulation <- pblapply(unique(BAModel2$SiteSpec), function(x){
    BAPop <- subset(BAModel, SiteSpec==x)
    BAPop$Year <- BAPop$Year - BAPop$STATUS_YR
    SlopeModel <- tryCatch(glm.nb(Count~Year + BA + BA*Year + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL})
    if(is.null(SlopeModel)){
      return(NULL)
    }
    SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff)
    SlopeOutput$Coef <- row.names(SlopeOutput)
    SlopeOutput$Model <- "Slope"
    
    MeanModel <- tryCatch(glm.nb(Count~ BA + offset(log(Hours)), link=log, data=BAPop), error=function(e){NULL})
    if(is.null(MeanModel)){
      return(NULL)
    }
    MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
    MeanOutput$Coef <- row.names(MeanOutput)
    MeanOutput$Model <- "Mean"
    
    ModelOutput <- rbind(SlopeOutput, MeanOutput)
    
    ModelOutput$SiteSpec <- x
    ModelOutput$Class <- unique(BAPop$BeforeSlope)
    return(ModelOutput)
  })
  
  if(length(BAByPopulation[!sapply(BAByPopulation, is.null)])/length(BAByPopulation)!=1){stop("There are still NUll models coming through")}
  BAByPopulation <- rbindlist(BAByPopulation)
  names(BAByPopulation) <- c("Estimate", "Error", "z", "P", "Coef", "Model", "SiteSpec", "Class")
  
  dir.create(file.path(MatchedFP, "ModelOutput/"), showWarnings = FALSE)
  write.csv(BAByPopulation, file=paste0(MatchedFP, "ModelOutput/BAMeanSlope.csv"), row.names=FALSE)
  
  BAMeanSlope <- subset(BAByPopulation, Model=="Mean" & Coef=="BA" | Coef=="Year:BA" & Model=="Slope")
  BAMeanSlope$Category <- ifelse(BAMeanSlope$P>0.05, "Insig", ifelse(BAMeanSlope$Estimate>0, "Pos", "Neg"))
  BAMeanSlope <- dcast(BAMeanSlope, SiteSpec~Model, value.var="Category")
  BAMeanSlope <- merge(BAMeanSlope, unique(BA[,c("SiteSpec", "BeforeSlope")]))
  BAMeanSlope$Outcome <- ifelse(BAMeanSlope$Mean =="Pos" | BAMeanSlope$Slope =="Pos", "Success", ifelse(BAMeanSlope$Mean=="Insig" & BAMeanSlope$Slope=="Insig", ifelse(BAMeanSlope$BeforeSlope=="-1", "NoChangeNeg", "NoChangePosStable"), "Failure"))
  
  BAModelledOutcomes <- rbind(BAMeanSlope[,c("SiteSpec", "Outcome")], NullModel[,c("SiteSpec", "Model")], use.names=FALSE)
  if(nrow(BAModelledOutcomes)!=length(unique(BAModel$SiteSpec))){stop("Some specmatch have been lost!")}
  
  BAModelledOutcomes$OutcomeSummary <- ifelse(BAModelledOutcomes$Outcome =="Failure", "Failure",
                                              ifelse(BAModelledOutcomes$Outcome == "Success", "Success",
                                                     ifelse(BAModelledOutcomes$Outcome=="Neutral" |BAModelledOutcomes$Outcome=="Null" | BAModelledOutcomes$Outcome =="NoChangeNeg" | BAModelledOutcomes$Outcome =="NoChangePosStable", "Neutral", "Excluded")))
  
  ### All zeroes
  BACat <- subset(BA, ModelCat=="Categorise")
  BACat <- unique(BACat[,c("SiteSpec", "CI", "BA", "AfterSlope", "BeforeSlope2")])
  BACatBefore <- subset(BACat, BA==0)[,c("SiteSpec", "CI", "BA", "BeforeSlope2")]
  names(BACatBefore) <- c("SiteSpec", "CI", "BA","Slope")
  BACatAfter <- subset(BACat, BA==1)[,c("SiteSpec", "CI", "BA", "AfterSlope")]
  names(BACatAfter) <- c("SiteSpec", "CI", "BA","Slope")
  BACat <- rbind(BACatBefore, BACatAfter)
  
  BACat <- dcast(as.data.table(BACat), SiteSpec~BA, value.var="Slope")
  names(BACat) <- c("SiteSpec", "Before", "After")
  BACat$Outcome <- ifelse(BACat$Before=="AllZero" & BACat$After=="Counts", "Colonised", "Extinct")
  
  BACat$OutcomeSummary <- ifelse(BACat$Outcome=="Colonised", "Success", "Failure")
  
  ### Combine and summarise
  BASummarise <- rbind(BAModelledOutcomes[,c("SiteSpec", "Outcome", "OutcomeSummary")], BACat[,c("SiteSpec", "Outcome", "OutcomeSummary")])
  if(nrow(BASummarise)!=length(unique(BASummarise$SiteSpec))){stop("There are duplicates")}
  if(nrow(BASummarise)!=length(unique(BA$SiteSpec))){stop("We've lost some specmatch")}
  
  BASummarise$Species <- str_split_fixed(BASummarise$SiteSpec, "[.]", 2)[,2]
  write.csv(BASummarise, file=paste0(MatchedFP, "ModelOutput/BACategories.csv"), row.names=FALSE)
  return("Done")
}
SuccessMap <- function(MapData, FileName){
  mapWorld <- borders("world", colour="gray70", fill="gray70")
  mapEurope <- borders(database="world", xlim = c(-11, 50), ylim = c(35, 75), colour="grey70", fill="grey70")
  
  World  <- ggplot() + mapWorld +
    geom_point(aes(x=MapData$Longitude, y=MapData$Latitude,  colour=MapData$SuccessPerc), size=2, shape=19, alpha=0.7)+
    scale_color_brewer(palette = "RdBu")+  
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
      plot.background=element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"))
  
  Europe  <- ggplot() + mapEurope +
    geom_point(aes(x=MapData$Longitude, y=MapData$Latitude,  colour=MapData$SuccessPerc), size=2, shape=19, alpha=0.7)+
    scale_color_brewer(palette = "RdBu", guide=FALSE)+ 
    coord_map(xlim = c(-11, 30), ylim = c(35, 60))+
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(size = 1, fill = NA, colour="black"),
      legend.justification = "top",
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.background=element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"))
  
  png(paste0(MatchedFP, "Results/", FileName, ".png"), 19, 9, units="in", res=300, bg="transparent") #Save as an image
  grid.newpage()
  v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
  v2<-viewport(width = 0.44, height = 0.39, x = -0.1, y = 0.02, just=c("left", "bottom")) #plot area for the inset map
  print(World,vp=v1) 
  print(Europe,vp=v2, mar=c(0,0,0,0))
  dev.off()
}

BlueCols <- colorRampPalette(c("#060437ff", "#3e6aabff", "#e0f3f8ff"))
RedCols <- colorRampPalette(c("#fffdbfff", "#d73027ff", "#680e18ff"))
x <- 5
plot(rep(1,x),col=BlueCols(x),pch=19,cex=3)
plot(rep(1,x),col=RedCols(x),pch=19,cex=3)


#### Define Scenario ####
Database <- "Waterbird" #LPI #Waterbird
CBCBuffer <- "CBCBuffer" #NoCBCbuffer #CBCBuffer
PABuffer <- "1km" #1km, 10km

DatabaseFP <- paste0(ResultsFP, Database,"_",CBCBuffer,"_PABuffer", PABuffer, "/")
dir.create(file.path(DatabaseFP), showWarnings = FALSE)

ZeroThresh <- 0.7
ZeroThreshFP <- 7 #6, 7, 8, 9, 10
TotalYearsBuffer <- 10 #5, 10, 20
MeasuredYearsBuffer <- 7 #3, 7, 15
ImputeFlag <- "ZeroesImputed" #NoImputedZeroes #ZeroesImputed
ParallelP <- 0.05
ParallelPFP <- "05"
InputFP <- paste0(DatabaseFP, "ZeroThresh", ZeroThreshFP, "_TimeSpan", TotalYearsBuffer, "_ParallelP", ParallelPFP,"_", ImputeFlag, "/")
dir.create(file.path(InputFP), showWarnings = FALSE)

StDiffThresh <- 0.25
StDiffThreshFP <- "025"
PropPairsThresh <- 0.8
PropPairsThreshFP <- "08"

MatchedFP <- paste0(InputFP, "StDiff", StDiffThreshFP, "_PropPairs", PropPairsThreshFP, "/")
dir.create(file.path(MatchedFP), showWarnings = FALSE)

#### Assess Matching ####
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
  dir.create(file.path(paste0(MatchedFP, "MatchSummaries/")), showWarnings = FALSE)
  dir.create(file.path(paste0(MatchedFP, "MatchData/")), showWarnings = FALSE)
  DeleteOldFiles <- c(list.files(paste0(MatchedFP, "MatchSummaries/"), full.names=TRUE), list.files(paste0(MatchedFP, "MatchData/"), full.names=TRUE))
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
      write.csv(Summary, paste0(MatchedFP, "MatchSummaries/Summary_", i, ".csv"), row.names=FALSE)
      write.csv(MatchingSub, paste0(MatchedFP, "MatchData/MatchData_", i, ".csv"), row.names=FALSE)
      return(i)
    } else {
      return(NULL)
    }
  }, mc.cores=ncores)
}
CompMatch <- CompareMatching(MatchingCovariates, StDiffThresh) #Run through and iteratively remove matches until we have the matched set for each species

SummariesFinal <- rbindlist(lapply(list.files(path=paste0(MatchedFP, "MatchSummaries/"), full.names=TRUE), fread)) #Read in that summary material
nrow(subset(SummariesFinal, PropPairs>PropPairsThresh))/nrow(SummariesFinal)
SummariesFinal <- subset(SummariesFinal, PropPairs>PropPairsThresh)
MatchDataFinal <- rbindlist(lapply(list.files(path=paste0(MatchedFP, "MatchData/"), full.names=TRUE), fread))
MatchingFinalCleaned <- MatchingFinal[MatchingFinal$SpecMatch %in% MatchDataFinal$SpecMatch,]
MatchingFinalCleaned <- MatchingFinalCleaned[MatchingFinalCleaned$Species %in% unique(SummariesFinal$Species),]

TheNumbers(subset(MatchingFinalCleaned, CI==1))
TheNumbers(subset(MatchingFinalCleaned, CI==0))

save(MatchingFinalCleaned, file=paste0(MatchedFP, "MatchingFinalCleaned.RData"))

#Extract species stats for comparison
load(file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))
load(file=paste0(InputFP, "UnprotectedCountsReadyForMatching.RData"))

TaxaData <- function(Dataset){
  AllSitesStats <- as.data.table(unique(Dataset[,c("Species","Genus", "Family", "Order", "SiteCode", "CI")]))
  AllSitesStatsFamSite <- unique(AllSitesStats[,c("Order", "Family", "SiteCode", "CI")])
  FamiliesSite <- dcast(AllSitesStatsFamSite, Order + Family~CI, length, value.var="SiteCode")
  AllSitesStatsFamSpec <- unique(AllSitesStats[,c("Order", "Family","Genus", "Species", "CI")])
  FamiliesSpec <- dcast(AllSitesStatsFamSpec, Order + Family~., length, value.var="Species")
  AllSitesStatsFamGen <- unique(AllSitesStats[,c("Order", "Family", "Genus", "CI")])
  FamiliesGen <- dcast(AllSitesStatsFamGen, Order + Family~., length, value.var="Genus")
  FamilyData <- merge(FamiliesSite, FamiliesGen, all=TRUE)
  names(FamilyData) <- c("Order", "Family", "UnprotectedSites", "ProtectedSites", "Genera")
  FamilyData <- merge(FamilyData, FamiliesSpec, all=TRUE)
  names(FamilyData) <- c("Order", "Family", "UnprotectedSites", "ProtectedSites", "Genera", "Species")
  NOrder <- length(unique(FamilyData$Order))
  NFamily <- length(unique(FamilyData$Family))
  NGenus <- length(unique(AllSitesStatsFamSpec$Genus))
  NSpecies <- length(unique(AllSitesStatsFamSpec$Species))
  TaxaList <- list(FamilyData, NOrder,NFamily,NGenus, NSpecies)
  names(TaxaList) <- c("FamilyData", "NumOrders", "NumFamilies", "NumGenera", "NumSpecies")
  return(TaxaList)
}

AllData <- rbind(unique(ProtectedCountData[,c("SiteCode", "Species", "Genus", "Family", "Order", "CI")]), unique(UnprotectedCountData[,c("SiteCode", "Species", "Genus", "Family", "Order", "CI")]))
AllTaxa <- TaxaData(AllData)
MatchTaxa <- TaxaData(MatchingFinalCleaned)

AllTaxa
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
                                                       paste0(TheNumbers(ProtectedCountData)[[1]], " (", TheNumbers(subset(MatchingFinalCleaned, CI==1))[[1]], ")"), 
                                                       paste0(TheNumbers(UnprotectedCountData)[[1]], " (", TheNumbers(subset(MatchingFinalCleaned, CI==0))[[1]], ")"))))
names(AllOrdAllFam) <- names(BothTaxa)
BothTaxa <- rbind(BothTaxa, AllOrdAllFam)
write.csv(BothTaxa, paste0(MatchedFP, "TaxaMatchStats.txt"), row.names=FALSE)

#### Clean and Model Data ####
load(file=paste0(MatchedFP, "MatchingFinalCleaned.RData"))
BACI <- CleanBACI(MatchingFinalCleaned, ZeroThresh, "None")
load(file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))
BA <- CleanBA(ProtectedCountData, ZeroThresh)

BACIPopulationCategorise(BACI, MatchedFP)
BAPopulationCategorise(BA, MatchedFP)


#### Stacked Bar Plots ####
#BACI
BACISummarise <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
BACISummariseForPlot <- melt(dcast(as.data.table(subset(BACISummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SpecMatch"), id.vars="Species")

cols <- colorRampPalette(c("lightslategrey", "aliceblue"))
ColourListPlot <- c(BlueCols(5)[1:3],cols(3)[[2]], "aliceblue", "lightgoldenrod1", "goldenrod1", RedCols(5)[3:5])
FactorLevels <- c("ColonisedIn", "ExtinctOut", "Success", "ColonisedBoth", "NoChangePosStable", "NoChangeNeg", "ExtinctBoth","Failure", "ColonisedOut", "ExtinctIn")
BACIResultsPlot <- BarPlotWidths(BACISummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Success", "ColonisedIn", "ExtinctOut"), MainResults=TRUE, aspectrat=0.5)
BACIResultsPlot[[1]]

dir.create(file.path(paste0(MatchedFP, "Results/")), showWarnings=FALSE)
ggsave(paste0(MatchedFP, "Results/BACICategoriesBarPlot.pdf"), BACIResultsPlot[[1]], device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf)

#BA
BASummarise <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv"))
BASummarise <- melt(dcast(as.data.table(subset(BASummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SiteSpec"), id.vars="Species")
FactorLevels <- c("Colonised", "Success", "NoChangePosStable", "NoChangeNeg", "Failure", "Extinct")
ColourListPlot <- c(BlueCols(5)[1:2],"azure3", "burlywood2", RedCols(5)[4:5])

BAResultsPlot <- BarPlotWidths(BASummarise, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Colonised", "Success"), MainResults=TRUE, aspectrat=0.5)
BAResultsPlot[[1]]
dir.create(file.path(paste0(MatchedFP, "Results/")), showWarnings = FALSE)
ggsave(paste0(MatchedFP, "Results/BACategoriesBarPlot.pdf"), BAResultsPlot[[1]], device="pdf", width = 400, height = 150, units = "mm") #Save as a pdf)

#### BA CI BACI Stacked Bar and Comparison ####
BACISummarise <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
BACISummarise$SiteSpec <- paste0(BACISummarise$SiteCode, ".", BACISummarise$Species)
BASummarise <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv"))
BASummarise <- BASummarise[BASummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, BASummarise$SiteSpec),]
BACISummarise <- BACISummarise[BACISummarise$SiteSpec %in% intersect(BACISummarise$SiteSpec, BASummarise$SiteSpec),]

#Get CI Data
CI <- subset(BACI[BACI$SpecMatch %in% BACISummarise$SpecMatch,], BA==1)
CICast <- dcast(as.data.table(unique(CI[,c("SpecMatch", "CI", "SiteCode", "AfterSlope")])), SpecMatch~AfterSlope, length, value.var="SiteCode")
CI <- CI[CI$SpecMatch %in% subset(CICast, AllZero==0)$SpecMatch,]

CIByPopulation <- pbmclapply(unique(CI$SpecMatch), function(x){
  BACIPop <- subset(CI, SpecMatch==x)
  SlopeModel <- tryCatch(glm.nb(Count~Year + CI*Year + offset(log(Hours)), link=log, data=BACIPop), error=function(e){NULL})
  if(is.null(SlopeModel)){
    return(NULL)
  }
  ModelOutput <- as.data.frame(summary(SlopeModel)$coeff)
  ModelOutput$Coef <- row.names(ModelOutput)
  
  ModelOutput$SpecMatch <- x
  return(ModelOutput)
}, mc.cores=ncores)

if(length(CIByPopulation[!sapply(CIByPopulation, is.null)])/length(CIByPopulation)!=1){stop("There are still NUll models coming through")}
CIByPopulation <- rbindlist(CIByPopulation)
names(CIByPopulation) <- c("Estimate", "Error", "z", "P", "Coef", "SpecMatch")
CISlope <- subset(CIByPopulation, Coef=="Year:CI")
CISlope <- merge(CISlope, unique(subset(CI, CI==1)[,c("SpecMatch", "SiteSpec")]))

CIAfterSlope <- CalculateBeforeSlopes(CI, ZeroThresh,parallelise=TRUE, keepinsig=FALSE)
names(CIAfterSlope) <- c("SiteSpec", "AfterSlope")

CISlope <- merge(CISlope, CIAfterSlope, by="SiteSpec")
CISlope$Outcome <- ifelse(CISlope$P<0.05, ifelse(CISlope$Estimate<0, "Failure", "Success"), ifelse(CISlope$AfterSlope=="-1", "NoChangeNeg", "NoChangePosStable"))
CISlope$OutcomeSummary <- ifelse(CISlope$Outcome=="Failure"|CISlope$Outcome=="Success", CISlope$Outcome, "Neutral")

#Plot BACI, BA and CI
BACISummarise <- BACISummarise[BACISummarise$SiteSpec %in% intersect(intersect(BACISummarise$SiteSpec,BASummarise$SiteSpec),CISlope$SiteSpec),]
BASummarise <- BASummarise[BASummarise$SiteSpec %in% intersect(intersect(BACISummarise$SiteSpec,BASummarise$SiteSpec),CISlope$SiteSpec),]
CISlope <- CISlope[CISlope$SiteSpec %in% intersect(intersect(BACISummarise$SiteSpec,BASummarise$SiteSpec),CISlope$SiteSpec),]

#Make a BACI plot
BACISummariseForPlot <- melt(dcast(as.data.table(subset(BACISummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SpecMatch"), id.vars="Species")
ColourListPlot <- c(BlueCols(5)[2:3],"azure3", "burlywood2", RedCols(5)[4])
FactorLevels <- c("ColonisedIn", "ExtinctOut", "Success", "ColonisedBoth", "NoChangePosStable", "NoChangeNeg", "ExtinctBoth","Failure", "ColonisedOut", "ExtinctIn")
BACIResultsPlot <- BarPlotWidths(BACISummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Success", "ColonisedIn", "ExtinctOut"), MainResults=TRUE, aspectrat=0.5)
WidthData <- BACIResultsPlot[[2]]

#Make BA Plot
BASummariseForPlot <- melt(dcast(as.data.table(subset(BASummarise, OutcomeSummary!="Excluded")), Species~Outcome, length, value.var="SiteSpec"), id.vars="Species")
BASummariseForPlot <- merge(BASummariseForPlot, unique(WidthData[,c("Species", "LogWidth", "Width", "Position")]), by="Species")
BASummarisePercent <- dcast(BASummariseForPlot, Species~., sum, value.var="value")
names(BASummarisePercent) <- c("Species", "Total")
BASummariseForPlot <- merge(BASummariseForPlot, BASummarisePercent)
BASummariseForPlot$PercentVal <- BASummariseForPlot$value/BASummariseForPlot$Total
BASummariseForPlot$variable <- factor(BASummariseForPlot$variable, levels=c("Colonised", "Success", "NoChangePosStable", "NoChangeNeg", "Failure", "Extinct"))
ColourListPlot <- c(BlueCols(5)[1:2],"azure3", "burlywood2", RedCols(5)[4:5])
BASummariseForPlotUnique <- unique(BASummariseForPlot[,c("Species", "LogWidth", "Position")])

BAStacked <- ggplot(BASummariseForPlot, aes(x=Position, y=PercentVal, fill=variable))+
  geom_bar(stat="identity", position="stack", width=BASummariseForPlot$LogWidth) +
  scale_fill_manual(values = ColourListPlot, name=("Category"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
  scale_x_continuous(expand = c(0, 0), breaks=BASummariseForPlotUnique$Position, labels=BASummariseForPlotUnique$Species) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Proportion")+
  xlab("Species")+
  theme(aspect.ratio=0.5, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_blank(),
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
BAStacked

#Make a CI plot
CISlope$Species <- str_split_fixed(CISlope$SiteSpec, "[.]", 2)[,2]
CISummariseForPlot <- melt(dcast(as.data.table(CISlope), Species~Outcome, length, value.var="SiteSpec"), id.vars="Species")
CISummariseForPlot <- merge(CISummariseForPlot, unique(WidthData[,c("Species", "LogWidth", "Width", "Position")]), by="Species")
CISummarisePercent <- dcast(CISummariseForPlot, Species~., sum, value.var="value")
names(CISummarisePercent) <- c("Species", "Total")
CISummariseForPlot <- merge(CISummariseForPlot, CISummarisePercent)
CISummariseForPlot$PercentVal <- CISummariseForPlot$value/CISummariseForPlot$Total
CISummariseForPlot$variable <- factor(CISummariseForPlot$variable, levels=c("Success", "NoChangePosStable", "NoChangeNeg", "Failure"))
ColourListPlot <- c(BlueCols(5)[2],"azure3", "burlywood2", RedCols(5)[4])
CISummariseForPlotUnique <- unique(CISummariseForPlot[,c("Species", "LogWidth", "Position")])

CIStacked <- ggplot(CISummariseForPlot, aes(x=Position, y=PercentVal, fill=variable))+
  geom_bar(stat="identity", position="stack", width=CISummariseForPlot$LogWidth) +
  scale_fill_manual(values = ColourListPlot, name=("Category"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
  scale_x_continuous(expand = c(0, 0), breaks=CISummariseForPlotUnique$Position, labels=CISummariseForPlotUnique$Species) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Proportion")+
  xlab("Species")+
  theme(aspect.ratio=0.5, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_blank(),
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
CIStacked

#Brind them together
library(cowplot)
AllPlot <- plot_grid(BACIResultsPlot[[1]], BAStacked,CIStacked, ncol=1, align="hv")
ggsave(paste0(MatchedFP, "Results/BACI_BA_CI_CompPlot.pdf"), AllPlot, device="pdf", width = 300, height = 450, units = "mm") #Save as a pdf)

###Compare the three
BACIBACI <- Reduce(merge, list(CISlope[,c("SiteSpec", "OutcomeSummary")], BACISummarise[,c("SiteSpec", "OutcomeSummary")], BASummarise[,c("SiteSpec", "OutcomeSummary")]))
names(BACIBACI) <- c("SiteSpec", "CI", "BACI", "BA")

BACIBACI$ChangeBA <- ifelse(BACIBACI$BA=="Success" & BACIBACI$BACI=="Failure", "Success to Failure",
                            ifelse(BACIBACI$BA=="Success" & BACIBACI$BACI=="Neutral", "Success to Neutral",
                                   ifelse(BACIBACI$BA=="Failure" & BACIBACI$BACI=="Success", "Failure to Success", 
                                          ifelse(BACIBACI$BA=="Failure" & BACIBACI$BACI=="Neutral", "Failure to Neutral",
                                                 ifelse(BACIBACI$BA=="Neutral" & BACIBACI$BACI=="Success", "Neutral to Success", 
                                                        ifelse(BACIBACI$BA=="Neutral" & BACIBACI$BACI=="Failure","Neutral to Failure",
                                                               ifelse(BACIBACI$BA=="Neutral" & BACIBACI$BACI=="Neutral","Neutral to Neutral",
                                                                      ifelse(BACIBACI$BA=="Success" & BACIBACI$BACI=="Success","Success to Success", "Failure to Failure"))))))))

BACIBACI$ChangeCI <- ifelse(BACIBACI$CI=="Success" & BACIBACI$BACI=="Failure", "Success to Failure",
                            ifelse(BACIBACI$CI=="Success" & BACIBACI$BACI=="Neutral", "Success to Neutral",
                                   ifelse(BACIBACI$CI=="Failure" & BACIBACI$BACI=="Success", "Failure to Success", 
                                          ifelse(BACIBACI$CI=="Failure" & BACIBACI$BACI=="Neutral", "Failure to Neutral",
                                                 ifelse(BACIBACI$CI=="Neutral" & BACIBACI$BACI=="Success", "Neutral to Success", 
                                                        ifelse(BACIBACI$CI=="Neutral" & BACIBACI$BACI=="Failure","Neutral to Failure",
                                                               ifelse(BACIBACI$CI=="Neutral" & BACIBACI$BACI=="Neutral","Neutral to Neutral",
                                                                      ifelse(BACIBACI$CI=="Success" & BACIBACI$BACI=="Success","Success to Success", "Failure to Failure"))))))))

BACIBACI$Species <- str_split_fixed(BACIBACI$SiteSpec, "[.]",2)[,2]

#CompareBA
CompareBASummariseForPlot <- melt(dcast(as.data.table(BACIBACI[,c("SiteSpec", "Species", "ChangeBA")]), Species~ChangeBA, length, value.var="SiteSpec"), id.vars="Species")
FactorLevels <- c("Failure to Success", "Neutral to Success", "Failure to Neutral", "Success to Success", "Neutral to Neutral", "Failure to Failure", "Success to Neutral", "Neutral to Failure", "Success to Failure")
ColourListPlot <- c(BlueCols(5)[c(1:3)],"antiquewhite", "antiquewhite3", "antiquewhite4", RedCols(5)[c(3:5)])

BACompPlot <- BarPlotWidths(CompareBASummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Success", "ColonisedIn", "ExtinctOut"), MainResults=TRUE, aspectrat=0.5)
BACompPlot[[1]]
ggsave(paste0(MatchedFP, "Results/BACI_BA_DifferencePlot.pdf"), BACompPlot[[1]], device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf)

#CompareCI
CompareCISummariseForPlot <- melt(dcast(as.data.table(BACIBACI[,c("SiteSpec", "Species", "ChangeCI")]), Species~ChangeCI, length, value.var="SiteSpec"), id.vars="Species")
CICompPlot <- BarPlotWidths(CompareCISummariseForPlot, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Success", "ColonisedIn", "ExtinctOut"), MainResults=TRUE, aspectrat=0.5)
CICompPlot[[1]]
ggsave(paste0(MatchedFP, "Results/BACI_CI_DifferencePlot.pdf"), CICompPlot[[1]], device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf)

#Straight comparison plots
#BACIBACI <- BACIBACISave
#BACIBACI <- subset(BACIBACI, BACI!="Neutral")
BACIBACISummarise1 <- dcast(BACIBACI, ChangeBA~., length, value.var="SiteSpec")
BACIBACISummarise2 <- dcast(BACIBACI, ChangeCI~., length, value.var="SiteSpec")
BACIBACISummarise <- merge(BACIBACISummarise1, BACIBACISummarise2, by.x="ChangeBA", by.y="ChangeCI")
names(BACIBACISummarise) <- c("Change", "BA", "CI")
BACIBACISummarise <- as.data.frame(melt(BACIBACISummarise, id.vars="Change"))
BACIBACISummarise$Change <- factor(BACIBACISummarise$Change, levels=rev(FactorLevels))

BA_CI_BACI_Diff <- ggplot(aes(x=Change, y=value), data=BACIBACISummarise)+
  geom_bar(stat="identity", aes(fill=Change))+
  facet_wrap(~variable)+
  coord_flip()+
  scale_fill_manual(values=rev(ColourListPlot), guide=FALSE)+
  xlab("Change from BA/CI to BACI")+
  ylab("Number of Populations")+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.ticks = element_line(colour="black"),
        legend.justification = "top")
BA_CI_BACI_Diff
ggsave(paste0(MatchedFP, "Results/BACI_BA_CI_DifferenceSum.pdf"), BA_CI_BACI_Diff, device="pdf", width = 200, height = 100, units = "mm") #Save as a pdf)

#### Plot MeanSlope ####
BACIMeanSlope <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACIMeanSlope.csv"))
BACIMeanSlope <- read.csv(file=paste0(MatchedFP, "ModelOutput/BAMeanSlope.csv"))
names(BACIMeanSlope)[names(BACIMeanSlope)=="SiteSpec"] <- "SpecMatch"
BACIMeanSlope <- as.data.table(subset(BACIMeanSlope, Coef=="Year:BA" & Model=="Slope" | Coef=="BA" & Model=="Mean"))

#comp level mean
MeanVsLevel <- subset(BACIMeanSlope, Coef=="BA")
right <- dcast(MeanVsLevel, SiteSpec~Model, value.var="Estimate")
plot(right$Mean, right$Slope)
cor(right$Mean, right$Slope)

SlopevsLevel <- subset(BACIMeanSlope, Model=="Slope")
right <- dcast(SlopevsLevel, SiteSpec~Coef, value.var="Estimate")
plot(right$BA, right$`Year:BA`)
cor(right$BA, right$`Year:BA`)

SlopevsMean <- as.data.table(subset(BACIMeanSlope, Coef=="Year:BA" & Model=="Slope" | Coef=="BA" & Model=="Mean"))
right <- dcast(SlopevsMean, SiteSpec~Model, value.var="Estimate")
plot(right$Mean, right$Slope)
cor(right$Mean, right$Slope)


BACIMeanSlopeEst <- dcast(BACIMeanSlope, SpecMatch + Class~Model, value.var="Estimate")
BACIMeanSlopeError <- dcast(BACIMeanSlope, SpecMatch + Class~Model, value.var="Error")
BACIMeanSlopeP <- dcast(BACIMeanSlope, SpecMatch + Class~Model, value.var="P")
BACIMeanSlope <- list(BACIMeanSlopeEst, BACIMeanSlopeError, BACIMeanSlopeP) %>% reduce(left_join, by = c("SpecMatch", "Class"))
names(BACIMeanSlope) <- c("SpecMatch", "BeforeSlope", "MeanEst", "SlopeEst", "MeanError", "SlopeError", "MeanP", "SlopeP")
BACIMeanSlope[,c(3:ncol(BACIMeanSlope))] <- apply(BACIMeanSlope[,c(3:ncol(BACIMeanSlope))], 2, as.numeric)

#Get Order of Species, add significance field
#BACIOrder <- unique(BACI[,c("Species", "Order")])
#BACIMeanSlope <- merge(BACIMeanSlope, BACIOrder, by="Species")
BACIMeanSlope$Significance <- ifelse(BACIMeanSlope$MeanP<0.05, ifelse(BACIMeanSlope$SlopeP<0.05, "Both Significant", "Mean Significant"), ifelse(BACIMeanSlope$SlopeP<0.05, "Slope Significant", "Neither Significant"))
BACIMeanSlope <- subset(BACIMeanSlope, Significance!="Neither Significant")
BACIMeanSlope <- BACIMeanSlope[complete.cases(BACIMeanSlope),]
BACIMeanSlope$Significance <- factor(BACIMeanSlope$Significance, levels=c("Both Significant", "Slope Significant", "Mean Significant", "Neither Significant"))
BACIMeanSlope2 <- subset(BACIMeanSlope, Significance!="Neither Significant")
BACIMeanSlope2 <- BACIMeanSlope2[order(BACIMeanSlope2$Significance, decreasing=TRUE),]
BACIMeanSlope2 <- merge(BACIMeanSlope2, AfterSlopeModelBACI, by="SpecMatch")

#Add in after slope
BACIMeanSlope2 <- subset(BACIMeanSlope2, Significance=="Both Significant")
#MeanSlopePlot <- 
ggplot(data=BACIMeanSlope2, aes(x=MeanEst, y=SlopeEst, colour=After))+ #, fill=Result
  geom_point(alpha=0.5)+ #aes(colour=Category,shape=Category)
  #facet_grid(~BeforeSlope)+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  #scale_colour_manual(values = c("black", "grey40", "grey80"), name=("Significance"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
  scale_colour_distiller(palette = "Spectral", limits=c(-1,1), direction=1)+
  xlab("Mean")+ #Proportional change in mean (before to after)
  ylab("Slope")+ #"Proportional change in slope (before to after)"
  #scale_shape_manual(values=c(0:3))+
  xlim(-5,5)+
  ylim(-3,3)+
  #geom_errorbar(aes(ymin=SlopeEst-SlopeError, ymax=SlopeEst+SlopeError), width=.2)+
  #geom_errorbarh(aes(xmin=MeanEst-MeanError, xmax=MeanEst+MeanError), height=.07)+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_text(size=plotfontsize, colour="black", angle = 90), #element_text(size=plotfontsize, colour="black", angle = 90)
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")

ggsave(paste0(MatchedFP, "Results/BACIMeanSlope.pdf"), MeanSlopePlot)

ok <- BACIMeanSlope2
ok <- subset(BACIMeanSlope2, Significance=="Both Significant")
nrow(subset(ok, MeanEst>0 & SlopeEst>0))
nrow(subset(ok, MeanEst>0 & SlopeEst<0))
nrow(subset(ok, MeanEst<0 & SlopeEst<0))
nrow(subset(ok, MeanEst<0 & SlopeEst>0))

#### Map Figure ####
#BA by category
BASummarise <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv"))
BASummarise$SiteCode <- str_split_fixed(BASummarise$SiteSpec, "[.]", 2)[,1]
BASummariseSuccess <- dcast(as.data.table(BASummarise), SiteCode~OutcomeSummary, length, value.var="SiteSpec")
BASummariseSuccess$SuccessPerc <- BASummariseSuccess$Success/rowSums(BASummariseSuccess[,c(2:5)])*100
BASummariseSuccess <- merge(BASummariseSuccess[,c("SiteCode", "SuccessPerc")], unique(BA[,c("SiteCode", "Latitude", "Longitude")]), by="SiteCode")
BASummariseSuccess <- BASummariseSuccess[order(BASummariseSuccess$SuccessPerc),]
BASummariseSuccess$SuccessPerc <- as.factor(round_any(BASummariseSuccess$SuccessPerc, 10, f=ceiling))
SuccessMap(BASummariseSuccess, "BACategoriesMap")

#BACI by category
BACISummarise <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
BACISummariseSuccess <- dcast(as.data.table(BACISummarise), SiteCode~OutcomeSummary, length, value.var="SpecMatch")
BACISummariseSuccess$SuccessPerc <- BACISummariseSuccess$Success/rowSums(BACISummariseSuccess[,c(2:5)])*100
BACISummariseSuccess <- merge(BACISummariseSuccess[,c("SiteCode", "SuccessPerc")], unique(BACI[,c("SiteCode", "Latitude", "Longitude")]), by="SiteCode")
BACISummariseSuccess <- BACISummariseSuccess[order(BACISummariseSuccess$SuccessPerc),]
BACISummariseSuccess$SuccessPerc <- as.factor(round_any(BACISummariseSuccess$SuccessPerc, 10, f=ceiling))
SuccessMap(BACISummariseSuccess, "BACICategoriesMap")

#### Absolute after slopes ####
BAMeanSlope <- read.csv(file=paste0(MatchedFP, "ModelOutput/BAMeanSlope.csv"))
BASlope <- subset(BAMeanSlope, Model=="Slope" & Coef=="Year" | Coef=="Year:BA")
#BASlope$Estimate <- ifelse(BASlope$P<0.05, BASlope$Estimate, 0)
AfterSlope <- dcast(as.data.table(BASlope), SiteSpec ~ Coef, value.var="Estimate")
AfterP <- dcast(as.data.table(BASlope), SiteSpec ~ Coef, value.var="P")
names(AfterP) <- c("SiteSpec", "YearP", "YearBAP")
AfterSlope <- merge(AfterSlope, AfterP)
AfterSlope$Year <- ifelse(AfterSlope$YearBAP>0.05, ifelse(AfterSlope$YearP<0.05, AfterSlope$Year, 0), AfterSlope$Year)
AfterSlope$`Year:BA` <- ifelse(AfterSlope$YearBAP<0.05, AfterSlope$`Year:BA`, 0)
AfterSlope$After <- AfterSlope$Year + AfterSlope$`Year:BA`

#If the interaction term is significant, we want the Year term regardless
#If the interaction term is insignificant, we want the Year term if it's significant

#Ok recalculate after slopes
BAAfter <- subset(BA, BA==1)
AfterSlopeModel <- rbindlist(pbmclapply(unique(BAAfter$SiteSpec), function(x){
  Pop <- subset(BAAfter, SiteSpec==x)
  PopModel <- tryCatch(glm.nb(Count~Year + offset(log(Hours)), link=log, data=Pop), error=function(e){NULL})
  if(is.null(PopModel)){
    return(data.frame("Estimate"=0, "StError"=0, "z"=0, "P"=0, "Coef"="Year", "SiteSpec"=x))
  }
  PopModelOutput <- as.data.frame(summary(PopModel)$coeff)
  PopModelOutput$Coef <- row.names(PopModelOutput)
  names(PopModelOutput) <- c("Estimate", "StError", "z", "P", "Coef")
  PopModelOutput <- subset(PopModelOutput, Coef=="Year")
  PopModelOutput$SiteSpec <- x
  return(PopModelOutput)
}, mc.cores=ncores))
AfterSlopeModel$Estimate2 <- ifelse(AfterSlopeModel$P<0.05, AfterSlopeModel$Estimate, 0)
AfterSlopeModel <- AfterSlopeModel[,c("Estimate", "SiteSpec")]
names(AfterSlopeModel) <- c("After2", "SiteSpec")

try <- merge(AfterSlopeModel, AfterSlope[,c("SiteSpec", "After")], by="SiteSpec")

### Now BACI
BACIAfter <- subset(BACI, BA==1 & CI==1)
AfterSlopeModelBACI <- rbindlist(pbmclapply(unique(BACIAfter$SpecMatch), function(x){
  Pop <- subset(BACIAfter, SpecMatch==x)
  PopModel <- tryCatch(glm.nb(Count~Year + offset(log(Hours)), link=log, data=Pop), error=function(e){NULL})
  if(is.null(PopModel)){
    return(data.frame("Estimate"=0, "StError"=0, "z"=0, "P"=0, "Coef"="Year", "SpecMatch"=x))
  }
  PopModelOutput <- as.data.frame(summary(PopModel)$coeff)
  PopModelOutput$Coef <- row.names(PopModelOutput)
  names(PopModelOutput) <- c("Estimate", "StError", "z", "P", "Coef")
  PopModelOutput <- subset(PopModelOutput, Coef=="Year")
  PopModelOutput$SpecMatch <- x
  return(PopModelOutput)
}, mc.cores=ncores))
AfterSlopeModelBACI <- AfterSlopeModelBACI[,c("Estimate", "P", "SpecMatch")]
names(AfterSlopeModelBACI) <- c("After", "AfterP", "SpecMatch")
#### Prop by Site or Species ####
BACISummarise <- as.data.table(read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv")))
BACISummarise$SiteCode <- str_split_fixed(BACISummarise$SiteSpec, "[.]",2)[,1]

SiteProp <- dcast(BACISummarise, SiteCode~OutcomeSummary, length, value.var="Species")
SiteProp <- as.data.table(cbind(as.character(SiteProp$SiteCode), prop.table(as.matrix(SiteProp[,c(2:5)]), margin=1)))
SiteMaxProp <- dcast(melt(SiteProp, id.vars="V1"), V1~., max, value.var="value")
names(SiteMaxProp) <- c("ID", "MaxProp")
SiteMaxProp$MaxProp <- as.numeric(SiteMaxProp$MaxProp)
mean(SiteMaxProp$MaxProp)
SiteMaxProp$Var <- "SiteCode"


SpeciesProp <- dcast(BACISummarise, Species~OutcomeSummary, length, value.var="SiteCode")
SpeciesProp <- as.data.table(cbind(as.character(SpeciesProp$Species), prop.table(as.matrix(SpeciesProp[,c(2:5)]), margin=1)))
SpeciesMaxProp <- dcast(melt(SpeciesProp, id.vars="V1"), V1~., max, value.var="value")
names(SpeciesMaxProp) <- c("ID", "MaxProp")
SpeciesMaxProp$MaxProp <- as.numeric(SpeciesMaxProp$MaxProp)
mean(SpeciesMaxProp$MaxProp)
SpeciesMaxProp$Var <- "Species"

MaxProp <- rbind(SiteMaxProp, SpeciesMaxProp)

ggplot(data=MaxProp, aes(x=Var, y=MaxProp))+
  geom_boxplot()+
  geom_beeswarm(size=0.05)

#### BACI Model by Species ####
BACIModel <- subset(BACI, ModelCat=="Model")
load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))

BACIBySpeciesWithCovariates <- rbindlist(lapply(c("1","0","-1"), function(BS){
  dir.create(file.path(paste0(MatchedFP, "SpeciesModels/")), showWarnings = FALSE)
  dir.create(file.path(paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/")), showWarnings = FALSE)
  Species <- rbindlist(pbmclapply(unique(BACIModel$Species), function(x){
    if(file.exists(paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/", x, "_", BS, ".csv"))){
      return(NULL)
    }
    write.csv(NULL, paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/", x, "_", BS, ".csv"))
    
    ModelOutput <- data.frame(matrix(ncol = 6, nrow = 1))
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model")
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- BS
    
    BACISpec <- subset(BACIModel, Species==x & BeforeSlope == BS)
    if(nrow(BACISpec)==0){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("NoRows")
      write.csv(ModelOutput, file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    if(length(unique(BACISpec$SpecMatch))<4){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("Lessthan4specmatch")
      write.csv(ModelOutput, file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    BACISpecVIF <- CollVIF(BACISpec, VariablesForMatchingByYear, scale=TRUE)
    BACISpec <- BACISpecVIF[[1]]
    BACISpec$YearScale <- BACISpec$Year-BACISpec$STATUS_YR
    SlopeFormula <- as.formula(paste("Count", paste(c("YearScale", "CI", "BA", "CI*BA", "YearScale*CI", "YearScale*BA", "YearScale*CI*BA", 
                                                      BACISpecVIF[[2]], "(1|SiteCode)", "offset(log(BACISpec$Hours))"), 
                                                    collapse = " + "), sep = " ~ "))
    SlopeModel <- tryCatch(glmmTMB(SlopeFormula, data=BACISpec, family=nbinom1(link="log")), error=function(e){NULL})
    
    MeanFormula <- as.formula(paste("Count", paste(c("CI", "BA","CI*BA", 
                                                     BACISpecVIF[[2]], "(1|SiteCode)", "offset(log(BACISpec$Hours))"), 
                                                   collapse = " + "), sep = " ~ "))
    MeanModel <- tryCatch(glmmTMB(MeanFormula, data=BACISpec, family=nbinom1(link="log")), error=function(e){NULL})
    
    if(is.null(SlopeModel) & is.null(MeanModel)){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("NullModels")
      write.csv(ModelOutput, file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    
    if(!is.null(SlopeModel)){
      SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff[[1]])
      SlopeOutput$Coef <- row.names(SlopeOutput)
      SlopeOutput$Model <- "Slope"
      if(!is.null(MeanModel)){
        MeanOutput <- as.data.frame(summary(MeanModel)$coeff[[1]])
        MeanOutput$Coef <- row.names(MeanOutput)
        MeanOutput$Model <- "Mean"
        ModelOutput <- rbind(SlopeOutput, MeanOutput)
      } else {
        ModelOutput <- SlopeOutput
      }
    } else {
      MeanOutput <- as.data.frame(summary(MeanModel)$coeff[[1]])
      MeanOutput$Coef <- row.names(MeanOutput)
      MeanOutput$Model <- "Mean"
      ModelOutput <- MeanOutput
    }
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- unique(BACISpec$BeforeSlope)
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model", "Species", "BeforeSlope")
    write.csv(ModelOutput, file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
    return(ModelOutput)
  }, mc.cores=ncores))
  return(Species)
}))
save(BACIBySpeciesWithCovariates, file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates.RData"))

BACIBySpeciesWithCovariates <- rbindlist(lapply(list.files(paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates/"), full.names=TRUE), fread))

BACIMeanSlope <- subset(BACIBySpeciesWithCovariates, Model=="Mean" & Coef=="CI:BA" | Coef=="YearScale:CI:BA" & Model=="Slope")
BACIMeanSlopeEst <- dcast(BACIMeanSlope, Species + BeforeSlope~Model, value.var="Estimate")
BACIMeanSlopeError <- dcast(BACIMeanSlope, Species + BeforeSlope~Model, value.var="StError")
BACIMeanSlopeP <- dcast(BACIMeanSlope, Species + BeforeSlope~Model, value.var="P")
BACIMeanSlope <- list(BACIMeanSlopeEst, BACIMeanSlopeError, BACIMeanSlopeP) %>% reduce(left_join, by = c("Species", "BeforeSlope"))
names(BACIMeanSlope) <- c("Species", "BeforeSlope", "MeanEst", "SlopeEst", "MeanError", "SlopeError", "MeanP", "SlopeP")
BACIMeanSlope[,c(3:ncol(BACIMeanSlope))] <- apply(BACIMeanSlope[,c(3:ncol(BACIMeanSlope))], 2, as.numeric)

#Add significance field
BACIMeanSlope$Significance <- ifelse(BACIMeanSlope$MeanP<0.05, ifelse(BACIMeanSlope$SlopeP<0.05, "Both Significant", "Mean Significant"), ifelse(BACIMeanSlope$SlopeP<0.05, "Slope Significant", "Neither Significant"))
BACIMeanSlope <- subset(BACIMeanSlope, Significance!="Neither Significant")
BACIMeanSlope <- BACIMeanSlope[complete.cases(BACIMeanSlope),]
BACIMeanSlope2 <- subset(BACIMeanSlope, Significance == "Both Significant")
MeanSlopePlot <- ggplot(data=BACIMeanSlope, aes(x=MeanEst, y=SlopeEst, colour=Significance))+ #, fill=Result
  geom_point(alpha=0.5)+ #aes(colour=Category,shape=Category)
  #facet_grid(~BeforeSlope)+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_brewer(palette="Set1")+
  xlab("Mean")+ #Proportional change in mean (before to after)
  ylab("Slope")+ #"Proportional change in slope (before to after)"
  #scale_shape_manual(values=c(0:3))+
  xlim(-2,2)+
  ylim(-0.4,0.4)+
  geom_errorbar(aes(ymin=SlopeEst-SlopeError, ymax=SlopeEst+SlopeError), width=.12)+
  geom_errorbarh(aes(xmin=MeanEst-MeanError, xmax=MeanEst+MeanError), height=.02)+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_text(size=plotfontsize, colour="black", angle = 90), #element_text(size=plotfontsize, colour="black", angle = 90)
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
MeanSlopePlot
ggsave(paste0(MatchedFP, "Results/BACIMeanSlope.pdf"), MeanSlopePlot)

SpeciesCovs <- unique(SpeciesSiteCovs[,c("Species", "MigStatus", "Family", "GenLength", "EOO", "RedList")])
BACIMeanSlopeModel <- merge(BACIMeanSlope, SpeciesCovs, by=c("Species"))

SlopeModel <- glmmTMB(SlopeEst~ GenLength + MigStatus + RedList + (1|Family), data=BACIMeanSlopeModel)
MeanModel <- glmmTMB(MeanEst~ GenLength + MigStatus + RedList + (1|Family), data=BACIMeanSlopeModel)

BACISpeciesCat <- BACIMeanSlopeModel
BACISpeciesCat$MeanEst <- ifelse(BACISpeciesCat$MeanP<0.05, BACISpeciesCat$MeanEst, "Insig")
BACISpeciesCat$SlopeEst <- ifelse(BACISpeciesCat$SlopeP<0.05, BACISpeciesCat$SlopeEst, "Insig")
BACISpeciesCat$Category <- ifelse(BACISpeciesCat$MeanEst=="Insig" & BACISpeciesCat$SlopeEst=="Insig", "Neutral",
                                  ifelse(BACISpeciesCat$MeanEst=="Insig" & BACISpeciesCat$SlopeEst>0, "Positive",
                                         ifelse(BACISpeciesCat$MeanEst=="Insig" & BACISpeciesCat$SlopeEst<0,"Negative", 
                                                ifelse(BACISpeciesCat$SlopeEst=="Insig" & BACISpeciesCat$MeanEst>0,"Positive",
                                                       ifelse(BACISpeciesCat$SlopeEst=="Insig" & BACISpeciesCat$MeanEst<0,"Negative",
                                                              ifelse(BACISpeciesCat$SlopeEst>0 | BACISpeciesCat$MeanEst>0,"Positive", "Negative"))))))
BACISpeciesCat$Category <- factor(BACISpeciesCat$Category, levels=c("Negative", "Neutral", "Positive"))
CategoryModel <- clmm(Category~MigStatus + RedList + GenLength + (1|Family), data=BACISpeciesCat)
#### BACI Model by Site ####
#BACI modelled by site (remove this, not helpful)
BACIModel <- subset(BACI, ModelCat=="Model")
BACIBySitewithCovariates

BACIBySiteWithCovariates <- rbindlist(lapply(c("1","0","-1"), function(BS){
  dir.create(file.path(paste0(MatchedFP, "SiteModels/")), showWarnings = FALSE)
  dir.create(file.path(paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/")), showWarnings = FALSE)
  ProtSites <- unique(subset(BACIModel, CI==1)$SiteCode)
  Sites <- rbindlist(pbmclapply(ProtSites, function(x){
    if(file.exists(paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/", x, "_", BS, ".csv"))){
      return(NULL)
    }
    write.csv(NULL, paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/", x, "_", BS, ".csv"))
    
    ModelOutput <- data.frame(matrix(ncol = 6, nrow = 1))
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model")
    ModelOutput$SiteCode <- x
    ModelOutput$BeforeSlope <- BS
    SpecMatchSite <- unique(subset(BACIModel, SiteCode==x & CI==1)$SpecMatch)
    BACISite <- subset(BACIModel[BACIModel$SpecMatch %in% SpecMatchSite,], BeforeSlope==BS)
    
    if(nrow(BACISite)==0){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("NoRows")
      write.csv(ModelOutput, file=paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    if(length(unique(BACISite$SpecMatch))<4){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("Lessthan4specmatch")
      write.csv(ModelOutput, file=paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    BACISiteVIF <- CollVIF(BACISite, VariablesForMatchingByYear, scale=TRUE)
    BACISite <- BACISiteVIF[[1]]
    BACISite$YearScale <- BACISite$Year-BACISite$STATUS_YR
    SlopeFormula <- as.formula(paste("Count", paste(c("YearScale", "CI", "BA", "CI*BA", "YearScale*CI", "YearScale*BA", "YearScale*CI*BA", 
                                                      BACISiteVIF[[2]], "(1|Species)", "offset(log(BACISite$Hours))"), 
                                                    collapse = " + "), sep = " ~ "))
    SlopeModel <- tryCatch(glmmTMB(SlopeFormula, data=BACISite, family=nbinom1(link="log")), error=function(e){NULL})
    
    MeanFormula <- as.formula(paste("Count", paste(c("CI", "BA","CI*BA", 
                                                     BACISiteVIF[[2]], "(1|Species)", "offset(log(BACISite$Hours))"), 
                                                   collapse = " + "), sep = " ~ "))
    MeanModel <- tryCatch(glmmTMB(MeanFormula, data=BACISite, family=nbinom1(link="log")), error=function(e){NULL})
    
    if(is.null(SlopeModel) & is.null(MeanModel)){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("NullModels")
      write.csv(ModelOutput, file=paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    
    if(!is.null(SlopeModel)){
      SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff[[1]])
      SlopeOutput$Coef <- row.names(SlopeOutput)
      SlopeOutput$Model <- "Slope"
      if(!is.null(MeanModel)){
        MeanOutput <- as.data.frame(summary(MeanModel)$coeff[[1]])
        MeanOutput$Coef <- row.names(MeanOutput)
        MeanOutput$Model <- "Mean"
        ModelOutput <- rbind(SlopeOutput, MeanOutput)
      } else {
        ModelOutput <- SlopeOutput
      }
    } else {
      MeanOutput <- as.data.frame(summary(MeanModel)$coeff[[1]])
      MeanOutput$Coef <- row.names(MeanOutput)
      MeanOutput$Model <- "Mean"
      ModelOutput <- MeanOutput
    }
    ModelOutput$SiteCode <- x
    ModelOutput$BeforeSlope <- unique(BACISite$BeforeSlope)
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model", "SiteCode", "BeforeSlope")
    write.csv(ModelOutput, file=paste0(MatchedFP, "SiteModels/BACIBySiteWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
    return(ModelOutput)
  }, mc.cores=ncores))
  return(Sites)
}))
save(BACIBySpeciesWithCovariates, file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates.RData"))


BACIMeanSlope <- subset(BACIBySiteWithCovariates, Model=="Mean" & Coef=="CI:BA" | Coef=="YearScale:CI:BA" & Model=="Slope")
BACIMeanSlopeEst <- dcast(BACIMeanSlope, SiteCode + BeforeSlope~Model, value.var="Estimate")
BACIMeanSlopeError <- dcast(BACIMeanSlope, SiteCode + BeforeSlope~Model, value.var="StError")
BACIMeanSlopeP <- dcast(BACIMeanSlope, SiteCode + BeforeSlope~Model, value.var="P")
BACIMeanSlope <- list(BACIMeanSlopeEst, BACIMeanSlopeError, BACIMeanSlopeP) %>% reduce(left_join, by = c("SiteCode", "BeforeSlope"))
names(BACIMeanSlope) <- c("SiteCode", "BeforeSlope", "MeanEst", "SlopeEst", "MeanError", "SlopeError", "MeanP", "SlopeP")
BACIMeanSlope[,c(3:ncol(BACIMeanSlope))] <- apply(BACIMeanSlope[,c(3:ncol(BACIMeanSlope))], 2, as.numeric)

#Add significance field
BACIMeanSlope$Significance <- ifelse(BACIMeanSlope$MeanP<0.05, ifelse(BACIMeanSlope$SlopeP<0.05, "Both Significant", "Mean Significant"), ifelse(BACIMeanSlope$SlopeP<0.05, "Slope Significant", "Neither Significant"))
BACIMeanSlope <- subset(BACIMeanSlope, Significance!="Neither Significant")
BACIMeanSlope <- BACIMeanSlope[complete.cases(BACIMeanSlope),]
BACIMeanSlope2 <- subset(BACIMeanSlope, Significance == "Both Significant")
MeanSlopePlot <- ggplot(data=BACIMeanSlope, aes(x=MeanEst, y=SlopeEst, colour=Significance))+ #, fill=Result
  geom_point(alpha=0.5)+ #aes(colour=Category,shape=Category)
  #facet_grid(~BeforeSlope)+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_brewer(palette="Set1")+
  xlab("Mean")+ #Proportional change in mean (before to after)
  ylab("Slope")+ #"Proportional change in slope (before to after)"
  #scale_shape_manual(values=c(0:3))+
  xlim(-2.5,2.5)+
  ylim(-0.6,0.6)+
  geom_errorbar(aes(ymin=SlopeEst-SlopeError, ymax=SlopeEst+SlopeError), width=.12)+
  geom_errorbarh(aes(xmin=MeanEst-MeanError, xmax=MeanEst+MeanError), height=.02)+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_text(size=plotfontsize, colour="black", angle = 90), #element_text(size=plotfontsize, colour="black", angle = 90)
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
MeanSlopePlot
ggsave(paste0(MatchedFP, "Results/BACIMeanSlopeBySite.pdf"), MeanSlopePlot)

SiteCovs <- unique(SpeciesSiteCovs[,c("SiteCode", "GovMean", "STATUS_YR", "ISO3", "PAAreaLog", "AnthMode", "Ramsar")])
BACIMeanSlopeModel <- merge(BACIMeanSlope, SiteCovs, by=c("SiteCode"))
BACIMeanSlopeModel$STATUS_YR <- BACIMeanSlopeModel$STATUS_YR-min(BACIMeanSlopeModel$STATUS_YR)

SlopeModel <- glmmTMB(SlopeEst~ GovMean + STATUS_YR + PAAreaLog + AnthMode + Ramsar + Ramsar*GovMean + PAAreaLog*GovMean + (1|ISO3), data=BACIMeanSlopeModel)
MeanModel <- glmmTMB(MeanEst~ GovMean + STATUS_YR + PAAreaLog + AnthMode + Ramsar + Ramsar*GovMean + PAAreaLog*GovMean + (1|ISO3), data=BACIMeanSlopeModel)

BACISpeciesCat <- BACIMeanSlopeModel
BACISpeciesCat$MeanEst <- ifelse(BACISpeciesCat$MeanP<0.05, BACISpeciesCat$MeanEst, "Insig")
BACISpeciesCat$SlopeEst <- ifelse(BACISpeciesCat$SlopeP<0.05, BACISpeciesCat$SlopeEst, "Insig")
BACISpeciesCat$Category <- ifelse(BACISpeciesCat$MeanEst=="Insig" & BACISpeciesCat$SlopeEst=="Insig", "Neutral",
                                  ifelse(BACISpeciesCat$MeanEst=="Insig" & BACISpeciesCat$SlopeEst>0, "Positive",
                                         ifelse(BACISpeciesCat$MeanEst=="Insig" & BACISpeciesCat$SlopeEst<0,"Negative", 
                                                ifelse(BACISpeciesCat$SlopeEst=="Insig" & BACISpeciesCat$MeanEst>0,"Positive",
                                                       ifelse(BACISpeciesCat$SlopeEst=="Insig" & BACISpeciesCat$MeanEst<0,"Negative",
                                                              ifelse(BACISpeciesCat$SlopeEst>0 | BACISpeciesCat$MeanEst>0,"Positive", "Negative"))))))
BACISpeciesCat$Category <- factor(BACISpeciesCat$Category, levels=c("Negative", "Neutral", "Positive"))
CategoryModel <- clmm(Category~MigStatus + RedList + GenLength + (1|Family), data=BACISpeciesCat)





#### BA Model by Species ####
BAModel <- subset(BA, ModelCat=="Model")
load(file=paste0(DatabaseFP, "VariablesForMatchingByYear.RData"))

BABySpeciesWithCovariates <- rbindlist(lapply(c(1,0,-1), function(BS){
  dir.create(file.path(paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/")), showWarnings = FALSE)
  Species <- rbindlist(pbmclapply(unique(BAModel$Species)[1:5], function(x){
    if(file.exists(paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/", x, "_", BS, ".csv"))){
      return(NULL)
    }
    write.csv(NULL, paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
    
    ModelOutput <- data.frame(matrix(ncol = 6, nrow = 1))
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model")
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- BS
    
    BASpec <- subset(BAModel, Species==x & BeforeSlope == BS)
    if(nrow(BASpec)==0){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("NoRows")
      write.csv(ModelOutput, paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    if(length(unique(BASpec$SiteSpec))<4){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("Lessthan4specmatch")
      write.csv(ModelOutput, paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    BASpecVIF <- CollVIF(BASpec, VariablesForMatchingByYear, scale=TRUE)
    BASpec <- BASpecVIF[[1]]
    BASpec$YearScale <- BASpec$Year-BASpec$STATUS_YR
    SlopeFormula <- as.formula(paste("Count", paste(c("YearScale", "BA", "YearScale*BA",
                                                      BASpecVIF[[2]], "(1|SiteCode)", "offset(log(BASpec$Hours))"), 
                                                    collapse = " + "), sep = " ~ "))
    SlopeModel <- tryCatch(glmer.nb(SlopeFormula, data=BASpec), error=function(e){NULL})
    MeanFormula <- as.formula(paste("Count", paste(c("BA",
                                                     BASpecVIF[[2]], "(1|SiteCode)", "offset(log(BASpec$Hours))"), 
                                                   collapse = " + "), sep = " ~ "))
    MeanModel <- tryCatch(glmer.nb(MeanFormula, data=BASpec), error=function(e){NULL})
    
    if(is.null(SlopeModel) & is.null(MeanModel)){
      ModelOutput[, c("Estimate", "StError", "Z", "P", "Coef", "Model")] <- c("NullModels")
      write.csv(ModelOutput, paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    
    if(!is.null(SlopeModel)){
      SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff)
      SlopeOutput$Coef <- row.names(SlopeOutput)
      SlopeOutput$Model <- "Slope"
      if(!is.null(MeanModel)){
        MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
        MeanOutput$Coef <- row.names(MeanOutput)
        MeanOutput$Model <- "Mean"
        ModelOutput <- rbind(SlopeOutput, MeanOutput)
      } else {
        ModelOutput <- SlopeOutput
      }
    } else {
      MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
      MeanOutput$Coef <- row.names(MeanOutput)
      MeanOutput$Model <- "Mean"
      ModelOutput <- MeanOutput
    }
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- unique(BASpec$BeforeSlope)
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model", "Species", "BeforeSlope")
    
    write.csv(ModelOutput, paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates/", x, "_", BS, ".csv"), row.names=FALSE)
    return(ModelOutput)
  }, mc.cores=ncores))
  return(Species)
}))
save(BABySpeciesWithCovariates, file=paste0(MatchedFP, "SpeciesModels/BABySpeciesWithCovariates.RData"))

#load(file=paste0(MatchedFP, "SpeciesModels/BACIBySpeciesWithCovariates.RData")) #UNHASH ONCE PROPER RUNS DONE

BAMeanSlope <- subset(BACIBySpeciesWithCovariates, Model=="Mean" & Coef=="BA" | Coef=="YearScale:BA" & Model=="Slope")
BAMeanSlopeEst <- dcast(BAMeanSlope, Species + BeforeSlope~Model, value.var="Estimate")
BAMeanSlopeError <- dcast(BAMeanSlope, Species + BeforeSlope~Model, value.var="StError")
BAMeanSlopeP <- dcast(BAMeanSlope, Species + BeforeSlope~Model, value.var="P")
BAMeanSlope <- list(BAMeanSlopeEst, BAMeanSlopeError, BAMeanSlopeP) %>% reduce(left_join, by = c("Species", "BeforeSlope"))
names(BAMeanSlope) <- c("Species", "BeforeSlope", "MeanEst", "SlopeEst", "MeanError", "SlopeError", "MeanP", "SlopeP")
BAMeanSlope[,c(3:ncol(BAMeanSlope))] <- apply(BAMeanSlope[,c(3:ncol(BAMeanSlope))], 2, as.numeric)

#Get Order of Species, add significance field
BAOrder <- unique(BA[,c("Species", "Order")])
BAMeanSlope <- merge(BAMeanSlope, BAOrder, by="Species")
BAMeanSlope$Significance <- ifelse(BAMeanSlope$MeanP<0.05, ifelse(BAMeanSlope$SlopeP<0.05, "Both Significant", "Mean Significant"), ifelse(BAMeanSlope$SlopeP<0.05, "Slope Significant", "Neither Significant"))
BAMeanSlope <- subset(BAMeanSlope, Significance!="Neither Significant")
BAMeanSlope <- BAMeanSlope[complete.cases(BAMeanSlope),]

MeanSlopePlot <- ggplot(data=BAMeanSlope, aes(x=MeanEst, y=SlopeEst, colour=Significance))+ #, fill=Result
  geom_point(alpha=0.5)+ #aes(colour=Category,shape=Category)
  #facet_grid(~BeforeSlope)+
  geom_vline(xintercept = 0, linetype="dashed")+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_color_brewer(palette="Set1")+
  xlab("Mean")+ #Proportional change in mean (before to after)
  ylab("Slope")+ #"Proportional change in slope (before to after)"
  #scale_shape_manual(values=c(0:3))+
  xlim(-2.5,2.5)+
  ylim(-0.75,0.75)+
  geom_errorbar(aes(ymin=SlopeEst-SlopeError, ymax=SlopeEst+SlopeError), width=.2)+
  geom_errorbarh(aes(xmin=MeanEst-MeanError, xmax=MeanEst+MeanError), height=.085)+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_text(size=plotfontsize, colour="black", angle = 90), #element_text(size=plotfontsize, colour="black", angle = 90)
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")

ggsave(paste0(MatchedFP, "Results/BAMeanSlope.pdf"), MeanSlopePlot)


#### Get PA Area ####
ProtectedSites <- fread(paste0(ResultsFP, "ProtectedSites_Uncleaned.csv"))
ProtectedSites <- ProtectedSites[ProtectedSites$SiteCode %in% BA$SiteCode]
ProtectedSites <- merge(ProtectedSites, unique(BA[,c("STATUS_YR", "SiteCode")]), by="SiteCode")
names(ProtectedSites)[names(ProtectedSites)=="STATUS_YR.y"] <- "UsedStatYr"
names(ProtectedSites)[names(ProtectedSites)=="STATUS_YR.x"] <- "STATUS_YR"

#Subset to cases where any other PAs were designated in the survey period
if(nrow(subset(ProtectedSites, STATUS_YR<UsedStatYr))>0){stop("There are cases of earlier status years - note this error will remain till you rerun matching")}
ProtectedSites <- subset(ProtectedSites, STATUS_YR>=UsedStatYr) #Remove me once matching is rerun!

ProtectedSites <- subset(ProtectedSites, STATUS_YR<(UsedStatYr+TotalYearsBuffer))

OriginalArea <- unique(subset(ProtectedSites, STATUS_YR==UsedStatYr)[,c("SiteCode", "GIS_AREA")])
OriginalArea <- dcast(OriginalArea, SiteCode~., max, value.var="GIS_AREA")
names(OriginalArea) <- c("SiteCode", "OriginalArea")
ProtectedSites <- merge(ProtectedSites, OriginalArea)
ProtectedSites$LargerThanOriginal <- ifelse(ProtectedSites$GIS_AREA>ProtectedSites$OriginalArea, 1, 0)

#Remove later designated PAs that are smaller than the original
ProtectedSites <- subset(ProtectedSites, STATUS_YR==UsedStatYr | GIS_AREA!=OriginalArea & LargerThanOriginal==1)

ProtectedSitesMean <- dcast(ProtectedSites, SiteCode~., mean, value.var="GIS_AREA")
names(ProtectedSitesMean) <- c("SiteCode", "MeanArea")
ProtectedSitesMean <- merge(ProtectedSitesMean, unique(ProtectedSites[,c("SiteCode", "OriginalArea")]))
ProtectedSitesMean$Diff <- abs(ProtectedSitesMean$MeanArea - ProtectedSitesMean$OriginalArea)/ProtectedSitesMean$OriginalArea
ProtectedSitesMean <- subset(ProtectedSitesMean, Diff<0.1)

PAArea <- ProtectedSitesMean[,c("SiteCode", "MeanArea")]
names(PAArea) <- c("SiteCode", "PAArea")

#### Get Covariates HAN DOUBLE CHECK THAT THE PA SYSTEM IS WORKING ####
#Family, Migstatus, Governance
load(file=paste0(InputFP, "ProtectedCountsReadyForMatching.RData"))
SpeciesSiteCovs <- unique(ProtectedCountData[,c("SiteCode", "Species", "MigStatus", "Family", "GovMean", "STATUS_YR", "ISO3")])

#Redlist
IUCN_REDLIST_KEY <- "412cdf89c6a8241d53e22f14d47c409c9a4baa87783ca3f3fc85dbc029228e2b"
RedListData <- rbindlist(pbmclapply(unique(SpeciesSiteCovs$Species), function(x) rl_search(name = x, key = IUCN_REDLIST_KEY)$result, mc.cores=ncores))
RedListData <- RedListData[,c("scientific_name", "category", "eoo_km2")]
names(RedListData) <- c("Species", "RedList", "EOO")
if(length(unique(RedListData$Species))!=length(unique(SpeciesSiteCovs$Species))){stop("SpeciesLost")}
SpeciesSiteCovs <- merge(SpeciesSiteCovs, RedListData)

#Gen Length
GenLength <- as.data.frame(fread(paste0(DataFP, "GenerationLength/genlengthbirdconsbiol2020table4.csv")))
GenLength <- GenLength[GenLength$Species %in% unique(SpeciesSiteCovs$Species),]
GenLength <- GenLength[,c("Species", "GenLength")]
if(length(unique(GenLength$Species))!=length(unique(SpeciesSiteCovs$Species))){stop("SpeciesLost")}
SpeciesSiteCovs <- merge(SpeciesSiteCovs, GenLength)

#PA Area
SpeciesSiteCovs <- merge(SpeciesSiteCovs, PAArea, all=T)

#Anthrome
BA$AnthRound <- as.numeric(as.character(BA$AnthRound))
AnthMode <- dcast(as.data.table(BA), SiteCode~., Mode, value.var="AnthRound")
names(AnthMode) <- c("SiteCode", "AnthMode")
if(length(unique(AnthMode$SiteCode))!=length(unique(SpeciesSiteCovs$SiteCode))){stop("SpeciesLost")}
SpeciesSiteCovs <- merge(SpeciesSiteCovs, AnthMode)

#Ramsar sites
ProtectedSites <- fread(paste0(ResultsFP, "ProtectedSites_Uncleaned.csv"))
Ramsar <- ProtectedSites[grep("\\Ramsar\\b", ProtectedSites$DESIG)]

SpeciesSiteCovs$Ramsar <- "Other"
SpeciesSiteCovs[SpeciesSiteCovs$SiteCode %in% Ramsar$SiteCode,]$Ramsar <- "Ramsar"

#Clean up
SpeciesSiteCovs$AnthMode <- as.character(SpeciesSiteCovs$AnthMode)
SpeciesSiteCovs$PAAreaLog <- log(SpeciesSiteCovs$PAArea)
SpeciesSiteCovs$RedList <- factor(SpeciesSiteCovs$RedList, levels=c("LC", "NT", "VU", "EN", "CR"))
SpeciesSiteCovs$MigStatus <- factor(SpeciesSiteCovs$MigStatus)
SpeciesSiteCovs$MigStatus <- recode_factor(SpeciesSiteCovs$MigStatus, "Resident"="Resident", "Breeding"="Migrant", "Non-breeding"="Migrant", "Passage"="Migrant")
SpeciesSiteCovs$Family <- factor(SpeciesSiteCovs$Family)
SpeciesSiteCovs$AnthMode <- factor(SpeciesSiteCovs$AnthMode, levels=c('10', '20', '30', '40', '50', '60'))
SpeciesSiteCovs$AnthMode <- recode_factor(SpeciesSiteCovs$AnthMode, '10'='Urban', '20'='Village', '30'='Croplands', '40'='Rangeland', '50'='Semi-natural', '60'='Wild')
SpeciesSiteCovs$Ramsar <- factor(SpeciesSiteCovs$Ramsar)
SpeciesSiteCovs$ISO3 <- factor(SpeciesSiteCovs$ISO3)

#### Covariate Models by Category ####
# Model Categorised Data
dir.create(file.path(paste0(MatchedFP, "Drivers/")), showWarnings = FALSE)

for(BABACI in c("BA", "BACI")){
  if(BABACI=="BA"){
    BAResultsCat <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv"))
    BAResultsCat$SiteCode <- str_split_fixed(BAResultsCat$SiteSpec, "[.]", 2)[,1]
    BAResultsCat <- merge(BAResultsCat, SpeciesSiteCovs, by=c("SiteCode", "Species"))
  } else {
    BAResultsCat <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
    BAResultsCat <- merge(BAResultsCat, SpeciesSiteCovs, by=c("SiteCode", "Species"))
  }
  ModelSpecsData <- function(ModelData, PAArea){
    ModelData$Note <- ""
    
    for(var in c("MigStatus", "RedList", "Ramsar", "Family")){
      if(min(subset(as.data.frame(table(ModelData[,c(var)])), Freq>0)$Freq)<4){
        ModelData$Note <- paste0(ModelData$Note, paste0(subset(as.data.frame(table(ModelData[,c(var)])), Freq==1)$Var1, collapse=","), " from ", var, " removed_")
        ModelData <- ModelData[!ModelData[,c(var)] %in% subset(as.data.frame(table(ModelData[,c(var)])), Freq==1)$Var1,]
      }
    }
    
    if(PAArea==FALSE){
      ModelData$PAAreaLog <- NULL
    }
    
    Formula <- as.formula(paste("OutcomeSummary", paste(c(paste0(c(names(ModelData)[names(ModelData) %in% c("MigStatus", "RedList", "GenLength", "GovMean", "AnthMode", "Ramsar", "STATUS_YR", "PAAreaLog", "Family")])), 
                                                    "(1|ISO3)", "(1|ISO3:SiteCode)", "(1|Species)"), 
                                                  collapse = " + "), sep = " ~ "))
    return(list(ModelData, Formula))
  }
  
  BAResultsCat$SiteCode <- as.factor(BAResultsCat$SiteCode)
  BAResultsCat$OutcomeSummary <- factor(BAResultsCat$OutcomeSummary, levels=c("Failure", "Neutral", "Success"))
  BAResultsCat$STATUS_YR <- BAResultsCat$STATUS_YR-min(BAResultsCat$STATUS_YR)
  write.csv(BAResultsCat, paste0(MatchedFP, "Drivers/", BABACI, "_CategoryModelsData.csv"), row.names = FALSE)
  
  #get just the overlapping sitespecs
  BACICategories <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
  BACICategories$SiteSpec <- paste0(BACICategories$SiteCode,".", BACICategories$Species)
  BACategories <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv"))
  BAResultsCat$SiteSpec <- paste0(BAResultsCat$SiteCode,".", BAResultsCat$Species)
  BAResultsCatOverlap <- BAResultsCat[BAResultsCat$SiteSpec %in% intersect(BACICategories$SiteSpec, BACategories$SiteSpec),]
  
  BAResultsCatOverlapModelPA <- ModelSpecsData(BAResultsCatOverlap[!is.na(BAResultsCatOverlap$PAAreaLog),], PAArea=TRUE)
  BAResultsCatOverlapModel <- ModelSpecsData(BAResultsCatOverlap, PAArea=FALSE)
  
  AllEffectsNoPAOverlap <- clmm(BAResultsCatOverlapModel[[2]], data=BAResultsCatOverlapModel[[1]])
  AllEffectsPAOverlap <- clmm(BAResultsCatOverlapModelPA[[2]], data=BAResultsCatOverlapModelPA[[1]])
  
  BAResultsCatModelPA <- ModelSpecsData(BAResultsCat[!is.na(BAResultsCat$PAAreaLog),], PAArea=TRUE)
  BAResultsCatModel <- ModelSpecsData(BAResultsCat, PAArea=FALSE)
  
  AllEffectsNoPA <- clmm(BAResultsCatModel[[2]], data=BAResultsCatModel[[1]]) #Family + 
  AllEffectsPA <- clmm(BAResultsCatModelPA[[2]], data=BAResultsCatModelPA[[1]]) #Family + 
  
  AllEffectsSummaryOverlap <- as.data.frame(summary(AllEffectsPAOverlap)[[1]])
  AllEffectsSummaryOverlap$Coef <- row.names(AllEffectsSummaryOverlap)
  AllEffectsSummaryOverlap$Model <- "AllOverlap"
  
  AllEffectsNoPASummaryOverlap <- as.data.frame(summary(AllEffectsPAOverlap)[[1]])
  AllEffectsNoPASummaryOverlap$Coef <- row.names(AllEffectsNoPASummaryOverlap)
  AllEffectsNoPASummaryOverlap$Model <- "NoPAOverlap"
  
  AllEffectsSummary <- as.data.frame(summary(AllEffectsPA)[[1]])
  AllEffectsSummary$Coef <- row.names(AllEffectsSummary)
  AllEffectsSummary$Model <- "All"
  
  AllEffectsNoPASummary <- as.data.frame(summary(AllEffectsNoPA)[[1]])
  AllEffectsNoPASummary$Coef <- row.names(AllEffectsNoPASummary)
  AllEffectsNoPASummary$Model <- "NoPA"
  
  AllEffectsSummary <- rbindlist(list(AllEffectsSummary, AllEffectsNoPASummary, AllEffectsSummaryOverlap, AllEffectsNoPASummaryOverlap))
  names(AllEffectsSummary) <- c("Estimate", "Error", "z", "p", "Coef", "Model")
  AllEffectsSummary$OverlapNotes <- unique(BAResultsCatOverlapModel[[1]]$Note)
  AllEffectsSummary$OverlapPANotes <- unique(BAResultsCatOverlapModelPA[[1]]$Note)
  AllEffectsSummary$Notes <- unique(BAResultsCatModel[[1]]$Note)
  AllEffectsSummary$PANotes <- unique(BAResultsCatModelPA[[1]]$Note)
  
  write.csv(AllEffectsSummary, paste0(MatchedFP, "Drivers/", BABACI, "_CategoryModelsOutput.csv"), row.names = FALSE)
}

hey <- read.csv(paste0(MatchedFP, "Drivers/BACI_CategoryModelsOutput.csv"))
hey$EstP <- paste0(round(hey$Estimate, 2), "_", round(hey$p, 2), ifelse(hey$p<0.05, "*", ""))

hey <- read.csv(paste0(MatchedFP, "Drivers/BA_CategoryModelsOutput.csv"))
hey$EstP <- paste0(round(hey$Estimate, 2), "_", round(hey$p, 2), ifelse(hey$p<0.05, "*", ""))

#### Plot Category Data by Covariate ####
### try again ###

AllEffectsPA <- clmm(OutcomeSummary ~ MigStatus + GovMean + STATUS_YR + RedList + 
                       GenLength + AnthMode + Ramsar + PAAreaLog + (1 | ISO3) + 
                       (1 | ISO3:SiteCode) + (1 | Species), data=BAResultsCatModelPA[[1]]) #Family + 
ModelData <- BAResultsCatModelPA[[1]]
try <- fitted(AllEffectsPA)
ordinal::predict.clm(AllEffectsPA)

library(MCMCglmm)
AllEffectsPA <- MCMCglmm(OutcomeSummary ~ MigStatus + Family + GovMean + STATUS_YR + RedList + 
                           GenLength + AnthMode + Ramsar + PAAreaLog, random=~ISO3 + Species, data=BAResultsCatModelPA[[1]], family="ordinal") #Family + 

ModelData <- BAResultsCatModelPA[[1]]
ModelData$GovMean <- min(ModelData$GovMean)
#ModelData$GenLength <- mean(ModelData$GenLength)
#ModelData$PAAreaLog <- mean(ModelData$PAAreaLog)
#ModelData$Family <- "Anatidae"
#ModelData$STATUS_YR <- median(ModelData$STATUS_YR)
#ModelData$STATUS_YR <- round_any(ModelData$STATUS_YR, 10)
#ModelData$MigStatus <- "Migrant"
#ModelData$RedList <- "LC"
#ModelData$AnthMode <- "Croplands"
#ModelData$Ramsar <- "Other"
ModelData <- ModelData[,c("OutcomeSummary", "ISO3", "Species", "MigStatus", "Family", "GovMean", "STATUS_YR", "RedList", "GenLength", "AnthMode", "Ramsar", "PAAreaLog")]
ModelData <- ModelData[complete.cases(ModelData),]
ModelData$OutcomeSummary <- 0
#ModelData$Residuals <- predict.MCMCglmm(AllEffectsPA, ModelData)+residuals.MCMCglmm(AllEffectsPA)
ModelData$Trend <- predict.MCMCglmm(AllEffectsPA, ModelData, marginal = AllEffectsPA$Random$formula)



predict.MCMCglmm(AllEffectsPA)

### Getting better at plots ###
#first, get partial residuals
library(ggeffects)
try <- ggpredict(AllEffectsPAOverlap, c("STATUS_YR"), type="re")
try <- ggpredict(AllEffectsPA, c("PAAreaLog", "Ramsar"), type="re")
plot(ggpredict(AllEffectsPA, c("PAAreaLog", "Ramsar"), type="re"))

#The groups - Species, ISO3, SiteCode
BAResultsCat <- BAResultsCat[!is.na(BAResultsCat$OutcomeSummary),]

BAResultsCat$PAAreaLogRound <- round_any(BAResultsCat$PAAreaLog, 5)

tryish <- dcast(as.data.table(BAResultsCat), Species + PAAreaLogRound + Ramsar~OutcomeSummary, length, value.var="SiteSpec")
tryish$FailureProp <- tryish$Failure/rowSums(tryish[,c(4:6)])
tryish$NeutralProp <- tryish$Neutral/rowSums(tryish[,c(4:6)])
tryish$SuccessProp <- tryish$Success/rowSums(tryish[,c(4:6)])

tryish <- tryish[,c("Species", "PAAreaLogRound", "Ramsar", "FailureProp", "NeutralProp", "SuccessProp")]
tryish <- melt(tryish, id.vars=c("Species", "PAAreaLogRound", "Ramsar"), variable.name="response.level", value.name="predicted")
names(tryish) <- c("Species", "x", "group", "response.level", "predicted")
tryish$Data <- "Raw"

PredictData <- as.data.frame(try)
PredictData$Data <- "Modelled"
PredictData <- rbind(PredictData, tryish, use.names=TRUE, fill=TRUE)
PredictData$response.level <- as.factor(PredictData$response.level)
PredictData$response.level <- recode_factor(PredictData$response.level, "1"="Failure", "2"="Neutral", "3"="Success", "FailureProp"="Failure", "NeutralProp" = "Neutral", "SuccessProp"= "Success")
PredictData <- subset(PredictData, predicted!="NaN")

ggplot()+
  geom_line(aes(x=x, y=predicted, group=Species, colour=Species), data=subset(PredictData, Data=="Raw"), alpha=0.5)+
  geom_line(aes(x=x, y=predicted, group=group, colour=group), data=subset(PredictData, Data=="Modelled"))+
  facet_wrap(~response.level)


BACICategorySigCovs <- ggpredict(AllEffectsPA, c("PAAreaLog", "Ramsar"), type="re")








try2 <- ggpredict(AllEffectsPA, c("Ramsar"), type="re")

Colours <- c("dodgerblue4", "grey90", "red4")
CategoryPlots <- function(BABACI, Colours){
  CategoryFigFP <- paste0(MatchedFP, "Drivers/", BABACI, "CategoryPlots/")
  dir.create(file.path(CategoryFigFP), showWarnings = FALSE)
  
  CategoryData <- fread(paste0(MatchedFP, "Drivers/", BABACI, "_CategoryModelsData.csv"))
  CategoryData <- CategoryData[complete.cases(CategoryData),]
  
  #Mig
  CategoryDataMig <- dcast(CategoryData, MigStatus~OutcomeSummary, length, value.var="SiteSpec")
  CategoryDataMig <- melt(CategoryDataMig, id.vars="MigStatus")
  CategoryDataMig$MigStatus <- factor(CategoryDataMig$MigStatus, levels=c("Resident", "Migrant"))
  MigWidths <- BarPlotWidths(CategoryDataMig, Log=FALSE, Cov="MigStatus", xlabel="", AngledxLabs=FALSE, ColourList=Colours, SuccessValues, MainResults=FALSE, aspectrat=0.5)
  ggsave(paste0(CategoryFigFP, "Mig.pdf"), MigWidths, device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf
  
  #Anthrome
  CategoryDataAnth <- dcast(CategoryData, AnthMode~OutcomeSummary, length, value.var="SiteSpec")
  CategoryDataAnth <- melt(CategoryDataAnth, id.vars="AnthMode")
  CategoryDataAnth$AnthMode <- factor(CategoryDataAnth$AnthMode, levels=c("Urban", "Village", "Croplands", "Rangeland", "Semi-natural", "Wild"))
  AnthPlot <- BarPlotWidths(CategoryDataAnth, Log=TRUE, Cov="AnthMode", xlabel="", AngledxLabs=FALSE, ColourList=Colours, SuccessValues, MainResults=FALSE, aspectrat=0.5)
  ggsave(paste0(CategoryFigFP, "Anthrome.pdf"), AnthPlot, device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf
  
  #Family
  CategoryDataFamily <- dcast(CategoryData, Family~OutcomeSummary, length, value.var="SiteSpec")
  CategoryDataFamily <- melt(CategoryDataFamily, id.vars="Family")
  CategoryDataFamily$Family <- factor(CategoryDataFamily$Family)
  FamilyPlot <- BarPlotWidths(CategoryDataFamily, Log=TRUE, Cov="Family", xlabel="", AngledxLabs=TRUE, ColourList=Colours, SuccessValues, MainResults=FALSE, aspectrat=0.5)
  ggsave(paste0(CategoryFigFP,"Family.pdf"), FamilyPlot, device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf
  
  #IUCN
  CategoryDataIUCN <- dcast(CategoryData, RedList~OutcomeSummary, length, value.var="SiteSpec")
  CategoryDataIUCN <- melt(CategoryDataIUCN, id.vars="RedList")
  CategoryDataIUCN$RedList <- factor(CategoryDataIUCN$RedList, levels=c('LC', 'NT', 'VU', 'EN', 'CR'))
  CategoryDataIUCN$RedList <- recode_factor(CategoryDataIUCN$RedList, 'LC'='Least Concern', 'NT'='Near Threatened', 'VU'='Vulnerable', 'EN'='Endangered', 'CR'='Critically Endangered')
  IUCNPlot <- BarPlotWidths(CategoryDataIUCN, Log=TRUE, Cov="RedList", xlabel="", AngledxLabs=FALSE, ColourList=Colours, SuccessValues, MainResults=FALSE, aspectrat=0.5)
  ggsave(paste0(CategoryFigFP,"IUCN.pdf"), IUCNPlot, device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf
  
  #Ramsar
  CategoryDataRamsar <- dcast(CategoryData, Ramsar~OutcomeSummary, length, value.var="SiteSpec")
  CategoryDataRamsar <- melt(CategoryDataRamsar, id.vars="Ramsar")
  CategoryDataRamsar$Ramsar <- factor(CategoryDataRamsar$Ramsar, levels=c("Other", "Ramsar"))
  RamsarWidths <- BarPlotWidths(CategoryDataRamsar, Log=FALSE, Cov="Ramsar", xlabel="", AngledxLabs=FALSE, ColourList=Colours, SuccessValues, MainResults=FALSE, aspectrat=0.5)
  ggsave(paste0(CategoryFigFP, "Ramsar.pdf"), RamsarWidths, device="pdf", width = 300, height = 150, units = "mm") #Save as a pdf
  
  #Category plots
  CategoryData$OutcomeSummary <- factor(CategoryData$OutcomeSummary, levels=c("Failure", "Neutral", "Success"))
  CategoryData <- as.data.frame(CategoryData)
  
  BoxPlotFigure <- function(BoxData, ylabel, Cov){
    BoxData$Cov <- BoxData[,c(Cov)]
    BoxData <- BoxData[,c("OutcomeSummary", "Cov")]
    ggplot(data=BoxData, aes(x=OutcomeSummary, y=Cov))+
      #geom_beeswarm(size=0.8, cex=0.2, colour="grey60", grouponX=FALSE)+
      geom_boxplot(fill=rev(Colours), alpha=0.8)+
      #coord_flip()+
      ylab(ylabel)+
      theme(
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=10, colour = "black"),
        axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=10, colour = "black"))
  }
  
  #Log Area
  AreaPlot <- BoxPlotFigure(CategoryData, ylabel="Log (Area) km2", Cov="PAAreaLog")
  ggsave(paste0(CategoryFigFP,"Area.pdf"), AreaPlot, device="pdf", width = 150, height = 150, units = "mm") #Save as a pdf
  
  #Gen
  GenPlot <- BoxPlotFigure(CategoryData, ylabel="Generation Length (years)", Cov="GenLength")
  ggsave(paste0(CategoryFigFP, "Gen.pdf"), GenPlot, device="pdf", width = 150, height = 150, units = "mm") #Save as a pdf
  
  #Gov
  GovernancePlot <- BoxPlotFigure(CategoryData, ylabel="Governance Index", Cov="GovMean")
  ggsave(paste0(CategoryFigFP, "Governance.pdf"), GovernancePlot, device="pdf", width = 150, height = 150, units = "mm") #Save as a pdf
}
CategoryPlots("BA", Colours)
CategoryPlots("BACI", Colours)

BACategories <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACategories.csv"))
BACategories$SiteCode <- str_split_fixed(BACategories$SiteSpec, "[.]", 2)[,1]

BACICategories <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
BACICategories$SiteSpec <- paste0(BACICategories$SiteCode,".", BACICategories$Species)
BACICategories <- BACICategories[,names(BACICategories) %in% names(BACategories)]

BACICategories$Counterfactual <- "BACI"
BACategories$Counterfactual <- "BA"

BACICategories <- BACICategories[BACICategories$SiteSpec %in% intersect(BACICategories$SiteSpec, BACategories$SiteSpec),]
BACategories <- BACategories[BACategories$SiteSpec %in% intersect(BACICategories$SiteSpec, BACategories$SiteSpec),]

Join <- rbind(BACICategories, BACategories)
Join <- subset(Join, OutcomeSummary!="Excluded")
JoinCast <- dcast(as.data.table(Join), SiteSpec~Counterfactual, value.var="OutcomeSummary")
JoinCast$Change <- ifelse(JoinCast$BA==JoinCast$BACI, "NoChange", ifelse(JoinCast$BA=="Success" & JoinCast$BACI=="Failure", "SuccesstoFailure",
                                                                         ifelse(JoinCast$BA=="Success" & JoinCast$BACI=="Neutral", "SuccesstoNeutral",
                                                                                ifelse(JoinCast$BA=="Failure" & JoinCast$BACI=="Success", "FailuretoSuccess", 
                                                                                       ifelse(JoinCast$BA=="Failure" & JoinCast$BACI=="Neutral", "FailuretoNeutral",
                                                                                              ifelse(JoinCast$BA=="Neutral" & JoinCast$BACI=="Success", "NeutraltoSuccess", "NeutraltoFailure"))))))

#### Covariate Models by Mean/Slope  ####
for(BABACI in c("BA", "BACI")){
  if(BABACI=="BA"){
    BAMeanSlopeModel <- read.csv(file=paste0(MatchedFP, "ModelOutput/BAMeanSlope.csv"))
    BAMeanSlopeModel$Species <- str_split_fixed(BAMeanSlopeModel$SiteSpec, "[.]", 2)[,2]
    BAMeanSlopeModel$SiteCode <- str_split_fixed(BAMeanSlopeModel$SiteSpec, "[.]", 2)[,1]
    MeanSlopeResults <- merge(BAMeanSlopeModel, SpeciesSiteCovs, by=c("SiteCode", "Species"))
  } else {
    BACIMeanSlopeModel <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACIMeanSlope.csv"))
    BACIMeanSlopeModel <- merge(BACIMeanSlopeModel, unique(subset(BACI, CI==1)[,c("SpecMatch", "SiteCode")]))
    BACIMeanSlopeModel$Species <- str_split_fixed(BACIMeanSlopeModel$SpecMatch, "[.]", 2)[,1]
    BACIMeanSlopeModel$SiteSpec <- paste0(BACIMeanSlopeModel$SiteCode, ".", BACIMeanSlopeModel$Species)
    check <- nrow(BACIMeanSlopeModel)
    BACIMeanSlopeModel <- merge(BACIMeanSlopeModel, SpeciesSiteCovs, by=c("SiteCode", "Species"))
    if(nrow(BACIMeanSlopeModel)/check<0.95){stop("SiteSpec lost!!!!")}
    MeanSlopeResults <- BACIMeanSlopeModel
  }
  
  MeanSlopeResults$SiteCode <- as.factor(MeanSlopeResults$SiteCode)
  MeanSlopeResults$STATUS_YR <- MeanSlopeResults$STATUS_YR-min(MeanSlopeResults$STATUS_YR)
  write.csv(MeanSlopeResults, paste0(MatchedFP, "Drivers/", BABACI, "_MeanSlopeModelsData.csv"), row.names = FALSE)
  
  PredictorModelAll <- function(MSResults, Response, Signif, overlap){
    if(Response=="Slope"){
      MSResults <- subset(MSResults, Model=="Slope" & P<Signif)
    } else {
      MSResults <- subset(MSResults, Model=="Mean" & P<Signif)
    }
    MSResults$ID <- row.names(MSResults)
    if(overlap==TRUE){
      BACIMeanSlopeModel <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACIMeanSlope.csv"))
      BACIMeanSlopeModel <- merge(BACIMeanSlopeModel, unique(subset(BACI, CI==1)[,c("SpecMatch", "SiteCode")]))
      BACIMeanSlopeModel$Species <- str_split_fixed(BACIMeanSlopeModel$SpecMatch, "[.]", 2)[,1]
      BACIMeanSlopeModel$SiteSpec <- paste0(BACIMeanSlopeModel$SiteCode, ".", BACIMeanSlopeModel$Species)
      
      BAMeanSlopeModel <- read.csv(file=paste0(MatchedFP, "ModelOutput/BAMeanSlope.csv"))
      
      MSResults <- MSResults[MSResults$SiteSpec %in% intersect(BACIMeanSlopeModel$SiteSpec, BAMeanSlopeModel$SiteSpec),]
    }
    ModelSpecsData <- function(ModelData, PAArea){
      ModelData$Note <- ""
      
      for(var in c("MigStatus", "RedList", "Ramsar", "Family")){
        if(min(subset(as.data.frame(table(ModelData[,c(var)])), Freq>0)$Freq)<4){
          ModelData$Note <- paste0(ModelData$Note, paste0(subset(as.data.frame(table(ModelData[,c(var)])), Freq<4)$Var1, collapse=","), " from ", var, " removed_")
          ModelData <- ModelData[!ModelData[,c(var)] %in% subset(as.data.frame(table(ModelData[,c(var)])), Freq<4)$Var1,]
        }
      }
      
      if(PAArea==FALSE){
        ModelData$PAAreaLog <- NULL
      }
      
      Formula <- as.formula(paste("Estimate", paste(c(paste0(c(names(ModelData)[names(ModelData) %in% c("MigStatus", "RedList", "GenLength", "GovMean", "AnthMode", "Ramsar", "STATUS_YR", "Family", "PAAreaLog")])), 
                                                      "(1|ISO3)", "(1|ISO3:SiteCode)", "(1|Species)"), 
                                                    collapse = " + "), sep = " ~ "))
      return(list(ModelData, Formula))
    }
    
    MSResults2 <- ModelSpecsData(MSResults[!is.na(MSResults$PAAreaLog),], PAArea=TRUE)
    PredModelPA <- glmmTMB(MSResults2[[2]], data=MSResults2[[1]])
    PredModelSummaryPA <- as.data.frame(summary(PredModelPA)[[6]]$cond)
    PredModelSummaryPA$Coef <- row.names(PredModelSummaryPA)
    PredModelSummaryPA$Model <- paste0(Response, "All")
    
    MSResults <- ModelSpecsData(MSResults, PAArea=FALSE)
    PredModelNoPA <- glmmTMB(MSResults[[2]], data=MSResults[[1]])
    PredModelSummary <- as.data.frame(summary(PredModelNoPA)[[6]]$cond)
    PredModelSummary$Coef <- row.names(PredModelSummary)
    PredModelSummary$Model <- paste0(Response, "NoPA")
    
    PredModelSummary <- rbindlist(list(PredModelSummaryPA, PredModelSummary))
    names(PredModelSummary) <- c("Estimate", "Error", "z", "P", "Coef", "Model")
    PredModelSummary$Overlap <- overlap
    PredModelSummary$Signif <- Signif
    MSResultsSummary <- MSResults[[1]]
    MSResultsSummary$FittedNoPA <- fitted.values(PredModelNoPA)
    MSResultsSummaryPA <- MSResults2[[1]]
    MSResultsSummaryPA$FittedPA <- fitted.values(PredModelPA)
    MSResultsSummary <- merge(MSResultsSummary, MSResultsSummaryPA[,c("ID", "FittedPA")], by="ID", all=T)
    MSResultsSummary$NotesPA <- unique(MSResults2[[1]]$Note)
    MSResultsSummary$Overlap <- overlap
    MSResultsSummary$Signif <- Signif
    return(list(PredModelSummary, MSResultsSummary))
  }
  
  MeanOutputOverlap <- PredictorModelAll(MeanSlopeResults, "Mean", 1, overlap=TRUE)
  SlopeOutputOverlap <- PredictorModelAll(MeanSlopeResults, "Slope", 1, overlap=TRUE)
  MeanOutput <- PredictorModelAll(MeanSlopeResults, "Mean", 1, overlap=FALSE)
  SlopeOutput <- PredictorModelAll(MeanSlopeResults, "Slope", 1, overlap=FALSE)
  MeanOutputOverlapP <- PredictorModelAll(MeanSlopeResults, "Mean", 0.05, overlap=TRUE)
  SlopeOutputOverlapP <- PredictorModelAll(MeanSlopeResults, "Slope", 0.05, overlap=TRUE)
  MeanOutputP <- PredictorModelAll(MeanSlopeResults, "Mean", 0.05, overlap=FALSE)
  SlopeOutputP <- PredictorModelAll(MeanSlopeResults, "Slope", 0.05, overlap=FALSE)
  
  
  MeanSlopeOutput <- rbindlist(list(MeanOutput[[1]], SlopeOutput[[1]], MeanOutputOverlap[[1]], SlopeOutputOverlap[[1]], MeanOutputOverlapP[[1]], SlopeOutputOverlapP[[1]], MeanOutputP[[1]], SlopeOutputP[[1]]))
  MeanData <- rbindlist(list(MeanOutput[[2]], MeanOutputOverlap[[2]], MeanOutputOverlapP[[2]], MeanOutputP[[2]]))
  SlopeData <- rbindlist(list(SlopeOutput[[2]], SlopeOutputOverlap[[2]], SlopeOutputOverlapP[[2]], SlopeOutputP[[2]]))
  
  write.csv(MeanSlopeOutput, paste0(MatchedFP, "Drivers/", BABACI, "_MeanSlopeModels.csv"), row.names = FALSE)
  write.csv(MeanData, paste0(MatchedFP, "Drivers/", BABACI, "_MeanData.csv"), row.names = FALSE)
  write.csv(SlopeData, paste0(MatchedFP, "Drivers/", BABACI, "_SlopeData.csv"), row.names = FALSE)
}

hey <- read.csv(paste0(ResultsFP, "Drivers/BACI_MeanSlopeModels.csv"))
hey$EstP <- paste0(round(hey$Estimate, 2), "_", round(hey$P, 2), ifelse(hey$P<0.05, "*", ""))

hey <- read.csv(paste0(ResultsFP, "Drivers/BA_MeanSlopeModels.csv"))
hey$EstP <- paste0(round(hey$Estimate, 2), "_", round(hey$P, 2), ifelse(hey$P<0.05, "*", ""))

#### Plot Mean/Slope by covariate UP TO HERE WITH UPDATING FILE PATHS####
#Testing

#CHECK WITH ALI IS IT RIGHT TO JUST USE RAW COEFFICIENT

MSResults2 <- ModelSpecsData(MSResults[!is.na(MSResults$PAAreaLog),], PAArea=TRUE)
PredModelPA <- glmmTMB(MSResults2[[2]], data=MSResults2[[1]])

PredModelPALME <- lmer(Estimate ~ MigStatus + Family + GovMean + STATUS_YR + RedList + 
                         GenLength + AnthMode + Ramsar + PAAreaLog + STATUS_YR + (1 | Species), data=MSResults2[[1]], na.action="na.fail")

PredModelPA <- glmmTMB(Estimate ~ MigStatus + GovMean + STATUS_YR + RedList + 
                         GenLength + AnthMode + Ramsar + PAAreaLog + (STATUS_YR | Family), data=MSResults2[[1]])

summary(PredModelPA)
#Anova(PredModelPALME)
#dredge(PredModelPALME)

ModelData <- MSResults2[[1]]

ModelData$GovMean <- mean(ModelData$GovMean)
ModelData$GenLength <- mean(ModelData$GenLength)
ModelData$PAAreaLog <- mean(ModelData$PAAreaLog)
#ModelData$Family <- "Anatidae"
#ModelData$STATUS_YR <- median(ModelData$STATUS_YR)
#ModelData$STATUS_YR <- round_any(ModelData$STATUS_YR, 10)
ModelData$MigStatus <- "Migrant"
ModelData$RedList <- "Not Threatened"
ModelData$AnthMode <- "Croplands"
ModelData$Ramsar <- "Other"

ModelData$Residuals <- predict(PredModelPA, ModelData)+residuals(PredModelPA)
ModelData$Trend <- predict(PredModelPA, ModelData)

right <- merge(BAMeanSlopeModel, SpeciesSiteCovs, by=c("SiteCode", "Species"))
right <- unique(right[,c("SiteCode", "STATUS_YR")])
names(right) <- c("SiteCode", "StatYr")
ModelData <- merge(ModelData, right, by="SiteCode")

ggplot(data=ModelData)+
  geom_point(aes(x=StatYr, y=Residuals, group=Family, colour=Family))+
  geom_line(aes(x=StatYr, y=Trend, group=Family, colour=Family))+
  #facet_grid(Ramsar~RedList)+
  geom_hline(yintercept=1, colour="black", linetype="dashed")+
  ylim(0,500)+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(), 
    text = element_text(size=10, colour = "black"),
    axis.text.x = element_text(size=10, colour = "black"), #
    axis.text.y = element_text(size=10, colour = "black"),
    panel.border = element_rect(size = 1, fill = NA),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, size=10, colour = "black"))


#### Gov and Status YR
PredModelPA <- glmmTMB(MSResults2[[2]], data=MSResults2[[1]])



PredModelPA <- glmmTMB(Estimate ~ MigStatus + Family + GovMean + StatYrRel + RedList + 
                         GenLength + AnthMode + Ramsar + PAAreaLog + GovMean*StatYrRel + (1 | Species), data=MSResults2[[1]])

summary(PredModelPA)

ModelData <- MSResults2[[1]]

ModelData$GovMean <- mean(ModelData$GovMean)
ModelData$GenLength <- mean(ModelData$GenLength)
ModelData$PAAreaLog <- mean(ModelData$PAAreaLog)
ModelData$Family <- "Anatidae"
ModelData$StatYrRel <- median(ModelData$StatYrRel)
#ModelData$STATUS_YR <- round_any(ModelData$STATUS_YR, 10)
#ModelData$StatYrRel <- round_any(ModelData$StatYrRel, 10)
ModelData$MigStatus <- "Migrant"
#ModelData$RedList <- "LC"
ModelData$AnthMode <- "Croplands"
#ModelData$Ramsar <- "Other"

ModelData$Residuals <- predict(PredModelPA, ModelData)+residuals(PredModelPA)
ModelData$Trend <- predict(PredModelPA, ModelData)

ggplot(data=ModelData)+
  geom_boxplot(aes(x=Ramsar, y=Residuals))+
  facet_wrap(~RedList)+
  geom_hline(yintercept=0, colour="black", linetype="dashed")+
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(), 
    text = element_text(size=10, colour = "black"),
    axis.text.x = element_text(size=10, colour = "black"), #
    axis.text.y = element_text(size=10, colour = "black"),
    panel.border = element_rect(size = 1, fill = NA),
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0, size=10, colour = "black"))



library(ggeffects)
plot(ggpredict(PredModelPA, c("PAAreaLog", "Ramsar", "Species")))



MeanSlopePlots <- function(BABACI, MeanSlope, ylabel){
  CategoryFigFP <- paste0(MatchedFP, "Drivers/", BABACI, "CategoryPlots/")
  dir.create(file.path(CategoryFigFP), showWarnings = FALSE)
  
  if(MeanSlope=="Slope"){
    MeanSlopeResultsPlotty <- fread(paste0(ResultsFP, "Drivers/", BABACI, "_SlopeData.csv"))
  } else {
    MeanSlopeResultsPlotty <- fread(paste0(ResultsFP, "Drivers/", BABACI, "_MeanData.csv"))
  }
  
  Boxplotplot <- function(MeanSlopeResultsPlottyPlot, Coef, xlabel, Family=FALSE){
    MeanSlopeResultsPlottyPlot <- as.data.frame(MeanSlopeResultsPlottyPlot)
    MeanSlopeResultsPlottyPlot$Coef <- MeanSlopeResultsPlottyPlot[,c(Coef)]
    MeanSlopeResultsPlottyPlot <- MeanSlopeResultsPlottyPlot[,c("Coef", "Estimate", "Fitted")]
    names(MeanSlopeResultsPlottyPlot) <- c("Coef", "Actual Values", "Fitted Values")
    MeanSlopeResultsPlottyPlot <- melt(as.data.table(MeanSlopeResultsPlottyPlot), id.vars="Coef")
    ggplot(data=MeanSlopeResultsPlottyPlot, aes(x=factor(Coef), y=log(value), fill=factor(variable)))+
      #geom_beeswarm(cex=0.2, size=0.2)+
      geom_boxplot(outlier.size=0.5, alpha=0.9)+
      #geom_boxplot(aes(y=log(Fitted)), outlier.size=0.1, fill="grey10", colour="grey10", alpha=0.5)+
      scale_fill_manual(values=c("#56B4E9", "#0072B2"), name="Data")+
      ylab(ylabel)+
      xlab(xlabel)+
      geom_hline(yintercept=0, colour="black", linetype="dashed")+
      theme(
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=10, colour = "black"),
        axis.text.x = element_text(size=10, colour = "black", angle = ifelse(Family==FALSE, 0, 45)), #
        axis.text.y = element_text(size=10, colour = "black"),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=10, colour = "black"))
  }
  ContPlot <- function(Coef, xlabel){
    MeanSlopeResultsPlottyPlot <- as.data.frame(MeanSlopeResultsPlotty)
    MeanSlopeResultsPlottyPlot$Coef <- MeanSlopeResultsPlottyPlot[,c(Coef)]
    
    if(Coef=="GenLength"){
      MeanSlopeResultsPlottyPlot$Species <- MeanSlopeResultsPlottyPlot$SiteCode
    }
    
    f <- as.formula(paste("Estimate", paste(c(Coef, "(1|Species)"), collapse = " + "), sep = " ~ "))
    PredModel <- lme4::lmer(f, data=MeanSlopeResultsPlottyPlot)
    MeanSlopeResultsPlottyPlot$Fitted <- fitted(PredModel)
    
    ggplot(data=MeanSlopeResultsPlottyPlot, aes(x=Coef, y=Estimate))+ #, fill=Result
      #geom_linerange(aes(ymin=log(MeanMean-SE), ymax=log(MeanMean+SE)), size=0.3, colour="grey60")+
      geom_point(colour="grey90")+
      geom_line(aes(x=Coef, y=Fitted, group=Species), alpha=0.5)+
      ylab(ylabel)+
      xlab(xlabel)+
      geom_hline(yintercept=0, colour="grey60", linetype="dashed")+
      #geom_line(aes(x=Coef, y=log(Model)), colour="red")+
      #geom_errorbarh(aes(xmin=Mean-MeanError, xmax=Mean+MeanError), width=.2)+
      theme(aspect.ratio=1, 
            panel.background = element_blank(),
            panel.grid = element_blank(), 
            text = element_text(size=plotfontsize),
            panel.border = element_rect(size = 1, fill = NA),
            strip.background = element_blank(),
            strip.text = element_text(hjust = 0, size=plotfontsize),
            legend.title=element_text(size=plotfontsize),
            axis.text.x = element_text(size=plotfontsize, colour="black"), #element_text(size=plotfontsize, colour="black", angle = 90)
            axis.ticks = element_line(colour="black"),
            legend.text = element_text(size=plotfontsize, colour="black"),
            legend.justification = "top")
  }
  
  #Mig
  MeanSlopeResultsPlottyMig <- MeanSlopeResultsPlotty
  MeanSlopeResultsPlottyMig$MigStatus <- factor(MeanSlopeResultsPlottyMig$MigStatus, levels=c("Resident", "Migrant"))
  MeanSlopeMig <- Boxplotplot(MeanSlopeResultsPlottyMig, "MigStatus", "Migratory Status")
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_Mig.pdf"), MeanSlopeMig, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
  
  #Anthrome
  MeanSlopeResultsPlottyAnth <- MeanSlopeResultsPlotty
  MeanSlopeResultsPlottyAnth$AnthMode <- factor(MeanSlopeResultsPlottyAnth$AnthMode, levels=c("Urban", "Village", "Croplands", "Rangeland", "Semi-natural", "Wild"))
  MeanSlopeAnth <- Boxplotplot(MeanSlopeResultsPlottyAnth, "AnthMode", "Land Category")
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_Anth.pdf"), MeanSlopeAnth, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
  
  #Family
  MeanSlopeFamily <- Boxplotplot(MeanSlopeResultsPlotty, "Family", "Family", Family=TRUE)
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_Family.pdf"), MeanSlopeFamily, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
  
  #IUCN
  MeanSlopeResultsPlottyIUCN <- MeanSlopeResultsPlotty
  MeanSlopeResultsPlottyIUCN$RedList <- factor(MeanSlopeResultsPlottyIUCN$RedList, levels=c('LC', 'NT', 'VU', 'EN', 'CR'))
  MeanSlopeResultsPlottyIUCN$RedList <- recode_factor(MeanSlopeResultsPlottyIUCN$RedList, 'LC'='Least Concern', 'NT'='Near Threatened', 'VU'='Vulnerable', 'EN'='Endangered', 'CR'='Critically Endangered')
  MeanSlopeIUCN <- Boxplotplot(MeanSlopeResultsPlottyIUCN, "RedList", "Red list category")
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_IUCN.pdf"), MeanSlopeIUCN, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
  
  #Area
  AreaPlot <- ContPlot("PAAreaLog", "Log (PA Area)")
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_Area.pdf"), AreaPlot, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
  
  #GenLength
  GenLengthPlot <- ContPlot("GenLength", "Generation length (years)")
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_GenLength.pdf"), GenLengthPlot, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
  
  #Governance
  GenLengthPlot <- ContPlot("GovMean", "Governance Index")
  ggsave(paste0(FiguresFP,BABACI, "/MeanSlope/", MeanSlope, "_Governance.pdf"), GenLengthPlot, device="pdf", width = 200, height = 150, units = "mm") #Save as a pdf
}

#### Covariate sensitivity summarise ####

#Get model output
BACIMeanSlope <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BACI_MeanSlopeModels.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  Dat$Scenario <- paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11])
  Dat$BABACI <- "BACI"
  Dat$Framework <- "MeanSlope"
  return(Dat)
}))
BACIMeanSlope <- subset(BACIMeanSlope, Overlap==FALSE & Signif=="0.05")
BAMeanSlope <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BA_MeanSlopeModels.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  Dat$Scenario <- paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11])
  Dat$BABACI <- "BA"
  Dat$Framework <- "MeanSlope"
  return(Dat)
}))
BAMeanSlope <- subset(BAMeanSlope, Overlap==FALSE &Signif=="0.05")
BACICategory <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BACI_CategoryModelsOutput.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  Dat$Scenario <- paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11])
  Dat$BABACI <- "BACI"
  names(Dat)[names(Dat)=="p"] <- "P"
  Dat$Framework <- "Category"
  return(Dat)
}))
BACategory <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BA_CategoryModelsOutput.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  Dat$Scenario <- paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11])
  Dat$BABACI <- "BA"
  names(Dat)[names(Dat)=="p"] <- "P"
  Dat$Framework <- "Category"
  return(Dat)
}))

#Get number of rows
BACIMeanSlopeData <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BACI_MeanSlopeModelsData.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  return(data.frame(Scenario=paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11]), nrow=nrow(Dat[complete.cases(Dat$PAAreaLog),]), BABACI="BACI", Framework="MeanSlope",
                    Ramsar=as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar))[2]/sum(as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar)))))
}))
BAMeanSlopeData <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BA_MeanSlopeModelsData.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  return(data.frame(Scenario=paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11]), nrow=nrow(Dat[complete.cases(Dat$PAAreaLog),]), BABACI="BA", Framework="MeanSlope",
                    Ramsar=as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar))[2]/sum(as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar)))))
}))
BACICategoryData <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BACI_CategoryModelsData.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  return(data.frame(Scenario=paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11]), nrow=nrow(Dat[complete.cases(Dat$PAAreaLog),]), BABACI="BACI", Framework="Category",
                    Ramsar=as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar))[2]/sum(as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar)))))
}))
BACategoryData <- rbindlist(lapply(list.files(path = DatabaseFP, pattern = "BA_CategoryModelsData.csv$", recursive = TRUE, full.names = TRUE), function(x){
  Dat <- fread(x)
  return(data.frame(Scenario=paste0(str_split_fixed(x, "[/]", 13)[,10], "_", str_split_fixed(x, "[/]", 13)[,11]), nrow=nrow(Dat[complete.cases(Dat$PAAreaLog),]), BABACI="BA", Framework="Category",
                    Ramsar=as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar))[2]/sum(as.vector(table(unique(Dat[,c("SiteCode", "Ramsar")])$Ramsar)))))
}))

AllModels <- rbindlist(list(BACIMeanSlope[,c("Estimate", "Error", "z", "P", "Coef", "Model", "Scenario", "BABACI", "Framework", "Signif")],
                            BAMeanSlope[,c("Estimate", "Error", "z", "P", "Coef", "Model", "Scenario", "BABACI", "Framework", "Signif")],
                            BACICategory[,c("Estimate", "Error", "z", "P", "Coef", "Model", "Scenario", "BABACI", "Framework")],
                            BACategory[,c("Estimate", "Error", "z", "P", "Coef", "Model", "Scenario", "BABACI", "Framework")]), fill=TRUE)
AllModelsData <- rbindlist(list(BACIMeanSlopeData, BAMeanSlopeData, BACICategoryData, BACategoryData))

AllModels$Significant <- ifelse(AllModels$P<0.05, ifelse(AllModels$Estimate>0, 1, -1), 0)

try <- dcast(AllModels, Scenario + Model + Signif + BABACI + Framework ~ Coef, value.var="Significant")
try <- subset(try, Model=="All" | Model=="MeanAll" | Model=="SlopeAll")
try[,c("Neutral|Success", "Failure|Neutral", "(Intercept)")] <- NULL

try$Framework <- ifelse(try$Model=="MeanAll", "Mean", ifelse(try$Model=="SlopeAll", "Slope", "Category"))

okidoke <- melt(try, id.vars=c("Scenario", "Model", "Signif", "BABACI", "Framework"))
#okidoke$nrowRound <- round_any(okidoke$nrow, 10000, ceiling)
okidoke <- dcast(okidoke, variable + BABACI + Framework~value, length, value.var="Scenario")
okidoke$Neg <- round(okidoke[,5]/rowSums(okidoke[,c(4:7)]), 2)
okidoke$Insig <- round(okidoke[,6]/rowSums(okidoke[,c(4:7)]), 2)
okidoke$Pos <- round(okidoke[,7]/rowSums(okidoke[,c(4:7)]), 2)
okidoke$Total <- rowSums(okidoke[,c(4:7)])
okidoke <- melt(okidoke[,c(1:3, 8:11)], id.vars=c("variable", "BABACI", "Framework", "Total"))
names(okidoke) <- c("Coef", "Counterfactual", "Framework","Total", "Outcome", "Proportion")

okidoke2 <- okidoke[grepl("Family", okidoke$Coef),]
okidoke2 <- subset(okidoke2, Coef!="MigStatusUncertain")
okidoke2$Outcome <- factor(okidoke2$Outcome, levels=c("Insig","Neg", "Pos"))
ggplot(okidoke2, aes(x=Coef, y=Proportion, fill=Outcome))+
  geom_bar(stat="identity", position="stack")+
  facet_grid(Counterfactual~Framework)+
  scale_fill_manual(values = c("gray82", "firebrick", "steelblue4"), name=("Outcome"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Proportion")+
  xlab("Coef")+
  theme(aspect.ratio=0.5, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.ticks = element_line(colour="black"),
        axis.text.x = element_text(angle=45),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")

#Check the num rows
Var <- "RamsarRamsar"

VarSub <- subset(AllModels, Coef==Var)
VarSub <- subset(VarSub, Model=="AllOverlap" | Model=="MeanAll" | Model=="SlopeAll")

VarSub <- merge(VarSub, AllModelsData, by=c("Scenario", "BABACI", "Framework"))
VarSub$Framework <- ifelse(VarSub$Model=="MeanAll", "Mean", ifelse(VarSub$Model=="SlopeAll", "Slope", "Category"))
VarSub$Signal <- ifelse(VarSub$Estimate>0, "Pos", "Neg")
VarSub$EstRound <- round(VarSub$Estimate, 2)
ggplot(VarSub, aes(x=nrow, y=P, colour=Signal))+
  geom_point()+
  geom_text(aes(label=EstRound),hjust= -0.03, vjust= -0.03)+
  #stat_smooth(method="lm")+
  #scale_colour_distiller(palette = "Spectral", limits=c(-0.1,0.1))+
  facet_wrap(BABACI~Framework, scales="free_x")+ #
  geom_hline(yintercept = 0.05, linetype="dashed")+
  #geom_vline(xintercept = 0.05, linetype="dashed")+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.ticks = element_line(colour="black"),
        axis.text.x = element_text(angle=45),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
  


okidoke <- melt(try, id.vars=c("Scenario", "Model", "Signif", "BABACI", "Framework"))
#okidoke$nrowRound <- round_any(okidoke$nrow, 10000, ceiling)


#We want a table with 



#### Case Study plots BACI ####
#HANNAH NOTE. TO COMBINE COEFFICIENTS PRE EXPONENTIALISING THEM YOU NEED TO ADD NOT MULTIPLE

BACIByPopulation <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACIMeanSlope.csv"))

hist(subset(BACIByPopulation, Coef=="Year:BA:CI" & P<0.05)$Estimate, breaks=50)
hist(subset(BACIByPopulation, Coef=="BA:CI" & Model=="Mean" & P<0.05)$Estimate, breaks=50)

BACIByPopulation$Coef <- as.factor(BACIByPopulation$Coef)
BACIByPopulation$Coef <- recode_factor(BACIByPopulation$Coef, "Year:BA:CI" = "SlopeChange", "BA:CI" = "MeanChange")
BACIByPopulationShrink <- BACIByPopulation[,c("SpecMatch", "Coef", "Estimate", "P")]
BACIByPopulationShrink <- merge(BACIByPopulationShrink, unique(BACI[,c("SpecMatch", "Dataset", "IUCN")]))

BACIByPopulationCBC <- subset(BACIByPopulationShrink, Dataset=="CBC")
BACIByPopulationIWC <- subset(BACIByPopulationShrink, Dataset=="IWC")

BACIByPopulationSigSlope <- subset(BACIByPopulationCBC, Coef=="SlopeChange" & P<0.05 & Estimate >0)
BACIByPopulationSigMean <- subset(BACIByPopulationCBC, Coef=="MeanChange" & P<0.05 & Estimate<0 | Coef=="MeanChange" & P>0.05)

BACIByPopulationMeanOrSlope <- BACIByPopulationCBC[BACIByPopulationCBC$SpecMatch %in% BACIByPopulationSigSlope$SpecMatch,] #intersect(as.character(BACIByPopulationSigSlope$SpecMatch), as.character(BACIByPopulationSigMean$SpecMatch))
BACIByPopulationMeanOrSlope <- merge(BACIByPopulationMeanOrSlope, unique(BACI[,c("SpecMatch", "SiteCode", "Latitude", "Longitude", "CI", "WDPAID")]))
PAs <- as.data.frame(str_split_fixed(BACIByPopulationMeanOrSlope$WDPAID, "[,]", 4))
BACIByPopulationMeanOrSlope2 <- cbind(BACIByPopulationMeanOrSlope[,c(1:ncol(BACIByPopulationMeanOrSlope)-1)], PAs)
BACIByPopulationMeanOrSlope2 <- melt(BACIByPopulationMeanOrSlope2, id.vars=c("SpecMatch", "Coef", "Estimate", "P", "Dataset", "IUCN", "SiteCode", "Latitude", "Longitude", "CI"), value.name="WDPAID")
#PPData <- read.dbf(paste0(DataFP, "ProtectedPlanet/WDPA_2020/WDPA_Jan2020-shapefile/WDPA_Jan2020-shapefile-polygons.dbf"))
BACIByPopulationMeanOrSlope2 <- merge(BACIByPopulationMeanOrSlope2, PPData, by="WDPAID")

BACIByPopulationMeanOrSlope2 <- subset(BACIByPopulationMeanOrSlope2, MARINE==0 | MARINE==1)
#BACIByPopulationMeanOrSlope2 <- BACIByPopulationMeanOrSlope2[grep("\\Ramsar\\b", BACIByPopulationMeanOrSlope2$DESIG),]

BACIByPopulationMeanOrSlope2
#Larus marinus.260
#Arenaria interpres.12
#Larus argentatus.417
Solo <- "Mareca strepera.612"
Solo <- "Cygnus olor.436"
library(gridExtra)
library(cowplot)

load(file=paste0(DatabaseFP, "ProtectedCountDataCollinearRemoved.RData"))
load(file=paste0(DatabaseFP, "UnprotectedCountDataCollinearRemoved.RData"))

Solo <- "Bucephala clangula.1149"


#  write.csv(BACISolo, paste0(ResultsFP, "TREE/ExamplePopData.csv"), row.names=FALSE)
hey <- subset(BACISolo, CI==1)
SoloModSlope <- glm.nb(Count~Year2 + BA + BA*Year2 + offset(log(Hours)), link=log, data=BACISolo)
summary(SoloModSlope)


#For CBC
lapply(unique(BACIByPopulationMeanOrSlope2$SpecMatch), function(Solo){
  print(Solo)
  BACISolo <- read.csv(paste0(ResultsFP, "TREE/ExamplePopData.csv"), header=TRUE)
  BACISolo$Year2 <- BACISolo$Year-BACISolo$STATUS_YR
  SoloModSlope <- glm.nb(Count~Year2 + BA + BA*Year2 + CI + CI*BA + CI*Year2 + BA*CI*Year2 + offset(log(Hours)), link=log, data=BACISolo)
  summary(SoloModSlope)
  SoloModMean <- glm.nb(Count~BA + CI + CI*BA + offset(log(Hours)), link=log, data=BACISolo)
  summary(SoloModMean)

  BACISoloPredict <- BACISolo[,c("Year2", "BA", "CI", "Hours")]
  BACISoloPredict$Hours <- mean(BACISoloPredict$Hours)
  AdjustedCounts <- exp(predict(SoloModSlope, BACISoloPredict) + residuals(SoloModSlope))
  AdjustedCounts <- BACISolo$Count/BACISolo$Hours
  
  ModelSummaries <- lapply(list(SoloModMean, SoloModSlope), function(Model){
    Summary <- as.data.frame(summary(Model)[[11]][,c(1,2,4)])
    Summary$Estimate <- paste0(sprintf('%.2f',Summary$Estimate), " (", sprintf('%.2f',Summary$`Std. Error`), ")", ifelse(Summary$`Pr(>|z|)`<0.001, "***", ifelse(Summary$`Pr(>|z|)`<0.00, "**", ifelse(Summary$`Pr(>|z|)`<0.05, "*", ""))))
    #Summary$p <- paste0(ifelse(Summary$`Pr(>|z|)`<0.001, "<0.001", sprintf('%.3f',Summary$`Pr(>|z|)`)))
    Summary$Coefficient <- row.names(Summary)
    Summary <- Summary[,c(4, 1)]
    return(Summary)
  })
  ModelSummaries <- merge(ModelSummaries[[1]], ModelSummaries[[2]], by="Coefficient", all=T)
  ModelSummaries[is.na(ModelSummaries)] <- ""
  ModelSummaries <- ModelSummaries[c(1,5,2,4,3,6,8,7),]
  ModelSummaries$Coefficient <- c("Intercept", "Year", "A", "I", "A:I", "Year:A", "Year:I", "Year:A:I")
  names(ModelSummaries) <- c("Coefficient", "Mean", "Level/Slope")
  write.csv(ModelSummaries, paste0(ResultsFP, "TREESummaryModel.csv"), row.names=FALSE)
  
  BACISolo$SlopeFitted <- exp(predict(SoloModSlope, BACISoloPredict))
  BACISolo$MeanFitted <- exp(predict(SoloModMean, BACISoloPredict))
  BACISolo$CountsAdjusted <- AdjustedCounts
  BACISoloCont <- subset(BACISolo, CI==0)
  BACISoloInt <- subset(BACISolo, CI==1)
  
  BACICols <- c("Unprotected" = "#CC79A7", "Protected" = "#0072B2")
  Logged <- ggplot()+
    geom_point(aes(x=BACISoloCont$Year, y=BACISoloCont$CountsAdjusted, color="Unprotected"), alpha=0.5, shape=16)+
    geom_point(aes(x=BACISoloInt$Year, y=BACISoloInt$CountsAdjusted, color="Protected"), alpha=0.5, shape=16)+
    geom_line(aes(x=BACISoloCont$Year, y=BACISoloCont$SlopeFitted, group=BACISoloCont$BA, color="Unprotected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$SlopeFitted, group=BACISoloInt$BA, color="Protected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloCont$Year,y=BACISoloCont$MeanFitted, group=BACISoloCont$BA, color="Unprotected"),  lty="dashed", show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$MeanFitted, group=BACISoloInt$BA, color="Protected"), lty="dashed", show.legend = FALSE)+
    geom_vline(xintercept=unique(BACISolo$STATUS_YR)+0.5, color="black", linetype="dashed")+
    scale_colour_manual(breaks=names(BACICols), values=BACICols, name="Groups")+
    scale_x_continuous(breaks=c(min(BACISolo$Year):max(BACISolo$Year)))+
    scale_y_continuous(expand=c(0, 0), limits=c(1,round_any(max(BACISolo$CountsAdjusted),1000, f=ceiling)+10), trans = "log10")+ #)+
    xlab("Year")+
    ylab("Count")+
    theme(aspect.ratio=0.55, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          #panel.border = element_rect(size = 1, fill = NA),
          axis.text.x = element_text(size=plotfontsize, colour="black"), 
          axis.text.y = element_text(size=plotfontsize, colour="black"), 
          axis.ticks = element_line(colour="black"),
          axis.line = element_line(size = 0.5),
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top", legend.title = element_blank(),
          legend.key=element_blank())
  Logged
  #ggsave(paste0("/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Chapter3/TREE/Figures/FinalFigures/Box2FigI.pdf"), Logged, width = 125, height = 55, units = "mm") #Save as a pdf
  
  #g_legend<-function(a.gplot){
  #  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  #  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  #  legend <- tmp$grobs[[leg]]
  #  return(legend)}
  #mylegend<-g_legend(Normal)
  #grid.arrange(arrangeGrob(Logged + theme(legend.position="none"),
  #                         Normal + theme(legend.position="none"),
  #                         nrow=2),
  #             mylegend, nrow=1)
  #Logged
  library(gridExtra)
  plot_grid(BACIResultsPlot[[1]], Plot, ncol=1, align="hv")
  
  #Plots <- plot_grid(Normal + theme(legend.position="none"), Logged + theme(legend.position="none"), 
  #                   labels=c("(A)", "(B)"), ncol = 1, nrow = 2, label_x=0, label_y=1, label_fontface="plain")
  #Plots
  #legend <- get_legend(Normal + theme(legend.box.margin = margin(0, 0, 0, 0)))
  #PlotsLegend <- plot_grid(Plots, legend, rel_widths = c(0.7, 0.2), axis = ("t"), align = "h")
  ggsave(paste0(ResultsFP, "PlotCheckDELETE/", Solo, ".pdf"), Logged, width = 125, height = 55, units = "mm") #Save as a pdf
  
  DataCheck <- BACISolo[,c("Count", "CountEffort", "BA", "CI", "Year")]
  DataCheck <- DataCheck[with(DataCheck, order(CI, Year)),]
  return("Done")
})

#For IWC
Solo <- "Gavia immer.68"
Solo <- "Numenius arquata.299"
Solo <- "Mareca strepera.612"
Solo <- "Spatula clypeata.472"
pblapply(unique(BACIByPopulationMeanOrSlope2$SpecMatch), function(Solo){
  print(Solo)
  BACISolo <- subset(BACI, SpecMatch==Solo & CI==1)[,c("Count", "Hours", "Year", "BA", "CI", "STATUS_YR")]
  ok <- subset(try, SiteName=="Fort Smith-Moffett" & Year<=max(BACISolo$Year) & Year>=min(BACISolo$Year))[,c("Count", "Hours", "Year")]
  ok$CI <- 0
  ok$BA <- ifelse(ok$Year>unique(BACISolo$STATUS_YR), 1,0)
  ok$STATUS_YR <- unique(BACISolo$STATUS_YR)
  
  BACISolo <- rbind(BACISolo, ok, use.names=TRUE, fill=TRUE)
  BACISolo$Year2 <- BACISolo$Year-BACISolo$STATUS_YR
  SoloModSlope <- glm.nb(Count~Year2 + BA + BA*Year2 + CI + CI*BA + CI*Year2 + BA*CI*Year2, link=log, data=BACISolo)
  summary(SoloModSlope)
  SoloModMean <- glm.nb(Count~BA + CI + CI*BA + offset(log(Hours)), link=log, data=BACISolo)
  summary(SoloModMean)
  
  SoloModSlopeNull <- glm.nb(Count~BA + Year + offset(log(Hours)), link=log, data=BACISolo)
  if(AIC(SoloModSlopeNull)<AIC(SoloModSlope)){return(NULL)}
  
  ModelSummaries <- lapply(list(SoloModMean, SoloModSlope), function(Model){
    Summary <- as.data.frame(summary(Model)[[11]][,c(1,2,4)])
    Summary$Estimate <- paste0(sprintf('%.2f',Summary$Estimate), " (", sprintf('%.2f',Summary$`Std. Error`), ")", ifelse(Summary$`Pr(>|z|)`<0.05, "*", ""))
    #Summary$p <- paste0(ifelse(Summary$`Pr(>|z|)`<0.001, "<0.001", sprintf('%.3f',Summary$`Pr(>|z|)`)))
    Summary$Coefficient <- row.names(Summary)
    Summary <- Summary[,c(4, 1)]
    return(Summary)
  })
  ModelSummaries <- merge(ModelSummaries[[1]], ModelSummaries[[2]], by="Coefficient", all=T)
  ModelSummaries[is.na(ModelSummaries)] <- ""
  ModelSummaries <- ModelSummaries[c(1,5,2,4,3,6,8,7),]
  ModelSummaries$Coefficient <- c("Intercept", "Year", "A", "I", "A:I", "Year:A", "Year:I", "Year:A:I")
  names(ModelSummaries) <- c("Coefficient", "Mean", "Level/Slope")
  write.csv(ModelSummaries, paste0(ResultsFP, "TREESummaryModel.csv"), row.names=FALSE)
  BACISolo$SlopeFitted <- fitted(SoloModSlope)
  BACISolo$MeanFitted <- fitted(SoloModMean)
  
  BACISoloCont <- subset(BACISolo, CI==0)
  BACISoloInt <- subset(BACISolo, CI==1)
  
  BACICols <- c("Unprotected" = "#CC79A7", "Protected" = "#0072B2")
  Logged <- ggplot()+
    geom_point(aes(x=BACISoloCont$Year, y=BACISoloCont$Count, color="Unprotected"), alpha=0.5, shape=16)+
    geom_point(aes(x=BACISoloInt$Year, y=BACISoloInt$Count, color="Protected"), alpha=0.5, shape=16)+
    geom_line(aes(x=BACISoloCont$Year, y=BACISoloCont$SlopeFitted, group=BACISoloCont$BA, color="Unprotected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$SlopeFitted, group=BACISoloInt$BA, color="Protected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloCont$Year,y=BACISoloCont$MeanFitted, group=BACISoloCont$BA, color="Unprotected"),  lty="dashed", show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$MeanFitted, group=BACISoloInt$BA, color="Protected"), lty="dashed", show.legend = FALSE)+
    geom_vline(xintercept=unique(BACISolo$STATUS_YR)+0.5, color="black", linetype="dashed")+
    scale_colour_manual(breaks=names(BACICols), values=BACICols, name="Groups")+
    #scale_x_continuous(breaks=c(min(BACISolo$Year):max(BACISolo$Year)))+
    scale_y_continuous(expand=c(0, 0), limits=c(1,round_any(max(BACISolo$Count),1000, f=ceiling)+10), trans = "log10")+ #)+
    xlab("Year")+
    ylab("Count")+
    theme(aspect.ratio=0.55, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          #panel.border = element_rect(size = 1, fill = NA),
          axis.text.x = element_text(size=plotfontsize, colour="black"), 
          axis.text.y = element_text(size=plotfontsize, colour="black"), 
          axis.ticks = element_line(colour="black"),
          axis.line = element_line(size = 0.5),
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top", legend.title = element_blank(),
          legend.key=element_blank())
  Logged
  #ggsave(paste0("/Users/hannahwauchope/OneDrive - University Of Cambridge/PhD/Chapter3/TREE/Figures/FinalFigures/Box2FigI.pdf"), Logged, width = 125, height = 55, units = "mm") #Save as a pdf
  ggsave(paste0(ResultsFP, "PlotCheckDELETE/", Solo, ".pdf"), Logged, width = 125, height = 55, units = "mm") #Save as a pdf
  
  DataCheck <- BACISolo[,c("Count", "BA", "CI", "Year")]
  DataCheck <- DataCheck[with(DataCheck, order(CI, Year)),]
  return("Done")
})


BACISpec <- subset(BACI, SpecMatch==Solo)
Sites <- unique(BACISpec[,c("SiteCode", "Latitude", "Longitude", "CI", "WDPAID")])
PPDataSite <- subset(PPData, WDPAID=="555587005")

BACISites <- BACI[BACI$SiteCode %in% BACISpec$SiteCode,]
BACISitesCheck <- dcast(as.data.table(unique(BACISites[,c("SiteCode", "SpecMatch")])), SpecMatch~., length, value.var="SiteCode")

BACISummaries <- read.csv(file=paste0(MatchedFP, "ModelOutput/BACICategories.csv"))
BACIOutcomeCheck <- BACISummaries[BACISummaries$SpecMatch %in% BACISitesCheck$SpecMatch,]

BACISpec2 <- subset(BACI, SpecMatch==subset(BACISitesCheck, .==2)[2,]$SpecMatch)
Sites <- unique(BACISpec2[,c("SiteCode", "Latitude", "Longitude", "CI", "WDPAID")])

PPData <- read.dbf(paste0(DataFP, "ProtectedPlanet/WDPA_2020/WDPA_Jan2020-shapefile/WDPA_Jan2020-shapefile-polygons.dbf"))
PPDataSite <- subset(PPData, WDPAID=="555586922" | WDPAID=="555586979" | WDPAID=="555586982")
PPDataSite <- subset(PPData, GIS_AREA<74243 & GIS_AREA>74242)

BACISummaries[BACISummaries$SpecMatch %in% subset(BACI, WDPAID=="555512051")$SpecMatch,]

#### Model by Species DEAD ####
BACISpec <- subset(BACIModel, Species==Spec)
BACIBySpecies <- lapply(unique(BACIModel$Species), function(x){
  Mod <- glmer.nb(Count~Year + BA + BA*Year + CI + CI*BA + CI*Year + BA*CI*Year + (1|SpecMatch) + offset(log(Hours)), data=BACISpec)
  Mod2 <- glmer.nb(log(Count+1)~Year + BA + BA*Year + CI + CI*BA + CI*Year + (1|SpecMatch) + offset(log(Hours)), data=BACISpec)
  Area <- anova(Mod, Mod2) 
})

ModelsBySpeciesBA <- function(ProtectedData, Sign){
  Models <- pbmclapply(unique(ProtectedData$Species), function(Spec){
    if(file.exists(paste0(ResultsFP, "ModelResults/BA/", Spec, "_", Sign, "CHECK.csv"))){
      return(NULL)
    }
    write.csv(NULL, paste0(ResultsFP, "ModelResults/BA/", Spec, "_", Sign, "CHECK.csv"))
    print(Spec)
    SpeciesSubset <- subset(ProtectedData, Species==Spec)
    SpeciesSubset$SiteSpec <- NULL
    
    SiteThresh <- 4
    if(length(unique(SpeciesSubset$SiteCode))<SiteThresh){
      write.csv(NULL, paste0(ResultsFP, "ModelResults/BA/", Spec, "_", Sign, ".csv"))
      return(NULL)
    }
    SpeciesSubset <- CollVIF(SpeciesSubset, "BA")
    f <- as.formula(paste("Count", paste(c("Year", "BA","Year*BA", 
                                           names(SpeciesSubset)[1:(ncol(SpeciesSubset)-10)], "(1|SiteCode)", "offset(log(SpeciesSubset$Hours))"), 
                                         collapse = " + "), sep = " ~ "))
    f2 <- as.formula(paste("Count", paste(c("BA", 
                                            names(SpeciesSubset)[1:(ncol(SpeciesSubset)-10)], "(1|SiteCode)", "offset(log(SpeciesSubset$Hours))"), 
                                          collapse = " + "), sep = " ~ "))
    SpeciesSubsetModel <- glmer.nb(f, data=SpeciesSubset)
    if(isSingular(SpeciesSubsetModel)==TRUE){stop(paste0("Slope Model is singular for ", Spec))}
    SpeciesSubsetModelMean <- glmer.nb(f2, data=SpeciesSubset)
    if(isSingular(SpeciesSubsetModelMean)==TRUE){stop(paste0("Mean Model is singular for ", Spec))}
    
    Output <- as.data.frame(summary(SpeciesSubsetModel)[[10]])
    Output$Variable <- row.names(Output)
    Output$Model <- "Slope"
    
    Output2 <- as.data.frame(summary(SpeciesSubsetModelMean)[[10]])
    Output2$Variable <- row.names(Output2)
    Output2$Model <- "Mean"
    
    Output3 <- rbind(Output, Output2)
    
    Output3$Species <- Spec
    write.csv(Output3, paste0(ResultsFP, "ModelResults/BA/", Spec, "_", Sign, ".csv"), row.names = FALSE)
    return(Output3)
  }, mc.cores=ncores)
  return(Models)
}

#### CONFUSION ####
hey <- subset(BACISolo, CI==1)
SoloModSlope <- glm.nb(Count~Year2 + BA + BA*Year2 + offset(log(Hours)), link=log, data=BACISolo)
summary(SoloModSlope)

head(BACIMeanSlope)
#For CBC
Solo <- "AL00003.Anas acuta"
BACISolo <- subset(BA, SiteSpec==Solo)[,c("Count", "Year", "Hours", "Dataset", "BA", "STATUS_YR")]
BACISolo$Year2 <- BACISolo$Year-BACISolo$STATUS_YR
SoloModSlope <- glm.nb(Count~Year2 + BA + BA*Year2 + offset(log(Hours)), link=log, data=BACISolo)
summary(SoloModSlope)

BACISolo$SlopeFitted <- exp(predict(SoloModSlope, BACISolo))
BACISolo <- BACISolo[order(BACISolo$Year),]

Logged <- ggplot()+
  geom_point(aes(x=BACISolo$Year, y=BACISolo$Count), alpha=0.5, shape=16)+
  geom_line(aes(x=BACISolo$Year,y=BACISolo$SlopeFitted, group=BACISolo$BA), show.legend = FALSE)+
  geom_vline(xintercept=unique(BACISolo$STATUS_YR)+0.5, color="black", linetype="dashed")+
  scale_x_continuous(breaks=c(min(BACISolo$Year):max(BACISolo$Year)))+
  scale_y_continuous(expand=c(0, 0), limits=c(1,round_any(max(BACISolo$Count),1000, f=ceiling)+10), trans = "log10")+ #)+
  xlab("Year")+
  ylab("Count")+
  theme(aspect.ratio=0.55, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        #panel.border = element_rect(size = 1, fill = NA),
        axis.text.x = element_text(size=plotfontsize, colour="black"), 
        axis.text.y = element_text(size=plotfontsize, colour="black"), 
        axis.ticks = element_line(colour="black"),
        axis.line = element_line(size = 0.5),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top", legend.title = element_blank(),
        legend.key=element_blank())
Logged
#### Graveyard ####
#For loop match matrix
for(i in colnames(MatchMatrix)){
  print(i)
  StatYr <- unique(subset(ProtectedCountsSpecies, SiteCode==i)$STATUS_YR)
  
  ProtectedCountsSpeciesSite <- subset(ProtectedCountsSpecies, SiteCode==i)
  
  UnprotectedCountsSpeciesSite <- subset(UnprotectedCountsSpecies, Year<= (StatYr+TotalYearsBuffer) & Year>(StatYr-TotalYearsBuffer))
  UnprotectedCountsSpeciesSite$BA <- ifelse(UnprotectedCountsSpeciesSite$Year>StatYr, 1, 0)
  UnprotectedMeasuredYears <- dcast(as.data.table(UnprotectedCountsSpeciesSite), SiteCode~BA, length, value.var="Count")
  UnprotectedMeasuredYears <- UnprotectedMeasuredYears[UnprotectedMeasuredYears$'0'>=MeasuredYearsBuffer & UnprotectedMeasuredYears$'1'>=MeasuredYearsBuffer,]
  UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% UnprotectedMeasuredYears$SiteCode,]
  
  if(nrow(UnprotectedCountsSpeciesSite)==0){
    MatchMatrix[,i] <- 10
    next
  }
  
  #Conduct exact matching
  ProtCovs <- GetMeanCovs(ProtectedCountsSpeciesSite)
  UnProtCovs <- GetMeanCovs(UnprotectedCountsSpeciesSite)
  
  Anth <- as.character(ProtCovs$AnthRound) #Protected anthrome
  Region <- as.character(ProtCovs$GeoRegion) #Protected region
  ProtBeforeSlope <- unique(subset(ProtectedCountsSpecies, SiteCode==i)$BeforeSlope) #Protected before slope
  UnProtCovsSub <- subset(UnProtCovs, AnthRound==Anth & GeoRegion==Region) #Subset to correct antrhome and region
  UnprotectedCountsSpeciesSite <- UnprotectedCountsSpeciesSite[UnprotectedCountsSpeciesSite$SiteCode %in% UnProtCovsSub$SiteCode,] #Subset to correct antrhome and region
  if(nrow(UnprotectedCountsSpeciesSite)==0){
    MatchMatrix[,i] <- 10
    next
  }
  
  if(length(unique(UnprotectedCountsSpeciesSite$SiteCode))!=length(unique(UnprotectedCountsSpeciesSite$SiteSpec))){stop("Something's gone wrong, unprotected sitecode ! = sitespec")}
  
  #Calculate unprotected before slopes
  BeforeSlopes <- CalculateBeforeSlopes(subset(UnprotectedCountsSpeciesSite, BA==0), parallelise=FALSE)
  BeforeSlopes <- subset(BeforeSlopes, BeforeSlope==ProtBeforeSlope)
  if(nrow(BeforeSlopes)==0){
    MatchMatrix[,i] <- 10
    next
  }
  
  print(c(ProtBeforeSlope, unique(BeforeSlopes$BeforeSlope)))
  
  FilteredSites <- str_split_fixed(BeforeSlopes$SiteSpec, "[.]",3)[,1]
  
  ### Get the distances and add to matrix
  DistanceValues <- MDYearList[[paste0(StatYr)]]
  if(ncol(DistanceValues)==1){
    DistanceValues <- DistanceValues[rownames(DistanceValues) %in% FilteredSites,]
    if(length(DistanceValues)==1){
      k <- paste0(FilteredSites)
      MatchMatrix[k, i] <- DistanceValues
      MatchMatrix[!rownames(MatchMatrix) %in% FilteredSites,i] <- 10
    } else {
      MatchMatrix[c(names(DistanceValues)), i] <- DistanceValues
      MatchMatrix[!rownames(MatchMatrix) %in% FilteredSites,i] <- 10
    }
  } else {
    DistanceValues <- DistanceValues[rownames(DistanceValues) %in% FilteredSites,]
    if(is.null(nrow(DistanceValues))){
      k <- paste0(FilteredSites)
      MatchMatrix[k, i] <- DistanceValues[[i]]
      MatchMatrix[!rownames(MatchMatrix) %in% FilteredSites,i] <- 10
    } else {
      MatchMatrix[c(names(DistanceValues[,i])), i] <- DistanceValues[,i]
      MatchMatrix[!rownames(MatchMatrix) %in% FilteredSites,i] <- 10
    }
  }
  if(length(rownames(MatchMatrix)[rownames(MatchMatrix) %in% i])==1){
    MatchMatrix[i,i] <- 10
  }
}

library(stargazer)
data <- as.data.frame(cbind(a = rnorm(30), b = rnorm(30)))
fit_lm <- lm(data, formula = a ~ b)
stargazer(SoloModMean, SoloModSlope, style = "all2", column.labels = c("Mean", "Level/Slope"), colnames = NULL,
          single.row=TRUE, column.sep.width = "20pt", digits=2, model.numbers=FALSE, omit.table.layout = "n",
          intercept.bottom=FALSE, no.space=FALSE, align = TRUE, dep.var.caption = "",dep.var.labels.include = FALSE,
          covariate.labels = c("Intercept", "Year", "A", "I", "Year:A", "A:I", "Year:I", "Year:A:I"), type = "html", out = paste0(ResultsFP, "fit_lm.html"))

#ajps is good, so is io

### Old plot code
lapply(unique(BACIByPopulationMeanOrSlope$SpecMatch), function(Solo){
  print(Solo)
  BACISolo <- subset(BACI, SpecMatch==Solo)[,c("Count", "Hours", "BA", "CI", "Year", "STATUS_YR")]
  BACISolo$YearCent <- BACISolo$Year-BACISolo$STATUS_YR
  
  #AIC Comp
  SlopeNull <- glm.nb(Count~YearCent + CI + CI*YearCent + offset(log(Hours)), link=log, data=BACISolo)
  MeanNull <- glm.nb(Count~CI + offset(log(Hours)), link=log, data=BACISolo)
  
  Slope <- glm.nb(Count~YearCent + BA + BA*YearCent + CI + CI*BA + CI*YearCent + BA*CI*YearCent + offset(log(Hours)), link=log, data=BACISolo)
  Mean <- glm.nb(Count~BA + CI + CI*BA + offset(log(Hours)), link=log, data=BACISolo)
  
  if(AIC(SlopeNull)<AIC(Slope)){return("Null model performs better")}
  
  #Get summary data
  summary(Slope)
  summary(Mean)
  
  GetSummaryData <- TRUE
  if(GetSummaryData == TRUE){
    ModelSummaries <- lapply(list(SoloModMean, SoloModSlope), function(Model){
      Summary <- as.data.frame(summary(Model)[[11]][,c(1,2,4)])
      Summary$Estimate <- paste0(sprintf('%.2f',Summary$Estimate), " (", sprintf('%.2f',Summary$`Std. Error`), ")", ifelse(Summary$`Pr(>|z|)`<0.05, "*", ""))
      #Summary$p <- paste0(ifelse(Summary$`Pr(>|z|)`<0.001, "<0.001", sprintf('%.3f',Summary$`Pr(>|z|)`)))
      Summary$Coefficient <- row.names(Summary)
      Summary <- Summary[,c(4, 1)]
      return(Summary)
    })
    ModelSummaries <- merge(ModelSummaries[[1]], ModelSummaries[[2]], by="Coefficient", all=T)
    ModelSummaries[is.na(ModelSummaries)] <- ""
    ModelSummaries <- ModelSummaries[c(1,5,2,4,3,6,8,7),]
    ModelSummaries$Coefficient <- c("Intercept", "Year", "A", "I", "A:I", "Year:A", "Year:I", "Year:A:I")
    names(ModelSummaries) <- c("Coefficient", "Mean", "Level/Slope")
    write.csv(ModelSummaries, paste0(ResultsFP, "TREESummaryModel.csv"), row.names=FALSE)
  }
  
  #Adjust for effort
  BACISolo$CountEffort <- (BACISolo$Count+log(BACISolo$Hours))+(coef(Slope)[[1]])+1
  if(min(BACISolo$CountEffort)<0){return("Next")}
  
  SlopeEffort <- glm.nb(CountEffort~YearCent + BA + BA*YearCent + CI + CI*BA + CI*YearCent + BA*CI*YearCent, link=log, data=BACISolo)
  MeanEffort <- glm.nb(CountEffort~BA + CI + CI*BA, link=log, data=BACISolo)
  
  BACISolo$SlopeFitted <- fitted(SlopeEffort)
  BACISolo$MeanFitted <- fitted(MeanEffort)
  
  BACISoloCont <- subset(BACISolo, CI==0)
  BACISoloInt <- subset(BACISolo, CI==1)
  
  BACICols <- c("Unprotected" = "dodgerblue4", "Protected" = "orangered4")
  
  Normal <- ggplot()+
    geom_point(aes(x=BACISoloCont$Year, y=BACISoloCont$CountEffort, color="Unprotected"))+
    geom_point(aes(x=BACISoloInt$Year, y=BACISoloInt$CountEffort, color="Protected"))+
    geom_line(aes(x=BACISoloCont$Year, y=BACISoloCont$SlopeFitted, group=BACISoloCont$BA, color="Unprotected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$SlopeFitted, group=BACISoloInt$BA, color="Protected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloCont$Year,y=BACISoloCont$MeanFitted, group=BACISoloCont$BA, color="Unprotected"),  lty="dashed", show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$MeanFitted, group=BACISoloInt$BA, color="Protected"), lty="dashed", show.legend = FALSE)+
    geom_vline(xintercept=unique(BACISolo$STATUS_YR)+0.5, color="grey70")+
    scale_colour_manual(breaks=names(BACICols), values=BACICols, name="Groups")+
    scale_y_continuous(expand=c(0, 0), limits=c(0,round_any(max(BACISolo$Count),10, f=ceiling)+10))+ #trans = "log10")+
    xlab("Year")+
    ylab("Count")+
    theme(aspect.ratio=0.5, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          axis.text.x = element_text(size=plotfontsize, colour="black"), 
          axis.text.y = element_text(size=plotfontsize, colour="black"), 
          axis.ticks = element_line(colour="black"),
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top", legend.title = element_blank(),
          legend.key=element_blank())
  Normal
  Logged <- ggplot()+
    geom_point(aes(x=BACISoloCont$Year, y=BACISoloCont$CountEffort, color="Unprotected"))+
    geom_point(aes(x=BACISoloInt$Year, y=BACISoloInt$CountEffort, color="Protected"))+
    geom_line(aes(x=BACISoloCont$Year, y=BACISoloCont$SlopeFitted, group=BACISoloCont$BA, color="Unprotected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$SlopeFitted, group=BACISoloInt$BA, color="Protected"), show.legend = FALSE)+
    geom_line(aes(x=BACISoloCont$Year,y=BACISoloCont$MeanFitted, group=BACISoloCont$BA, color="Unprotected"),  lty="dashed", show.legend = FALSE)+
    geom_line(aes(x=BACISoloInt$Year,y=BACISoloInt$MeanFitted, group=BACISoloInt$BA, color="Protected"), lty="dashed", show.legend = FALSE)+
    geom_vline(xintercept=unique(BACISolo$STATUS_YR)+0.5, color="grey70")+
    scale_colour_manual(breaks=names(BACICols), values=BACICols, name="Groups")+
    scale_y_continuous(expand=c(0, 0), limits=c(1,round_any(max(BACISolo$CountEffort),100, f=ceiling)), trans = "log10")+ #)+
    xlab("Year")+
    ylab("Count")+
    theme(aspect.ratio=0.5, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          axis.text.x = element_text(size=plotfontsize, colour="black"), 
          axis.text.y = element_text(size=plotfontsize, colour="black"), 
          axis.ticks = element_line(colour="black"),
          legend.text = element_text(size=plotfontsize, colour="black"),
          legend.justification = "top", legend.title = element_blank(),
          legend.key=element_blank())
  Logged
  
  BothPlots <- FALSE
  if(BothPlots==TRUE){
    library(cowplot)
    Plots <- plot_grid(Normal + theme(legend.position="none"), Logged + theme(legend.position="none"), 
                       labels=c("(A)", "(B)"), ncol = 1, nrow = 2, label_x=0, label_y=1, label_fontface="plain")
    legend <- get_legend(Normal + theme(legend.box.margin = margin(0, 0, 0, 0)))
    Logged <- plot_grid(Plots, legend, rel_widths = c(0.7, 0.2), axis = ("t"), align = "h")
  }
  ggsave(paste0(ResultsFP, "PlotCheckDELETE/", Solo, ".pdf"), Logged, width = 200, height = 150, units = "mm") #Save as a pdf
  
  return("Done")
})

BABySpecies <- rbindlist(lapply(c(1,0,-1), function(BS){
  dir.create(file.path(paste0(paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/"))), showWarnings = FALSE)
  Species <- rbindlist(pbmclapply(unique(BAModel$Species), function(x){
    if(file.exists(paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"))){
      return(NULL)
    }
    write.csv(NULL, paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
    
    ModelOutput <- data.frame(matrix(ncol = 8, nrow = 1))
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- BS
    names(ModelOutput) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model", "Species", "BeforeSlope")
    
    BASpec <- subset(BAModel, Species==x & BeforeSlope == BS)
    if(nrow(BASpec)==0){
      ModelOutput[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model")] <- c("NoRows")
      write.csv(ModelOutput, paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    if(length(unique(BASpec$SiteSpec))<4){
      ModelOutput[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model")] <- c("Lessthan4specmatch")
      write.csv(ModelOutput, paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    BASpec$YearScale <- BASpec$Year-BASpec$STATUS_YR
    SlopeModel <- tryCatch(glmer.nb(Count~YearScale + BA + BA*YearScale + (1|SiteCode) + offset(log(Hours)), data=BASpec), error=function(e){NULL})
    MeanModel <- tryCatch(glmer.nb(Count~BA + (1|SiteCode) + offset(log(Hours)), data=BASpec), error=function(e){NULL})
    if(is.null(SlopeModel) & is.null(MeanModel)){
      ModelOutput[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model")] <- c("NullModels")
      write.csv(ModelOutput, paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    
    if(!is.null(SlopeModel)){
      SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff)
      SlopeOutput$Coef <- row.names(SlopeOutput)
      SlopeOutput$Model <- "Slope"
      if(!is.null(MeanModel)){
        MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
        MeanOutput$Coef <- row.names(MeanOutput)
        MeanOutput$Model <- "Mean"
        ModelOutput <- rbind(SlopeOutput, MeanOutput)
      } else {
        ModelOutput <- SlopeOutput
      }
    } else {
      MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
      MeanOutput$Coef <- row.names(MeanOutput)
      MeanOutput$Model <- "Mean"
      ModelOutput <- MeanOutput
    }
    
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- unique(BASpec$BeforeSlope)
    names(ModelOutput) <- c("Estimate", "StError", "Z", "P", "Coef", "Model", "Species", "BeforeSlope")
    
    write.csv(ModelOutput, paste0(ResultsFP, "SpeciesModels/BABySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
    return(ModelOutput)
  }, mc.cores=ncores))
  return(Species)
}))
save(BABySpecies, file=paste0(ResultsFP, "SpeciesModels/BABySpecies.RData"))

BACIBySpecies <- rbindlist(lapply(c(1,0,-1), function(BS){
  dir.create(file.path(paste0(paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/"))), showWarnings = FALSE)
  Species <- rbindlist(pbmclapply(unique(BACIModel$Species), function(x){
    if(file.exists(paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"))){
      return(NULL)
    }
    write.csv(NULL, paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"))
    
    ModelOutput <- data.frame(matrix(ncol = 8, nrow = 1))
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- BS
    names(ModelOutput) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model", "Species", "BeforeSlope")
    
    BACISpec <- subset(BACIModel, Species==x & BeforeSlope == BS)
    if(nrow(BACISpec)==0){
      ModelOutput[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model")] <- c("NoRows")
      write.csv(ModelOutput, file=paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    if(length(unique(BACISpec$SpecMatch))<4){
      ModelOutput[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model")] <- c("Lessthan4specmatch")
      write.csv(ModelOutput, file=paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    BACISpec$YearScale <- BACISpec$Year-BACISpec$STATUS_YR
    SlopeModel <- tryCatch(glmer.nb(Count~YearScale + BA + BA*YearScale + CI + CI*BA + CI*YearScale + BA*CI*YearScale + (1|SiteCode) + offset(log(Hours)), data=BACISpec), error=function(e){NULL})
    MeanModel <- tryCatch(glmer.nb(Count~BA + CI + CI*BA + (1|SiteCode) + offset(log(Hours)), data=BACISpec), error=function(e){NULL})
    if(is.null(SlopeModel) & is.null(MeanModel)){
      ModelOutput[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)", "Coef", "Model")] <- c("NullModels")
      write.csv(ModelOutput, file=paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
      return(ModelOutput)
    }
    
    if(!is.null(SlopeModel)){
      SlopeOutput <- as.data.frame(summary(SlopeModel)$coeff)
      SlopeOutput$Coef <- row.names(SlopeOutput)
      SlopeOutput$Model <- "Slope"
      if(!is.null(MeanModel)){
        MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
        MeanOutput$Coef <- row.names(MeanOutput)
        MeanOutput$Model <- "Mean"
        ModelOutput <- rbind(SlopeOutput, MeanOutput)
      } else {
        ModelOutput <- SlopeOutput
      }
    } else {
      MeanOutput <- as.data.frame(summary(MeanModel)$coeff)
      MeanOutput$Coef <- row.names(MeanOutput)
      MeanOutput$Model <- "Mean"
      ModelOutput <- MeanOutput
    }
    
    ModelOutput$Species <- x
    ModelOutput$BeforeSlope <- unique(BACISpec$BeforeSlope)
    names(Species) <- c("Estimate", "StError", "Z", "P", "Coef", "Model", "Species", "BeforeSlope")
    write.csv(ModelOutput, file=paste0(ResultsFP, "SpeciesModels/BACIBySpecies_",Database, "_", ZeroThreshFP, "ZeroThresh_", TotalYearsBuffer, "Buffer_", ImputeFlag, "/", x, "_", BS, ".csv"), row.names=FALSE)
    return(ModelOutput)
  }, mc.cores=ncores))
  return(Species)
}))
save(BACIBySpecies, file=paste0(MatchedFP, "SpeciesModels/BACIBySpecies.RData"))

### Remove colonisations and extinctions
#BACISummariseForPlot <- subset(BACISummariseForPlot, variable!="ExtinctIn" & variable!="ExtinctBoth" & variable!="ColonisedOut")
#FactorLevels <- c("ColonisedIn", "ExtinctOut", "Success", "ColonisedBoth", "NoChangePosStable", "NoChangeNeg", "Failure")
#ColourListPlot <- c(BlueCols(5)[1:3],"azure3", "azure2", "#fffdbfff", RedCols(5)[5])

#Now BACI after slope plot
AfterSlopes2 <- CalculateBeforeSlopes(subset(BACI, BA==1), ZeroThresh, parallelise = TRUE, keepinsig = TRUE)
AfterSlopes2 <- merge(AfterSlopes2, unique(BACI[,c("SiteSpec", "SpecMatch", "CI")]))
names(AfterSlopes2) <- c("SiteSpec", "AfterSlope", "SpecMatch", "CI")
BACISummarise2 <- merge(BACISummarise, subset(AfterSlopes2, CI==1)[,c("SpecMatch", "AfterSlope")])
BACISummariseForPlot2 <- melt(dcast(as.data.table(subset(BACISummarise2, OutcomeSummary!="Excluded")), Species~AfterSlope, length, value.var="SpecMatch"), id.vars="Species")
BACISummariseForPlot2$variable <- factor(BACISummariseForPlot2$variable)
BACISummariseForPlot2$variable <- recode_factor(BACISummariseForPlot2$variable, "-1" = "Declining", "0" = "Stable", "1" = "Increasing", "AllZero" = "NotPresent")
BACISummariseForPlot2$variable <- factor(BACISummariseForPlot2$variable, levels=c("Increasing", "Stable", "Declining", "NotPresent"))

ColourListPlot <- c(BlueCols(5)[1],"grey50", RedCols(5)[5], "grey10")
FactorLevels <- c("Increasing", "Stable", "Declining", "NotPresent")
BACISummariseForPlot2 <- subset(BACISummariseForPlot2, variable!="NotPresent")
FactorLevels <- c("Increasing", "Stable", "Declining")
BACIResultsPlotAfter <- BarPlotWidths(BACISummariseForPlot2, Log=TRUE, Cov="Species", FactorLevels = FactorLevels, xlabel="Species", AngledxLabs=FALSE, ColourList=ColourListPlot, SuccessValues=c("Increasing"), MainResults="MainWithSpecies", aspectrat=0.5)
plot_grid(BACIResultsPlot[[1]], BACIResultsPlotAfter[[1]], ncol=1, align="hv")

WidthData <- BACIResultsPlot[[2]]

BACISummariseForPlot2 <- merge(BACISummariseForPlot2, unique(WidthData[,c("Species", "LogWidth", "Width", "Position")]), by="Species")
BACISummariseForPlot2Percent <- dcast(BACISummariseForPlot2, Species~., sum, value.var="value")
names(BACISummariseForPlot2Percent) <- c("Species", "Total")
BACISummariseForPlot2 <- merge(BACISummariseForPlot2, BACISummariseForPlot2Percent)
BACISummariseForPlot2$PercentVal <- BACISummariseForPlot2$value/BACISummariseForPlot2$Total
ColourList2 <- c(BlueCols(5)[1],"grey50", RedCols(5)[5], "grey10")
BACISummariseForPlot2Unique <- unique(BACISummariseForPlot2[,c("Species", "LogWidth", "Position")])

Plot <- ggplot(BACISummariseForPlot2, aes(x=Position, y=PercentVal, fill=variable))+
  geom_bar(stat="identity", position="stack", width=BACISummariseForPlot2$LogWidth) +
  scale_fill_manual(values = ColourList2, name=("Category"))+ #, "grey60", rev(brewer.pal(11, "RdYlBu")[7:11]))
  scale_x_continuous(expand = c(0, 0), breaks=BACISummariseForPlot2Unique$Position, labels=BACISummariseForPlot2Unique$Species) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Proportion")+
  xlab(xlabel)+
  theme(aspect.ratio=aspectrat, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=plotfontsize),
        legend.title=element_text(size=plotfontsize),
        axis.text.x = element_text(size=6, colour="black", angle=90),
        axis.ticks = element_line(colour="black"),
        legend.text = element_text(size=plotfontsize, colour="black"),
        legend.justification = "top")
Plot
plot_grid(BACIResultsPlot[[1]], Plot, ncol=1, align="hv")

### Covariate plots
BAResultsCat$OutcomeSummary <- factor(BAResultsCat$OutcomeSummary, levels=c("Failure", "Neutral", "Success"))
BAResultsCat$OutcomeSummary <- recode_factor(BAResultsCat$OutcomeSummary, "Failure" = "1", "Neutral" = "2", "Success" = "3")
BAResultsCat$OutcomeSummary <- as.numeric(as.character(BAResultsCat$OutcomeSummary))
BAResultsCat <- BAResultsCat[complete.cases(BAResultsCat),]
AllEffects2 <- brm(OutcomeSummary~MigStatus + RedList + GenLength + Family + GovMean + AnthMode + PAAreaLog + Ramsar + (1|SiteCode) + (1|Species), BAResultsCat, family='cumulative')

sjPlot::plot_model(AllEffects, type="re", terms="PAAreaLog")[[1]]
