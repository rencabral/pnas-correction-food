#Escapement generator_PNAS correction
#This consider's Hilborn and Ovando et al.'s suggestions

gc()
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)

#library(ramlegacy)

#LOAD THE MOST LATEST RAM Data v.4.491
#I've downloaded this file from the RAM website and saved in my local directory
load("/Users/ren/Documents/GitHub/pnas-correction-food/DBdata[asmt][v4.491].RData")
RAMDATA2 <- timeseries_values_views
head(RAMDATA2)
colnames(RAMDATA2)
RAMDATA3 <- RAMDATA2 %>% dplyr::select(stockid,year,F,FdivFmsy)
head(RAMDATA3)

#Filter all the data then merge everything! I will use megadata as merger file.
F_RAM<-RAMDATA3 %>% dplyr::select(stockid,year,F) %>% filter(!is.na(F)) %>% group_by(stockid) %>% slice(which.max(year)) %>% dplyr::select(-year)
FdivFmsy_RAM <- RAMDATA3 %>% dplyr::select(stockid,year,FdivFmsy) %>% filter(!is.na(FdivFmsy)) %>% group_by(stockid) %>% slice(which.max(year)) %>% dplyr::select(-year)

#load Fmsy from RAM
RAM_bioparams <- bioparams_values_views %>% dplyr::select(stockid,Fmsy) 
head(RAM_bioparams)

#load MegaData
MegaData <- readRDS(file = "/Users/ren/Documents/GitHub/pnas-correction-food/MegaData.rds")
head(MegaData)
Mng2_MegaData <- MegaData %>% dplyr::select("SpeciesID","Manage","stockid","r_fin","Efin_BAU1","MSYfin") %>% filter(Manage==1)
head(Mng2_MegaData)

MatchedFile <- left_join(Mng2_MegaData, RAM_bioparams, by='stockid') %>%
  left_join(., F_RAM, by='stockid') %>% left_join(., FdivFmsy_RAM, by='stockid')
head(MatchedFile)
dim(MatchedFile)
sortedMatchedFile <- MatchedFile[order(- MatchedFile$MSYfin),]

sum(sortedMatchedFile$FdivFmsy>=0,na.rm=T)

#FILL GAP: GET COSTELLO et al. (2016) F/Fmsy for the assessed stocks
CostelloData <- read.csv("/Users/ren/Documents/CODES/FoodProvision/Aquamaps/ProjectionData.csv", stringsAsFactors = FALSE)
Costello2012 <- CostelloData %>% filter(Year=="2012")

##Check Trachurus trachurus
Trachurus <- Costello2012 %>% filter(SciName == "Trachurus trachurus")
sum(Trachurus$MSY)
Trachurus_MegaData <- MegaData %>% filter(SciName == "Trachurus trachurus")
sum(Trachurus_MegaData$MSYfin)

Costello2012_ram <- Costello2012 %>% filter(Dbase=="RAM") %>% dplyr::select(IdOrig,FvFmsy)
head(Costello2012_ram) #ok, match this with our species names PNAS

CostelloFvFmsy<-Costello2012_ram %>% separate(IdOrig, c("A", "stockid_temp","C","D","E")) 
CostelloFvFmsy$stockid <- CostelloFvFmsy$stockid_temp

#manually checked - recovering stockid
for (i in c(44,45,46,62:81,83:85,87,106:113,143,144,184:191)){
  CostelloFvFmsy$stockid[i] <- CostelloFvFmsy$C[i]
}

for (i in c(228,229,362,371,378,391,392)){
  CostelloFvFmsy$stockid[i] <- paste0(CostelloFvFmsy$stockid_temp[i],"-",CostelloFvFmsy$C[i])
}

#384 is a special case
CostelloFvFmsy$stockid[384] <- paste0(CostelloFvFmsy$stockid_temp[384],"-",CostelloFvFmsy$C[384],"-",CostelloFvFmsy$D[384])
CostelloFvFmsy <- CostelloFvFmsy %>% dplyr::select(stockid,FvFmsy)

names(CostelloFvFmsy)[names(CostelloFvFmsy) == "FvFmsy"] <- "FvFmsy_Costello"
head(CostelloFvFmsy)

##MERGE COSTELLO's F/Fmsy with our main database
sortedMatchedFile_V2 <- left_join(sortedMatchedFile,CostelloFvFmsy,by="stockid")
sum(sortedMatchedFile_V2$FvFmsy_Costello>=0,na.rm=T) #315 entries

#0k, now combine F and Fmsy RAM values
sortedMatchedFile_V2$FcombFmsy <- sortedMatchedFile_V2$F/sortedMatchedFile_V2$Fmsy

#now, FvFmsy_Fin
sortedMatchedFile_V2$FvFmsy_Ray <- sortedMatchedFile_V2$FdivFmsy

sum(is.na(sortedMatchedFile_V2$FvFmsy_Ray)==F) #210

for (i in 1:dim(sortedMatchedFile_V2)[1]){
  if (is.na(sortedMatchedFile_V2$FvFmsy_Ray[i])==T){
    sortedMatchedFile_V2$FvFmsy_Ray[i] <- sortedMatchedFile_V2$FcombFmsy[i]
  }
}

sum(is.na(sortedMatchedFile_V2$FvFmsy_Ray)==F)

#now, add the Costello et al. F/Fmsy
for (i in 1:dim(sortedMatchedFile_V2)[1]){
  if (is.na(sortedMatchedFile_V2$FvFmsy_Ray[i])==T){
    sortedMatchedFile_V2$FvFmsy_Ray[i] <- sortedMatchedFile_V2$FvFmsy_Costello[i]
  }
}

sum(is.na(sortedMatchedFile_V2$FvFmsy_Ray)==F)
sum(sortedMatchedFile_V2$FvFmsy_Ray>=0,na.rm=T) #ok, 340 stocks.

#let us add the new data to MagaData. Make others ER=0. MPAs will do nothing to the stocks that we do not have data on.
head(MegaData)
FinalFvFmsy<-sortedMatchedFile_V2 %>% dplyr::select(stockid,FvFmsy_Ray)
MegaData_Ray <- left_join(MegaData,FinalFvFmsy,by="stockid")

#indicate that stocks included in our analysis
MegaData_Ray$INCLUDE <- 1
MegaData_Ray$INCLUDE[is.na(MegaData_Ray$FvFmsy_Ray) & MegaData_Ray$Manage==1] <- 0

##Update Efin_BAU1
MegaData_Ray <- MegaData_Ray %>% mutate(Efin_BAU1_Ray = 1-(r_fin*FvFmsy_Ray/2))

#if negative, make it zero
MegaData_Ray$Efin_BAU1_Ray[MegaData_Ray$Efin_BAU1_Ray<0] <- 0

for (i in 1:dim(MegaData_Ray)[1]){
  if (MegaData_Ray$Manage[i]==0){
    MegaData_Ray$Efin_BAU1_Ray[i] <- MegaData_Ray$Efin_BAU1[i]
  }
}

#for no RAM/Costello et al. data (RAM stocks), make escapement = 1 (meaning, they are unfished)
for (i in 1:dim(MegaData_Ray)[1]){
  if (MegaData_Ray$Manage[i]==1 & is.na(MegaData_Ray$Efin_BAU1_Ray[i])==T){
    MegaData_Ray$Efin_BAU1_Ray[i]<-1
  }
}

#make Escapement of unassessed Trachurus trachurus = 1 (meaning, no fishing for that species)
for (i in 1:dim(MegaData_Ray)[1]){
  if (MegaData_Ray$SciName[i]=="Trachurus trachurus" & MegaData_Ray$Manage[i]==0){
    MegaData_Ray$Efin_BAU1_Ray[i]<-1
    MegaData_Ray$INCLUDE[i] <- 0
  }
}

#for the remaining trachurus trachurus, use Costello MSY of 50,352 MMT.
for (i in 1:dim(MegaData_Ray)[1]){
  if (MegaData_Ray$stockid[i]=="HMACKWA"){
    MegaData_Ray$MSYfin[i]<-50352
    MegaData_Ray$Kfin[i]<-4*50352/MegaData_Ray$r_fin[i]
  }
}

#if Exploitation rate is > r, make Exploitation rate == r. At Exploitation rate = r, biomass and catch will be zero 
#if catch is < 0, we set the catch to zero anyway in our code
MegaData_Ray <- MegaData_Ray %>% mutate(Efin_BAU1_Ray= (Efin_BAU1_Ray*((1-Efin_BAU1_Ray)<=r_fin)) + ((1-r_fin)*((1-Efin_BAU1_Ray)>r_fin)))

#Ensure that min=0 and max=1
min(MegaData_Ray$Efin_BAU1_Ray)
max(MegaData_Ray$Efin_BAU1_Ray)

MegaData_Ray$ExploitationRate_BAU1_Ray<-1-MegaData_Ray$Efin_BAU1_Ray

#how many stocks?
sum(MegaData_Ray$INCLUDE)

MegaData_Ray <- MegaData_Ray %>% select(SpeciesID,Manage,stockid,SciName,r_fin,m_fin,MSYfin,Kfin,FvFmsy_Ray,INCLUDE,Efin_BAU1_Ray,ExploitationRate_BAU1_Ray)
head(MegaData_Ray)
saveRDS(MegaData_Ray,file="/Users/ren/Documents/GitHub/pnas-correction-food/MegaData_Ray.rds")