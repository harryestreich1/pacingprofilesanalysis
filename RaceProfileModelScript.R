### Load libraries
library(fda)
library(tidyverse)
library(depmixS4)
library(ggpubr)

source("RaceProfileModelFormulas.R")
### Load data
menData = read.csv("MK11000m.csv", encoding="UTF-8", stringsAsFactors=FALSE,check.names=FALSE)
womenData = read.csv("WK1500m.csv", encoding="UTF-8", stringsAsFactors=FALSE,check.names=FALSE)

## Variables (can be changed)
allEvents = FALSE ## set FALSE to only do Domestic Finals and All International Races
minimumRaces = 10 ## set to minimum required races to make plot

# Plot Variables
runCareerStates = TRUE ## set FALSE to ignore career state plot
runPC = TRUE ## set FALSE to ignore PC tracking plot

# Variables (do not change)
convergenceIterations = 200 # default 200, reducing with make model less accurate but will run quicker
menBasis = 20 ## do not change, default 20
womenBasis = 10 ## do not change, default 10
minDate = min(min(menData$Year), min(womenData$Year))
maxDate = max(max(menData$Year), max(womenData$Year))




## Create PCA Objects
menPCA = dataToPCA(menData, menBasis)
womenPCA = dataToPCA(womenData, womenBasis)

## Create HMM Objects
menDataT = dataTransform(menData, menPCA, allEvents, menBasis)
menHMM = PCAtoHMM(menDataT, convergenceIterations)
womenDataT = dataTransform(womenData, womenPCA, allEvents, womenBasis)
womenHMM = PCAtoHMM(womenDataT, convergenceIterations)

## Create Career State Summary Plots
if(runCareerStates == TRUE)
{
  exportCareerStatePlot(menDataT, menHMM, minDate, maxDate, minimumRaces, "Men")
  exportCareerStatePlot(womenDataT, womenHMM, minDate, maxDate, minimumRaces, "Women")
}

## Create PC Summary Plots
if(runPC == TRUE)
{
  exportPCPlot(menDataT, menHMM, minDate, maxDate, minimumRaces, "Men")
  exportPCPlot(womenDataT, womenHMM, minDate, maxDate, minimumRaces, "Women")
}


