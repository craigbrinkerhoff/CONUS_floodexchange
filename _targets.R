# _targets.R file

#necessary global packages for targets and pipelining. ALl other packages are hardcoded
library(targets)
library(tarchetypes)
library(dplyr)
library(tidyr)
library(crew)
library(ggplot2)
library(tibble)

source("src/functions_floodplain.R")
source("src/validation.R")
source('src/utils.R')
source('src/figures.R')

## usgs dem vertical accuracy: https://www.mdpi.com/2072-4292/14/4/940

## PARAMETERS
gageRecordStart <- '1983-01-01'
gageRecordEnd <- '2023-12-31'
minRecordLength <- 10 #[yrs] minimum number of years on record for a gage to be included

minADCPMeas <- 20 #min depth stage adcp measurements for AHG
minAHGr2 <- 0.30 #min depth AHG fit

usgs_maps <- sf::st_read('data/path_to_data/CONUS_connectivity_data/USGS_models/USGSmodels_area.shp')

#distributed computing across cpu-node
# tar_option_set(
#   controller = crew_controller_local(workers = 45)
# )

## MAP PIPELINE OVER STATIC BRANCHES OF GAGES
gageAnalysis <- tar_map(
  values = tibble( # Setup static branching
    method_function = rlang::syms(c("getBasinGages")),
    huc4 = c('0108', '0109') #0110
  ),
  names='huc4',

  ## PREP GAGES
  tar_target(gagesBasin, method_function(huc4, gageRecordStart, gageRecordEnd)),
  tar_target(gage, prepGage(gagesBasin),
            pattern=map(gagesBasin),
            iteration='list'),
  tar_target(gageRecord, prepFlowRecord(gage, gageRecordStart, gageRecordEnd, minRecordLength),
            pattern=map(gage),
            iteration='list'),
  tar_target(gages_fixed, unlistGages(gage)),

  ## PREP BASIN DATA
  tar_target(basinData, buildBasinDataPackage(huc4)),

  ## BUILD DEPTH AHG RELATION
  tar_target(depthAHG, buildDepthAHG(gage, minADCPMeas),
            pattern=map(gage),
            iteration='list'),

  ## UPSCALE TO RIVER NETWORK (12/20/24: to do: add resources arguments so that x number of cores are passed to each basin Analysis)
  tar_target(upscalingModel, buildUpscalingModel(huc4, gageRecord, gage, BHGmodel, depthAHG, minAHGr2, gageRecordStart, gageRecordEnd)),
  tar_target(basinModel, buildNetworkModel(huc4, upscalingModel, BHGmodel, basinData), deployment='main'),
  tar_target(basinAnalysis, runNetworkModel(huc4, basinData, basinModel), deployment='main'),

  ## HORTON SCALING TO CAPTURE INUNDATION IN STREAMS < 10M WIDE
  tar_target(hortonResults, hortonScaling(basinAnalysis, huc4)),

  ## VALIDATE AT REACHES WITH USGS INUNDATION MODELS
  tar_target(val_FEMA, valModelFEMA(huc4, preppedFEMA, basinAnalysis, basinData)),
  tar_target(val_USGS, valModelUSGS(huc4, basinAnalysis, usgs_maps)),

  ## BASIN SPECIFIC FIGURES
  tar_target(upscalingPlot, upscalingFig(upscalingModel, huc4))
)


list(
  ## PREP BANKFULL MODELs AND DATA
  tar_target(BHGdata, dataBHG()),
  tar_target(BHGmodel, modelsBHG()),

  ## PREP VALIDATION MAPS
  tar_target(preppedFEMA, prepFEMA()),
  
  ## RUN GAGE ANALYSIS (STATIC BRANCHED)
  gageAnalysis,

  ## COMBINE
  tar_combine(gage_combined, gageAnalysis$gages_fixed, command = dplyr::bind_rows(!!!.x)),
  tar_combine(val_FEMA_combined, gageAnalysis$val_FEMA, command = dplyr::bind_rows(!!!.x)),
  tar_combine(val_USGS_combined, gageAnalysis$val_USGS, command = dplyr::bind_rows(!!!.x)),
  tar_combine(hortonResults_combined, gageAnalysis$hortonResults, command = dplyr::bind_rows(!!!.x)),

  ## FIGURES
  tar_target(valFEMAFig, makeValFEMA(val_FEMA_combined, gage_combined, BHGdata)),
  tar_target(valUSGSFig, makeValUSGS(val_USGS_combined, BHGdata))
)