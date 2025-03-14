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

snapping_thresh <- 5000 #[m]
area_thresh_perc <- 0.05

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

  ## UPSCALE TO RIVER NETWORK (2/18/25: to do: add resources arguments so that x number of cores are passed to each basin Analysis and regulateFlooding, which runs in parallel too)
  tar_target(upscalingModel, buildUpscalingModel(huc4, gageRecord, gage, BHGmodel, depthAHG, minAHGr2, gageRecordStart, gageRecordEnd)),
  tar_target(basinModel, buildNetworkModel(huc4, upscalingModel, BHGmodel, basinData), deployment='main'),
  tar_target(basinAnalysis, runNetworkModel(huc4, basinData, basinModel), deployment='main'),

  ## CALCULATE UPSTREAM FLOW REGULATION BY LAKES/RESERVOIRS
  tar_target(basinAnalysis_fin, regulateFlooding(basinAnalysis, GWD, area_thresh_perc, huc4, snapping_thresh)),

  ## HORTON SCALING TO CAPTURE INUNDATION IN STREAMS < 10M WIDE
  tar_target(hortonResults, hortonScaling(basinAnalysis, huc4, 'barrier')),
 # tar_target(hortonResults_nobarrier, hortonScaling(basinAnalysis_barriers, huc4, 'no_barrier')),

  ## VALIDATE AT REACHES WITH USGS INUNDATION MODELS
  tar_target(val_FEMA, valModelFEMA(huc4, preppedFEMA, basinAnalysis, basinData)),
  tar_target(val_USGS, valModelUSGS(huc4, basinAnalysis, usgs_maps)),
  tar_target(val_USGSvols, valModelUSGSvols(huc4, basinAnalysis, val_USGS, c('q0_2', 'q0_5', 'q1', 'q2', 'q4', 'q10', 'q20','q50', 'q80', 'q90', 'q96', 'q98', 'q99', 'q99_5', 'q99_8'))),

  ## BASIN SPECIFIC FIGURES
  tar_target(upscalingPlot, upscalingFig(upscalingModel, huc4))
)


list(
  ## PREP BANKFULL MODELs AND DATA
  tar_target(BHGdata, dataBHG()),
  tar_target(BHGmodel, modelsBHG()),

  ## PREP FLOW BARRIER DATASET
  tar_target(GWD, prepGWD()),

  ## PREP VALIDATION MAPS
  tar_target(preppedFEMA, prepFEMA()),
  
  ## RUN GAGE ANALYSIS (STATIC BRANCHED)
  gageAnalysis,

  ## COMBINE
  tar_combine(gage_combined, gageAnalysis$gages_fixed, command = dplyr::bind_rows(!!!.x)),
  tar_combine(val_FEMA_combined, gageAnalysis$val_FEMA, command = dplyr::bind_rows(!!!.x)),
  tar_combine(val_USGS_combined, gageAnalysis$val_USGS, command = dplyr::bind_rows(!!!.x)),
  tar_combine(val_USGSvols_combined, gageAnalysis$val_USGSvols, command = dplyr::bind_rows(!!!.x)),
  
  tar_combine(model_combined, gageAnalysis$basinAnalysis_fin, command = dplyr::bind_rows(!!!.x)),
  
  tar_combine(hortonResults_combined, gageAnalysis$hortonResults, command = dplyr::bind_rows(!!!.x)),
  #tar_combine(hortonResults_combined_nobarrier, gageAnalysis$hortonResults_nobarrier, command = dplyr::bind_rows(!!!.x)),

  ## FIGURES
  tar_target(valFEMAFig, makeValFEMA(val_FEMA_combined, gage_combined, BHGdata)),
  tar_target(valUSGSareaFig, makeValUSGSarea(val_USGS_combined)),
  tar_target(valUSGSvolFig, makeValUSGSvol(val_USGSvols_combined)),
  tar_target(streamOrderFig, makeHortonScalingFig(hortonResults_combined, model_combined)),
  tar_target(regulationFig, makeRegulationFig(model_combined)),
  tar_target(mapFig, makeMapFig(usgs_maps))
)