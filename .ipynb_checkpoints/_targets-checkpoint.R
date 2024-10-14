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
source('src/utils.R')

## usgs dem vertical accuracy: https://www.mdpi.com/2072-4292/14/4/940

benchmark_index <- 3 #just a dummy parameter to force the pipeline to rerun entirely for benchmarking
print(paste0('benchmark test ', benchmark_index))

## PARAMETERS
path_to_data <- '/nas/cee-water/cjgleason/craig/CONUS_ephemeral_data'
gageRecordStart <- '2013-01-01'
gageRecordEnd <- '2023-12-31'
mcSamples <- 1
minEventDuration <- 0 #[hr] minimum flood event length.
minRecordLength <- 10 #[yrs] minimum number of years on record for a gage to be included

minADCPMeas <- 20 #min depth stage adcp measurements used to convert bankfull stage to bankfull depth

#distributed computing across cpu-node
# tar_option_set(
#   controller = crew_controller_local(workers = 32)
# )

## MAP PIPELINE OVER STATIC BRANCHES OF GAGES
gageAnalysis <- tar_map(
  values = tibble( # Setup static branching
    method_function = rlang::syms(c("getBasinGages")),
    huc4 = c('0108')#, '0107', '0502', '1710', '0204')
  ),
  names='huc4',
  
  ## PREP GAGES
  tar_target(gagesBasin, method_function(path_to_data, huc4, gageRecordStart, gageRecordEnd)),
  tar_target(gage, prepGage(path_to_data, gagesBasin),
             pattern=map(gagesBasin),
             iteration='list'),
  tar_target(gageRecord, prepFlowRecord(gage, gageRecordStart, gageRecordEnd, minRecordLength),
             pattern=map(gage),
             iteration='list'),
  
  ## PREP BASIN DATA
  tar_target(basinData, buildBasinDataPackage(path_to_data, huc4)),

  ## BUILD BANKFULL STAGE MODEL
  tar_target(bankfullModel, buildGageBankfullModel(gage, minADCPMeas),
             pattern=map(gage),
             iteration='list'),
  tar_target(bankfull_mc, setupMonteCarlo(mcSamples, bankfullModel, 'justone'),
             pattern=map(bankfullModel),
             iteration='list'), #(see function for the third parameter)
  
  ## BUILD FLOOD EVENT ANALYSIS (INC. FLOODPLAIN PROFILES AND INUNDATED AREAS)
  tar_target(floodEvents, buildEventAnalysis(path_to_data, huc4, gageRecord, gage, bankfull_mc, BHGmodel, minEventDuration),
             pattern=map(gage, gageRecord, bankfull_mc), #branch over mc realization of bankfull stage
             iteration='list'),
  
  ## BUILD FLOODWATER PROFILES FROM GAGED FLOODS (for validation)
  # tar_target(floodwaterProfiles, buildFloodwaterProfile(path_to_data, huc4, basinData, gage, floodEvents),
  #            pattern=map(gage, floodEvents),
  #            iteration='list'),
  
  ## SUMMARISE MONTE CARLO SIMULATIONS INTO MEAN AND SIGMA
  tar_target(basinFloodRecord, buildEventAnalysis_full(floodEvents, gage)),
  
  ## UPSCALE TO RIVER NETWORK
 # tar_target(upscalingModel_returnperiod, buildUpscalingModel_returnperiod(huc4, basinFloodRecord, gageRecordStart, gageRecordEnd)),
  tar_target(upscalingModel, buildUpscalingModel(huc4, basinFloodRecord, gageRecordStart, gageRecordEnd)),
  tar_target(basinModel, buildNetworkModel(path_to_data, huc4, upscalingModel, BHGmodel, basinData)),
  tar_target(basinAnalysis, runNetworkModel(path_to_data, huc4, basinData, basinModel, benchmark_index))
)


list(
  ## PREP BANKFULL MODELs AND DATA
  tar_target(BHGmodel, modelsBHG()),
  tar_target(BHGdata, dataBHG()),
  
  ## RUN GAGE ANALYSIS (STATIC BRANCHED)
  gageAnalysis,
  
  ## VALIDATION
  tar_combine(combined_bankfullModel, gageAnalysis$bankfullModel, command = c(!!!.x)),
  tar_target(bankfullModelValidation, validateBankfullModel(combined_bankfullModel, BHGdata))
)