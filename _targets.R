# _targets.R file
# Craig Brinkerhoff
# Spring 2025
# Master pipeline for river-floodplain exchange model

#necessary global packages for pipelining.
library(targets)
library(tarchetypes)
library(dplyr)
library(tidyr)
library(crew)
library(crew.cluster)
library(ggplot2)
library(tibble)

#load functions
#source("src/functions_floodplain.R")
source('src/functions_v3.R')
#source("src/validation.R")
source('src/utils.R')
source('src/figures.R')
#source('src/validation_tau.R')

##BIG TO DOS
#(3/21/25): Probably go back and make sure we are filtering and QAQC'ing the nhd properly..
# (3/21/25): Think about log transform bias in horton scaling (is it realll.....?)

# usgs dem vertical accuracy: https://www.mdpi.com/2072-4292/14/4/940


###### PARAMETERS ######
gageRecordStart <- '1983-01-01'
gageRecordEnd <- '2023-12-31'
minRecordLength <- 20 #[yrs] minimum number of years on record for a gage to be included
minADCPMeas <- 20 #min depth stage adcp measurements for AHG
minAHGr2 <- 0.30 #min depth AHG fit
#snapping_thresh <- 5000 #[m]
#area_thresh_perc <- 0.10
#minWidth <- 10 #m, minimum bankful river width to map inundation at gages (remember dem res is 10m)
#maxProbDiff <- 0.01 #max tolerance when matching flow events to probabilities (%)
#usgs_maps <- sf::st_read('data/path_to_data/CONUS_connectivity_data/USGS_models/USGSmodels_area.shp')


###### DISTRIBUTED COMPUTING SETUP ######
# tar_option_set(
#   controller = crew_controller_local(workers = 6)
# )

#local setup (used for all but the big jobs with internal parallelism)
# controller_local <- crew_controller_local(
#   name = "my_local_controller",
#   workers = 3)

# #slurm setup (to submit external, big jobs for the heft targets with internal parallelism)
# controller_slurm <- crew_controller_slurm(
#   name = "my_slurm_controller",
#   workers = 1,
#   options_cluster = crew.cluster::crew_options_slurm(#script_lines = "module load R",
#                                                     partition = 'day',
#                                                     memory_gigabytes_required = 180,
#                                                     cpus_per_task = 35,
#                                                     time_minutes = 60*12, #12hr time limit
#                                                     log_output = "logs/crew_out_%A.txt",
#                                                     log_error = "logs/crew_err_%A.txt"))

# #define default settings
# tar_option_set(
#   controller = crew_controller_group(controller_local, controller_slurm),
#   resources = tar_resources(
#     crew = tar_resources_crew(controller = "my_local_controller")
#   )
# )



###### BUILD TRAINING SET BY HUC4 BASIN ######
gageAnalysis <- tar_map(
  values = tibble( # Setup static branching
    method_function = rlang::syms(c("getBasinGages")),
    huc4 = c('0101', '0102', '0103', '0104', '0105', '0106', '0107', '0108', '0109', '0110')
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
  tar_target(depthAHG, buildDepthAHG(gage, minADCPMeas),
            pattern=map(gage),
            iteration='list'),
  tar_target(gagePrepped, prepBankfullHydraulics('deploy', gageRecord, gage, BHGmodel_jacknife, BHGmodel, minAHGr2, gageRecordStart, gageRecordEnd)),

  ## RUN GAGE FLOOD VOLUME & DISCHARGE
  tar_target(gageQexc, calc_Qexc('deploy', gageRecord, gagePrepped, depthAHG),
            pattern=map(gageRecord),
            iteration='list'),
  tar_target(gageFlux, buildGageFloodFunctions(huc4, BHGmodel, gageQexc, 'average of event means')), #the specifics of this calculation are hardcoded in gageQexc, but we add the description here
  tar_target(gageVolume, runDEMModel(huc4, gageFlux)),

  ## VALIDATE VOLUME BATHTUB MODEL
  tar_target(reaches_val, collectValReaches(huc4)),#grab reaches joined a priori to gage network (using gages as proxy for USGS volume model mainstems (b/c they are calibrated to these specific gages))
  tar_target(depths_val, wrangleDepthGrids(huc4, reaches_val)),
  tar_target(gageFlux_val, buildGageFloodFunctions_volumeval(huc4, BHGmodel, depths_val)), #also passes along the observed volumes, compared to the normal function above
  tar_target(gageVolume_val, runDEMModel(huc4, gageFlux_val))
)


###### RUN PIPELINE ######
list(
  ## PREP BANKFULL MODELS AND DATA
  tar_target(BHGdata, dataBHG()),
  tar_target(BHGmodel, modelsBHG()),

  ## PREP FLOW BARRIER DATASET
  tar_target(GWD, prepGWD()),
  
  ## RUN HUC4 ANALYSIS
  gageAnalysis,

  ## COMBINE HUC4 OBJECTS
  tar_combine(gageVolume_val_combined , gageAnalysis$gageVolume_val, command = dplyr::bind_rows(!!!.x)),

  #EXCHANGE TIME VALIDATION VIA JACKNIFE REGRESSION (more or less a LOOCV for the regression models)
  tar_target(BHGmodel_jacknife, modelsJacknifeBHG()),
  tar_target(gagesBasin_val, getBasinGagesVal(BHGmodel_jacknife)),
  tar_target(gage_val, prepGage(gagesBasin_val),
            pattern=map(gagesBasin_val),
            iteration='list'),
  tar_target(gageRecord_val, prepFlowRecord(gage_val, gageRecordStart, gageRecordEnd, minRecordLength),
            pattern=map(gage_val),
            iteration='list'),
  tar_target(depthAHG_val, buildDepthAHG(gage_val, minADCPMeas),
          pattern=map(gage_val),
          iteration='list'),
  tar_target(gagePrepped_val, prepBankfullHydraulics('val', gageRecord_val, gage_val, BHGmodel_jacknife, BHGmodel, minAHGr2, gageRecordStart, gageRecordEnd)),
  tar_target(gageQexc_val, calc_Qexc('val', gageRecord_val, gagePrepped_val, depthAHG_val),
          pattern=map(gageRecord_val),
          iteration='list'),

  ##FIGURES
  tar_target(fig_validationCalculation, makeCalculationValFig(gageVolume_val_combined, gageQexc_val))


  ## PREP VALIDATION MAPS
  # tar_target(preppedFEMA, prepFEMA('1')),
  
  ## BUILD HUC2 UPSCALING MODELS
  # buildUpscalingModels,

    # ## COMBINE RESULTS ACROSS REGIONS
  # tar_combine(gage_combined, list(basinAnalysis_01$gages_fixed, basinAnalysis_02$gages_fixed), command = dplyr::bind_rows(!!!.x)),
  # tar_combine(val_FEMA_combined, list(basinAnalysis_01$val_FEMA, basinAnalysis_02$val_FEMA), command = dplyr::bind_rows(!!!.x)),
  # tar_combine(val_USGS_combined, list(basinAnalysis_01$val_USGS, basinAnalysis_02$val_USGS), command = dplyr::bind_rows(!!!.x)),
  # tar_combine(val_USGSvols_combined, list(basinAnalysis_01$val_USGSvols, basinAnalysis_02$val_USGS), command = dplyr::bind_rows(!!!.x)),

  # ## COMBINE RESULTS WITHIN REGIONS
  # tar_combine(model_combined_01, basinAnalysis_01$basinAnalysis_fin_df, command = dplyr::bind_rows(!!!.x)),
  # tar_combine(model_combined_02, basinAnalysis_02$basinAnalysis_fin_df, command = dplyr::bind_rows(!!!.x)),
  
  # tar_combine(hortonResults_combined_01, basinAnalysis_01$hortonResults, command = dplyr::bind_rows(!!!.x)),
  # tar_combine(hortonResults_combined_02, basinAnalysis_02$hortonResults, command = dplyr::bind_rows(!!!.x)),

  # ## FIGURES
  # tar_target(valFEMAFig, makeValFEMA(val_FEMA_combined, gage_combined, BHGdata)),
  # tar_target(valUSGSareaFig, makeValUSGSarea(val_USGS_combined)),
  # tar_target(valUSGSvolFig, makeValUSGSvol(val_USGSvols_combined)),
  # #tar_target(streamOrderFig, makeHortonScalingFig(hortonResults_combined, model_combined)),
  # tar_target(hrt_test, residenceTimeFig(basinAnalysis_hrt_0107)),
  # tar_target(regulationFig_01, makeRegulationFig(model_combined_02)),
  # tar_target(mapFig, makeMapFig(usgs_maps)) #UPDATE


  ##COMBINE OUTPUTS
  # tar_combine(output_combined, buildGageDataset$output, command=dplyr::bind_rows(!!!.x)),
  # tar_combine(gageTau_combined, buildGageDataset$gageTau_summ, command=dplyr::bind_rows(!!!.x)),
)





  # ## PREP GAGES
  # tar_target(gagesBasin, method_function(huc4, gageRecordStart, gageRecordEnd)),
  # tar_target(gage, prepGage(gagesBasin),
  #           pattern=map(gagesBasin),
  #           iteration='list'),
  # tar_target(gageRecord, prepFlowRecord(gage, gageRecordStart, gageRecordEnd, minRecordLength),
  #           pattern=map(gage),
  #           iteration='list'),
  # tar_target(gages_fixed, unlistGages(gage)),

  # ## BUILD DEPTH AHG RELATION
  # tar_target(depthAHG, buildDepthAHG(gage, minADCPMeas),
  #           pattern=map(gage),
  #           iteration='list'),



  ## PREP BASIN DATA
  # tar_target(basinData, buildBasinDataPackage(huc4), resources= tar_resources(crew = tar_resources_crew(controller = 'my_slurm_controller'))),
  #tar_target(gages_fixed, unlistGages(gage_01)),

  ## CALCULATE UPSTREAM FLOW REGULATION BY LAKES/RESERVOIRS
  # tar_target(basinAnalysis_fin, regulateFlooding(basinAnalysis_hrt, GWD, area_thresh_perc, huc4, snapping_thresh), resources= tar_resources(crew = tar_resources_crew(controller = 'my_slurm_controller'))),
  # tar_target(basinAnalysis_fin_df, unshapefify(basinAnalysis_fin)),


  # ## HORTON SCALING TO CAPTURE INUNDATION IN STREAMS < 10M WIDE
  # tar_target(hortonResults, hortonScaling(basinAnalysis_fin, huc4)),

  # ## VALIDATE AT REACHES WITH USGS INUNDATION MODELS (unfortunately, the flow probabilities are hardcoded throughout, not just at htis target...)
  # tar_target(val_FEMA, valModelFEMA(huc4, preppedFEMA, basinAnalysis_fin, basinData)),
  # tar_target(val_USGS, valModelUSGS(huc4, basinAnalysis_fin, usgs_maps)),
  # tar_target(val_USGSvols, valModelUSGSvols(huc4, basinAnalysis_fin, val_USGS, c('q0_2', 'q0_5', 'q1', 'q2', 'q4', 'q10', 'q20','q50', 'q80', 'q90', 'q96', 'q98', 'q99', 'q99_5', 'q99_8')))




###### HUC2 WATER LEVEL RELATIONSHIPS RECIPE ######
# buildUpscalingModels <- tar_map(
#   values = tibble( # Setup static branching
#     method_function = rlang::syms(c("getBasinGages")),
#     huc2 = c('01', '02', '03')
#   ),
#   names='huc2',

#   ## PREP GAGES
#   # tar_target(gagesBasin, method_function(huc2, gageRecordStart, gageRecordEnd))
#   # tar_target(gage, prepGage(gagesBasin),
#   #           pattern=map(gagesBasin),
#   #           iteration='list'),
#   # tar_target(gageRecord, prepFlowRecord(gage, gageRecordStart, gageRecordEnd, minRecordLength),
#   #           pattern=map(gage),
#   #           iteration='list'),

#   # ## BUILD DEPTH AHG RELATION
#   # tar_target(depthAHG, buildDepthAHG(gage, minADCPMeas),
#   #           pattern=map(gage),
#   #           iteration='list')
  
#   # ## BUILD UPSCALING MODELS
#   # tar_target(upscalingModel, buildUpscalingModel(huc2, gageRecord, gage, BHGmodel, depthAHG, minAHGr2, gageRecordStart, gageRecordEnd)),
  
#   # ## BUILD REGION SPECIFIC UPSCALING FIGURES
#   # tar_target(upscalingPlot, upscalingFig(upscalingModel, huc2))
# )