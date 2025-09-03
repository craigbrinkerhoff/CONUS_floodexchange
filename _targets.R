# _targets.R file
# Craig Brinkerhoff
# Summer 2025
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
source('src/functions_v3.R')
source('src/utils.R')
source('src/figures.R')

# usgs dem vertical accuracy: https://www.mdpi.com/2072-4292/14/4/940

###### PARAMETERS ######
gageRecordStart <- '1983-01-01'
gageRecordEnd <- '2023-12-31'
minRecordLength <- 20 #[yrs] minimum number of years on record for a gage to be included
minADCPMeas <- 20 #min depth stage adcp measurements for AHG
minAHGr2 <- 0.30 #min depth AHG fit

#ml parameters
nInnerFolds <- 10
nOuterFolds <- 5
numGrid <- 50
numRepeats <- 1

###### DISTRIBUTED COMPUTING SETUP ######
# tar_option_set(
#   controller = crew_controller_local(workers = 5)
# )

vol_grids <- list.files('data/path_to_data/CONUS_connectivity_data/volume_validation/', pattern = "\\.tif$", full.names=TRUE)



###### BUILD TRAINING SET BY HUC4 BASIN ######
gageAnalysis <- tar_map(
  values = tibble( # Setup static branching
    method_function = rlang::syms(c("getBasinGages")),
    huc4 = c('0101', '0102', '0103', '0104', '0105', '0106', '0107', '0108', '0109', '0110',
          '0202', '0203', '0204', '0205', '0206', '0207', '0208',
          '0301', '0302', '0303', '0304', '0305', '0306', '0307', '0308', '0309', '0310', '0311', '0312', '0313', '0314', '0315', '0316', '0317', '0318',
          '0401', '0402', '0403', '0404', '0405', '0406', '0407', '0408', '0409', '0410', '0411', '0412', '0413', '0414', '0420', '0427', '0429', '0430', '0431',
          '0501', '0502', '0503', '0504', '0505', '0506', '0507', '0508', '0509', '0510', '0511', '0512', '0513', '0514',
          '0601', '0602', '0603', '0604',
          '0701', '0702','0703', '0704', '0705', '0706', '0707', '0708', '0709', '0710', '0711', '0712', '0713', '0714',
          '0801', '0802', '0803', '0804', '0805', '0806', '0807', '0808', '0809',
          '0901', '0902', '0903', '0904',
          '1002', '1003', '1004', '1005', '1006', '1007', '1008', '1009', '1010', '1011', '1012', '1013', '1014', '1015', '1016', '1017', '1018', '1019', '1020', '1021', '1022', '1023', '1024', '1025', '1026', '1027', '1028', '1029', '1030',
          '1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108', '1109', '1110', '1111', '1112', '1113', '1114',
          '1201', '1202', '1203', '1204', '1205', '1206', '1207', '1208', '1209','1210','1211',
          '1301','1302','1303','1304','1305','1306','1307','1308','1309',
          '1401','1402','1403','1404','1405','1406','1407','1408',
          '1501','1502','1503','1504','1505','1506','1507','1508',
          '1601','1602','1603','1604','1605','1606',
          '1701','1702','1703','1704','1705','1706','1707','1708','1709','1710','1711','1712',
          '1801','1802','1803','1804','1805','1806','1807','1808','1809','1810')
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
  tar_target(gageFlux, buildGageFloodFunctions(huc4, BHGmodel, gageQexc)), #the specifics of this calculation are hardcoded in gageQexc, but we add the description here
  #tar_target(gageFluxFlowProb, selectDepthandQ(gageFlux)), #USE THIS FUNCTION TO PICK THE FLOW PROBABILITY YOU'RE MODELING (see src/utils.R)
  tar_target(gageVolume, runDEMModel(huc4, gageFlux)),

  ## PREP FOR ML
  tar_target(gageForModel, addOtherNHDFeatures(gageVolume, huc4)),
  tar_target(conusForModel, buildCONUSnetwork(huc4, BHGmodel)),

  ## PREP FOR MAPPING AND SUMMARIZING
  tar_target(gage_df, makeGageDF(gage, gageForModel, huc4)),

  ## VALIDATE VOLUME BATHTUB MODEL per HUC4
  tar_target(reaches_val, collectValReaches(huc4)),#grab reaches joined a priori to gage network (using gages as proxy for USGS volume model mainstems (b/c they are calibrated to these specific gages))
  tar_target(depths_val, wrangleDepthGrids(huc4, reaches_val, volVal)),
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

  ## ASSIGN HUC4 TO VOLUME VALIDATION DATA
  tar_target(volVal, assignVolVals(vol_grids)),

  ## RUN HUC4 ANALYSIS
  gageAnalysis,

  ## COMBINE HUC4 OBJECTS
  tar_combine(gageVolume_val_combined, gageAnalysis$gageVolume_val, command = dplyr::bind_rows(!!!.x)),
  tar_combine(gageForModel_combined, gageAnalysis$gageForModel, command = dplyr::bind_rows(!!!.x)), #gages for training
  tar_target(modelDF, cleanUpDF(gageForModel_combined)), #gages for training
  tar_combine(gages_df_combined, gageAnalysis$gage_df, command=dplyr::bind_rows(!!!.x)), #gages for map
  tar_target(gagesDF, cleanUpGages(gages_df_combined, modelDF)), #gages for map
  tar_combine(conusDF, gageAnalysis$conusForModel, command = dplyr::bind_rows(!!!.x)), #conus for deploy

  ## TRAIN ML MODELS
  tar_target(model_V_eval, trainModelEval_V(modelDF, nInnerFolds, nOuterFolds, numGrid, numRepeats)),
  tar_target(model_V, trainModelFin_V(modelDF, nInnerFolds, numGrid)),
  tar_target(model_Q_eval, trainModelEval_Q(modelDF, nInnerFolds, nOuterFolds, numGrid, numRepeats)),
  tar_target(model_Q, trainModelFin_Q(modelDF, nInnerFolds, numGrid)),

  ## PREDICT MEAN MONTHLY TAU ACROSS UNITED STATES
  tar_target(conus_fin_1, deployModel(conusDF, model_Q, model_V, 1)), #Jan
  tar_target(conus_fin_2, deployModel(conusDF, model_Q, model_V, 2)), #Feb
  tar_target(conus_fin_3, deployModel(conusDF, model_Q, model_V, 3)), #Mar
  tar_target(conus_fin_4, deployModel(conusDF, model_Q, model_V, 4)), #Apr
  tar_target(conus_fin_5, deployModel(conusDF, model_Q, model_V, 5)), #May
  tar_target(conus_fin_6, deployModel(conusDF, model_Q, model_V, 6)), #Jun
  tar_target(conus_fin_7, deployModel(conusDF, model_Q, model_V, 7)), #Jul
  tar_target(conus_fin_8, deployModel(conusDF, model_Q, model_V, 8)), #Aug
  tar_target(conus_fin_9, deployModel(conusDF, model_Q, model_V, 9)), #Sep
  tar_target(conus_fin_10, deployModel(conusDF, model_Q, model_V, 10)), #Oct
  tar_target(conus_fin_11, deployModel(conusDF, model_Q, model_V, 11)), #Nov
  tar_target(conus_fin_12, deployModel(conusDF, model_Q, model_V, 12)), #Dec

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
  tar_target(fig_map, makeMap(conus_fin_1)),
  #tar_target(fig_map_storage, makeMap_storage(conus_fin)),
  tar_target(fig_SO, makeStreamOrderFig(conus_fin)),
  tar_target(fig_validationCalculation, makeCalculationValFig(gageVolume_val_combined, gageQexc_val)),
  tar_target(fig_validationML, makeMLValFig(model_Q_eval, model_V_eval)),
  tar_target(fig_gageMap, makeGageMap(gagesDF)),
  tar_target(fig_VIP, makeVIPPlot(model_Q, model_V))
)










#snapping_thresh <- 5000 #[m]
#area_thresh_perc <- 0.10
#minWidth <- 10 #m, minimum bankful river width to map inundation at gages (remember dem res is 10m)
#maxProbDiff <- 0.01 #max tolerance when matching flow events to probabilities (%)
#usgs_maps <- sf::st_read('data/path_to_data/CONUS_connectivity_data/USGS_models/USGSmodels_area.shp')


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
