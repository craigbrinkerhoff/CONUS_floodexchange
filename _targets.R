# _targets.R file
# Craig Brinkerhoff
# Winter 2026
# Master pipeline for floodwater residence time study

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
source('src/functions.R')
source('src/utils.R')
source('src/figures.R')

###### PARAMETERS ######
gageRecordStart <- '1983-01-01'
gageRecordEnd <- '2023-12-31'
minRecordLength <- 20 #[yrs] minimum number of years on record for a gage to be included
minADCPMeas <- 20 #min depth~stage in situ flow measurements for AHG
minAHGr2 <- 0.30 #min depth AHG fit
min_floods <- 3 #min number of floods on record to compute probabilities

#ml parameters
nInnerFolds <- 10
nOuterFolds <- 5
numGrid <- 50
numRepeats <- 1

###### DISTRIBUTED COMPUTING SETUP ######
# tar_option_set(
#   controller = crew_controller_local(workers = 9)
# )

vol_grids <- list.files('data/path_to_data/CONUS_connectivity_data/volume_validation/', pattern = "\\.tif$", full.names=TRUE)

###### BATCH PROCESSING BY HUC4 BASIN ######
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
  tar_target(allGages, grabAllGages(huc4)),

  ## CALCULATE FLOOD VOLUME & DISCHARGE AT GAGE
  tar_target(gageQexc, calc_Qexc('deploy', gageRecord, gagePrepped, depthAHG),
            pattern=map(gageRecord),
            iteration='list'),
  tar_target(gageFlux, buildGageFloodFunctions(huc4, BHGmodel, gageQexc, min_floods)),
  tar_target(gageVolume, runDEMModel(huc4, gageFlux)),

  ## PREP FOR ML
  tar_target(gageForModel, addOtherNHDFeatures(gageVolume, huc4)),
  tar_target(conusForModel, buildCONUSnetwork(huc4, BHGmodel)),

  ## APPLY MODEL TO BASIN AND SUMMARIZE BY STREAM ORDER
  tar_target(basinPredictions, predictBasin(huc4, conusForModel, model_Qf, model_V, model_Q)),
  tar_target(basinSummarySO, summarizeBasinSO(huc4, basinPredictions)),
  tar_target(basinSummary, summarizeBasin(huc4, basinPredictions)),
  tar_target(nReaches, tallyReaches(basinPredictions)),

  ## PREP FOR MAPPING
  tar_target(gage_df, makeGageDF(gage, gageForModel, huc4)),
  tar_target(mapTau_02, makeMapBasin(basinPredictions, 0.02)),
  tar_target(mapTau_10, makeMapBasin(basinPredictions, 0.10)),
  tar_target(mapTau_20, makeMapBasin(basinPredictions, 0.20)),
  tar_target(mapTau_50, makeMapBasin(basinPredictions, 0.50)),

  ## VALIDATE VOLUME MODEL per BASIN
  tar_target(reaches_val, collectValReaches(huc4)),
  tar_target(depths_val, wrangleDepthGrids(huc4, reaches_val, volVal)),
  tar_target(gageFlux_val, buildGageFloodFunctions_volumeval(huc4, BHGmodel, depths_val)),
  tar_target(gageVolume_val, runDEMModel(huc4, gageFlux_val))
)


###### RUN PIPELINE ######
list(
  ## PREP BANKFULL MODELS AND DATA
  tar_target(BHGdata, dataBHG()),
  tar_target(BHGmodel, modelsBHG()),

  ## PREP VALIDATION DATA
  tar_target(volVal, assignVolVals(vol_grids)),

  ## BATCH BASINS (see above)
  gageAnalysis,

  ## COMBINE BASIN OBJECTS FOR MODEL TRAINING
  tar_combine(gageVolume_val_combined, gageAnalysis$gageVolume_val, command = dplyr::bind_rows(!!!.x)),
  tar_combine(gageForModel_combined, gageAnalysis$gageForModel, command = dplyr::bind_rows(!!!.x)),
  tar_target(modelDF, cleanUpDF(gageForModel_combined)),
  tar_combine(gages_df_combined, gageAnalysis$gage_df, command=dplyr::bind_rows(!!!.x)),
  tar_target(gagesDF, cleanUpGages(gages_df_combined, modelDF)),
  tar_combine(allGages_combined, gageAnalysis$allGages, command = dplyr::bind_rows(!!!.x)),
  tar_combine(nReaches_combined, gageAnalysis$nReaches, command = vctrs::vec_c(!!!.x)),

  ## TRAIN ML MODELS
  tar_target(model_V_eval, trainModelEval_V(modelDF, nInnerFolds, nOuterFolds, numGrid, numRepeats)),
  tar_target(model_V, trainModelFin_V(modelDF, nInnerFolds, numGrid)),
  tar_target(model_Qf_eval, trainModelEval_Qf(modelDF, nInnerFolds, nOuterFolds, numGrid, numRepeats)),
  tar_target(model_Qf, trainModelFin_Qf(modelDF, nInnerFolds, numGrid)),
  tar_target(model_Q_eval, trainModelEval_Q(modelDF, nInnerFolds, nOuterFolds, numGrid, numRepeats)),
  tar_target(model_Q, trainModelFin_Q(modelDF, nInnerFolds, numGrid)),

  ## COMBINE BASIN OBJECTS POST-MODEL MAPPING & PREDICTION
  tar_combine(basinsList_02, gageAnalysis$mapTau_02, command = list(!!!.x)),
  tar_combine(basinsList_10, gageAnalysis$mapTau_10, command = list(!!!.x)),
  tar_combine(basinsList_20, gageAnalysis$mapTau_20, command = list(!!!.x)),
  tar_combine(basinsList_50, gageAnalysis$mapTau_50, command = list(!!!.x)),
  tar_combine(combined_basinSummarySO, gageAnalysis$basinSummarySO, command = dplyr::bind_rows(!!!.x)),
  tar_combine(combined_basinSummary, gageAnalysis$basinSummary, command = dplyr::bind_rows(!!!.x)),

  #QEXC COMPARISON VIA JACKNIFE REGRESSION & IN SITU QB
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

  ##PAPER FIGURES
  tar_target(fig_validationML, makeMLValFig(model_Qf_eval, model_V_eval, modelDF)),
  tar_target(fig_reachTauMap, makeReachMap(basinsList_02, basinsList_10, basinsList_20, basinsList_50)),
  tar_target(fig_reachTauMap_insets, makeReachMapInset(mapTau_50_0205, mapTau_50_0206, mapTau_50_0207, mapTau_50_0208, mapTau_50_0502, mapTau_50_0503, mapTau_50_0501, mapTau_50_0505)),
  tar_target(fig_streamorder, makeReachBoxplotsFig(combined_basinSummarySO)),
  tar_target(fig_basinExchange, makeBasinExchangeMap(combined_basinSummary)),

  ## EXTENDED FIGURES
  tar_target(fig_validationCalculation, makeCalculationValFig(gageVolume_val_combined, gageQexc_val)),
  tar_target(fig_gageMap, makeGageMap(gagesDF)),
  tar_target(fig_totalQVal, makeMLQFig(model_Q_eval, modelDF)),
  tar_target(fig_huc4s, makeHuc4Map())
)