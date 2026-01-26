## Utility functions
## Craig Brinkerhoff
## Fall 2025


#' grabAllGages
#'
#' Grabs all streamgages from the NHD
#'
#' @param huc4 basin ID code
#'
#' @return streamgages joined a priori to the NHD
grabAllGages <- function(huc4){
  huc2 <- substr(huc4, 1, 2)

  network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
  network_gages <- network_gages %>%
    dplyr::select('GageID')
  return(network_gages)
}


#' cleanUpDF
#'
#' Makes sure that same training data are used for V_flood, Q_flood, and Q_total machine learning models
#'
#' @param df model training dataframe
#'
#' @return filtered model training dataframe
cleanUpDF <- function(df){
  out <- df %>%
    tidyr::drop_na(V_m3) %>%
    tidyr::drop_na(Qexc_m3dy) %>%
    tidyr::drop_na(Q_m3dy)

  return(out)
}


#' cleanUpDF_timing
#'
#' Makes sure that good training data used for timing model
#'
#' @param df model training dataframe
#'
#' @return filtered model training dataframe
cleanUpDF_timing <- function(df){
  out <- df %>%
    tidyr::drop_na(meanMonth)

  return(out)
}


#' cleanUpGages
#'
#' Ensure that model training dataframe and dataframe for plotting gages agree
#'
#' @param gages_df_combined Dataframe of all streamgages
#' @param modelDF Dataframe of all streamgages used for model training
#'
#' @return filtered gage dataframe for mapping
cleanUpGages <- function(gages_df_combined, modelDF){
  out <- gages_df_combined %>%
    dplyr::inner_join(modelDF, by=c('site_no'='GageID')) %>%
    dplyr::distinct() # a few duplicate rows, presumably from spatial joining of physio upstream. Not a problem in the modelDF

  return(out)
}


#' getBasinGages
#'
#' Grabs streamgages that meet our QAQC filters
#'
#' @param huc4id NHD basin code
#' @param gageRecordStart Start of time period for gage record
#' @param gageRecordEnd End of time period for gage record
#'
#' @return streamgage IDs that pass our quality control
getBasinGages <- function(huc4id, gageRecordStart, gageRecordEnd){
  set.seed(435)

  #grab all  huc8 ids to loop through
  huc2 <- substr(huc4id, 1, 2)
  huc8s <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU8.shp')) %>%
    dplyr::filter(substr(huc8, 1, 4)==huc4id)

  huc8ids <- huc8s$huc8

  #loop through hub8 basins, grabbing gages that pass initial QAQC parameters
  gages_fin <- data.frame()
  for(i in huc8ids){
    gages <- tryCatch({dataRetrieval::whatNWISdata(huc = i,
                                        parameterCd = c('00065', '00060'), #discharge parameter code [cfs]
                                        startDt = gageRecordStart, #only checks the gagues were active between these dates, need to pull the gage record to calucalte actual periods of redcord, etc. (see below)
                                        endDt = gageRecordEnd)},
                      error=function(x){data.frame('site_no'=i,
                                                  'parm_cd'=NA,
                                                  'data_type_cd'=NA)}) #trycatch for huc8s without any gages can't return anything

    #remove 'water quality' discharge measurements tagged under qw
    gages <- gages %>%
      dplyr::select(c('site_no', 'parm_cd', 'data_type_cd')) %>%
      dplyr::filter((parm_cd == '00060' & data_type_cd == 'uv')) %>% 
      dplyr::select(c('site_no', 'parm_cd'))

    gages_fin <- rbind(gages_fin, gages)
  }

  if(nrow(gages_fin)==0){
    return(NA)
  }

  return(gages_fin$site_no)
}


#' getBasinGagesVal
#'
#' Grabs all ids EXCEPT the one withheld for jacknife regression
#'
#' @param BHGmodel_jacknife Bankfull hydraulics models dataframe
#'
#' @return Streamgage IDs sans the withheld streamgage
getBasinGagesVal <- function(BHGmodel_jacknife){

    gages <- BHGmodel_jacknife$GageID
    gages <- gages[gages != 'ungaged']
    return(gages)
}


#' makeGageDF
#'
#' Makes a dataframe from list of streamgages
#'
#' @param gage Gage record dataframe
#' @param model_gages model river network
#' @param huc4 basin code
#'
#' @return Dataframe of streamgages in a basin that are also in the model
makeGageDF <- function(gage, model_gages, huc4){
  if(nrow(gage %>% dplyr::bind_rows())==0) {return(data.frame())}
  out <- gage %>%
    dplyr::bind_rows() %>%
    dplyr::filter(site_no %in% model_gages$GageID) %>%
    dplyr::mutate(huc4 = huc4) %>%
    dplyr::relocate(huc4)

  return(out)
}


#' tallyReaches
#'
#' count number of modeled reaches in a basin
#'
#' @param basinPredictions hydrography dataframe of modeled predictions
#'
#' @return number of modeled reaches in a basin
tallyReaches <- function(basinPredictions){
    out <- nrow(basinPredictions)/4

    return(out)
}


#' fixGeometries
#'
#' corrects geometry errors in sf object hydrography (if existing)
#'
#' @param rivnet hydrography as sf object
#'
#' @return sf object hydrography with fixed geometries
fixGeometries <- function(rivnet){
  curveLines <- dplyr::filter(rivnet, sf::st_geometry_type(rivnet) == 'MULTICURVE')
  if(nrow(curveLines) > 0){ #if saved as a curve, recast geometry as a line
    rivnet <- sf::st_cast(rivnet, 'MULTILINESTRING')
  }
  
  return(rivnet)
}
