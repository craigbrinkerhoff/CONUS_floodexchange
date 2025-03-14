## Utility functions
## Craig Brinkerhoff
## Summer 2024


getBasinGages <- function(huc4id, gageRecordStart, gageRecordEnd, benchmark_index){
  huc2 <- substr(huc4id, 1, 2)
  huc8s <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU8.shp')) %>%
    dplyr::filter(substr(huc8, 1, 4)==huc4id)
  
  
  huc8ids <- huc8s$huc8

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
  
  return(gages_fin$site_no)
}





modelsBHG <- function(){
  dataset <- readr::read_csv('data/bhg_us_database_bieger_2015.csv') %>% #available by searching for paper at https://swat.tamu.edu/search
    dplyr::select(c('Physiographic Division', '...9', '...11','...13')) #some necessary manual munging for colnames from dataset
  
  colnames(dataset) <- c('DIVISION', 'DA_km2', 'Qb_cms','Wb_m')
  
  dataset$Wb_m <- as.numeric(dataset$Wb_m)
  dataset$Qb_cms <- as.numeric(dataset$Qb_cms)
  dataset$DA_km2 <- as.numeric(dataset$DA_km2)
  
  dataset <- tidyr::drop_na(dataset)
  
  division <- toupper(sort(unique(dataset$DIVISION))) #make lowercase and sort
  
  #build models, grouped by physiographic region
  models <- dplyr::group_by(dataset, DIVISION) %>%
    dplyr::do(model_Wb = lm(log10(Wb_m)~log10(DA_km2), data=.),
              model_Qb = lm(log10(Qb_cms)~log10(DA_km2), data=.)) %>% #fit models by physiographic regions
    dplyr::summarise(a_Wb = 10^(model_Wb$coef[1]), #model intercept
                     b_Wb = model_Wb$coef[2], #model exponent
                     r2_Wb = summary(model_Wb)$r.squared, #model performance
                     mean_residual_Wb = mean(model_Wb$residuals, na.rm=T),
                     sd_residual_Wb = sd(model_Wb$residuals, na.rm=T),
                     see_Wb = sd(model_Wb$residuals, na.rm=T),
                     a_Qb = 10^(model_Qb$coef[1]),
                     b_Qb = model_Qb$coef[2],
                     r2_Qb = summary(model_Qb)$r.squared) %>%
    dplyr::mutate(division = division)
  
  
  models[models$division == "INTERMONTANE PLATEAU",]$division <- "INTERMONTANE PLATEAUS" #make sure names line up
  
  return(models)
}




dataBHG <- function(){
  dataset <- readr::read_csv('data/bhg_us_database_bieger_2015.csv') %>% #available by searching for paper at https://swat.tamu.edu/search
    dplyr::select(c('USGS Station No.', 'Physiographic Division', '...9', '...13')) #some necessary manual munging for colnames from dataset
  
  colnames(dataset) <- c('GageID', 'DIVISION', 'DA_gage_skm','Wb_m_obs')
  
  dataset$Wb_m_obs <- as.numeric(dataset$Wb_m_obs)
  dataset$DA_gage_skm <- as.numeric(dataset$DA_gage_skm)
  
  dataset <- tidyr::drop_na(dataset) #only keep those with a usgs site and a non-NA Qb value (some data didn't report Qb)
  
  return(dataset)
}






fixGeometries <- function(rivnet){
  curveLines <- dplyr::filter(rivnet, sf::st_geometry_type(rivnet) == 'MULTICURVE')
  if(nrow(curveLines) > 0){ #if saved as a curve, recast geometry as a line
    rivnet <- sf::st_cast(rivnet, 'MULTILINESTRING')
  }
  
  return(rivnet)
}





# setupMonteCarlo <- function(m,bankfullModel, variable){
#   set.seed(654)
  
#   if(variable == 'justone'){
#     out <- rnorm(1000,bankfullModel$bankfull_stage_mu_m, bankfullModel$bankfull_stage_se_m)
#     out <- mean(out) #no MC, just take the median of the simulated distribution 
#   }
#   else{ #normal MC
#     out <- rnorm(m,bankfullModel$bankfull_stage_mu_m, bankfullModel$bankfull_stage_se_m)
#   }
  
#   return(list('stage_m'=out,
#               'depth_m'=bankfullModel$bankfull_depth_mu_m))
# }



#' Note: Sf catchments are handled internally (using sql queries) because for some reason dplyr::filter won't work on these objects passed betweenn functions... ANd the sql is likely faster anyway
buildBasinDataPackage <- function(huc4) {
  huc2 <- substr(huc4, 1, 2)
  
  #setup basin shapefiles
  dem <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/elev_cm.tif')) #[cm]
  d8 <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/fdr.tif')) #[cm]

  prelist <- list('dem'=dem,
       'd8'=d8)
  outlist <- lapply(prelist, terra::wrap) #terra holds C++ pointers in memory so making lists of objects (for distributed computing won't work- we need to 'wrap' the objects into a packed state that can be sent over a serialized connection)
  
  
  #build waterbody type lookup table
  # sf::sf_use_s2(FALSE)
  # waterbodies <-  sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
  #                             layer = 'NHDWaterbody',
  #                             quiet = TRUE) %>%
  #   sf::st_zm()
  # waterbodies <- fixGeometries(waterbodies)
  # flowlines <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
  #                             layer = 'NHDFlowline',
  #                             quiet = TRUE) %>%
  #   sf::st_zm()
  # flowlines <- fixGeometries(flowlines)
  # waterbodies_fin <- sf::st_join(waterbodies, flowlines, join=sf::st_contains, largest=TRUE)
  
  # waterbodyLookUp <- data.frame('NHDPlusID_reach'=waterbodies_fin$NHDPlusID.y,
  #                               'NHDPlusID_polygon'=waterbodies_fin$NHDPlusID.x,
  #                               'waterbody_type'=as.character(waterbodies_fin$FType.x))
  
  # return(list('terra'=outlist,
  #             'waterbodyLookUp'=waterbodyLookUp))
  return('terra'=outlist)
}



unlistGages <- function(gage_list){
  df <- dplyr::bind_rows(gage_list)
  return(df)
}




prepGWD <- function(){
    #prep barriers dataset
    states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
    states <- dplyr::filter(states, !(NAME %in% c('Alaska',
                                                'American Samoa',
                                                'Commonwealth of the Northern Mariana Islands',
                                                'Guam',
                                                'District of Columbia',
                                                'Puerto Rico',
                                                'United States Virgin Islands',
                                                'Hawaii'))) #remove non CONUS states/territories

    states <- sf::st_union(states) %>%
        sf::st_transform(crs=sf::st_crs(4326))

    #append dam drainage areas via the GDW (Global Dam Watch database- https://www.globaldamwatch.org/database)
    barriers <- sf::st_read('data/path_to_data/CONUS_sediment_data/GDW_barriers_v1_0.shp') %>%
        sf::st_intersection(states) %>%
        dplyr::filter(INSTREAM == 'Instream') %>% #remove offstream barriers not on hydrosheds/hydrobasins framework
        dplyr::filter(DAM_TYPE %in% c('Dam', 'Lake Control Dam', 'Low Permeable Dam')) %>% #remove locks, which don't trap sediment in the way we care about here
        dplyr::mutate(RESV_CATCH_SKM = CATCH_SKM)%>%
        dplyr::select(c('HYRIV_ID', 'RESV_CATCH_SKM'))
    
    return(barriers)
}