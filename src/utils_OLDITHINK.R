makeValGageRecord <- function(gageRecord, BHGmodel_jacknife){

  ids <- c(do.call(rbind, (lapply(gageRecord, function(x) nrow(x) > 0))))
  ids <- which(ids == TRUE)
  gageRecord <- gageRecord[ids]

  ids2 <- c(do.call(rbind, (lapply(gageRecord, function(x) x[1,]$site_no %in% BHGmodel_jacknife$GageID))))
  gageRecord <- gageRecord[ids2]

  return(gageRecord)
}

#' Note: Sf catchments are handled internally (using sql queries) because for some reason dplyr::filter won't work on these objects passed betweenn functions... ANd the sql is likely faster anyway
buildBasinDataPackage <- function(huc4) {
  huc2 <- substr(huc4, 1, 2)
  
  #setup basin shapefiles
  dem <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/elev_cm.tif')) #[cm]
  d8 <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/fdr.tif')) #[cm]

  prelist <- list('dem'=dem,
       'd8'=d8)
  outlist <- lapply(prelist, terra::wrap) #terra holds C++ pointers in memory so making lists of objects (for distributed computing won't work- we need to 'wrap' the objects into a packed state that can be sent over a serialized connection)
  
  
  return('terra'=outlist)
}


unlistGages <- function(gage_list){
  df <- dplyr::bind_rows(gage_list)
  return(df)
}



unshapefify <- function(shp){
  shp <- shp %>% 
    sf::st_drop_geometry()
  
  return(shp)
}