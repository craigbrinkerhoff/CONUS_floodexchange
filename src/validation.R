## Validation functions
## Craig Brinkerhoff
## Winter 2024




valModelUSGS <- function(huc4, basinData, basinModel){
  shp <- sf::st_read('data/path_to_data/CONUS_connectivity_data/deerfield_nhd_reaches.shp')

  basinModel <- basinModel %>%
    dplyr::filter(NHDPlusID %in% shp$NHDPlusID)

  #dummy setup for now
  analysis <- runNetworkModel(huc4, basinData, basinModel)

  return(analysis)
}





prepFEMA <- function(){
    #read in fema
    states <- list.files('data/path_to_data/CONUS_connectivity_data/CONUS_FEMA')

    fema_shp <- data.frame()
    for(i in states){
        geodatabase <- list.files(paste0('data/path_to_data/CONUS_connectivity_data/CONUS_FEMA/', i), pattern = "\\.gdb$")
        
        #load in fema maps and filter for 100-yr flood polygons
        fema_temp <- sf::st_read(paste0('data/path_to_data/CONUS_connectivity_data/CONUS_FEMA/', i, '/', geodatabase),
                                    layer = 'S_FLD_HAZ_AR',
                                    quiet = TRUE) %>%
            dplyr::filter(FLD_ZONE %in% c('A', 'AE', 'AH', 'AO')) %>% #A zones correspond to 100-yr flood maps, built with a variety of methods. See README.
            dplyr::mutate(STATE = i) %>%
            dplyr::select(c('FLD_AR_ID', 'STATE', 'FLD_ZONE'))
        
        fema_shp <- rbind(fema_shp, fema_temp)
    }

    return(fema_shp)
}






valModelFEMA <- function(huc4id, preppedFEMA, basinModel, basinAnalysis){
    sf::sf_use_s2(FALSE)

    #prep
    huc2 <- substr(huc4id, 1, 2)

    #read in huc4 basin
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE)

    #preemptively filter fema map to our basin (saves time)
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)

    fema_shp <- preppedFEMA %>%
        sf::st_filter(basin)

    #only keep unit catchments that overlap the fema maps
    catch_keep <- unit_catchments %>%
        sf::st_intersects(fema_shp)
    catch_keep <- unique(unlist(catch_keep))

    basinAnalysis <- basinAnalysis[catch_keep,]

    #filter model for 100-yr AEP and reaches in FEMA map
    basinAnalysis <- basinAnalysis %>%
        dplyr::filter(!(is.na(Af_qFEMA_m2))) %>%
        dplyr::mutate(A_qFEMA_m2 = Af_qFEMA_m2 + (Wb_m*LengthKM*1000)) %>% #add the channel area back into the estimate (for validation)
        dplyr::select(c('NHDPlusID', 'GageID', 'StreamCalc', 'A_qFEMA_m2'))
    
    #remove catchments with narrow rivers from fema_shps (then calcualate inundated area)   
    unit_catchments <- unit_catchments %>%
        dplyr::filter(NHDPlusID %in% basinAnalysis$NHDPlusID)

    fema_shp <- fema_shp %>%
        sf::st_intersection(unit_catchments) %>%
        dplyr::mutate(A_femaAEP_m2 = sf::st_area(.))

    #build validation table
    out <- data.frame('huc4'=huc4id,
                    'n_catchments'=nrow(basinAnalysis),
                    'A_qFEMA_km2'=1e-6 * sum(as.numeric(basinAnalysis$A_qFEMA_m2), na.rm=T),
                    'A_femaAEP_km2'=1e-6 * sum(as.numeric(fema_shp$A_femaAEP_m2), na.rm=T))
    
    return(out)
}