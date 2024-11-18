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






valModelFEMA <- function(huc4id, preppedFEMA, basinAnalysis){
    sf::sf_use_s2(FALSE)

    #prep
    huc2 <- substr(huc4id, 1, 2)
    # huc8s <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU8.shp')) %>%
    #     dplyr::select(c('huc8'))

    #read in huc4 basin
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE)

    #preemptively filter fema map to our basin (saves time)
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)

    #to speed up below calculations
    fema_shp <- preppedFEMA %>%
        sf::st_filter(basin)

    #filter model for 100-yr AEP and reaches with modeled results (non-NA)
    basinAnalysis <- basinAnalysis %>%
        dplyr::filter(!(is.na(A_qFEMA_m2))) %>%
        dplyr::select(c('NHDPlusID', 'GageID', 'StreamCalc', 'A_qFEMA_m2'))
    
    #only keep model reaches that align with the fema maps
    catch_keep <- fema_shp %>%
        sf::st_intersects(basinAnalysis)
    catch_keep <- unique(unlist(catch_keep))

    basinAnalysis <- basinAnalysis[catch_keep,]
    
    #remove catchments with narrow rivers from fema_shps (then calcualate inundated area)   
    unit_catchments <- unit_catchments %>%
        dplyr::filter(NHDPlusID %in% basinAnalysis$NHDPlusID)
    
    basinAnalysis <- basinAnalysis %>% #unit_catchments %>%
    #    sf::st_join(huc8s, join=sf::st_intersects, largest=TRUE) %>% #if in multiple huc8s, assign to the one it's mostly in (this should handle just a few edge cases)
        dplyr::group_by(NHDPlusID) %>%
        dplyr::summarise(A_qFEMA_m2 = sum(A_qFEMA_m2, na.rm=T))

    fema_shp <- fema_shp %>%
        sf::st_intersection(unit_catchments) %>%  #make sure fema polygons don't align with this subset of modeled catchments. It's unfair to include FEMA maps outside of the subset of our model with predictions (that overlaps FEMA) and its unfair to include FEMA maps that extend beyond our catchment scheme
        sf::st_transform('+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs') %>% #equal area projection to match model
        dplyr::mutate(A_femaAEP_m2 = sf::st_area(.)) %>%
      #  sf::st_transform(sf::st_crs(huc8s)) %>% #put back for spatial join
      #  sf::st_join(huc8s, join=sf::st_intersects, largest=TRUE) %>%
        dplyr::group_by(NHDPlusID) %>%
        dplyr::summarise(A_femaAEP_km2 = 1e-6 * sum(as.numeric(A_femaAEP_m2), na.rm=T)) %>%
        dplyr::select(c('NHDPlusID', 'A_femaAEP_km2'))

    #build validation table
    out <- data.frame('huc4'=huc4id,
                    'NHDPlusID'=basinAnalysis$NHDPlusID,
                    'A_qFEMA_km2'=basinAnalysis$A_qFEMA_m2*1e-6)

    out <- out %>%
        dplyr::left_join(fema_shp, by='NHDPlusID') %>%
        tidyr::drop_na() #a few edge cases where  fema shp is 0 km2, but just brushes against model catchment and intersects --> NA fema area. SO remove them.
    
    return(out)
}