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
    fema_elev_shp <- data.frame()
    for(i in states){
        geodatabase <- list.files(paste0('data/path_to_data/CONUS_connectivity_data/CONUS_FEMA/', i), pattern = "\\.gdb$")
        
        #load in fema maps and filter for 100-yr flood polygons
        fema_temp <- sf::st_read(paste0('data/path_to_data/CONUS_connectivity_data/CONUS_FEMA/', i, '/', geodatabase),
                                    layer = 'S_FLD_HAZ_AR',
                                    quiet = TRUE) %>%
            dplyr::filter(FLD_ZONE %in% c('A', 'AE', 'AH', 'AO')) %>% #A zones correspond to 100-yr flood maps, built with a variety of methods. See README.
            dplyr::mutate(STATE = i) %>%
            dplyr::select(c('FLD_AR_ID', 'STATE', 'FLD_ZONE'))
        
        #grab 100yr flood wse
        fema_elev <- sf::st_read(paste0('data/path_to_data/CONUS_connectivity_data/CONUS_FEMA/', i, '/', geodatabase),
                                    layer = 'S_BFE',
                                    quiet = TRUE) %>%
            dplyr::mutate(elev_m = ELEV * 0.3048) %>% #ft to m
            dplyr::select(c('BFE_LN_ID', 'elev_m'))
        
        fema_shp <- rbind(fema_shp, fema_temp)
        fema_elev_shp <- rbind(fema_elev_shp, fema_elev)
    }

    return(list('fema_shp'=fema_shp,
                'fema_elev_shp'=fema_elev_shp))
}







valModelFEMA <- function(huc4id, preppedFEMA, basinAnalysis, basinData){
    sf::sf_use_s2(FALSE)
    
    preppedFEMA_depths <- preppedFEMA$fema_elev_shp
    preppedFEMA <- preppedFEMA$fema_shp
    dem <- terra::unwrap(basinData$terra$dem) #[cm] need to unwrap the packed terra package (for distributed computing)

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
        dplyr::select(c('NHDPlusID', 'GageID', 'Wb_m', 'LengthKM', 'A_qFEMA_m2', 'V_qFEMA_m2', 'Htf_qFEMA_m'))
    
    #only keep model reaches that align with the fema maps
    basinAnalysis <- basinAnalysis %>%
        sf::st_filter(fema_shp, .predicate = sf::st_within) #Only keep model reaches that fall within the fema maps. FEMA maps are inconsistent in space, so we can really only validate in places where we know the entire model reach has been mapped by FEMA (to avoid conflating real false positives with incomplete FEMA maps)

    #keep gage IDs for later
    gageLookup <- basinAnalysis %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(NHDPlusID) %>%
        dplyr::summarise(GageID = dplyr::first(GageID)) #just pass through
    
    #remove catchments with narrow rivers from fema_shps (then calcualate inundated area)   
    unit_catchments <- unit_catchments %>%
        dplyr::filter(NHDPlusID %in% basinAnalysis$NHDPlusID)
    
    basinAnalysis <- basinAnalysis %>% #unit_catchments %>%
    #    sf::st_join(huc8s, join=sf::st_intersects, largest=TRUE) %>% #if in multiple huc8s, assign to the one it's mostly in (this should handle just a few edge cases)
        sf::st_drop_geometry() %>%
        dplyr::group_by(NHDPlusID) %>%
        dplyr::summarise(A_qFEMA_m2 = sum(A_qFEMA_m2, na.rm=T),
                        V_qFEMA_m3 = sum(V_qFEMA_m2, na.rm=T),
                        Htf_qFEMA_m = mean(Htf_qFEMA_m, na.rm=T),
                        channel_area_km2=mean(((Wb_m/1000)*LengthKM)), #mean to pass through groupby
                        LengthKM = mean(LengthKM)) #mean to pass through groupby

    fema_shp <- fema_shp %>%
        sf::st_intersection(unit_catchments) %>%  #make sure fema polygons don't align with this subset of modeled catchments. It's unfair to include FEMA maps outside of the subset of our model with predictions (that overlaps FEMA) and its unfair to include FEMA maps that extend beyond our catchment scheme
        sf::st_transform('+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +type=crs') %>% #equal area projection to match model
        dplyr::mutate(A_femaAEP_m2 = sf::st_area(.)) %>%
      #  sf::st_transform(sf::st_crs(huc8s)) %>% #put back for spatial join
      #  sf::st_join(huc8s, join=sf::st_intersects, largest=TRUE) %>%
        dplyr::group_by(NHDPlusID) %>%
        dplyr::summarise(A_femaAEP_km2 = 1e-6 * sum(as.numeric(A_femaAEP_m2), na.rm=T)) %>%
        dplyr::select(c('NHDPlusID', 'A_femaAEP_km2'))

    #get average fema ground elevation by nhd catchment
    fema_terra <- terra::vect(preppedFEMA_depths)
    ground_evals <- terra::extract(dem, fema_terra, fun=min) #thalweg elevation
    preppedFEMA_depths$groundElevLookup_m <- ground_evals$elev_cm * 0.01 #[cm to m]

    preppedFEMA_depths <- unit_catchments %>%
        sf::st_join(preppedFEMA_depths, join=sf::st_intersects) %>%
        sf::st_drop_geometry() %>%
        dplyr::left_join(basinAnalysis, by='NHDPlusID') %>%
        dplyr::select(c('NHDPlusID', 'elev_m', 'groundElevLookup_m', 'Htf_qFEMA_m'))

    #build validation table
    out <- data.frame('huc4'=huc4id,
                    'NHDPlusID'=basinAnalysis$NHDPlusID,
                    'A_qFEMA_km2'=basinAnalysis$A_qFEMA_m2*1e-6,
                    'V_qFEMA_km3'=basinAnalysis$V_qFEMA_m3*1e-9,
                    'channel_area_km2'=basinAnalysis$channel_area_km2,
                    'LengthKM'=basinAnalysis$LengthKM)

    out <- out %>%
        dplyr::left_join(fema_shp, by='NHDPlusID') %>%
        dplyr::left_join(gageLookup, by='NHDPlusID') %>%
        dplyr::relocate(GageID, .after=NHDPlusID) %>%
        tidyr::drop_na(A_femaAEP_km2) #a few edge cases where  fema shp is 0 km2, but just brushes against model catchment and intersects --> NA fema area. SO remove them.
    
    return(list('for_area'=out,
                'for_depth'=preppedFEMA_depths))
}