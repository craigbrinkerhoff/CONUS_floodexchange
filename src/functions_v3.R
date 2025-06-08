



buildGageFloodFunctions_volumeval <- function(huc4, BHGmodel, gageDF) {
    #dummy dataframe to get rest of function to run
    gageDF <- data.frame(site_no=gageDF$NHDPlusID,
                        id = 'volume_validation',
                        Htf_m=gageDF$depth_m,
                        V_usgs_m3 = gageDF$V_usgs_m3)

    huc2 <- substr(huc4, 1, 2)

    network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        sf::st_zm() %>%
        dplyr::left_join(network_VAA, by='NHDPlusID') %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

    #add reach type
    waterbodies <-  sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                            layer = 'NHDWaterbody',
                            quiet = TRUE) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(c('Permanent_Identifier', 'FType', 'AreaSqKm'))
        colnames(waterbodies) <- c('Permanent_Identifier', 'FType', 'WBAreaSqKm')

    network <- network %>%
        dplyr::left_join(waterbodies, by=c('WBArea_Permanent_Identifier'='Permanent_Identifier')) %>%
        dplyr::mutate(waterbody_type = ifelse(FType %in% c('436', '390') & WBAreaSqKm > 0, 'lake/reservoir', 'river')) %>% #to be a waterbody, must have an area > 0
        dplyr::relocate(WBArea_Permanent_Identifier, .after=waterbody_type) %>%
        dplyr::relocate(WBAreaSqKm, .after=WBArea_Permanent_Identifier) %>%
        dplyr::select(!FType)

    network <- network %>%
        dplyr::filter(NHDPlusID %in% gageDF$site_no)

    #remove impossible flowlines with no drainage area or divergent starting reaches (streamalc == 0; see pg. 45 at https://pubs.usgs.gov/of/2019/1096/ofr20191096.pdf)
    network <- network %>%
        dplyr::filter(AreaSqKm > 0 & TotDASqKm > 0 & StreamCalc > 0) 

    huc4id <- huc4
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)
    basin <- fixGeometries(basin)

    #setup bankfull geometry
    sf::sf_use_s2(FALSE)
    regions <- sf::st_read('data/physio.shp') #physiographic regions
    regions <- fixGeometries(regions)
    shp_temp <- sf::st_join(basin, regions, largest=TRUE) #take the physiographic region that the basin is mostly in (dominant spatial intersection)

    physio_region <- shp_temp$DIVISION

    #add bankfull hydraulics
    network$a_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Wb)})
    network$b_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Wb)})
    network$mean_residual_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$mean_residual_Wb)})
    network$Wb_m <- 10^(network$a_Wb + network$b_Wb*log10(network$TotDASqKm)) * network$mean_residual_Wb #see other notes on bias correction

    #Add flood stages
    network <- network %>%
        dplyr::left_join(gageDF, by=c('NHDPlusID'='site_no'))

    network <- network %>%
        dplyr::filter(AreaSqKm > 0) %>%
        dplyr::relocate(Shape, .after=tidyselect::last_col())

    return(network)
}




buildGageFloodFunctions <- function(huc4, BHGmodel, gageDF, code) {
    #prep gages for filtering
    gageDF <- dplyr::bind_rows(gageDF) %>%
        dplyr::filter(Qexc_m3dy > 0) %>% # integration necessitates storms be 2+ days long to calculate a flux
        dplyr::group_by(site_no) %>%
        dplyr::summarise(id = code,
                        length_dys = mean(length_dy),
                        Htf_m = mean(Htf_m, na.rm=T),
                        Qexc_m3dy = mean(Qexc_m3dy, na.rm=T))

    huc2 <- substr(huc4, 1, 2)

    network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        sf::st_zm() %>%
        dplyr::left_join(network_VAA, by='NHDPlusID') %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

    #add reach type
    waterbodies <-  sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                            layer = 'NHDWaterbody',
                            quiet = TRUE) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(c('Permanent_Identifier', 'FType', 'AreaSqKm'))
        colnames(waterbodies) <- c('Permanent_Identifier', 'FType', 'WBAreaSqKm')

    network <- network %>%
        dplyr::left_join(waterbodies, by=c('WBArea_Permanent_Identifier'='Permanent_Identifier')) %>%
        dplyr::mutate(waterbody_type = ifelse(FType %in% c('436', '390') & WBAreaSqKm > 0, 'lake/reservoir', 'river')) %>% #to be a waterbody, must have an area > 0
        dplyr::relocate(WBArea_Permanent_Identifier, .after=waterbody_type) %>%
        dplyr::relocate(WBAreaSqKm, .after=WBArea_Permanent_Identifier) %>%
        dplyr::select(!FType)

    network <- network %>%
        dplyr::filter(GageID %in% gageDF$site_no)

    #remove impossible flowlines with no drainage area or divergent starting reaches (streamalc == 0; see pg. 45 at https://pubs.usgs.gov/of/2019/1096/ofr20191096.pdf)
    network <- network %>%
        dplyr::filter(AreaSqKm > 0 & TotDASqKm > 0 & StreamCalc > 0) 

    huc4id <- huc4
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)
    basin <- fixGeometries(basin)

    #setup bankfull geometry
    sf::sf_use_s2(FALSE)
    regions <- sf::st_read('data/physio.shp') #physiographic regions
    regions <- fixGeometries(regions)
    shp_temp <- sf::st_join(basin, regions, largest=TRUE) #take the physiographic region that the basin is mostly in (dominant spatial intersection)

    physio_region <- shp_temp$DIVISION

    #add bankfull hydraulics
    network$a_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Wb)})
    network$b_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Wb)})
    network$mean_residual_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$mean_residual_Wb)})
    network$Wb_m <- 10^(network$a_Wb + network$b_Wb*log10(network$TotDASqKm)) * network$mean_residual_Wb #see other notes on bias correction

    network$a_Qb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Qb)})
    network$b_Qb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Qb)})
    network$mean_residual_Qb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$mean_residual_Qb)})
    network$Qb_cms <- 10^(network$a_Qb + network$b_Qb*log10(network$TotDASqKm)) * network$mean_residual_Qb #see other notes on bias correction

    #Add flood stages
    network <- network %>%
        dplyr::left_join(gageDF, by=c('GageID'='site_no'))

    network <- network %>%
        dplyr::filter(AreaSqKm > 0) %>%
        dplyr::relocate(Shape, .after=tidyselect::last_col())

    return(network)
}





runDEMModel <- function(huc4, network){
    #setup basin shapefiles
    huc2 <- substr(huc4, 1, 2)
    dem <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/elev_cm.tif')) #[cm]
    d8 <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/fdr.tif')) #[cm]

    #run embarresingly parallel inundation model
    inundationWrapper <- function(reach) {
        #lower bound on dem resolution (also, these small rivers don't reallllly have floodplains, right?)
        if(reach$Wb_m <= 10 | reach$waterbody_type == 'lake/reservoir'){
            reach$A_m2 <-NA
            reach$V_m3 <- NA
            return(reach)
        }

        #run inundation model
        inundationDataPackage <- prepInundationData(huc4, reach$NHDPlusID, reach$Wb_m, dem, d8)
        if(inundationDataPackage$flag == 'no pour points'){
            reach$A_m2 <-NA
            reach$V_m3 <- NA
            return(reach)
        }

        #loop through flood sizes and model inundation
        areaWrapper <- function(depth){
            floodMap <- modelInundation(inundationDataPackage, depth)

            #convert map to flood area
            flooded_pixels <- floodMap$binary #terra::project(floodMap$binary, "epsg:4269") #NAD83 unprojected, mathce FEMA flood maps
            flooded_pixels[flooded_pixels != 1] <- NA
            flooded_pixels_cellsize <- terra::cellSize(flooded_pixels, mask=TRUE, transform=FALSE) #uses the equal-area projection of the mask (gets same result if transformed to gcs)
            flood_area_m2 <- terra::global(flooded_pixels_cellsize, 'sum', na.rm=T)
            flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

            flood_area_m2 <- flood_area_m2[[1]]

            return(flood_area_m2)
        }

        volWrapper <- function(depth){
            floodMap <- modelInundation(inundationDataPackage, depth)

            #convert map to flood area
            flooded_pixels <- floodMap$binary #terra::project(floodMap$binary, "epsg:4269") #NAD83 unprojected, mathce FEMA flood maps
            flooded_pixels[flooded_pixels != 1] <- NA
            flooded_pixels_cellsize <- terra::cellSize(flooded_pixels, mask=TRUE, transform=FALSE)
            flood_area_m2 <- terra::global(flooded_pixels_cellsize, 'sum', na.rm=T)
            flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

            #convert flood map to volume
            if(flood_area_m2 > 0) { #only if there's any actual inundated area at 10m resolution
                floodmask <- floodMap$binary
                floodmask[floodmask != 1] = 0
                flood_depths <- floodMap$depths * floodmask

                flood_vols <- flood_depths * flooded_pixels_cellsize #m3
                flood_vol_m3 <- terra::global(flood_vols, fun="sum", na.rm=T)

                flood_vol_m3 <- flood_vol_m3[[1]]
            }
            else{
                flood_vol_m3 <- 0
            }

        return(flood_vol_m3)
        }

        reach$A_m2 <- sapply(reach$Htf_m, areaWrapper)
        reach$V_m3 <- sapply(reach$Htf_m, volWrapper)

        return(reach)
    }

    library(plyr)
    library(doParallel)
    registerDoParallel(cores=detectCores()-2)

    network2 <- sf::st_drop_geometry(network)
    network_list <- setNames(split(network2, seq(nrow(network2))), rownames(network2))
    result <- llply(network_list, inundationWrapper, .parallel=TRUE)
    network_fin <- dplyr::bind_rows(result)

    #turn back into shapefile and prep for output
    network <- network %>%
        dplyr::select(c('NHDPlusID', 'Htf_m','Shape')) %>%
        dplyr::left_join(network_fin, by=c('NHDPlusID', 'Htf_m')) %>%
        dplyr::relocate(Shape, .after = tidyselect::last_col())

    return(network)
}



collectValReaches <- function(huc4id){#, usgs_maps){
    sf::sf_use_s2(FALSE)

    #grab maps within our huc4
    usgs_maps <- usgs_maps %>%
        sf::st_transform(crs=sf::st_crs(4269)) %>% #just to make sure
        dplyr::filter(huc4 == huc4id)

    #prep
    huc2 <- substr(huc4id, 1, 2)

    #read in huc4 basin
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)

    unit_catchments <- unit_catchments %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::filter(is.na(GageID) == 0)# %>%
       # sf::st_filter(usgs_maps, .predicate = sf::st_within) #Only keep model reaches that fall within the usgs maps. These maps are inconsistent in space, so we can really only validate in places where we know the entire model reach has been mapped by USGS (to avoid conflating real false positives with incomplete USGS maps)

    #check for basins with no usgs models
    if(nrow(unit_catchments)==0){
        return(data.frame())
    }

    return(unit_catchments$NHDPlusID)
}



wrangleDepthGrids <- function(huc4, reaches_val){
    huc2 <- substr(huc4,1, 2)

    #get nhd
    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        sf::st_zm() %>%
        dplyr::filter(NHDPlusID %in% reaches_val) %>%
        sf::st_transform()
    
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE) %>%
        dplyr::filter(NHDPlusID %in% network$NHDPlusID)
    
    unit_catchments_terra <- terra::vect(unit_catchments)
    
    #get depth grids
    grids <- list.files(paste0('data/path_to_data/CONUS_connectivity_data/USGS_models/depth_grids/',huc4), pattern = "\\.tif$", full.names=TRUE)
    grid_names <- list.files(paste0('data/path_to_data/CONUS_connectivity_data/USGS_models/depth_grids/',huc4), pattern = "\\.tif$", full.names=FALSE)

    #loop trhough and extract mean thalweg depth
    k <- 1
    out <- data.frame()
    for(i in grids){
        grid <- terra::rast(i)
        name <- grid_names[k]

        #get depth at flowline
        depths <- terra::extract(grid, unit_catchments, fun=function(x){mean(x, na.rm=T)})
        flowline_depth_m <- depths[,2]*0.3048 #ft to m (assuming all are ft....)

        #get actual volume
        grid_d <- grid
        res <- terra::res(grid_d)[1]
        grid_d <- terra::project(grid_d, "epsg:4269")
        grid_d <- grid_d * 0.3048 #ft to m
        grid_a <- terra::cellSize(grid_d, unit="m", transform=TRUE)

        grid_v <- grid_d * grid_a #m3
        vols <- terra::extract(grid_v, unit_catchments_terra, fun=function(x){sum(x, na.rm=T)})
        vols <- vols[,2]

        temp <- data.frame('grid'= substr(name, 1, nchar(name)-4),
                        'NHDPlusID'=network$NHDPlusID,
                        'depth_m'=flowline_depth_m)
        temp <- temp %>%
            dplyr::mutate(V_usgs_m3 = vols) %>%
            dplyr::filter(depth_m > 0 & V_usgs_m3 > 0)
        
        out <- rbind(out, temp)

        k <- k + 1
    }
    return(out)
}




getBasinGagesVal <- function(BHGmodel_jacknife){

    gages <- BHGmodel_jacknife$GageID
    gages <- gages[gages != 'ungaged']
    return(gages)
}


# combineTau <- function(gageTau){
#     out <- dplyr::bind_rows(gageTau)
#     return(out)
# }





# prepForTauValidation <- function(gageTau_val, gageTau){
#     df_obs <- dplyr::bind_rows(gageTau_val)
#     df_model <- gageTau

#     colnames(df_obs) <- c('Qprob', 'site_no', 'avg_event_Q_flood_cms_OBS','avg_event_exc_dys_OBS', 'avg_ann_exc_dys_OBS', 'avg_ann_n_events_OBS')

#     #join bankful dataset with gages
#     df <- df_obs %>%
#         dplyr::inner_join(df_model, by=c('site_no', 'Qprob'))

#     return(df)
# }



modelsJacknifeBHG <- function(){
  dataset <- readr::read_csv('data/bhg_us_database_bieger_2015.csv') %>% #available by searching for paper at https://swat.tamu.edu/search
    dplyr::select(c('USGS Station No.','Physiographic Division', '...9', '...11','...13','...17')) #some necessary manual munging for colnames from dataset

  colnames(dataset) <- c('GageID', 'DIVISION', 'DA_km2', 'Qb_cms','Wb_m','Ab_m2')

  dataset$Wb_m <- as.numeric(dataset$Wb_m)
  dataset$Ab_m2 <- as.numeric(dataset$Ab_m2)
  dataset$Qb_cms <- as.numeric(dataset$Qb_cms)
  dataset$DA_km2 <- as.numeric(dataset$DA_km2)

  dataset <- tidyr::drop_na(dataset)

  division <- toupper(sort(unique(dataset$DIVISION))) #make lowercase and sort

  #prep for jacknife, only running the cidivison necessary
  dataset$id <- 1:nrow(dataset)

  out <- data.frame()
  for(i in dataset$id){
    dataset$division_i <- dataset[which(dataset$id == i)[1],]$DIVISION

    #build models, grouped by physiographic region
    models <- dataset %>%
      dplyr::filter(id != i & DIVISION == division_i) %>%
      dplyr::mutate(Wb_log10 = log10(Wb_m),
                    Qb_log10 = log10(Qb_cms),
                    DA_log10 = log10(DA_km2),
                    Ab_log10 = log10(Ab_m2)) %>%
      dplyr::group_by(DIVISION) %>%
      dplyr::do(model_Wb = lm(Wb_log10~DA_log10, data=.),
                model_Qb = lm(Qb_log10~DA_log10, data=.), #fit models by physiographic regions
                model_Ab = lm(Ab_log10~DA_log10, data=.)) %>%
      dplyr::summarise(a_Wb = model_Wb$coef[1], #model intercept
                      b_Wb = model_Wb$coef[2], #model exponent
                      r2_Wb = summary(model_Wb)$r.squared, #model performance
                      mean_residual_Wb = mean(10^(model_Wb$residuals), na.rm=T),
                      see_Wb = summary(model_Wb)$sigma,
                      a_Qb = model_Qb$coef[1],
                      b_Qb = model_Qb$coef[2],
                      r2_Qb = summary(model_Qb)$r.squared,
                      mean_residual_Qb = mean(10^(model_Qb$residuals), na.rm=T),
                      see_Qb = summary(model_Qb)$sigma,
                      a_Ab = model_Ab$coef[1],
                      b_Ab = model_Ab$coef[2],
                      r2_Ab = summary(model_Ab)$r.squared,
                      mean_residual_Ab = mean(10^(model_Ab$residuals), na.rm=T),
                      see_Ab = summary(model_Ab)$sigma)
    
    temp <- data.frame('GageID'=dataset[which(dataset$id == i),]$GageID,
                      'division'=dataset[which(dataset$id == i),]$division_i,
                      'Wb_m_OBS'=dataset[which(dataset$id == i),]$Wb_m,
                      'Ab_m2_OBS'=dataset[which(dataset$id == i),]$Ab_m2,
                      'Qb_cms_OBS'=dataset[which(dataset$id == i),]$Qb_cms)
    temp <- cbind(temp, models)
    out <- rbind(out, temp)
  }

  out <- out %>%
    dplyr::distinct()

  return(out)
}





prepBankfullHydraulics <- function(mode, gageRecord, gage, BHGmodel_jacknife, BHGmodel, minAHGr2, gageRecordStart, gageRecordEnd){
    gageRecord <- dplyr::bind_rows(gageRecord) %>%
        dplyr::select(c('site_no', 'fema100yrflood_cms', 'Q_cms', 'exceed_prob')) #date_hr

    #filter for flows above bankfull, to calculate exceedance probabilities for flood events (aside from the max-annual AEP used to validate against FEMA)
    gage <- sf::st_drop_geometry(gage) %>% 
        dplyr::bind_rows() %>% 
        dplyr::select(c('site_no', 'DA_skm', 'physio_region')) %>% 
        sf::st_drop_geometry() %>%
        dplyr::filter(!(is.na(physio_region))) %>%
        dplyr::distinct()

    BHGmodel_jacknife <- BHGmodel_jacknife %>%
        dplyr::filter(GageID %in% gage$site_no) %>%
        dplyr::select(c('GageID', 'a_Qb', 'b_Qb', 'a_Wb', 'b_Wb', 'a_Ab', 'b_Ab','Wb_m_OBS', 'Ab_m2_OBS', 'Qb_cms_OBS','mean_residual_Qb', 'mean_residual_Wb', 'mean_residual_Ab', 'see_Qb'))
    
    BHGmodel <- BHGmodel %>%
        dplyr::select(c('division', 'a_Qb', 'b_Qb', 'a_Wb', 'b_Wb', 'a_Ab', 'b_Ab', 'mean_residual_Qb', 'mean_residual_Wb', 'mean_residual_Ab', 'see_Qb'))

    if(mode == 'val'){
        gage <- gage %>%
            dplyr::left_join(BHGmodel_jacknife, by=c('site_no'='GageID')) %>%
            tidyr::drop_na()
    }
    else if(mode == 'deploy'){
        gage <- gage %>%
            dplyr::left_join(BHGmodel, by=c('physio_region'='division')) %>%
            tidyr::drop_na()
    }
    else{
        return('error! Please set mode to val or deploy')
    }

    gage$Qb_cms_log10 <- gage$a_Qb + gage$b_Qb*log10(gage$DA_skm)

    gage$Wb_m_log10 <- gage$a_Wb + gage$b_Wb*log10(gage$DA_skm)

    gage$Ab_m2_log10 <- gage$a_Ab + gage$b_Ab*log10(gage$DA_skm)

    if(mode == 'val'){
        gage <- gage %>%
            dplyr::mutate(Qb_cms_log10_low = Qb_cms_log10 - 1.96*see_Qb,
                        Qb_cms_log10_high = Qb_cms_log10 + 1.96*see_Qb) %>%
            dplyr::mutate(Qb_cms = 10^(Qb_cms_log10) * mean_residual_Qb,
                        Qb_cms_low = 10^(Qb_cms_log10_low) * mean_residual_Qb,
                        Qb_cms_high = 10^(Qb_cms_log10_high) * mean_residual_Qb,
                        Wb_m = 10^(Wb_m_log10) * mean_residual_Wb,
                        Ab_m2 = 10^(Ab_m2_log10) * mean_residual_Ab) %>% #using the mean log residual, in transformed to ntural space, a la eq 9.24 from https://pubs.usgs.gov/tm/04/a03/tm4a3.pdf (also see https://nrtwq.usgs.gov/co/methods/)
            dplyr::mutate(Ub_ms = Qb_cms / Ab_m2,
                        Ub_ms_OBS = Qb_cms_OBS / Ab_m2_OBS,
                        Hb_m = Ab_m2 / Wb_m,
                        Hb_m_OBS = Ab_m2_OBS / Wb_m_OBS)
    }

    else if(mode == 'deploy'){
        gage <- gage %>%
            dplyr::mutate(Qb_cms_log10_low = Qb_cms_log10 - 1.96*see_Qb,
                        Qb_cms_log10_high = Qb_cms_log10 + 1.96*see_Qb) %>%
            dplyr::mutate(Qb_cms = 10^(Qb_cms_log10) * mean_residual_Qb,
                        Qb_cms_low = 10^(Qb_cms_log10_low) * mean_residual_Qb,
                        Qb_cms_high = 10^(Qb_cms_log10_high) * mean_residual_Qb,
                        Wb_m = 10^(Wb_m_log10) * mean_residual_Wb,
                        Ab_m2 = 10^(Ab_m2_log10) * mean_residual_Ab) %>% #using the mean log residual, in transformed to ntural space, a la eq 9.24 from https://pubs.usgs.gov/tm/04/a03/tm4a3.pdf (also see https://nrtwq.usgs.gov/co/methods/)
            dplyr::mutate(Ub_ms = Qb_cms / Ab_m2,
                        Hb_m = Ab_m2 / Wb_m)
    }
    else{
        return('error! Please set mode to val or deploy')
    }

    out <- gage %>%
        tidyr::drop_na() %>%
        dplyr::distinct() #in basin 0306 there were duplicate enteries for one gage, so just confirm this isn't happening elsewhere

    return(out)
}








calc_Qexc <- function(mode, gageRecord, gagePrepped, depAHG){
    if(nrow(gageRecord)==0){
        return(data.frame())
    }

    if(mode == 'val'){
        probs <- gagePrepped %>%
            dplyr::filter(site_no == gageRecord[1,]$site_no) %>%
            dplyr::filter(Qb_cms > 0 & Qb_cms_OBS > 0 & Ub_ms > 0 & Ub_ms_OBS > 0 & Wb_m > 0 & Wb_m_OBS > 0)
    }

    if(mode == 'deploy'){
        probs <- gagePrepped %>%
            dplyr::filter(site_no == gageRecord[1,]$site_no) %>%
            dplyr::filter(Qb_cms > 0 & Ub_ms > 0  & Wb_m > 0)
    }

    if(nrow(probs)==0){
        return(data.frame())
    }

    depAHG <- dplyr::bind_rows(depAHG) %>%
        dplyr::filter(site_no == gageRecord[1,]$site_no)

    if(nrow(depAHG)==0){
        return(data.frame())
    }

    gageRecord <- gageRecord %>%
        dplyr::filter(Q_cms > 0) %>%
        dplyr::left_join(depAHG, by='site_no') %>%
        dplyr::mutate(Htf_m_notcorrected = 10^(lm_depth_a + (lm_depth_b * log10(Q_cms)))) %>%
        dplyr::mutate(Htf_m = Htf_m_notcorrected * lm_depth_bias_correct) %>%
        dplyr::filter(Htf_m > 0)

    if(nrow(gageRecord)==0){
        return(data.frame())
    }

    #identify floods
    #validation version
    if(mode == 'val'){
        gageRecord$flood <- ifelse((gageRecord$Q_cms > probs$Qb_cms) & (gageRecord$Q_cms > (probs$Ub_ms*probs$Wb_m*gageRecord$Htf_m)), 'yes', 'no')
        gageRecord$flood_OBS <- ifelse((gageRecord$Q_cms > probs$Qb_cms_OBS) & (gageRecord$Q_cms > (probs$Ub_ms_OBS*probs$Wb_m_OBS*gageRecord$Htf_m)), 'yes', 'no')
        gageRecord$flood_id <- NA
        gageRecord$flood_id_OBS <- NA
        ticker <- 1
        ticker_OBS <- 1
        for(i in 1:nrow(gageRecord)){
            if(gageRecord[i,]$flood == 'yes') {
                gageRecord[i,]$flood_id <- ticker

                if(is.na(gageRecord[i+1,]$site_no) == FALSE & gageRecord[i+1,]$flood == 'no'){
                    ticker <- ticker + 1
                }
            }
            if(gageRecord[i,]$flood_OBS == 'yes') {
                gageRecord[i,]$flood_id_OBS <- ticker_OBS

                if(is.na(gageRecord[i+1,]$site_no) == FALSE & gageRecord[i+1,]$flood_OBS == 'no'){
                    ticker_OBS <- ticker_OBS + 1
                }
            }
        }

        years <- max(lubridate::year(gageRecord$dateTime))-min(lubridate::year(gageRecord$dateTime))

    #https://chemenggcalc.com/residence-time-distribution-in-cstr-and-pfr/
    # E_function <- function(t, hrt){ #assume the floodplain is a well-mixed, homegenous reactor (CTSR/continuous stirred-tank reactor- these are generally treated as having exponential distributions of parcle residence times)
    #     return(1/hrt * exp(-t/hrt))
    # }

        events <- gageRecord %>%
            dplyr::filter(flood_id > 0) %>%
            dplyr::mutate(Ub_ms = probs[1,]$Ub_ms,
                        Wb_m = probs[1,]$Wb_m,
                        dy = dplyr::row_number()) %>%
            dplyr::group_by(flood_id)%>%
            dplyr::summarise(flag='model',
                            n_yrs=years,
                            Qexc_m3dy = (pracma::trapz(dy, Q_cms*86400) - pracma::trapz(dy, (Ub_ms*Wb_m*Htf_m*86400)))/n(), #event floodplain flux (m3/dy)
                            length_dy = n(), #length of event (dys)
                            Htf_m = max(Htf_m)) 

        if(nrow(events)==0){
            return(data.frame())
        }

    # events <- events %>%
    #     dplyr::mutate(event_length_dys = mapply(function(z,x) {integrate(f=function(t) t*E_function(t,z), lower=0, upper=x, subdivisions=1000)$value}, z=hrt_dy, x=length_dy)) %>%#sum((Q_cms*86400) >= event_avg_Q)) %>%
    #     dplyr::mutate(flag = 'model',
    #                 n_yrs = years)

        events_obs <- gageRecord %>%
            dplyr::filter(flood_id_OBS > 0) %>%
            dplyr::mutate(Ub_ms = probs[1,]$Ub_ms_OBS,
                        Wb_m = probs[1,]$Wb_m_OBS,
                        dy = dplyr::row_number()) %>%
            dplyr::group_by(flood_id_OBS)%>%
            dplyr::summarise(flag='obs',
                            n_yrs=years,
                            Qexc_m3dy = (pracma::trapz(dy, Q_cms*86400) - pracma::trapz(dy, (Ub_ms*Wb_m*Htf_m*86400)))/n(), #event floodplain flux (m3/dy)
                            length_dy = n(), #length of event (dys)
                            Htf_m = max(Htf_m)) 

        if(nrow(events_obs)==0){
            return(data.frame())
        }

    # events_obs <- events_obs %>%
    #     dplyr::mutate(event_length_dys = mapply(function(z,x) {integrate(f=function(t) t*E_function(t,z), lower=0, upper=x, subdivisions=1000)$value}, z=hrt_dy, x=length_dy)) %>%#sum((Q_cms*86400) >= event_avg_Q)) %>%
    #     dplyr::mutate(flag = 'obs',
    #                 n_yrs = years)

        colnames(events_obs) <- c('flood_id', 'flag','n_yrs','Qexc_m3dy','length_dy', 'Htf_m')

        out <- rbind(events, events_obs)
        out$site_no <- probs[1,]$site_no
        out$Qb_obs <- probs[1,]$Qb_cms_OBS
    }

    #dploy version
    else if(mode == 'deploy'){
        gageRecord$flood <- ifelse((gageRecord$Q_cms > probs$Qb_cms) & (gageRecord$Q_cms > (probs$Ub_ms*probs$Wb_m*gageRecord$Htf_m)), 'yes', 'no')
        gageRecord$flood_id <- NA
        ticker <- 1

        for(i in 1:nrow(gageRecord)){
            if(gageRecord[i,]$flood == 'yes') {
                gageRecord[i,]$flood_id <- ticker

                if(is.na(gageRecord[i+1,]$site_no) == FALSE & gageRecord[i+1,]$flood == 'no'){
                    ticker <- ticker + 1
                }
            }
        }

        years <- max(lubridate::year(gageRecord$dateTime))-min(lubridate::year(gageRecord$dateTime))

    #https://chemenggcalc.com/residence-time-distribution-in-cstr-and-pfr/
    # E_function <- function(t, hrt){ #assume the floodplain is a well-mixed, homegenous reactor (CTSR/continuous stirred-tank reactor- these are generally treated as having exponential distributions of parcle residence times)
    #     return(1/hrt * exp(-t/hrt))
    # }

        events <- gageRecord %>%
            dplyr::filter(flood_id > 0) %>%
            dplyr::mutate(Ub_ms = probs[1,]$Ub_ms,
                        Wb_m = probs[1,]$Wb_m,
                        dy = dplyr::row_number()) %>%
            dplyr::group_by(flood_id)%>%
            dplyr::summarise(n_yrs=years,
                            Qexc_m3dy = (pracma::trapz(dy, Q_cms*86400) - pracma::trapz(dy, (Ub_ms*Wb_m*Htf_m*86400)))/n(), #event floodplain flux (m3/dy)
                            length_dy = n(), #length of event (dys)
                            Htf_m = max(Htf_m)) #event stage (i.e. crest height) (m))

        if(nrow(events)==0){
            return(data.frame())
        }

    # events <- events %>%
    #     dplyr::mutate(event_length_dys = mapply(function(z,x) {integrate(f=function(t) t*E_function(t,z), lower=0, upper=x, subdivisions=1000)$value}, z=hrt_dy, x=length_dy)) %>%#sum((Q_cms*86400) >= event_avg_Q)) %>%
    #     dplyr::mutate(flag = 'model',
    #                 n_yrs = years)

        out <- events
        out$site_no <- probs[1,]$site_no
    }
    else{
        return('error! Please set mode to val or deploy')
    }

    out <- out %>%
        dplyr::relocate(site_no)

    return(out)
}



prepGage <- function(gageID){
  site <- dataRetrieval::readNWISsite(siteNumbers = gageID) %>%
    dplyr::mutate(lat = dec_lat_va,
                  lon = dec_long_va,
                  DA_skm = drain_area_va * 2.58999) %>% #mi2 to km2
    dplyr::select(c('site_no', 'lat', 'lon', 'DA_skm'))

  #get physiograhic region (for banfull depth model)
  sf::sf_use_s2(FALSE)
  regions <- sf::st_read('data/physio.shp') #physiographic regions
  regions <- fixGeometries(regions)
  site_shp <- tryCatch(sf::st_as_sf(site, coords=c('lon', 'lat'), crs=sf::st_crs(4269)),
                      error = function(m){return(data.frame())})
  if(nrow(site_shp)==0) {return(data.frame())} #if no gage data, move on
  
  shp_temp <- sf::st_join(site_shp, regions) #take the physiographic region that the basin is mostly in (dominant spatial intersection)
  site_shp$physio_region <- shp_temp$DIVISION
  
  site <- site_shp %>%
    sf::st_drop_geometry()

  return(site)
}





prepFlowRecord <- function(gage, gageRecordStart, gageRecordEnd, minRecordLength){
  #prep
  gageID <- gage$site_no
  
  #grab sub-daily discharge data
  gagedata <- tryCatch(dataRetrieval::readNWISdata(sites = gageID,
                                          service='dv',
                                          parameterCd = '00060', #discharge parameter code [cfs]
                                          startDate = gageRecordStart,
                                          endDate = gageRecordEnd),
                      error = function(m){return(data.frame())})

  if(nrow(gagedata)==0) {return(data.frame())} #if no gage data, move on
  if(!("X_00060_00003" %in% colnames(gagedata))){return(data.frame())} #if incorrect discharge column name (happens sometimes), move on

  #convert to metric
  gagedata$Q_cms <- gagedata$X_00060_00003 * 0.0283 #cfs to cms
  
  #only keep flowing events
  gagedata <- gagedata %>%
    dplyr::filter(Q_cms > 0)

  # convert to date (workaround to handle known date class bug with midnight- https://github.com/tidyverse/lubridate/issues/1124)
  gagedata$date <- lubridate::ymd_hms(format(as.POSIXct(gagedata$dateTime), format = "%Y-%m-%d %T %Z"))

  # gagedata_og <- gagedata

  #downsample to hourly (just to be consistent and reduce volume of data)
  # gagedata <- gagedata_og %>%
  #   dplyr::mutate(date_hr = lubridate::ymd_h(paste0(lubridate::year(date), '-', lubridate::month(date),'-',lubridate::day(date), '-', lubridate::hour(date)))) %>%
  #   dplyr::group_by(date_hr) %>%
  #   dplyr::summarise(site_no = dplyr::first(site_no), #pass through the group by
  #                   Q_cms = mean(Q_cms, na.rm=T))

  #get max annual floods for calculating the FEMA 100yr AEP
  gagedata_aep <- gagedata %>%
    dplyr::mutate(year = lubridate::year(date)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(site_no = dplyr::first(site_no), #pass through the group by
                  Q_cms = max(Q_cms, na.rm=T))

  gagedata_aep$rank <- rank(-gagedata_aep$Q_cms) #annual flood data
  gagedata_aep$exceed_prob <- (gagedata_aep$rank/(nrow(gagedata_aep) + 1))

  #get number of years on recordf
  start <- min(gagedata_aep$year)
  end <- max(gagedata_aep$year)
  spread <- end - start

  if(spread < minRecordLength){return(data.frame())}
  if(nrow(gagedata)== 0){return(data.frame())}
  
  gagedata$spread <- spread


  #get exceedance probs for each flow value
  gagedata$rank <- rank(-gagedata$Q_cms) #hourly flow data
  gagedata$exceed_prob <- (gagedata$rank/(nrow(gagedata) + 1)) #exceed prob for a given hour over 'spread' years

  #grab rating table to convert to stage
  #source('src/utils.R')
  ratingTable <- tryCatch(dataRetrieval::readNWISrating(gageID, type='exsa'), #readNWISrating_CRAIG
                      error = function(m){return(data.frame())})
  ratingTable$stage_m <- (ratingTable$INDEP  + ratingTable$SHIFT) * 0.3048 #ft to m
  ratingTable$Q_cms <- ratingTable$DEP * 0.0283

  #convert historical streamflow record to stages using rating table
  gagedata$stage_m <- sapply(gagedata$Q_cms, function(i){ratingTable[which.min(abs(ratingTable$Q_cms - i)),]$stage_m})

  if(is.list(gagedata$stage_m)){return(data.frame())} #sometimes you get an empty list b/c the rating table is empty

  #add flag if gagedata is outside of the rating table
  gagedata$ratingflag <- ifelse(gagedata$Q_cms > max(ratingTable$Q_cms) | gagedata$Q_cms < min(ratingTable$Q_cms), 1, 0)

  #add fema 100 yr aep
  fema100yrflood_cms <- gagedata_aep[which.min(abs(gagedata_aep$exceed_prob - 0.01)),]$Q_cms

  gagedata$fema100yrflood_cms <- fema100yrflood_cms

  return(gagedata)
}




prepInundationData <- function(huc4, reachID, bankfullWidth, dem, d8){
  #prep
  huc2 <- substr(huc4, 1, 2)
  
  sql <- paste0('SELECT * FROM NHDPlusCatchment WHERE NHDPlusID = ', reachID)
  unit_catchment <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    query = sql, #sql query to only load the catchment we want efficiently (see above)
                                    quiet = TRUE)
  
  unit_catchment <- terra::vect(unit_catchment)
  unit_catchment <- terra::project(unit_catchment, terra::crs(dem))

  #read in the nhd centerline
  sql <- paste0('SELECT * FROM NHDFlowline WHERE NHDPlusID = ', reachID)
  flowline <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                    layer = 'NHDFlowline',
                                    query = sql, #sql query to only load the catchment we want efficiently (see above)
                                    quiet = TRUE) %>%
    dplyr::filter(NHDPlusID == reachID) %>%
    sf::st_transform(crs=sf::st_crs(terra::crs(dem))) %>%
    sf::st_zm()

  flowline_terra <- terra::vect(flowline)

  #math the nhdplus unit catchment 
  dem_clipped <- terra::crop(dem, unit_catchment)
  dem_clipped <- terra::mask(dem_clipped, unit_catchment)
  
  #fill and breach the dem to calculate flow directions
  r <- terra::rast(flowline_terra, resolution=10, extent=terra::ext(dem_clipped))
  flowline_rast <- terra::rasterize(flowline_terra, r)
  flowline_rast[flowline_rast == 1] = -5000 #50m following https://doi.org/10.1016/j.rse.2024.114333
  flowline_rast[is.na(flowline_rast)] = 0
  flowline_rast <- terra::mask(flowline_rast, dem_clipped)

  dem_burned <- dem_clipped + flowline_rast
  dem_burned <- flowdem::breach(dem_burned) #https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.10648
  hydrodem <- flowdem::fill(dem_burned, epsilon=TRUE) #https://www.sciencedirect.com/science/article/pii/S0098300413001337
  d8_clipped <- terra::terrain(hydrodem, 'flowdir')#flowdem::dirs(hydrodem, mode='d8')

  #snap pour point to the largest flow accumulation cell (to avoid incorrect snapping based on nhd topology)
  flowaccum <- terra::flowAccumulation(d8_clipped) #flowdem::accum(d8_clipped, mode='d8')
  flowaccum <- terra::mask(flowaccum, unit_catchment)

  #extract dem
  dem_clipped <- dem_clipped*0.01 #cm to m
  
  #use HAND method to get relative elevations for inundation mapping
  #https://nhess.copernicus.org/articles/19/2405/2019/nhess-19-2405-2019.html
  #https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022WR032039
  #https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12661
  #get thalweg mask
  thalweg_pour_points <- terra::extract(flowaccum, flowline_terra, xy=TRUE)

  #get pixel unit catchments and HAND (in parallel)
  #8/14/2024: NOTE: I know there are teeny tiny catchments getting wrapped up in here because the
  #nhd flowline technically intersects them just a bit, but they aren't the actual accumulated flow river pixels.
  #I think it's ok? In theory those catchments should be draining into other pixels but the watersheds get recalculated
  #for every pixel, so as flow accum goes up, it should just include these erroneous catchments, which were run way at the beginning
  #(since it's sorted by DA). Meaning these little catchments are wrapped into the correct per-pixel unit catchments, and when they are taken
  #care of at the top, they don't overlap the downstream catchments so they don't influence the per-pixel unit catchment calculation.
  # I think this makes sense, but worth another think through.
  pixel_pour_points <- na.omit(thalweg_pour_points[order(thalweg_pour_points$flowdir),]) #terra doesn't rename the bands, this is actually flow accum from previous function

  #first, ensure the pour point is on the dem... (there are some edge cases where dem and flowline don't align perfectly). In these cases FOR NOW, i"ll just pass an error term
  if(nrow(pixel_pour_points) == 0){ #terra doesn't rename the bands, this is actually flow accum from previous function
    return(list('hand'=NA,
              'bankfullMask'=NA,
              'flag'='no pour points'))
  }
  
  for(i in 1:nrow(pixel_pour_points)){
    pixel_pour_point <- terra::vect(pixel_pour_points[i,], c('x','y'))
    pixel_watershed <- terra::watershed(x=d8_clipped, pourpoint=cbind(pixel_pour_points[i,]$x, pixel_pour_points[i,]$y)) #can't pass a spatVector in dev version... needs to be list of point coords
    pixel_watershed[pixel_watershed == 0] <- NA
    pixel_dem <- terra::mask(dem_clipped, pixel_watershed)
    
    #if most upstream catchment, there's nothing upstream to remove
    if(i == 1){
      #get HAND
      elev <- terra::extract(pixel_dem, pixel_pour_point)
      rel_elev <- pixel_dem - elev[1,]$elev_cm #band wasn't renamed, this is already in m
      next
    }
    
    #remove upstream pixel_watersheds (so this is just the unit pixel watershed)
    pixel_watershed_up <- terra::watershed(x=d8_clipped, pourpoint=cbind(pixel_pour_points[i-1,]$x, pixel_pour_points[i-1,]$y))
    pixel_watershed_up[pixel_watershed_up == 0] <- NA
    pixel_mask <- terra::app(terra::sds(pixel_watershed, pixel_watershed_up), fun="sum", na.rm=T)
    pixel_mask[pixel_mask == 2] <- NA
    
    #get pixel per-pixel unit catchment dem
    pixel_dem <- terra::mask(pixel_dem, pixel_mask)
    
    #get HAND
    elev <- terra::extract(pixel_dem, pixel_pour_point)
    rel_elev_temp <- pixel_dem - elev[1,]$elev_cm #band wasn't renamed, this is already in m
    rel_elev <- terra::merge(rel_elev, rel_elev_temp)
  }
  
  hand <- rel_elev
  hand[hand < 0] = 0
  
  #return what's needed for inundation modeling
  return(list('hand'=hand,
              'flag'='has pour points'))
}




modelInundation <- function(inundationData, floodDepth){
  library(terra)
  hand <- inundationData$hand
  
  #actual inundation calculation
  floodDepths <- floodDepth - hand
  floodMap <- floodDepths
  floodMap[floodMap < 0] <- 0
  floodMap[floodMap > 0] <- 1
  
  floodDepths[floodDepths < 0] <- 0 # if depth is negative (artifact of models and hand), just set it to 0

  return(list('binary'=floodMap,
              'depths'=floodDepths))
}





buildDepthAHG <- function(gage, minADCPMeas){
  if(nrow(gage)==0){return(data.frame())}
  #prep
  gageID <- gage$site_no

  #prep custom rating table from available adcp measurements
  surfaceData <- tryCatch(dataRetrieval::readNWISmeas(gageID, expanded = TRUE),
                    error=function(m){return(data.frame('site_no'=gageID,
                                                      'lm_r2'=NA,
                                                      'lm_depth_a'=NA,
                                                      'lm_depth_b'=NA))})
  
  if(nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
                                                      'lm_r2'=NA,
                                                      'lm_depth_a'=NA,
                                                      'lm_depth_b'=NA))}
  surfaceData <- surfaceData %>%
    dplyr::filter(measured_rating_diff %in% c('Fair','Good', 'Excellent')) %>% #Fair are < 8% flow diff, Good are < 5% flow diff, and Excellent are <2% flow diff
    dplyr::mutate('Q_cms'=chan_discharge* 0.0283,
                  'width_m'=chan_width*0.3048,
                  'area_m2'=chan_area*0.092903,
                  'velocity_ms'=chan_velocity*0.3048) %>%
    dplyr::select(c('site_no', 'Q_cms', 'width_m', 'area_m2', 'velocity_ms')) %>%
    dplyr::mutate(depth_m=area_m2 / width_m) %>%
    dplyr::filter(depth_m > 0 & Q_cms > 0 & is.finite(depth_m) & is.finite(Q_cms))

  #remove likely erronous outlier depth measurements
  iqr <- IQR(surfaceData$depth_m, na.rm=T)
  
  surfaceData <- surfaceData %>%
    dplyr::filter(depth_m < (quantile(depth_m, 0.75, na.rm=T) + 1.5*iqr)) %>%
    dplyr::filter(depth_m > (quantile(depth_m, 0.25, na.rm=T) - 1.5*iqr))
  
  if(nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
                                                      'lm_r2'=NA,
                                                      'lm_depth_a'=NA,
                                                      'lm_depth_b'=NA))}
  
  #get bankfull depth
  lm_depth <- lm(log10(depth_m)~log10(Q_cms), data=surfaceData)
  lm_depth_a <- coef(lm_depth)[1]
  lm_depth_b <- coef(lm_depth)[2]
  lm_depth_bias <- mean(10^(lm_depth$residuals), na.rm=T)
  
  #setup site info for later export
  site_info <- data.frame('site_no'=gageID,
                          'lm_r2'=summary(lm_depth)$r.squared,
                          'lm_depth_a'=lm_depth_a,
                          'lm_depth_b'=lm_depth_b,
                          'lm_depth_bias_correct'=lm_depth_bias)

  rownames(site_info) <- NULL
  
  return(site_info)
}




# hortonScaling <- function(basinAnalysis, huc4){

#   scaledBasin <- basinAnalysis %>%
#     sf::st_drop_geometry() %>%
#     dplyr::group_by(StreamCalc) %>% 
#     dplyr::summarise(Af_by_order = sum(A_q1_m2 - (Wb_m*LengthKM*1000), na.rm=T),
#                     A_by_order = sum(A_q1_m2, na.rm=T),
#                     Vf_by_order = sum(V_q1_m3 - (Wb_m*LengthKM*1000*Htf_q1_m), na.rm=T),
#                     V_by_order = sum(V_q1_m3, na.rm=T),
#                     frac = round((sum(Wb_m > 10, na.rm=T)/n()),2),
#                     lengthKM_by_order = sum(LengthKM, na.rm=T)) %>%
#     dplyr::mutate(frac2 = ifelse(frac < 1, NA, frac))

#   #fit Horton equations for floodplain water (only using the orders with 100% model coverage (i.e. Wb > 10m))
#   area_model <- lm(log(Af_by_order*frac2)~log(StreamCalc), data=scaledBasin)
#   vol_model <- lm(log(Vf_by_order*frac2)~log(StreamCalc), data=scaledBasin)
#   vol_model_all <- lm(log(V_by_order*frac2)~log(StreamCalc), data=scaledBasin)

#   #also fit length
#   length_model <- lm(log(lengthKM_by_order)~log(StreamCalc), data=scaledBasin)

#   #use horton eqs to calculate inundated areas and volumes for orders with frac < 100%
#   scaledBasin$Af_by_order_km2_fin <- exp(predict(area_model, scaledBasin)) * 1e-6 #[m2 to km2]
#   scaledBasin$Vf_by_order_km3_fin <- exp(predict(vol_model, scaledBasin)) * 1e-9 #[m3 to km3]
#   scaledBasin$V_by_order_km3_fin <- exp(predict(vol_model_all, scaledBasin)) * 1e-9 #[m3 to km3]

#   #prep output
#   scaledBasin <- scaledBasin %>%
#     dplyr::mutate('huc4'=huc4,
#                 'Hf_avg_by_order_cm_fin'=(Vf_by_order_km3_fin / Af_by_order_km2_fin / (Af_by_order_km2_fin*1e6/100))*100000, #100m2 cell size
#                 'area_model_r2'=summary(area_model)$r.squared,
#                 'vol_model_r2'=summary(vol_model)$r.squared,
#                 'vol_all_model_r2'=summary(vol_model_all)$r.squared) %>%
#     dplyr::select(c('huc4', 'StreamCalc', 'frac', 'area_model_r2', 'vol_model_r2', 'vol_all_model_r2', 'Af_by_order_km2_fin', 'Vf_by_order_km3_fin', 'V_by_order_km3_fin', 'Hf_avg_by_order_cm_fin'))

#   return(scaledBasin)
# }




# regulateFlooding <- function(basinAnalysis, barriers, area_thresh_perc, huc4id, snapping_thresh){
#   huc2 <- substr(huc4id, 1, 2)

#   network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE) %>%
#     dplyr::select(c('NHDPlusID', 'TotDASqKm', 'ToNode', 'FromNode'))

#   basins <- basinAnalysis

#   basins_all <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
#     sf::st_zm() %>%
#     dplyr::left_join(network_VAA, by='NHDPlusID')

#   basins <- basins %>%
#     dplyr::select(!'TotDASqKm') %>%
#     dplyr::left_join(network_VAA, by='NHDPlusID')

#   huc_shp <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
#      dplyr::filter(huc4 == huc4id)

#   #join Global Dam Watch Database to nhd
#   barriers <- barriers %>% #sf::st_read('data/path_to_data/CONUS_sediment_data/GDW_barriers_v1_0.shp') %>%
#     sf::st_transform(crs=sf::st_crs(basins_all)) %>%
#     sf::st_crop(huc_shp)

#   buffered_barriers <- barriers %>%
#     sf::st_buffer(snapping_thresh)
  
#   barriers_init <- buffered_barriers %>%
#     sf::st_join(basins_all) %>%
#     sf::st_drop_geometry() %>%
#     dplyr::filter(abs(RESV_CATCH_SKM-TotDASqKm)/RESV_CATCH_SKM <= area_thresh_perc & is.na(RESV_CATCH_SKM)==0) %>%
#     dplyr::select(c('HYRIV_ID', 'NHDPlusID', 'RESV_CATCH_SKM'))
  
#   for_dist <- barriers %>%
#     dplyr::left_join(barriers_init, by='HYRIV_ID') %>%
#     dplyr::filter(is.na(NHDPlusID)==0)

#   #loop through dams and find the nhd reach closest to the barrier (after filtering for only reaches within 1km and within 5% drainage area agreeement)
#   for_dist2 <- for_dist %>%
#     dplyr::group_by(HYRIV_ID) %>%
#     dplyr::summarise(n=n())
  
#   out <- NA
#   for(i in for_dist2$HYRIV_ID){
#       ids <- for_dist[for_dist$HYRIV_ID == i,]$NHDPlusID
#       check_reaches <- basins_all[basins_all$NHDPlusID %in% ids,]

#       nearest_reach <- check_reaches[sf::st_nearest_feature(for_dist2[for_dist2$HYRIV_ID == i,], check_reaches),]
#       out <- c(out, nearest_reach$NHDPlusID)
#   }

#   barriers_fin <- barriers_init %>%
#     dplyr::filter(NHDPlusID %in% out) %>%
#     dplyr::select(c('NHDPlusID', 'RESV_CATCH_SKM'))
    
    
#     # dplyr::group_by(NHDPlusID) %>% 
#     # #dplyr::slice_min(abs(RESV_CATCH_SKM-TotDASqKm)/RESV_CATCH_SKM) %>%
#     # dplyr::ungroup() %>%
#     # dplyr::select(c('NHDPlusID', 'RESV_CATCH_SKM'))

#   #runs asynchrounsouly for reach i using doparallel
#   effContrbWrapper <- function(basin, basins_all){
#     #crawl upstream to find reservoir drainage areas
#     ids_vec <- basin$FromNode
#     lookup <- data.frame()
#     while(length(ids_vec)> 0){
#       up_basins <- basins_all %>%
#         sf::st_drop_geometry() %>%
#         dplyr::filter(ToNode %in% ids_vec) %>%
#         dplyr::left_join(barriers_fin, by='NHDPlusID') %>%
#         dplyr::group_by(NHDPlusID) %>% #can have multiple dams on a reach, so just keep the most downstream one (see next line)
#         dplyr::slice_max(TotDASqKm) %>%
#         dplyr::ungroup() %>%
#         dplyr::mutate(RESV_CATCH_SKM = ifelse(RESV_CATCH_SKM > 0, TotDASqKm, NA)) #for consistency, use the nhd drainage area

#       #if a reservoir is found, store that catchment area away and stop crawling that part of the network
#       if(sum(up_basins$RESV_CATCH_SKM, na.rm=T) > 0){
#           temp <- up_basins %>%
#             dplyr::filter(RESV_CATCH_SKM > 0 & is.na(RESV_CATCH_SKM)==0) %>%
#             dplyr::select(c('NHDPlusID', 'FromNode', 'RESV_CATCH_SKM'))
#           lookup <- rbind(lookup, temp)

#           #remove that basin from the list of crawling basins
#           up_basins <- up_basins %>%
#             dplyr::filter(!(NHDPlusID %in% temp$NHDPlusID))
#         }
#       ids_vec <- up_basins$FromNode
#     }
#     lookup <- lookup[!duplicated(lookup),]
#   print(lookup)
#     basin$regulatedDA_skm <- sum(lookup$RESV_CATCH_SKM, na.rm=T)

#     return(basin)
#   }

#     library(plyr)
#     library(doParallel)
#     registerDoParallel(cores=detectCores()-2)

#     basin_list <- setNames(split(basins, seq(nrow(basins))), rownames(basins))
#     result <- llply(basin_list, effContrbWrapper, basins_all, .parallel=TRUE)
#     basins_fin <- dplyr::bind_rows(result)

#   return(basins_fin)
# }
