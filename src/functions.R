## Main functions
## Craig Brinkerhoff
## Spring 2026


#' buildGageFloodFunctions
#'
#' Build df with gauge Qflood
#' Add all auxillary variables needed, i.e. bankful hydraulics
#'
#' @param HUC4 basin ID code
#' @param BHGmodel Bankfull hydraulics models
#' @param gageDF gauge metadata dataframe
#' @param min_floods Minimum number of floods on record we allow (see _targets.R)
#'
#' @return dataframe of gauges with aux info for computing Q_flood
buildGageFloodFunctions <- function(huc4, BHGmodel, gageDF, min_floods) {
    if(nrow(gageDF %>% dplyr::bind_rows())==0){
        return(data.frame())
    }

    #prep gauges for filtering
    gageDF_02 <- dplyr::bind_rows(gageDF) %>%
        dplyr::filter(Qexc_m3dy > 0) %>% # note that integration necessitates storms be 2+ days long to calculate a flux. Fine for the large floods we focus on.
        dplyr::mutate(prob = 0.02)
    
    gageDF_10 <- dplyr::bind_rows(gageDF) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.10)

    gageDF_50 <- dplyr::bind_rows(gageDF) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.50)
    
    gageDF_20 <- dplyr::bind_rows(gageDF) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.20)

    gageDF <- rbind(gageDF_02, gageDF_10, gageDF_50, gageDF_20)
    gageDF <- gageDF %>%
        dplyr::group_by(site_no, prob)%>%
        dplyr::summarise(Htf_m = quantile(Htf_m, 1-prob, na.rm=T),
                        Qexc_m3dy = quantile(Qexc_m3dy, 1-prob, na.rm=T),
                        Q_m3dy = quantile(Q_m3dy, 1-prob, na.rm=T),
                        n_floods=n()) %>%
        dplyr::filter(n_floods >= min_floods) %>%
        dplyr::distinct()

    #join hydrography attributes
    huc2 <- substr(huc4, 1, 2)
    network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        sf::st_zm() %>%
        dplyr::left_join(network_VAA, by='NHDPlusID') %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

    #load Indiana fixed basins (see manuscript)
    indiana <- readr::read_csv('data/indiana.csv')
    if(huc4 %in% indiana$huc4){
        network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/indiana/indiana_fixed_',huc4,'.shp'),quiet = TRUE) %>%
            dplyr::select(!path) %>%
            sf::st_zm() %>%
            dplyr::left_join(network_VAA, by='NHDPlusID') %>%
            dplyr::left_join(network_gages, by='NHDPlusID') %>%
            dplyr::mutate(WBArea_Permanent_Identifier = WBArea_Per)

        sf::st_geometry(network) <- "Shape"

        network <- network %>%
            dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))
    }

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

    #keep gauged reaches only
    network <- network %>%
        dplyr::filter(GageID %in% gageDF$site_no)

    #remove impossible flowlines with no drainage area or divergent starting reaches (streamalc == 0; see pg. 45 at https://pubs.usgs.gov/of/2019/1096/ofr20191096.pdf)
    network <- network %>%
        dplyr::filter(AreaSqKm > 0 & TotDASqKm > 0 & StreamCalc > 0)
    
    if(nrow(network)==0){
        return(data.frame())
    }

    #fix geometeries (if necessary)
    huc4id <- huc4
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)
    basin <- fixGeometries(basin)

    #setup bankfull geometry
    sf::sf_use_s2(FALSE)
    regions <- sf::st_read('data/physio.shp') #physiographic regions
    regions <- fixGeometries(regions)
    shp_temp <- sf::st_join(basin, regions, largest=TRUE) #take the physiographic region that the basin is primarily overlapping (dominant spatial intersection)

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
        dplyr::mutate(huc4=huc4) %>%
        dplyr::relocate(Shape, .after=tidyselect::last_col()) %>%
        dplyr::relocate(huc4)

    return(network)
}


#' buildGageFloodFunctions_volumeval
#'
#' Build df with gauge Qflood for validation sites (gauges with in situ Qb values)
#' Add all auxiliary variables needed, i.e. bankfull hydraulics
#'
#' @param HUC4 basin ID code
#' @param BHGmodel Bankfull hydraulics models
#' @param gageDF gauge metadata dataframe
#'
#' @return dataframe of gauges with aux info for computing Q_flood, filtered for those with bankfull observations
buildGageFloodFunctions_volumeval <- function(huc4, BHGmodel, gageDF) {
    if(nrow(gageDF)==0){
        return(data.frame())
    }

    #dummy dataframe to get rest of function to run
    gageDF <- data.frame(site_no=gageDF$NHDPlusID,
                        id = 'volume_validation',
                        Htf_m=gageDF$depth_m,
                        V_usgs_m3 = gageDF$V_usgs_m3)

    #load hydrography attributes
    huc2 <- substr(huc4, 1, 2)
    network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        sf::st_zm() %>%
        dplyr::left_join(network_VAA, by='NHDPlusID') %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

    #load Indiana fixed basins (see manuscript)
    indiana <- readr::read_csv('data/indiana.csv')
    if(huc4 %in% indiana$huc4){
        network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/indiana/indiana_fixed_',huc4,'.shp'),quiet = TRUE) %>%
            dplyr::select(!path) %>%
            sf::st_zm() %>%
            dplyr::left_join(network_VAA, by='NHDPlusID') %>%
            dplyr::left_join(network_gages, by='NHDPlusID') %>%
            dplyr::mutate(WBArea_Permanent_Identifier = WBArea_Per)

        sf::st_geometry(network) <- "Shape"

        network <- network %>%
            dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

    }

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

    #onlyy keep gauged reaches
    network <- network %>%
        dplyr::filter(NHDPlusID %in% gageDF$site_no)

    #remove impossible flowlines with no drainage area or divergent starting reaches (streamalc == 0; see pg. 45 at https://pubs.usgs.gov/of/2019/1096/ofr20191096.pdf)
    network <- network %>%
        dplyr::filter(AreaSqKm > 0 & TotDASqKm > 0 & StreamCalc > 0) 

    #fix geometries (if applicable)
    huc4id <- huc4
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)
    basin <- fixGeometries(basin)

    #setup bankfull geometry via physiographic regions
    sf::sf_use_s2(FALSE)
    regions <- sf::st_read('data/physio.shp')
    regions <- fixGeometries(regions)
    shp_temp <- sf::st_join(basin, regions, largest=TRUE) #take the physiographic region that the basin is primarily in (dominant spatial intersection)

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
        dplyr::mutate(huc4=huc4) %>%
        dplyr::relocate(Shape, .after=tidyselect::last_col()) %>%
        dplyr::relocate(huc4)

    return(network)
}


#' runDEMModel
#'
#' Run inundation model and add Vflood to gauge dataframe
#'
#' @param huc4 basin ID code
#' @param network hydrography dataframe for inundation mapping
#'
#' @return hydrography with inundation estimates
runDEMModel <- function(huc4, network){
    if(nrow(network)==0){
        return(data.frame())
    }

    #setup basin shapefiles
    huc2 <- substr(huc4, 1, 2)
    dem <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/elev_cm.tif')) #[cm]
    d8 <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/fdr.tif')) #[cm]

    #wrapper for running inundation model in parallel
    inundationWrapper <- function(reach) {
        #lower bound on dem resolution
        if(reach$Wb_m <= 10 | reach$waterbody_type == 'lake/reservoir'){
            reach$A_m2 <-NA
            reach$V_m3 <- NA
            return(reach)
        }

        #run inundation model
        inundationDataPackage <- prepInundationData(huc4, reach$NHDPlusID, reach$Wb_m, dem, d8) #see ~src/functions.R
        if(inundationDataPackage$flag == 'no pour points'){
            reach$A_m2 <-NA
            reach$V_m3 <- NA
            return(reach)
        }

        #loop through flood sizes and model inundation
        areaWrapper <- function(depth){
            floodMap <- modelInundation(inundationDataPackage, depth)

            #convert map to flood area
            flooded_pixels <- floodMap$binary
            flooded_pixels[flooded_pixels != 1] <- NA
            flooded_pixels_cellsize <- terra::cellSize(flooded_pixels, mask=TRUE, transform=FALSE) #uses conus equal-area projection of the mask
            flood_area_m2 <- terra::global(flooded_pixels_cellsize, 'sum', na.rm=T)
            flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

            flood_area_m2 <- flood_area_m2[[1]]

            return(flood_area_m2)
        }

        volWrapper <- function(depth){
            floodMap <- modelInundation(inundationDataPackage, depth)

            #convert map to flood area
            flooded_pixels <- floodMap$binary
            flooded_pixels[flooded_pixels != 1] <- NA
            flooded_pixels_cellsize <- terra::cellSize(flooded_pixels, mask=TRUE, transform=FALSE)
            flood_area_m2 <- terra::global(flooded_pixels_cellsize, 'sum', na.rm=T)
            flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

            #convert flood map to volume
            if(flood_area_m2 > 0) {
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

    #run inundation in parallel (within a single machine)
    library(plyr)
    library(doParallel)
    registerDoParallel(cores=detectCores()-2)

    #run in parallel
    network2 <- sf::st_drop_geometry(network)
    network_list <- setNames(split(network2, seq(nrow(network2))), rownames(network2))
    result <- llply(network_list, inundationWrapper, .parallel=TRUE)
    network_fin <- dplyr::bind_rows(result)

    #turn back into shapefile and prep for output
    network <- network %>%
        sf::st_drop_geometry() %>%
        dplyr::select(c('NHDPlusID', 'Htf_m')) %>%
        dplyr::left_join(network_fin, by=c('NHDPlusID', 'Htf_m')) %>%
        dplyr::mutate(huc4 = huc4) %>%
        dplyr::relocate(huc4)

    return(network)
}


#' addOtherNHDFeatures
#'
#' Add NHD features for ML model training on gauges
#'
#' @param network hydrography network dataframe
#' @param huc4id basin ID code
#'
#' @return hydrography dataframe with added features
addOtherNHDFeatures <- function(network, huc4id){
    if(nrow(network)==0){
        return(data.frame())
    }

    huc2 <- substr(huc4id,1,2)

    #load precip and temperature datasets from NHD
    flow_df <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusEROMMA', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'QBMA')) #naturalized flow

    precip_01 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM01', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM01')) #precip jan

    precip_02 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM02', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM02')) #precip feb
    
    precip_03 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM03', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM03')) #precip mar
    
    precip_04 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM04', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM04')) #precip apr
    
    precip_05 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM05', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM05')) #precip may
    
    precip_06 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM06', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM06')) #precip jun
    
    precip_07 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM07', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM07')) #precip jul
    
    precip_08 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM08', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM08')) #precip aug
    
    precip_09 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM09', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM09')) #precip sep
    
    precip_10 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM10', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM10')) #precip oct
    
    precip_11 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM11', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM11')) #precip nov
    
    precip_12 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM12', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM12')) #precip dec

    temp_01 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM01', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM01')) #temp jan

    temp_02 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM02', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM02')) #temp feb
    
    temp_03 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM03', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM03')) #temp mar
    
    temp_04 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM04', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM04')) #temp apr
    
    temp_05 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM05', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM05')) #temp may
    
    temp_06 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM06', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM06')) #temp jun
    
    temp_07 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM07', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM07')) #temp jul
    
    temp_08 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM08', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempVMM08')) #temp aug
    
    temp_09 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM09', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM09')) #temp sep
    
    temp_10 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM10', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM10')) #temp oct

    temp_11 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM11', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM11')) #temp nov
    
    temp_12 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM12', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM12')) #temp dec
    
    #join everything into single data frame
    network <- network %>%
        dplyr::filter(waterbody_type == 'river') %>%
        dplyr::select(c('NHDPlusID', 'GageID', 'StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'prob','Q_m3dy', 'Qexc_m3dy', 'V_m3')) %>%
        dplyr::left_join(flow_df, by='NHDPlusID') %>%
        dplyr::left_join(precip_01, by='NHDPlusID') %>%
        dplyr::left_join(precip_02, by='NHDPlusID') %>%
        dplyr::left_join(precip_03, by='NHDPlusID') %>%
        dplyr::left_join(precip_04, by='NHDPlusID') %>%
        dplyr::left_join(precip_05, by='NHDPlusID') %>%
        dplyr::left_join(precip_06, by='NHDPlusID') %>%
        dplyr::left_join(precip_07, by='NHDPlusID') %>%
        dplyr::left_join(precip_08, by='NHDPlusID') %>%
        dplyr::left_join(precip_09, by='NHDPlusID') %>%
        dplyr::left_join(precip_10, by='NHDPlusID') %>%
        dplyr::left_join(precip_11, by='NHDPlusID') %>%
        dplyr::left_join(precip_12, by='NHDPlusID') %>%
        dplyr::left_join(temp_01, by='NHDPlusID') %>%
        dplyr::left_join(temp_02, by='NHDPlusID') %>%
        dplyr::left_join(temp_03, by='NHDPlusID') %>%
        dplyr::left_join(temp_04, by='NHDPlusID') %>%
        dplyr::left_join(temp_05, by='NHDPlusID') %>%
        dplyr::left_join(temp_06, by='NHDPlusID') %>%
        dplyr::left_join(temp_07, by='NHDPlusID') %>%
        dplyr::left_join(temp_08, by='NHDPlusID') %>%
        dplyr::left_join(temp_09, by='NHDPlusID') %>%
        dplyr::left_join(temp_10, by='NHDPlusID') %>%
        dplyr::left_join(temp_11, by='NHDPlusID') %>%
        dplyr::left_join(temp_12, by='NHDPlusID') %>%
        dplyr::filter(!(is.na(GageID))) %>% #there is a single reach that somehow got joined with no gauge and no data, just remove here
        dplyr::distinct()
    
    return(network)
}


#' buildCONUSnetwork
#'
#' Build hydrography network dataframe for ungauged prediction
#'
#' @param huc4 basin ID code
#' @param BHGmodel Bankfull hydraulic geometry models
#'
#' @return hydrography dataframe
buildCONUSnetwork <- function(huc4, BHGmodel){
    huc2 <- substr(huc4,1,2)

    #read in nhd
    network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        sf::st_zm() %>%
        dplyr::left_join(network_VAA, by='NHDPlusID') %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::select(c('NHDPlusID', 'GageID', 'WBArea_Permanent_Identifier', 'StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

    #load Indiana fixed basins (see manuscript)
    indiana <- readr::read_csv('data/indiana.csv')
    if(huc4 %in% indiana$huc4){
        network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/indiana/indiana_fixed_',huc4,'.shp'),quiet = TRUE) %>%
            dplyr::select(!path) %>%
            sf::st_zm() %>%
            dplyr::left_join(network_VAA, by='NHDPlusID') %>%
            dplyr::left_join(network_gages, by='NHDPlusID') %>%
            dplyr::mutate(WBArea_Permanent_Identifier = WBArea_Per)

        sf::st_geometry(network) <- "Shape"

        network <- network %>%
            dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))
    }

    #fix geometries (if applicable)
    network <- fixGeometries(network)

    #clip to CONUS land boundary (remove international streams and coastal estuaries like Long Island Sound and Chesapeake Bay)
    # CONUS boundary
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
        sf::st_transform(crs=sf::st_crs(network))
    
    network <- network %>%
        sf::st_intersection(states)

    #add reach type
    waterbodies <-  sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                            layer = 'NHDWaterbody',
                            quiet = TRUE) %>%
        sf::st_drop_geometry() %>%
        dplyr::select(c('Permanent_Identifier', 'FType', 'AreaSqKm'))
        colnames(waterbodies) <- c('Permanent_Identifier', 'FType', 'WBAreaSqKm')

    #join waterbodies/reach types
    network <- network %>%
        dplyr::left_join(waterbodies, by=c('WBArea_Permanent_Identifier'='Permanent_Identifier')) %>%
        dplyr::mutate(waterbody_type = ifelse(FType %in% c('436', '390') & WBAreaSqKm > 0, 'lake/reservoir', 'river')) %>% #to be a waterbody, must have an area > 0
        dplyr::relocate(WBArea_Permanent_Identifier, .after=waterbody_type) %>%
        dplyr::relocate(WBAreaSqKm, .after=WBArea_Permanent_Identifier) %>%
        dplyr::select(!FType)

    #remove impossible flowlines with no drainage area or divergent starting reaches (streamalc == 0; see pg. 45 at https://pubs.usgs.gov/of/2019/1096/ofr20191096.pdf)
    network <- network %>%
        dplyr::filter(AreaSqKm > 0 & TotDASqKm > 0 & StreamCalc > 0)  

    huc4id <- huc4
    basin <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
        dplyr::filter(huc4 == huc4id)
    basin <- fixGeometries(basin)

    #setup bankfull geometry
    sf::sf_use_s2(FALSE)
    regions <- sf::st_read('data/physio.shp')
    regions <- fixGeometries(regions)
    shp_temp <- sf::st_join(basin, regions, largest=TRUE) #take the physiographic region that the basin is primarily (dominant spatial intersection)

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

    network$a_Ab <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Ab)})
    network$b_Ab <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Ab)})
    network$mean_residual_Ab <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$mean_residual_Ab)})
    network$Ab_m2 <- 10^(network$a_Ab + network$b_Ab*log10(network$TotDASqKm)) * network$mean_residual_Ab #see other notes on bias correction

    network$Hb_m <- network$Ab_m2 / network$Wb_m

    # filter for rivers with bankfull width > 10m (matching inundation simulation)
    network <- network %>%
        dplyr::filter(AreaSqKm > 0 & Wb_m > 10) %>%
        dplyr::mutate(huc4=huc4)

    #load precip and temperature datasets from NHD
    flow_df <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusEROMMA', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'QBMA')) #naturalized flow

    precip_01 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM01', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM01')) #precip jan

    precip_02 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM02', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM02')) #precip feb
    
    precip_03 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM03', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM03')) #precip mar
    
    precip_04 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM04', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM04')) #precip apr
    
    precip_05 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM05', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM05')) #precip may
    
    precip_06 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM06', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM06')) #precip jun
    
    precip_07 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM07', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM07')) #precip jul
    
    precip_08 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM08', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM08')) #precip aug
    
    precip_09 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM09', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM09')) #precip sep
    
    precip_10 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM10', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM10')) #precip oct
    
    precip_11 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM11', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM11')) #precip nov
    
    precip_12 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrPrecipMM12', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'PrecipMM12')) #precip dec

    temp_01 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM01', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM01')) #temp jan

    temp_02 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM02', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM02')) #temp feb
    
    temp_03 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM03', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM03')) #temp mar
    
    temp_04 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM04', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM04')) #temp apr
    
    temp_05 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM05', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM05')) #temp may
    
    temp_06 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM06', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM06')) #temp jun
    
    temp_07 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM07', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM07')) #temp jul
    
    temp_08 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM08', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempVMM08')) #temp aug
    
    temp_09 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM09', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM09')) #temp sep
    
    temp_10 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM10', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM10')) #temp oct

    temp_11 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM11', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM11')) #temp nov
    
    temp_12 <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusIncrTempMM12', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'TempMM12')) #temp dec

    #filter add flow probs
    network <- network %>%
        dplyr::mutate(prob = purrr::map(.x = NHDPlusID, ~c(0.02, 0.10, 0.20, 0.50))) %>%
        tidyr::unnest(prob)

    #wrangle into single dataframe
    network <- network %>%
        dplyr::filter(waterbody_type == 'river') %>%
        dplyr::select(c('huc4', 'GageID', 'NHDPlusID', 'prob', 'Wb_m','Hb_m', 'Qb_cms', 'StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope')) %>%
        dplyr::left_join(flow_df, by='NHDPlusID') %>%
        dplyr::left_join(precip_01, by='NHDPlusID') %>%
        dplyr::left_join(precip_02, by='NHDPlusID') %>%
        dplyr::left_join(precip_03, by='NHDPlusID') %>%
        dplyr::left_join(precip_04, by='NHDPlusID') %>%
        dplyr::left_join(precip_05, by='NHDPlusID') %>%
        dplyr::left_join(precip_06, by='NHDPlusID') %>%
        dplyr::left_join(precip_07, by='NHDPlusID') %>%
        dplyr::left_join(precip_08, by='NHDPlusID') %>%
        dplyr::left_join(precip_09, by='NHDPlusID') %>%
        dplyr::left_join(precip_10, by='NHDPlusID') %>%
        dplyr::left_join(precip_11, by='NHDPlusID') %>%
        dplyr::left_join(precip_12, by='NHDPlusID') %>%
        dplyr::left_join(temp_01, by='NHDPlusID') %>%
        dplyr::left_join(temp_02, by='NHDPlusID') %>%
        dplyr::left_join(temp_03, by='NHDPlusID') %>%
        dplyr::left_join(temp_04, by='NHDPlusID') %>%
        dplyr::left_join(temp_05, by='NHDPlusID') %>%
        dplyr::left_join(temp_06, by='NHDPlusID') %>%
        dplyr::left_join(temp_07, by='NHDPlusID') %>%
        dplyr::left_join(temp_08, by='NHDPlusID') %>%
        dplyr::left_join(temp_09, by='NHDPlusID') %>%
        dplyr::left_join(temp_10, by='NHDPlusID') %>%
        dplyr::left_join(temp_11, by='NHDPlusID') %>%
        dplyr::left_join(temp_12, by='NHDPlusID')
    
    return(network)
}


#' trainModelEval_V
#'
#' Evaluate ml performance for predicting V_flood
#'
#' @param gageForModel_combined gauge dataframe for model training
#' @param nInnerFolds number inner layer folds for nested resampling cross-validation
#' @param nOuterFolds number outer layer folds for nested resampling cross-validation
#' @param numGrid Hyperparameter search space grid search parameter for nested resampling cross-validation
#' @param numRepeats number of repeats, to account for stochastic differences from nested resampling cross-validation
#'
#' @return model evaluation across nOuterFolds
trainModelEval_V <- function(gageForModel_combined, nInnerFolds, nOuterFolds, numGrid, numRepeats){
    library(tidymodels)
    library(parsnip)
    library(bonsai)

    #make sure all training sets comprise the same gauges
    data_fin <- gageForModel_combined %>%
        dplyr::mutate(V_m3 = log10(V_m3)) %>%
        dplyr::select(!Qexc_m3dy) %>%
        dplyr::select(!Q_m3dy) %>%
        dplyr::select(!NHDPlusID) %>%
        dplyr::select(!GageID)

    #tune using nested resampling procedure
    set.seed(76)
    folds_out <- vfold_cv(data_fin, repeats=numRepeats, v=nOuterFolds)

    #wrapper function for tuning across folds
    nestedWrapper <- function(object){
        set.seed(86)
        split <- object
        split_train <- training(split)
        split_test  <- testing(split)
        folds_in <- vfold_cv(split_train, repeats=1, v=nInnerFolds)

        #recipe
        recipe_full <- recipe(V_m3 ~ ., data=split_train) %>%
            step_cut(StreamCalc, breaks=seq(1,max(data_fin$StreamCalc),1)) %>% #one hot encoding of stream order
            step_dummy(StreamCalc) %>%
            step_normalize(all_numeric_predictors()) #normalize all features to mean 0, sd 1

        #model
        lightgbmboost_model <- boost_tree(tree_depth=tune(), min_n=tune(), trees=tune(), learn_rate=0.1) %>% #https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html and a default learning rate of 0.1 (learning rate and ntrees are inverses of one another so only one should be tuned)
            set_engine('lightgbm', num_leaves=tune()) %>% #added to bonsai package ad hoc rather than in tidymodels
            set_mode("regression")

        #workflow
        workflow <- workflow() %>%
            add_model(lightgbmboost_model) %>%
            add_recipe(recipe_full)

        #tune via grid search
        set.seed(835)
        tuned_workflow <- tune_grid(workflow, resamples = folds_in, grid=numGrid, verbose=TRUE)

        #evalualte hyperparameter tuning
        bestModel <- tuned_workflow %>%
            select_best(metric='rsq')

        final_wf <- workflow %>%
            finalize_workflow(bestModel)

        #fit the final, best model (across all folds) and evaluate on the withheld test set
        fit_fin <- final_wf %>%    
            last_fit(split)

        return(list('summary'=fit_fin,
            'model_fin'=extract_workflow(fit_fin)))
    }

    #run across folds
    result <- plyr::llply(folds_out$splits, nestedWrapper, .parallel=FALSE)

    return(result)
}


#' trainModelFin_V
#'
#' Train ml model on all data to predict V_flood from NHD
#'
#' @param gageForModel_combined gauge dataframe for model training
#' @param nInnerFolds number inner layer folds for nested resampling cross-validation
#'
#' @return model for predicting V_flood
trainModelFin_V <- function(gageForModel_combined, nInnerFolds, numGrid){
    library(tidymodels)
    library(parsnip)
    library(bonsai)

    #make sure all training sets comprise the same gauges
    data_fin <- gageForModel_combined %>%
        dplyr::mutate(V_m3 = log10(V_m3)) %>%
        dplyr::select(!Qexc_m3dy) %>%
        dplyr::select(!Q_m3dy) %>%
        dplyr::select(!NHDPlusID) %>%
        dplyr::select(!GageID)

    #tune using nested resampling procedure
    set.seed(76)

    #retrain workflow on all data (using cv for hyperparameter tuning)
    folds <- vfold_cv(data_fin, repeats=1, v=nInnerFolds)

    #recipe
    recipe_full <- recipe(V_m3 ~ ., data=data_fin) %>%
        step_cut(StreamCalc, breaks=seq(1,max(data_fin$StreamCalc),1)) %>% #one hot encoding of stream order
        step_dummy(StreamCalc) %>%
        step_normalize(all_numeric_predictors()) #normalize all features to mean 0, sd 1

    #model
    lightgbmboost_model <- boost_tree(tree_depth=tune(), min_n=tune(), trees=tune(), learn_rate=0.1) %>% #https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html and a default learning rate of 0.1 (learning rate and ntrees are inverses of one another so only one should be tuned)
        set_engine('lightgbm', num_leaves=tune()) %>% #added to bonsai package ad hoc rather than in tidymodels
        set_mode("regression")

    #workflow
    workflow <- workflow() %>%
        add_model(lightgbmboost_model) %>%
        add_recipe(recipe_full)

    #tune via grid search
    set.seed(875)
    tuned_workflow <- tune_grid(workflow, resamples = folds, grid=numGrid)

    #pick best hyperparameter set after tuning
    bestModel <- tuned_workflow %>%
        select_best(metric='rsq')

    final_wf <- workflow %>%
        finalize_workflow(bestModel)

    #finally, fit to all data
    set.seed(295)
    final_model <- final_wf %>% 
        fit(data_fin)

    return(final_model)
}


#' trainModelEval_Qf
#'
#' Evaluate ml performance for predicting Q_flood
#'
#' @param gageForModel_combined gauge dataframe for model training
#' @param nInnerFolds number inner layer folds for nested resampling cross-validation
#' @param nOuterFolds number outer layer folds for nested resampling cross-validation
#' @param numGrid Hyperparameter search space grid search parameter for nested resampling cross-validation
#' @param numRepeats number of repeats, to account for stochastic differences from nested resampling cross-validation
#'
#' @return model evaluation across nOuterFolds
trainModelEval_Qf <- function(gageForModel_combined, nInnerFolds, nOuterFolds, numGrid, numRepeats){
    library(tidymodels)
    library(parsnip)
    library(bonsai)

    #make sure all training sets comprise the same gauges
    data_fin <- gageForModel_combined %>%
        dplyr::mutate(Qexc_m3dy = log10(Qexc_m3dy)) %>%
        dplyr::select(!V_m3) %>%
        dplyr::select(!Q_m3dy) %>%
        dplyr::select(!NHDPlusID) %>%
        dplyr::select(!GageID)

    #tune using nested resampling procedure
    set.seed(76)
    folds_out <- vfold_cv(data_fin, repeats=numRepeats, v=nOuterFolds)

    #wrapper function for tuning across folds
    nestedWrapper <- function(object){
        set.seed(86)
        split <- object
        split_train <- training(split)
        split_test  <- testing(split)
        folds_in <- vfold_cv(split_train, repeats=1, v=nInnerFolds)

        #recipe
        recipe_full <- recipe(Qexc_m3dy ~ ., data=split_train) %>%
            step_cut(StreamCalc, breaks=seq(1,max(data_fin$StreamCalc),1)) %>% #one hot encoding of stream order
            step_dummy(StreamCalc) %>%
            step_normalize(all_numeric_predictors()) #normalize all features to mean 0, sd 1

        #model
        lightgbmboost_model <- boost_tree(tree_depth=tune(), min_n=tune(), trees=tune(), learn_rate=0.1) %>% #https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html and a default learning rate of 0.1 (learning rate and ntrees are inverses of one another so only one should be tuned)
            set_engine('lightgbm', num_leaves=tune()) %>% #added to bonsai package ad hoc rather than tidymodels
            set_mode("regression")

        #workflow
        workflow <- workflow() %>%
            add_model(lightgbmboost_model) %>%
            add_recipe(recipe_full)

        #tune via grid search
        set.seed(835)
        tuned_workflow <- tune_grid(workflow, resamples = folds_in, grid=numGrid, verbose=TRUE)

        #evaluate hyperparameter tuning
        bestModel <- tuned_workflow %>%
            select_best(metric='rsq')

        final_wf <- workflow %>%
            finalize_workflow(bestModel)

        #fit the final, best model (across all folds) and evaluate on the withheld test set
        fit_fin <- final_wf %>%    
            last_fit(split)

        return(list('summary'=fit_fin,
            'model_fin'=extract_workflow(fit_fin)))
    }

    #run across folds
    result <- plyr::llply(folds_out$splits, nestedWrapper, .parallel=FALSE)

    return(result)
}


#' trainModelFin_Qf
#'
#' Train ml model on all data to predict Q_flood from NHD
#'
#' @param gageForModel_combined gauge dataframe for model training
#' @param nInnerFolds number inner layer folds for nested resampling cross-validation
#'
#' @return model for predicitng Q_flood
trainModelFin_Qf <- function(gageForModel_combined, nInnerFolds, numGrid){
    library(tidymodels)
    library(parsnip)
    library(bonsai)

    #make sure all training sets comprise the same gauges
    data_fin <- gageForModel_combined %>%
        dplyr::mutate(Qexc_m3dy = log10(Qexc_m3dy)) %>%
        dplyr::select(!V_m3) %>%
        dplyr::select(!Q_m3dy) %>%
        dplyr::select(!NHDPlusID) %>%
        dplyr::select(!GageID)

    #tune using nested resampling procedure
    set.seed(76)

    #retrain workflow on all data (using cv for hyperparameter tuning)
    folds <- vfold_cv(data_fin, repeats=1, v=nInnerFolds)

    #recipe
    recipe_full <- recipe(Qexc_m3dy ~ ., data=data_fin) %>%
        step_cut(StreamCalc, breaks=seq(1,max(data_fin$StreamCalc),1)) %>% #one hot encoding of stream order
        step_dummy(StreamCalc) %>%
        step_normalize(all_numeric_predictors()) #normalize all features to mean 0, sd 1

    #model
    lightgbmboost_model <- boost_tree(tree_depth=tune(), min_n=tune(), trees=tune(), learn_rate=0.1) %>% #https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html and a default learning rate of 0.1 (learning rate and ntrees are inverses of one another so only one should be tuned)
        set_engine('lightgbm', num_leaves=tune()) %>% #added to bonsai package ad hoc rather than tidymodels
        set_mode("regression")

    #workflow
    workflow <- workflow() %>%
        add_model(lightgbmboost_model) %>%
        add_recipe(recipe_full)

    #tune via grid search
    set.seed(875)
    tuned_workflow <- tune_grid(workflow, resamples = folds, grid=numGrid)

    #pick best hyperparameter set after tuning
    bestModel <- tuned_workflow %>%
        select_best(metric='rsq')

    final_wf <- workflow %>%
        finalize_workflow(bestModel)

    #finally, fit to all data
    set.seed(295)
    final_model <- final_wf %>% 
        fit(data_fin)

    return(final_model)
}


#' trainModelEval_Q
#'
#' Evaluate ml performance for predicitng Q_total
#'
#' @param gageForModel_combined gauge dataframe for model training
#' @param nInnerFolds number inner layer folds for nested resampling cross-validation
#' @param nOuterFolds number outer layer folds for nested resampling cross-validation
#' @param numGrid Hyperparameter search space grid search parameter for nested resampling cross-validation
#' @param numRepeats number of repeats, to account for stochastic differences from nested resampling cross-validation
#'
#' @return model evaluation across nOuterFolds
trainModelEval_Q <- function(gageForModel_combined, nInnerFolds, nOuterFolds, numGrid, numRepeats){
    library(tidymodels)
    library(parsnip)
    library(bonsai)

    #make sure all training sets comprise the same gauges
    data_fin <- gageForModel_combined %>%
        dplyr::mutate(Q_m3dy = log10(Q_m3dy)) %>%
        dplyr::select(!V_m3) %>%
        dplyr::select(!Qexc_m3dy) %>%
        dplyr::select(!NHDPlusID) %>%
        dplyr::select(!GageID)

    #tune using nested resampling procedure
    set.seed(76)
    folds_out <- vfold_cv(data_fin, repeats=numRepeats, v=nOuterFolds)

    nestedWrapper <- function(object){
        set.seed(86)
        split <- object
        split_train <- training(split)
        split_test  <- testing(split)
        folds_in <- vfold_cv(split_train, repeats=1, v=nInnerFolds)

        #recipe
        recipe_full <- recipe(Q_m3dy ~ ., data=split_train) %>%
            step_cut(StreamCalc, breaks=seq(1,max(data_fin$StreamCalc),1)) %>% #one hot encoding of stream order
            step_dummy(StreamCalc) %>%
            step_normalize(all_numeric_predictors()) #normalize all features to mean 0, sd 1

        #model
        lightgbmboost_model <- boost_tree(tree_depth=tune(), min_n=tune(), trees=tune(), learn_rate=0.1) %>% #https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html and a default learning rate of 0.1 (learning rate and ntrees are inverses of one another so only one should be tuned)
            set_engine('lightgbm', num_leaves=tune()) %>% #added to bonsai package ad hoc rather than tidymodels
            set_mode("regression")

        #workflow
        workflow <- workflow() %>%
            add_model(lightgbmboost_model) %>%
            add_recipe(recipe_full)

        #tune via grid search
        set.seed(835)
        tuned_workflow <- tune_grid(workflow, resamples = folds_in, grid=numGrid, verbose=TRUE)

        #evaluate hyperparameter tuning
        bestModel <- tuned_workflow %>%
            select_best(metric='rsq')

        final_wf <- workflow %>%
            finalize_workflow(bestModel)

        #fit the final, best model (across all folds) and evaluate on the withheld test set
        fit_fin <- final_wf %>%    
            last_fit(split)

        return(list('summary'=fit_fin,
            'model_fin'=extract_workflow(fit_fin)))
    }

    #run aross folds
    result <- plyr::llply(folds_out$splits, nestedWrapper, .parallel=FALSE)

    return(result)
}


#' trainModelFin_Q
#'
#' Train ml model on all data to predict Q_total from NHD
#'
#' @param gageForModel_combined gauge dataframe for model training
#' @param nInnerFolds number inner layer folds for nested resampling cross-validation
#'
#' @return model for predicitng Q_total
trainModelFin_Q <- function(gageForModel_combined, nInnerFolds, numGrid){
    library(tidymodels)
    library(parsnip)
    library(bonsai)

    #make sure all training sets comprise the same gauges
    data_fin <- gageForModel_combined %>%
        dplyr::mutate(Q_m3dy = log10(Q_m3dy)) %>%
        dplyr::select(!V_m3) %>%
        dplyr::select(!Qexc_m3dy) %>%
        dplyr::select(!NHDPlusID) %>%
        dplyr::select(!GageID)

    #tune using nested resampling procedure
    set.seed(76)

    #retrain workflow on all data (using cv for hyperparameter tuning)
    folds <- vfold_cv(data_fin, repeats=1, v=nInnerFolds)

    #recipe
    recipe_full <- recipe(Q_m3dy ~ ., data=data_fin) %>%
        step_cut(StreamCalc, breaks=seq(1,max(data_fin$StreamCalc),1)) %>% #one hot encoding of stream order
        step_dummy(StreamCalc) %>%
        step_normalize(all_numeric_predictors()) #normalize all features to mean 0, sd 1

    #model
    lightgbmboost_model <- boost_tree(tree_depth=tune(), min_n=tune(), trees=tune(), learn_rate=0.1) %>% #https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html and a default learning rate of 0.1 (learning rate and ntrees are inverses of one another so only one should be tuned)
        set_engine('lightgbm', num_leaves=tune()) %>% #added to bonsai package ad hoc rather than tidymodels
        set_mode("regression")

    #workflow
    workflow <- workflow() %>%
        add_model(lightgbmboost_model) %>%
        add_recipe(recipe_full)

    #tune via grid search
    set.seed(875)
    tuned_workflow <- tune_grid(workflow, resamples = folds, grid=numGrid)

    #pick best hyperparameter set after tuning
    bestModel <- tuned_workflow %>%
        select_best(metric='rsq')

    final_wf <- workflow %>%
        finalize_workflow(bestModel)

    #finally, fit to all data
    set.seed(295)
    final_model <- final_wf %>% 
        fit(data_fin)

    return(final_model)
}


#' predictBasin
#'
#' Run three ML models on hydrography, predicitng Q_flood, V_flood, and Q_total
#'
#' @param huc4 basin ID code
#' @param conusDF basin hydrography dataframe for model prediction
#' @param model_Qf trained Q_flood model
#' @param model_V trained V_flood model
#' @param model_Q trained Q_total model
#'
#' @return hydrography dataframe with predictions
predictBasin <- function(huc4, conusDF, model_Qf, model_V, model_Q){
    if(nrow(conusDF)==0){
        return(data.frame())
    }

    #prep for models
    sf::sf_use_s2(FALSE)
    library(tidymodels)

    forPredict <- conusDF %>%
        sf::st_drop_geometry()

    nhdIDs <- forPredict$NHDPlusID

    forPredict <- forPredict %>%
        dplyr::select(!c('huc4', 'GageID', 'NHDPlusID', 'Wb_m', 'Hb_m', 'Qb_cms',))

    #predict river-floodplain exchange terms
    forPredict$log10V_m3 <- predict(model_V, forPredict)$.pred
    forPredict$log10Qexc_m3dy <- predict(model_Qf, forPredict)$.pred
    forPredict$log10Q_m3dy <- predict(model_Q, forPredict)$.pred
    forPredict$NHDPlusID <- nhdIDs #add ids back to prediction df

    #wrangle dataframe
    forPredict_fin <- forPredict %>%
        dplyr::select(c('NHDPlusID','prob', 'log10V_m3', 'log10Qexc_m3dy', 'log10Q_m3dy'))

    #additional calculations for residence time
    conusDF <- conusDF %>%
        dplyr::left_join(forPredict_fin, by=c('NHDPlusID', 'prob')) %>%
        dplyr::mutate(log10tau_hr = (log10V_m3 - log10Qexc_m3dy) + log10(24), #no channel or bankfull excess, just what's in the fp
                    log10tau_channel_hr = log10(Wb_m*Hb_m*LengthKM*1000) - log10(Qb_cms*86400) + log10(24)) %>% #no bankfull excess, just the channel
        dplyr::select(c('huc4', 'NHDPlusID', 'LengthKM', 'StreamCalc', 'AreaSqKm', 'TotDASqKm', 'prob', 'log10V_m3', 'log10Qexc_m3dy', 'log10Q_m3dy', 'log10tau_hr', 'log10tau_channel_hr', 'Shape'))

    return(conusDF)
}


#' summarizeBasinSO
#'
#' summarize model predictions by basin, flood size, and streamorder
#'
#' @param huc4 basin ID code
#' @param basinPredictions sf dataframe with predictions
#'
#' @return dataframe summarized by flood size and streamorder
summarizeBasinSO <- function(huc4, basinPredictions){
    if(nrow(basinPredictions)==0){
        return(data.frame())
    }

    #wrangle by stream order and flood size
    out <- basinPredictions %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(prob, StreamCalc) %>%
        dplyr::summarize(
                    median_tau_channel_hr_km = median(10^(log10tau_channel_hr)/LengthKM, na.rm=T),
                    median_tau_hr_km = median(10^(log10tau_hr)/LengthKM, na.rm=T)) %>%
        dplyr::mutate(huc4 = huc4)

    return(out)
}


#' summarizeBasin
#'
#' summarize model predictions by basin and flood size
#'
#' @param huc4 basin ID code
#' @param basinPredictions sf dataframe with predictions
#'
#' @return dataframe summarized by flood size
summarizeBasin <- function(huc4, basinPredictions){
    if(nrow(basinPredictions)==0){
        return(data.frame())
    }

    #wrangle by flood size
    out <- basinPredictions %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(prob) %>%
        dplyr::summarize(
                    median_tau_channel_hr_km = median(10^(log10tau_channel_hr)/LengthKM, na.rm=T),
                    median_tau_hr_km = median(10^(log10tau_hr)/LengthKM, na.rm=T),
                    total_exchange_frac = sum(10^(log10Qexc_m3dy), na.rm=T)/sum(10^(log10Q_m3dy), na.rm=T)) %>%
        dplyr::mutate(huc4 = huc4) %>%
        dplyr::mutate(total_exchange_frac = ifelse(total_exchange_frac > 1, 1, total_exchange_frac)) #a handful of basins are slightly above 1 because the numbers come from independent models, so we just round them off

    return(out)
}


#' makeMapBasin
#'
#' prep basin sf object for mapping
#'
#' @param basinPredictions sf dataframe with predictions
#' @param chosen_prob flood size
#'
#' @return basin sf object with mapping bin generated, given a flood size
makeMapBasin <- function(basinPredictions, chosen_prob){
    if(nrow(basinPredictions)==0){
        return(data.frame())
    }

    library(sf)

    #bin model results for mapping
	basinPredictions <- basinPredictions %>%
        dplyr::filter(prob == chosen_prob) %>%
        dplyr::mutate(tau_col = dplyr::case_when(
            log10(10^(log10tau_hr)/LengthKM) <= -0.5 ~ '-0.5'
            ,log10(10^(log10tau_hr)/LengthKM) <= 0 ~ '0'
            ,log10(10^(log10tau_hr)/LengthKM) <= 0.5 ~ '0.5'
            ,log10(10^(log10tau_hr)/LengthKM) <= 1 ~ '1'
            ,TRUE ~ '1+'
        ))

    #convert to factor
    basinPredictions$tau_col <- factor(basinPredictions$tau_col, levels = c('-0.5', '0', '0.5', '1', '1+'))

    return(basinPredictions)
}


#' assignVolVals
#'
#' Assign basin codes to hydrodynamic inundation rasters
#'
#' @param vol_grids List of hydrodynamic volume rasters used for inundation verification
#'
#' @return lookup table between inundation rasters and basin ID codes
assignVolVals <- function(vol_grids){

    huc4s <- terra::vect('data/path_to_data/CONUS_connectivity_data/HUC4s.shp')

    #loop through usgs models and assign huc4 basin codes (for matchups)
    out <- data.frame()
    for(i in vol_grids){
        rast <- terra::rast(i)
        res_og <- terra::res(rast)[1]
        unit_og <- terra::linearUnits(rast)
        rast <- terra::project(rast, terra::crs(huc4s)) #use the huc4 projection

        ints <- terra::intersect(huc4s, terra::ext(rast)) #find regional basin that the usgs model is within

        temp <- as.data.frame(ints) %>%
            dplyr::select('huc4') %>%
            dplyr::mutate('rast'=i,
                        'res_og_m'=res_og*unit_og)

        out <- rbind(out, temp)
    }

    return(out)
}


#' collectValReaches
#'
#' Collect all modeled reaches to test hydrodynamic inundation simulations
#' Using gauges as proxy for USGS volume model mainstems b/c they are calibrated to these specific gauges on the mainstem
#' This is a useful proxy for ensuring we only validate at reaches where the flood model actually passes through the entire catchment (bc their associatied gauges are not on the edge of the model domains but on the mainstems)
#'
#' @param huc4id basin ID code
#'
#' @return hydrography reach IDs
collectValReaches <- function(huc4id){
    sf::sf_use_s2(FALSE)

    #prep
    huc2 <- substr(huc4id, 1, 2)

    #read in huc4 basin
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE)
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)

    #join to gauged reaches
    unit_catchments <- unit_catchments %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        dplyr::filter(is.na(GageID) == 0)

    return(unit_catchments$NHDPlusID)
}


#' wrangleDepthGrids
#'
#' Get estimated volumes from local USGS hydrodynamic inundation simulations given an observed depth
#' Also, prep for Q_flood comparison to jacknife regression estimates
#' 
#' @param huc4id basin ID code
#' @param reaches_val NHD reaches that may overlap with hydrodynamic inundation simulatios
#' @param volVal All the hydrodynamic inundation simulations
#'
#' @return V_flood comparsion result
wrangleDepthGrids <- function(huc4id, reaches_val, volVal){
    huc2 <- substr(huc4id,1, 2)

    #grab model of interest
    volVal <- volVal %>%
        dplyr::filter(huc4 %in% huc4id)

    if(nrow(volVal) == 0){
        return(data.frame())
    }

    #get nhd hydrography for basin
    network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'GageID'))

    network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
        dplyr::left_join(network_gages, by='NHDPlusID') %>%
        sf::st_zm() %>%
        dplyr::filter(NHDPlusID %in% reaches_val)

    #load Indiana fixed basins (see manuscript)
    indiana <- readr::read_csv('data/indiana.csv')
    if(huc4id %in% indiana$huc4){
        network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/indiana/indiana_fixed_',huc4id,'.shp'),quiet = TRUE) %>%
            dplyr::select(!path) %>%
            sf::st_zm() %>%
            dplyr::left_join(network_gages, by='NHDPlusID') %>%
            dplyr::filter(NHDPlusID %in% reaches_val)

        sf::st_geometry(network) <- "Shape"
    }

    if(nrow(network) == 0){
        return(data.frame())
    }
    
    #pull in the unit catchments used for reach-by-reach comparisons
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE) %>%
        dplyr::filter(NHDPlusID %in% network$NHDPlusID)
    
    unit_catchments_terra <- terra::vect(unit_catchments)
    
    #get depth grids from usgs models
    grids <- volVal$rast
    grid_names <- basename(volVal$rast)
    if(length(grids) == 0){
        return(data.frame())
    }

    #make sure verything is listed in the same order for the merging of the extract
    network <- network[order(network$NHDPlusID),]
    unit_catchments <- unit_catchments[order(unit_catchments$NHDPlusID),]
    unit_catchments_terra <- unit_catchments_terra[order(unit_catchments_terra$NHDPlusID),]

    #loop through and extract modeled volumes from usgs models
    k <- 1
    out <- data.frame()
    for(i in grids){
        grid <- terra::rast(i)
        grid <- terra::project(grid, "epsg:4269") #NAD1983 unprojected to match below
        name <- grid_names[k]

        #get depth at nhd flowline
        depths <- terra::extract(grid, unit_catchments_terra, fun=function(x){mean(x, na.rm=T)})
        flowline_depth_m <- depths[,2]*0.3048 #ft to m

        #get actual volume
        grid_d <- grid
        res <- terra::res(grid_d)[1]
        grid_d <- terra::project(grid_d, "epsg:4269")
        grid_d <- grid_d * 0.3048 #ft to m
        grid_a <- terra::cellSize(grid_d, unit="m", transform=TRUE)

        grid_v <- grid_d * grid_a #m3
        vols <- terra::extract(grid_v, unit_catchments_terra, fun=function(x){sum(x, na.rm=T)})
        vols <- vols[,2]

        #prep output
        temp <- data.frame('grid'= substr(name, 1, nchar(name)-4),
                        'NHDPlusID'=network$NHDPlusID,
                        'GageID'=network$GageID,
                        'depth_m'=flowline_depth_m)

        temp <- temp %>%
            dplyr::mutate(V_usgs_m3 = vols) %>%
            dplyr::filter(depth_m > 0 & V_usgs_m3 > 0)
        
        if(nrow(temp)==0){next}
        
        out <- rbind(out, temp)

        k <- k + 1
    }

    if(nrow(out)==0){
        return(data.frame())
    }

    #final wrangle
    out$huc4 <- huc4id
    out <- out %>%
        dplyr::relocate(huc4)

    return(out)
}


#' modelsJacknifeBHG
#'
#' Fit jacknife regression models for bankful hydraulic geometry
#' 
#' 
#' @return dataframe of model results
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

    #prep for jacknife
    dataset$id <- 1:nrow(dataset)

    #loop through data, witholding i data point to fit jacknife regressions
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


#' prepBankfullHydraulics
#'
#' Computes bankfull hydraulic geometry for gauge (to calculate Q_flood)
#' 
#' @param mode 'deploy' or 'val'
#' @param gageRecord streamgauge record dataframe
#' @param gage streamgauge metadata dataframe
#' @param BHGmodel_jacknife bankfull hydraulic models fit via jacknife regression (for val mode)
#' @param BHGmodel bankfull hydraulic models fit (for deploy mode) 
#' @param gageRecordStart starting date for streamflow records
#' @param gageRecordEnd ending data for streamflow records 
#'
#' @return dataframe with bankfull hydraulics added to gauge record
prepBankfullHydraulics <- function(mode, gageRecord, gage, BHGmodel_jacknife, BHGmodel){
    if(nrow(gage %>% dplyr::bind_rows())==0){return(data.frame())}
    if(nrow(gageRecord %>% dplyr::bind_rows())==0){return(data.frame())}

    gageRecord <- dplyr::bind_rows(gageRecord) %>%
        dplyr::select(c('site_no', 'Q_cms', 'exceed_prob'))

    #filter for flows above bankfull, to calculate exceedance probabilities for flood events
    gage <- sf::st_drop_geometry(gage) %>% 
        dplyr::bind_rows() %>% 
        dplyr::select(c('site_no', 'DA_skm', 'physio_region')) %>% 
        dplyr::filter(site_no %in% gageRecord$site_no) %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(!(is.na(physio_region))) %>%
        dplyr::distinct()

    #wrangling
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

    #fit banfull hydraulics models
    gage$Qb_cms_log10 <- gage$a_Qb + gage$b_Qb*log10(gage$DA_skm)
    gage$Wb_m_log10 <- gage$a_Wb + gage$b_Wb*log10(gage$DA_skm)
    gage$Ab_m2_log10 <- gage$a_Ab + gage$b_Ab*log10(gage$DA_skm)

    #handle transform bias
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

    #handle transform bias
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


#' calc_Qexc
#'
#' Calculates Q_flood_event for floods in streamgauge record
#' 
#' @param mode 'deploy' or 'val'
#' @param gageRecord streamgauge record dataframe
#' @param gagePrepped streamgauge metadata dataframe
#' @param depAHG depth AHG model parameters 
#' @param minAHGr2 minimum allowable r2 for depth AHG fit (if smaller than minAHGr2, remove)
#'
#' @return dataframe of flood event variables
calc_Qexc <- function(mode, gageRecord, gagePrepped, depAHG, minAHGr2){
    if(nrow(gageRecord)==0){
        return(data.frame())
    }

    #prep
    if(mode == 'val'){
        probs <- gagePrepped %>%
            dplyr::filter(site_no == gageRecord[1,]$site_no) %>% #make sure dataset matches gageRecord dataset
            dplyr::filter(Qb_cms > 0 & Qb_cms_OBS > 0 & Ub_ms > 0 & Ub_ms_OBS > 0 & Wb_m > 0 & Wb_m_OBS > 0)
    }

    if(mode == 'deploy'){
        probs <- gagePrepped %>%
            dplyr::filter(site_no == gageRecord[1,]$site_no) %>% #make sure dataset matches gageRecord dataset
            dplyr::filter(Qb_cms > 0 & Ub_ms > 0  & Wb_m > 0)
    }

    if(nrow(probs)==0){
        return(data.frame())
    }

    depAHG <- dplyr::bind_rows(depAHG) %>%
        dplyr::filter(site_no == gageRecord[1,]$site_no) %>%  #make sure dataset matches gageRecord dataset
        dplyr::filter(lm_r2 >= minAHGr2)

    if(nrow(depAHG)==0){
        return(data.frame())
    }

    #compute flow depth using AHG
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

        #loop through and assign flood ids to multi-day floods
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

        #calculate metrics for each identified flood (using modeled Qb)
        events <- gageRecord %>%
            dplyr::filter(flood_id > 0) %>%
            dplyr::mutate(Ub_ms = probs[1,]$Ub_ms,
                        Wb_m = probs[1,]$Wb_m,
                        dy = dplyr::row_number()) %>%
            dplyr::group_by(flood_id)%>%
            dplyr::summarise(flag='model',
                            n_yrs=years,
                            Qexc_m3dy = (pracma::trapz(dy, Q_cms*86400) - pracma::trapz(dy, (Ub_ms*Wb_m*Htf_m*86400)))/n(), #event floodplain flux (m3/dy)
                            Q_m3dy = pracma::trapz(dy, Q_cms*86400)/n(), #total flow
                            Htf_m = mean(Htf_m)) 

        if(nrow(events)==0){
            return(data.frame())
        }

        #calculate metrics for each identified flood (using observed Qb)
        events_obs <- gageRecord %>%
            dplyr::filter(flood_id_OBS > 0) %>%
            dplyr::mutate(Ub_ms = probs[1,]$Ub_ms_OBS,
                        Wb_m = probs[1,]$Wb_m_OBS,
                        dy = dplyr::row_number()) %>%
            dplyr::group_by(flood_id_OBS)%>%
            dplyr::summarise(flag='obs',
                            n_yrs=years,
                            Qexc_m3dy = (pracma::trapz(dy, Q_cms*86400) - pracma::trapz(dy, (Ub_ms*Wb_m*Htf_m*86400)))/n(), #event floodplain flux (m3/dy)
                            Q_m3dy = pracma::trapz(dy, Q_cms*86400)/n(), #total flow
                            Htf_m = mean(Htf_m)) 

        if(nrow(events_obs)==0){
            return(data.frame())
        }

        colnames(events_obs) <- c('flood_id', 'flag','n_yrs','Qexc_m3dy', 'Q_m3dy', 'Htf_m')

        out <- rbind(events, events_obs)
        out$site_no <- probs[1,]$site_no
        out$Qb_obs <- probs[1,]$Qb_cms_OBS
    }

    #identify floods
    #deploy version
    else if(mode == 'deploy'){
        gageRecord$flood <- ifelse((gageRecord$Q_cms > probs$Qb_cms) & (gageRecord$Q_cms > (probs$Ub_ms*probs$Wb_m*gageRecord$Htf_m)), 'yes', 'no')
        gageRecord$flood_id <- NA

        #loop through and assign flood ids to multi-day floods
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

        #calculate metrics for each identified flood (using modeled Qb)
        events <- gageRecord %>%
            dplyr::filter(flood_id > 0) %>%
            dplyr::mutate(Ub_ms = probs[1,]$Ub_ms,
                        Wb_m = probs[1,]$Wb_m,
                        dy = dplyr::row_number()) %>%
            dplyr::group_by(flood_id)%>%
            dplyr::summarise(date = dplyr::first(dateTime),
                            n_yrs=years,
                            Qexc_m3dy = (pracma::trapz(dy, Q_cms*86400) - pracma::trapz(dy, (Ub_ms*Wb_m*Htf_m*86400)))/n(), #event floodplain flux (m3/dy)
                            Q_m3dy = pracma::trapz(dy, Q_cms*86400)/n(), #total flow
                            Htf_m = mean(Htf_m)) #event stage (i.e. crest height) (m))

        if(nrow(events)==0){
            return(data.frame())
        }

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


#' prepGage
#'
#' Get streamgauge metadata
#' 
#' @param gageID streamgauge ID code
#'
#' @return dataframe gauge metadata
prepGage <- function(gageID){
    if(is.na(gageID)){
        return(data.frame())
    }

    #grab gauge metadata
    site <- dataRetrieval::readNWISsite(siteNumbers = gageID) %>%
        dplyr::mutate(lat = dec_lat_va,
                    lon = dec_long_va,
                    DA_skm = drain_area_va * 2.58999) %>% #mi2 to km2
        dplyr::select(c('site_no', 'lat', 'lon', 'DA_skm'))

    #get physiograhic region
    sf::sf_use_s2(FALSE)
    regions <- sf::st_read('data/physio.shp')
    regions <- fixGeometries(regions)
    site_shp <- tryCatch(sf::st_as_sf(site, coords=c('lon', 'lat'), crs=sf::st_crs(4269)),
                        error = function(m){return(data.frame())})
    if(nrow(site_shp)==0) {return(data.frame())} #if no gage data, move on

    shp_temp <- sf::st_join(site_shp, regions) #take the physiographic region that the basin is primarily in (dominant spatial intersection)
    site_shp$physio_region <- shp_temp$DIVISION

    #convert to df
    site <- site_shp %>%
        dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                    lat = sf::st_coordinates(.)[,2]) %>%
        sf::st_drop_geometry()

    return(site)
}


#' prepFlowRecord
#'
#' Build gauge flow record
#' 
#' @param gage dataframe of gauge metadata
#' @param gageRecordStart Start of time period for gauge record
#' @param gageRecordEnd End of time period for gauge record
#' @param minRecordLength minimum length of years
#'
#' @return streamgauge flow record dataframe
prepFlowRecord <- function(gage, gageRecordStart, gageRecordEnd, minRecordLength){
    set.seed(435)

    if(nrow(gage)==0){
        return(data.frame())
    }

    #prep
    gageID <- gage$site_no

    #grab daily discharge data
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

    #only keep flowing days
    gagedata <- gagedata %>%
        dplyr::filter(Q_cms > 0)

    # convert to date (workaround to handle known date class bug with midnight- https://github.com/tidyverse/lubridate/issues/1124)
    gagedata$date <- lubridate::ymd_hms(format(as.POSIXct(gagedata$dateTime), format = "%Y-%m-%d %T %Z"))
    gagedata$year <- lubridate::year(gagedata$date)

    #get number of years on record
    start <- min(gagedata$year)
    end <- max(gagedata$year)
    spread <- end - start

    if(spread < minRecordLength){return(data.frame())}
    if(nrow(gagedata)== 0){return(data.frame())}

    gagedata$spread <- spread

    #get exceedance probs for each flow value
    gagedata$rank <- rank(-gagedata$Q_cms)
    gagedata$exceed_prob <- (gagedata$rank/(nrow(gagedata) + 1)) #exceed prob for a given day over 'spread' years

    gagedata <- gagedata %>%
        dplyr::filter(spread >= minRecordLength & date >= gageRecordStart & date <= gageRecordEnd) #ensure previous filtering worked

    if(nrow(gagedata)== 0){return(data.frame())}

    return(gagedata)
}


#' prepInundationData
#'
#' Prep for DEM inundation
#' Compute HAND rasters
#' Run for a single model reach (reachID), used within the indunation wrapper function
#' 
#' @param huc4 basin ID code
#' @param reachID Model reach ID
#' @param bankfullWidth reach bankfull width
#' @param dem dem raster
#' @param d8 flow direction raster
#'
#' @return hand raster
prepInundationData <- function(huc4, reachID, bankfullWidth, dem, d8){
    #prep
    huc2 <- substr(huc4, 1, 2)

    #read in nhd unit catchment
    sql <- paste0('SELECT * FROM NHDPlusCatchment WHERE NHDPlusID = ', reachID)
    unit_catchment <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                        layer = 'NHDPlusCatchment',
                                        query = sql, #sql query for efficient loading of desired catchment
                                        quiet = TRUE)

    unit_catchment <- terra::vect(unit_catchment)
    unit_catchment <- terra::project(unit_catchment, terra::crs(dem))

    #read in the nhd centerline
    sql <- paste0('SELECT * FROM NHDFlowline WHERE NHDPlusID = ', reachID)
    flowline <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                    layer = 'NHDFlowline',
                                    query = sql, #sql query for efficient loading of desired catchment
                                    quiet = TRUE) %>%
        dplyr::filter(NHDPlusID == reachID) %>%
        sf::st_transform(crs=sf::st_crs(terra::crs(dem))) %>%
        sf::st_zm()

    flowline_terra <- terra::vect(flowline)

    #match the nhdplus unit catchment 
    dem_clipped <- terra::crop(dem, unit_catchment)
    dem_clipped <- terra::mask(dem_clipped, unit_catchment)

    #burn levee heights into dem (using national levee database)
    nld <- sf::st_read('data/path_to_data/CONUS_connectivity_data/NLD_floodwalls.shp')
    nld$levee_height_cm <- nld$WALL_HEIGH * 30.48 #ft to cm
    nld <- terra::vect(nld)
    nld <- terra::project(nld, terra::crs(dem_clipped))
    raster_template = terra::rast(terra::ext(dem_clipped), resolution = 10,crs = terra::crs(dem_clipped))
    nld_raster <- terra::rasterize(nld, raster_template, field = "levee_height_cm")
    nld_raster <- terra::subst(nld_raster, NA, 0)
    dem_clipped <- dem_clipped + nld_raster

    #fill and breach the dem to calculate flow directions
    r <- terra::rast(flowline_terra, resolution=10, extent=terra::ext(dem_clipped))
    flowline_rast <- terra::rasterize(flowline_terra, r)
    flowline_rast[flowline_rast == 1] = -5000 #50m following https://doi.org/10.1016/j.rse.2024.114333
    flowline_rast[is.na(flowline_rast)] = 0
    flowline_rast <- terra::mask(flowline_rast, dem_clipped)

    dem_burned <- dem_clipped + flowline_rast
    dem_burned <- flowdem::breach(dem_burned) #https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.10648
    hydrodem <- flowdem::fill(dem_burned, epsilon=TRUE) #https://www.sciencedirect.com/science/article/pii/S0098300413001337
    d8_clipped <- terra::terrain(hydrodem, 'flowdir')

    #snap pour point to the largest flow accumulation cell (to avoid incorrect snapping based on nhd topology)
    flowaccum <- terra::flowAccumulation(d8_clipped) 
    flowaccum <- terra::mask(flowaccum, unit_catchment)

    #extract dem
    dem_clipped <- dem_clipped*0.01 #cm to m

    #use HAND method to get relative elevations for inundation mapping
    #get thalweg mask
    thalweg_pour_points <- terra::extract(flowaccum, flowline_terra, xy=TRUE)

    #get pixel unit catchments and HAND
    pixel_pour_points <- na.omit(thalweg_pour_points[order(thalweg_pour_points$flowdir),]) #terra doesn't automatically rename the bands, this is actually flow accum from previous function

    #first, ensure the pour point is on the dem... (there are some edge cases where dem and flowline don't align perfectly). In these cases pass an error
    if(nrow(pixel_pour_points) == 0){
        return(list('hand'=NA,
                'bankfullMask'=NA,
                'flag'='no pour points'))
    }

    #do pour points sequentially, moving downstream and computing HAND
    for(i in 1:nrow(pixel_pour_points)){
        pixel_pour_point <- terra::vect(pixel_pour_points[i,], c('x','y'))
        pixel_watershed <- terra::watershed(x=d8_clipped, pourpoint=cbind(pixel_pour_points[i,]$x, pixel_pour_points[i,]$y))
        pixel_watershed[pixel_watershed == 0] <- NA
        pixel_dem <- terra::mask(dem_clipped, pixel_watershed)

        #if most upstream catchment, there's nothing upstream to remove
        if(i == 1){
            #get HAND
            elev <- terra::extract(pixel_dem, pixel_pour_point)
            rel_elev <- pixel_dem - elev[1,]$elev_cm #band wasn't renamed, this was previously converted to m
            next
        }

        #remove upstream pixel_watersheds, i.e. isolate the unit pixel watershed
        pixel_watershed_up <- terra::watershed(x=d8_clipped, pourpoint=cbind(pixel_pour_points[i-1,]$x, pixel_pour_points[i-1,]$y))
        pixel_watershed_up[pixel_watershed_up == 0] <- NA
        pixel_mask <- terra::app(terra::sds(pixel_watershed, pixel_watershed_up), fun="sum", na.rm=T)
        pixel_mask[pixel_mask == 2] <- NA

        #get per-pixel unit catchment dem
        pixel_dem <- terra::mask(pixel_dem, pixel_mask)

        #get HAND
        elev <- terra::extract(pixel_dem, pixel_pour_point)
        rel_elev_temp <- pixel_dem - elev[1,]$elev_cm #band wasn't renamed, this was previously converted to m
        rel_elev <- terra::merge(rel_elev, rel_elev_temp)
    }

    hand <- rel_elev
    hand[hand < 0] = 0

    #return what's needed for inundation modeling
    return(list('hand'=hand,
                'flag'='has pour points'))
}


#' modelInundation
#'
#' Maps inundation
#' 
#' @param inundationData inundation package
#' @param floodDepth flood stage for specific event
#'
#' @return inundation map
modelInundation <- function(inundationData, floodDepth){
    library(terra)
    hand <- inundationData$hand

    #inundation calculation
    floodDepths <- floodDepth - hand
    floodMap <- floodDepths
    floodMap[floodMap < 0] <- 0
    floodMap[floodMap > 0] <- 1

    floodDepths[floodDepths < 0] <- 0 # if depth is negative (rare artifact of models and HAND), just set it to 0

    return(list('binary'=floodMap,
                'depths'=floodDepths))
}


#' buildDepthAHG
#'
#' Fit a depth AHG curve and estimate flow depth for flood events
#' 
#' @param gage streamgauge metadata dataframe
#' @param minADCPMeas minimum allowed field measurements to fit depth AHG (if smaller, skip gauge)
#'
#' @return depth AHG parameters and fit statistics
buildDepthAHG <- function(gage, minADCPMeas){
    if(nrow(gage)==0){return(data.frame())}

    #prep
    gageID <- gage$site_no

    #read in situ depth~flow measurements
    surfaceData <- tryCatch(dataRetrieval::readNWISmeas(gageID, expanded = TRUE),
                                    error=function(m){return(data.frame('site_no'=gageID,
                                                        'lm_r2'=NA,
                                                        'lm_depth_a'=NA,
                                                        'lm_depth_b'=NA))})

    if(nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
                                                        'lm_r2'=NA,
                                                        'lm_depth_a'=NA,
                                                        'lm_depth_b'=NA))}

    #wrangle depth~flow paired in situ measurements
    surfaceData <- surfaceData %>%
        dplyr::filter(measured_rating_diff %in% c('Fair','Good', 'Excellent')) %>% #Fair are < 8% flow diff, Good are < 5% flow diff, and Excellent are <2% flow diff
        dplyr::mutate('Q_cms'=chan_discharge* 0.0283, #imperial to metric
                    'width_m'=chan_width*0.3048,
                    'area_m2'=chan_area*0.092903,
                    'velocity_ms'=chan_velocity*0.3048) %>%
        dplyr::select(c('site_no', 'Q_cms', 'width_m', 'area_m2', 'velocity_ms')) %>%
        dplyr::mutate(depth_m=area_m2 / width_m) %>%
        dplyr::filter(depth_m > 0 & Q_cms > 0 & is.finite(depth_m) & is.finite(Q_cms))

    #remove likely erronous outlier depth measurements that can overlay sway rating curve estimates
    iqr <- IQR(surfaceData$depth_m, na.rm=T)

    surfaceData <- surfaceData %>%
        dplyr::filter(depth_m < (quantile(depth_m, 0.75, na.rm=T) + 1.5*iqr)) %>%
        dplyr::filter(depth_m > (quantile(depth_m, 0.25, na.rm=T) - 1.5*iqr))

    if(nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
                                                        'lm_r2'=NA,
                                                        'lm_depth_a'=NA,
                                                        'lm_depth_b'=NA))}

    #fit depth ahg
    lm_depth <- lm(log10(depth_m)~log10(Q_cms), data=surfaceData)
    lm_depth_a <- coef(lm_depth)[1]
    lm_depth_b <- coef(lm_depth)[2]
    lm_depth_bias <- mean(10^(lm_depth$residuals), na.rm=T)

    #wrangle metadata export
    site_info <- data.frame('site_no'=gageID,
                            'lm_r2'=summary(lm_depth)$r.squared,
                            'lm_depth_a'=lm_depth_a,
                            'lm_depth_b'=lm_depth_b,
                            'lm_depth_bias_correct'=lm_depth_bias)

    rownames(site_info) <- NULL

    return(site_info)
}


#' modelsBHG
#'
#' Fit bakfull hydraulics models using Bieger et al. (2015) data
#'
#' @return list of models
modelsBHG <- function(){
    dataset <- readr::read_csv('data/bhg_us_database_bieger_2015.csv') %>% #available by searching for paper at https://swat.tamu.edu/search
        dplyr::select(c('Physiographic Division', '...9', '...11','...13','...17')) #some necessary manual munging for colnames from dataset

    colnames(dataset) <- c('DIVISION', 'DA_km2', 'Qb_cms','Wb_m','Ab_m2')

    dataset$Wb_m <- as.numeric(dataset$Wb_m)
    dataset$Ab_m2 <- as.numeric(dataset$Ab_m2)
    dataset$Qb_cms <- as.numeric(dataset$Qb_cms)
    dataset$DA_km2 <- as.numeric(dataset$DA_km2)

    dataset <- tidyr::drop_na(dataset)

    division <- toupper(sort(unique(dataset$DIVISION))) #make lowercase and sort

    #build models, grouped by physiographic region
    models <- dataset %>%
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
                    see_Ab = summary(model_Ab)$sigma) %>%
        dplyr::mutate(division = division)

    models[models$division == "INTERMONTANE PLATEAU",]$division <- "INTERMONTANE PLATEAUS" #make sure names line up

    #we found that the WB_m model for the Interior Highlands, i.e. mountainious region in Arkansas, Missouri, Oklahoma, and Ilinois, had really poor fit (r2: 0.27)
    #and a very high model intercept. Likely due to extremely small sample size (7 sites). So, we swap that model for the region that is most similar, the Appalachian Highlands.
    models[models$division == 'INTERIOR HIGHLANDS',]$a_Wb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$a_Wb
    models[models$division == 'INTERIOR HIGHLANDS',]$b_Wb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$b_Wb
    models[models$division == 'INTERIOR HIGHLANDS',]$r2_Wb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$r2_Wb
    models[models$division == 'INTERIOR HIGHLANDS',]$mean_residual_Wb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$mean_residual_Wb
    models[models$division == 'INTERIOR HIGHLANDS',]$see_Wb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$see_Wb

    models[models$division == 'INTERIOR HIGHLANDS',]$a_Qb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$a_Qb
    models[models$division == 'INTERIOR HIGHLANDS',]$b_Qb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$b_Qb
    models[models$division == 'INTERIOR HIGHLANDS',]$r2_Qb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$r2_Qb
    models[models$division == 'INTERIOR HIGHLANDS',]$mean_residual_Qb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$mean_residual_Qb
    models[models$division == 'INTERIOR HIGHLANDS',]$see_Qb <- models[models$division == 'APPALACHIAN HIGHLANDS',]$see_Qb

    models[models$division == 'INTERIOR HIGHLANDS',]$a_Ab <- models[models$division == 'APPALACHIAN HIGHLANDS',]$a_Ab
    models[models$division == 'INTERIOR HIGHLANDS',]$b_Ab <- models[models$division == 'APPALACHIAN HIGHLANDS',]$b_Ab
    models[models$division == 'INTERIOR HIGHLANDS',]$r2_Ab <- models[models$division == 'APPALACHIAN HIGHLANDS',]$r2_Ab
    models[models$division == 'INTERIOR HIGHLANDS',]$mean_residual_Ab <- models[models$division == 'APPALACHIAN HIGHLANDS',]$mean_residual_Ab
    models[models$division == 'INTERIOR HIGHLANDS',]$see_Ab <- models[models$division == 'APPALACHIAN HIGHLANDS',]$see_Ab

    return(models)
}


#' dataBHG
#'
#' Wrangle Bieger et al. (2015) data
#'
#' @return cleaned dataset as dataframe
dataBHG <- function(){
    dataset <- readr::read_csv('data/bhg_us_database_bieger_2015.csv') %>% #available by searching for paper at https://swat.tamu.edu/search
        dplyr::select(c('USGS Station No.', 'Physiographic Division', '...9', '...11','...13')) #some necessary manual munging for colnames from dataset

    colnames(dataset) <- c('GageID', 'DIVISION', 'DA_gage_skm','Qb_cms','Wb_m_obs')

    dataset$Qb_cms <- as.numeric(dataset$Qb_cms)
    dataset$Wb_m_obs <- as.numeric(dataset$Wb_m_obs)
    dataset$DA_gage_skm <- as.numeric(dataset$DA_gage_skm)

    dataset <- tidyr::drop_na(dataset) #only keep those with a usgs site and a non-NA Qb value (some Bieger data didn't report Qb)

    return(dataset)
}


#' assignRegion
#'
#' Assign physiographic region to huc4 basin. Use the most overlapping region for a given huc4
#'
#' @return dataframe lookup table
assignRegion <- function(){
    library(sf)

    #loop through huc2 basins, running spatial join for each subset of huc4 basins
    basin <- data.frame()
    for(i in c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18')){
        basin_temp <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', i, '/WBD_', i, '_HU2_Shape/Shape/WBDHU4.shp'))
        basin_temp <- fixGeometries(basin_temp)

        sf::sf_use_s2(FALSE)
        regions <- sf::st_read('data/physio.shp')
        regions <- fixGeometries(regions)
        shp_temp <- sf::st_join(basin_temp, regions, largest=TRUE) #take the physiographic region that the basin is primarily (dominant spatial intersection)

        temp <- shp_temp %>%
            sf::st_drop_geometry() %>%
            dplyr::select(c('huc4', 'DIVISION', 'PROVINCE', 'SECTION')) 

        basin <- rbind(basin, temp)
    }

    #only keep the 205 huc4 basins used in this study
    hucs <- c('0101', '0102', '0103', '0104', '0105', '0106', '0107', '0108', '0109', '0110',
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
    
    out <- basin %>%
        dplyr::filter(huc4 %in% hucs)
    
    return(out)
}


#' findRegulatedReaches
#'
#' Find regulated river reaches in huc4 basin, i.e. reaches immediately downstream of a dammed reach
#'
#' @param huc4 basin code
#' @param basinPredictions sf dataframe of model predictions for huc4 basin
#' @param perc_thresh threshold for drainage area match to join dam to reach
#' @param buffer_dist radius of buffer around dam to search for 'best matching' river reaches to join to.
#'
#' @return list of 2 dataframes: unregulated basin reaches and regulated basin reaches
findRegulatedReaches <- function(huc4, basinPredictions, perc_thresh, buffer_dist){
    library(sf)
    huc2 <- substr(huc4, 1, 2)

    if(nrow(basinPredictions) == 0){return(data.frame())}

    #load network topology df for huc4 basin
    network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE) %>%
        dplyr::select(c('NHDPlusID', 'FromNode', 'ToNode'))

    #load global dam watch database and filter for US
    gdw <- sf::st_read('data/path_to_data/CONUS_connectivity_data/GDW_barriers_v1_0.shp') %>%
        dplyr::filter(COUNTRY == 'United States') %>%
        sf::st_transform(crs=sf::st_crs(basinPredictions))

    #sptially join dams to modeled reaches using quality control hyperparameters in _targets.R
    dammed_reaches <- sf::st_join(basinPredictions, gdw, join = sf::st_is_within_distance, dist = buffer_dist)
    dammed_reaches <- dammed_reaches %>%
        dplyr::mutate(diff = abs(CATCH_SKM - TotDASqKm)/TotDASqKm) %>%
        dplyr::filter(!(is.na(diff))) %>%
        dplyr::group_by(NHDPlusID, prob) %>%
        dplyr::slice_min(diff) %>%
        dplyr::filter(diff <= perc_thresh) %>%
        dplyr::left_join(network_VAA, by='NHDPlusID')
    
    regulated_reaches <- basinPredictions %>%
        dplyr::left_join(network_VAA, by='NHDPlusID') %>%
        dplyr::filter(FromNode %in% dammed_reaches$ToNode) #filter for reaches immedately downstream of dammed reach
    
    other_reaches <- basinPredictions %>%
        dplyr::filter(!(NHDPlusID %in% c(regulated_reaches$NHDPlusID))) %>%
        dplyr::filter(!(NHDPlusID %in% c(dammed_reaches$NHDPlusID))) #filter for reaches that are (1) undammed and (2) not immedately downstream of a dammed reach

    return(list('other_reaches'=other_reaches,
                'regulated_reaches'=regulated_reaches))
}


#' writeToFile
#'
#' Write model results to file for export
#'
#' @param huc4 basin code
#' @param basinPredictions sf dataframe of model predictions for huc4 basin
#'
#' @return text confirming whether written to file
writeToFile <- function(huc4, basinPredictions){
    if(nrow(basinPredictions)==0){
        return('empty df')
    }

    out <- basinPredictions %>%
        sf::st_drop_geometry() %>%
        dplyr::select(c('huc4', 'NHDPlusID', 'prob','LengthKM','log10tau_hr','log10Qexc_m3dy'))

    readr::write_csv(out, paste0('cache/results/results_',huc4,'.csv'))
    return('printed to file')
}