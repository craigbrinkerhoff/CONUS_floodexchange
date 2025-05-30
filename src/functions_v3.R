



buildGageFloodFunctions_val <- function(huc4, BHGmodel, gageDF) {
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




buildGageFloodFunctions <- function(mode, huc4, BHGmodel, gageDF, code) {
    if(mode == 'deploy'){
        #prep gages for filtering
        gageDF <- dplyr::bind_rows(gageDF) %>%
            dplyr::filter(Qexc_m3dy > 0) %>% # integration necessitates storms be 2+ days long to calculate a flux
            dplyr::group_by(site_no) %>%
            dplyr::summarise(id = code,
                            length_dys = mean(length_dy),
                            Htf_m = mean(Htf_m, na.rm=T),
                            Qexc_m3dy = mean(Qexc_m3dy, na.rm=T))
    }
    else if(mode == 'val'){
        #dummy dataframe to get rest of function to run
        gageDF <- data.frame(site_no=gageDF$NHDPlusID,
                            id = 'volume_validation',
                            length_dys=NA,
                            Htf_m=gageDF$depth_m,
                            Qexc_m3dy=NA)
    }
    else{
        return('error! Please set mode to val or deploy')
    }

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

    if(mode == 'deploy'){
        network <- network %>%
            dplyr::filter(GageID %in% gageDF$site_no)
    }
    else if(mode == 'val'){
        network <- network %>%
            dplyr::filter(NHDPlusID %in% gageDF$site_no)
    }
    else{
        return('error! Please set mode to val or deploy')
    }

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
    if(mode == 'deploy'){
        network <- network %>%
            dplyr::left_join(gageDF, by=c('GageID'='site_no'))
    }
    if(mode == 'val'){
        network <- network %>%
            dplyr::left_join(gageDF, by=c('NHDPlusID'='site_no'))
    }

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
        dplyr::select(c('NHDPlusID', 'Shape')) %>%
        dplyr::left_join(network_fin, by='NHDPlusID') %>%
        dplyr::relocate(Shape, .after = tidyselect::last_col())

    return(network)
}



collectValReaches <- function(huc4id, usgs_maps){
    sf::sf_use_s2(FALSE)

    #grab maps within our huc4
    usgs_maps <- usgs_maps %>%
        sf::st_transform(crs=sf::st_crs(4269)) %>% #just to make sure
        dplyr::mutate(key_exdprob = keyProb) %>%
        dplyr::filter(huc4 == huc4id)

    #prep
    huc2 <- substr(huc4id, 1, 2)

    #read in huc4 basin
    unit_catchments <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4id,'_HU4_GDB/NHDPLUS_H_',huc4id,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    quiet = TRUE)

    unit_catchments <- unit_catchments %>%
        sf::st_filter(usgs_maps, .predicate = sf::st_within) #Only keep model reaches that fall within the usgs maps. These maps are inconsistent in space, so we can really only validate in places where we know the entire model reach has been mapped by USGS (to avoid conflating real false positives with incomplete USGS maps)

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
        dplyr::filter(NHDPlusID %in% reaches_val)
    
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
        depths <- terra::extract(grid, network, fun=function(x){mean(x, na.rm=T)}) #SHOULD REALLY BE MAX I THINK...
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
            tidyr::drop_na()
        
        out <- rbind(out, temp)

        k <- k + 1
    }
    return(out)
}