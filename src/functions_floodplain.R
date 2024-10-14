## Functions for floodplain model
## Craig Brinkerhoff
## Summer 2024


prepGage <- function(path_to_data, gageID){
  site <- dataRetrieval::readNWISsite(siteNumbers = gageID) %>%
    dplyr::mutate(lat = dec_lat_va,
                  lon = dec_long_va,
                  DA_skm = drain_area_va * 2.58999) %>% #mi2 to km2
    dplyr::select(c('site_no', 'lat', 'lon', 'DA_skm'))
  
  #get physiograhic region (for banfull depth model)
  sf::sf_use_s2(FALSE)
  regions <- sf::st_read('data/physio.shp') #physiographic regions
  regions <- fixGeometries(regions)
  site_shp <- sf::st_as_sf(site, coords=c('lon', 'lat'), crs=sf::st_crs(4269))
  
  shp_temp <- sf::st_join(site_shp, regions) #take the physiographic region that the basin is mostly in (dominant spatial intersection)
  site_shp$physio_region <- shp_temp$DIVISION
  
  return(site_shp)
}






prepFlowRecord <- function(gage, gageRecordStart, gageRecordEnd, minRecordLength){
  #prep
  gageID <- gage$site_no
  
  #grab sub-daily discharge data
  gagedata <- tryCatch(dataRetrieval::readNWISdata(sites = gageID,
                                          service='uv', #15' or 30' stage and flow
                                          parameterCd = c('00065', '00060'), #discharge and stage parameter codes [cfs and ft, respectively]
                                          startDate = gageRecordStart,
                                          endDate = gageRecordEnd),
                      error = function(m){return(data.frame())})
  
  if(nrow(gagedata)==0) {return(data.frame())} #if no gage data, move on
  if(!("X_00065_00000" %in% colnames(gagedata))){return(data.frame())} #if incorrect stage column name (happens sometimes), move on
  if(!("X_00060_00000" %in% colnames(gagedata))){return(data.frame())} #if incorrect discharge column name (happens sometimes), move on
  
  #hardcoded (for now at least) downscampling of sensor data to 30'
  #if sensor records 15' data, downsample to 30' data to match other gages
  if(as.numeric(gagedata$dateTime[2] - gagedata$dateTime[1]) == 15){ #15'
    gagedata <- gagedata[seq(1,nrow(gagedata),2),]
  }
  
  #if sensor records 5' data, downsample to 30' data to match other gages
  if(as.numeric(gagedata$dateTime[2] - gagedata$dateTime[1]) == 5){ #5'
    gagedata <- gagedata[seq(1,nrow(gagedata),6),]
  }
  
  #convert to metric
  gagedata$Q_cms <- gagedata$X_00060_00000 * 0.0283 #cfs to cms
  gagedata$log10_Q_cms <- log10(gagedata$Q_cms)
  gagedata$stage_m <- gagedata$X_00065_00000 * 0.3048 #ft to m
  
  #only keep flowing events
  gagedata <- gagedata %>%
    dplyr::filter(stage_m > 0 & Q_cms > 0)

  # get flood event duration from flow timeseries (workaround to handle known date class bug with midnight- https://github.com/tidyverse/lubridate/issues/1124)
  gagedata$date <- lubridate::ymd_hms(format(as.POSIXct(gagedata$dateTime), format = "%Y-%m-%d %T %Z"))

  #get number of years on record
  start <- lubridate::year(min(gagedata$date))
  end <- lubridate::year(max(gagedata$date))
  spread <- end - start
  
  if(spread < minRecordLength){return(data.frame())}

  gagedata$spread <- spread

  #get exceedance probs for each flow value
  probs <- quantile(gagedata$Q_cms, 1 - seq(0,1,0.00001)) #calc exdeeance probs down to 0.001%
  names(probs) <- NULL
  gagedata$exceed_prob <- sapply(gagedata$Q_cms, function(x){which.min(abs(probs - x)) / length(probs)})

  return(gagedata)
}




prepInundationData <- function(path_to_data, huc4, reachID, bankfullWidth, dem, d8){
  #prep
  huc2 <- substr(huc4, 1, 2)
  
  sql <- paste0('SELECT * FROM NHDPlusCatchment WHERE NHDPlusID = ', reachID)
  gage_unit_catchment <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                     layer = 'NHDPlusCatchment',
                                     query = sql, #sql query to only load the catchment we want efficiently (see above)
                                     quiet = TRUE)
  
  gage_unit_catchment <- terra::vect(gage_unit_catchment)
  gage_unit_catchment <- terra::project(gage_unit_catchment, terra::crs(dem))

  #read in the nhd centerline
  sql <- paste0('SELECT * FROM NHDFlowline WHERE NHDPlusID = ', reachID)
  flowline <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                     layer = 'NHDFlowline',
                                     query = sql, #sql query to only load the catchment we want efficiently (see above)
                                     quiet = TRUE) %>%
    dplyr::filter(NHDPlusID == reachID) %>%
    sf::st_transform(crs=sf::st_crs(terra::crs(dem))) %>%
    sf::st_zm()

  flowline_terra <- terra::vect(flowline)

  #math the nhdplus unit catchment 
  dem_clipped <- terra::crop(dem, gage_unit_catchment)
  dem_clipped <- terra::mask(dem_clipped, gage_unit_catchment)
  
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
  flowaccum <- terra::mask(flowaccum, gage_unit_catchment)

  #extract gage's dem
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
      gage_elev <- terra::extract(pixel_dem, pixel_pour_point)
      rel_elev <- pixel_dem - gage_elev[1,]$elev_cm #band wasn't renamed, this is already in m
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
    gage_elev <- terra::extract(pixel_dem, pixel_pour_point)
    rel_elev_temp <- pixel_dem - gage_elev[1,]$elev_cm #band wasn't renamed, this is already in m
    rel_elev <- terra::merge(rel_elev, rel_elev_temp)
  }
  
  hand <- rel_elev
  hand[hand < 0] = 0
  
  #mask out the bankfull channel (given bankful width)
  flowline <- terra::vect(flowline)
  bankfullMask <- terra::buffer(flowline, bankfullWidth) #don't let it go below the spatial resolution of the dem
  r <- terra::rast(bankfullMask, resolution=10, extent=terra::ext(hand))
  bankfullMask <- terra::rasterize(bankfullMask, r)
  bankfullMask[is.na(bankfullMask[])] <- 0
  bankfullMask <- terra::subst(bankfullMask, 1, 2) #reacast mask as value '2'
  
  reachWBAreaPermanentIdentifier <- flowline$WBArea_Permanent_Identifier
  
  #return what's needed for inundation modeling
  return(list('hand'=hand,
              'bankfullMask'=bankfullMask,
              'flag'='has pour points'))
}




modelInundation <- function(inundationData, floodDepth){
  hand <- inundationData$hand
  bankfullMask <- inundationData$bankfullMask
  
  #actual inundation calculation
  floodDepths <- floodDepth - hand
  floodMap <- floodDepths
  floodMap[floodMap < 0] <- 0
  floodMap[floodMap > 0] <- 1
  floodMap <- floodMap + bankfullMask #add bankfull channel mask to flood map

  floodMap <- terra::subst(floodMap, 3, 2) #recast 'bankfullMask + flooded' as river, to copacetic with 'river only'
  
  floodDepths[floodDepths < 0] <- 0

  return(list('binary'=floodMap,
              'depths'=floodDepths))
}






# buildFloodwaterProfile <- function(path_to_data, huc4, basinData, gage, floodEvents){
#   if(nrow(floodEvents)==0){return(data.frame())} #if gage didn't run under floodEvents, don't run now
  
#   huc2 <- substr(huc4, 1, 2)
  
#   #grab all necessary parameters from floodevent record
#   gageID <- floodEvents[1,]$site_no
#   bankfullWidth <- floodEvents[1,]$Wb_m
#   floodDepths <- floodEvents$Htf_eventmean_m
  
#   #setup basin shapefiles
#   dem <- terra::unwrap(basinData$terra$dem) #[cm] need to unwrap the packed terra package (for distributed computing)
#   d8 <- terra::unwrap(basinData$terra$d8) #[cm] need to unwrap the packed terra package (for distributed computing)
  
#   #just for the gage analysis, we need to get the NHD reach IDs for the gages a priori (repeats what's done within prepInundationData() out of necessity)
#   reachID <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE) %>%
#     dplyr::filter(GageID == gageID)
#   if(nrow(reachID)==0){return(data.frame())}
#   reachID <- reachID$NHDPlusID
  
#   #grab waterbody type (if any)
#   lookup <- basinData$waterbodyLookUp %>%
#     dplyr::filter(!(is.na(NHDPlusID_reach)))
#   waterbody_type <- lookup[as.character(lookup$NHDPlusID_reach) == reachID,]$waterbody_type
#   waterbody_type <- ifelse(length(waterbody_type) == 0, 'floodplain',
#                            ifelse(waterbody_type == '466', 'wetland',
#                                   ifelse(waterbody_type == '436', 'reservoir',
#                                          ifelse(waterbody_type == '390', 'lake', 'other'))))

#   #build necessary datasets for inundation modeling
#   inundationDataPackage <- prepInundationData(path_to_data, huc4, reachID, bankfullWidth, dem, d8)

#   #setup lookup table
#   floodplainLookupTable <- data.frame()
#   img_ticker <- 0.001 #arbitrary number to name imgs and keep them in order...
  
#   #loop through floods and build floodplina profile lookup table
#   for(i in floodDepths) {
    
#     floodMap <- modelInundation(inundationDataPackage, i)
    
#     #save map
#     png(paste0('cache/for_gif/',img_ticker,'.png'))
#     terra::plot(floodMap, box=FALSE)
#     dev.off()
    
#     flooded_pixels <- terra::freq(floodMap)
#     flood_area <- flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
#     flood_area <- ifelse(length(flood_area)==0, 0, flood_area)
    
#     floodMap <- data.frame('Htf_m'=i, #cm to m
#                        'Af_m2'=flood_area)
#     floodplainLookupTable <- rbind(floodMap, floodplainLookupTable)
    
#     img_ticker <- img_ticker + 0.001
#   }
  
#   #make gif using magick package
#   imgs <- list.files('cache/for_gif/')
#   img_list <- lapply(paste0('cache/for_gif/',imgs), magick::image_read)
#   img_joined <- magick::image_join(img_list)
#   img_animated <- magick::image_animate(img_joined, fps=2)
#   magick::image_write(img_animated, paste0("cache/gagePlots/", gageID, '_', waterbody_type, '_gif.gif'))
#   unlink('cache/for_gif/*')
  
#   #make floodplain profile plot
#   theme_set(theme_classic())
  
#   floodplain_elevation_profile <- ggplot(floodplainLookupTable, aes(x=Af_m2, y=Htf_m)) +
#     geom_point(size=5, color='#218f8d') +
#     ylab(bquote(bold("Flood level ["*m*"]"))) +
#     xlab(bquote(bold("Inundated area ["*m^2*"]"))) +
#     theme(axis.title = element_text(face = 'bold', size=18),
#           axis.text = element_text(size=15))
#   ggsave(paste0('cache/gagePlots/', gageID, '_floodwater_profile.png'), floodplain_elevation_profile, width=7, height=7)
  
#   #add flow back to lookuptable so we can validate later
#   floodplainLookupTable$Q_eventmean_cms <- floodEvents$Q_eventmean_cms
#   floodplainLookupTable$site_no <- floodEvents$site_no
  
#   floodplainLookupTable <- floodplainLookupTable %>%
#     dplyr::relocate(site_no, Q_eventmean_cms)
  
#   return(floodplainLookupTable)
# }






buildGageBankfullModel <- function(gage, minADCPMeas){
  #prep
  gageID <- gage$site_no
  
  #prep rating table
  ratingTable <- dataRetrieval::readNWISrating(gageID, "exsa")
  ratingTable$stage_m <- ratingTable$INDEP * 0.3048# + ratingTable$SHIFT #ft to m
  ratingTable$Q_cms <- ratingTable$DEP * 0.0283
  ratingTable$log10_Q_cms <- log10(ratingTable$Q_cms)
  
  #prep custom rating table from available adcp measurements
  surfaceData <- dataRetrieval::readNWISmeas(gageID, expanded = TRUE) %>%
    dplyr::filter(measured_rating_diff %in% c('Good', 'Excellent')) %>%
    dplyr::mutate('stage_m'=gage_height_va * 0.3048,
                  'Q_cms'=chan_discharge* 0.0283,
                  'width_m'=chan_width*0.3048,
                  'area_m2'=chan_area*0.092903,
                  'velocity_ms'=chan_velocity * 0.3048) %>%
    dplyr::select(c('site_no', 'stage_m', 'Q_cms', 'width_m', 'area_m2', 'velocity_ms')) %>%
    dplyr::mutate(depth_m=Q_cms / (width_m*velocity_ms),
                  log10_Q_cms = log10(Q_cms))
  
  #only keep flowing events
  ratingTable <- ratingTable %>%
    dplyr::filter(stage_m > 0 & Q_cms > 0)
  
  surfaceData <- surfaceData %>%
    dplyr::filter(stage_m > 0 & Q_cms > 0)
  
  #remove likely erronous outlier depth measurements
  iqr <- IQR(surfaceData$depth_m, na.rm=T)
  
  surfaceData <- surfaceData %>%
    dplyr::filter(depth_m < (quantile(depth_m, 0.75, na.rm=T) + 1.5*iqr)) %>%
    dplyr::filter(depth_m > (quantile(depth_m, 0.25, na.rm=T) - 1.5*iqr))
  
  if(nrow(ratingTable)==0 | nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
                                                                    'DA_skm'=NA,
                                                                    'seg_r2'=NA,
                                                                    'seg_nbreaks'=NA,
                                                                    'depthstage_r2'=NA,
                                                                    'bankfull_stage_mu_m'=NA,
                                                                    'bankfull_stage_se_m'=NA,
                                                                    'bankfull_depth_mu_m'=NA,
                                                                    'bankfull_depth_lower'=NA,
                                                                    'bankfull_depth_upper'=NA,
                                                                    'bankfull_Q_mu_m'=NA,
                                                                    'bankfull_Q_lower'=NA,
                                                                    'bankfull_Q_upper'=NA))}
  
  #segmented regression
  #linear regression for segmented
  lm <- lm(log10_Q_cms~stage_m, data=ratingTable)

  #fit breakpoint models with 0-2 breakpoints, keeping the one with the best bic
  seg <- segmented::selgmented(lm, Kmax=2, type='bic')

  bankfull_stage_mu_m <- seg$psi[nrow(seg$psi),2] #always use higher of the breakpoints
  bankfull_stage_se_m <- seg$psi[nrow(seg$psi),3] #always use higher of the breakpoints
  
  #get bankfull Q
  breakpoint_Q_mean <- (ratingTable[which(ratingTable$INDEP == round(bankfull_stage_mu_m*3.28084,2)),]$DEP)*0.0283
  breakpoint_Q_lower <- (ratingTable[which(ratingTable$INDEP == round((bankfull_stage_mu_m-bankfull_stage_se_m)*3.28084,2)),]$DEP)*0.0283
  breakpoint_Q_upper <- (ratingTable[which(ratingTable$INDEP == round((bankfull_stage_mu_m+bankfull_stage_se_m)*3.28084,2)),]$DEP)*0.0283
  
  #get bankfull depth
  lm_depth <- lm((depth_m)~(stage_m), data=surfaceData)
  lm_depth_a <- coef(lm_depth)[1]
  lm_depth_b <- coef(lm_depth)[2]
  
  bankfull_depth_mean <- lm_depth_a + lm_depth_b*bankfull_stage_mu_m
  bankfull_depth_lower <- lm_depth_a + lm_depth_b*(bankfull_stage_mu_m-bankfull_stage_se_m) 
  bankfull_depth_upper <- lm_depth_a + lm_depth_b*(bankfull_stage_mu_m+bankfull_stage_se_m)
  
  theme_set(theme_classic())
  plot <- ggplot(data=surfaceData, aes(x=stage_m, y=depth_m)) +
    geom_point(size=2) +
    geom_smooth(method='lm', se=F, size=1.5)+
    xlab(bquote(bold("Stage ["*m*"]"))) +
    ylab(bquote(bold("Depth ["*m*"]"))) +
    theme(axis.title = element_text(face = 'bold', size=18),
          axis.text = element_text(size=15))
  ggsave(paste0('cache/gagePlots/', gageID, '_stagedepth.png'), width=8, height=8)
  
  #save a plot of the bankfull regression model
  theme_set(theme_classic())
  
  temp <- data.frame('breakpoint_stage_mean' = bankfull_stage_mu_m,
                     'breakpoint_stage_sd' = bankfull_stage_se_m,
                     'breakpoint_Q_mean' = breakpoint_Q_mean,
                     'breakpoint_Q_lower' = breakpoint_Q_lower,
                     'breakpoint_Q_upper' = breakpoint_Q_upper)

  plot <- ggplot() +
    geom_point(data=ratingTable, aes(x=stage_m, y=Q_cms), alpha=0.3, size=1.5) +
    geom_pointrange(data=temp, aes(x=breakpoint_stage_mean, y=breakpoint_Q_mean, xmin = breakpoint_stage_mean-breakpoint_stage_sd, xmax = breakpoint_stage_mean + breakpoint_stage_sd, ymin = breakpoint_Q_lower, ymax = breakpoint_Q_upper), color='darkred', size=1)+
    scale_y_log10() +
    xlab(bquote(bold("Stage ["*m*"]"))) +
    ylab(bquote(bold("Discharge ["*m^3*s^-1*"]"))) +
    theme(axis.title = element_text(face = 'bold', size=18),
          axis.text = element_text(size=15))
  ggsave(paste0('cache/gagePlots/', gageID, '_bankfullModel.png'), width=8, height=8)

  #setup site info for later export
  site_info <- data.frame('site_no'=gageID,
                          'DA_skm'=gage$DA_skm,
                          'seg_r2'=summary(seg)$r.squared,
                          'seg_nbreaks'=nrow(seg$psi),
                          'depthstage_r2'=summary(lm_depth)$r.squared,
                          'depthstage_nObs'=nrow(surfaceData),
                          'bankfull_stage_mu_m'=bankfull_stage_mu_m,
                          'bankfull_stage_se_m'=bankfull_stage_se_m,
                          'bankfull_depth_mu_m'=bankfull_depth_mean,
                          'bankfull_depth_lower'=bankfull_depth_lower,
                          'bankfull_depth_upper'=bankfull_depth_upper,
                          'bankfull_Q_mu_m'=breakpoint_Q_mean,
                          'bankfull_Q_lower'=breakpoint_Q_lower,
                          'bankfull_Q_upper'=breakpoint_Q_upper)
  return(site_info)
}









validateBankfullModel <- function(bankfullModel, BHGdata){
  library(ggplot2)
  theme_set(theme_classic())
  
  df <- dplyr::bind_rows(bankfullModel)

  
  bankfullValidationSet_Qb <- df %>%
    dplyr::select(c('site_no', 'bankfull_Q_mu_m')) %>%
    dplyr::left_join(BHGdata, by='site_no') %>%
    dplyr::filter(!(is.na(Qb_cms)))
  
  bankfullValidationSet_Hb <- df %>%
    dplyr::select(c('site_no', 'bankfull_depth_mu_m')) %>%
    dplyr::left_join(BHGdata, by='site_no') %>%
    dplyr::mutate(Hb_m = as.numeric(Hb_m)) %>%
    dplyr::filter(!(is.na(Hb_m)))

  plot_Qb <- ggplot(bankfullValidationSet_Qb, aes(x=Qb_cms, y=bankfull_Q_mu_m)) +
    geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
    geom_point(size=5) +
    geom_smooth(method='lm', se=F) +
    annotate("text", x = 10^1, y = 10^3, label = bquote(.(length(unique(bankfullValidationSet_Qb$site_no)))~'gages'), size=7)+
    scale_x_log10(labels = scales::label_log(digits = 2))+#, limits=c(10^1, 10^8))+
    scale_y_log10(labels = scales::label_log(digits = 2))+#, limits=c(10^1, 10^8))+
    xlab(bquote(bold('In situ bankfull discharge ['*m^3~s^-1*']')))+
    ylab(bquote(bold('Estimated bankfull discharge ['*m^3~s^-1*']')))+
    theme(axis.title=element_text(face='bold', size=18),
          axis.text = element_text(size=15),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size=15))
  
  plot_Hb <- ggplot(bankfullValidationSet_Hb, aes(x=Hb_m, y=bankfull_depth_mu_m)) +
    geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
    geom_point(size=5) +
    geom_smooth(method='lm', se=F) +
    annotate("text", x = 0.75, y = 3, label = bquote(.(length(unique(bankfullValidationSet_Hb$site_no)))~'gages'), size=7)+
 #   scale_x_log10(labels = scales::label_log(digits = 2))+#, limits=c(10^1, 10^8))+
#    scale_y_log10(labels = scales::label_log(digits = 2))+#, limits=c(10^1, 10^8))+
    xlab(bquote(bold('In situ bankfull depth ['*m*']')))+
    ylab(bquote(bold('Estimated bankfull depth ['*m*']')))+
    theme(axis.title=element_text(face='bold', size=18),
          axis.text = element_text(size=15),
          legend.position = c(0.8, 0.2),
          legend.text = element_text(size=15))
  
  design <- "
  AB
  "
  
  comboPlot <- patchwork::wrap_plots(A=plot_Qb, B=plot_Hb, design = design)
  
  ggsave('cache/bankfullValidation.png', comboPlot, width=12, height=7)
  
  
  return(list('Qb'=bankfullValidationSet_Qb,
              'Hb'=bankfullValidationSet_Hb))
}









buildEventAnalysis <- function(path_to_data, huc4, gagedata, gage, bankfullModel, BHGmodel, minDuration){
  bankfull_stage_m <- bankfullModel$stage_m
  bankfull_depth_m <- bankfullModel$depth_m
  
  #if gageRecord is empty, just skip
  if(nrow(gagedata) == 0 | is.na(bankfull_stage_m) | is.na(gage$DA_skm) | gage$DA_skm == 0) { #return empty if no gagerecord, NA stage (because of missing rating curve) or NA drainage area (because missing from gage metadata)
    return(data.frame())
  }
  
  #prep
  gageID <- gage$site_no

  #calculate bankfull discharge from the published shift-adjusted rating tables
  ratingTable <- dataRetrieval::readNWISrating(gageID, "exsa")
  ratingTable$INDEP <- ratingTable$INDEP# + ratingTable$SHIFT
  
  #only keep flowing events
  ratingTable$stage_m <- ratingTable$INDEP * 0.3048# + ratingTable$SHIFT #ft to m
  ratingTable$Q_cms <- ratingTable$DEP * 0.0283
  ratingTable <- ratingTable %>%
    dplyr::filter(stage_m > 0 & Q_cms > 0)
  
  #flood event analysis
  bankfull_discharge_cms <- (ratingTable[which(ratingTable$INDEP == round(bankfull_stage_m*3.28084,2)),]$DEP)*0.0283

  # get flood event duration from flow timeseries (workaround to handle known date class bug with midnight- https://github.com/tidyverse/lubridate/issues/1124)
  gagedata$date <- lubridate::ymd_hms(format(as.POSIXct(gagedata$dateTime), format = "%Y-%m-%d %T %Z"))
  
  floods_mc <- gagedata %>%
    dplyr::mutate(S_m = round(stage_m,3)) %>%
    dplyr::select(c('site_no', 'date', 'S_m', 'Q_cms', 'exceed_prob')) %>%
    dplyr::mutate(Sb_m = round(bankfull_stage_m,3),
                  Qb_cms = bankfull_discharge_cms) %>%
    dplyr::filter(S_m > Sb_m) %>% #both are rounded to 2 decimals (~1/2 an inch) to avoid over precision when thresholding floods
    tidyr::drop_na() %>% #handle NA dates in the flow record (errors in USGS data)
    dplyr::mutate(act_next = dplyr::lead(date)) %>%
    tidyr::drop_na() #remove final timestep, where 'next' is NA
  
  if(nrow(floods_mc)==0){return(data.frame())} #if no floods, move on
  
  #calculate the 'expected' next sample timestamp, given the full flow record
  floods_mc$index <- sapply(floods_mc$date, function(x){return(which(gagedata$date == x)+1)})
  floods_mc$exp_next <- gagedata[floods_mc$index,]$date
  
  #loop through flow record, finding starting and ending points of flood events
  floods_mc$event_id <- 0
  ticker <- 1
  for(i in 1:nrow(floods_mc)) { #already sorted by date
    floods_mc[i,]$event_id <- ifelse(floods_mc[i,]$exp_next == floods_mc[i,]$act_next, ticker, 0)
    
    if(floods_mc[i,]$event_id == 0) {
      floods_mc[i,]$event_id <- ticker
      ticker <- ticker + 1
      next
    }
    else {
      next
    }
  }
  
  floods_mc$event_id <- as.character(floods_mc$event_id)
  
  floods_mc <- floods_mc %>%
    dplyr::select(c('site_no', 'date', 'Sb_m', 'Qb_cms', 'S_m', 'Q_cms', 'exceed_prob', 'event_id'))
  
  #get physiographic region (and thus bankful depth)
  floods_mc$Hb_m <- bankfull_depth_m
  
  floods_mc$a_Wb <- BHGmodel[BHGmodel$division == gage$physio_region,]$a_Wb
  floods_mc$b_Wb <- BHGmodel[BHGmodel$division == gage$physio_region,]$b_Wb
  floods_mc$Wb_m <- floods_mc$a_Wb * (gage$DA_skm)^floods_mc$b_Wb
  
  # calculate exchange metrics for overbank flood events
  event_metrics <- floods_mc %>%
    dplyr::mutate(date=lubridate::ymd_hms(date))%>%
    dplyr::group_by(event_id) %>%
    dplyr::summarise(date_start = min(date),
                     date_end = max(date),
                     duration_min = lubridate::interval(min(date),max(date)) %/% lubridate::minutes(1),
                     spread = mean(spread, nna.rm=T), #mean to pass constant through
                     Hb_m = mean(Hb_m, na.rm=T), #mean to pass constant through
                     Wb_m = mean(Wb_m, na.rm=T), #mean to pass constant through
                     Hf_eventmean_m = mean(S_m - Sb_m, na.rm=T),
                     S_eventmean_m = mean(S_m, na.rm=T),
                     Q_eventmean_cms = mean(Q_cms, na.rm=T),
                     Q_eventpeak_cms = max(Q_cms, na.rm=T),
                     exceed_prob_eventpeak = min(exceed_prob, na.rm=T),
                     exceed_prob_eventmean = mean(exceed_prob, na.rm=T)) %>% #exceedance prob associated with the peak discharge (i.e. lowest prob)
    dplyr::mutate(Htf_eventmean_m = Hb_m + Hf_eventmean_m,
                  site_no = gage$site_no) %>%
    dplyr::relocate(site_no) %>%
    dplyr::filter((duration_min/60) > minDuration) #remove 'events' that are too short (see input parameter)

  event_metrics$event_id <- as.numeric(event_metrics$event_id)

  event_metrics <- event_metrics[order(event_metrics$event_id),]
  
  return(event_metrics)
}







buildEventAnalysis_full <- function(event_metrics, gages) {
  gages_fin <- dplyr::bind_rows(gages) %>%
    sf::st_drop_geometry()
  
  #summarise across monte carlo simulations
  event_metrics_fin <- dplyr::bind_rows(event_metrics) %>%
    dplyr::group_by(site_no, event_id) %>% #to summarise across MC simuls
    dplyr::summarise_at(c('Wb_m', 'Hb_m', 'date_start', 'duration_min','exceed_prob_eventmean', 'Q_eventmean_cms', 'Htf_eventmean_m', 'Hf_eventmean_m'), list(mean=mean, sd=sd)) %>%
    dplyr::rename(date_start = date_start_mean, #these cols don't vary across the MC simuls, so just recast back to static variables
                  duration_min = duration_min_mean,
                  Wb_m = Wb_m_mean,
                  Hb_m = Hb_m_mean) %>%
    dplyr::select(!c('duration_min_sd', 'Wb_m_sd', 'Hb_m_sd')) %>%
    dplyr::left_join(gages_fin, by='site_no') %>%
    dplyr::relocate(site_no, physio_region, DA_skm, Wb_m, Hb_m) %>%
    dplyr::rename(region=physio_region)
  
  return(event_metrics_fin)
}






buildUpscalingModel <- function(huc4, floodEvents, gageRecordStart, gageRecordEnd){
  huc2 <- substr(huc4, 1, 2)
  
  #get duration of entire gage record
  duration_dys <- lubridate::ymd(gageRecordEnd) - lubridate::ymd(gageRecordStart)
  
  #get mean flood depths across flood events
  forModel_mean <- floodEvents %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(meanflood_Htf_m = mean(Htf_eventmean_m_mean, na.rm=T),
                     DA_skm = mean(DA_skm, na.rm=T)) #mean to carry through summarise
  
  #Q0.5 flood
  q05_forModel <- floodEvents %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_eventmean_mean - 0.005)) %>% 
    dplyr::slice_min(diff) %>%
    dplyr::ungroup() %>%
    dplyr::filter(diff < 0.01) %>% #make sure the floods are actually in the record, i.e. flow must be within 1% of expected exceedance probability
    dplyr::mutate(q05flood_Htf_m = Htf_eventmean_m_mean) %>%
    dplyr::select(c('site_no','q05flood_Htf_m'))


  #Q1 flood
  q1_forModel <- floodEvents %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_eventmean_mean - 0.01)/0.01) %>% 
    dplyr::slice_min(diff) %>%
    dplyr::ungroup() %>%
    dplyr::filter(diff < 0.20) %>% #make sure the floods are actually in the record, i.e. flow must be within 1% of expected exceedance probability
    dplyr::mutate(q1flood_Htf_m = Htf_eventmean_m_mean) %>%
    dplyr::select(c('site_no','q1flood_Htf_m'))
  
  #Q5 flood
  q5_forModel <- floodEvents %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_eventmean_mean - 0.05)/0.05) %>% 
    dplyr::slice_min(diff) %>%
    dplyr::ungroup() %>%
    dplyr::filter(diff < 0.20) %>% #make sure the floods are actually in the record, i.e. flow must be within 1% of expected exceedance probability
    dplyr::mutate(q5flood_Htf_m = Htf_eventmean_m_mean) %>%
    dplyr::select(c('site_no','q5flood_Htf_m'))

  #Q10 flood
  q10_forModel <- floodEvents %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_eventmean_mean - 0.1)/0.1) %>% 
    dplyr::slice_min(diff) %>%
    dplyr::ungroup() %>%
    dplyr::filter(diff < 0.2) %>% #make sure the floods are actually in the record, i.e. flow must be within 1% of expected exceedance probability
    dplyr::mutate(q10flood_Htf_m = Htf_eventmean_m_mean) %>%
    dplyr::select(c('site_no','q10flood_Htf_m'))

  #Q25 flood
  q25_forModel <- floodEvents %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_eventmean_mean - 0.25)/0.25) %>% 
    dplyr::slice_min(diff) %>%
    dplyr::ungroup() %>%
    dplyr::filter(diff < 0.20) %>% #make sure the floods are actually in the record, i.e. flow must be within 1% of expected exceedance probability
    dplyr::mutate(q25flood_Htf_m = Htf_eventmean_m_mean) %>%
    dplyr::select(c('site_no','q25flood_Htf_m'))

  forModel <- forModel_mean %>%
    dplyr::left_join(q05_forModel, by='site_no') %>%
    dplyr::left_join(q1_forModel, by='site_no') %>%
    dplyr::left_join(q5_forModel, by='site_no') %>%
    dplyr::left_join(q10_forModel, by='site_no') %>%
    dplyr::left_join(q25_forModel, by='site_no') %>%
    tidyr::pivot_longer(c('meanflood_Htf_m', 'q05flood_Htf_m', 'q1flood_Htf_m', 'q5flood_Htf_m', 'q10flood_Htf_m', 'q25flood_Htf_m'))

  #fit Htf models
  models_r2 <- forModel %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(broom::glance(lm(log(value) ~ log(DA_skm)))) %>%
    dplyr::summarise(n_gages=sum(!(is.na(value))),
                    rsq = mean(adj.r.squared, na.rm=T)) %>% #mean to pass the constant through
    dplyr::select(c(name, n_gages, rsq))

  models <- forModel %>%
    dplyr::group_by(name) %>%
    dplyr::group_modify(~ broom::tidy(lm(log(value) ~ log(DA_skm), data = .x))) %>%
    dplyr::select(name:std.error) %>%
    dplyr::left_join(models_r2, by='name') %>%
    tidyr::pivot_wider(names_from=term, values_from=c(estimate, std.error, n_gages, rsq)) %>%
    dplyr::select(c('name', 'estimate_(Intercept)','estimate_log(DA_skm)', 'std.error_(Intercept)', 'std.error_log(DA_skm)', 'n_gages_(Intercept)', 'rsq_(Intercept)'))
  
  colnames(models) <- c('name', 'coef','exp', 'stderr_coef', 'stderr_exp', 'n_gages', 'rsq')
  
  theme_set(theme_classic())
  
  plot <- ggplot(forModel, aes(x=DA_skm, y=value)) +
    geom_point(size=5, color='#153131')+
    geom_smooth(method='lm', se=F, color='black',linewidth=1.25) +
    scale_x_log10()+
    scale_y_log10() +
    ylab(bquote(bold("Flood level ["*m*"]"))) +
    xlab(bquote(bold("Drainage Area ["*km^2*"]"))) +
    theme(axis.title = element_text(face = 'bold', size=18),
          axis.text = element_text(size=15),
          plot.title = element_text(face = 'bold', size=22)) +
    facet_wrap(~factor(name, levels=c('meanflood_Htf_m', 'q05flood_Htf_m', 'q1flood_Htf_m', 'q5flood_Htf_m', 'q10flood_Htf_m', 'q25flood_Htf_m'),
                          labels=c('Mean flood', '0.5% flow', '1% flow', '5% flow', '10% flow', '25% flow'))) +
    geom_text(data = models_r2, aes(x = 10^1.5, y = 5, label = paste0('rsq: ',round(rsq,2))), size=6) +
    geom_text(data = models_r2, aes(x = 10^1.5, y = 3.5, label = paste0(n_gages, ' gages')), size=6) +
    theme(strip.text.x = element_text(size = 20))
  
  ggsave(paste0('cache/upscalingModel.png'), plot, width=12, height=12)
  
  return(models)
}








buildNetworkModel <- function(path_to_data, huc4, upscalingModel, BHGmodel, basinData) {
  huc2 <- substr(huc4, 1, 2)
  
  network_VAA <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
  network_gages <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
  network <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
    sf::st_zm() %>%
    dplyr::left_join(network_VAA, by='NHDPlusID') %>%
    dplyr::left_join(network_gages, by='NHDPlusID') %>%
    dplyr::select(c('NHDPlusID', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))
  
  huc4id <- huc4
  basin <- sf::st_read(paste0(path_to_data, '/HUC2_', huc2, '/WBD_', huc2, '_HU2_Shape/Shape/WBDHU4.shp')) %>%
    dplyr::filter(huc4 == huc4id)
  basin <- fixGeometries(basin)

  #setup bankfull geometry
  sf::sf_use_s2(FALSE)
  regions <- sf::st_read('data/physio.shp') #physiographic regions
  regions <- fixGeometries(regions)
  shp_temp <- sf::st_join(basin, regions, largest=TRUE) #take the physiographic region that the basin is mostly in (dominant spatial intersection)

  physio_region <- shp_temp$DIVISION
  
  network$a_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Wb)})
  network$b_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Wb)})
  network$Wb_m <- network$a_Wb * (network$TotDASqKm)^network$b_Wb

  #upscale mean flood depths by drainage area
  network$Htf_avg_m <- exp(upscalingModel[upscalingModel$name == 'meanflood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'meanflood_Htf_m',]$exp
  network$Htf_q05_m <- exp(upscalingModel[upscalingModel$name == 'q05flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q05flood_Htf_m',]$exp
  network$Htf_q1_m <- exp(upscalingModel[upscalingModel$name == 'q1flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q1flood_Htf_m',]$exp
  network$Htf_q5_m <- exp(upscalingModel[upscalingModel$name == 'q5flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q5flood_Htf_m',]$exp
  network$Htf_q10_m <- exp(upscalingModel[upscalingModel$name == 'q10flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q10flood_Htf_m',]$exp
  network$Htf_q25_m <- exp(upscalingModel[upscalingModel$name == 'q25flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q25flood_Htf_m',]$exp

  network <- network %>%
    dplyr::filter(AreaSqKm > 0) %>%
    dplyr::relocate(Shape, .after=tidyselect::last_col())
  
  #add reach type
  lookup <- basinData$waterbodyLookUp %>%
    dplyr::filter(!(is.na(NHDPlusID_reach))) %>%
    dplyr::select(c('NHDPlusID_reach', 'waterbody_type'))

  network <- network %>%
    dplyr::left_join(lookup, by=c('NHDPlusID'='NHDPlusID_reach'))
  network$waterbody_type <- ifelse(is.na(network$waterbody_type), 'floodplain',
                                ifelse(network$waterbody_type == '466', 'wetland',
                                    ifelse(network$waterbody_type == '436', 'reservoir',
                                         ifelse(network$waterbody_type == '390', 'lake', 'other'))))
  
  return(network)
}




runNetworkModel <- function(path_to_data, huc4, basinData, network){
  #setup basin shapefiles
  dem <- terra::unwrap(basinData$terra$dem) #[cm] need to unwrap the packed terra package (for distributed computing)
  d8 <- terra::unwrap(basinData$terra$d8) #[cm] need to unwrap the packed terra package (for distributed computing)

  #run embarresingly parallel inundation model
  inundationWrapper <- function(reach) {
    #lower bound on dem resolution (also, these small rivers don't reallllly have floodplains, right?)
    if(reach$Wb_m <= 10 | reach$waterbody_type == 'lake'){
      reach[,c('Af_avg_m2', 'Af_q05_m2', 'Af_q1_m2', 'Af_q5_m2', 'Af_q10_m2', 'Af_q25_m2')] <-NA
      reach[,c('Vol_avg_m3', 'Vol_q05_m3', 'Vol_q1_m3', 'Vol_q5_m3', 'Vol_q10_m3', 'Vol_q25_m3')] <- NA
      return(reach)
    }
    
    #run inundation model
    inundationDataPackage <- prepInundationData(path_to_data, huc4, reach$NHDPlusID, reach$Wb_m, dem, d8)
    if(inundationDataPackage$flag == 'no pour points'){
      reach[,c('Af_avg_m2', 'Af_q05_m2', 'Af_q1_m2', 'Af_q5_m2', 'Af_q10_m2', 'Af_q25_m2')] <- NA
      reach[,c('Vol_avg_m3', 'Vol_q05_m3', 'Vol_q1_m3', 'Vol_q5_m3', 'Vol_q10_m3', 'Vol_q25_m3')] <- NA
      return(reach)
    }

    #loop through flood sizes and model inundation
    areaWrapper <- function(depth){
      floodMap <- modelInundation(inundationDataPackage, depth)

      #convert map to flood area
      flooded_pixels <- terra::freq(floodMap$binary)
      flood_area_m2 <- flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
      flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

      return(flood_area_m2)
    }

    volWrapper <- function(depth){
      floodMap <- modelInundation(inundationDataPackage, depth)

      #convert map to flood area
      flooded_pixels <- terra::freq(floodMap$binary)
      flood_area_m2 <- flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
      flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

      #convert flood map to volume
      if(flood_area_m2 > 0) { #only if there's any actual inundated area at 10m resolution
        floodmask <- floodMap$binary
        floodmask[floodmask != 1] = 0
        flood_depths <- floodMap$depths * floodmask
        flood_depths <- terra::global(flood_depths, fun="sum", na.rm=T)
        flood_vol_m3 <- flood_depths * flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
        flood_vol_m3 <- flood_vol_m3[[1]]
      }
      else{
        flood_vol_m3 <- 0
      }
    
    return(flood_vol_m3)
    }
    reach[,c('Af_avg_m2', 'Af_q05_m2', 'Af_q1_m2', 'Af_q5_m2', 'Af_q10_m2', 'Af_q25_m2')] <- lapply(reach[,c('Htf_avg_m', 'Htf_q05_m', 'Htf_q1_m', 'Htf_q5_m', 'Htf_q10_m', 'Htf_q25_m')], areaWrapper)
    reach[,c('Vol_avg_m3', 'Vol_q05_m3', 'Vol_q1_m3', 'Vol_q5_m3', 'Vol_q10_m3', 'Vol_q25_m3')] <- lapply(reach[,c('Htf_avg_m', 'Htf_q05_m', 'Htf_q1_m', 'Htf_q5_m', 'Htf_q10_m', 'Htf_q25_m')], volWrapper)

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




hortonScaling <- function(basinAnalysis){

  scaledBasin <- basinAnalysis_0108 %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(StreamCalc > 0)  %>%  #RERUN WHOLE ANALYSIS WITH THIS REMOVED
    dplyr::group_by(StreamCalc) %>% 
    dplyr::summarise(Af_by_order = sum(Af_avg_m2, na.rm=T),
                    Vol_by_order = sum(Vol_avg_m3, na.rm=T),
                    frac = round((sum(Wb_m > 10, na.rm=T)/n()),2),
                    lengthKM_by_order = sum(LengthKM, na.rm=T)) %>%
    dplyr::mutate(frac = ifelse(frac == 0, NA, frac))
  
  scaledBasin$Af_by_order_scaled <- scaledBasin$Af_by_order * 1/scaledBasin$frac #for orders without 100% of rivers > 10m, this number neds to be inflated by the fraction of the order that is 'visible', i.e. Wb > 10m
  scaledBasin$Vol_by_order_scaled <- scaledBasin$Vol_by_order * 1/scaledBasin$frac #for orders without 100% of rivers > 10m, this number neds to be inflated by the fraction of the order that is 'visible', i.e. Wb > 10m

  #fit Horton equations
  area_model <- lm(log(Af_by_order_scaled)~log(StreamCalc), data=scaledBasin)
  vol_model <- lm(log(Vol_by_order_scaled)~log(StreamCalc), data=scaledBasin)
  length_model <- lm(log(lengthKM_by_order)~log(StreamCalc), data=scaledBasin)

  #use horton eqs to calculate inundated areas and volumes for orders with frac < 100%
  scaledBasin$Af_by_order_km2_fin <- exp(predict(area_model, scaledBasin)) * 1e-6 #[m2 to km2]
  scaledBasin$frac <- ifelse(is.na(scaledBasin$frac), 0, scaledBasin$frac)
  scaledBasin$Af_by_order_km2_fin <- ifelse(scaledBasin$frac == 1, scaledBasin$Af_by_order *1e-6, scaledBasin$Af_by_order_km2_fin)

  scaledBasin$Vol_by_order_km3_fin <- exp(predict(vol_model, scaledBasin)) * 1e-9 #[m3 to km3]
  scaledBasin$frac <- ifelse(is.na(scaledBasin$frac), 0, scaledBasin$frac)
  scaledBasin$Vol_by_order_km3_fin <- ifelse(scaledBasin$frac == 1, scaledBasin$Vol_by_order *1e-9, scaledBasin$Vol_by_order_km3_fin)

  #prep output
  scaledBasin <- scaledBasin %>%
    dplyr::select(c('StreamCalc', 'Af_by_order_km2_fin', 'Vol_by_order_km3_fin'))

  return(list('scaling_results'=scaledBasin,
              'area_model_summary'=summary(area_model),
              'vol_model_summary'=summary(vol_model)))
}