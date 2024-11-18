## Functions for floodplain model
## Craig Brinkerhoff
## Summer 2024


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
                                          service='uv', #15' or 30' stage and flow (depending on gage)
                                          parameterCd = c('00060'), #discharge parameter code [cfs]
                                          startDate = gageRecordStart,
                                          endDate = gageRecordEnd),
                      error = function(m){return(data.frame())})

  if(nrow(gagedata)==0) {return(data.frame())} #if no gage data, move on
  if(!("X_00060_00000" %in% colnames(gagedata))){return(data.frame())} #if incorrect discharge column name (happens sometimes), move on

  #convert to metric
  gagedata$Q_cms <- gagedata$X_00060_00000 * 0.0283 #cfs to cms
  
  #only keep flowing events
  gagedata <- gagedata %>%
    dplyr::filter(Q_cms > 0)

  # convert to date (workaround to handle known date class bug with midnight- https://github.com/tidyverse/lubridate/issues/1124)
  gagedata$date <- lubridate::ymd_hms(format(as.POSIXct(gagedata$dateTime), format = "%Y-%m-%d %T %Z"))

  gagedata_og <- gagedata

  #downsample to hourly (just to be consistent and reduce volume of data)
  gagedata <- gagedata_og %>%
    dplyr::mutate(date_hr = lubridate::ymd_h(paste0(lubridate::year(date), '-', lubridate::month(date),'-',lubridate::day(date), '-', lubridate::hour(date)))) %>%
    dplyr::group_by(date_hr) %>%
    dplyr::summarise(site_no = dplyr::first(site_no), #pass through the group by
                    Q_cms = mean(Q_cms, na.rm=T))

  #get max annual floods for calculating the FEMA 100yr AEP
  gagedata_aep <- gagedata_og %>%
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
  ratingTable <- tryCatch(dataRetrieval::readNWISrating(gageID, "exsa"),
                      error = function(m){return(data.frame())})
  ratingTable$stage_m <- (ratingTable$INDEP  + ratingTable$SHIFT) * 0.3048 #ft to m
  ratingTable$Q_cms <- ratingTable$DEP * 0.0283

  #convert historical streamflow record to stages using rating table
  gagedata$stage_m <- sapply(gagedata$Q_cms, function(i){ratingTable[which.min(abs(ratingTable$Q_cms - i)),]$stage_m})

  if(is.list(gagedata$stage_m)){return(data.frame())} #sometimes you get an empty list, presumbaly b/c the rating table is empty

  #add flag if gagedata is outside of the rating table
  gagedata$ratingflag <- ifelse(gagedata$Q_cms > max(ratingTable$Q_cms) | gagedata$Q_cms < min(ratingTable$Q_cms), 1, 0)

  #add fema 100 yr aep
  fema100yrflood_cms <- gagedata_aep[which.min(abs(gagedata_aep$exceed_prob - 0.01)),]$Q_cms
  bankfull_1_5yr_cms <- gagedata_aep[which.min(abs(gagedata_aep$exceed_prob - 0.667)),]$Q_cms

  gagedata$fema100yrflood_cms <- fema100yrflood_cms
  gagedata$bankfull_1_5yr_cms <- bankfull_1_5yr_cms

  return(gagedata)
}




prepInundationData <- function(huc4, reachID, bankfullWidth, dem, d8){
  #prep
  huc2 <- substr(huc4, 1, 2)
  
  sql <- paste0('SELECT * FROM NHDPlusCatchment WHERE NHDPlusID = ', reachID)
  gage_unit_catchment <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),
                                    layer = 'NHDPlusCatchment',
                                    query = sql, #sql query to only load the catchment we want efficiently (see above)
                                    quiet = TRUE)
  
  gage_unit_catchment <- terra::vect(gage_unit_catchment)
  gage_unit_catchment <- terra::project(gage_unit_catchment, terra::crs(dem))

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
  # flowline <- terra::vect(flowline)
  # bankfullMask <- terra::buffer(flowline, bankfullWidth) #don't let it go below the spatial resolution of the dem
  # r <- terra::rast(bankfullMask, resolution=10, extent=terra::ext(hand))
  # bankfullMask <- terra::rasterize(bankfullMask, r)
  # bankfullMask[is.na(bankfullMask[])] <- 0
  # bankfullMask <- terra::subst(bankfullMask, 1, 2) #reacast mask as value '2'
  
  # reachWBAreaPermanentIdentifier <- flowline$WBArea_Permanent_Identifier
  
  #return what's needed for inundation modeling
  return(list('hand'=hand,
            #  'bankfullMask'=bankfullMask,
              'flag'='has pour points'))
}




modelInundation <- function(inundationData, floodDepth){
  hand <- inundationData$hand
  #bankfullMask <- inundationData$bankfullMask
  
  #actual inundation calculation
  floodDepths <- floodDepth - hand
  floodMap <- floodDepths
  floodMap[floodMap < 0] <- 0
  floodMap[floodMap > 0] <- 1
  #floodMap <- floodMap + bankfullMask #add bankfull channel mask to flood map

  #floodMap <- terra::subst(floodMap, 3, 2) #recast 'bankfullMask + flooded' as river, to copacetic with 'river only'
  
  floodDepths[floodDepths < 0] <- 0

  return(list('binary'=floodMap,
              'depths'=floodDepths))
}






# buildStageDepthRelationship <- function(gage, minADCPMeas){
#   #prep
#   gageID <- gage$site_no

#   #prep custom rating table from available adcp measurements
#   surfaceData <- tryCatch(dataRetrieval::readNWISmeas(gageID, expanded = TRUE),
#                     error=function(m){return(data.frame('site_no'=gageID,
#                                                       'lm_r2'=NA,
#                                                       'lm_depth_a'=NA,
#                                                       'lm_depth_b'=NA))})
  
#   if(nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
#                                                       'lm_r2'=NA,
#                                                       'lm_depth_a'=NA,
#                                                       'lm_depth_b'=NA))}
#   surfaceData <- surfaceData %>%
#     dplyr::filter(measured_rating_diff %in% c('Good', 'Excellent')) %>%
#     dplyr::mutate('stage_m'=gage_height_va * 0.3048,
#                   'Q_cms'=chan_discharge* 0.0283,
#                   'width_m'=chan_width*0.3048,
#                   'area_m2'=chan_area*0.092903,
#                   'velocity_ms'=chan_velocity * 0.3048) %>%
#     dplyr::select(c('site_no', 'stage_m', 'Q_cms', 'width_m', 'area_m2', 'velocity_ms')) %>%
#     dplyr::mutate(depth_m=Q_cms / (width_m*velocity_ms),
#                   log10_Q_cms = log10(Q_cms)) %>%
#     dplyr::filter(stage_m > 0 & Q_cms > 0)

#   #remove likely erronous outlier depth measurements
#   iqr <- IQR(surfaceData$depth_m, na.rm=T)
  
#   surfaceData <- surfaceData %>%
#     dplyr::filter(depth_m < (quantile(depth_m, 0.75, na.rm=T) + 1.5*iqr)) %>%
#     dplyr::filter(depth_m > (quantile(depth_m, 0.25, na.rm=T) - 1.5*iqr))
  
#   if(nrow(surfaceData)< minADCPMeas){return(data.frame('site_no'=gageID,
#                                                       'lm_r2'=NA,
#                                                       'lm_depth_a'=NA,
#                                                       'lm_depth_b'=NA))}
  
#   #get bankfull depth
#   lm_depth <- lm((stage_m)~(depth_m), data=surfaceData)
#   lm_depth_a <- coef(lm_depth)[1]
#   lm_depth_b <- coef(lm_depth)[2]
  
#   #setup site info for later export
#   site_info <- data.frame('site_no'=gageID,
#                           'lm_r2'=summary(lm_depth)$r.squared,
#                           'lm_depth_a'=lm_depth_a,
#                           'lm_depth_b'=lm_depth_b)

#   rownames(site_info) <- NULL
  
#   return(site_info)
# }




buildDepthAHG <- function(gage, minADCPMeas){
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
  lm_depth <- lm(log(depth_m)~log(Q_cms), data=surfaceData)
  lm_depth_a <- coef(lm_depth)[1]
  lm_depth_b <- coef(lm_depth)[2]
  
  #setup site info for later export
  site_info <- data.frame('site_no'=gageID,
                          'lm_r2'=summary(lm_depth)$r.squared,
                          'lm_depth_a'=lm_depth_a,
                          'lm_depth_b'=lm_depth_b)

  rownames(site_info) <- NULL
  
  return(site_info)
}






buildUpscalingModel <- function(huc4, gageRecord, gage, BHGmodel, depAHG, minAHGr2, gageRecordStart, gageRecordEnd){
  huc2 <- substr(huc4, 1, 2)

  #setup
  depAHG <- dplyr::bind_rows(depAHG)

  gageRecord <- dplyr::bind_rows(gageRecord) %>%
    dplyr::select(c('site_no', 'date_hr', 'fema100yrflood_cms', 'bankfull_1_5yr_cms', 'Q_cms', 'exceed_prob'))

  #filter for flows above bankfull, to calculate exceedance probabilities for flood events
  gage <- sf::st_drop_geometry(gage) %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(c('site_no', 'DA_skm', 'physio_region')) %>% 
    sf::st_drop_geometry()
  
  gage$a_Qb <- sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Qb)})
  gage$b_Qb <- sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Qb)})
  gage$Qb_m <- gage$a_Qb * (gage$DA_skm)^gage$b_Qb
  
  #bring everything together
  upscaling_df <- gageRecord %>%
    dplyr::left_join(gage, by='site_no') %>%
    dplyr::left_join(depAHG, by='site_no') %>%
    dplyr::filter(!(is.na(DA_skm))) %>%
    dplyr::filter(lm_r2 > minAHGr2)

  #filter for flows above bankfull and calculate exceed prob from that
  num <- upscaling_df %>%
    dplyr::filter(Q_cms > Qb_m) %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(n=n())
  
  exceed_prob_df <- upscaling_df %>%
    dplyr::left_join(num, by='site_no') %>%
    dplyr::filter(Q_cms > Qb_m) %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(rank = rank(-Q_cms)) %>%
    dplyr::mutate(exceed_prob_bankfull = rank / (n+1)) %>%
    dplyr::select(c('site_no', 'Q_cms','exceed_prob_bankfull'))

  upscaling_df <- upscaling_df %>%
    dplyr::left_join(exceed_prob_df, by=c('site_no', 'Q_cms')) %>%
    dplyr::filter(!(is.na(exceed_prob_bankfull)))

  #get duration of entire gage record
  duration_dys <- lubridate::ymd(gageRecordEnd) - lubridate::ymd(gageRecordStart)

  #calculate water level, relative to channel bottom
  upscaling_df$Htf_m <- exp(upscaling_df$lm_depth_a)*upscaling_df$Q_cms^upscaling_df$lm_depth_b
  upscaling_df$femaHtf_m <- exp(upscaling_df$lm_depth_a)*upscaling_df$fema100yrflood_cms^upscaling_df$lm_depth_b

  #100yr FEMA AEP flood
  qFEMA_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(qFEMAflood_Htf_m = dplyr::first(femaHtf_m)) %>% #it's a constant so
    dplyr::select(c('site_no', 'qFEMAflood_Htf_m'))

  #Q0.1 flood
  q01_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.001)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q01flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no', 'q01flood_Htf_m'))

  #Q1 flood
  q1_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.01)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q1flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q1flood_Htf_m'))
  
  #Q2.5 flood
  q205_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.025)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q205flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q205flood_Htf_m'))
  
  #Q5 flood
  q5_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.05)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q5flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q5flood_Htf_m'))

  #Q10 flood
  q10_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.1)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q10flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q10flood_Htf_m'))

  #Q15 flood
  q15_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.15)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q15flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q15flood_Htf_m'))

  #Q20 flood
  q20_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.20)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q20flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q20flood_Htf_m'))
  
  #Q25 flood
  q25_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.25)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(q25flood_Htf_m = Htf_m) %>%
    dplyr::select(c('site_no','q25flood_Htf_m'))

  forModel <- gage %>%
    dplyr::left_join(qFEMA_forModel, by='site_no') %>%
    dplyr::left_join(q01_forModel, by='site_no') %>%
    dplyr::left_join(q1_forModel, by='site_no') %>%
    dplyr::left_join(q205_forModel, by='site_no') %>%
    dplyr::left_join(q5_forModel, by='site_no') %>%
    dplyr::left_join(q10_forModel, by='site_no') %>%
    dplyr::left_join(q15_forModel, by='site_no') %>%
    dplyr::left_join(q20_forModel, by='site_no') %>%
    dplyr::left_join(q25_forModel, by='site_no') %>%
    tidyr::pivot_longer(c('qFEMAflood_Htf_m', 'q01flood_Htf_m', 'q1flood_Htf_m', 'q205flood_Htf_m', 'q5flood_Htf_m', 'q10flood_Htf_m', 'q15flood_Htf_m', 'q20flood_Htf_m', 'q25flood_Htf_m'))

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
    ylab(bquote(bold("Water level ["*m*"]"))) +
    xlab(bquote(bold("Drainage Area ["*km^2*"]"))) +
    theme(axis.title = element_text(face = 'bold', size=18),
          axis.text = element_text(size=15),
          plot.title = element_text(face = 'bold', size=22)) +
    facet_wrap(~factor(name, levels=c('qFEMAflood_Htf_m', 'q01flood_Htf_m', 'q1flood_Htf_m', 'q205flood_Htf_m', 'q5flood_Htf_m', 'q10flood_Htf_m', 'q15flood_Htf_m', 'q20flood_Htf_m', 'q25flood_Htf_m'),
                          labels=c('100yr AEP flood','0.1% flood', '1% flood', '2.5% flood','5% flood', '10% flood', '15% flood','20% flood', '25% flood'))) +
    geom_text(data = models_r2, aes(x = 10^0.8, y = 5, label = paste0('rsq: ',round(rsq,2))), size=6) +
    geom_text(data = models_r2, aes(x = 10^0.8, y = 3.5, label = paste0(n_gages, ' gages')), size=6) +
    theme(strip.text.x = element_text(size = 20))
  
  ggsave(paste0('cache/upscalingModel.png'), plot, width=12, height=12)
  
  return(models)
}





buildNetworkModel <- function(huc4, upscalingModel, BHGmodel, basinData) {
  huc2 <- substr(huc4, 1, 2)
  
  network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
  network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
  network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
    sf::st_zm() %>%
    dplyr::left_join(network_VAA, by='NHDPlusID') %>%
    dplyr::left_join(network_gages, by='NHDPlusID') %>%
    dplyr::select(c('NHDPlusID', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))
  
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
  
  network$a_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Wb)})
  network$b_Wb <- sapply(physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Wb)})
  network$Wb_m <- network$a_Wb * (network$TotDASqKm)^network$b_Wb

  #upscale mean flood depths by drainage area
  #network$Htf_qFEMA_m <- exp(upscalingModel[upscalingModel$name == 'qFEMAflood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'qFEMAflood_Htf_m',]$exp
  # network$Htf_q01_m <- exp(upscalingModel[upscalingModel$name == 'q01flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q01flood_Htf_m',]$exp
  network$Htf_q1_m <- exp(upscalingModel[upscalingModel$name == 'q1flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q1flood_Htf_m',]$exp
  # network$Htf_q205_m <- exp(upscalingModel[upscalingModel$name == 'q205flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q205flood_Htf_m',]$exp
  # network$Htf_q5_m <- exp(upscalingModel[upscalingModel$name == 'q5flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q5flood_Htf_m',]$exp
  # network$Htf_q10_m <- exp(upscalingModel[upscalingModel$name == 'q10flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q10flood_Htf_m',]$exp
  # network$Htf_q15_m <- exp(upscalingModel[upscalingModel$name == 'q15flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q15flood_Htf_m',]$exp
  # network$Htf_q20_m <- exp(upscalingModel[upscalingModel$name == 'q20flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q20flood_Htf_m',]$exp
  # network$Htf_q25_m <- exp(upscalingModel[upscalingModel$name == 'q25flood_Htf_m',]$coef) * network$TotDASqKm^upscalingModel[upscalingModel$name == 'q25flood_Htf_m',]$exp

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




runNetworkModel <- function(huc4, basinData, network){
  #setup basin shapefiles
  dem <- terra::unwrap(basinData$terra$dem) #[cm] need to unwrap the packed terra package (for distributed computing)
  d8 <- terra::unwrap(basinData$terra$d8) #[cm] need to unwrap the packed terra package (for distributed computing)

  #run embarresingly parallel inundation model
  inundationWrapper <- function(reach) {
    #lower bound on dem resolution (also, these small rivers don't reallllly have floodplains, right?)
    if(reach$Wb_m <= 10 | reach$waterbody_type == 'lake'){
      reach[,c('A_qFEMA_m2')] <-NA
      reach[,c('Vol_qFEMA_m3')] <- NA
      # reach[,c('Af_q01_m2', 'Af_q02_m2', 'Af_q1_m2', 'Af_q205_m2','Af_q5_m2', 'Af_q10_m2', 'Af_q15_m2', 'Af_q20_m2', 'Af_q25_m2')] <-NA
      # reach[,c('Vol_q01_m3', 'Vol_q02_m3', 'Vol_q1_m3', 'Vol_q205_m3', 'Vol_q5_m3', 'Vol_q10_m3', 'Vol_q15_m3', 'Vol_q20_m3', 'Vol_q25_m3')] <- NA
      return(reach)
    }
    
    #run inundation model
    inundationDataPackage <- prepInundationData(huc4, reach$NHDPlusID, reach$Wb_m, dem, d8)
    if(inundationDataPackage$flag == 'no pour points'){
      reach[,c('A_qFEMA_m2')] <-NA
      reach[,c('Vol_qFEMA_m3')] <- NA
      # reach[,c('Af_q01_m2', 'Af_q02_m2', 'Af_q1_m2', 'Af_q205_m2','Af_q5_m2', 'Af_q10_m2', 'Af_q15_m2', 'Af_q20_m2', 'Af_q25_m2')] <-NA
      # reach[,c('Vol_q01_m3', 'Vol_q02_m3', 'Vol_q1_m3', 'Vol_q205_m3', 'Vol_q5_m3', 'Vol_q10_m3', 'Vol_q15_m3', 'Vol_q20_m3', 'Vol_q25_m3')] <- NA
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
      # flooded_pixels <- terra::freq(floodMap$binary)
      # flood_area_m2 <- flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
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
      # flooded_pixels <- terra::freq(floodMap$binary)
      # flood_area_m2 <- flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
      flood_area_m2 <- ifelse(length(flood_area_m2)==0, 0, flood_area_m2)

      #convert flood map to volume
      if(flood_area_m2 > 0) { #only if there's any actual inundated area at 10m resolution
        floodmask <- floodMap$binary
        floodmask[floodmask != 1] = 0
        flood_depths <- floodMap$depths * floodmask
        flood_depths <- terra::global(flood_depths, fun="sum", na.rm=T)
        flood_vol_m3 <- flood_depths * flood_area_m2
        # flood_vol_m3 <- flood_depths * flooded_pixels[flooded_pixels$value == 1,]$count * 10^2 #10m pixels
        flood_vol_m3 <- flood_vol_m3[[1]]
      }
      else{
        flood_vol_m3 <- 0
      }
    
    return(flood_vol_m3)
    }
    reach[,c('A_qFEMA_m2')] <- lapply(reach[,c('Htf_q1_m')], areaWrapper)
    reach[,c('Vol_qFEMA_m3')] <- lapply(reach[,c('Htf_q1_m')], volWrapper)

    # reach[,c('Af_q01_m2', 'Af_q02_m2', 'Af_q1_m2', 'Af_q205_m2','Af_q5_m2', 'Af_q10_m2', 'Af_q15_m2', 'Af_q20_m2', 'Af_q25_m2')] <- lapply(reach[,c('Htf_q01_m', 'Htf_q02_m', 'Htf_q1_m', 'Htf_q205_m', 'Htf_q5_m', 'Htf_q10_m', 'Htf_q15_m', 'Htf_q20_m','Htf_q25_m')], areaWrapper)
    # reach[,c('Vol_q01_m3', 'Vol_q02_m3', 'Vol_q1_m3', 'Vol_q205_m3', 'Vol_q5_m3', 'Vol_q10_m3', 'Vol_q15_m3', 'Vol_q20_m3', 'Vol_q25_m3')] <- lapply(reach[,c('Htf_q01_m', 'Htf_q02_m', 'Htf_q1_m', 'Htf_q205_m', 'Htf_q5_m', 'Htf_q10_m', 'Htf_q15_m', 'Htf_q20_m','Htf_q25_m')], volWrapper)

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
  scaledBasin <- basinAnalysis %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(StreamCalc > 0)  %>%  #RERUN WHOLE ANALYSIS WITH THIS REMOVED
    dplyr::group_by(StreamCalc) %>% 
    dplyr::summarise(Af_by_order = sum(Af_qFEMA_m2, na.rm=T),
                    Vol_by_order = sum(Vol_qFEMA_m3, na.rm=T),
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