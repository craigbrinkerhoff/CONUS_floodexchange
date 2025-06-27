## Utility functions
## Craig Brinkerhoff
## Summer 2025





getBasinGages <- function(huc4id, gageRecordStart, gageRecordEnd){
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

  return(models)
}








dataBHG <- function(){
  dataset <- readr::read_csv('data/bhg_us_database_bieger_2015.csv') %>% #available by searching for paper at https://swat.tamu.edu/search
    dplyr::select(c('USGS Station No.', 'Physiographic Division', '...9', '...11','...13')) #some necessary manual munging for colnames from dataset
  
  colnames(dataset) <- c('GageID', 'DIVISION', 'DA_gage_skm','Qb_cms','Wb_m_obs')
  
  dataset$Qb_cms <- as.numeric(dataset$Qb_cms)
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



#Our bandaid version of readNWISrating that manually constructs our own exsa url
readNWISrating_CRAIG <- function(siteNumber, type='base', convertType = TRUE) {

  # No rating xml service
  #BROKEN AS OF 3/21/25. NWIS seems to have changed internal urls for rating tables, so dataRetrieval::readNWISrating() breaks because dataRetrieval::constructNWISURL() breaks. A band-aid is to directly constructing our own url using a custom version of readNWISrating() (in ~/src/utils.R), circumventing readNWISdata.
  #url <- constructNWISURL(siteNumber, service = "rating", ratingType = type)
  url <- paste0('https://waterdata.usgs.gov/nwisweb/data/ratings/', type, '/USGS.', siteNumber, '.', type, '.rdb')

  data <- dataRetrieval::importRDB1(url, asDateTime = FALSE, convertType = convertType)

  if ("current_rating_nu" %in% names(data)) {
    intColumns <- intColumns[!("current_rating_nu" %in% names(data)[intColumns])]
    data$current_rating_nu <- gsub(" ", "", data$current_rating_nu)
  }

  if (nrow(data) > 0) {
    if (type == "base") {
      Rat <- grep("//RATING ", comment(data), value = TRUE, fixed = TRUE)
      Rat <- sub("# //RATING ", "", Rat)
      Rat <- scan(text = Rat, sep = " ", what = "")
      attr(data, "RATING") <- Rat
    }

    siteInfo <- suppressMessages(dataRetrieval::readNWISsite(siteNumbers = siteNumber))

    attr(data, "siteInfo") <- siteInfo
    attr(data, "variableInfo") <- NULL
    attr(data, "statisticInfo") <- NULL
  }

  return(data)
}




makeGageDF <- function(gage, model_gages, huc4){
  out <- gage %>%
    dplyr::bind_rows() %>%
    dplyr::filter(site_no %in% model_gages$GageID) %>%
    dplyr::mutate(huc4 = huc4) %>%
    dplyr::relocate(huc4)

  return(out)
}