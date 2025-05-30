## Craig Brinkerhoff
## Summer 2025
## Functions for ML model




summariseInn <- function(gageInundation, gageTau){
    gageInundation <- gageInundation %>%
        dplyr::mutate(Af_q0_2_m2 = A_q0_2_m2 - (Wb_m*LengthKM*1000),
                    Af_q0_5_m2 = A_q0_5_m2 - (Wb_m*LengthKM*1000),
                    Af_q1_m2 = A_q1_m2 - (Wb_m*LengthKM*1000),
                    Af_q2_m2 = A_q2_m2 - (Wb_m*LengthKM*1000),
                    Af_q4_m2 = A_q4_m2 - (Wb_m*LengthKM*1000),
                    Af_q10_m2 = A_q10_m2 - (Wb_m*LengthKM*1000),
                    Af_q20_m2 = A_q20_m2 - (Wb_m*LengthKM*1000),
                    Af_q50_m2 = A_q50_m2 - (Wb_m*LengthKM*1000),
                    Af_q80_m2 = A_q80_m2 - (Wb_m*LengthKM*1000),
                    Af_q90_m2 = A_q90_m2 - (Wb_m*LengthKM*1000),
                    Vf_q0_2_m3 = V_q0_2_m3 - (Wb_m*LengthKM*1000*Htf_q0_2_m),
                    Vf_q0_5_m3 = V_q0_5_m3 - (Wb_m*LengthKM*1000*Htf_q0_5_m),
                    Vf_q1_m3 = V_q1_m3 - (Wb_m*LengthKM*1000*Htf_q1_m),
                    Vf_q2_m3 = V_q2_m3 - (Wb_m*LengthKM*1000*Htf_q2_m),
                    Vf_q4_m3 = V_q4_m3 - (Wb_m*LengthKM*1000*Htf_q4_m),
                    Vf_q10_m3 = V_q10_m3 - (Wb_m*LengthKM*1000*Htf_q10_m),
                    Vf_q20_m3 = V_q20_m3 - (Wb_m*LengthKM*1000*Htf_q20_m),
                    Vf_q50_m3 = V_q50_m3 - (Wb_m*LengthKM*1000*Htf_q50_m),
                    Vf_q80_m3 = V_q80_m3 - (Wb_m*LengthKM*1000*Htf_q80_m),
                    Vf_q90_m3 = V_q90_m3 - (Wb_m*LengthKM*1000*Htf_q90_m)) %>%
        dplyr::mutate(Af_q0_2_m2 = ifelse(Af_q0_2_m2 < 0, 0, Af_q0_2_m2),
                    Af_q0_5_m2 = ifelse(Af_q0_5_m2 < 0, 0, Af_q0_5_m2),
                    Af_q1_m2 = ifelse(Af_q1_m2 < 0, 0, Af_q1_m2),
                    Af_q2_m2 = ifelse(Af_q2_m2 < 0, 0, Af_q2_m2),
                    Af_q4_m2 = ifelse(Af_q4_m2 < 0, 0, Af_q4_m2),
                    Af_q10_m2 = ifelse(Af_q10_m2 < 0, 0, Af_q10_m2),
                    Af_q20_m2 = ifelse(Af_q20_m2 < 0, 0, Af_q20_m2),
                    Af_q50_m2 = ifelse(Af_q50_m2 < 0, 0, Af_q50_m2),
                    Af_q80_m2 = ifelse(Af_q80_m2 < 0, 0, Af_q80_m2),
                    Af_q90_m2 = ifelse(Af_q90_m2 < 0, 0, Af_q90_m2),
                    Vf_q0_2_m3 = ifelse(Af_q0_2_m2 < 0, 0, Af_q0_2_m2),
                    Vf_q0_5_m3 = ifelse(Af_q0_5_m2 < 0, 0, Af_q0_5_m2),
                    Vf_q1_m3 = ifelse(Af_q1_m2 < 0, 0, Af_q1_m2),
                    Vf_q2_m3 = ifelse(Af_q2_m2 < 0, 0, Af_q2_m2),
                    Vf_q4_m3 = ifelse(Af_q4_m2 < 0, 0, Af_q4_m2),
                    Vf_q10_m3 = ifelse(Af_q10_m2 < 0, 0, Af_q10_m2),
                    Vf_q20_m3 = ifelse(Af_q20_m2 < 0, 0, Af_q20_m2),
                    Vf_q50_m3 = ifelse(Af_q50_m2 < 0, 0, Af_q50_m2),
                    Vf_q80_m3 = ifelse(Af_q80_m2 < 0, 0, Af_q80_m2),
                    Vf_q90_m3 = ifelse(Af_q90_m2 < 0, 0, Af_q90_m2))
    
    
    gageTau <- gageTau %>%
        dplyr::bind_rows() %>%
        dplyr::distinct() %>% #occasionally, i.e. in huc 0306, there are duplicate gages from dataRetrieval that mess things up (and are unncessary). So remove those rows
        tidyr::pivot_wider(id_cols=site_no, names_from=Qprob, values_from=c(avg_event_exc_dys, avg_ann_exc_dys, avg_ann_n_events))

    out <- gageInundation %>%
        dplyr::left_join(gageTau, by=c('GageID'='site_no')) %>%
        sf::st_drop_geometry()
    
    return(out)
}




calcGageTau <- function(gageRecord, preppedForInn_probs, maxDiff){
    if(nrow(gageRecord)==0){
        return(data.frame())
    }

    probs <- preppedForInn_probs %>%
        dplyr::filter(site_no == gageRecord[1,]$site_no)
    
    if(nrow(probs)==0){
        return(data.frame())
    }

    years <- max(lubridate::year(gageRecord$dateTime))-min(lubridate::year(gageRecord$dateTime))

    #Q0.2 flood
    events_02 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q0_2_cms)/probs$Q_q0_2_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q0_2_dys')
    
    #Q0.5 flood
    events_05 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q0_5_cms)/probs$Q_q0_5_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q0_5_dys')
    
    #Q1 flood
    events_1 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q1_cms)/probs$Q_q1_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q1_dys')
    
    #Q2 flood
    events_2 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q2_cms)/probs$Q_q2_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q2_dys')

    #Q4 flood
    events_4 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q4_cms)/probs$Q_q4_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q4_dys')

    #Q10 flood
    events_10 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q10_cms)/probs$Q_q10_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q10_dys')

    #Q20 flood
    events_20 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q20_cms)/probs$Q_q20_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q20_dys')

    #Q50 flood
    events_50 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q50_cms)/probs$Q_q50_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q50_dys')

    #Q80 flood
    events_80 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q80_cms)/probs$Q_q80_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q80_dys')

    #Q90 flood
    events_90 <- gageRecord %>%
        dplyr::mutate(diff = abs(Q_cms - probs$Q_q90_cms)/probs$Q_q90_cms) %>%
        dplyr::filter(diff < maxDiff) %>%
        dplyr::mutate(Qprob = 'tau_q90_dys')

    events <- rbind(events_02, events_05, events_1, events_2, events_4, events_10, events_20, events_50, events_80, events_90)

    if(nrow(events)==0){
        return(data.frame())
    }

    out <- data.frame()
    
    for(m in 1:nrow(events)){
        event_row <- which(gageRecord$dateTime == events[m,]$dateTime)

        flag <- 0
        i <- 1
        k <- 1
        event_tally <- events[m,] %>%
                dplyr::select(c('site_no', 'dateTime', 'Q_cms', 'stage_m'))

        #crawl forward and backward in time to find event length
        while(flag == 0){
            forward <- gageRecord[event_row+i,] %>% #note, going beyond the end of df will just result in NA df
                dplyr::select(c('site_no', 'dateTime', 'Q_cms', 'stage_m'))
            
            if((event_row-k) < 1){ #can't index 0, so manually create NA df
              backward <- data.frame('site_no'=NA,'dateTime'=NA, 'Q_cms'=NA, 'stage_m'=NA)
            }
            else{
              backward <- gageRecord[event_row-k,]%>%
                dplyr::select(c('site_no', 'dateTime', 'Q_cms', 'stage_m'))
            }

            if(sum(is.na(forward)) == 0 & forward$Q_cms > probs$Qb_cms & forward$Q_cms <= events[m,]$Q_cms){
              event_tally <- rbind(event_tally, forward)
              i <- i + 1
            }
            if(sum(is.na(backward)) == 0 & backward$Q_cms > probs$Qb_cms& backward$Q_cms <= events[m,]$Q_cms){
                event_tally <- rbind(event_tally, backward)
                k <- k + 1
            }
          else{
                flag <- 1
          }
        }

        temp <- data.frame('site_no'=events[m,]$site_no,
                        'Qprob'=events[m,]$Qprob,
                        'tau_dys'=nrow(event_tally),
                        'Q_flood_cms'=sum(event_tally$Q_cms) - sum(probs$Ub_ms*probs$Wb_m*(event_tally$stage_m-probs$Hb_m)),
                        'years'=years)
        out <- rbind(out, temp)
    }

    #average across events
    out <- out %>%
        dplyr::group_by(Qprob) %>%
        dplyr::summarise(site_no = first(site_no), #dummy to pass through
                        avg_event_Q_flood_cms = mean(Q_flood_cms),
                        avg_event_exc_dys=mean(tau_dys),
                        avg_ann_exc_dys = sum(tau_dys)/mean(years), #mean of constant to pass through
                        avg_ann_n_events=n()/mean(years)) #mean of constant to pass through
    
    return(out)
}









getBasinGageshuc4 <- function(huc4id, gageRecordStart, gageRecordEnd, benchmark_index){
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



prepForInn <- function(huc4){
  huc2 <- substr(huc4, 1, 2)
  
  #setup basin shapefiles
  dem <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/elev_cm.tif')) #[cm]
  d8 <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/fdr.tif')) #[cm]

  prelist <- list('dem'=dem,
       'd8'=d8)
  outlist <- lapply(prelist, terra::wrap) #terra holds C++ pointers in memory so making lists of objects (for distributed computing) won't work- we need to 'wrap' the objects into a packed state that can be sent over a serialized connection)
  
  
  return('terra'=outlist)
}





calcFlowProbs <- function(huc4, gageRecord, gage, BHGmodel, depAHG, minAHGr2, gageRecordStart, gageRecordEnd){

  #setup
  depAHG <- dplyr::bind_rows(depAHG)

  gageRecord <- dplyr::bind_rows(gageRecord) %>%
    dplyr::select(c('site_no', 'fema100yrflood_cms', 'Q_cms', 'exceed_prob')) #date_hr

  #filter for flows above bankfull, to calculate exceedance probabilities for flood events (aside from the max-annual AEP used to validate against FEMA)
  gage <- sf::st_drop_geometry(gage) %>% 
    dplyr::bind_rows() %>% 
    dplyr::select(c('site_no', 'DA_skm', 'physio_region')) %>% 
    sf::st_drop_geometry() %>%
    dplyr::filter(!(is.na(physio_region)))
  
  gage$a_Qb <- as.numeric(sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Qb)}))
  gage$b_Qb <- as.numeric(sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Qb)}))
  gage$Qb_cms_log10 <- gage$a_Qb + gage$b_Qb*log10(gage$DA_skm)

  gage$a_Wb <- as.numeric(sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Wb)}))
  gage$b_Wb <- as.numeric(sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Wb)}))
  gage$Wb_m_log10 <- gage$a_Wb + gage$b_Wb*log10(gage$DA_skm)

  gage$a_Ab <- as.numeric(sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$a_Ab)}))
  gage$b_Ab <- as.numeric(sapply(gage$physio_region, function(x){return(BHGmodel[BHGmodel$division == x,]$b_Ab)}))
  gage$Ab_m2_log10 <- gage$a_Ab + gage$b_Ab*log10(gage$DA_skm)

  #adjust for back-transformed bias using the mean residual
  BHGmodel <- BHGmodel %>%
    dplyr::select('division', 'mean_residual_Qb', 'mean_residual_Wb', 'mean_residual_Ab')
  
  gage <- gage %>%
    dplyr::left_join(BHGmodel, by=c('physio_region'='division')) %>%
    dplyr::mutate(Qb_cms = 10^(Qb_cms_log10) * mean_residual_Qb,
                  Wb_m = 10^(Wb_m_log10) * mean_residual_Wb,
                  Ab_m2 = 10^(Ab_m2_log10) * mean_residual_Ab) %>% #using the mean log residual, in transformed to ntural space, a la eq 9.24 from https://pubs.usgs.gov/tm/04/a03/tm4a3.pdf (also see https://nrtwq.usgs.gov/co/methods/)
    dplyr::mutate(Ub_ms = Qb_cms / Ab_m2)
  
  #bring everything together
  upscaling_df <- gageRecord %>%
    dplyr::left_join(gage, by='site_no') %>%
    dplyr::left_join(depAHG, by='site_no') %>%
    dplyr::filter(!(is.na(DA_skm))) %>%
    dplyr::filter(lm_r2 > minAHGr2)

  #filter for flows above bankfull and calculate exceed prob from that
  num <- upscaling_df %>%
    dplyr::filter(Q_cms > Qb_cms) %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(n=n())
  
  exceed_prob_df <- upscaling_df %>%
    dplyr::left_join(num, by='site_no') %>%
    dplyr::filter(Q_cms > Qb_cms) %>%
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
  upscaling_df$Htf_m_notcorrected <- 10^(upscaling_df$lm_depth_a + (upscaling_df$lm_depth_b * log10(upscaling_df$Q_cms)))
  upscaling_df$Htf_m <- upscaling_df$Htf_m_notcorrected * upscaling_df$lm_depth_bias_correct

  upscaling_df$femaHtf_m_notcorrected <- 10^(upscaling_df$lm_depth_a + (upscaling_df$lm_depth_b * log10(upscaling_df$fema100yrflood_cms)))
  upscaling_df$femaHtf_m <- upscaling_df$femaHtf_m_notcorrected * upscaling_df$lm_depth_bias_correct

  #Bankful depth
  qBankful_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(Qb_cms - Q_cms)) %>%
    dplyr::slice_min(diff, with_ties=FALSE) %>%
    dplyr::mutate(Hb_m = Htf_m) %>%
    dplyr::select(c('site_no','Hb_m'))

  #100yr FEMA AEP flood
  qFEMA_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::summarise(Htf_qFEMA_m = dplyr::first(femaHtf_m)) %>% #it's a constant so just use the first
    dplyr::select(c('site_no', 'Htf_qFEMA_m'))

  #Q0.2 flood
  q0_2_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.002)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q0_2_m = Htf_m,
                Q_q0_2_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q0_2_m', 'Q_q0_2_cms'))

  #Q0.5 flood
  q0_5_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.005)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q0_5_m = Htf_m,
                Q_q0_5_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q0_5_m', 'Q_q0_5_cms'))

  #Q1 flood
  q1_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.01)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q1_m = Htf_m,
                Q_q1_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q1_m', 'Q_q1_cms'))

  #Q2 flood
  q2_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.02)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q2_m = Htf_m,
                Q_q2_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q2_m', 'Q_q2_cms'))
  
  #Q4 flood
  q4_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.04)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q4_m = Htf_m,
                Q_q4_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q4_m', 'Q_q4_cms'))

  #Q10 flood
  q10_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.10)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q10_m = Htf_m,
                Q_q10_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q10_m', 'Q_q10_cms'))
  
  #Q20 flood
  q20_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.20)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q20_m = Htf_m,
                Q_q20_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q20_m', 'Q_q20_cms'))
  
  #Q50 flood
  q50_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.50)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q50_m = Htf_m,
                Q_q50_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q50_m', 'Q_q50_cms'))
  
  #Q80 flood
  q80_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.80)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q80_m = Htf_m,
                Q_q80_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q80_m', 'Q_q80_cms'))

  #Q90 flood
  q90_forModel <- upscaling_df %>%
    dplyr::group_by(site_no) %>%
    dplyr::mutate(diff = abs(exceed_prob_bankfull - 0.90)) %>% 
    dplyr::slice_min(diff, with_ties=FALSE) %>% #some exceed prob ties have stage that are ~0.01m apart, just take the first
    dplyr::mutate(Htf_q90_m = Htf_m,
                Q_q90_cms = Q_cms) %>%
    dplyr::select(c('site_no','Htf_q90_m', 'Q_q90_cms'))
  
  out <- gage %>%
    dplyr::left_join(qFEMA_forModel, by='site_no') %>%
    dplyr::left_join(q0_2_forModel, by='site_no') %>%
    dplyr::left_join(q0_5_forModel, by='site_no') %>%
    dplyr::left_join(q1_forModel, by='site_no') %>%
    dplyr::left_join(q2_forModel, by='site_no') %>%
    dplyr::left_join(q4_forModel, by='site_no') %>%
    dplyr::left_join(q10_forModel, by='site_no') %>%
    dplyr::left_join(q20_forModel, by='site_no') %>%
    dplyr::left_join(q50_forModel, by='site_no') %>%
    dplyr::left_join(q80_forModel, by='site_no') %>%
    dplyr::left_join(q90_forModel, by='site_no') %>%
    dplyr::left_join(qBankful_forModel, by='site_no') %>%
    tidyr::drop_na() %>%
    dplyr::distinct() #in basin 0306 there were duplicate enteries for one gage, so just confirm this isn't happening elsewhere

  return(out)
}




prepForInn_gages <- function(huc4, depthAHG, BHGmodel, preppedForInn_probs) {
  huc2 <- substr(huc4, 1, 2)

  gage_list <- dplyr::bind_rows(depthAHG)
  
  network_VAA <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusFlowlineVAA', quiet=TRUE)
  network_gages <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'), layer='NHDPlusEROMQAMA', quiet=TRUE)
  network <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/NHDPLUS_H_',huc4,'_HU4_GDB.gdb'),layer = 'NHDFlowline',quiet = TRUE) %>%
    sf::st_zm() %>%
    dplyr::left_join(network_VAA, by='NHDPlusID') %>%
    dplyr::left_join(network_gages, by='NHDPlusID') %>%
    dplyr::select(c('NHDPlusID', 'WBArea_Permanent_Identifier', 'GageID','StreamCalc', 'AreaSqKm', 'TotDASqKm', 'LengthKM', 'Slope', 'Shape'))

  #remove impossible flowlines with no drainage area or divergent starting reaches (streamalc == 0; see pg. 45 at https://pubs.usgs.gov/of/2019/1096/ofr20191096.pdf)
  network <- network %>%
    dplyr::filter(AreaSqKm > 0 & TotDASqKm > 0 & StreamCalc > 0) #%>%
    #dplyr::filter(GageID %in% gage_list$site_no)
  
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


  network <- network %>%
    dplyr::left_join(preppedForInn_probs, by=c('GageID'='site_no')) %>%
    tidyr::drop_na() %>%
    dplyr::select(!c('WBArea_Permanent_Identifier', 'a_Wb', 'b_Wb', 'mean_residual_Wb', 'a_Qb', 'b_Qb', 'mean_residual_Qb'))
  
  return(network)
}





calcGageInundation <- function(huc4, network, minWidth){
  #setup basin shapefiles
  # dem <- terra::unwrap(basinData$dem) #[cm] need to unwrap the packed terra package (for distributed computing)
  # d8 <- terra::unwrap(basinData$d8) #[cm] need to unwrap the packed terra package (for distributed computing)
  huc2 <- substr(huc4, 1, 2)
  dem <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/elev_cm.tif')) #[cm]
  d8 <- terra::rast(paste0('data/path_to_data/CONUS_ephemeral_data/HUC2_', huc2, '/NHDPLUS_H_',huc4,'_HU4_GDB/fdr.tif')) #[cm]

  #run embarresingly parallel inundation model
  inundationWrapper <- function(reach) {
    #lower bound on dem resolution (also, these small rivers don't reallllly have floodplains, right?)
    if(reach$Wb_m <= minWidth){
      reach[,c('A_qFEMA_m2', 'A_q0_2_m2', 'A_q0_5_m2', 'A_q1_m2', 'A_q2_m2', 'A_q4_m2', 'A_q10_m2', 'A_q20_m2', 'A_q50_m2', 'A_q80_m2', 'A_q90_m2')] <-NA
      reach[,c('V_qFEMA_m3', 'V_q0_2_m3', 'V_q0_5_m3', 'V_q1_m3', 'V_q2_m3', 'V_q4_m3', 'V_q10_m3', 'V_q20_m3', 'V_q50_m3', 'V_q80_m3', 'V_q90_m3')] <- NA
      return(reach)
    }
    
    #run inundation model
    inundationDataPackage <- prepInundationData(huc4, reach$NHDPlusID, reach$Wb_m, dem, d8)
    if(inundationDataPackage$flag == 'no pour points'){
      reach[,c('A_qFEMA_m2', 'A_q0_2_m2', 'A_q0_5_m2', 'A_q1_m2', 'A_q2_m2', 'A_q4_m2', 'A_q10_m2', 'A_q20_m2', 'A_q50_m2', 'A_q80_m2', 'A_q90_m2')] <-NA
      reach[,c('V_qFEMA_m3', 'V_q0_2_m3', 'V_q0_5_m3', 'V_q1_m3', 'V_q2_m3', 'V_q4_m3', 'V_q10_m3', 'V_q20_m3', 'V_q50_m3', 'V_q80_m3', 'V_q90_m3')] <- NA
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

    reach[,c('A_qFEMA_m2', 'A_q0_2_m2', 'A_q0_5_m2', 'A_q1_m2', 'A_q2_m2', 'A_q4_m2', 'A_q10_m2', 'A_q20_m2', 'A_q50_m2', 'A_q80_m2', 'A_q90_m2')] <- lapply(reach[,c('Htf_qFEMA_m', 'Htf_q0_2_m', 'Htf_q0_5_m', 'Htf_q1_m', 'Htf_q2_m', 'Htf_q4_m', 'Htf_q10_m', 'Htf_q20_m', 'Htf_q50_m', 'Htf_q80_m', 'Htf_q90_m')], areaWrapper)
    reach[,c('V_qFEMA_m3', 'V_q0_2_m3', 'V_q0_5_m3', 'V_q1_m3', 'V_q2_m3', 'V_q4_m3', 'V_q10_m3', 'V_q20_m3', 'V_q50_m3', 'V_q80_m3', 'V_q90_m3')] <- lapply(reach[,c('Htf_qFEMA_m', 'Htf_q0_2_m', 'Htf_q0_5_m', 'Htf_q1_m', 'Htf_q2_m', 'Htf_q4_m', 'Htf_q10_m', 'Htf_q20_m', 'Htf_q50_m', 'Htf_q80_m', 'Htf_q90_m')], volWrapper)

    return(reach)
  }

  library(plyr)
#   library(doParallel)
#   registerDoParallel(cores=detectCores()-2)

  network2 <- sf::st_drop_geometry(network)
  network_list <- setNames(split(network2, seq(nrow(network2))), rownames(network2))
  result <- llply(network_list, inundationWrapper, .parallel=FALSE)
  network_fin <- dplyr::bind_rows(result)

  #turn back into shapefile and prep for output
  network <- network %>%
    dplyr::select(c('NHDPlusID', 'Shape')) %>%
    dplyr::left_join(network_fin, by='NHDPlusID') %>%
    dplyr::relocate(Shape, .after = tidyselect::last_col()) %>%
    tidyr::drop_na() #drop reaches too narrow for model
  
  return(network)
}