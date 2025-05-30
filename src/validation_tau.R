## Craig Brinkerhoff
## Summer 2025

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








calc_Qexc <- function(mode, gageRecord, gagePrepped, depAHG, maxDiff){
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


