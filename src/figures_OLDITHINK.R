makeValFEMA <- function(val_FEMA_combined, gage_combined, BHGdata){
    library(ggplot2)
    theme_set(theme_classic())

    #make inundated area val df
    val_FEMA_area <- val_FEMA_combined %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(A_femaAEP_km2 > 0 & A_qFEMA_km2 > 0)    

    #make floodplain only val df
    val_FEMA_area_fp <- val_FEMA_area %>%
        dplyr::left_join(BHGdata, by='GageID') %>%
        dplyr::left_join(gage_combined, by=c('GageID'='site_no')) %>%
        dplyr::filter(!(is.na(Wb_m_obs))) %>%
        dplyr::mutate(physio_region = stringr::str_to_title(physio_region))

    #refit Wb equations while witholding the sites we can validate (WB obs)
    Wb_model_val <- BHGdata %>%
        dplyr::filter(!(GageID %in% val_FEMA_area_fp$GageID)) %>%
        dplyr::group_by(DIVISION) %>%
        dplyr::do(model_Wb = lm(log10(Wb_m_obs)~log10(DA_gage_skm), data=.)) %>%
        dplyr::mutate(a_Wb = 10^(model_Wb$coef[1]), #model intercept
                    b_Wb = model_Wb$coef[2], #model exponent
                    r2_Wb = summary(model_Wb)$r.squared) #model performance
    
    Wb_model_val[Wb_model_val$DIVISION == "Intermontane Plateau",]$DIVISION <- "Intermontane Plateaus" #make sure names line up

    val_FEMA_area_fp$a_valmodel_Wb <- sapply(val_FEMA_area_fp$physio_region, function(x){return(Wb_model_val[Wb_model_val$DIVISION == x,]$a_Wb)})
    val_FEMA_area_fp$b_valmodel_Wb <- sapply(val_FEMA_area_fp$physio_region, function(x){return(Wb_model_val[Wb_model_val$DIVISION == x,]$b_Wb)})
    val_FEMA_area_fp$Wb_valmodel_m <- val_FEMA_area_fp$a_valmodel_Wb * (val_FEMA_area_fp$DA_gage_skm)^val_FEMA_area_fp$b_valmodel_Wb

    val_FEMA_area_fp <- val_FEMA_area_fp %>%
        dplyr::mutate(channel_area_obs_km2 = (Wb_m_obs/1000) * LengthKM,
                    channel_area_valmodel_km2 = (Wb_valmodel_m/1000) * LengthKM) %>%
        dplyr::mutate(Af_femaAEP_km2 = A_femaAEP_km2 - channel_area_obs_km2,
                    Af_qFEMA_valmodel_km2 = A_qFEMA_km2 - channel_area_valmodel_km2) %>%
        dplyr::filter(Af_femaAEP_km2 > 0 & Af_qFEMA_valmodel_km2 > 0)

    lm_inn_a <- lm(log(A_qFEMA_km2)~log(A_femaAEP_km2), data=val_FEMA_area)
    lm_flood <- lm(log(Af_qFEMA_valmodel_km2)~log(Af_femaAEP_km2), data=val_FEMA_area_fp)

    scatter_area <- ggplot(val_FEMA_area, aes(x=A_femaAEP_km2*1e6, y=A_qFEMA_km2*1e6))+
        geom_hex(bins=20)+
        #geom_point(size=6, alpha=0.3, color='#588157') +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        scale_fill_distiller(palette='YlOrBr', name='Count', direction=1)+
       # geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks",limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab(bquote(bold("FEMA 100yr flood ["*m^2*"]"))) +
        ylab(bquote(bold("Modeled 100yr flood ["*m^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
           # axis.text.x=element_blank(),
           # axis.ticks.x=element_blank(),
            axis.line = element_blank(),
            legend.position=c(0.9,0.175),
            plot.title = element_text(size=22,hjust = 0.5),
            panel.border = element_rect(colour = "black", fill = NA)) +
        labs(title='Total area')+
        annotate("text", x = 1e2, y = 10^6.5, label = bquote(.(nrow(val_FEMA_area))*' reaches'), size=6)+
        annotate("text", x = 1e2, y = 10^6, label = bquote(R^2*'= '*.(round(summary(lm_inn_a)$r.squared,2))), size=6)

    scatter_area_fp <- ggplot(val_FEMA_area_fp, aes(x=Af_femaAEP_km2*1e6, y=Af_qFEMA_valmodel_km2*1e6))+
        geom_point(size=6, color='#588157') +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab(bquote(bold("FEMA 100yr inundated area ["*m^2*"]"))) +
        ylab(bquote(bold("Modeled 100yr inundated area ["*m^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none',
            axis.line = element_blank(),
            plot.title = element_text(size=22, hjust = 0.5),
            panel.border = element_rect(colour = "black", fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        labs(title='Floodplain area')+
        annotate("text", x = 1e2, y = 10^6.5, label = bquote(.(nrow(val_FEMA_area_fp))*' reaches'), size=6)+
        annotate("text", x = 1e2, y = 10^6, label = bquote(R^2*'= '*.(round(summary(lm_flood)$r.squared,2))), size=6)
    
    layout <- "
        A
        B
    "

    comboPlot <- patchwork::wrap_plots(A=scatter_area, B=scatter_area_fp, design=layout)

    ggsave('cache/validationFEMA.png', comboPlot, width=6, height=12)
    ggsave('cache/validationFEMA_jobtalk.png', scatter_area, width=7, height=7)

    return(val_FEMA_combined)
}




makeValUSGSarea <- function(val_USGS_combined){
    library(ggplot2)
    theme_set(theme_classic())

    #make inundated area val df
    val_USGS_area <- val_USGS_combined %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(A_usgs_km2 > 0 & A_model_km2 > 0) %>%
        dplyr::mutate(plot_var = factor(key_exdprob, c('A_q0_2', 'A_q0_5', 'A_q1', 'A_q2', 'A_q4', 'A_q10', 'A_q20', 'A_q50'), labels=c('0.2% flood', '0.5% flood', '1% flood', '2% flood', '4% flood', '10% flood', '20% flood', '50% flood')))

    lm_inn_a <- lm(log(A_model_km2)~log(A_usgs_km2), data=val_USGS_area)

    #get r2 and num reaches for each characteristic flood
    df_r2 <- val_USGS_area %>%
        dplyr::group_by(plot_var) %>%
        dplyr::group_modify(~ broom::glance(lm(log10(A_model_km2) ~ log10(A_usgs_km2), data = .x))) %>% 
        dplyr::select(c('plot_var', 'r.squared'))
    
    df_num <- val_USGS_area %>%
        dplyr::group_by(plot_var) %>%
        dplyr::summarise(n=n())
    
    lookup <- df_r2 %>%
        dplyr::left_join(df_num, by='plot_var') %>%
        dplyr::mutate(r2fin = sprintf("italic(R^2) == %.2f", r.squared),
                    nfin = paste0(n, ' reaches'))

    #plot
    scatter_area <- ggplot(val_USGS_area, aes(x=A_usgs_km2*1e6, y=A_model_km2*1e6))+
        geom_point(size=5, color='#9a8c98') +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab(bquote(bold("USGS-hydrodynamic inundated area ["*m^2*"]"))) +
        ylab(bquote(bold("Modeled inundated area ["*m^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            strip.text = element_text(size=22),
            axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA)) +
        facet_wrap(vars(plot_var), nrow=2) +
        geom_text(x=2, y=6, aes(label=r2fin), data=lookup, parse=TRUE, size=5)+
        geom_text(x=2, y=6.5, aes(label=nfin), data=lookup, size=5)

    ggsave('cache/validationUSGS_area.png', scatter_area, width=16, height=10)#width=6, height=12)

    return(val_USGS_combined)
}






makeValUSGSvol <- function(val_USGSvols_combined){    
    library(ggplot2)
    theme_set(theme_classic())

    # bias_correct <- val_USGSvols_combined %>%
    #     dplyr::filter(V_model_km3 > 0 & V_usgs_km3 > 0) #%>%
    #     dplyr::group_by(key_exdprob) %>%
    #     dplyr::summarise(bias = Metrics::bias(log(V_usgs_km3), log(V_model_km3)),
    #                     n = n())

    #make inundated area val df
    val_USGS_vol <- val_USGSvols_combined %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(V_usgs_km3 > 0 & V_model_km3 > 0) %>%
        # dplyr::left_join(bias_correct, by='key_exdprob') %>%
        # dplyr::mutate(V_model_km3_biascorrect = exp(log(V_model_km3) + bias)) %>%
        # tidyr::gather(key=key_bias, value=value, c('V_model_km3', 'V_model_km3_biascorrect')) %>%
        dplyr::mutate(plot_var = factor(key_exdprob, c('V_q0_2', 'V_q0_5', 'V_q1', 'V_q2', 'V_q4', 'V_q10', 'V_q20', 'V_q50'), labels=c('0.2% flood', '0.5% flood', '1% flood', '2% flood', '4% flood', '10% flood', '20% flood', '50% flood')))

    #get r2 and num reaches for each characteristic flood
    df_r2 <- val_USGS_vol %>%
        # dplyr::filter(key_bias == 'V_model_km3_biascorrect') %>%
        dplyr::group_by(plot_var) %>%
        dplyr::group_modify(~ broom::glance(lm(log10(V_model_km3) ~ log10(V_usgs_km3), data = .x))) %>% 
        dplyr::select(c('plot_var', 'r.squared'))
    
    df_num <- val_USGS_vol %>%
        # dplyr::filter(key_bias == 'V_model_km3_biascorrect') %>%
        dplyr::group_by(plot_var) %>%
        dplyr::summarise(n=n())
    
    lookup <- df_r2 %>%
        dplyr::left_join(df_num, by='plot_var') %>%
        dplyr::mutate(r2fin = sprintf("italic(R^2) == %.2f", r.squared),
                    nfin = paste0(n, ' reaches'))

    #plot
    scatter_vol <- ggplot(val_USGS_vol, aes(x=V_usgs_km3*1e9, y=V_model_km3*1e9))+#, color=key_bias))+
        geom_point(size=5, color='#a39171') +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        # scale_color_manual(name='', labels=c('Model', 'Bias Corrected Model'), values=c('#81b29a', '#3d405b'))+
        scale_x_log10(guide = "axis_logticks", limits=c(1e2, 1e7), breaks=c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e2, 1e7), breaks=c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab(bquote(bold("USGS-hydrodynamic inundated volume ["*m^3*"]"))) +
        ylab(bquote(bold("Modeled inundated volume ["*m^3*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            strip.text = element_text(size=22),
            axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            legend.position='bottom',
            legend.text=element_text(size=15)) +
        facet_wrap(vars(plot_var), nrow=2)+
        geom_text(x=3, y=6, aes(label=r2fin), data=lookup, parse=TRUE, size=5)+
        geom_text(x=3, y=6.5, aes(label=nfin), data=lookup, size=5)

    ggsave('cache/validationUSGS_vol.png', scatter_vol, width=16, height=10)

    return(val_USGSvols_combined)
}







upscalingFig <- function(df, huc2){
    library(ggplot2)

    theme_set(theme_classic())

    forLabels <- df$models %>%
        dplyr::filter(stringr::str_detect(name, 'Htf_m')) %>%
        dplyr::mutate(r2fin = sprintf("italic(R^2) == %.2f", rsq),
                    nfin = paste0(n_gages, ' gages'))
  
    forPlot <- df$df %>%
        dplyr::filter(stringr::str_detect(name, 'Htf_m'))

    plot <- ggplot(forPlot, aes(x=DA_skm, y=value)) +
        geom_point(size=5, color='#153131')+
        geom_smooth(method='lm', se=F, color='black',linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks")+
        scale_y_log10(guide = "axis_logticks") +
        ylab(bquote(bold("Water level ["*m*"]"))) +
        xlab(bquote(bold("Drainage Area ["*km^2*"]"))) +
        theme(strip.text.x = element_text(size = 20),
            axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            strip.text = element_text(size=22),
            axis.line = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA)) +
        facet_wrap(~factor(name, levels=c('qFEMAflood_Htf_m', 'q0_2flood_Htf_m', 'q0_5flood_Htf_m', 'q1flood_Htf_m', 'q2flood_Htf_m', 'q4flood_Htf_m', 'q10flood_Htf_m', 'q20flood_Htf_m', 'q50flood_Htf_m', 'q80flood_Htf_m', 'q90flood_Htf_m', 'q96flood_Htf_m', 'q98flood_Htf_m', 'q99flood_Htf_m', 'q99_5flood_Htf_m', 'q99_8flood_Htf_m'),
                                labels=c('100yr AEP flood', '0.2% flood', '0.5% flood', '1% flood', '2% flood', '4% flood', '10% flood', '20% flood', '50% flood', '80% flood', '90% flood', '96% flood', '98% flood', '99% flood', '99.5% flood', '99.8% flood'))) +
        geom_text(x=0.05, y=0.95, data = forLabels, aes(label = r2fin), size=6, parse=TRUE) +
        geom_text(x=0.05, y=0.8, data = forLabels, aes(label = nfin), size=6)
  
  ggsave(paste0('cache/upscalingModel_', huc2, '.png'), plot, width=14, height=14)

  return(paste0('cache/upscalingModel_', huc2, '.png'))
}





# residenceTimeFig <- function(df){
#     library(ggplot2)
#     theme_set(theme_classic())

#     forPlot_f <- df %>%
#         sf::st_drop_geometry() %>%
#         tidyr::gather(key=key, value=value, c('HRTf_q0_5_s','HRTf_q1_s','HRTf_q10_s','HRTf_q20_s','HRTf_q50_s','HRTf_q80_s','HRTf_q90_s','HRTf_q99_s','HRTf_q99_5_s')) %>%
#         dplyr::filter(!(is.na(value))) %>%
#         dplyr::mutate(key = forcats::fct_rev(factor(key, levels=c('HRTf_q0_5_s','HRTf_q1_s','HRTf_q10_s','HRTf_q20_s','HRTf_q50_s','HRTf_q80_s','HRTf_q90_s','HRTf_q99_s','HRTf_q99_5_s'),
#                                 labels=c('0.5%', '1%', '10%', '20%', '50%', '80%', '90%','99%', '99.5'))))
    
#     forPlot_c <- df %>%
#         sf::st_drop_geometry() %>%
#         tidyr::gather(key=key, value=value, c('HRTc_q0_5_s','HRTc_q1_s','HRTc_q10_s','HRTc_q20_s','HRTc_q50_s','HRTc_q80_s','HRTc_q90_s','HRTc_q99_s','HRTc_q99_5_s')) %>%
#         dplyr::filter(!(is.na(value))) %>%
#         dplyr::mutate(key = forcats::fct_rev(factor(key, levels=c('HRTc_q0_5_s','HRTc_q1_s','HRTc_q10_s','HRTc_q20_s','HRTc_q50_s','HRTc_q80_s','HRTc_q90_s','HRTc_q99_s','HRTc_q99_5_s'),
#                                 labels=c('0.5%', '1%', '10%', '20%', '50%', '80%', '90%','99%', '99.5'))))

#     forPlot_tf <- df %>%
#         sf::st_drop_geometry() %>%
#         tidyr::gather(key=key, value=value, c('HRTtf_q0_5_s','HRTtf_q1_s','HRTtf_q10_s','HRTtf_q20_s','HRTtf_q50_s','HRTtf_q80_s','HRTtf_q90_s','HRTtf_q99_s','HRTtf_q99_5_s')) %>%
#         dplyr::filter(!(is.na(value))) %>%
#         dplyr::mutate(key = forcats::fct_rev(factor(key, levels=c('HRTtf_q0_5_s','HRTtf_q1_s','HRTtf_q10_s','HRTtf_q20_s','HRTtf_q50_s','HRTtf_q80_s','HRTtf_q90_s','HRTtf_q99_s','HRTtf_q99_5_s'),
#                                 labels=c('0.5%', '1%', '10%', '20%', '50%', '80%', '90%','99%', '99.5'))))

#     plot <- ggplot()+
#         stat_summary(data=forPlot_f, aes(x=key, y=value/60/60/24), color='darkgreen', fun = "median",  size = 8, geom = "point") +
#         stat_summary(data=forPlot_c, aes(x=key, y=value/60/60/24), fun = "median",  size = 8, geom = "point") +
#         stat_summary(data=forPlot_tf, aes(x=key, y=value/60/60/24), color='orange', fun = "median",  size = 8, geom = "point") +
#         ylab('HRT [dys]')+
#         xlab('Flow') +
#         scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10))
    
#     ggsave('cache/hrt_test.png', plot, width=15, height=10)
# }




makeHortonScalingFig <- function(hortonResults, model_combined){
    library(ggplot2)
    theme_set(theme_classic())

   # hortonResults <- rbind(hortonResults_barrier, hortonResults_nobarrier)
    
    #sum across huc basins
    hortonResults <- hortonResults %>%
        dplyr::group_by(StreamCalc) %>%
        dplyr::summarise(Vf_by_order_km3_fin = sum(Vf_by_order_km3_fin, na.rm=T))

    cumm_volume <- ggplot(hortonResults, aes(x=factor(StreamCalc), y=Vf_by_order_km3_fin)) +
        geom_point(size=8, color='#a39171') +
        #scale_color_manual(values=c('#9a8c98', '#a39171'))+
        scale_y_log10(labels = scales::label_log(base=10))+
        ylab(bquote(bold("Combined 1% flood volume ["*km^3*"]"))) +
        xlab('Stream Order')+
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none')


    model_combined <- model_combined %>%
        dplyr::mutate(Vf_q1_m3 = V_q1_m3 - (Wb_m*LengthKM*1000*Htf_q1_m))
    reach_volume <- ggplot(model_combined, aes(x=factor(StreamCalc), y=Vf_q1_m3)) +
        geom_boxplot(fill='#a39171') +
        scale_y_log10(labels = scales::label_log(base=10)) +
      #  scale_fill_manual(name='', labels=c('Barrier', 'No barrier'), values=c('#9a8c98', '#a39171'))+
        xlab('Stream Order') +
        ylab(bquote(bold("Individual reach 1% flood storage ["*m^3*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position=c(0.1, 0.9))
    
        layout <- "
        A
        B
    "

    comboPlot <- patchwork::wrap_plots(A=cumm_volume, B=reach_volume, design=layout)
    
    ggsave('cache/streamOrderPlot.png', comboPlot, width=8, height=15)
}




makeMapFig <- function(usgs_maps) {
    library(ggplot2)
    theme_set(theme_classic())

        #filter for the US only
    # CONUS boundary
    states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
    # states <- dplyr::filter(states, !(NAME %in% c('Alaska',
    #                                             'American Samoa',
    #                                             'Commonwealth of the Northern Mariana Islands',
    #                                             'Guam',
    #                                             'District of Columbia',
    #                                             'Puerto Rico',
    #                                             'United States Virgin Islands',
    #                                             'Hawaii'))) #remove non CONUS states/territories
        states <- dplyr::filter(states, NAME %in% c('Massachusetts', 'Rhode Island', 'New Hampshire', 'Maine', 'Vermont', 'Connecticut')) #remove non CONUS states/territories

    states <- sf::st_union(states) %>%
        sf::st_transform(crs=sf::st_crs(4326))

    codes_huc02 <- c('01')
    #read in all HUC4 basins and make a single shapefile
    basins_overall <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data', '/HUC2_', codes_huc02[1], '/NHDPLUS_H_0109_HU4_GDB/NHDPLUS_H_0109_HU4_GDB.gdb'), layer='WBDHU4') %>% dplyr::select('HUC4')
    # for(i in codes_huc02[-1]){
    #   basins <- sf::st_read(paste0('data/path_to_data/CONUS_ephemeral_data', '/HUC2_', i, '/WBD_', i, '_HU2_Shape/Shape/WBDHU', '4', '.shp')) %>% dplyr::select('huc4', 'name')
    #   basins_overall <- rbind(basins_overall, basins)
    # }

    basins_overall <- basins_overall %>%
        dplyr::filter(HUC4 %in% c('0108', '0109'))

#  basins_overall <- fixGeometries(basins_overall)
    
    map <- ggplot(usgs_maps) +
        geom_sf(data=states,
            color='black',
            size=1.25,
            alpha=0)+
        geom_sf(data=basins_overall,
            color='darkgreen',
            size=1.25,
            alpha=0)+
        geom_sf(aes(color='black'),
            size=5) +
        # annotate("text", x = -119, y = 25.1, label = bquote(.(nrow(usgs_maps))~"validation points"), size=6)+
        theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
        theme(legend.position='none')+
        theme(axis.text = element_text(family = "Futura-Medium", size=15),
            legend.title = element_text(face = "bold", size = 18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))+
        xlab('')+
        ylab('')
    ggsave('cache/validationMap.png', map, width=8, height=8)
}





makeRegulationFig <- function(model_combined, flow_perc){
    library(ggplot2)
    theme_set(theme_classic())

    model_combined <- model_combined %>%
        sf::st_drop_geometry() %>% 
        dplyr::mutate(Vf_q1_m3 = V_q1_m3 - (Wb_m*LengthKM*1000*Htf_q1_m),
                    Vf_q10_m3 = V_q10_m3 - (Wb_m*LengthKM*1000*Htf_q10_m),
                    Vf_q90_m3 = V_q90_m3 - (Wb_m*LengthKM*1000*Htf_q90_m),
                    flag = ifelse(regulatedDA_skm > 0, 'regulated', 'unregulated')) %>%
        dplyr::mutate(Vf_q1_m3 = ifelse(Vf_q1_m3 < 0, 0, Vf_q1_m3),
                    Vf_q10_m3 = ifelse(Vf_q10_m3 < 0, 0, Vf_q10_m3),
                    Vf_q90_m3 = ifelse(Vf_q90_m3 < 0, 0, Vf_q90_m3))
    
    num <- model_combined %>% 
        dplyr::group_by(StreamCalc, flag) %>% 
        dplyr::summarise(n=n()) %>%
        dplyr::mutate(filt = paste0(StreamCalc, '_', flag)) %>%
        dplyr::filter(n >= 100)

    model_combined <- model_combined %>%
        dplyr::mutate(filt = paste0(StreamCalc, '_', flag)) %>%
        dplyr::filter(filt %in% num$filt)

    orders_q1 <- ggplot(model_combined, aes(x=factor(StreamCalc), y=Vf_q1_m3/V_q1_m3, linetype=flag)) +
        geom_boxplot(alpha=0.3, fill='#66c2a5') +
        stat_summary(fun = mean,geom = 'point', aes(group=flag, pch=flag), color='#66c2a5', size=8, show.legend = F, position=position_dodge(width=0.5))+
        stat_summary(fun = mean,geom = 'line', aes(group=flag), color='#66c2a5', linewidth=1.25, show.legend = F, position=position_dodge(width=0.5))+
        scale_linetype(name='1% flood', labels=c(c('Upstream flow regulation', 'No upstream flow regulation')))+
        xlab('') +
        ylab('') +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position=c(0.8, 0.9),
            legend.title=element_text(size=20))
    
    orders_q10 <- ggplot(model_combined, aes(x=factor(StreamCalc), y=Vf_q10_m3/V_q10_m3, linetype=flag)) +
        geom_boxplot(alpha=0.3, fill='#fc8d62') +
        stat_summary(fun = mean,geom = 'point', aes(group=flag, pch=flag), color='#fc8d62', size=8, show.legend = F, position=position_dodge(width=0.5))+
        stat_summary(fun = mean,geom = 'line', aes(group=flag), color='#fc8d62', linewidth=1.25, show.legend = F, position=position_dodge(width=0.5))+
        scale_linetype(name='10% flood', labels=c(c('Upstream flow regulation', 'No upstream flow regulation')))+
        xlab('Stream Order') +
        ylab('Percent water in the floodplain') +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position=c(0.8, 0.9),
            legend.title=element_text(size=20))

    orders_q90 <- ggplot(model_combined, aes(x=factor(StreamCalc), y=Vf_q90_m3/V_q90_m3, linetype=flag)) +
        geom_boxplot(alpha=0.3, fill='#8da0cb') +
        stat_summary(fun = mean,geom = 'point', aes(group=flag, pch=flag), color='#8da0cb', size=8, show.legend = F, position=position_dodge(width=0.5))+
        stat_summary(fun = mean,geom = 'line', aes(group=flag), color='#8da0cb', linewidth=1.25, show.legend = F, position=position_dodge(width=0.5))+
        scale_linetype(name='90% flood', labels=c(c('Upstream flow regulation', 'No upstream flow regulation')))+
        xlab('Stream Order') +
        ylab('') +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position=c(0.8, 0.9),
            legend.title=element_text(size=20))

    forPlot <- model_combined %>%
        dplyr::mutate(q1_frac = Vf_q1_m3/V_q1_m3,
                    q10_frac = Vf_q10_m3/V_q10_m3,
                    q90_frac = Vf_q90_m3/V_q90_m3) %>%
        tidyr::gather(key=key, value=value, c(q1_frac, q10_frac, q90_frac)) %>%
        dplyr::group_by(StreamCalc, flag, key) %>%
        dplyr::summarise(value_fin = mean(value, na.rm=T))

    orders_qALL <- ggplot(forPlot, aes(x=factor(StreamCalc), y=value_fin, color=key, pch=flag, linetype=flag, group=interaction(key, flag))) +
        geom_line(linewidth=1.25)+
        geom_point(size=8) +
        scale_color_brewer('', palette='Set2', labels=c('1% flood', '10% flood', '90% flood'))+
        scale_shape('', labels=c('Upstream flow regulation', 'No upstream flow regulation'))+
        scale_linetype('', labels=c('Upstream flow regulation', 'No upstream flow regulation'))+
        ylim(0,1)+
        xlab('') +
        ylab('Percent water in floodplain') +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position=c(0.8, 0.8))

    layout <- "
        AB
        CD
    "

    comboPlot <- patchwork::wrap_plots(A=orders_qALL, B=orders_q1, C=orders_q10, D=orders_q90, design=layout)

    ggsave('cache/regulationPlot.png', comboPlot, width=14, height=14)
}









# makeValUSGSvolCombined <- function(val_USGSvols_combined){    
#     library(ggplot2)
#     theme_set(theme_classic())

#     #make inundated area val df
#     val_USGS_vol <- val_USGSvols_combined %>%
#         sf::st_drop_geometry() %>%
#         dplyr::filter(V_usgs_km3 > 0 & V_model_km3 > 0) %>%
#         dplyr::mutate(plot_var = factor(key_exdprob, c('V_q0_2', 'V_q0_5', 'V_q1', 'V_q2', 'V_q4', 'V_q10', 'V_q20', 'V_q50'), labels=c('0.2% flood', '0.5% flood', '1% flood', '2% flood', '4% flood', '10% flood', '20% flood', '50% flood')))

#     #get r2 and num reaches for each characteristic flood
#     df_r2 <- val_USGS_vol %>%
#         dplyr::group_by(plot_var) %>%
#         dplyr::group_modify(~ broom::glance(lm(log10(V_model_km3) ~ log10(V_usgs_km3), data = .x))) %>% 
#         dplyr::select(c('plot_var', 'r.squared'))
    
#     df_num <- val_USGS_vol %>%
#         dplyr::group_by(plot_var) %>%
#         dplyr::summarise(n=n())
    
#     lookup <- df_r2 %>%
#         dplyr::left_join(df_num, by='plot_var') %>%
#         dplyr::mutate(r2fin = sprintf("italic(R^2) == %.2f", r.squared),
#                     nfin = paste0(n, ' reaches'))

#     #plot
#     scatter_vol <- ggplot(val_USGS_vol, aes(x=V_usgs_km3*1e9, y=V_model_km3*1e9))+#, color=key_bias))+
#         geom_point(size=5, color='#a39171') +
#         geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
#         geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
#         # scale_color_manual(name='', labels=c('Model', 'Bias Corrected Model'), values=c('#81b29a', '#3d405b'))+
#         scale_x_log10(guide = "axis_logticks", limits=c(1e2, 1e7), breaks=c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
#         scale_y_log10(guide = "axis_logticks", limits=c(1e2, 1e7), breaks=c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
#         xlab(bquote(bold("USGS-hydrodynamic inundated volume ["*m^3*"]"))) +
#         ylab(bquote(bold("Modeled inundated volume ["*m^3*"]"))) +
#         theme(axis.title = element_text(face = 'bold', size=18),
#             axis.text = element_text(size=15),
#             strip.text = element_text(size=22),
#             axis.line = element_blank(),
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             strip.background = element_blank(),
#             panel.border = element_rect(colour = "black", fill = NA),
#             legend.position='bottom',
#             legend.text=element_text(size=15)) +
#         facet_wrap(vars(plot_var), nrow=2)+
#         geom_text(x=3, y=6, aes(label=r2fin), data=lookup, parse=TRUE, size=5)+
#         geom_text(x=3, y=6.5, aes(label=nfin), data=lookup, size=5)

#     ggsave('cache/validationUSGS_vol.png', scatter_vol, width=16, height=10)

#     return(val_USGSvols_combined)
# }