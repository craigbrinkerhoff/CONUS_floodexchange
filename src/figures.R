## Figures functions
## Craig Brinkerhoff
## Winter 2024

makeValFEMA <- function(val_FEMA_combined, gage_combined, BHGdata){
    library(ggplot2)
    theme_set(theme_classic())

    #make inundated area val df
    val_FEMA_area <- val_FEMA_combined$for_area %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(A_femaAEP_km2 > 0 & A_qFEMA_km2 > 0)

    val_FEMA_wse <- val_FEMA_combined$for_depth %>%
        dplyr::mutate(H_qFEMA_m = (Htf_qFEMA_m+groundElevLookup_m)) %>%
        dplyr::filter(H_qFEMA_m > 0 & elev_m > 0)

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
    lm_inn_v <- lm(log(H_qFEMA_m)~log(elev_m), data=val_FEMA_wse)
    lm_flood <- lm(log(Af_qFEMA_valmodel_km2)~log(Af_femaAEP_km2), data=val_FEMA_area_fp)

    scatter_area <- ggplot(val_FEMA_area, aes(x=A_femaAEP_km2*1e6, y=A_qFEMA_km2*1e6, color=huc4))+
        geom_point(size=6, alpha=0.2) +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        scale_color_brewer(palette='Dark2')+
        geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab('') +
        ylab(bquote(bold("Modeled 100yr inundated area ["*m^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none',
            plot.title = element_text(size=22)) +
        labs(title='A')+
        annotate("text", x = 1e5, y = 10^1.5, label = bquote(.(nrow(val_FEMA_area))*' river reaches'), size=6)+
        annotate("text", x = 1e5, y = 10^1, label = bquote(r^2*': '*.(round(summary(lm_inn_a)$r.squared,2))), size=6)
    
    scatter_wse <- ggplot(val_FEMA_wse, aes(x=elev_m, y=H_qFEMA_m))+
        geom_point(size=6, alpha=0.2) +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        xlab(bquote(bold("FEMA 100yr wse ["*m*"]"))) +
        ylab(bquote(bold("Modeled 100yr wse ["*m*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none',
            plot.title = element_text(size=22))+
        labs(title='C')+
        annotate("text", x = 350, y = 50, label = bquote(.(nrow(val_FEMA_wse))*' river cross-sections'), size=6)+
        annotate("text", x = 350, y = 10, label = bquote(r^2*': '*.(round(summary(lm_inn_v)$r.squared,2))), size=6)
    

    scatter_area_fp <- ggplot(val_FEMA_area_fp, aes(x=Af_femaAEP_km2*1e6, y=Af_qFEMA_valmodel_km2*1e6, color=huc4))+
        geom_point(size=6) +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        scale_color_brewer(palette='Dark2')+
        geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab(bquote(bold("FEMA 100yr inundated area ["*m^2*"]"))) +
        ylab(bquote(bold("Modeled 100yr inundated area ["*m^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none',
            plot.title = element_text(size=22)) +
        labs(title='B')+
        annotate("text", x = 1e5, y = 10^1.5, label = bquote(.(nrow(val_FEMA_area_fp))*' river reaches'), size=6)+
        annotate("text", x = 1e5, y = 10^1, label = bquote(r^2*': '*.(round(summary(lm_flood)$r.squared,2))), size=6)
    
    layout <- "
        AC
        BD
    "

    comboPlot <- patchwork::wrap_plots(A=scatter_area, B=scatter_area_fp, C= scatter_wse, design=layout)

    ggsave('cache/validationFEMA.png', comboPlot, width=12, height=12)

    return(val_FEMA_combined)
}





upscalingFig <- function(df, huc4){
    library(ggplot2)

    theme_set(theme_classic())
  
    plot <- ggplot(df$df, aes(x=DA_skm, y=value)) +
        geom_point(size=5, color='#153131')+
        geom_smooth(method='lm', se=F, color='black',linewidth=1.25) +
        scale_x_log10()+
        scale_y_log10() +
        ylab(bquote(bold("Water level ["*m*"]"))) +
        xlab(bquote(bold("Drainage Area ["*km^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            plot.title = element_text(face = 'bold', size=22)) +
        facet_wrap(~factor(name, levels=c('qFEMAflood_Htf_m', 'q1flood_Htf_m', 'q10flood_Htf_m', 'q25flood_Htf_m', 'q50flood_Htf_m', 'q75flood_Htf_m'),
                                labels=c('100yr AEP flood', '1% flood', '10% flood', '25% flood', '50% flood', '75% flood'))) +
        geom_text(data = df$models_r2, aes(x = 10^0.8, y = 5, label = paste0('rsq: ',round(df$models$rsq,2))), size=6) +
        geom_text(data = df$models_r2, aes(x = 10^0.8, y = 3.5, label = paste0(df$models$n_gages, ' gages')), size=6) +
        theme(strip.text.x = element_text(size = 20))
  
  ggsave(paste0('cache/upscalingModel_', huc4, '.png'), plot, width=12, height=10)

  return(paste0('cache/upscalingModel_', huc4, '.png'))
}





# hortonScalingFig <- function(hortonResults){
#     library(ggplot2)
#     theme_set(theme_classic())

#     plot_volume <- ggplot(hortonResults, aes(x=StreamCalc, y=Vf_by_order_km3_fin)) +
#         geom_point(size=5) +
#         scale_y_log10(labels = scales::label_log(base=10))+
#         scale_x_log10()+
#         ylab(bquote(bold("Floodplain inundated volume ["*km^3*"]"))) +
#         xlab('Stream Order')+
#         theme(axis.title = element_text(face = 'bold', size=18),
#             axis.text = element_text(size=15),
#             legend.position='none')
    
#     plot_area <- ggplot(hortonResults, aes(x=StreamCalc, y=Af_by_order_km2_fin)) +
#         geom_point(size=5) +
#         scale_y_log10(labels = scales::label_log(base=10))+
#         scale_x_log10()+
#         ylab(bquote(bold("Floodplain inundated area ["*km^2*"]"))) +
#         xlab('')+
#         theme(axis.title = element_text(face = 'bold', size=18),
#             axis.text = element_text(size=15),
#             legend.position='none')
    

#     layout <- "
#         A
#         B
#     "

#     comboPlot <- patchwork::wrap_plots(A=plot_area, B=plot_volume, design=layout)
    
#     ggsave('cache/scalingPlot.png', comboPlot, width=8, height=15)
# }