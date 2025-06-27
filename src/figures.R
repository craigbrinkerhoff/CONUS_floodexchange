## Figures functions
## Craig Brinkerhoff
## Winter 2024




makeCalculationValFig <- function(gageVolume_val_combined, gageQexc_val){    
    library(ggplot2)
    theme_set(theme_classic())

    df <- gageQexc_val %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(site_no, flag) %>%
        dplyr::summarise(Qexc_m3dy = mean(Qexc_m3dy, na.rm=T)) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        tidyr::pivot_wider(id_cols=site_no, names_from=flag, values_from=c('Qexc_m3dy'))

    lm <- lm(log10(model)~log10(obs), data=df)
    plot1 <- ggplot(df, aes(x=obs, y= model)) +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_point(fill='#c1121f', pch=21, size=10, color='black') +
        geom_smooth(method='lm', se=F, color='black', linewidth=2)+
        annotate("text", x = 10^4.5, y = 10^7.5, label = bquote(R^2*'= '*.(round(summary(lm)$r.squared,2))), size=6)+
        annotate("text", x = 10^4.5, y = 10^7.25, label = bquote(n*'= '*.(nrow(df))*' gages'), size=6)+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(10^3.75, 10^7.5))+
        scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10), limits=c(10^3.75, 10^7.5))+
        ylab(bquote(bold('Estimated mean floodplain flux ['~m^3/dy~']')))+
        xlab(bquote(bold('Observed mean floodplain flux ['~m^3/dy~']')))+
        labs(tag='A')+
        theme(axis.title = element_text(size=20, face='bold'))+
        theme(legend.position='none',
            panel.border = element_rect(colour = "black", fill = NA),
            axis.line = element_blank())+
        theme(axis.text = element_text(family = "Futura-Medium", size=18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))

    forPlot <- gageVolume_val_combined %>%
        dplyr::mutate(V_usgs_m3 = ifelse(V_usgs_m3 == 0, NA, V_usgs_m3),
                    V_model_m3 = V_m3) %>%
        tidyr::drop_na(V_model_m3, V_usgs_m3)

    #plot
    lm <- lm(log10(V_model_m3)~log10(V_usgs_m3), data=forPlot)
    plot2 <- ggplot(forPlot, aes(x=V_usgs_m3, y=V_model_m3))+
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_point(pch=21, size=10, color='black', fill='#669bbc') +
        geom_smooth(color='black', method='lm', se=F, linewidth=2) +
        annotate("text", x = 10^4, y = 10^8, label = bquote(R^2*'= '*.(round(summary(lm)$r.squared,2))), size=6)+
        annotate("text", x = 10^4, y = 10^7.65, label = bquote(n*'= '*.(nrow(forPlot))*' reaches'), size=6)+
        scale_x_log10(guide = "axis_logticks", limits=c(1e3, 1e8), breaks=c(1e3, 1e4, 1e5, 1e6, 1e7,1e8), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e3, 1e8), breaks=c(1e3, 1e4, 1e5, 1e6, 1e7,1e8), labels = scales::label_log(base=10))+
        xlab(bquote(bold("USGS-hydrodynamic inundated volume ["*m^3*"]"))) +
        ylab(bquote(bold("Estimated inundated volume ["*m^3*"]"))) +
        labs(tag='B')+
        theme(axis.title = element_text(size=20, face='bold'))+
        theme(legend.position='none',
            panel.border = element_rect(colour = "black", fill = NA),
            axis.line = element_blank())+
        theme(axis.text = element_text(family = "Futura-Medium", size=18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))
    
    layout <- "
    AB
    "

    combo_plot <- patchwork::wrap_plots(A=plot1, B=plot2, design=layout)

    ggsave('cache/validation_calculation.png', combo_plot, width=16, height=8)

    return(forPlot)
}






makeMLValFig <- function(trainedModel_Q, trainedModel_V){
    library(tidymodels)
    library(ggplot2)
    theme_set(theme_classic())

    # Q model
    full_predictions_Q <- data.frame()
    for(i in 1:length(trainedModel_Q)){
        results <- trainedModel_Q[[i]]$summary %>% collect_predictions()
        results$fold <- i
        results <- results %>%
            dplyr::relocate(fold)
        full_predictions_Q <- rbind(full_predictions_Q, results)
    }

    #V model
    full_predictions_V <- data.frame()
    for(i in 1:length(trainedModel_V)){
        results <- trainedModel_V[[i]]$summary %>% collect_predictions()
        results$fold <- i
        results <- results %>%
            dplyr::relocate(fold)
        full_predictions_V <- rbind(full_predictions_V, results)
    }

    r2_Q <- round(summary(lm(.pred~Qexc_m3dy, data=full_predictions_Q))$r.squared,2)
    scatter_Q <- ggplot(full_predictions_Q, aes(x=10^(Qexc_m3dy), y=10^(.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(pch=21, size=8, color='black', fill='#a3b18a') +
        geom_smooth(method='lm', se=F, linewidth=2, color='black')+
        annotate("text", x = 1e5, y = 1e10, label = bquote(bold(n:~.(nrow(full_predictions_Q))~'gages')), size=8, color='#a3b18a')+
        annotate("text", x = 1e5, y = 10^9.5, label = bquote(bold(r^2*':'~.(r2_Q))), size=8, color='#a3b18a')+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e3, 1e10))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e3, 1e10))+
        xlab(bquote(bold('Observed '~Q[flood]~'['~m^3*dy^-1~']')))+
        ylab(bquote(bold('Modeled '~Q[flood]~'['~m^3*dy^-1~']')))+
        theme(axis.title=element_text(face='bold', size=18),
            axis.text = element_text(size=15),
            panel.border = element_rect(colour = "black", fill = NA),
            plot.tag = element_text(size=24,
                                    face='bold'))

    r2_V <- round(summary(lm(.pred~V_m3, data=full_predictions_V))$r.squared,2)
    scatter_V <- ggplot(full_predictions_V, aes(x=10^(V_m3), y=10^(.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(pch=21, size=8, color='black', fill='#588157') +
        geom_smooth(method='lm', se=F, linewidth=2, color='black')+
        annotate("text", x = 1e3, y = 1e9, label = bquote(bold(n:~.(nrow(full_predictions_V))~'gages')), size=8, color='#588157')+
        annotate("text", x = 1e3, y = 10^8.5, label = bquote(bold(r^2*':'~.(r2_V))), size=8, color='#588157')+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e1, 1e9))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e1, 1e9))+
        xlab(bquote(bold('Observed '~V[flood]~'['~m^3~']')))+
        ylab(bquote(bold('Modeled '~V[flood]~'['~m^3~']')))+
        theme(axis.title=element_text(face='bold', size=18),
            axis.text = element_text(size=15),
            panel.border = element_rect(colour = "black", fill = NA),
            plot.tag = element_text(size=24,
                                    face='bold'))

    #combine (using full join for twice repated so when brought together ~4 predictions at each site)
    full_predictions <- full_predictions_Q %>%
        dplyr::full_join(full_predictions_V, by='.row') %>%
        dplyr::mutate(tau_dy = 10^(V_m3 / Qexc_m3dy),
                    log10tau_dy.pred = (.pred.y / .pred.x)) %>%
        dplyr::mutate(transform_bias =  mean(10^(log10(tau_dy)-log10tau_dy.pred), na.rm=T)) %>% #eq 9.24 from https://pubs.usgs.gov/tm/04/a03/tm4a3.pdf (also see https://nrtwq.usgs.gov/co/methods/)
        dplyr::mutate(tau_dy.pred = 10^(log10tau_dy.pred) * transform_bias) %>%
        tidyr::drop_na()

    r2 <- round(summary(lm(tau_dy.pred~tau_dy, data=full_predictions))$r.squared,2)

    scatter_tau <- ggplot(full_predictions, aes(x=(tau_dy), y=(tau_dy.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(pch=21, size=8, color='black', fill='#3a5a40') +
        geom_smooth(method='lm', se=F, linewidth=2, color='black')+
        coord_cartesian(xlim=c(1,30), ylim=c(1,30))+
        annotate("text", x = 7.5, y = 30, label = bquote(bold(n:~.(nrow(full_predictions))~'gages')), size=8, color='#3a5a40')+
        annotate("text", x = 7.5, y = 28.5, label = bquote(bold(r^2*':'~.(r2))), size=8, color='#3a5a40')+
        xlab(bquote(bold('Observed'~tau[flood]~'[dys]')))+
        ylab(bquote(bold('Modeled'~tau[flood]~'[dys]')))+
        theme(axis.title=element_text(face='bold', size=18),
            axis.text = element_text(size=15),
            panel.border = element_rect(colour = "black", fill = NA),
            plot.tag = element_text(size=24,
                                    face='bold'))

    layout <- "
        ABC
    "

    comboPlot <- patchwork::wrap_plots(A=scatter_Q, B=scatter_V, C=scatter_tau, design=layout)

    ggsave('cache/validationML.png', comboPlot, width=20, height=9)

    return('cache/validationML.png')
}




makeVIPPlot <- function(final_model_Q, final_model_V){
    library(tidymodels)

    theme_set(theme_classic())

    scores_Q <- vip::vi(final_model_Q, method="model", type='gain') #for lightgbm, this is the coefficient magnitude for each feature, re-scaled to relative percentages. Gain "represents fractional contribution of each feature to the model based on the total gain of this feature's splits"
            #https://koalaverse.github.io/vip/reference/vi_model.html
    scores_Q <- scores_Q[1:15,]

    scores_V <- vip::vi(final_model_V, method="model", type='gain') #for lightgbm, this is the coefficient magnitude for each feature, re-scaled to relative percentages. Gain "represents fractional contribution of each feature to the model based on the total gain of this feature's splits"
            #https://koalaverse.github.io/vip/reference/vi_model.html
    scores_V <- scores_V[1:15,]
    
    plot_Q <- ggplot(scores_Q, aes(x=Importance, y=forcats::fct_reorder(Variable, Importance)))+
        geom_col(fill='#c1121f', color='black', size=1.2) +
        scale_fill_brewer(palette='Set2', name='')+
        xlab(bquote(bold(Q[flood]~'model feature importance'))) +
        ylab('') +
        coord_cartesian(xlim=c(0,1))+
        theme(axis.title = element_text(size=22, face='bold'),
            axis.text = element_text(size=18,face='bold'),
            legend.position=c(0.75, 0.1),
            legend.text = element_text(size=22))

    plot_V <- ggplot(scores_V, aes(x=Importance, y=forcats::fct_reorder(Variable, Importance)))+
        geom_col(fill='#669bbc', color='black', size=1.2) +
        scale_fill_brewer(palette='Set2', name='')+
        xlab(bquote(bold(V[flood]~'model feature importance'))) +
        ylab('') +
        coord_cartesian(xlim=c(0,1))+
        theme(axis.title = element_text(size=22, face='bold'),
            axis.text = element_text(size=18,face='bold'),
            legend.position=c(0.75, 0.1),
            legend.text = element_text(size=22))
    
    layout <- "
        AB
    "

    comboPlot <- patchwork::wrap_plots(A=plot_Q, B=plot_V, design=layout)

    ggsave('cache/vip.png', comboPlot, width=20, height=10)

    return('cache/vip.png')
}




makeGageMap <- function(gage_df) {
    sf::sf_use_s2(FALSE)
    library(ggplot2)
    theme_set(theme_classic())

    #filter for the US only
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
        sf::st_transform(crs=sf::st_crs(4326))

    #build gage shapefile
    gage_shp <- sf::st_as_sf(gage_df, coords=c('lon', 'lat'), crs='epsg:4326')

    map <- ggplot(gage_shp) +
        geom_sf(data=states,
            color='black',
            size=1.25,
            alpha=0)+
        geom_sf(fill='#3a5a40',
            color='black',
            pch=21,
            size=3) +
        annotate("text", x = -119, y = 25.1, label = bquote(.(nrow(gage_shp))~"gages"), size=6)+
        theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
        theme(legend.position='bottom')+
        theme(text = element_text(family = "Futura-Medium"),
            legend.title = element_text(face = "bold", size = 18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))+
        xlab('')+
        ylab('')
    ggsave('cache/gageMap.png', map, width=10, height=7)

    return(gage_shp)
}