## Figures functions
## Craig Brinkerhoff
## Fall 2025



makeMap <- function(conus_fin, month){
    library(sf)
    library(ggplot2)
    theme_set(theme_classic())

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

    #floodplain
	conus_fin <- conus_fin %>%
        dplyr::mutate(tau_col = dplyr::case_when(
            log10tau_hr <= -2 ~ '-2'
            ,log10tau_hr <= -1 ~ '-1'
            ,log10tau_hr <= 0 ~ '0'
            ,log10tau_hr <= 1 ~ '1'
            ,log10tau_hr <= 2 ~ '2'
            ,TRUE ~ '3'))

    conus_fin$tau_col <- factor(conus_fin$tau_col, levels = c('-2', '-1', '0', '1', '2', '3'))

    map_fp <- ggplot(conus_fin) +
        geom_sf(data=states,
            color='black',
            size=1,
            alpha=0)+
        geom_sf(aes(color= tau_col),#log10(TotDASqKm)),
            linewidth = 0.2,
            show.legend = "line") +
        scale_color_manual(values=c('#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51'),#palette='YlGnBu',
                        labels=c(bquote(10^-2), bquote(10^-1), bquote(10^0), bquote(10^1), bquote(10^2), bquote(10^3)),
                        name=bquote(bold(Mean~.(month)~tau[flood]~'[hr]')))+
        theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
        theme(legend.position=c(0.9, 0.2))+
        theme(text = element_text(family = "Futura-Medium"),
            legend.title = element_text(face = "bold", size = 18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))+
        labs(tag='A')+
        xlab('')+
        ylab('')

    #channel
	conus_fin <- conus_fin %>%
        dplyr::mutate(tau_col = dplyr::case_when(
            log10tau_channel_hr <= -2 ~ '-2'
            ,log10tau_channel_hr <= -1 ~ '-1'
            ,log10tau_channel_hr <= 0 ~ '0'
            ,log10tau_channel_hr <= 1 ~ '1'
            ,log10tau_channel_hr <= 2 ~ '2'
            ,TRUE ~ '3'))

    conus_fin$tau_col <- factor(conus_fin$tau_col, levels = c('-2', '-1', '0', '1', '2', '3'))

    map_ch <- ggplot(conus_fin) +
        geom_sf(data=states,
            color='black',
            size=1,
            alpha=0)+
        geom_sf(aes(color= tau_col),# linewidth=log10(TotDASqKm)),
            linewidth=0.2,
            show.legend = "line") +
        scale_linewidth(breaks = c(log10(10), log10(1000), log10(100000)), range=c(0.15, 1), guide='none')+
        scale_color_manual(values=c('#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51'), #palette='YlGnBu',
                        labels=c(bquote(10^-2), bquote(10^-1), bquote(10^0), bquote(10^1), bquote(10^2), bquote(10^3)),
                        name=bquote(bold(tau[channel]~'[hr]')))+
        theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
        theme(legend.position=c(0.9, 0.2))+
        labs(tag='B')+
        theme(text = element_text(family = "Futura-Medium"),
            legend.title = element_text(face = "bold", size = 18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))+
        xlab('')+
        ylab('')

    layout <- "
    A
    B
    "

    combo_plot <- patchwork::wrap_plots(A=map_fp, B=map_ch, design=layout)

    ggsave(paste0('cache/tau_fp_v_ch_',month,'.png'), combo_plot, width=13, height=15)
}




makeBasinTauMap <- function(basinShapefile){
    library(ggplot2)
    theme_set(theme_classic())

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

    basinShapefile$prob_lab <- ifelse(basinShapefile$prob == 0.02, "2% flood",
                                    ifelse(basinShapefile$prob == 0.1, "10% flood",
                                        ifelse(basinShapefile$prob == 0.5, "50% flood", NA)))
    basinShapefile$prob_lab <- factor(basinShapefile$prob_lab, levels=c('2% flood', '10% flood', '50% flood'))

    map <- ggplot(basinShapefile) +
        geom_sf(aes(fill=tau_hr_km), linewidth=0) +
        geom_sf(data=states,
            color='black',
            linewidth=0.75,
            alpha=0)+
        scale_fill_gradient(low='white', high='#723d46', name=bquote(tau[flood]~'['*hr^1*km^-1*']')) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position='bottom',
            axis.text=element_text(size=18),
            legend.text=element_text(size=16),
            legend.title=element_text(size=18),
            legend.key.spacing.x = unit(0.5, "cm"),
            strip.text.x = element_text(size = 22, color = "black", face = "bold.italic"),
            strip.text.y = element_text(size = 22, color = "black", face = "bold.italic"))+
        facet_wrap(vars(prob_lab), nrow=2)

    
    ggsave('cache/basinTauMap.png', map, width=12, height=10)

}



makeBasinRemovalMap <- function(basinShapefile){
    library(ggplot2)
    theme_set(theme_classic())

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

    basinShapefile$prob_lab <- ifelse(basinShapefile$prob == 0.02, "2% flood",
                                    ifelse(basinShapefile$prob == 0.1, "10% flood",
                                        ifelse(basinShapefile$prob == 0.5, "50% flood", NA)))
    basinShapefile$prob_lab <- factor(basinShapefile$prob_lab, levels=c('2% flood', '10% flood', '50% flood'))

    map <- ggplot(basinShapefile) +
        geom_sf(aes(fill=trap_eff_totalload), linewidth=0) +
        geom_sf(data=states,
            color='black',
            linewidth=0.75,
            alpha=0)+
        scale_fill_gradient(low='white', high='#344e41',name='% suspended sediment\ndeposited on floodplain', limits=c(0,15)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.position='bottom',
            axis.text=element_text(size=18),
            legend.text=element_text(size=16),
            legend.title=element_text(size=18),
            legend.spacing.x = unit(5, 'cm'),
            strip.text.x = element_text(size = 22, color = "black", face = "bold.italic"),
            strip.text.y = element_text(size = 22, color = "black", face = "bold.italic"))+
        facet_wrap(vars(prob_lab), nrow=2)

    
    ggsave('cache/basinRemovalMap.png', map, width=12, height=10)

}




makeReachBoxplotsFig <- function(combined_basinSummarySO){
    library(ggplot2)
    theme_set(theme_classic())

    forPlot <- combined_basinSummarySO %>%
        tidyr::gather(key=key, value=value, c('median_tau_hr_km', 'median_tau_channel_hr_km'))

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

    forPlot$prob_lab <- ifelse(forPlot$prob == 0.02, "2% flood",
                                    ifelse(forPlot$prob == 0.1, "10% flood",
                                        ifelse(forPlot$prob == 0.5, "50% flood", NA)))
    forPlot$prob_lab <- factor(forPlot$prob_lab, levels=c('2% flood', '10% flood', '50% flood'))

    boxplot <- ggplot(forPlot, aes(x=factor(StreamCalc), y=value, fill=key), color='black')+
        geom_boxplot() +
        scale_fill_manual(values=c('#86bbd8', '#b1b695'), name='', labels=c('Channel flow', 'Flood water')) +
        scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10))+
        xlab(bquote(bold(Stream~Order)))+
        ylab(bquote(bold(tau~'['*hr^1*km^-1*']'))) +
        theme(legend.position='bottom',
            axis.text=element_text(size=18),
            axis.title = element_text(size=20, face='bold'),
            legend.text=element_text(size=16),
            legend.title=element_text(size=18),
          #  legend.spacing.x = unit(5, 'cm'),
            strip.text.x = element_text(size = 22, color = "black", face = "bold.italic"),
            strip.text.y = element_text(size = 22, color = "black", face = "bold.italic"))+
        facet_wrap(vars(prob_lab), nrow=2)

    ggsave('cache/streamOrderFig.png', boxplot, width=12, height=10)
}







# makeMapProposal <- function(conus_fin, month){
#     library(sf)
#     library(ggplot2)
#     theme_set(theme_classic())

#     # CONUS boundary
#     states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
#     states <- dplyr::filter(states, !(NAME %in% c('Alaska',
#                                                 'American Samoa',
#                                                 'Commonwealth of the Northern Mariana Islands',
#                                                 'Guam',
#                                                 'District of Columbia',
#                                                 'Puerto Rico',
#                                                 'United States Virgin Islands',
#                                                 'Hawaii'))) #remove non CONUS states/territories

#     states <- sf::st_union(states) %>%
#         sf::st_transform(crs=sf::st_crs(4326))

#     #floodplain
# 	conus_fin <- conus_fin %>%
#         dplyr::mutate(tau_col = dplyr::case_when(
#             log10tau_hr <= -2 ~ '-2'
#             ,log10tau_hr <= -1 ~ '-1'
#             ,log10tau_hr <= 0 ~ '0'
#             ,log10tau_hr <= 1 ~ '1'
#             ,log10tau_hr <= 2 ~ '2'
#             ,TRUE ~ '3'))

#     conus_fin$tau_col <- factor(conus_fin$tau_col, levels = c('-2', '-1', '0', '1', '2', '3'))

#     map_fp <- ggplot(conus_fin) +
#         geom_sf(data=states,
#             color='black',
#             size=1,
#             alpha=0)+
#         geom_sf(aes(color= tau_col),#log10(TotDASqKm)),
#             linewidth = 0.2,
#             show.legend = "line") +
#         scale_color_manual(values=c('#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51'),#palette='YlGnBu',
#                         labels=c(bquote(10^-2), bquote(10^-1), bquote(10^0), bquote(10^1), bquote(10^2), bquote(10^3)),
#                         name=bquote(bold(Mean~.(month)~tau[flood]~'[hr]')))+
#         theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position=c(0.9, 0.2))+
#         theme(text = element_text(family = "Futura-Medium"),
#             legend.title = element_text(face = "bold", size = 18),
#             legend.text = element_text(family = "Futura-Medium", size = 18),
#             plot.tag = element_text(size=26,
#                                     face='bold'))+
#         labs(tag='A')+
#         xlab('')+
#         ylab('')

#     ggsave(paste0('cache/PROPOSAL_tau_fp_v_ch_',month,'.png'), map_fp, width=20, height=15)
# }








# makeMapSI <- function(conus_fin, month){
#     library(sf)
#     library(ggplot2)
#     theme_set(theme_classic())

#     # CONUS boundary
#     states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
#     states <- dplyr::filter(states, !(NAME %in% c('Alaska',
#                                                 'American Samoa',
#                                                 'Commonwealth of the Northern Mariana Islands',
#                                                 'Guam',
#                                                 'District of Columbia',
#                                                 'Puerto Rico',
#                                                 'United States Virgin Islands',
#                                                 'Hawaii'))) #remove non CONUS states/territories

#     states <- sf::st_union(states) %>%
#         sf::st_transform(crs=sf::st_crs(4326))

#     #floodplain
# 	conus_fin <- conus_fin %>%
#         dplyr::mutate(tau_col = dplyr::case_when(
#             log10tau_hr <= -2 ~ '-2'
#             ,log10tau_hr <= -1 ~ '-1'
#             ,log10tau_hr <= 0 ~ '0'
#             ,log10tau_hr <= 1 ~ '1'
#             ,log10tau_hr <= 2 ~ '2'
#             ,TRUE ~ '3'))

#     conus_fin$tau_col <- factor(conus_fin$tau_col, levels = c('-2', '-1', '0', '1', '2', '3'))

#     map <- ggplot(conus_fin) +
#         geom_sf(data=states,
#             color='black',
#             size=1,
#             alpha=0)+
#         geom_sf(aes(color= tau_col),#log10(TotDASqKm)),
#             linewidth = 0.2,
#             show.legend = "line") +
#         scale_color_manual(values=c('#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51'),#palette='YlGnBu',
#                         labels=c(bquote(10^-2), bquote(10^-1), bquote(10^0), bquote(10^1), bquote(10^2), bquote(10^3)),
#                         name=bquote(bold(Mean~.(month)~tau[flood]~'[hr]')))+
#         theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position='bottom')+
#         theme(plot.title = element_text(face = "italic", size = 26),
#             axis.text = element_text(size = 22),
#             plot.tag = element_text(size=26, face='bold'),
#             legend.position=c(0.9,0.15),
#             legend.text = element_text(size=18),
#             legend.title = element_text(size=22,face="bold"),
#             legend.spacing.y = unit(0.1, 'cm'))+
#         xlab('')+
#         ylab('')

#     ggsave(paste0('cache/tau_fp_v_ch_',month,'.png'), map, width=13, height=8)
# }





# makeTrappingMap <- function(conus_fin, month){
#     library(sf)
#     library(ggplot2)
#     theme_set(theme_classic())

#     # CONUS boundary
#     states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
#     states <- dplyr::filter(states, !(NAME %in% c('Alaska',
#                                                 'American Samoa',
#                                                 'Commonwealth of the Northern Mariana Islands',
#                                                 'Guam',
#                                                 'District of Columbia',
#                                                 'Puerto Rico',
#                                                 'United States Virgin Islands',
#                                                 'Hawaii'))) #remove non CONUS states/territories

#     states <- sf::st_union(states) %>%
#         sf::st_transform(crs=sf::st_crs(4326))

#     #Springtime (May)
# 	conus_fin <- conus_fin %>%
#         dplyr::mutate(frac_col = dplyr::case_when(
#             trap_eff <= 2.5 ~ '2.5'
#             ,trap_eff <= 5 ~ '5'
#             ,trap_eff <= 10 ~ '10'
#             ,trap_eff <= 20 ~ '20'
#             ,trap_eff <= 40, ~ '40'
#             ,trap_eff <= 100 ~ '100'
#             ,TRUE ~ NA))

#     conus_fin$frac_col <- factor(conus_fin$frac_col, levels = c('2.5', '5', '10', '20','40', '100', NA))

#     map_fp <- ggplot(conus_fin) +
#         geom_sf(data=states,
#             color='black',
#             size=1,
#             alpha=0)+
#         geom_sf(aes(color= frac_col),
#             linewidth = 0.2,
#             show.legend = "line") +
#         scale_color_manual(values=c('#dad7cd', '#a3b18a', '#588157', '#3a5a40', '#344e41', 'black'),
#                         labels=c('2.5%', '5%', '10%', '20%', '40%', '100%'),
#                         name = paste0(month, " flood sediment\ntrapping efficiency"))+
#         theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position=c(0.88, 0.2))+
#         theme(text = element_text(family = "Futura-Medium"),
#             legend.title = element_text(face = "bold", size = 18),
#             legend.text = element_text(family = "Futura-Medium", size = 18),
#             plot.tag = element_text(size=26,
#                                     face='bold'))+
#         labs(tag='A')+
#         xlab('')+
#         ylab('')

#     # #channel
# 	# conus_fin <- conus_fin %>%
#     #     dplyr::mutate(tau_col = dplyr::case_when(
#     #         log10tau_channel_hr <= -2 ~ '-2'
#     #         ,log10tau_channel_hr <= -1 ~ '-1'
#     #         ,log10tau_channel_hr <= 0 ~ '0'
#     #         ,log10tau_channel_hr <= 1 ~ '1'
#     #         ,log10tau_channel_hr <= 2 ~ '2'
#     #         ,TRUE ~ '3'))

#     # conus_fin$tau_col <- factor(conus_fin$tau_col, levels = c('-2', '-1', '0', '1', '2', '3'))

#     # map_ch <- ggplot(conus_fin) +
#     #     geom_sf(data=states,
#     #         color='black',
#     #         size=1,
#     #         alpha=0)+
#     #     geom_sf(aes(color= tau_col),# linewidth=log10(TotDASqKm)),
#     #         linewidth=0.2,
#     #         show.legend = "line") +
#     #     scale_linewidth(breaks = c(log10(10), log10(1000), log10(100000)), range=c(0.15, 1), guide='none')+
#     #     scale_color_manual(values=c('#264653', '#2a9d8f', '#8ab17d', '#e9c46a', '#f4a261', '#e76f51'), #palette='YlGnBu',
#     #                     labels=c(bquote(10^-2), bquote(10^-1), bquote(10^0), bquote(10^1), bquote(10^2), bquote(10^3)),
#     #                     name=bquote(bold(tau[channel]~'[hr]')))+
#     #     theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
#     #     theme(legend.position=c(0.9, 0.3))+
#     #     labs(tag='B')+
#     #     theme(text = element_text(family = "Futura-Medium"),
#     #         legend.title = element_text(face = "bold", size = 18),
#     #         legend.text = element_text(family = "Futura-Medium", size = 18),
#     #         plot.tag = element_text(size=26,
#     #                                 face='bold'))+
#     #     xlab('')+
#     #     ylab('')
    
#     # layout <- "
#     # A
#     # B
#     # "

#     # combo_plot <- patchwork::wrap_plots(A=map_fp, B=map_ch, design=layout)

#     ggsave(paste0('cache/trapping_', month, '.png'), map_fp, width=13, height=15)
# }







# makeSeasonalityFig <- function(seasonality_1, seasonality_5){
#     library(ggplot2)
#     theme_set(theme_classic())

#     forPlot <- rbind(seasonality_1, seasonality_5)
#     forPlot <- tidyr::gather(forPlot, key=key, value=value, c('log10tau_hr', 'log10tau_channel_hr'))

#     boxplot <- ggplot(forPlot, aes(x=factor(month), y=10^value, fill=key))+
#         geom_boxplot(color='black', linewidth=1) +
#         scale_fill_brewer(palette='Set2') +
#         scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10))+
#         theme(axis.title = element_text(size=22, face='bold'),
#             axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position='bottom')+
#         theme(text = element_text(family = "Futura-Medium"),
#             legend.title = element_text(face = "bold", size = 18),
#             legend.text = element_text(family = "Futura-Medium", size = 18),
#             plot.tag = element_text(size=26,
#                                     face='bold'))+
#         xlab(bquote(bold(tau~'[hr]')))+
#         ylab('Month')
    
#     ggsave('cache/seasonality.png', boxplot, width=9, height=9)
# }








# makeStreamOrderFig <- function(conus_fin){
#     library(ggplot2)
#     theme_set(theme_classic())

#     #stream order boxplots
#     forPlot <- conus_fin %>%
#         tidyr::gather(key=key, value=value, c('log10tau_hr', 'log10tau_channel_hr'))
#     boxplot <- ggplot(forPlot, aes(x=factor(StreamCalc), y=10^value, fill=key), color='black')+
#         geom_boxplot() +
#         scale_fill_manual(values=c('#86bbd8', '#b1b695'), name='', labels=c('Channel', 'Flood')) +
#         scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10))+
#         theme(axis.title = element_text(size=22, face='bold'),
#             axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position='bottom')+
#         theme(text = element_text(family = "Futura-Medium"),
#             legend.title = element_text(face = "bold", size = 18),
#             legend.text = element_text(family = "Futura-Medium", size = 18),
#             plot.tag = element_text(size=26,
#                                     face='bold'))+
#         xlab(bquote(bold(Stream~Order)))+
#         ylab(bquote(bold(tau~'[hr]'))) +
#         labs(tag="A")
    
#     #flood tau histogram
#     histFlood <- ggplot(conus_fin, aes(x=10^log10tau_hr)) +
#         geom_histogram(color='black', linewidth=1, fill='#b1b695', bins=50) +
#         scale_x_log10(guide = "axis_logticks",  labels = scales::label_log(base=10), limits=c(10^-5, 10^4))+
#         scale_y_continuous(breaks=c(50000, 150000, 250000), labels=c(50000, 150000, 250000))+
#         coord_cartesian(ylim=c(0, 300000))+
#         theme(axis.title = element_text(size=22, face='bold'),
#             axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position='bottom')+
#         theme(text = element_text(family = "Futura-Medium"),
#             legend.title = element_text(face = "bold", size = 18),
#             legend.text = element_text(family = "Futura-Medium", size = 18),
#             plot.tag = element_text(size=26,
#                                     face='bold'))+
#         ylab('')+
#         xlab(bquote(bold(tau[flood]~'[hr]'))) +
#         labs(tag='C')

#     #channel tau histogram
#     histChannel <- ggplot(conus_fin, aes(x=10^log10tau_channel_hr)) +
#         geom_histogram(color='black', linewidth=1, fill='#86bbd8', bins=50) +
#         scale_x_log10(guide = "axis_logticks",  labels = scales::label_log(base=10),  limits=c(10^-5, 10^4))+
#         scale_y_continuous(breaks=c(50000, 150000, 250000), labels=c(50000, 150000, 250000))+
#         coord_cartesian(ylim=c(0, 300000))+
#         theme(axis.title = element_text(size=22, face='bold'),
#             axis.text = element_text(family="Futura-Medium", size=20))+
#         theme(legend.position='bottom')+
#         theme(text = element_text(family = "Futura-Medium"),
#             legend.title = element_text(face = "bold", size = 18),
#             legend.text = element_text(family = "Futura-Medium", size = 18),
#             plot.tag = element_text(size=26,
#                                     face='bold'))+
#         ylab(bquote(bold(Count)))+
#         xlab(bquote(bold(tau[channel]~'[hr]'))) +
#         labs(tag='B')

#     #combo plot
#     layout <- "
#     AA
#     AA
#     BC
#     "

#     combo_plot <- patchwork::wrap_plots(A=boxplot, B=histChannel, C=histFlood, design=layout)

#     ggsave('cache/distributions.png', combo_plot, height=12, width=10)
# }





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
        annotate("text", x = 1e5, y = 1e10, label = bquote(bold(n:~.(nrow(full_predictions_Q)))), size=8, color='#a3b18a')+
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
        annotate("text", x = 1e3, y = 1e9, label = bquote(bold(n:~.(nrow(full_predictions_V)))), size=8, color='#588157')+
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

    #combine
    full_predictions <- full_predictions_Q %>%
        dplyr::left_join(full_predictions_V, by='.row') %>%
        dplyr::mutate(log10tau_hr = (V_m3 - Qexc_m3dy) + log10(24),
                    log10tau_hr.pred = (.pred.y - .pred.x) + log10(24)) %>%
        tidyr::drop_na()

    r2 <- round(summary(lm(log10tau_hr.pred~log10tau_hr, data=full_predictions))$r.squared,2)

    scatter_tau <- ggplot(full_predictions, aes(x=10^(log10tau_hr), y=10^(log10tau_hr.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(pch=21, size=8, color='black', fill='#3a5a40') +
        geom_smooth(method='lm', se=F, linewidth=2, color='black')+
        # coord_cartesian(xlim=c(1,20), ylim=c(1,20))+
        annotate("text", x = 10^-2, y = 10^5, label = bquote(bold(n:~.(nrow(full_predictions)))), size=8, color='#3a5a40')+ #7.5 30
        annotate("text", x = 10^-2, y = 10^4.5, label = bquote(bold(r^2*':'~.(r2))), size=8, color='#3a5a40')+ #7.5 28.5
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e-4, 1e5))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e-4, 1e5))+
        xlab(bquote(bold('Observed'~tau[flood]~'[hr]')))+
        ylab(bquote(bold('Modeled'~tau[flood]~'[hr]')))+
        theme(axis.title=element_text(face='bold', size=18),
            axis.text = element_text(size=15),
            panel.border = element_rect(colour = "black", fill = NA),
            plot.tag = element_text(size=24,
                                    face='bold'))

    layout <- "
        ABC
    "

    comboPlot <- patchwork::wrap_plots(A=scatter_Q, B=scatter_V, C=scatter_tau, design=layout)

    ggsave('cache/validationML.png', comboPlot, width=20, height=8)

    return(full_predictions)
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
    gage_shp <- gage_df %>%
        sf::st_as_sf(coords=c('lon', 'lat'), crs='epsg:4326')

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








# experimentDOC <- function(raymondSaiers_data, conus_fin){
#     library(ggplot2)
#     theme_set(theme_classic())

#     conus_fin <- conus_fin %>%
#         sf::st_drop_geometry()%>%
#         dplyr::inner_join(raymondSaiers_data, by=c('GageID'='site_no')) %>%
#         dplyr::filter(!(is.na(GageID))) %>%
#         dplyr::select(c('log10tau_hr', 'mean_DOC_mgL'))

#     # raymondSaiers_data <- raymondSaiers_data %>%
#     #     dplyr::group_by(site_no) %>%
#     #     dplyr::summarise(tau_dy = mean(V_m3)/mean(Qexc_m3dy),
#     #                     mean_DOC_mgL = mean(mean_DOC_mgL))


#     plot <- ggplot(conus_fin, aes(x=10^(log10tau_hr), y=mean_DOC_mgL))+
#         geom_point(pch=21, color='black', fill = 'green', size=8, linewidth=1) +
#         scale_x_log10(guide = "axis_logticks",  labels = scales::label_log(base=10))+
#         scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10))+
#         ylab(bquote(DOC~'['*mg^1*L^-1*']')) +
#         xlab(bquote(tau[flood]~'[hr]'))

#     ggsave('cache/DOC.png', plot)
# }





# makeVIPPlot <- function(final_model_Q, final_model_V){
#     library(tidymodels)

#     theme_set(theme_classic())

#     scores_Q <- vip::vi(final_model_Q, method="model", type='gain') #for lightgbm, this is the coefficient magnitude for each feature, re-scaled to relative percentages. Gain "represents fractional contribution of each feature to the model based on the total gain of this feature's splits"
#             #https://koalaverse.github.io/vip/reference/vi_model.html
#     scores_Q <- scores_Q[1:15,]

#     scores_V <- vip::vi(final_model_V, method="model", type='gain') #for lightgbm, this is the coefficient magnitude for each feature, re-scaled to relative percentages. Gain "represents fractional contribution of each feature to the model based on the total gain of this feature's splits"
#             #https://koalaverse.github.io/vip/reference/vi_model.html
#     scores_V <- scores_V[1:15,]
    
#     plot_Q <- ggplot(scores_Q, aes(x=Importance, y=forcats::fct_reorder(Variable, Importance)))+
#         geom_col(fill='#c1121f', color='black', size=1.2) +
#         scale_fill_brewer(palette='Set2', name='')+
#         xlab(bquote(bold(Q[flood]~'model feature importance'))) +
#         ylab('') +
#         coord_cartesian(xlim=c(0,1))+
#         theme(axis.title = element_text(size=22, face='bold'),
#             axis.text = element_text(size=18,face='bold'),
#             legend.position=c(0.75, 0.1),
#             legend.text = element_text(size=22))

#     plot_V <- ggplot(scores_V, aes(x=Importance, y=forcats::fct_reorder(Variable, Importance)))+
#         geom_col(fill='#669bbc', color='black', size=1.2) +
#         scale_fill_brewer(palette='Set2', name='')+
#         xlab(bquote(bold(V[flood]~'model feature importance'))) +
#         ylab('') +
#         coord_cartesian(xlim=c(0,1))+
#         theme(axis.title = element_text(size=22, face='bold'),
#             axis.text = element_text(size=18,face='bold'),
#             legend.position=c(0.75, 0.1),
#             legend.text = element_text(size=22))
    
#     layout <- "
#         AB
#     "

#     comboPlot <- patchwork::wrap_plots(A=plot_Q, B=plot_V, design=layout)

#     ggsave('cache/vip.png', comboPlot, width=20, height=10)

#     return('cache/vip.png')
# }
