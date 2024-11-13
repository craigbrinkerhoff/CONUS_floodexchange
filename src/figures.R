## Figures functions
## Craig Brinkerhoff
## Winter 2024

makeValFEMA <- function(val_FEMA_combined){
    library(ggplot2)
    theme_set(theme_classic())

    val_FEMA_combined <- val_FEMA_combined %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(A_femaAEP_km2 > 0 & A_qFEMA_km2 > 0) %>%
        dplyr::mutate(huc4 = substr(huc8, 1, 4))

    scatter <- ggplot(val_FEMA_combined, aes(x=A_femaAEP_km2, y=A_qFEMA_km2, color=huc4))+
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_point(size=6) +
        scale_color_brewer(palette='Dark2')+
        #geom_smooth(color='black', method='lm', se=F, linewidth=2) +
        scale_x_continuous(breaks=c(2,5,8,11), limits=c(0,13))+
        scale_y_continuous(breaks=c(2,5,8,11), limits=c(0,13))+
        ylab(bquote(bold("FEMA 100yr inundated area ["*km^2*"]"))) +
        xlab(bquote(bold("Modeled 100yr inundated area ["*km^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none')+
        annotate("text", x = 2.5, y = 9, label = paste0(nrow(val_FEMA_combined), ' HUC8 basins'), size=6)
    
    ggsave('cache/validationFEMA.png', scatter, width=8, height=8)

    return(val_FEMA_combined)
}