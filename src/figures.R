## Figures functions
## Craig Brinkerhoff
## Winter 2024

makeValFEMA <- function(val_FEMA_combined){
    library(ggplot2)
    theme_set(theme_classic())

    val_FEMA_combined <- val_FEMA_combined %>%
        sf::st_drop_geometry() %>%
        dplyr::filter(A_femaAEP_km2 > 0 & A_qFEMA_km2 > 0)
    
    lm <- lm(log(A_qFEMA_km2)~log(A_femaAEP_km2), data=val_FEMA_combined)

    scatter <- ggplot(val_FEMA_combined, aes(x=A_femaAEP_km2*1e6, y=A_qFEMA_km2*1e6, color=huc4))+
        geom_point(size=6, alpha=0.2) +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        scale_color_brewer(palette='Dark2')+
        geom_smooth(color='black', method='lm', se=F, linewidth=1.25) +
        scale_x_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e1, 1e7), breaks=c(1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels = scales::label_log(base=10))+
        xlab(bquote(bold("FEMA 100yr inundated area ["*m^2*"]"))) +
        ylab(bquote(bold("Modeled 100yr inundated area ["*m^2*"]"))) +
        theme(axis.title = element_text(face = 'bold', size=18),
            axis.text = element_text(size=15),
            legend.position='none')+
        annotate("text", x = 1e6, y = 10^1.5, label = bquote(.(nrow(val_FEMA_combined))*' rivers'), size=6)+
        annotate("text", x = 1e6, y = 10^1, label = bquote(r^2*': '*.(round(summary(lm)$r.squared,2))), size=6)
    
    ggsave('cache/validationFEMA.png', scatter, width=8, height=8)

    return(val_FEMA_combined)
}