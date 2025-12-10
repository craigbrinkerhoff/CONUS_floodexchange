## Figures functions
## Craig Brinkerhoff
## Fall 2025




#' makeBasinExchangeMap
#'
#' Makes basin flood discharge map
#' 
#' @param combined_basinSummary dataframe of basin Qflood and total Q and fraction of flow in floodplain
#'
#' @return figure written to file
makeBasinExchangeMap <- function(combined_basinSummary){
    library(sf)
    library(ggplot2)
    theme_set(theme_classic())

    #CONUS boundary
    states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
    states <- dplyr::filter(states, !(NAME %in% c('Alaska',
                                                'American Samoa',
                                                'Commonwealth of the Northern Mariana Islands',
                                                'Guam',
                                                'District of Columbia',
                                                'Puerto Rico',
                                                'United States Virgin Islands',
                                                'Hawaii'))) #remove non CONUS states/territories

    #regional basin ids
    hucs <- c('0101', '0102', '0103', '0104', '0105', '0106', '0107', '0108', '0109', '0110',
        '0202', '0203', '0204', '0205', '0206', '0207', '0208',
        '0301', '0302', '0303', '0304', '0305', '0306', '0307', '0308', '0309', '0310', '0311', '0312', '0313', '0314', '0315', '0316', '0317', '0318',
        '0401', '0402', '0403', '0404', '0405', '0406', '0407', '0408', '0409', '0410', '0411', '0412', '0413', '0414', '0420', '0427', '0429', '0430', '0431',
        '0501', '0502', '0503', '0504', '0505', '0506', '0507', '0508', '0509', '0510', '0511', '0512', '0513', '0514',
        '0601', '0602', '0603', '0604',
        '0701', '0702','0703', '0704', '0705', '0706', '0707', '0708', '0709', '0710', '0711', '0712', '0713', '0714',
        '0801', '0802', '0803', '0804', '0805', '0806', '0807', '0808', '0809',
        '0901', '0902', '0903', '0904',
        '1002', '1003', '1004', '1005', '1006', '1007', '1008', '1009', '1010', '1011', '1012', '1013', '1014', '1015', '1016', '1017', '1018', '1019', '1020', '1021', '1022', '1023', '1024', '1025', '1026', '1027', '1028', '1029', '1030',
        '1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108', '1109', '1110', '1111', '1112', '1113', '1114',
        '1201', '1202', '1203', '1204', '1205', '1206', '1207', '1208', '1209','1210','1211',
        '1301','1302','1303','1304','1305','1306','1307','1308','1309',
        '1401','1402','1403','1404','1405','1406','1407','1408',
        '1501','1502','1503','1504','1505','1506','1507','1508',
        '1601','1602','1603','1604','1605','1606',
        '1701','1702','1703','1704','1705','1706','1707','1708','1709','1710','1711','1712',
        '1801','1802','1803','1804','1805','1806','1807','1808','1809','1810')

    #wrangle basins
    basins <- sf::st_read('data/HUC4s.shp') %>%
        dplyr::filter(huc4 %in% hucs) %>%
        dplyr::filter(name != 'Lake Erie')

    states <- sf::st_union(states) %>%
        sf::st_transform(crs=sf::st_crs(basins))

    basins <- basins %>%
        sf::st_intersection(states)

    #add model results
    combined_basinSummary <- combined_basinSummary$out %>%
        dplyr::select(c('huc4','prob', 'total_exchange_frac')) %>%
        tidyr::drop_na(total_exchange_frac) %>%
        dplyr::mutate(total_exchange_frac = ifelse(total_exchange_frac > 1, 1, total_exchange_frac)) #a handful of basins are just above 1 because the numbers come from independent models, so we just round them off

    basins <- basins %>%
        dplyr::left_join(combined_basinSummary, by='huc4') %>%
        tidyr::drop_na(total_exchange_frac)
    
    basins$prob_lab <- ifelse(basins$prob == 0.02, "2% flood",
                                    ifelse(basins$prob == 0.1, "10% flood",
                                        ifelse(basins$prob == 0.2, "20% flood",
                                            ifelse(basins$prob == 0.5, "50% flood", NA))))
    basins$prob_lab <- factor(basins$prob_lab, levels=c('2% flood', '10% flood', '20% flood', '50% flood'))

    #plot
    map <- ggplot(basins) +
        geom_sf(aes(fill=total_exchange_frac*100), color='black', linewidth=0.5) +
        geom_sf(data=states,color='black',linewidth=1,alpha=0) +
        scale_fill_gradient(name='% of flow exchanged\nwith floodplain', low='white', high='#003049', limits=c(0,100))+
        theme(legend.position='bottom',
            strip.text.x = element_text(size = 18),
            strip.text.y = element_text(size = 16),
            axis.text = element_text(size=14),
            legend.title=element_text(size=18),
            legend.text = element_text(size=14))+
        guides(color='none')+
        facet_wrap(vars(prob_lab), nrow=2) +
        theme(strip.background = element_blank(),
            panel.border=element_rect(colour="black",size=1, fill=NA))

    #write to file
    ggsave('cache/basinExchange.png', map, width=10.5, height=8.5)
}

#' makeReachMap
#'
#' Makes figure of CONUS model results
#' 
#' @param basinsList_02 list of basin sf mapping objects for 2% flood
#' @param basinsList_10 list of basin sf mapping objects for 10% flood
#' @param basinsList_20 list of basin sf mapping objects for 20% flood
#' @param basinsList_50 list of basin sf mapping objects for 50% flood
#'
#' @return figure written to file
makeReachMap <- function(basinsList_02, basinsList_10, basinsList_20, basinsList_50){
    library(ggplot2)
    theme_set(theme_classic())

    #CONUS boundary
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

    basinsList_02 <- Filter(function(x) dim(x)[1] > 0, basinsList_02)
    basinsList_10 <- Filter(function(x) dim(x)[1] > 0, basinsList_10)
    basinsList_20 <- Filter(function(x) dim(x)[1] > 0, basinsList_20)
    basinsList_50 <- Filter(function(x) dim(x)[1] > 0, basinsList_50)

    #plot
    map_02 <- ggplot() +
            geom_sf(data=basinsList_02[[1]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[2]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[3]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[4]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[5]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[6]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[7]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[8]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[9]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[10]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[11]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[12]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[13]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[14]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[15]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[16]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[17]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[18]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[19]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[20]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[21]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[22]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[23]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[24]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[25]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[26]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[27]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[28]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[29]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[30]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[31]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[32]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[33]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[34]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[35]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[36]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[37]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[38]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[39]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[40]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[41]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[42]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[43]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[44]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[45]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[46]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[47]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[48]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[49]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[50]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[51]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[52]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[53]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[54]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[55]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[56]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[57]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[58]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[59]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[60]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[61]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[62]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[63]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[64]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[65]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[66]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[67]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[68]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[69]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[70]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[71]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[72]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[73]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[74]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[75]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[76]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[77]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[78]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[79]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[80]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[81]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[82]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[83]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[84]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[85]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[86]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[87]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[88]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[89]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[90]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[91]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[92]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[93]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[94]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[95]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[96]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[97]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[98]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[99]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[100]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[101]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[102]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[103]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[104]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[105]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[106]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[107]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[108]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[109]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[110]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[111]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[112]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[113]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[114]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[115]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[116]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[117]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[118]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[119]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[120]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[121]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[122]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[123]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[124]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[125]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[126]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[127]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[128]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[129]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[130]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[131]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[132]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[133]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[134]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[135]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[136]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[137]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[138]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[139]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[140]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[141]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[142]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[143]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[144]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[145]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[146]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[147]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[148]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[149]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[150]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[151]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[152]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[153]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[154]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[155]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[156]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[157]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[158]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[159]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[160]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[161]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[162]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[163]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[164]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[165]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[166]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[167]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[168]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[169]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[170]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[171]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[172]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[173]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[174]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[175]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[176]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[177]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[178]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[179]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[180]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[181]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[182]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[183]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[184]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[185]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[186]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[187]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[188]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[189]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[190]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[191]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[192]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[193]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[194]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[195]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[196]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[197]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[198]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[199]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[201]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[202]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[203]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[204]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_02[[205]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            ggtitle('2% flood')+
            geom_sf(data=states,color='black',size=1.5,alpha=0)+  
            scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
            xlab('')+
            ylab('')+
            theme(axis.text = element_text(size=22))+
            theme(legend.position='bottom')+
            theme(legend.title = element_text(size = 28),
                legend.text = element_text(size = 24),
                plot.title = element_text(hjust = 0.5, size=30),
                panel.border=element_rect(colour="black",size=1, fill=NA))+
            guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 8), title.position = "top", title.hjust = 0.5))

        map_10 <- ggplot() +
            geom_sf(data=basinsList_10[[1]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[2]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[3]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[4]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[5]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[6]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[7]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[8]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[9]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[10]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[11]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[12]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[13]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[14]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[15]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[16]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[17]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[18]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[19]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[20]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[21]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[22]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[23]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[24]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[25]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[26]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[27]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[28]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[29]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[30]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[31]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[32]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[33]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[34]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[35]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[36]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[37]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[38]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[39]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[40]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[41]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[42]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[43]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[44]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[45]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[46]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[47]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[48]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[49]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[50]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[51]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[52]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[53]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[54]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[55]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[56]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[57]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[58]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[59]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[60]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[61]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[62]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[63]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[64]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[65]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[66]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[67]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[68]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[69]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[70]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[71]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[72]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[73]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[74]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[75]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[76]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[77]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[78]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[79]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[80]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[81]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[82]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[83]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[84]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[85]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[86]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[87]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[88]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[89]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[90]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[91]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[92]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[93]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[94]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[95]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[96]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[97]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[98]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[99]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[100]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[101]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[102]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[103]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[104]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[105]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[106]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[107]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[108]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[109]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[110]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[111]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[112]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[113]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[114]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[115]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[116]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[117]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[118]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[119]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[120]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[121]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[122]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[123]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[124]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[125]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[126]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[127]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[128]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[129]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[130]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[131]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[132]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[133]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[134]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[135]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[136]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[137]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[138]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[139]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[140]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[141]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[142]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[143]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[144]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[145]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[146]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[147]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[148]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[149]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[150]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[151]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[152]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[153]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[154]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[155]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[156]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[157]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[158]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[159]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[160]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[161]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[162]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[163]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[164]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[165]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[166]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[167]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[168]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[169]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[170]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[171]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[172]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[173]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[174]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[175]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[176]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[177]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[178]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[179]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[180]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[181]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[182]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[183]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[184]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[185]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[186]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[187]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[188]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[189]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[190]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[191]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[192]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[193]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[194]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[195]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[196]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[197]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[198]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[199]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[201]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[202]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[203]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[204]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_10[[205]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            ggtitle('10% flood')+
            geom_sf(data=states,color='black',size=1.5,alpha=0)+
            scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
            xlab('')+
            ylab('')+
            theme(axis.text = element_text(size=22))+
            theme(legend.position='bottom')+
            theme(legend.title = element_text(size = 28),
                legend.text = element_text(size = 24),
                plot.title = element_text(hjust = 0.5, size=30),
                panel.border=element_rect(colour="black",size=1, fill=NA))+
            guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 8), title.position = "top", title.hjust = 0.5))

        map_50 <- ggplot() +
            geom_sf(data=basinsList_50[[1]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[2]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[3]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[4]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[5]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[6]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[7]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[8]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[9]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[10]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[11]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[12]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[13]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[14]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[15]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[16]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[17]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[18]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[19]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[20]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[21]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[22]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[23]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[24]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[25]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[26]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[27]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[28]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[29]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[30]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[31]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[32]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[33]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[34]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[35]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[36]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[37]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[38]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[39]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[40]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[41]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[42]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[43]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[44]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[45]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[46]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[47]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[48]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[49]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[50]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[51]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[52]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[53]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[54]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[55]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[56]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[57]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[58]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[59]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[60]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[61]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[62]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[63]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[64]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[65]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[66]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[67]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[68]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[69]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[70]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[71]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[72]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[73]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[74]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[75]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[76]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[77]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[78]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[79]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[80]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[81]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[82]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[83]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[84]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[85]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[86]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[87]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[88]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[89]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[90]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[91]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[92]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[93]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[94]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[95]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[96]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[97]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[98]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[99]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[100]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[101]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[102]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[103]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[104]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[105]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[106]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[107]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[108]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[109]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[110]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[111]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[112]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[113]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[114]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[115]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[116]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[117]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[118]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[119]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[120]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[121]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[122]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[123]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[124]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[125]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[126]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[127]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[128]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[129]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[130]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[131]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[132]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[133]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[134]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[135]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[136]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[137]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[138]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[139]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[140]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[141]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[142]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[143]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[144]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[145]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[146]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[147]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[148]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[149]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[150]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[151]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[152]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[153]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[154]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[155]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[156]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[157]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[158]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[159]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[160]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[161]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[162]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[163]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[164]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[165]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[166]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[167]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[168]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[169]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[170]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[171]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[172]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[173]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[174]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[175]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[176]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[177]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[178]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[179]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[180]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[181]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[182]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[183]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[184]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[185]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[186]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[187]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[188]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[189]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[190]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[191]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[192]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[193]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[194]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[195]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[196]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[197]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[198]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[199]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[201]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[202]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[203]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[204]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_50[[205]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            ggtitle('50% flood')+
            geom_sf(data=states,color='black',size=1.5,alpha=0)+
            scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
            xlab('')+
            ylab('')+
            theme(axis.text = element_text(size=22))+
            theme(legend.position='bottom')+
            theme(legend.title = element_text(size = 28),
                legend.text = element_text(size = 24),
                plot.title = element_text(hjust = 0.5, size=30),
                panel.border=element_rect(colour="black",size=1, fill=NA))+
            guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 8), title.position = "top", title.hjust = 0.5))


    map_20 <- ggplot() +
            geom_sf(data=basinsList_20[[1]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[2]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[3]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[4]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[5]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[6]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[7]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[8]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[9]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[10]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[11]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[12]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[13]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[14]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[15]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[16]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[17]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[18]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[19]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[20]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[21]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[22]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[23]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[24]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[25]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[26]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[27]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[28]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[29]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[30]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[31]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[32]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[33]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[34]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[35]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[36]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[37]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[38]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[39]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[40]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[41]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[42]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[43]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[44]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[45]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[46]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[47]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[48]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[49]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[50]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[51]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[52]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[53]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[54]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[55]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[56]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[57]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[58]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[59]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[60]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[61]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[62]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[63]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[64]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[65]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[66]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[67]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[68]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[69]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[70]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[71]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[72]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[73]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[74]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[75]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[76]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[77]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[78]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[79]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[80]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[81]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[82]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[83]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[84]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[85]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[86]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[87]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[88]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[89]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[90]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[91]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[92]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[93]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[94]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[95]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[96]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[97]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[98]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[99]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[100]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[101]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[102]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[103]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[104]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[105]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[106]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[107]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[108]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[109]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[110]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[111]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[112]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[113]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[114]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[115]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[116]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[117]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[118]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[119]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[120]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[121]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[122]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[123]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[124]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[125]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[126]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[127]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[128]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[129]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[130]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[131]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[132]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[133]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[134]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[135]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[136]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[137]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[138]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[139]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[140]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[141]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[142]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[143]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[144]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[145]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[146]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[147]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[148]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[149]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[150]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[151]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[152]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[153]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[154]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[155]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[156]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[157]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[158]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[159]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[160]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[161]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[162]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[163]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[164]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[165]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[166]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[167]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[168]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[169]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[170]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[171]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[172]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[173]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[174]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[175]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[176]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[177]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[178]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[179]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[180]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[181]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[182]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[183]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[184]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[185]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[186]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[187]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[188]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[189]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[190]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[191]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[192]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[193]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[194]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[195]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[196]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[197]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[198]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[199]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[201]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[202]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[203]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[204]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            geom_sf(data=basinsList_20[[205]], aes(color= tau_col), linewidth = 0.25,show.legend = "line") +
            ggtitle('20% flood')+
            geom_sf(data=states,color='black',size=1.5,alpha=0)+  
            scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
            xlab('')+
            ylab('')+
            theme(axis.text = element_text(size=22))+
            theme(legend.position='bottom')+
            theme(legend.title = element_text(size = 28),
                legend.text = element_text(size = 24),
                plot.title = element_text(hjust = 0.5, size=30),
                panel.border=element_rect(colour="black",size=1, fill=NA))+
            guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 8), title.position = "top", title.hjust = 0.5))

    #setup inset box    
    zoom_to <- c(-79.5000, 39.1399)

    #zoom bounds
    zoom_level <- 3
    lon_span <- 360 / 5^zoom_level
    lat_span <- 360 / 5^zoom_level
    lon_bounds_1 <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
    lat_bounds_1 <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

    #inset box
    df <- data.frame(lon_bounds_1, lat_bounds_1)
    box_1 <- df %>% 
        sf::st_as_sf(coords = c("lon_bounds_1", "lat_bounds_1"), crs = 4326) %>% 
        sf::st_bbox() %>% 
        sf::st_as_sfc()

    # add inset 1 bounding box to map
    map_50 <- map_50 +
        geom_sf(data=box_1,
                color='#fca311',
                linewidth=3,
                alpha=0)

    #write to file
    layout <- "
        AB
        CD
        "

    combo_plot <- patchwork::wrap_plots(A=map_02, B=map_10, C=map_20, D=map_50, design=layout, guides='collect') &
        theme(legend.position='bottom',
            legend.title = element_text(size = 28))
    
    ggsave('cache/reachTauMap.png', combo_plot, width=18, height=12)
}

makeReachMapInset <- function(map_0205, map_0206,map_0207,map_0208,map_0502, map_0503, map_0501, map_0505){
    library(ggplot2)
	theme_set(theme_classic())

	#setup inset box   
    #centering
    zoom_to <- c(-79.5000, 39.1399)

	# INSET 1
    #zoom bounds
    zoom_level <- 3
    lon_span <- 360 / 5^zoom_level
    lat_span <- 360 / 5^zoom_level
    lon_bounds_1 <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
    lat_bounds_1 <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

    #set up inset box
    df <- data.frame(lon_bounds_1, lat_bounds_1)
    box_1 <- df %>% 
        sf::st_as_sf(coords = c("lon_bounds_1", "lat_bounds_1"), crs = 4326) %>% 
        sf::st_bbox() %>% 
        sf::st_as_sfc()

	#INSET 2
    #zoom bounds
    zoom_level <- 4
	lon_span <- 360 / 5^zoom_level
	lat_span <- 360 / 5^zoom_level
	lon_bounds_2 <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
	lat_bounds_2 <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

	#inset box
    df <- data.frame(lon_bounds_2, lat_bounds_2)
	box_2 <- df %>% 
        sf::st_as_sf(coords = c("lon_bounds_2", "lat_bounds_2"), crs = 4326) %>% 
        sf::st_bbox() %>% 
        sf::st_as_sfc()

	#INSET 3
    #zoom bounds
    zoom_level <- 5
	lon_span <- 360 / 5^zoom_level
	lat_span <- 360 / 5^zoom_level
	lon_bounds_3 <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
	lat_bounds_3 <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

	#inset box
    df <- data.frame(lon_bounds_3, lat_bounds_3)
	box_3 <- df %>% 
        sf::st_as_sf(coords = c("lon_bounds_3", "lat_bounds_3"), crs = 4326) %>% 
        sf::st_bbox() %>% 
        sf::st_as_sfc()

	#INSET 4
    #zoom bounds
    zoom_level <- 6
	lon_span <- 360 / 5^zoom_level
	lat_span <- 360 / 5^zoom_level
	lon_bounds_4 <- c(zoom_to[1] - lon_span / 2, zoom_to[1] + lon_span / 2)
	lat_bounds_4 <- c(zoom_to[2] - lat_span / 2, zoom_to[2] + lat_span / 2)

	#inset box
    df <- data.frame(lon_bounds_4, lat_bounds_4)
	box_4 <- df %>% 
        sf::st_as_sf(coords = c("lon_bounds_4", "lat_bounds_4"), crs = 4326) %>% 
        sf::st_bbox() %>% 
        sf::st_as_sfc()  		

	#build complete inset river network
    insetNet <- rbind(map_0205, map_0206,map_0207,map_0208,map_0502, map_0503, map_0501, map_0505)

    #plot inset 1
    insetShp1 <- sf::st_crop(insetNet, xmin=lon_bounds_1[1], xmax=lon_bounds_1[2], ymin=lat_bounds_1[1], ymax=lat_bounds_1[2])
    inset1 <- ggplot() +
        geom_sf(data = insetShp1, aes(color = tau_col), linewidth=0.5) +  		
        geom_sf(data=box_2,
            color='#fca311',
            linewidth=3,
            alpha=0) +
        scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
        ggsn::scalebar(data=insetShp1,location='bottomleft', dist = 50, dist_unit = "km",transform = TRUE, model = "WGS84", box.fill=c('white','red'),st.color='white',st.dist=0.05,border.size=0.05) +  
        xlab('')+
        ylab('')+
        coord_sf(datum=NA) +
        theme(legend.position='none',
            panel.background = element_rect(fill = "black"))

    #plot inset 2
    insetShp2 <- sf::st_crop(insetNet, xmin=lon_bounds_2[1], xmax=lon_bounds_2[2], ymin=lat_bounds_2[1], ymax=lat_bounds_2[2])    
    inset2 <- ggplot() +
        geom_sf(data = insetShp2, aes(color = tau_col), linewidth=0.5) + 		
        geom_sf(data=box_3,
            color='#fca311',
            linewidth=3,
            alpha=0) +
        scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
        ggsn::scalebar(data=insetShp2,location='bottomleft', dist = 10, dist_unit = "km",transform = TRUE, model = "WGS84", box.fill=c('white','red'),st.color='white',st.dist=0.05,border.size=0.05) +  
        xlab('')+
        ylab('')+  		
        coord_sf(datum=NA) +
        theme(legend.position='none',
            panel.background = element_rect(fill = "black"))

    #plot inset 3
    insetShp3 <- sf::st_crop(insetNet, xmin=lon_bounds_3[1], xmax=lon_bounds_3[2], ymin=lat_bounds_3[1], ymax=lat_bounds_3[2])
    inset3 <- ggplot() +
    geom_sf(data = insetShp3, aes(color = tau_col), linewidth=0.5) + 		
        scale_color_manual(values=c('#a7bea9', '#84a98c', '#52796f', '#354f52', '#2f3e46'),
                        name=bquote(bold(tau[flood]~'['*hr^1*km^-1*']')),
                        labels=c(bquote(10^-0.5), bquote(10^0), bquote(10^0.5), bquote(10^1), bquote('Over'~10^1)),
                        drop=FALSE)+
        ggsn::scalebar(data=insetShp3,location='bottomleft', dist = 2, dist_unit = "km",transform = TRUE, model = "WGS84", box.fill=c('white','red'),st.color='white',st.dist=0.05,border.size=0.05) +
        xlab('')+
        ylab('')+
        coord_sf(datum=NA) +						
        theme(legend.position='none',
            panel.background = element_rect(fill = "black"))

    # write to file
    design <- "
    CDE
    "

    comboPlot <- patchwork::wrap_plots(C=inset3, D=inset2, E=inset1, design=design)

	ggsave(filename="cache/reachTauMap_inset.png",plot=comboPlot,width=20,height=5)

	return('written to file')
}


#' makeReachBoxplotsFig
#'
#' Makes figure of CONUS model results by stream order
#' 
#' @param combined_basinSummarySO basin stream order summary dataframe
#'
#' @return figure written to file
makeReachBoxplotsFig <- function(combined_basinSummarySO){
    library(ggplot2)
    theme_set(theme_classic())

    #wrangle
    forPlot <- combined_basinSummarySO %>%
        dplyr::select(!median_tau_channel_hr_km)
    colnames(forPlot) <- c('prob', 'StreamCalc', 'median_tau_hr_km', 'huc4')

    channel <- combined_basinSummarySO %>%
        dplyr::filter(prob == 0.02) %>% #for bankfull, idential acros sprobs so just pick one
        dplyr::select(!median_tau_hr_km) %>%
        dplyr::mutate(prob = -999)
    colnames(channel) <- c('prob', 'StreamCalc', 'median_tau_hr_km', 'huc4')

    forPlot <- rbind(forPlot, channel)

    forPlot$prob_lab <- ifelse(forPlot$prob == 0.02, "2% flood",
                                    ifelse(forPlot$prob == 0.1, "10% flood",
                                        ifelse(forPlot$prob == 0.2, "20% flood",
                                            ifelse(forPlot$prob == 0.5, "50% flood",
                                                ifelse(forPlot$prob == -999, "Channel", NA)))))
    forPlot$prob_lab <- factor(forPlot$prob_lab, levels=c('Channel', '2% flood', '10% flood', '20% flood', '50% flood'))

    boxplot_all <- ggplot()+
        geom_boxplot(data=forPlot, aes(x=factor(StreamCalc), y=median_tau_hr_km, fill=prob_lab), color='black') +
        scale_fill_manual(values=c('#ad6a6c', '#faf3dd', '#c8d5b9', '#8fc0a9', '#68b0ab'), name='') +
        scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10), breaks=c(10^-2, 10^-1, 10^0, 10^1))+
        xlab(bquote(bold(Stream~Order)))+
        ylab(bquote(bold(tau~'['*hr^1*km^-1*']'))) +
        theme(legend.position='bottom',
            axis.text=element_text(size=20),
            axis.title = element_text(size=22, face='bold'),
            legend.text=element_text(size=20))
    
    ggsave('cache/streamOrderFig.png', boxplot_all, width=10, height=8)
}

#' makeCalculationValFig
#'
#' Makes figure of gage calculation experiment
#' 
#' @param gageVolume_val_combined hydrodynamic inundation simulation experiment results dataframe
#' @param gageQexc_val Q_flood calculation experiment results dataframe
#'
#' @return figure written to file
makeCalculationValFig <- function(gageVolume_val_combined, gageQexc_val){    
    library(ggplot2)
    theme_set(theme_classic())

    #wrangle
    gageDF_02 <- dplyr::bind_rows(gageQexc_val) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.02)
    
    gageDF_10 <- dplyr::bind_rows(gageQexc_val) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.10)

    gageDF_20 <- dplyr::bind_rows(gageQexc_val) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.20)

    gageDF_50 <- dplyr::bind_rows(gageQexc_val) %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        dplyr::mutate(prob = 0.50)

    gageDF <- rbind(gageDF_02, gageDF_10, gageDF_20, gageDF_50)

    #join obs and model datasets
    df <- gageDF %>%
        dplyr::group_by(site_no, prob, flag)%>%
        dplyr::summarise(Qexc_m3dy = quantile(Qexc_m3dy, 1-prob, na.rm=T)) %>%
        dplyr::distinct() %>%
        dplyr::filter(Qexc_m3dy > 0) %>%
        tidyr::pivot_wider(id_cols=c(site_no, prob), names_from=flag, values_from=c('Qexc_m3dy')) %>%
        tidyr::drop_na() #remove gages with no modeled flows (and a handful of 'gages' with no flow within our time record)

    #plot
    lm <- lm(log10(model)~log10(obs), data=df)
    df$prob <- factor(df$prob, levels=c(0.50, 0.20, 0.10, 0.02))
    plot1 <- ggplot(df, aes(x=obs, y= model)) +
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_point(aes( fill=prob), pch=21, size=10, color='black') +
        geom_smooth(method='lm', se=F, color='black', linewidth=2)+
        annotate("text", x = 10^4, y = 10^9, label = bquote(R^2*'= '*.(round(summary(lm)$r.squared,2))), size=6)+
        annotate("text", x = 10^4, y = 10^8.65, label = bquote(n*'= '*.(nrow(df)/4)*' gages'), size=6)+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(10^3, 10^9), breaks=c(1e3, 1e4, 1e5, 1e6, 1e7,1e8,1e9))+
        scale_y_log10(guide = "axis_logticks",  labels = scales::label_log(base=10), limits=c(10^3, 10^9), breaks=c(1e3, 1e4, 1e5, 1e6, 1e7,1e8,1e9))+
        scale_fill_brewer(palette='Set2', name='', labels=c('50%', '20%', '10%', '2%'))+
        ylab(bquote(bold('Estimated floodplain flux ['~m^3/dy~']')))+
        xlab(bquote(bold('Observed floodplain flux ['~m^3/dy~']')))+
        labs(tag='A')+
        theme(axis.title = element_text(size=20, face='bold'))+
        theme(legend.position=c(0.8,0.2),
            panel.border=element_rect(colour="black",size=1, fill=NA),
            axis.line = element_blank())+
        theme(axis.text = element_text(family = "Futura-Medium", size=18),
            legend.text = element_text(family = "Futura-Medium", size = 16),
            plot.tag = element_text(size=26,
                                    face='bold'))

    #this reach should be removed because of erroneous linking b/w usgs model and our model
    #There happens to be a gage on an adjacent reach, so it gets picked up. But USGS only simulates ~100m of a 12km reach, so complete mismatch. See qgis file for specifics.
    forPlot <- gageVolume_val_combined[gageVolume_val_combined$NHDPlusID != '60002600017895',]

    #plot
    lm <- lm(log10(V_m3)~log10(V_usgs_m3), data=forPlot)
    plot2 <- ggplot(forPlot, aes(x=V_usgs_m3, y=V_m3))+
        geom_abline(linewidth=2, color='darkgrey', linetype='dashed') +
        geom_point(pch=21, size=10, color='black', fill='#669bbc') +
        geom_smooth(color='black', method='lm', se=F, linewidth=2) +
        annotate("text", x = 10^3, y = 10^8, label = bquote(R^2*'= '*.(round(summary(lm)$r.squared,2))), size=6)+
        annotate("text", x = 10^3, y = 10^7.65, label = bquote(n*'= '*.(nrow(forPlot))*' reaches'), size=6)+
        scale_x_log10(guide = "axis_logticks", limits=c(1e2, 1e8), breaks=c(1e2,1e3, 1e4, 1e5, 1e6, 1e7,1e8), labels = scales::label_log(base=10))+
        scale_y_log10(guide = "axis_logticks", limits=c(1e2, 1e8), breaks=c(1e2,1e3, 1e4, 1e5, 1e6, 1e7,1e8), labels = scales::label_log(base=10))+
        xlab(bquote(bold("Observed inundated volume ["*m^3*"]"))) +
        ylab(bquote(bold("Estimated inundated volume ["*m^3*"]"))) +
        labs(tag='B')+
        theme(axis.title = element_text(size=20, face='bold'))+
        theme(legend.position='none',
            panel.border=element_rect(colour="black",size=1, fill=NA),
            axis.line = element_blank())+
        theme(axis.text = element_text(family = "Futura-Medium", size=18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))

    #write to file
    layout <- "
    AB
    "

    combo_plot <- patchwork::wrap_plots(A=plot1, B=plot2, design=layout)

    ggsave('cache/validation_calculation.png', combo_plot, width=16, height=8)

    return(forPlot)
}

#' makeMLValFig
#'
#' Makes figure of CONUS model validation
#' 
#' @param trainedModel_Q trained model for Q_flood
#' @param trainedModel_V trained model for V_flood
#' @param modelDF dataframe of gages for model training
#'
#' @return figure written to file
makeMLValFig <- function(trainedModel_Q, trainedModel_V, modelDF){
    library(tidymodels)
    library(ggplot2)
    theme_set(theme_classic())

    # Q model wrangle
    full_predictions_Q <- data.frame()
    for(i in 1:length(trainedModel_Q)){
        results <- trainedModel_Q[[i]]$summary %>% collect_predictions()
        results$fold <- i
        results <- results %>%
            dplyr::relocate(fold)
        full_predictions_Q <- rbind(full_predictions_Q, results)
    }

    #V model wrangle
    full_predictions_V <- data.frame()
    for(i in 1:length(trainedModel_V)){
        results <- trainedModel_V[[i]]$summary %>% collect_predictions()
        results$fold <- i
        results <- results %>%
            dplyr::relocate(fold)
        full_predictions_V <- rbind(full_predictions_V, results)
    }

    #add flow prob for plotting
    modelDF$.row <- row_number(modelDF)
    modelDF <- modelDF %>%
        dplyr::select(c('.row', 'prob'))

    full_predictions_Q <- full_predictions_Q %>%
        dplyr::left_join(modelDF, by='.row') %>%
        dplyr::filter(is.na(prob)==0)
    
    full_predictions_V <- full_predictions_V %>%
        dplyr::left_join(modelDF, by='.row') %>%
        dplyr::filter(is.na(prob)==0)

    #plot
    full_predictions_Q$prob <- factor(full_predictions_Q$prob, levels=c(0.50, 0.20, 0.10, 0.02))
    r2_Q <- round(summary(lm(.pred~Qexc_m3dy, data=full_predictions_Q))$r.squared,2)
    scatter_Q <- ggplot(full_predictions_Q, aes(x=10^(Qexc_m3dy), y=10^(.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(aes(fill=factor(prob)), pch=21, size=8, color='black')+
        geom_smooth(method='lm', se=F, linewidth=2, color='black', show.legend = FALSE)+
        annotate("text", x = 10^4.5, y = 1e10, label = bquote(bold(n:~.(nrow(full_predictions_Q)))), size=8, color='black')+
        annotate("text", x = 10^4.5, y = 10^9.5, label = bquote(bold(r^2*':'~.(r2_Q))), size=8, color='black')+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e3, 1e10))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e3, 1e10))+
        scale_fill_brewer(palette='Set2', name='Flood Probability', labels=c('50%', '20%', '10%', '2%'))+
        xlab(bquote(bold('Observed '~Q[flood]~'['~m^3*dy^-1~']')))+
        ylab(bquote(bold('Modeled '~Q[flood]~'['~m^3*dy^-1~']')))+
        theme(legend.position='bottom',
            legend.title=element_text(size=26),
            legend.text=element_text(size=22),
            axis.title=element_text(size=24),
            axis.text = element_text(size=22),
            panel.border=element_rect(colour="black",size=1, fill=NA))

    full_predictions_V$prob <- factor(full_predictions_V$prob, levels=c(0.50, 0.20, 0.10, 0.02))
    r2_V <- round(summary(lm(.pred~V_m3, data=full_predictions_V))$r.squared,2)
    scatter_V <- ggplot(full_predictions_V, aes(x=10^(V_m3), y=10^(.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(aes(fill=factor(prob)),pch=21, size=8, color='black')+
        geom_smooth(method='lm', se=F, linewidth=2, color='black', show.legend = FALSE)+
        annotate("text", x = 10^2.5, y = 1e9, label = bquote(bold(n:~.(nrow(full_predictions_V)))), size=8, color='black')+
        annotate("text", x = 10^2.5, y = 10^8.4, label = bquote(bold(r^2*':'~.(r2_V))), size=8, color='black')+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e1, 1e9))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e1, 1e9))+
        scale_fill_brewer(palette='Set2', name='Flood Probability', labels=c('50%', '20%', '10%', '2%'))+
        xlab(bquote(bold('Observed '~V[flood]~'['~m^3~']')))+
        ylab(bquote(bold('Modeled '~V[flood]~'['~m^3~']')))+
        theme(legend.position='bottom',
            legend.title=element_text(size=26),
            legend.text=element_text(size=22),
            axis.title=element_text(size=24),
            axis.text = element_text(size=22),
            panel.border=element_rect(colour="black",size=1, fill=NA))

    #combine for tau
    full_predictions <- full_predictions_Q %>%
        dplyr::select(!prob)%>%
        dplyr::left_join(full_predictions_V, by='.row') %>%
        dplyr::mutate(log10tau_hr = (V_m3 - Qexc_m3dy) + log10(24),
                    log10tau_hr.pred = (.pred.y - .pred.x) + log10(24)) %>%
        tidyr::drop_na()

    r2 <- round(summary(lm(log10tau_hr.pred~log10tau_hr, data=full_predictions))$r.squared,2)
    full_predictions$prob <- factor(full_predictions$prob, levels=c(0.50, 0.20, 0.10, 0.02))
    scatter_tau <- ggplot(full_predictions, aes(x=10^(log10tau_hr), y=10^(log10tau_hr.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(aes(fill=factor(prob)),pch=21, size=8, color='black')+
        geom_smooth(method='lm', se=F, linewidth=2, color='black',show.legend = FALSE)+
        annotate("text", x = 10^-3.5, y = 10^4, label = bquote(bold(n:~.(nrow(full_predictions)))), size=8, color='black')+
        annotate("text", x = 10^-3.5, y = 10^3.3, label = bquote(bold(r^2*':'~.(r2))), size=8, color='black')+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e-5, 1e4))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e-5, 1e4))+
        scale_fill_brewer(palette='Set2', name='Flood Probability', labels=c('50%', '20%', '10%', '2%'))+
        xlab(bquote(bold('Observed'~tau[flood]~'[hr]')))+
        ylab(bquote(bold('Modeled'~tau[flood]~'[hr]')))+
        theme(legend.position='bottom',
            legend.title=element_text(size=26),
            legend.text=element_text(size=22),
            axis.title=element_text(size=24),
            axis.text = element_text(size=22),
            panel.border=element_rect(colour="black",size=1, fill=NA))

    #write to file
    layout <- "
        ABC
    "

    comboPlot <- patchwork::wrap_plots(A=scatter_Q, B=scatter_V, C=scatter_tau, design=layout , guides='collect') &
        theme(legend.position='bottom')

    ggsave('cache/validationML.png', comboPlot, width=20, height=8)

    return(list('r2_Q'=r2_Q,
                'r2_V'=r2_V,
                'r2_tau'=r2))
}

#' makeMLQFig
#'
#' Makes figure of CONUS Q_total model validation
#' 
#' @param trainedModel_Q trained model for Q_total
#' @param modelDF dataframe of gages for model training
#'
#' @return figure written to file
makeMLQFig <- function(trainedModel_Q, modelDF){
    library(tidymodels)
    library(ggplot2)
    theme_set(theme_classic())

    # Q model wrangle
    full_predictions_Q <- data.frame()
    for(i in 1:length(trainedModel_Q)){
        results <- trainedModel_Q[[i]]$summary %>% collect_predictions()
        results$fold <- i
        results <- results %>%
            dplyr::relocate(fold)
        full_predictions_Q <- rbind(full_predictions_Q, results)
    }

    #add flow prob for plotting
    modelDF$.row <- row_number(modelDF)
    modelDF <- modelDF %>%
        dplyr::select(c('.row', 'prob'))

    full_predictions_Q <- full_predictions_Q %>%
        dplyr::left_join(modelDF, by='.row') %>%
        dplyr::filter(is.na(prob)==0)

    #plot
    full_predictions_Q$prob <- factor(full_predictions_Q$prob, levels=c(0.50, 0.20, 0.10, 0.02))
    r2_Q <- round(summary(lm(.pred~Q_m3dy, data=full_predictions_Q))$r.squared,2)
    scatter_Q <- ggplot(full_predictions_Q, aes(x=10^(Q_m3dy), y=10^(.pred))) + #natural space so we can use scale_log10
        geom_abline(linetype='dashed', color='darkgrey', linewidth=1.75) +
        geom_point(aes(fill=factor(prob)), pch=21, size=8, color='black')+
        geom_smooth(method='lm', se=F, linewidth=2, color='black', show.legend = FALSE)+
        annotate("text", x = 10^4.5, y = 1e10, label = bquote(bold(n:~.(nrow(full_predictions_Q)))), size=8, color='black')+
        annotate("text", x = 10^4.5, y = 10^9.5, label = bquote(bold(r^2*':'~.(r2_Q))), size=8, color='black')+
        scale_x_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e4, 1e10))+
        scale_y_log10(guide = "axis_logticks", labels = scales::label_log(base=10), limits=c(1e4, 1e10))+
        scale_fill_brewer(palette='Set2', name='Flood Probability', labels=c('50%', '20%', '10%', '2%'))+
        xlab(bquote(bold('Observed '~Q[total]~'['~m^3*dy^-1~']')))+
        ylab(bquote(bold('Modeled '~Q[total]~'['~m^3*dy^-1~']')))+
        theme(legend.position='bottom',
            legend.title=element_text(size=20),
            legend.text=element_text(size=15),
            axis.title=element_text(face='bold', size=18),
            axis.text = element_text(size=15),
            panel.border=element_rect(colour="black",size=1, fill=NA))

    #write to file
    ggsave('cache/validationTotalQML.png', scatter_Q, height=7, width=7)

    return('cache/validationTotalQML.png')
}

#' makeGageMap
#'
#' Makes figure of CONUS streamgages used in study
#' 
#' @param gage_df dataframe of gages for model training
#'
#' @return figure written to file
makeGageMap <- function(gage_df) {
    sf::sf_use_s2(FALSE)
    library(ggplot2)
    theme_set(theme_classic())

    #filter for the US only
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

    #plot
    map <- ggplot(gage_shp) +
        geom_sf(data=states,
            color='black',
            size=2,
            alpha=0)+
        geom_sf(fill='#3a5a40',
            color='black',
            pch=21,
            size=4) +
        annotate("text", x = -119, y = 25.1, label = bquote(.(nrow(gage_shp)/4)~"gages"), size=6)+
        theme(axis.title = element_text(size=26, face='bold'),axis.text = element_text(family="Futura-Medium", size=20))+
        theme(legend.position='bottom')+
        theme(text = element_text(family = "Futura-Medium"),
            legend.title = element_text(face = "bold", size = 18),
            legend.text = element_text(family = "Futura-Medium", size = 18),
            plot.tag = element_text(size=26,
                                    face='bold'))+
        xlab('')+
        ylab('')

    #write to file
    ggsave('cache/gageMap.png', map, width=10, height=7)

    return(gage_shp)
}

#' makeHuc4Map
#'
#' Makes figure of CONUS basins used in study
#'
#' @return figure written to file
makeHuc4Map <- function(){
    library(sf)
    library(ggplot2)
    theme_set(theme_classic())

    #CONUS boundary
    states <- sf::st_read('data/path_to_data/CONUS_sediment_data/cb_2018_us_state_5m.shp')
    states <- dplyr::filter(states, !(NAME %in% c('Alaska',
                                                'American Samoa',
                                                'Commonwealth of the Northern Mariana Islands',
                                                'Guam',
                                                'District of Columbia',
                                                'Puerto Rico',
                                                'United States Virgin Islands',
                                                'Hawaii'))) #remove non CONUS states/territories

    #regional basin ids
    hucs <- c('0101', '0102', '0103', '0104', '0105', '0106', '0107', '0108', '0109', '0110',
        '0202', '0203', '0204', '0205', '0206', '0207', '0208',
        '0301', '0302', '0303', '0304', '0305', '0306', '0307', '0308', '0309', '0310', '0311', '0312', '0313', '0314', '0315', '0316', '0317', '0318',
        '0401', '0402', '0403', '0404', '0405', '0406', '0407', '0408', '0409', '0410', '0411', '0412', '0413', '0414', '0420', '0427', '0429', '0430', '0431',
        '0501', '0502', '0503', '0504', '0505', '0506', '0507', '0508', '0509', '0510', '0511', '0512', '0513', '0514',
        '0601', '0602', '0603', '0604',
        '0701', '0702','0703', '0704', '0705', '0706', '0707', '0708', '0709', '0710', '0711', '0712', '0713', '0714',
        '0801', '0802', '0803', '0804', '0805', '0806', '0807', '0808', '0809',
        '0901', '0902', '0903', '0904',
        '1002', '1003', '1004', '1005', '1006', '1007', '1008', '1009', '1010', '1011', '1012', '1013', '1014', '1015', '1016', '1017', '1018', '1019', '1020', '1021', '1022', '1023', '1024', '1025', '1026', '1027', '1028', '1029', '1030',
        '1101', '1102', '1103', '1104', '1105', '1106', '1107', '1108', '1109', '1110', '1111', '1112', '1113', '1114',
        '1201', '1202', '1203', '1204', '1205', '1206', '1207', '1208', '1209','1210','1211',
        '1301','1302','1303','1304','1305','1306','1307','1308','1309',
        '1401','1402','1403','1404','1405','1406','1407','1408',
        '1501','1502','1503','1504','1505','1506','1507','1508',
        '1601','1602','1603','1604','1605','1606',
        '1701','1702','1703','1704','1705','1706','1707','1708','1709','1710','1711','1712',
        '1801','1802','1803','1804','1805','1806','1807','1808','1809','1810')

    #wrangle basins
    basins <- sf::st_read('data/HUC4s.shp') %>%
        dplyr::filter(huc4 %in% hucs) %>%
        dplyr::filter(name != 'Lake Erie')

    states <- sf::st_union(states) %>%
        sf::st_transform(crs=sf::st_crs(basins))

    basins <- basins %>%
        sf::st_intersection(states)

    #plot
    map <- ggplot(basins) +
        geom_sf(fill='#a4ac86', color='#333d29', linewidth=0.75) +
        geom_sf(data=states,color='black',linewidth=1.5,alpha=0) +
        theme(axis.text = element_text(size=22))

    #write to file
    ggsave('cache/huc4s.png', map, width=12, height=10)
}