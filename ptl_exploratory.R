#PTL - Exploratory

require(pROC)
require(dlnm)
require(SpatialEpi)
require(grid)
require(tsModel)
require(INLA)
library(RColorBrewer)
library(quantmod)
require(maps)
require(caret)
require(leafem)
require(htmlwidgets)



#Figure 1 (Left) Visualisation of study area  -----
piura_tumbes_lambayeque <- c("Piura", "Tumbes", "Lambayeque")
ptl_cols <- hex
peru_cities <- subset(world.cities, country.etc == "Peru")
ptl_cities <- subset(peru_cities, name %in% piura_tumbes_lambayeque)
capitals_geo <- st_as_sf(ptl_cities, coords = c("long", "lat"),
                         crs= "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
capitals_geo

#Map --
mapviewOptions(legend.pos =  "topright")
ptl_map <- subset(peru_adm_boundaries, shapeName %in% piura_tumbes_lambayeque)
ptl_map$rgns <- as.factor(ptl_map$shapeName)
tmp_map <- mapview(peru_adm_boundaries,alpha.regions = 0.1, map.type = c("CartoDB.Voyager"),layer.name = "Regions", legend = FALSE,color = "black",
                   col.regions = "orange")+  mapview(capitals_geo,legend = FALSE,col.regions = "white",cex = 2.75)+
  mapview(ptl_map, layer.name = "Department",  zcol = "rgns" ,col.regions = ptl_cols,legend = TRUE,
          alpha.regions = 0.7)
tmp_map
mapshot(tmp_map,  file = file.path(peru.exploratory.output,"tmp_map_coloured.pdf"))
mapshot(tmp_map,  file = file.path(peru.exploratory.output,"tmp_map_coloured.html"))



#Data.base for Piura, Tumbes, and Lambayeque
ptl_dt <- copy(monthly_peru_cases)
ptl_dt <- subset(ptl_dt, select = c("REGION", "TIME","MONTH", "YEAR", "ym_cases", 
                                                      "DIR","POP", "NATURAL_REGION", "POP_DENSITY",
                                                      "POP_INTERPOLATED", "DIR_POP_INTERP"))
ptl_dt[, POP_OFFSET_INTERP:= POP_INTERPOLATED/100000]
# ptl_dt <- merge(ptl_dt, unique(subset(areas_dt, select = c("REGION", "REGION_AREA_KM2"))), by = "REGION")
ptl_dt[, POP_DENSITY_INTERP:= POP_INTERPOLATED/REGION_AREA_KM2]
ptl_dt[, POP_OFFSET:= POP/100000]
setnames(ptl_dt, "ym_cases", "CASES")
#Remove first four months as need to use climate lags and RSI (based on 3 months of DIR)
ptl_dt <- ptl_dt[which(TIME > 4)]
#Start counting of months at 1 in May 2010
ptl_dt[, TIME:= TIME - 4]



setkeyv(ptl_dt, c( "YEAR", "MONTH", "REGION"))

#Convert to Numeric class
ptl_dt[, CASES:= as.numeric(CASES)]
ptl_dt[, DIR:= as.numeric(DIR)]
ptl_dt[, POP_OFFSET:= as.numeric(POP_OFFSET)]
ptl_dt[, MONTH:= as.numeric(MONTH)]
ptl_dt[, RGN_IND:= as.numeric(RGN_IND)]

#Year indexing starting from 1 in 2010
ptl_dt[, YEAR:= as.numeric(YEAR)-2009]

#Setting up binary seasonality indicator
summer_months <- c(12, seq(1, 4))
ptl_dt[, SEASON:= if_else(MONTH %in% summer_months, 1, 0)]
ptl_dt <- subset(ptl_dt, REGION %in% piura_tumbes_lambayeque)


#RSI - Relative Strength Index
tmp <- subset(monthly_peru_cases , REGION %in% piura_tumbes_lambayeque)
setkeyv(tmp, c("REGION", "TIME"))
tmp2 <- tmp[, list(RSI_DIR = RSI(DIR, n = 3)), by = c("REGION")]
tmp[, RSI_DIR:= tmp2$RSI_DIR]
tmp <- subset(tmp, select = c("RSI_DIR", "REGION", "TIME"))
tmp[, TIME:= TIME - 4]
ptl_dt <- merge(ptl_dt, tmp, by = c("REGION", "TIME"))
tmp[, RSI_DIR_LAG:= Lag(RSI_DIR, 1)]
tmp[, RSI_DIR_LAG_2:= Lag(RSI_DIR, 2)]
tmp
tmp_lag <- subset(tmp, select = c("RSI_DIR_LAG", "REGION", "TIME"))
tmp_lag2 <- subset(tmp, select = c("RSI_DIR_LAG_2", "REGION", "TIME"))

ptl_dt <- merge(ptl_dt, tmp_lag, by = c("REGION", "TIME"))
ptl_dt <- merge(ptl_dt, tmp_lag2, by = c("REGION", "TIME"))

#RSI DIR pop interpolated----
tmp <- subset(monthly_peru_cases , REGION %in% piura_tumbes_lambayeque)
tmp <- na.omit(tmp)
setkeyv(tmp, c("REGION", "TIME"))
tmp2 <- tmp[, list(RSI_DIR_POP_INTERP = RSI(DIR_POP_INTERP, n = 3)), by = c("REGION")]
tmp[, RSI_DIR_POP_INTERP:= tmp2$RSI_DIR_POP_INTERP]
tmp <- subset(tmp, select = c("RSI_DIR_POP_INTERP", "REGION", "TIME"))
tmp[, TIME:= TIME - 4]
ptl_dt <- merge(ptl_dt, tmp, by = c("REGION", "TIME"))
tmp[, RSI_DIR_POP_INTERP_LAG:= Lag(RSI_DIR_POP_INTERP, 1)]
tmp[, RSI_DIR_POP_INTERP_LAG_2:= Lag(RSI_DIR_POP_INTERP, 2)]
tmp_lag <- subset(tmp, select = c("RSI_DIR_POP_INTERP_LAG", "REGION", "TIME"))
tmp_lag2 <- subset(tmp, select = c("RSI_DIR_POP_INTERP_LAG_2", "REGION", "TIME"))
ptl_dt <- merge(ptl_dt, tmp_lag, by = c("REGION", "TIME"))
ptl_dt <- merge(ptl_dt, tmp_lag2, by = c("REGION", "TIME"))

setkeyv(ptl_dt, c("TIME", "REGION"))
tmp2 <- unique(ptl_dt$REGION)
tmp <- data.table(REGION = tmp2, RGN_IND = c(1:length(tmp2)))
ptl_dt <- merge(ptl_dt, tmp, by = "REGION")
setkeyv(ptl_dt, c("TIME", "REGION"))
ptl_dt[,RGN_IND:= as.numeric(RGN_IND)]


#Exclude 2022 Cases as climate data not yet available
ptl_dt <- ptl_dt[which(YEAR < 13), ]


#Climate data.table specifically ----
#Need separate dt for DLNM specification (aside from dt with case data)
ptl_climate_dt <- merge(tmax_2010_2019, tmin_2010_2019, by = c("REGION", "YEAR", "TIME", "MONTH"))
ptl_climate_dt <- merge(ptl_climate_dt, prec_2010_2019, by = c("REGION", "YEAR", "TIME", "MONTH"))
ptl_climate_dt <- merge(ptl_climate_dt, spi_dt, by = c("REGION", "YEAR","MONTH"))
ptl_climate_dt <- merge(ptl_climate_dt, oni_data, by = c("YEAR","MONTH"))
ptl_climate_dt <- merge(ptl_climate_dt, icen_data, by = c("YEAR","MONTH"))

ptl_climate_dt <- ptl_climate_dt[which(REGION %in% piura_tumbes_lambayeque)]
ptl_climate_dt[, tmax_prec:= tmax * prec]

ptl_climate_20_21 <- merge(tmax_2020_2021, prec_2020_2021, by = c("REGION", "YEAR", "TIME", "MONTH"))
ptl_climate_20_21 <- merge(ptl_climate_20_21, tmin_2020_2021, by = c("REGION", "YEAR", "TIME", "MONTH"))
ptl_climate_20_21 <- merge(ptl_climate_20_21, icen_data, by = c("YEAR",  "MONTH"))
ptl_climate_20_21 <- merge(ptl_climate_20_21, oni_data, by = c("YEAR",  "MONTH"))
ptl_climate_20_21 <- merge(ptl_climate_20_21, spi_dt, by = c("REGION", "YEAR",  "MONTH"))
ptl_climate_20_21[, tmax_prec:= tmax * prec]
setkeyv(ptl_climate_20_21, c("TIME",  "REGION"))
ptl_climate_20_21[, TIME:= TIME + 120]
ptl_climate_20_21 <- subset(ptl_climate_20_21, REGION %in% piura_tumbes_lambayeque)
setkeyv(ptl_climate_20_21, c("TIME", "REGION"))
ptl_climate_dt <- rbind(ptl_climate_dt, ptl_climate_20_21)
