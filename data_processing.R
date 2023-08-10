#Script to process inputs and perform exploratory analysis


#Load in Required Packages
require(ncf)
require(scales)
require(readxl)
require(ISOweek)
require(stringr)
require(geodata)
require(weathermetrics)
require(reshape2)
require(ncdf4)
require(here)
require(raster)
require(ggplot2)
require(dplyr)
require(viridis)
require(bfast)
require(rgeoboundaries)
require(sf)
require(data.table)
require(zoo)
require(mapview)
require(spdep)
require(ggpubr)
require(tidyverse)


# 
# peru.data.in.dir <- "C:/Users/mills/Documents/peru_dengue"
# peru.inla.data.in.dir <- "C:/Users/mills/Documents/peru_dengue/INLA/Input"
# peru.spi.data.in.dir <- "C:/Users/mills/Documents/peru_dengue/spi_data"
# peru.inla.data.out.dir <- "C:/Users/mills/Documents/peru_dengue/INLA/Output"
# peru.case_data.in.dir <- "C:/Users/mills/Documents/peru_dengue/case_data"
# peru.exploratory.output <- "C:/Users/mills/Documents/peru_dengue/exploratory/Output"


peru <- rgdal::readOGR("LandCoverData/peru.shp")


#Neighbourhood Matirix---

shp <- st_read("per_admbnda_adm1_ign_20200714.shp")
neighbour_list <- poly2nb(st_make_valid(shp))
neighbour_list
nb2INLA(file.path(peru.inla.data.in.dir, "nbr_mat.graph"),
        neighbour_list)


#Administrative Boundaries ----
peru_adm_boundaries <- geoboundaries("Peru", "adm1")


#Department Area Data ----
tmp2 <- st_as_sf(peru_adm_boundaries)
tmp2 <- as_Spatial(tmp2)
areas_dt <- data.table(REGION = tmp2$shapeName, REGION_AREA_KM2 = area(tmp2)/ 1000000)


#Population Data----
peru_pops <- data.table(read_excel(file.path(peru.data.in.dir, "peru_pop.xlsx")))
peru_pops <- melt(peru_pops, id.var = "REGION",
                  variable.name = "YEAR", value.name = "POP")
peru_pops[, YEAR:= as.numeric(as.character(YEAR))]
peru_pops[, POP:= as.numeric(POP)]
peru_pops[, REGION:= str_to_title(REGION)]

peru_pops[, MONTH:= rep(7, nrow(peru_pops))]

fill_pops_dt <- function(data){
  filled_pops <- data.table(YEAR = rep(unique(data$YEAR, each = 12 * length(unique(data$REGION)))),
                            MONTH = rep(seq(1, 12), length(unique(data$YEAR))*length(unique(data$REGION))))
  return(filled_pops)
}
peru_pops_filled <- fill_pops_dt(peru_pops)
setkeyv(peru_pops_filled, c("YEAR", "MONTH"))
peru_pops_filled[, REGION:= rep(unique(peru_pops$REGION), nrow(peru_pops_filled)/length(unique(peru_pops$REGION)))]
peru_pops_filled <- merge(peru_pops_filled, peru_pops, by = c("YEAR", "MONTH", "REGION"), all.x = TRUE)
peru_pops_filled[, TIME:= rep(seq(1, 276), each = length(unique(REGION)))]
peru_pops_filled <- peru_pops_filled[which(TIME>= 7)]
peru_pops_filled
for(i in 1:length(unique(peru_pops_filled$REGION))){
  rgn_in_q <- unique(peru_pops_filled$REGION)[i]
  tmp <- subset(peru_pops_filled, REGION == rgn_in_q) 
  tmp2 <- na.approx(tmp$POP, xout = 1:270)
  peru_pops_filled[which(REGION == rgn_in_q & TIME <= 271), POP:= tmp2]
}
peru_pops_filled

setnames(peru_pops_filled, "POP", "POP_INTERPOLATED")
peru_pops_filled[, TIME:= NULL]

#Monthly Reported Incidence ----
raw_peru_cases <- data.table(read.csv(file.path(peru.case_data.in.dir, "2010_2019_cases_full_data.csv")))
peru_cases <- raw_peru_cases[, list(TOTAL_CASES = length(Provincia)), by = c("Departamnento", "Semana", "Ano")]
setnames(peru_cases, colnames(peru_cases), c("REGION", "WEEK", "YEAR", "TOTAL_CASES"))
peru_cases[, REGION:= str_to_title(REGION)]

missing_regions_from_cases <- unique(peru_cases$REGION)[which(!(unique(peru_cases$REGION) %in% spatial_adm_boundaries$shapeName))]
peru_cases[which(REGION == "Madre De Dios"), REGION:= "Madre de Dios"]
peru_cases[which(REGION == "San Martin"), REGION:= "San Martín"]
peru_cases[which(REGION == "Junin"), REGION:= "Junín"]
peru_cases[which(REGION == "Huanuco"), REGION:= "Huánuco"]


peru_pops[which(REGION == "Áncash"), REGION:= "Ancash"]
peru_pops[which(REGION == "Madre De Dios"), REGION:= "Madre de Dios"]
peru_pops[which(REGION == "Madre De Dios"), REGION:= "Madre de Dios"]
merge(peru_cases, peru_pops, by = c("REGION", "YEAR"))
subset(peru_cases, YEAR == 2023) #360 observations in 2023 without pop estimate
peru_cases <- merge(peru_cases, peru_pops, by = c("REGION", "YEAR"))
peru_cases[, YEAR:= as.numeric(as.character(YEAR))]
peru_cases[, POP:= as.numeric(POP)]
peru_cases[, TOTAL_CASES:= as.numeric(TOTAL_CASES)]

peru_cases <- merge(peru_cases, peru_pops, by = c("REGION", "YEAR"))
peru_cases[, YEAR:= as.numeric(as.character(YEAR))]
peru_cases[, POP:= as.numeric(POP)]
peru_cases[, TOTAL_CASES:= as.numeric(TOTAL_CASES)]

peru_cases[, REGION:= str_to_title(REGION)]

setkeyv(peru_cases, c("REGION", "YEAR", "WEEK"))

#Fill out surveillance data to include regions with zero cases per week
fill_out_dt <- function(raw_table){
  incomplete_table <- copy(raw_table)
  keycol <- c("YEAR", "WEEK")
  setkeyv(incomplete_table, keycol)
  num_weeks <- length(unique(incomplete_table$WEEK))
  num_years <- length(unique(incomplete_table$YEAR)) 
  complete_table <- data.table(REGION = rep(unique(incomplete_table$REGION), each = num_weeks * num_years),
                               YEAR = rep(unique(incomplete_table$YEAR), num_weeks * length(unique(incomplete_table$REGION))),
                               WEEK = rep(rep(unique(incomplete_table$WEEK), num_years), length(unique(incomplete_table$REGION))))
  complete_table[, TOTAL_CASES:= rep(0, nrow(complete_table))]
  setkeyv(complete_table, c("REGION", "YEAR", "WEEK"))
  for(i in 1:length(unique(complete_table$REGION))){
    region_in_q <- unique(complete_table$REGION)[i]
    for(j in 1:num_years){
      year_in_q <- unique(complete_table$YEAR)[j]
      for(k in 1:num_weeks){
        week_in_q <- unique(complete_table$WEEK)[k]
        tmp <- subset(incomplete_table, REGION== region_in_q & YEAR == year_in_q & WEEK == week_in_q)
        if(nrow(tmp) > 0){
          complete_table[which(REGION== region_in_q & YEAR == year_in_q & WEEK == week_in_q),
                         TOTAL_CASES:= tmp$TOTAL_CASES]
        }
      }
    }
  }
  return(complete_table)
}


complete_peru_cases <- fill_out_dt(peru_cases)
complete_peru_cases

#Ensuring correct region-week case totals
complete_peru_cases <- complete_peru_cases[,list(TOTAL_CASES = sum(TOTAL_CASES)), by = c("REGION", "YEAR", "WEEK")]
complete_peru_cases[, YEAR_WEEK:= paste0(YEAR, "-W", WEEK, "-1")]
complete_peru_cases[which(WEEK < 10), YEAR_WEEK:= paste0(YEAR, "-W0", WEEK, "-1")]
complete_peru_cases[, START_DATE:= ISOweek2date(complete_peru_cases$YEAR_WEEK)]

#Setting up MONTHLY REGIONAL ANALYSIS----
complete_peru_cases[, MONTH:= substr(START_DATE, 1, 7)]
monthly_peru_cases <- complete_peru_cases %>%
  group_by(REGION, MONTH) %>%
  summarize(ym_cases = sum(TOTAL_CASES))
monthly_peru_cases <- as.data.table(monthly_peru_cases)
monthly_peru_cases[, m:=  substr(MONTH, nchar(MONTH) - 1, nchar(MONTH))]
monthly_peru_cases[, YEAR:= substr(MONTH, 1, 4)]
monthly_peru_cases[, m:= as.numeric(m)]
monthly_peru_cases[, MONTH:= NULL]
setnames(monthly_peru_cases, "m", "MONTH")
monthly_peru_cases[, YEAR:= as.numeric(YEAR)]
monthly_peru_cases <- merge(monthly_peru_cases, areas_dt, by = c("REGION"))
monthly_peru_cases[, POP_DENSITY:= POP/REGION_AREA_KM2]

monthly_peru_cases[which(REGION == "Madre De Dios"), REGION:= "Madre de Dios"]
monthly_peru_cases <- subset(monthly_peru_cases, YEAR != 2023)
setkeyv(monthly_peru_cases, c("REGION", "YEAR", "MONTH"))
monthly_peru_cases
# 156/12 = 13 years
monthly_peru_cases[, TIME:= rep(seq(1, 156),length(unique(monthly_peru_cases$REGION)))]

#Check every region-year has 12 monthly totals
table(monthly_peru_cases[, length(MONTH), by = c("REGION", "YEAR")]$V1)

#Merge in population data (both interpolated and yearly)
monthly_peru_cases <- merge(monthly_peru_cases, peru_pops, by = c("REGION", "YEAR"))
monthly_peru_cases <- merge(monthly_peru_cases, peru_pops_filled, by = c("MONTH", "YEAR", "REGION"))

monthly_peru_cases[, DIR:= ym_cases/POP * 100000]
monthly_peru_cases[, DIR_POP_INTERP:= ym_cases/POP_INTERPOLATED * 100000]



#SPI 6 Data --- 
tmp2 <- st_as_sf(peru_adm_boundaries) 
tmp2 <- as_Spatial(tmp2)
spi_file_name_func <- function(year){
  spi_file_name <- file.path(peru.spi.data.in.dir, paste0("spg06_m_wld_",year,"0101_",year,"1201_m.nc"))
  return(spi_file_name)
}
# tmp <- copy(peru_adm_boundaries)
tmp2[, ID:=NULL]
setnames(tmp2, colnames(tmp2)[1], "REGION" )

#Set up data.table for SPI-6

# spi_val_dt <- subset(tmp2, select = c("REGION",colnames(tmp2)[2]))
# spi_dt <- data.table(YEAR = rep(seq(2001, 2022, by = 1), each = length(peru_adm_boundaries$shapeName)*12),
#                      MONTH = rep( seq(1, 12), length(seq(2001, 2022, by = 1))*length(unique(peru_adm_boundaries$shapeName))))
# spi_dt[, REGION:= rep(peru_adm_boundaries$shapeName, each = nrow(spi_dt)/length(peru_adm_boundaries$shapeName))]
# 
# 

spi_data_table_func <- function()
{
  spi_dt <- data.table(YEAR = rep(seq(2001, 2022, by = 1), each = length(peru_adm_boundaries$shapeName)*12),
                       MONTH = rep( seq(1, 12), length(seq(2001, 2022, by = 1))*length(unique(peru_adm_boundaries$shapeName))))
  spi_dt[, REGION:= rep(rep(peru_adm_boundaries$shapeName, each = 12), length(seq(2001, 2022, by = 1)))]
  for(i in 1:length(unique(spi_dt$YEAR))){
    year <- unique(spi_dt$YEAR)[i]
    spi_file_name<- spi_file_name_func(year)
    r <- rast(spi_file_name_func(year))
    r <- crop(r, extent(peru_adm_boundaries))
    tmp <- copy(peru_adm_boundaries)
    tmp$spi <- data.table(raster::extract(r, tmp, fun="mean", weights = TRUE, na.rm = TRUE))
    tmp2 <- data.table(cbind(tmp$shapeName, tmp$spi))
    tmp2[, ID:= NULL]
    setnames(tmp2, colnames(tmp2)[1],c("REGION"))
    for(j in 1:12){
      spi_val_dt <- subset(tmp2, select = c("REGION",colnames(tmp2)[j+1]))
      for(k in 1:length(unique(tmp$shapeName))){
        tmp_region <- unique(tmp$shapeName)[k]
        spi_val_dt2 <- subset(spi_val_dt, REGION == tmp_region)
        spi_val <- subset(spi_val_dt2, select = colnames(spi_val_dt2)[2])
        spi_dt[which(REGION == tmp_region & MONTH == j & YEAR == year), SPI_6:= spi_val]
      }
    }
  }
  return(spi_dt)
}


spi_dt <- spi_data_table_func()
spi_dt <- data.table(spi_dt)


spi_dt[, MONTH:= as.integer(MONTH)]
spi_dt[, YEAR:= as.numeric(YEAR)]

#ENSO Indices ----
#ICEN--
icen_data <- data.table(read.table(file.path(peru.data.in.dir, "icen.txt"), 
                                   header = TRUE))
setnames(icen_data, colnames(icen_data), toupper(colnames(icen_data)))
monthly_peru_cases <- merge(monthly_peru_cases, icen_data, by = c("MONTH", "YEAR"))


#prec
prec_data <- data.table(read.table("prec.ascii.txt", header = TRUE))
prec_data[, SEAS:= as.factor(SEAS)]
prec_data <- subset(prec_data, YR>= 1999)
prec_data[, YR:= as.numeric(YR)]
prec_data[, TOTAL:= as.numeric(TOTAL)]
prec_data[, ANOM:= as.numeric(ANOM)]
setnames(prec_data, colnames(prec_data)[1:2], c("MONTH", "YEAR"))
months_chars <- (unique(prec_data$SEAS))
prec_data[, MONTH:= as.numeric(factor(prec_data$MONTH, levels = months_chars))]
prec_data[, YEAR:= as.numeric(YEAR)]
prec_data[, MONTH:= as.numeric(MONTH)]
monthly_peru_cases <- merge(monthly_peru_cases, prec_data, by = c("MONTH", "YEAR"))


#Worldclim Data ----
#Function to list worldclim 2010-2019 files for climatic variables
list_worldclim_variable_tifs <- function(years_sequence, climate_variable){
  climate_file_list <- 0
  for(i in 1:length(years_sequence)){
    file_year <- years_sequence[i]
    for(j in seq(1,12)){
      if(j < 10)
      {
        j<- paste0("0",j)
      }
      file_month <- j
      if(climate_variable == "tmax"){
        climate_file <- unname(file.path(peru.data.in.dir, "tmax_2010-2019",
                                         paste0("wc2.1_2.5m_",climate_variable,"_",
                                                file_year,"-", file_month,".tif")))
      }
      if(climate_variable == "tmin"){
        climate_file <- unname(file.path(peru.data.in.dir, "tmin_2010-2019",
                                         paste0("wc2.1_2.5m_",climate_variable,"_",
                                                file_year,"-", file_month,".tif")))
      }
      if(climate_variable == "prec"){
        climate_file <- unname(file.path(peru.data.in.dir, "prec_2010-2019",
                                         paste0("wc2.1_2.5m_",climate_variable,"_",
                                                file_year,"-", file_month,".tif")))
      }
      print(climate_file)
      climate_file_list <- c(climate_file_list, climate_file)
    }
  }
  return(climate_file_list[2:length(climate_file_list)])
}

#Function to list worldclim 2020-2021 files for climatic variables
list_worldclim_variable_tifs_20_21 <- function(years_sequence, climate_variable){
  climate_file_list <- 0
  for(i in 1:length(years_sequence)){
    file_year <- years_sequence[i]
    for(j in seq(1,12)){
      if(j < 10)
      {
        j<- paste0("0",j)
      }
      file_month <- j
      if(climate_variable == "tmax"){
        climate_file <- unname(file.path(peru.data.in.dir, "tmax_2020-2021",
                                         paste0("wc2.1_2.5m_",climate_variable,"_",
                                                file_year,"-", file_month,".tif")))
      }
      if(climate_variable == "tmin"){
        climate_file <- unname(file.path(peru.data.in.dir, "tmin_2020-2021",
                                         paste0("wc2.1_2.5m_",climate_variable,"_",
                                                file_year,"-", file_month,".tif")))
      }
      if(climate_variable == "prec"){
        climate_file <- unname(file.path(peru.data.in.dir, "prec_2020-2021",
                                         paste0("wc2.1_2.5m_",climate_variable,"_",
                                                file_year,"-", file_month,".tif")))
      }
      print(climate_file)
      climate_file_list <- c(climate_file_list, climate_file)
    }
  }
  return(climate_file_list[2:length(climate_file_list)])
}

#Function to extract 2010-2019 data
extract_worldclim_variable_regions <- function(years_sequence, climate_variable){
  worldclim_variable_dt <- data.table(YEAR = rep(years_sequence, each = 12*length(peru_adm_boundaries$shapeName)),
                                      REGION = rep(peru_adm_boundaries$shapeName, length(years_sequence)))
  worldclim_variable_dt[, MONTH:= rep(rep(seq(1, 12), each = length(peru_adm_boundaries$shapeName)), length(unique(years_sequence)))]
  worldclim_variable_dt[, TIME:= rep(seq(1,12*length(unique(YEAR))), each = length(unique(REGION)))]
  worldclim_variable_dt[, VALUE:= rep(0, nrow(worldclim_variable_dt))]
  for(i in 1:length(years_sequence)){
    year <- years_sequence[i]
    climate_variable_monthly_files <- list_worldclim_variable_tifs(years_sequence[i], climate_variable)
    for(j in 1:12){
      tmp_climate_var_peru <- mask(raster(climate_variable_monthly_files[j]), as_Spatial(peru_sf))
      tmp_climate_regions <- (raster::extract(tmp_climate_var_peru, peru_adm_boundaries, fun = "mean", weights = TRUE, na.rm= TRUE))
      worldclim_variable_dt[which(YEAR == year & MONTH==j), VALUE:= tmp_climate_regions]
    }
  }
  setnames(worldclim_variable_dt, "VALUE", climate_variable)
  return(worldclim_variable_dt)
}
seq_2010_2019 <- seq(2010, 2019, by = 1)
tmax_2010_2019 <- extract_worldclim_variable_regions(seq_2010_2019, "tmax")
tmin_2010_2019 <- extract_worldclim_variable_regions(seq_2010_2019, "tmin")
prec_2010_2019 <- extract_worldclim_variable_regions(seq_2010_2019, "prec")

#Function to extract 2020-2021 data
extract_worldclim_variable_regions_20_21 <- function(years_sequence, climate_variable){
  worldclim_variable_dt <- data.table(YEAR = rep(years_sequence, each = 12*length(peru_adm_boundaries$shapeName)),
                                      REGION = rep(peru_adm_boundaries$shapeName, length(years_sequence)))
  worldclim_variable_dt[, MONTH:= rep(rep(seq(1, 12), each = length(peru_adm_boundaries$shapeName)), length(unique(years_sequence)))]
  worldclim_variable_dt[, TIME:= rep(seq(1,12*length(unique(YEAR))), each = length(unique(REGION)))]
  worldclim_variable_dt[, VALUE:= rep(0, nrow(worldclim_variable_dt))]
  for(i in 1:length(years_sequence)){
    year <- years_sequence[i]
    climate_variable_monthly_files <- list_worldclim_variable_tifs_20_21(years_sequence[i], climate_variable)
    print(climate_variable_monthly_files)
    for(j in 1:12){
      tmp_climate_var_peru <- mask(raster(climate_variable_monthly_files[j]), as_Spatial(peru_sf))
      tmp_climate_regions <- (raster::extract(tmp_climate_var_peru, peru_adm_boundaries, fun = "mean", weights = TRUE, na.rm= TRUE))
      worldclim_variable_dt[which(YEAR == year & MONTH==j), VALUE:= tmp_climate_regions]
    }
  }
  setnames(worldclim_variable_dt, "VALUE", climate_variable)
  return(worldclim_variable_dt)
}


seq_2020_2021 <- seq(2020, 2021, by = 1)
tmax_2020_2021 <- extract_worldclim_variable_regions_20_21(seq_2020_2021, "tmax")
tmin_2020_2021 <- extract_worldclim_variable_regions_20_21(seq_2020_2021, "tmin")
prec_2020_2021 <- extract_worldclim_variable_regions_20_21(seq_2020_2021, "prec")
