
# ------------------------------------------------------------------------------- #
# Table calculations
# ------------------------------------------------------------------------------- #
library(tidyverse)
library(fst)
library(raster)
df <- fst::read.fst(path = '//catalogue/Workspace14/WFP_ClimateRiskPr/vihiga.fst')
head(df)
ID <- 311463
df <- df[df$id == ID,] # All raster
sl <- fst::read.fst(path = '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/soil/KEN/soilcp_data.fst')
scp <- raster::rasterFromXYZ(xyz = sl[,c('x','y','scp')])
ssat <- raster::rasterFromXYZ(xyz = sl[,c('x','y','ssat')])
soilcp <- raster::extract(x = scp, y = unique(df[,c('x','y')]))
soilst <- raster::extract(x = ssat, y = unique(df[,c('x','y')]))
source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')

watbal_loc <- watbal_wrapper(out_all = df, soilcp = soilcp, soilsat = soilst)
plot(watbal_loc$AVAIL, ty = 'l')

season <- 7:10
cnd <- lubridate::month(df$Date) %in% season
yrs_dts <<- split(df$Date[cnd],cumsum(c(1,diff(df$Date[cnd])!=1)))
yrs_dts <- c(yrs_dts[[1]],yrs_dts[[2]],yrs_dts[[3]])

df2 <- df[df$Date %in% yrs_dts,]
watbal_loc2 <- watbal_wrapper(out_all = df2, soilcp = soilcp, soilsat = soilst)

plot(x = watbal_loc$AVAIL[df$Date %in% yrs_dts], y = watbal_loc2$AVAIL, pch = 20)
abline(0, 1)
plot(x = watbal_loc$ERATIO[df$Date %in% yrs_dts], y = watbal_loc2$ERATIO, pch = 20)
abline(0, 1)
# ------------------------------------------------------------------------------- #
# Raster calculations
# ------------------------------------------------------------------------------- #

# Raster template
tmp <- terra::rast('//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif')
shp <- terra::vect(shp_fl)
tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
tmp[!is.na(tmp)] <- 1

# Loading data
tmn <- terra::rast(tmn_fls[tmn_dts %in% yrs_dts[[i]]])
tmn <- tmn %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
tmn <- tmn - 273.15
tmn <- tmn %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
tav <- terra::rast(tav_fls[tav_dts %in% yrs_dts[[i]]])
tav <- tav %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
tav <- tav - 273.15
tav <- tav %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
tmx <- tmx %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
tmx <- tmx - 273.15
tmx <- tmx %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
srd <- terra::rast(srd_fls[srd_dts %in% yrs_dts[[i]]])
srd <- srd %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
srd <- srd/1000000
srd <- srd %>% terra::resample(x = ., y = tmp) %>% terra::mask(tmp)
prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
prc <- prc %>% terra::crop(terra::ext(tmp)) %>% terra::mask(tmp)
prc[prc == -9999] <- 0

# Maximum evapotranspiration
ETMAX <- terra::lapp(x = terra::sds(srd,tmn,tav,tmx), fun = peest)

# Soil data
scp <- terra::rast('D:/india/soil/soilcp.tif')
sst <- terra::rast('D:/india/soil/soilsat.tif')
scp <- scp %>% terra::resample(tmp) %>% terra::mask(tmp) # Soil water capacity
sst <- sst %>% terra::resample(tmp) %>% terra::mask(tmp) # Soil water saturation point

# Compute water balance model
AVAIL <- tmp
AVAIL[!is.na(AVAIL)] <- 0
watbal <- 1:terra::nlyr(ETMAX) %>%
  purrr::map(.f = function(i){
    water_balance <- eabyep_calc(soilcp  = scp,
                                 soilsat = ssat,
                                 avail   = AVAIL[[terra::nlyr(AVAIL)]],
                                 rain    = prc[[i]],
                                 evap    = ETMAX[[i]])
    AVAIL <<- water_balance$Availability
    return(water_balance)
  })
ERATIO  <- watbal %>% purrr::map('Eratio') %>% terra::rast()
LOGGING <- watbal %>% purrr::map('Logging') %>% terra::rast()
IRR     <- ETMAX - prc
GDAY    <- terra::lapp(x = terra::sds(tav, ERATIO), fun = function(TAV, ERATIO){ifelse(TAV >= 6 & ERATIO >= 0.35, 1, 0)})

NDWS    <- terra::app(x = ERATIO, fun = function(ERATIO){ifelse(ERATIO < 0.5, 1, 0)}) %>% sum()
NWLD    <- terra::app(x = LOGGING, fun = function(LOGGING){ifelse(LOGGING > 0, 1, 0)}) %>% sum()
