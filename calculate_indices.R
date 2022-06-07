# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
.rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, gtools, sf, furrr, future))

source('https://raw.githubusercontent.com/CIAT-DAPA/agro-clim-indices/main/_main_functions.R')
calc_SHIMP <- function(TMAX, RH){SHI <- ifelse(TMAX >= 29 & RH > 50, 1, 0); return(SHI)}; calc_SHIMP <- compiler::cmpfun(calc_SHIMP)

seasons <- list(Kharif = 7:10, Rabi = c(10:12,1:2), Zaid = 3:6)
shp_fl <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/shps/india/ind_regions/ind_regions_d.shp'

# Function to compute Agro-climatic indices
calc_AgrClm <- function(season = season, shp_fl = shp_fl){
  
  ## ROI: regions of interest
  shp <- terra::vect(shp_fl)
  
  ## Daily files
  # Precipitation
  chr_pth <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/Chirps'
  chr_fls <- gtools::mixedsort(list.files(chr_pth, pattern = '*.tif$', full.names = T))
  chr_dts <- strsplit(x = chr_fls, split = 'chirps-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
  chr_dts <- strsplit(x = chr_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
  chr_dts <- as.Date(gsub('.', '-', chr_dts, fixed = T))
  
  # Tmax
  era5Dir <- '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5'
  tmx_pth <- paste0(era5Dir,'/2m_temperature-24_hour_maximum')
  tmx_fls <- gtools::mixedsort(list.files(tmx_pth, pattern = '*.nc$', full.names = T))
  tmx_dts <- strsplit(x = tmx_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  tmx_dts <- strsplit(x = tmx_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  tmx_dts <- as.Date(tmx_dts, "%Y%m%d")
  
  # Tmin
  tmn_pth <- paste0(era5Dir,'/2m_temperature-24_hour_minimum')
  tmn_fls <- gtools::mixedsort(list.files(tmn_pth, pattern = '*.nc$', full.names = T))
  tmn_dts <- strsplit(x = tmn_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  tmn_dts <- strsplit(x = tmn_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  tmn_dts <- as.Date(tmn_dts, "%Y%m%d")
  
  # Tmean
  tav_pth <- paste0(era5Dir,'/2m_temperature-24_hour_mean')
  tav_fls <- gtools::mixedsort(list.files(tav_pth, pattern = '*.nc$', full.names = T))
  tav_dts <- strsplit(x = tav_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  tav_dts <- strsplit(x = tav_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  tav_dts <- as.Date(tav_dts, "%Y%m%d")
  
  # Solar radiation
  srd_pth <- paste0(era5Dir,'/solar_radiation_flux')
  srd_fls <- gtools::mixedsort(list.files(srd_pth, pattern = '*.nc$', full.names = T))
  srd_dts <- strsplit(x = srd_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  srd_dts <- strsplit(x = srd_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  srd_dts <- as.Date(srd_dts, "%Y%m%d")
  
  # Relative humidity
  rhy_pth <- paste0(era5Dir,'/2m_relative_humidity')
  rhy_fls <- gtools::mixedsort(list.files(rhy_pth, pattern = '*.nc$', full.names = T))
  rhy_dts <- strsplit(x = rhy_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
  rhy_dts <- strsplit(x = rhy_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
  rhy_dts <- as.Date(rhy_dts, "%Y%m%d")
  
  # Filtering days within the season
  yrs <- lubridate::year(tmx_dts)
  yrs <- names(table(yrs)[table(yrs) %in% 365:366])
  tmx_dts <- tmx_dts[lubridate::year(tmx_dts) %in% yrs]
  cnd <- lubridate::month(tmx_dts) %in% season # Days within the season
  yrs_dts <- split(tmx_dts[cnd],cumsum(c(1,diff(tmx_dts[cnd])!=1)))
  yrs_dts <<- yrs_dts[-length(yrs_dts)]
  
  cat('..... Computing: Average temperature.\n')
  AT <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tav <- terra::rast(tav_fls[tav_dts %in% yrs_dts[[i]]])
      tav <- tav %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tav <- tav - 273.15
      AT <- terra::app(x = tav, fun = function(x){ y = mean(x, na.rm = T); return(y) })
      names(AT) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(AT)
    }) %>% terra::rast()
  AT <- AT %>% terra::mask(shp)
  
  cat('..... Computing: Total rainfall.\n')
  TR <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      TR <- terra::app(x = prc, fun = function(x){ y = sum(x, na.rm = T); return(y) })
      names(TR) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(TR)
    }) %>% terra::rast()
  TR <- TR %>% terra::mask(shp)
  
  cat('..... Computing: Number of days with Tmax above 35C.\n')
  NTx35 <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      NTx35 <- terra::app(x = tmx, fun = function(x){ y = calc_htsCMP(tmax = x, t_thresh = 35); return(y) })
      names(NTx35) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(NTx35)
    }) %>% terra::rast()
  NTx35 <- NTx35 %>% terra::mask(shp)
  
  cat('..... Computing: CDD: consecutive dry days.\n')
  CDD <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      CDD <- terra::app(x = prc, fun = function(x){ y = calc_cddCMP(PREC = x, p_thresh = 1); return(y) })
      names(CDD) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(CDD)
    }) %>% terra::rast()
  CDD <- CDD %>% terra::mask(shp)
  
  cat('..... Computing: P5D: rolling average of 5 days.\n')
  P5D <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      P5D <- terra::app(x = prc, fun = function(x){ y = calc_p5dCMP(PREC = x); return(y) })
      names(P5D) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(P5D)
    }) %>% terra::rast()
  P5D <- P5D %>% terra::mask(shp)
  
  cat('..... Computing: P95: percentile 95% of daily precipitation.\n')
  P95 <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      prc <- terra::rast(chr_fls[chr_dts %in% yrs_dts[[i]]])
      prc <- prc %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      prc[prc == -9999] <- NA
      P95 <- terra::app(x = prc, fun = function(x){ y = calc_p95CMP(PREC = x); return(y) })
      names(P95) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(P95)
    }) %>% terra::rast()
  P95 <- P95 %>% terra::mask(shp)
  
  cat('..... Computing: CSDI Cold spell duration Index.\n')
  CSDI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmn <- terra::rast(tmn_fls[tmn_dts %in% yrs_dts[[i]]])
      tmn <- tmn %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmn <- tmn - 273.15
      CSDI <- terra::app(x = tmn, fun = function(x){ y = calc_csdiMP(TMIN = x); return(y) })
      names(CSDI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(CSDI)
    }) %>% terra::rast()
  CSDI <- CSDI %>% terra::mask(shp)
  
  cat('..... Computing: SHI: number of days with maximum temperatures > 29C and relative humidity > 50%.\n')
  SHI <- 1:length(yrs_dts) %>%
    purrr::map(.f = function(i){
      tmx <- terra::rast(tmx_fls[tmx_dts %in% yrs_dts[[i]]])
      tmx <- tmx %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      tmx <- tmx - 273.15
      rhy <- terra::rast(rhy_fls[rhy_dts %in% yrs_dts[[i]]])
      rhy <- rhy %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
      SHI <- terra::lapp(x = terra::sds(tmx, rhy), fun = calc_SHIMP)
      SHI <- sum(SHI)
      names(SHI) <- lubridate::year(yrs_dts[[i]]) %>% unique() %>% paste0(collapse = '-')
      return(SHI)
    }) %>% terra::rast()
  SHI <- SHI %>% terra::mask(shp)
  
  cat('..... End.\n')
  return(list(AT    = AT,
              TR    = TR,
              NTx35 = NTx35,
              CDD   = CDD,
              P5D   = P5D,
              P95   = P95,
              CSDI  = CSDI,
              SHI   = SHI))
  
}

# Loop through seasons
1:length(seasons) %>%
  purrr::map(.f = function(s){
    cat(paste0('Processing season ',names(seasons)[s],':\n'))
    # Indices calculation
    indices <- calc_AgrClm(seasons[[s]], shp_fl)
    # Load 5 km raster template
    tmp <- terra::rast('//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif')
    shp <- terra::vect(shp_fl)
    tmp <- tmp %>% terra::crop(terra::ext(shp)) %>% terra::mask(shp)
    tmp[!is.na(tmp)] <- 1
    # Indices resampling
    indices <- indices %>% purrr::map(.f = function(r){r <- r %>% terra::resample(x = ., y = tmp) %>% terra::mask(shp); return(r)})
    # Saving results
    out <- paste0('D:/india/results/',names(seasons)[s]); if(!dir.exists(out)){dir.create(out,F,T)}
    1:length(names(indices)) %>%
      purrr::map(.f = function(j){
        terra::writeRaster(x = indices[[j]], filename = paste0(out,'/',names(indices)[j],'.tif'), overwrite = T)
      })
    return(cat('Process finished successfully!\n'))
  })
