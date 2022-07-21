# R options
g <- gc(reset = T); rm(list = ls()) # Empty garbage collector
# .rs.restartR()                      # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, terra, lubridate))

root <- '//CATALOGUE/Workspace14/WFP_ClimateRiskPr/1.Data/climate/CMIP6'

gcm <- 'ACCESS-ESM1-5'
var <- 'tasmin' # pr, tasmax, tasmin
prd <- c(2021, 2040)
iso <- 'MLI'
out <- '//catalogue/Workspace14/WFP_ClimateRiskPr/0.Project_Documents'

bc_delta <- function(gcm, var, prd, iso, out){
  
  cat(paste0('Processing: model ',gcm,', variable ',var,', future period ',prd[1],'-',prd[2],' in ',iso,'\n'))
  
  ## --------------------------------------------------------------------------------- ##
  ## Baseline
  ## --------------------------------------------------------------------------------- ##
  
  # Lists all the historical files
  bsl <- c(1995, 2014)
  fls <- list.files(path = paste0(root,'/download_data'), pattern = gcm, full.names = T)
  fls_his <- grep(pattern = 'historical', x = fls, value = T)
  fls_his <- grep(pattern = var, x = fls_his, value = T)
  
  # Identify the right files
  dts <- strsplit(x = fls_his, split = '_gn_') %>% purrr::map(2) %>% as.character()
  dts <- gsub(pattern = '.nc', replacement = '', x = dts, fixed = T)
  dts <- strsplit(x = dts, split = '-')
  ini <- dts %>% purrr::map(.f = function(vct){
    vct <- as.Date(vct, "%Y%m%d")
    chck <- as.Date('1995-01-01') %within% lubridate::interval(vct[1], vct[2])
    return(chck)
  }) %>% unlist() %>% which()
  end <- dts %>% purrr::map(.f = function(vct){
    vct <- as.Date(vct, "%Y%m%d")
    chck <- as.Date('2014-12-31') %within% lubridate::interval(vct[1], vct[2])
    return(chck)
  }) %>% unlist() %>% which()
  
  # Load the right files
  fls_his <- fls_his[ini:end]
  gcm_his <- fls_his %>%
    purrr::map(.f = function(fl){
      gcm <- terra::rast(fl)
      gcm <- gcm[[lubridate::year(terra::time(gcm)) %in% 1995:2014]]
      return(gcm)
    }) %>% terra::rast()
  
  ## --------------------------------------------------------------------------------- ##
  ## Future period
  ## --------------------------------------------------------------------------------- ##
  
  # Lists all the future files
  fls <- list.files(path = paste0(root,'/download_data'), pattern = gcm, full.names = T)
  fls_fut <- grep(pattern = 'ssp585', x = fls, value = T)
  fls_fut <- grep(pattern = var, x = fls_fut, value = T)
  
  # Identify the right files
  dts <- strsplit(x = fls_fut, split = '_gn_') %>% purrr::map(2) %>% as.character()
  dts <- gsub(pattern = '.nc', replacement = '', x = dts, fixed = T)
  dts <- strsplit(x = dts, split = '-')
  ini <- dts %>% purrr::map(.f = function(vct){
    vct <- as.Date(vct, "%Y%m%d")
    chck <- as.Date(paste0(prd[1],'-01-01')) %within% lubridate::interval(vct[1], vct[2])
    return(chck)
  }) %>% unlist() %>% which()
  end <- dts %>% purrr::map(.f = function(vct){
    vct <- as.Date(vct, "%Y%m%d")
    chck <- as.Date(paste0(prd[2],'-01-01')) %within% lubridate::interval(vct[1], vct[2])
    return(chck)
  }) %>% unlist() %>% which()
  
  # Load the right files
  fls_fut <- fls_fut[ini:end]
  gcm_fut <- fls_fut %>%
    purrr::map(.f = function(fl){
      gcm <- terra::rast(fl)
      gcm <- gcm[[lubridate::year(terra::time(gcm)) %in% prd[1]:prd[2]]]
      return(gcm)
    }) %>% terra::rast()
  
  # Country mask
  shp <- geodata::gadm(country = iso, level = 0, path = base::tempdir())
  
  result <- 1:12 %>%
    purrr::map(.f = function(mnth){
      
      cat(paste0('Processing month: ',mnth,'\n'))
      
      ## Compute monthly statistics for historical period
      # Days within the month
      cnd <- lubridate::month(terra::time(gcm_his)) == mnth
      yrs_dts <<- split(terra::time(gcm_his)[cnd],cumsum(c(1,diff(terra::time(gcm_his)[cnd])!=1)))
      
      his_stack <- yrs_dts %>%
        purrr::map(.f = function(flt){
          if(var %in% c('tasmax','tasmin')){
            avg <- mean(gcm_his[[terra::time(gcm_his) %in% flt]])
          } else {
            avg <- sum(gcm_his[[terra::time(gcm_his) %in% flt]])
          }
          return(avg)
        })
      his_stack <- terra::rast(his_stack)
      avg_his <- mean(his_stack)
      
      ## Compute monthly statistics for future period
      # Days within the month
      cnd <- lubridate::month(terra::time(gcm_fut)) == mnth # Days within the month
      yrs_dts <- split(terra::time(gcm_fut)[cnd],cumsum(c(1,diff(terra::time(gcm_fut)[cnd])!=1)))
      
      fut_stack <- yrs_dts %>%
        purrr::map(.f = function(flt){
          if(var %in% c('tasmax','tasmin')){
            avg <- mean(gcm_fut[[terra::time(gcm_fut) %in% flt]])
          } else {
            avg <- sum(gcm_fut[[terra::time(gcm_fut) %in% flt]])
          }
          return(avg)
        })
      fut_stack <- terra::rast(fut_stack)
      avg_fut <- mean(fut_stack)
      
      if(var %in% c('tasmax','tasmin')){
        avg_his <- avg_his - 273.15
        avg_fut <- avg_fut - 273.15
      } else {
        avg_his <- avg_his * 86400
        avg_fut <- avg_fut * 86400
      }
      
      avg_his <- terra::rotate(avg_his)
      avg_fut <- terra::rotate(avg_fut)
      
      if(var %in% c('tasmax','tasmin')){
        anom <- avg_fut - avg_his
      } else {
        anom <- (avg_fut - avg_his)/avg_his
      }
      
      anom <- anom %>% terra::crop(terra::ext(shp), snap = 'out')#  %>% terra::mask(shp, touches = T)
      
      crds <- anom %>% terra::as.data.frame(xy = T)
      
      empty <- anom
      terra::values(empty) <- NA
      
      anomalies_values <- unique(crds[,'mean'])
      
      tps  <- fields::Tps(x = crds[,c('x','y')], Y = crds[,'mean'])
      
      ref <- terra::rast("//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/chirps-v2.0.2020.01.01.tif")
      ref <- ref %>% terra::crop(terra::ext(shp), snap = 'out') %>% terra::mask(shp, touches = T)
      
      intp <- raster::interpolate(raster::raster(ref), tps)
      intp <- terra::rast(intp) %>% terra::mask(mask = ref)
      
      era5Dir <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/ERA5'
      chrpDir <- '//catalogue/Workspace14/WFP_ClimateRiskPr/1.Data/Chirps'
      
      if(var == 'tasmin'){
        obs_pth <- paste0(era5Dir,'/2m_temperature-24_hour_minimum')
      } else {
        if(var == 'tasmax'){
          obs_pth <- paste0(era5Dir,'/2m_temperature-24_hour_maximum')
        } else {
          if(var == 'pr'){
            obs_pth <- chrpDir
          }
        }
      }
      
      if(var %in% c('tasmin','tasmax')){
        obs_fls <- gtools::mixedsort(list.files(obs_pth, pattern = '*.nc$', full.names = T))
        obs_dts <- strsplit(x = obs_fls, split = 'glob-agric_AgERA5_', fixed = T) %>% purrr::map(2) %>% unlist()
        obs_dts <- strsplit(x = obs_dts, split = '_final-v1.0.nc', fixed = T) %>% purrr::map(1) %>% unlist()
        obs_dts <- as.Date(obs_dts, "%Y%m%d")
        obs_fls <- obs_fls[lubridate::year(obs_dts) %in% 1995:2014]
        obs_dts <- obs_dts[lubridate::year(obs_dts) %in% 1995:2014]
      } else {
        obs_fls <- gtools::mixedsort(list.files(obs_pth, pattern = '*.tif$', full.names = T))
        obs_dts <- strsplit(x = obs_fls, split = 'chirps-v2.0.', fixed = T) %>% purrr::map(2) %>% unlist()
        obs_dts <- strsplit(x = obs_dts, split = '.tif', fixed = T) %>% purrr::map(1) %>% unlist()
        obs_dts <- as.Date(gsub('.', '-', obs_dts, fixed = T))
      }
      
      # Days within the season
      cnd <- lubridate::month(obs_dts) == mnth
      yrs_dts <- split(obs_dts[cnd],cumsum(c(1,diff(obs_dts[cnd])!=1)))
      
      fut_ds_bc <- yrs_dts %>%
        purrr::map(.f = function(flt){
          
          obs <- terra::rast(obs_fls[obs_dts %in% flt])
          obs <- obs %>% terra::crop(terra::ext(shp), snap = 'out') %>% terra::mask(shp, touches = T)
          if(var %in% c('tasmin','tasmax')){
            obs <- obs - 273.15
          }
          obs <- terra::resample(x = obs, y = ref)
          
          if(var %in% c('tasmin','tasmax')){
            fut_ds_bc <- obs + intp
          } else {
            fut_ds_bc <- obs * (1 + intp)
          }
          return(fut_ds_bc)
          
        }) %>% terra::rast()
      
      return(fut_ds_bc)
      
    }) %>%
    terra::rast()
  result <- result[[order(terra::time(result))]]
  
  terra::writeRaster(x = result, filename = paste0(out,'/',iso,'_',var,'_',prd[1],'-',prd[2],'.tif'), overwrite = T)
  
}
bc_delta(gcm, var, prd, iso, out)