
#------------------------------------------------------------------------------#
# Extreme values -- ADAPTA 
# H. Achicanoy 
# Alliance Bioversity-CIAT, 2022
#------------------------------------------------------------------------------#
# main function
#------------------------------------------------------------------------------#
ExtVal <- function(cdd1 = infile, outfile = outfile){
  cdd1_ExtVal_rp2  <- terra::app(x = infile, fun = function(x){
    tryCatch(expr = {
      if(sum(is.na(x)) == length(x)){
        out <- NA
      } else {
        thr <- mean(x, na.rm = T)
        fit <- extRemes::fevd(x = x, threshold = thr, type = "GP", span = length(x), time.units = "years", period.basis = "year", verbose = F, na.action = na.exclude)
        out <- return.level(fit, 2) %>% as.numeric()
      }
    },
    error = function(e){
      cat("Modeling process failed\n")
      return("Done\n")
    })
    if(exists('out')){return(out)} else {out <- NA}
  })
  cdd1_ExtVal_rp5  <- terra::app(x = infile, fun = function(x){
    tryCatch(expr = {
      if(sum(is.na(x)) == length(x)){
        out <- NA
      } else {
        thr <- mean(x, na.rm = T)
        fit <- extRemes::fevd(x = x, threshold = thr, type = "GP", span = length(x), time.units = "years", period.basis = "year", verbose = F, na.action = na.exclude)
        out <- return.level(fit, 5) %>% as.numeric()
      }
    },
    error = function(e){
      cat("Modeling process failed\n")
      return("Done\n")
    })
    if(exists('out')){return(out)} else {out <- NA}
  })
  cdd1_ExtVal_rp10 <- terra::app(x = infile, fun = function(x){
    tryCatch(expr = {
      if(sum(is.na(x)) == length(x)){
        out <- NA
      } else {
        thr <- mean(x, na.rm = T)
        fit <- extRemes::fevd(x = x, threshold = thr, type = "GP", span = length(x), time.units = "years", period.basis = "year", verbose = F, na.action = na.exclude)
        out <- return.level(fit, 10) %>% as.numeric()
      }
    },
    error = function(e){
      cat("Modeling process failed\n")
      return("Done\n")
    })
    if(exists('out')){return(out)} else {out <- NA}
  })
  
  cdd1_ExtVal <- terra::rast(list(cdd1_ExtVal_rp2, cdd1_ExtVal_rp5, cdd1_ExtVal_rp10))
  names(cdd1_ExtVal) <- paste0('Return.period.',c(2,5,10),'y')
  
  terra::writeRaster(x = cdd1_ExtVal, outfile, overwrite = T)
  return(cdd1_ExtVal)
}
#------------------------------------------------------------------------------#