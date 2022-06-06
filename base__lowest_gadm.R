lowest_gadm <- function(iso = 'KEN', out = NULL){
  suppressMessages(if(!require(pacman)){install.packages('pacman');library(pacman)})
  suppressMessages(pacman::p_load(geodata,terra,sf))
  levels <- 5:1
  for(i in 1:length(levels)){
    tryCatch(expr = {
      shp <- geodata::gadm(country = iso, level = levels[i], path = tempdir(), version = '4.0')
      break
    },
    error = function(e){
      cat(paste0("Getting GADM level ",levels[i]," failed... Trying a higher level\n"))
      return("\n")
    })
  }
  if(!is.null(out)){
    shp <- sf::st_as_sf(shp)
    sf::st_write(shp, paste0(out,'/',iso,'.gpkg'))
  }
  shp <- as(shp, 'Spatial')
  return(shp)
}