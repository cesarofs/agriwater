#SINGLE VALUES OF AGROMETEOROLOGICAL DATA
#SENTINEL-2

#'
#' Crop coefficient (ETa / ET0) using Sentinel-2 images with single agrometeorological data.
#'
#' @param doy   is the Day of Year (DOY)
#' @param RG  is the global solar radiation
#' @param Ta  is the average air temperature
#' @param a  is one of the regression coefficients of SAFER algorithm
#' @param b  is one of the regression coefficients of SAFER algorithm
#' @export
#' @import terra
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), Crop Coefficient ("kc") and net radiation ("Rn_MJ").
#' @examples
#' library(agriwater)
#'
#' # dependencies of package 'agriwater'
#' library(terra)
#'
#' # Using a temporary folder to run example
#' wd <- tempdir()
#' initial = getwd()
#' setwd(wd)
#'
#' # creating raster which simulate Sentinel-2 reflectances - for using
#' # real data, please download:
#' # https://drive.google.com/open?id=14E1wHNLxG7_Dh4I-GqNYakj8YJDgKLzk
#'
#' xy <- matrix(rnorm(4, mean = 0.07, sd = 0.01), 2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B2.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B3.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.03, sd = 0.018),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B4.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B8.tif"),filetype = "GTiff", overwrite=TRUE)
#  # creating mask of study area
#' mask <- as.polygons(rast)
#' writeVector(mask, file.path(getwd(),"mask.shp"), overwrite=TRUE)
#'
#' # using "agriwater"
#' kc_s2(doy = 134, RG = 17.6, Ta = 27.9, a = 1.8, b = -0.008)
#'
#' #Exiting temporary folder and returning to previous workspace
#' setwd(initial)

kc_s2 = function(doy, RG, Ta, a, b){

  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b8 <- rast("B8.tif")

  mask <- vect("mask.shp")


  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b8_crop <- crop(b8, ext(mask)[1:4])
  b8_mascara <- mask(b8_crop, mask)

  b2_mascara <- b2_mascara/10000
  b3_mascara <- b3_mascara/10000
  b4_mascara <- b4_mascara/10000
  b8_mascara <- b8_mascara/10000

  Alb_Top = b2_mascara*0.32+b3_mascara*0.26+b4_mascara*0.25+b8_mascara*0.17

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149


  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)

  NDVI =(b8_mascara-b4_mascara)/(b8_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", filetype = "GTiff", overwrite=TRUE)

  lati <- long <- b2_mascara
  xy <- crds(b2_mascara)
  long[] <- xy[, 1]
  long <- crop(long, ext(mask)[1:4])
  lati[] <- xy[, 2]
  lati <- crop(lati, ext(mask)[1:4])


  map1 <- (long/long)*((2*pi)/365)*(doy-1)

  Et <- (0.000075+0.001868*cos(map1)-0.032077*sin(map1)-0.014615*cos(2*map1)-0.04089*sin(2*map1))

  LAT <- (13+(4*long/60)+(Et/60))

  Dec <- 0.006918-0.399912*cos(map1)+0.070257*sin(map1)+0.006758*cos(2*map1)+0.000907*sin(2*map1)-0.002697*cos(3*map1)+0.00148*sin(3*map1)

  W <- 15*(LAT-12)*(pi/180)

  cos_zwn <- sin(lati*pi/180)*sin(Dec)+cos(lati*pi/180)*cos(Dec)*cos(W)

  E0 <- (1.00011+0.034221*cos(map1)+0.00128*sin(map1)+0.000719*cos(2*map1)+0.000077*sin(2*map1))

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b2_mascara, method="bilinear")

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*Ta-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", filetype = "GTiff", overwrite=TRUE)

  slope =(4098*(0.6108*exp((17.27*(Ta))/((Ta)+237.3)))/((Ta)+237.3)^2)

  LEeq = (slope*Rn)/(slope+0.066)

  rm(b2, b3, b4, b8, b2_mascara, b3_mascara, b4_mascara, b8_mascara, slope, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

  RR =Alb_24*RG

  Emiss_atm = 0.9364*(((-1)*log(Transm))^0.1135)

  Emiss_atm[Emiss_atm > 1] <- 1

  RLdown_wm2 =(Emiss_atm*5.67*(10^(-8))*((Ta +273.15)^4))

  RL_down =RLdown_wm2/11.6

  RL_up =(RG-RR+RL_down-Rn_MJ)

  Esurf_r1 <- NDVI

  Esurf_r1[NDVI < 0] <- 1

  Esurf_r1[NDVI >= 0] <- NA

  Esurf_r2 <- 1.0035+0.0589*log(NDVI)

  Esurf <- merge(Esurf_r1, Esurf_r2)

  TS24 =((RL_up*11.6)/((Esurf*5.67)*(10^(-8))))^(0.25)

  TS24[TS24 < 273.15] = NA

  writeRaster(TS24, "LST", filetype = "GTiff", overwrite=TRUE)

  NDVI[NDVI <= 0] = NA

  kc=exp((a)+(b*((TS24-273.15)/(Alb_24*NDVI))))

  writeRaster(kc, "kc", filetype = "GTiff", overwrite=TRUE)
}

#'
#' Actual evapotranspiration (ETa) using Sentinel-2 images with single agrometeorological data.
#' @param doy  is the Day of Year (DOY)
#' @param RG is the global solar radiation
#' @param Ta is the average air temperature
#' @param ET0 is the reference evapotranspiration
#' @param a  is one of the regression coefficients of SAFER algorithm
#' @param b is one of the regression coefficients of SAFER algorithm
#' @export
#' @import terra
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), net radiation ("Rn_MJ"), Crop Coefficient ("kc") and Actual Evapotranspiration (evapo).
#' @examples
#' library(agriwater)
#'
#' # dependencies of package 'agriwater'
#' library(terra)
#'
#'
#' # Using a temporary folder to run example
#' wd <- tempdir()
#' initial = getwd()
#' setwd(wd)
#'
#' # creating raster which simulate Sentinel-2 reflectances - for using
#' # real data, please download:
#' # https://drive.google.com/open?id=14E1wHNLxG7_Dh4I-GqNYakj8YJDgKLzk
#'
#' xy <- matrix(rnorm(4, mean = 0.07, sd = 0.01), 2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B2.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B3.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.03, sd = 0.018),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B4.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B8.tif"),filetype = "GTiff", overwrite=TRUE)
#  # creating mask of study area
#' mask <- as.polygons(rast)
#' writeVector(mask, file.path(getwd(),"mask.shp"), overwrite=TRUE)
#'
#' # using "agriwater"
#' evapo_s2(doy = 134, RG = 17.6, Ta = 27.9, ET0 = 3.8, a = 1.8, b = -0.008)
#'
#' #Exiting temporary folder and returning to previous workspace
#' setwd(initial)

evapo_s2 = function(doy, RG, Ta, ET0, a, b){

  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b8 <- rast("B8.tif")

  mask <- vect("mask.shp")


  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b8_crop <- crop(b8, ext(mask)[1:4])
  b8_mascara <- mask(b8_crop, mask)

  b2_mascara <- b2_mascara/10000
  b3_mascara <- b3_mascara/10000
  b4_mascara <- b4_mascara/10000
  b8_mascara <- b8_mascara/10000

  Alb_Top = b2_mascara*0.32+b3_mascara*0.26+b4_mascara*0.25+b8_mascara*0.17

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149


  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)

  NDVI =(b8_mascara-b4_mascara)/(b8_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", filetype = "GTiff", overwrite=TRUE)


  lati <- long <- b2_mascara
  xy <- crds(b2_mascara)
  long[] <- xy[, 1]
  long <- crop(long, ext(mask)[1:4])
  lati[] <- xy[, 2]
  lati <- crop(lati, ext(mask)[1:4])


  map1 <- (long/long)*((2*pi)/365)*(doy-1)

  Et <- (0.000075+0.001868*cos(map1)-0.032077*sin(map1)-0.014615*cos(2*map1)-0.04089*sin(2*map1))

  LAT <- (13+(4*long/60)+(Et/60))

  Dec <- 0.006918-0.399912*cos(map1)+0.070257*sin(map1)+0.006758*cos(2*map1)+0.000907*sin(2*map1)-0.002697*cos(3*map1)+0.00148*sin(3*map1)

  W <- 15*(LAT-12)*(pi/180)

  cos_zwn <- sin(lati*pi/180)*sin(Dec)+cos(lati*pi/180)*cos(Dec)*cos(W)

  E0 <- (1.00011+0.034221*cos(map1)+0.00128*sin(map1)+0.000719*cos(2*map1)+0.000077*sin(2*map1))

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b2_mascara, method="bilinear")

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*Ta-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", filetype = "GTiff", overwrite=TRUE)

  slope =(4098*(0.6108*exp((17.27*(Ta))/((Ta)+237.3)))/((Ta)+237.3)^2)

  LEeq = (slope*Rn)/(slope+0.066)

  rm(b2, b3, b4, b8, b2_mascara, b3_mascara, b4_mascara, b8_mascara, slope, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

  RR =Alb_24*RG

  Emiss_atm = 0.9364*(((-1)*log(Transm))^0.1135)

  Emiss_atm[Emiss_atm > 1] <- 1

  RLdown_wm2 =(Emiss_atm*5.67*(10^(-8))*((Ta +273.15)^4))

  RL_down =RLdown_wm2/11.6

  RL_up =(RG-RR+RL_down-Rn_MJ)

  Esurf_r1 <- NDVI

  Esurf_r1[NDVI < 0] <- 1

  Esurf_r1[NDVI >= 0] <- NA

  Esurf_r2 <- 1.0035+0.0589*log(NDVI)

  Esurf <- merge(Esurf_r1, Esurf_r2)

  TS24 =((RL_up*11.6)/((Esurf*5.67)*(10^(-8))))^(0.25)

  TS24[TS24 < 273.15] = NA

  writeRaster(TS24, "LST", filetype = "GTiff", overwrite=TRUE)

  NDVI[NDVI <= 0] = NA

  kc=exp((a)+(b*((TS24-273.15)/(Alb_24*NDVI))))

  writeRaster(kc, "kc", filetype = "GTiff", overwrite=TRUE)

  ET=kc*ET0

  writeRaster(ET, "evapo", filetype = "GTiff", overwrite=TRUE)
}

#'Energy balance using Sentinel-2 images with single agrometeorological data.
#'@param doy is the Day of Year (DOY)
#'@param RG is the global solar radiation
#'@param Ta is the average air temperature
#'@param ET0  is the reference evapotranspiration
#'@param a is one of the regression coefficients of SAFER algorithm
#'@param b is one of the regression coefficients of SAFER algorithm
#'@export
#' @import terra
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), Crop Coefficient ("kc"), Actual Evapotranspiration (evapo), latent heat flux "LE_MJ"), net radiation ("Rn_MJ"), ground heat flux ("G_MJ") and the sensible heat flux ("H_MJ").
#' @examples
#' library(agriwater)
#'
#' # dependencies of package 'agriwater'
#' library(terra)
#'
#'
#' # Using a temporary folder to run example
#' wd <- tempdir()
#' initial = getwd()
#' setwd(wd)
#'
#' # creating raster which simulate Sentinel-2 reflectances - for using
#' # real data, please download:
#' # https://drive.google.com/open?id=14E1wHNLxG7_Dh4I-GqNYakj8YJDgKLzk
#'
#' xy <- matrix(rnorm(4, mean = 0.07, sd = 0.01), 2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B2.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B3.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.03, sd = 0.018),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B4.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B8.tif"),filetype = "GTiff", overwrite=TRUE)
#  # creating mask of study area
#' mask <- as.polygons(rast)
#' writeVector(mask, file.path(getwd(),"mask.shp"), overwrite=TRUE)
#'
#' # using "agriwater"
#' radiation_s2(doy = 134, RG = 17.6, Ta = 27.9, ET0 = 3.8, a = 1.8, b = -0.008)
#'
#' #Exiting temporary folder and returning to previous workspace
#' setwd(initial)

radiation_s2 =  function(doy, RG, Ta, ET0, a, b){


  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b8 <- rast("B8.tif")

  mask <- vect("mask.shp")


  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b8_crop <- crop(b8, ext(mask)[1:4])
  b8_mascara <- mask(b8_crop, mask)

  b2_mascara <- b2_mascara/10000
  b3_mascara <- b3_mascara/10000
  b4_mascara <- b4_mascara/10000
  b8_mascara <- b8_mascara/10000

  Alb_Top = b2_mascara*0.32+b3_mascara*0.26+b4_mascara*0.25+b8_mascara*0.17

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149

  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)

  NDVI =(b8_mascara-b4_mascara)/(b8_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", filetype = "GTiff", overwrite=TRUE)


  lati <- long <- b2_mascara
  xy <- crds(b2_mascara)
  long[] <- xy[, 1]
  long <- crop(long, ext(mask)[1:4])
  lati[] <- xy[, 2]
  lati <- crop(lati, ext(mask)[1:4])


  map1 <- (long/long)*((2*pi)/365)*(doy-1)

  Et <- (0.000075+0.001868*cos(map1)-0.032077*sin(map1)-0.014615*cos(2*map1)-0.04089*sin(2*map1))

  LAT <- (13+(4*long/60)+(Et/60))

  Dec <- 0.006918-0.399912*cos(map1)+0.070257*sin(map1)+0.006758*cos(2*map1)+0.000907*sin(2*map1)-0.002697*cos(3*map1)+0.00148*sin(3*map1)

  W <- 15*(LAT-12)*(pi/180)

  cos_zwn <- sin(lati*pi/180)*sin(Dec)+cos(lati*pi/180)*cos(Dec)*cos(W)

  E0 <- (1.00011+0.034221*cos(map1)+0.00128*sin(map1)+0.000719*cos(2*map1)+0.000077*sin(2*map1))

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b2_mascara, method="bilinear")

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*Ta-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", filetype = "GTiff", overwrite=TRUE)

  slope =(4098*(0.6108*exp((17.27*(Ta))/((Ta)+237.3)))/((Ta)+237.3)^2)

  LEeq = (slope*Rn)/(slope+0.066)

  rm(b2, b3, b4, b8, b2_mascara, b3_mascara, b4_mascara, b8_mascara, slope, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

  RR =Alb_24*RG

  Emiss_atm = 0.9364*(((-1)*log(Transm))^0.1135)

  Emiss_atm[Emiss_atm > 1] <- 1

  RLdown_wm2 =(Emiss_atm*5.67*(10^(-8))*((Ta +273.15)^4))

  RL_down =RLdown_wm2/11.6

  RL_up =(RG-RR+RL_down-Rn_MJ)

  Esurf_r1 <- NDVI

  Esurf_r1[NDVI < 0] <- 1

  Esurf_r1[NDVI >= 0] <- NA

  Esurf_r2 <- 1.0035+0.0589*log(NDVI)

  Esurf <- merge(Esurf_r1, Esurf_r2)

  TS24 =((RL_up*11.6)/((Esurf*5.67)*(10^(-8))))^(0.25)

  TS24[TS24 < 273.15] = NA

  writeRaster(TS24, "LST", filetype = "GTiff", overwrite=TRUE)

  NDVI[NDVI <= 0] = NA

  kc=exp((a)+(b*((TS24-273.15)/(Alb_24*NDVI))))

  ET=kc*ET0

  LE_MJ =ET*2.45

  writeRaster(LE_MJ, "LE_MJ", filetype = "GTiff", overwrite=TRUE)

  G_Rn =3.98*exp(-25.47*Alb_24)

  G_MJ =G_Rn*Rn_MJ

  writeRaster(G_MJ, "G_MJ", filetype = "GTiff", overwrite=TRUE)

  H_MJ =Rn_MJ-LE_MJ-G_MJ

  writeRaster(H_MJ, "H_MJ", filetype = "GTiff", overwrite=TRUE)

}

#' Surface Albedo using Sentinel-2 images.
#' @export
#' @import terra
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24").
#' @examples
#' library(agriwater)
#'
#' # dependencies of package 'agriwater'
#' library(terra)
#'
#' # Using a temporary folder to run example
#' wd <- tempdir()
#' initial = getwd()
#' setwd(wd)
#'
#' # creating raster which simulate Sentinel-2 reflectances - for using
#' # real data, please download:
#' # https://drive.google.com/open?id=14E1wHNLxG7_Dh4I-GqNYakj8YJDgKLzk
#'
#' xy <- matrix(rnorm(4, mean = 0.07, sd = 0.01), 2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B2.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B3.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.03, sd = 0.018),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B4.tif"),filetype = "GTiff", overwrite=TRUE)
#' xy <- matrix(rnorm(4, mean = 0.05, sd = 0.015),2, 2)
#' rast <- rast(xy, crs="+proj=longlat +datum=WGS84")
#' ext(rast) <- c(-40.5,-40.45,-9.5,-9.45)
#' writeRaster(rast, file.path(wd, "B8.tif"),filetype = "GTiff", overwrite=TRUE)
#  # creating mask of study area
#' mask <- as.polygons(rast)
#' writeVector(mask, file.path(getwd(),"mask.shp"), overwrite=TRUE)
#'
#' # using "agriwater"
#' albedo_s2()
#'
#' #Exiting temporary folder and returning to previous workspace
#' setwd(initial)

albedo_s2 = function(){

  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b8 <- rast("B8.tif")

  mask <- vect("mask.shp")


  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b8_crop <- crop(b8, ext(mask)[1:4])
  b8_mascara <- mask(b8_crop, mask)

  b2_mascara <- b2_mascara/10000
  b3_mascara <- b3_mascara/10000
  b4_mascara <- b4_mascara/10000
  b8_mascara <- b8_mascara/10000

  Alb_Top = b2_mascara*0.32+b3_mascara*0.26+b4_mascara*0.25+b8_mascara*0.17

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149


  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)
}
