#GRID OF AGROMETEOROLOGICAL DATA
#LANDSAT-8

#' Crop coefficient (ETa / ET0) using Landsat-8 images with a grid of agrometeorological data.
#'
#' @param doy   is the Day of Year (DOY)
#' @param a  is one of the regression coefficients of SAFER algorithm
#' @param b  is one of the regression coefficients of SAFER algorithm
#' @export
#' @import terra
#' @importFrom utils read.csv
#'
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), Crop Coefficient ("kc") and net radiation ("Rn_MJ").

kc_l8_grid  = function(doy, a, b){

  b1 <- rast("B1.tif")
  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b5 <- rast("B5.tif")
  b6 <- rast("B6.tif")
  b7 <- rast("B7.tif")
  RG <- rast("RG.tif")
  Ta <- rast("Ta.tif")

  mask <- vect("mask.shp")

  b1_crop <- crop(b1, ext(mask)[1:4])
  b1_mascara <- mask(b1_crop, mask)
  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b5_crop <- crop(b5, ext(mask)[1:4])
  b5_mascara <- mask(b5_crop, mask)
  b6_crop <- crop(b6, ext(mask)[1:4])
  b6_mascara <- mask(b6_crop, mask)
  b7_crop <- crop(b7, ext(mask)[1:4])
  b7_mascara <- mask(b7_crop, mask)
  RG_crop <- crop(RG, ext(mask)[1:4])
  RG <- mask(RG_crop, mask)
  Ta_crop <- crop(Ta, ext(mask)[1:4])
  Ta <- mask(Ta_crop, mask)

  metadata <- list.files(pattern = "txt")
  m <- read.csv(metadata, header = T)
  metadata <- paste( unlist(m), collapse='')
  rm(m)

  RADIANCE_MAXIMUM_BAND_1 <- "^.*RADIANCE_MAXIMUM_BAND_1 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_1 <- gsub(RADIANCE_MAXIMUM_BAND_1, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_1 <- as.numeric(RADIANCE_MAXIMUM_BAND_1)

  RADIANCE_MAXIMUM_BAND_2 <- "^.*RADIANCE_MAXIMUM_BAND_2 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_2 <- gsub(RADIANCE_MAXIMUM_BAND_2, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_2 <- as.numeric(RADIANCE_MAXIMUM_BAND_2)

  RADIANCE_MAXIMUM_BAND_3 <- "^.*RADIANCE_MAXIMUM_BAND_3 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_3 <- gsub(RADIANCE_MAXIMUM_BAND_3, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_3 <- as.numeric(RADIANCE_MAXIMUM_BAND_3)

  RADIANCE_MAXIMUM_BAND_4 <- "^.*RADIANCE_MAXIMUM_BAND_4 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_4 <- gsub(RADIANCE_MAXIMUM_BAND_4, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_4 <- as.numeric(RADIANCE_MAXIMUM_BAND_4)

  RADIANCE_MAXIMUM_BAND_5 <- "^.*RADIANCE_MAXIMUM_BAND_5 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_5 <- gsub(RADIANCE_MAXIMUM_BAND_5, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_5 <- as.numeric(RADIANCE_MAXIMUM_BAND_5)

  RADIANCE_MAXIMUM_BAND_6 <- "^.*RADIANCE_MAXIMUM_BAND_6 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_6 <- gsub(RADIANCE_MAXIMUM_BAND_6, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_6 <- as.numeric(RADIANCE_MAXIMUM_BAND_6)

  RADIANCE_MAXIMUM_BAND_7 <- "^.*RADIANCE_MAXIMUM_BAND_7 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_7 <- gsub(RADIANCE_MAXIMUM_BAND_7, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_7 <- as.numeric(RADIANCE_MAXIMUM_BAND_7)

  RADIANCE_MINIMUM_BAND_1 <- "^.*RADIANCE_MINIMUM_BAND_1 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_1 <- gsub(RADIANCE_MINIMUM_BAND_1, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_1 <- as.numeric(RADIANCE_MINIMUM_BAND_1)

  RADIANCE_MINIMUM_BAND_2 <- "^.*RADIANCE_MINIMUM_BAND_2 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_2 <- gsub(RADIANCE_MINIMUM_BAND_2, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_2 <- as.numeric(RADIANCE_MINIMUM_BAND_2)

  RADIANCE_MINIMUM_BAND_3 <- "^.*RADIANCE_MINIMUM_BAND_3 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_3 <- gsub(RADIANCE_MINIMUM_BAND_3, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_3 <- as.numeric(RADIANCE_MINIMUM_BAND_3)

  RADIANCE_MINIMUM_BAND_4 <- "^.*RADIANCE_MINIMUM_BAND_4 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_4 <- gsub(RADIANCE_MINIMUM_BAND_4, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_4 <- as.numeric(RADIANCE_MINIMUM_BAND_4)

  RADIANCE_MINIMUM_BAND_5 <- "^.*RADIANCE_MINIMUM_BAND_5 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_5 <- gsub(RADIANCE_MINIMUM_BAND_5, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_5 <- as.numeric(RADIANCE_MINIMUM_BAND_5)

  RADIANCE_MINIMUM_BAND_6 <- "^.*RADIANCE_MINIMUM_BAND_6 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_6 <- gsub(RADIANCE_MINIMUM_BAND_6, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_6 <- as.numeric(RADIANCE_MINIMUM_BAND_6)

  RADIANCE_MINIMUM_BAND_7 <- "^.*RADIANCE_MINIMUM_BAND_7 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_7 <- gsub(RADIANCE_MINIMUM_BAND_7, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_7 <- as.numeric(RADIANCE_MINIMUM_BAND_7)

  b1_mascara = ((RADIANCE_MAXIMUM_BAND_1+RADIANCE_MINIMUM_BAND_1)/65535)*b1_mascara-(RADIANCE_MINIMUM_BAND_1)

  b2_mascara = ((RADIANCE_MAXIMUM_BAND_2+RADIANCE_MINIMUM_BAND_2)/65535)*b2_mascara-(RADIANCE_MINIMUM_BAND_2)

  b3_mascara = ((RADIANCE_MAXIMUM_BAND_3+RADIANCE_MINIMUM_BAND_3)/65535)*b3_mascara-(RADIANCE_MINIMUM_BAND_3)

  b4_mascara = ((RADIANCE_MAXIMUM_BAND_4+RADIANCE_MINIMUM_BAND_4)/65535)*b4_mascara-(RADIANCE_MINIMUM_BAND_4)

  b5_mascara = ((RADIANCE_MAXIMUM_BAND_5+RADIANCE_MINIMUM_BAND_5)/65535)*b5_mascara-(RADIANCE_MINIMUM_BAND_5)

  b6_mascara = ((RADIANCE_MAXIMUM_BAND_6+RADIANCE_MINIMUM_BAND_6)/65535)*b6_mascara-(RADIANCE_MINIMUM_BAND_6)

  b7_mascara = ((RADIANCE_MAXIMUM_BAND_7+RADIANCE_MINIMUM_BAND_7)/65535)*b7_mascara-(RADIANCE_MINIMUM_BAND_7)


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

  b1_mascara=(pi*b1_mascara)/(1718.75*cos_zwn*E0)

  b1_mascara[b1_mascara <= 0] = 0

  b2_mascara=(pi*b2_mascara)/(1810.42*cos_zwn*E0)

  b2_mascara[b2_mascara <= 0] = 0

  b3_mascara=(pi*b3_mascara)/(1741.67*cos_zwn*E0)

  b3_mascara[b3_mascara <= 0] = 0

  b4_mascara=(pi*b4_mascara)/(1558.33*cos_zwn*E0)

  b4_mascara[b4_mascara <= 0] = 0

  b5_mascara=(pi*b5_mascara)/(962.5*cos_zwn*E0)

  b5_mascara[b5_mascara <= 0] = 0

  b6_mascara=(pi*b6_mascara)/(206.25*cos_zwn*E0)

  b6_mascara[b6_mascara <= 0] = 0

  b7_mascara=(pi*b7_mascara)/(68.75*cos_zwn*E0)

  b7_mascara[b7_mascara <= 0] = 0

  Alb_Top = b1_mascara*0.1+b2_mascara*0.31+b3_mascara*0.30+b4_mascara*0.13+b5_mascara*0.08+b6_mascara*0.05+b7_mascara*0.04

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149


  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)

  NDVI =(b5_mascara-b4_mascara)/(b5_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", filetype = "GTiff", overwrite=TRUE)

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b7_mascara, method="bilinear")

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*Ta-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", filetype = "GTiff", overwrite=TRUE)

  slope =(4098*(0.6108*exp((17.27*(Ta))/((Ta)+237.3)))/((Ta)+237.3)^2)

  LEeq = (slope*Rn)/(slope+0.066)

  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, slope, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

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
#' Actual evapotranspiration (ETa) using Landsat-8 images with a grid of agrometeorological data.
#' @param doy  is the Day of Year (DOY)
#' @param a  is one of the regression coefficients of SAFER algorithm
#' @param b is one of the regression coefficients of SAFER algorithm
#' @export
#' @import terra
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), net radiation ("Rn_MJ"), Crop Coefficient ("kc") and Actual Evapotranspiration (evapo).


evapo_l8_grid  = function(doy, a, b){

  b1 <- rast("B1.tif")
  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b5 <- rast("B5.tif")
  b6 <- rast("B6.tif")
  b7 <- rast("B7.tif")
  RG <- rast("RG.tif")
  Ta <- rast("Ta.tif")
  ET0 <- rast("ET0.tif")

  mask <- vect("mask.shp")

  b1_crop <- crop(b1, ext(mask)[1:4])
  b1_mascara <- mask(b1_crop, mask)
  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b5_crop <- crop(b5, ext(mask)[1:4])
  b5_mascara <- mask(b5_crop, mask)
  b6_crop <- crop(b6, ext(mask)[1:4])
  b6_mascara <- mask(b6_crop, mask)
  b7_crop <- crop(b7, ext(mask)[1:4])
  b7_mascara <- mask(b7_crop, mask)
  RG_crop <- crop(RG, ext(mask)[1:4])
  RG <- mask(RG_crop, mask)
  Ta_crop <- crop(Ta, ext(mask)[1:4])
  Ta <- mask(Ta_crop, mask)
  ET0_crop <- crop(ET0, ext(mask)[1:4])
  ET0 <- mask(ET0_crop, mask)

  metadata <- list.files(pattern = "txt")
  m <- read.csv(metadata, header = T)
  metadata <- paste( unlist(m), collapse='')
  rm(m)

  RADIANCE_MAXIMUM_BAND_1 <- "^.*RADIANCE_MAXIMUM_BAND_1 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_1 <- gsub(RADIANCE_MAXIMUM_BAND_1, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_1 <- as.numeric(RADIANCE_MAXIMUM_BAND_1)

  RADIANCE_MAXIMUM_BAND_2 <- "^.*RADIANCE_MAXIMUM_BAND_2 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_2 <- gsub(RADIANCE_MAXIMUM_BAND_2, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_2 <- as.numeric(RADIANCE_MAXIMUM_BAND_2)

  RADIANCE_MAXIMUM_BAND_3 <- "^.*RADIANCE_MAXIMUM_BAND_3 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_3 <- gsub(RADIANCE_MAXIMUM_BAND_3, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_3 <- as.numeric(RADIANCE_MAXIMUM_BAND_3)

  RADIANCE_MAXIMUM_BAND_4 <- "^.*RADIANCE_MAXIMUM_BAND_4 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_4 <- gsub(RADIANCE_MAXIMUM_BAND_4, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_4 <- as.numeric(RADIANCE_MAXIMUM_BAND_4)

  RADIANCE_MAXIMUM_BAND_5 <- "^.*RADIANCE_MAXIMUM_BAND_5 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_5 <- gsub(RADIANCE_MAXIMUM_BAND_5, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_5 <- as.numeric(RADIANCE_MAXIMUM_BAND_5)

  RADIANCE_MAXIMUM_BAND_6 <- "^.*RADIANCE_MAXIMUM_BAND_6 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_6 <- gsub(RADIANCE_MAXIMUM_BAND_6, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_6 <- as.numeric(RADIANCE_MAXIMUM_BAND_6)

  RADIANCE_MAXIMUM_BAND_7 <- "^.*RADIANCE_MAXIMUM_BAND_7 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_7 <- gsub(RADIANCE_MAXIMUM_BAND_7, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_7 <- as.numeric(RADIANCE_MAXIMUM_BAND_7)

  RADIANCE_MINIMUM_BAND_1 <- "^.*RADIANCE_MINIMUM_BAND_1 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_1 <- gsub(RADIANCE_MINIMUM_BAND_1, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_1 <- as.numeric(RADIANCE_MINIMUM_BAND_1)

  RADIANCE_MINIMUM_BAND_2 <- "^.*RADIANCE_MINIMUM_BAND_2 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_2 <- gsub(RADIANCE_MINIMUM_BAND_2, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_2 <- as.numeric(RADIANCE_MINIMUM_BAND_2)

  RADIANCE_MINIMUM_BAND_3 <- "^.*RADIANCE_MINIMUM_BAND_3 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_3 <- gsub(RADIANCE_MINIMUM_BAND_3, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_3 <- as.numeric(RADIANCE_MINIMUM_BAND_3)

  RADIANCE_MINIMUM_BAND_4 <- "^.*RADIANCE_MINIMUM_BAND_4 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_4 <- gsub(RADIANCE_MINIMUM_BAND_4, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_4 <- as.numeric(RADIANCE_MINIMUM_BAND_4)

  RADIANCE_MINIMUM_BAND_5 <- "^.*RADIANCE_MINIMUM_BAND_5 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_5 <- gsub(RADIANCE_MINIMUM_BAND_5, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_5 <- as.numeric(RADIANCE_MINIMUM_BAND_5)

  RADIANCE_MINIMUM_BAND_6 <- "^.*RADIANCE_MINIMUM_BAND_6 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_6 <- gsub(RADIANCE_MINIMUM_BAND_6, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_6 <- as.numeric(RADIANCE_MINIMUM_BAND_6)

  RADIANCE_MINIMUM_BAND_7 <- "^.*RADIANCE_MINIMUM_BAND_7 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_7 <- gsub(RADIANCE_MINIMUM_BAND_7, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_7 <- as.numeric(RADIANCE_MINIMUM_BAND_7)

  b1_mascara = ((RADIANCE_MAXIMUM_BAND_1+RADIANCE_MINIMUM_BAND_1)/65535)*b1_mascara-(RADIANCE_MINIMUM_BAND_1)

  b2_mascara = ((RADIANCE_MAXIMUM_BAND_2+RADIANCE_MINIMUM_BAND_2)/65535)*b2_mascara-(RADIANCE_MINIMUM_BAND_2)

  b3_mascara = ((RADIANCE_MAXIMUM_BAND_3+RADIANCE_MINIMUM_BAND_3)/65535)*b3_mascara-(RADIANCE_MINIMUM_BAND_3)

  b4_mascara = ((RADIANCE_MAXIMUM_BAND_4+RADIANCE_MINIMUM_BAND_4)/65535)*b4_mascara-(RADIANCE_MINIMUM_BAND_4)

  b5_mascara = ((RADIANCE_MAXIMUM_BAND_5+RADIANCE_MINIMUM_BAND_5)/65535)*b5_mascara-(RADIANCE_MINIMUM_BAND_5)

  b6_mascara = ((RADIANCE_MAXIMUM_BAND_6+RADIANCE_MINIMUM_BAND_6)/65535)*b6_mascara-(RADIANCE_MINIMUM_BAND_6)

  b7_mascara = ((RADIANCE_MAXIMUM_BAND_7+RADIANCE_MINIMUM_BAND_7)/65535)*b7_mascara-(RADIANCE_MINIMUM_BAND_7)

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

  b1_mascara=(pi*b1_mascara)/(1718.75*cos_zwn*E0)

  b1_mascara[b1_mascara <= 0] = 0

  b2_mascara=(pi*b2_mascara)/(1810.42*cos_zwn*E0)

  b2_mascara[b2_mascara <= 0] = 0

  b3_mascara=(pi*b3_mascara)/(1741.67*cos_zwn*E0)

  b3_mascara[b3_mascara <= 0] = 0

  b4_mascara=(pi*b4_mascara)/(1558.33*cos_zwn*E0)

  b4_mascara[b4_mascara <= 0] = 0

  b5_mascara=(pi*b5_mascara)/(962.5*cos_zwn*E0)

  b5_mascara[b5_mascara <= 0] = 0

  b6_mascara=(pi*b6_mascara)/(206.25*cos_zwn*E0)

  b6_mascara[b6_mascara <= 0] = 0

  b7_mascara=(pi*b7_mascara)/(68.75*cos_zwn*E0)

  b7_mascara[b7_mascara <= 0] = 0

  Alb_Top = b1_mascara*0.1+b2_mascara*0.31+b3_mascara*0.30+b4_mascara*0.13+b5_mascara*0.08+b6_mascara*0.05+b7_mascara*0.04

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149


  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)

  NDVI =(b5_mascara-b4_mascara)/(b5_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", filetype = "GTiff", overwrite=TRUE)

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b7_mascara, method="bilinear")

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*Ta-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", filetype = "GTiff", overwrite=TRUE)

  slope =(4098*(0.6108*exp((17.27*(Ta))/((Ta)+237.3)))/((Ta)+237.3)^2)

  LEeq = (slope*Rn)/(slope+0.066)

  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, slope, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

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

#'Energy balance using Landsat-8 images with a grid of agrometeorological data.
#'@param doy is the Day of Year (DOY)
#'@param a is one of the regression coefficients of SAFER algorithm
#'@param b is one of the regression coefficients of SAFER algorithm
#'@export
#' @import terra
#' @importFrom utils read.csv
#'
#'@return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), Crop Coefficient ("kc"), Actual Evapotranspiration (evapo), latent heat flux "LE_MJ"), net radiation ("Rn_MJ"), ground heat flux ("G_MJ") and the sensible heat flux ("H_MJ").


radiation_l8_grid  =  function(doy, a, b){


  b1 <- rast("B1.tif")
  b2 <- rast("B2.tif")
  b3 <- rast("B3.tif")
  b4 <- rast("B4.tif")
  b5 <- rast("B5.tif")
  b6 <- rast("B6.tif")
  b7 <- rast("B7.tif")
  RG <- rast("RG.tif")
  Ta <- rast("Ta.tif")
  ET0 <- rast("ET0.tif")

  mask <- vect("mask.shp")

  b1_crop <- crop(b1, ext(mask)[1:4])
  b1_mascara <- mask(b1_crop, mask)
  b2_crop <- crop(b2, ext(mask)[1:4])
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, ext(mask)[1:4])
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, ext(mask)[1:4])
  b4_mascara <- mask(b4_crop, mask)
  b5_crop <- crop(b5, ext(mask)[1:4])
  b5_mascara <- mask(b5_crop, mask)
  b6_crop <- crop(b6, ext(mask)[1:4])
  b6_mascara <- mask(b6_crop, mask)
  b7_crop <- crop(b7, ext(mask)[1:4])
  b7_mascara <- mask(b7_crop, mask)
  RG_crop <- crop(RG, ext(mask)[1:4])
  RG <- mask(RG_crop, mask)
  Ta_crop <- crop(Ta, ext(mask)[1:4])
  Ta <- mask(Ta_crop, mask)
  ET0_crop <- crop(ET0, ext(mask)[1:4])
  ET0 <- mask(ET0_crop, mask)

  metadata <- list.files(pattern = "txt")
  m <- read.csv(metadata, header = T)
  metadata <- paste( unlist(m), collapse='')
  rm(m)

  RADIANCE_MAXIMUM_BAND_1 <- "^.*RADIANCE_MAXIMUM_BAND_1 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_1 <- gsub(RADIANCE_MAXIMUM_BAND_1, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_1 <- as.numeric(RADIANCE_MAXIMUM_BAND_1)

  RADIANCE_MAXIMUM_BAND_2 <- "^.*RADIANCE_MAXIMUM_BAND_2 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_2 <- gsub(RADIANCE_MAXIMUM_BAND_2, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_2 <- as.numeric(RADIANCE_MAXIMUM_BAND_2)

  RADIANCE_MAXIMUM_BAND_3 <- "^.*RADIANCE_MAXIMUM_BAND_3 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_3 <- gsub(RADIANCE_MAXIMUM_BAND_3, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_3 <- as.numeric(RADIANCE_MAXIMUM_BAND_3)

  RADIANCE_MAXIMUM_BAND_4 <- "^.*RADIANCE_MAXIMUM_BAND_4 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_4 <- gsub(RADIANCE_MAXIMUM_BAND_4, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_4 <- as.numeric(RADIANCE_MAXIMUM_BAND_4)

  RADIANCE_MAXIMUM_BAND_5 <- "^.*RADIANCE_MAXIMUM_BAND_5 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_5 <- gsub(RADIANCE_MAXIMUM_BAND_5, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_5 <- as.numeric(RADIANCE_MAXIMUM_BAND_5)

  RADIANCE_MAXIMUM_BAND_6 <- "^.*RADIANCE_MAXIMUM_BAND_6 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_6 <- gsub(RADIANCE_MAXIMUM_BAND_6, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_6 <- as.numeric(RADIANCE_MAXIMUM_BAND_6)

  RADIANCE_MAXIMUM_BAND_7 <- "^.*RADIANCE_MAXIMUM_BAND_7 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MAXIMUM_BAND_7 <- gsub(RADIANCE_MAXIMUM_BAND_7, "\\1", metadata)
  RADIANCE_MAXIMUM_BAND_7 <- as.numeric(RADIANCE_MAXIMUM_BAND_7)

  RADIANCE_MINIMUM_BAND_1 <- "^.*RADIANCE_MINIMUM_BAND_1 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_1 <- gsub(RADIANCE_MINIMUM_BAND_1, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_1 <- as.numeric(RADIANCE_MINIMUM_BAND_1)

  RADIANCE_MINIMUM_BAND_2 <- "^.*RADIANCE_MINIMUM_BAND_2 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_2 <- gsub(RADIANCE_MINIMUM_BAND_2, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_2 <- as.numeric(RADIANCE_MINIMUM_BAND_2)

  RADIANCE_MINIMUM_BAND_3 <- "^.*RADIANCE_MINIMUM_BAND_3 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_3 <- gsub(RADIANCE_MINIMUM_BAND_3, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_3 <- as.numeric(RADIANCE_MINIMUM_BAND_3)

  RADIANCE_MINIMUM_BAND_4 <- "^.*RADIANCE_MINIMUM_BAND_4 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_4 <- gsub(RADIANCE_MINIMUM_BAND_4, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_4 <- as.numeric(RADIANCE_MINIMUM_BAND_4)

  RADIANCE_MINIMUM_BAND_5 <- "^.*RADIANCE_MINIMUM_BAND_5 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_5 <- gsub(RADIANCE_MINIMUM_BAND_5, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_5 <- as.numeric(RADIANCE_MINIMUM_BAND_5)

  RADIANCE_MINIMUM_BAND_6 <- "^.*RADIANCE_MINIMUM_BAND_6 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_6 <- gsub(RADIANCE_MINIMUM_BAND_6, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_6 <- as.numeric(RADIANCE_MINIMUM_BAND_6)

  RADIANCE_MINIMUM_BAND_7 <- "^.*RADIANCE_MINIMUM_BAND_7 = *?[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+).*"
  RADIANCE_MINIMUM_BAND_7 <- gsub(RADIANCE_MINIMUM_BAND_7, "\\1", metadata)
  RADIANCE_MINIMUM_BAND_7 <- as.numeric(RADIANCE_MINIMUM_BAND_7)

  b1_mascara = ((RADIANCE_MAXIMUM_BAND_1+RADIANCE_MINIMUM_BAND_1)/65535)*b1_mascara-(RADIANCE_MINIMUM_BAND_1)

  b2_mascara = ((RADIANCE_MAXIMUM_BAND_2+RADIANCE_MINIMUM_BAND_2)/65535)*b2_mascara-(RADIANCE_MINIMUM_BAND_2)

  b3_mascara = ((RADIANCE_MAXIMUM_BAND_3+RADIANCE_MINIMUM_BAND_3)/65535)*b3_mascara-(RADIANCE_MINIMUM_BAND_3)

  b4_mascara = ((RADIANCE_MAXIMUM_BAND_4+RADIANCE_MINIMUM_BAND_4)/65535)*b4_mascara-(RADIANCE_MINIMUM_BAND_4)

  b5_mascara = ((RADIANCE_MAXIMUM_BAND_5+RADIANCE_MINIMUM_BAND_5)/65535)*b5_mascara-(RADIANCE_MINIMUM_BAND_5)

  b6_mascara = ((RADIANCE_MAXIMUM_BAND_6+RADIANCE_MINIMUM_BAND_6)/65535)*b6_mascara-(RADIANCE_MINIMUM_BAND_6)

  b7_mascara = ((RADIANCE_MAXIMUM_BAND_7+RADIANCE_MINIMUM_BAND_7)/65535)*b7_mascara-(RADIANCE_MINIMUM_BAND_7)


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

  b1_mascara=(pi*b1_mascara)/(1718.75*cos_zwn*E0)

  b1_mascara[b1_mascara <= 0] = 0

  b2_mascara=(pi*b2_mascara)/(1810.42*cos_zwn*E0)

  b2_mascara[b2_mascara <= 0] = 0

  b3_mascara=(pi*b3_mascara)/(1741.67*cos_zwn*E0)

  b3_mascara[b3_mascara <= 0] = 0

  b4_mascara=(pi*b4_mascara)/(1558.33*cos_zwn*E0)

  b4_mascara[b4_mascara <= 0] = 0

  b5_mascara=(pi*b5_mascara)/(962.5*cos_zwn*E0)

  b5_mascara[b5_mascara <= 0] = 0

  b6_mascara=(pi*b6_mascara)/(206.25*cos_zwn*E0)

  b6_mascara[b6_mascara <= 0] = 0

  b7_mascara=(pi*b7_mascara)/(68.75*cos_zwn*E0)

  b7_mascara[b7_mascara <= 0] = 0

  Alb_Top = b1_mascara*0.1+b2_mascara*0.31+b3_mascara*0.30+b4_mascara*0.13+b5_mascara*0.08+b6_mascara*0.05+b7_mascara*0.04

  Alb_sur = 0.6054*Alb_Top + 0.0797

  Alb_24 =  1.0223*Alb_sur + 0.0149


  writeRaster(Alb_24, "Alb_24", filetype = "GTiff", overwrite=TRUE)

  NDVI =(b5_mascara-b4_mascara)/(b5_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", filetype = "GTiff", overwrite=TRUE)

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b7_mascara, method="bilinear")

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*Ta-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", filetype = "GTiff", overwrite=TRUE)

  slope =(4098*(0.6108*exp((17.27*(Ta))/((Ta)+237.3)))/((Ta)+237.3)^2)

  LEeq = (slope*Rn)/(slope+0.066)

  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, slope, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

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
