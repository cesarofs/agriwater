#SINGLE VALUES OF AGROMETEOROLOGICAL DATA
#LANDSAT-8 com thermal

#'Crop coefficient (ETa / ET0) using Landsat-8 images (including thermal bands) with single agrometeorological data.
#'
#' @param doy   is the Day of Year (DOY)
#' @param RG  is the global solar radiation
#' @param Ta  is the average air temperature
#' @param a  is one of the regression coefficients of SAFER algorithm
#' @param b  is one of the regression coefficients of SAFER algorithm
#' @export
#' @import raster
#' @import sp
#' @import rgdal
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), Crop Coefficient ("kc") and net radiation ("Rn_MJ").

kc_l8t = function(doy, RG, Ta, a, b){

  b1 <- raster("B1.tif")
  b2 <- raster("B2.tif")
  b3 <- raster("B3.tif")
  b4 <- raster("B4.tif")
  b5 <- raster("B5.tif")
  b6 <- raster("B6.tif")
  b7 <- raster("B7.tif")
  b10 <- raster("B10.tif")
  b11 <- raster("B11.tif")

  mask <- readOGR("mask.shp")

  b1_crop <- crop(b1, extent(mask))
  b1_mascara <- mask(b1_crop, mask)
  b2_crop <- crop(b2, extent(mask))
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, extent(mask))
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, extent(mask))
  b4_mascara <- mask(b4_crop, mask)
  b5_crop <- crop(b5, extent(mask))
  b5_mascara <- mask(b5_crop, mask)
  b6_crop <- crop(b6, extent(mask))
  b6_mascara <- mask(b6_crop, mask)
  b7_crop <- crop(b7, extent(mask))
  b7_mascara <- mask(b7_crop, mask)
  b10_crop <- crop(b10, extent(mask))
  b10_mascara <- mask(b10_crop, mask)
  b11_crop <- crop(b11, extent(mask))
  b11_mascara <- mask(b11_crop, mask)

  metadata <- list.files(pattern = "txt")
  m <- read.csv(metadata, header = T)

  radiance = m[c(82:103),]
  radiance = substr(radiance, 31,41)
  radiance = as.numeric(radiance)

  b1_mascara = ((radiance[1]


                 -radiance[2])/65535)*b1_mascara+(radiance[2])

  b2_mascara = ((radiance[3]-radiance[4])/65535)*b2_mascara+(radiance[4])

  b3_mascara = ((radiance[5]-radiance[6]
  )/65535)*b3_mascara+(radiance[6])

  b4_mascara = ((radiance[7]
                 -radiance[8])/65535)*b4_mascara+(radiance[8])

  b5_mascara = ((radiance[9]-radiance[10])/65535)*b5_mascara+(radiance[10])

  b6_mascara = ((radiance[11]-radiance[12])/65535)*b6_mascara+(radiance[12]
  )

  b7_mascara = ((radiance[13]-radiance[14])/65535)*b7_mascara+(radiance[14]
  )

  thermal = m[c(82:103),]
  thermal = substr(thermal, 32,42)
  thermal = as.numeric(thermal)

  b10_mascara = ((thermal[19]-thermal[20])/65535)*b10_mascara+(thermal[20]
  )

  b11_mascara = ((thermal[21]-thermal[22])/65535)*b10_mascara+(thermal[22]
  )


  lati <- long <- b2_mascara
  xy <- coordinates(b2_mascara)
  long[] <- xy[, 1]
  long <- crop(long, extent(mask))
  lati[] <- xy[, 2]
  lati <- crop(lati, extent(mask))

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


  writeRaster(Alb_24, "Alb_24", format = "GTiff", overwrite=TRUE)

  NDVI =(b5_mascara-b4_mascara)/(b5_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", format = "GTiff", overwrite=TRUE)

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b7_mascara, method="bilinear")



  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara,  R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

  thermal = m[c(192:195),]
  thermal = substr(thermal, 27, 34)
  thermal = as.numeric(thermal)

  Tbright_10 = thermal[3]/log((thermal[1]/(b10_mascara+1)))

  Tbright_11 = thermal[4]/log((thermal[2]/(b11_mascara+1)))

  Tbright = ((Tbright_10+Tbright_11)/2)

  TS24 = 1.0694*Tbright-20.173

  writeRaster(TS24, "LST", format = "GTiff", overwrite=TRUE)

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*(Ta-273.15)-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", format = "GTiff", overwrite=TRUE)

  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)


  NDVI[NDVI <= 0] = NA

  kc=exp((a)+(b*((TS24-273.15)/(Alb_sur*NDVI))))

  writeRaster(kc, "kc", format = "GTiff", overwrite=TRUE)
}


#'Actual evapotranspiration (ETa) using Landsat-8 (including thermal bands) images with single agrometeorological data.
#' @param doy  is the Day of Year (DOY)
#' @param RG is the global solar radiation
#' @param Ta  is the average air temperature
#' @param ET0 is the reference evapotranspiration
#' @param a  is one of the regression coefficients of SAFER algorithm
#' @param b is one of the regression coefficients of SAFER algorithm
#' @export
#' @import raster
#' @import sp
#' @import rgdal
#' @importFrom utils read.csv
#'
#' @return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), net radiation ("Rn_MJ"), Crop Coefficient ("kc") and Actual Evapotranspiration (evapo).

evapo_l8t = function(doy, RG, Ta, ET0, a, b){

  b1 <- raster("B1.tif")
  b2 <- raster("B2.tif")
  b3 <- raster("B3.tif")
  b4 <- raster("B4.tif")
  b5 <- raster("B5.tif")
  b6 <- raster("B6.tif")
  b7 <- raster("B7.tif")
  b10 <- raster("B10.tif")
  b11 <- raster("B11.tif")

  mask <- readOGR("mask.shp")

  b1_crop <- crop(b1, extent(mask))
  b1_mascara <- mask(b1_crop, mask)
  b2_crop <- crop(b2, extent(mask))
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, extent(mask))
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, extent(mask))
  b4_mascara <- mask(b4_crop, mask)
  b5_crop <- crop(b5, extent(mask))
  b5_mascara <- mask(b5_crop, mask)
  b6_crop <- crop(b6, extent(mask))
  b6_mascara <- mask(b6_crop, mask)
  b7_crop <- crop(b7, extent(mask))
  b7_mascara <- mask(b7_crop, mask)
  b10_crop <- crop(b10, extent(mask))
  b10_mascara <- mask(b10_crop, mask)
  b11_crop <- crop(b11, extent(mask))
  b11_mascara <- mask(b11_crop, mask)

  metadata <- list.files(pattern = "txt")
  m <- read.csv(metadata, header = T)

  radiance = m[c(82:103),]
  radiance = substr(radiance, 31,41)
  radiance = as.numeric(radiance)

  b1_mascara = ((radiance[1]


                 -radiance[2])/65535)*b1_mascara+(radiance[2])

  b2_mascara = ((radiance[3]-radiance[4])/65535)*b2_mascara+(radiance[4])

  b3_mascara = ((radiance[5]-radiance[6]
  )/65535)*b3_mascara+(radiance[6])

  b4_mascara = ((radiance[7]
                 -radiance[8])/65535)*b4_mascara+(radiance[8])

  b5_mascara = ((radiance[9]-radiance[10])/65535)*b5_mascara+(radiance[10])

  b6_mascara = ((radiance[11]-radiance[12])/65535)*b6_mascara+(radiance[12]
  )

  b7_mascara = ((radiance[13]-radiance[14])/65535)*b7_mascara+(radiance[14]
  )

  thermal = m[c(82:103),]
  thermal = substr(thermal, 32,42)
  thermal = as.numeric(thermal)

  b10_mascara = ((thermal[19]-thermal[20])/65535)*b10_mascara+(thermal[20]
  )

  b11_mascara = ((thermal[21]-thermal[22])/65535)*b11_mascara+(thermal[22]
  )


  lati <- long <- b2_mascara
  xy <- coordinates(b2_mascara)
  long[] <- xy[, 1]
  long <- crop(long, extent(mask))
  lati[] <- xy[, 2]
  lati <- crop(lati, extent(mask))

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


  writeRaster(Alb_24, "Alb_24", format = "GTiff", overwrite=TRUE)

  NDVI =(b5_mascara-b4_mascara)/(b5_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", format = "GTiff", overwrite=TRUE)

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b7_mascara, method="bilinear")




  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, Rn_coeff, RsTOP, RsTOP_aux, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)

  thermal = m[c(192:195),]
  thermal = substr(thermal, 27, 34)
  thermal = as.numeric(thermal)

  Tbright_10 = thermal[3]/log((thermal[1]/(b10_mascara+1)))

  Tbright_11 = thermal[4]/log((thermal[2]/(b11_mascara+1)))

  Tbright = ((Tbright_10+Tbright_11)/2)

  TS24 = 1.0694*Tbright-20.173

  writeRaster(TS24, "LST", format = "GTiff", overwrite=TRUE)

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*(Ta-273.15)-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", format = "GTiff", overwrite=TRUE)

  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)


  NDVI[NDVI <= 0] = NA

  kc=exp((a)+(b*((TS24-273.15)/(Alb_sur*NDVI))))

  writeRaster(kc, "kc", format = "GTiff", overwrite=TRUE)

  ET=kc*ET0

  writeRaster(ET, "evapo", format = "GTiff", overwrite=TRUE)
}


#'Energy balance using Landsat-8 images (including thermal bands) with single agrometeorological data.
#'@param doy is the Day of Year (DOY)
#'@param RG is the global solar radiation
#'@param Ta  is the average air temperature
#'@param ET0  is the reference evapotranspiration
#'@param a is one of the regression coefficients of SAFER algorithm
#'@param b is one of the regression coefficients of SAFER algorithm
#'@export
#'@import raster
#' @import sp
#' @import rgdal
#' @importFrom utils read.csv
#'
#'@return It returns in raster format (.tif) the Surface Albedo at 24h scale ("Alb_24"), NDVI, Surface Temperature ("LST"), Crop Coefficient ("kc"), Actual Evapotranspiration (evapo), latent heat flux "LE_MJ"), net radiation ("Rn_MJ"), ground heat flux ("G_MJ") and the sensible heat flux ("H_MJ").


radiation_l8t =  function(doy, RG, Ta, ET0, a, b){


  b1 <- raster("B1.tif")
  b2 <- raster("B2.tif")
  b3 <- raster("B3.tif")
  b4 <- raster("B4.tif")
  b5 <- raster("B5.tif")
  b6 <- raster("B6.tif")
  b7 <- raster("B7.tif")
  b10 <- raster("B10.tif")
  b11 <- raster("B11.tif")

  mask <- readOGR("mask.shp")

  b1_crop <- crop(b1, extent(mask))
  b1_mascara <- mask(b1_crop, mask)
  b2_crop <- crop(b2, extent(mask))
  b2_mascara <- mask(b2_crop, mask)
  b3_crop <- crop(b3, extent(mask))
  b3_mascara <- mask(b3_crop, mask)
  b4_crop <- crop(b4, extent(mask))
  b4_mascara <- mask(b4_crop, mask)
  b5_crop <- crop(b5, extent(mask))
  b5_mascara <- mask(b5_crop, mask)
  b6_crop <- crop(b6, extent(mask))
  b6_mascara <- mask(b6_crop, mask)
  b7_crop <- crop(b7, extent(mask))
  b7_mascara <- mask(b7_crop, mask)
  b10_crop <- crop(b10, extent(mask))
  b10_mascara <- mask(b10_crop, mask)
  b11_crop <- crop(b11, extent(mask))
  b11_mascara <- mask(b11_crop, mask)

  metadata <- list.files(pattern = "txt")
  m <- read.csv(metadata, header = T)

  radiance = m[c(82:103),]
  radiance = substr(radiance, 31,41)
  radiance = as.numeric(radiance)

  b1_mascara = ((radiance[1]


                 -radiance[2])/65535)*b1_mascara+(radiance[2])

  b2_mascara = ((radiance[3]-radiance[4])/65535)*b2_mascara+(radiance[4])

  b3_mascara = ((radiance[5]-radiance[6]
  )/65535)*b3_mascara+(radiance[6])

  b4_mascara = ((radiance[7]
                 -radiance[8])/65535)*b4_mascara+(radiance[8])

  b5_mascara = ((radiance[9]-radiance[10])/65535)*b5_mascara+(radiance[10])

  b6_mascara = ((radiance[11]-radiance[12])/65535)*b6_mascara+(radiance[12]
  )

  b7_mascara = ((radiance[13]-radiance[14])/65535)*b7_mascara+(radiance[14]
  )

  thermal = m[c(82:103),]
  thermal = substr(thermal, 32,42)
  thermal = as.numeric(thermal)

  b10_mascara = ((thermal[19]-thermal[20])/65535)*b10_mascara+(thermal[20]
  )

  b11_mascara = ((thermal[21]-thermal[22])/65535)*b11_mascara+(thermal[22]
  )


  lati <- long <- b2_mascara
  xy <- coordinates(b2_mascara)
  long[] <- xy[, 1]
  long <- crop(long, extent(mask))
  lati[] <- xy[, 2]
  lati <- crop(lati, extent(mask))

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


  writeRaster(Alb_24, "Alb_24", format = "GTiff", overwrite=TRUE)

  NDVI =(b5_mascara-b4_mascara)/(b5_mascara+b4_mascara)

  writeRaster(NDVI, "NDVI", format = "GTiff", overwrite=TRUE)

  Ws = acos(((-1)*tan(lati*pi/180))*tan(Dec))

  R =(Ws*sin(lati*pi/180)*sin(Dec))+(cos(lati*pi/180)*cos(Dec)*sin(Ws))

  RsTOP_aux =(1367/pi)*E0*R

  RsTOP = resample(RsTOP_aux, b7_mascara, method="bilinear")


  thermal = m[c(192:195),]
  thermal = substr(thermal, 27, 34)
  thermal = as.numeric(thermal)

  Tbright_10 = thermal[3]/log((thermal[1]/(b10_mascara+1)))

  Tbright_11 = thermal[4]/log((thermal[2]/(b11_mascara+1)))

  Tbright = ((Tbright_10+Tbright_11)/2)

  TS24 = 1.0694*Tbright-20.173

  Transm =(RG*11.6)/RsTOP

  Rn_coeff =6.99*(Ta-273.15)-39.99

  Rn =((1-Alb_24)*(RG*11.6))-(Rn_coeff*Transm)

  Rn_MJ =Rn/11.6

  writeRaster(Rn_MJ, "Rn_MJ", format = "GTiff", overwrite=TRUE)

  rm(b1, b2, b3, b4, b5, b6, b7, b1_mascara, b2_mascara, b3_mascara, b4_mascara, b5_mascara,b6_mascara,b7_mascara, R, Ws, E0, cos_zwn, W, Dec, LAT, Et, map1, lati, long)


  writeRaster(TS24, "LST", format = "GTiff", overwrite=TRUE)

  NDVI[NDVI <= 0] = NA

  kc=exp((a)+(b*((TS24-273.15)/(Alb_sur*NDVI))))

  ET=kc*ET0

  LE_MJ =ET*2.45

  writeRaster(LE_MJ, "LE_MJ", format = "GTiff", overwrite=TRUE)

  G_Rn =3.98*exp(-25.47*Alb_24)

  G_MJ =G_Rn*Rn_MJ

  writeRaster(G_MJ, "G_MJ", format = "GTiff", overwrite=TRUE)

  H_MJ =Rn_MJ-LE_MJ-G_MJ

  writeRaster(H_MJ, "H_MJ", format = "GTiff", overwrite=TRUE)

}

