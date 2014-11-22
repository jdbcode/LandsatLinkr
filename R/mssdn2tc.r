
#' Transform original MSS bands to Tasseled Cap
#'
#' Transform original MSS bands to Tasseled Cap brightness, greenness, and angle (angle between brightness and greenness). TC yellowness is currently excluded.
#' @param file The full path name of the *archv file
#' @param Parameter indicating which outputs to provide: !list options!
#' @import raster
#' @export
#' @details Transformation based on coefficents and equation presented in: Pflugmacher, D., Cohen, W. B., & E. Kennedy, R. (2012). Using Landsat-derived disturbance history (1972-2010) to predict current forest structure. Remote Sensing of Environment, 122, 146-165. doi:10.1016/j.rse.2011.09.025


mssdn2tc = function(file, tcset){
  
  #1 = angle
  #2 = brightness, greenness, angle
  #3 = brightness, greenness
  #4 = brightness, greenness, yellowness
  #5 = brightness, greenness, yellowness, angle
  
  b = brick(file)
  img = as.array(b)
  
  sensor = substr(basename(file), 3, 3)
  
  #coefficients
  if (sensor == "1") {
    bcoef = c(0.2031, 0.352, 0.4627, 0.0848, -12.7421) #bias = -12.7421
    gcoef = c(-0.1731, -0.3852, 0.3953, 0.1253, -9.5016) #bias = -9.5016
    if(tcset >= 4){ycoef = c(-0.5497, 0.2499, 0.052, -0.0132, 5.9304)} #bias = 5.9304
  }
  
  if (sensor == "2") {
    bcoef = c(0.2334, 0.3679, 0.3975, 0.065, -11.6690) #bias = -11.6690
    gcoef = c(-0.1989, -0.4026, 0.3396, 0.096, -0.6579) #bias = -0.6579
    if(tcset >= 4){ycoef = c(-0.6318, 0.2612, 0.0446, -0.0101, -2.4038)} #bias = -2.4038
  }
  
  if (sensor == "3") {
    bcoef = c(0.2336, 0.3504, 0.39, 0.0642, 7.3048) #bias = 7.3048
    gcoef = c(-0.1991, -0.3835, 0.3332, 0.0948, 0.6085) #bias = 0.6085
    if(tcset >= 4){ycoef = c(-0.6323, 0.2488, 0.0438, -0.0100, 1.9141)} #bias = 1.9141
  }
  
  if (sensor == "4") {
    bcoef = c(0.2153, 0.3774, 0.3418, 0.0682 ,6.0451) #bias = 6.0451
    gcoef = c(-0.1835, -0.4130, 0.292, 0.1007, 1.3184) #bias = 1.3184
    if(tcset >= 4){ycoef = c(-0.5829, 0.268, 0.0384, -0.0106, 3.5011)} #bias = 3.5011
  }
  
  if (sensor == "5") {
    bcoef = c(0.2208, 0.3492, 0.3688, 0.0611, 6.0382) #bias = 6.0382
    gcoef = c(-0.1881, -0.3822, 0.3151, 0.0903, 2.2989) #bias = 2.2989
    if(tcset >= 4){ycoef = c(-0.5976, 0.248, 0.0414, -0.0095, 2.7263)} #bias = 2.7263
  }
    
  bright =                (((img[,,1]*bcoef[1])+(img[,,2]*bcoef[2])+(img[,,3]*bcoef[3])+(img[,,4]*bcoef[4])) + bcoef[5]) *100 #+ bcoef[5]
  green =                 (((img[,,1]*gcoef[1])+(img[,,2]*gcoef[2])+(img[,,3]*gcoef[3])+(img[,,4]*gcoef[4])) + gcoef[5]) *100 #+ gcoef[5]
  if(tcset >= 4){yellow = (((img[,,1]*ycoef[1])+(img[,,2]*ycoef[2])+(img[,,3]*ycoef[3])+(img[,,4]*ycoef[4])) + ycoef[5]) *100} #+ ycoef[5]
  
  tcfile = sub("archv.tif","tc.tif", file)
  tcafile = sub("archv.tif","tca.tif", file)
  
  #calc tc and convert to a raster
  if(tcset >= 2){
    tcb = setValues(raster(file),bright)
    tcg = setValues(raster(file),green)
    if(tcset >= 4){
      tcy = setValues(raster(file),yellow)
      tc = stack(tcb,tcg,tcy)
    } else {tc = stack(tcb,tcg)}
    projection(tc) = set_projection(file)
    tc = as(tc, "SpatialGridDataFrame")
    writeGDAL(tc, tcfile, drivername = "GTiff", type = "Int16", mvFlag = -32768)
  }
  
  #calc tc angle and convert to a raster
  if(tcset <= 2 | tcset == 5){
    tca = atan(green/bright) * (180/pi) * 100 #need to multiply by (180/pi) to get degrees because atan() returns radians, 100 is a scalar to preserve two decimal places
    tca = setValues(raster(file),tca)
    projection(tca) = set_projection(file)
    tca = as(tca, "SpatialGridDataFrame")
    writeGDAL(tca, tcafile, drivername = "GTiff", type = "Int16", mvFlag = -32768)
  }
}