#' Convert DN values to toa reflectance
#'
#' Convert DN values to toa reflectance
#' @param file The full path name of the *archv file
#' @import raster
#' @import rgdal
#' @export


mssdn2refl = function(file){
  
  #link to the equations to convert DN to TOA and TOA to SR
  #http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
  
  #read in the image metadata  
  info = get_metadata(file)
  
  #define esun values for mss (chander et al 2009 summary of current radiometric calibration coefficients... RSE 113)
  if(info$sensor == "LANDSAT_1"){esun = c(1823,1559,1276,880.1)}
  if(info$sensor == "LANDSAT_2"){esun = c(1829,1539,1268,886.6)}
  if(info$sensor == "LANDSAT_3"){esun = c(1839,1555,1291,887.9)}
  if(info$sensor == "LANDSAT_4"){esun = c(1827,1569,1260,866.4)}
  if(info$sensor == "LANDSAT_5"){esun = c(1824,1570,1249,853.4)}
  
  #define the earth sun distance 
  #d = 1.01497  #Earth-Sun Distance in Astronomical Units at DOY 213
  doy = as.numeric(substr(basename(file),14,16))
  d = eudist(doy)
  
  #read in the DN image data  
  b = brick(file)
  img = as.array(b) 
  
  img[,,1] = ((info$b1gain*img[,,1])+info$b1bias)
  img[,,2] = ((info$b2gain*img[,,2])+info$b2bias)
  img[,,3] = ((info$b3gain*img[,,3])+info$b3bias)
  img[,,4] = ((info$b4gain*img[,,4])+info$b4bias)
  
  img[img<0] = 0
  
  #convert solar zenith angle to radians for cos calculation (r expects radians)
  sunzenith = info$sunzen/57.2958 
  
  #calulate SR from DOS TOA
  img[,,1] = (pi * img[,,1] * (d^2))/(esun[1] * cos(sunzenith))
  img[,,2] = (pi * img[,,2] * (d^2))/(esun[2] * cos(sunzenith))
  img[,,3] = (pi * img[,,3] * (d^2))/(esun[3] * cos(sunzenith))
  img[,,4] = (pi * img[,,4] * (d^2))/(esun[4] * cos(sunzenith))
  
  #scale
  img = round(img * 10000)
  
  img = setValues(b,img)
  dataType(img) = "INT2S"
  projection(img) = set_projection(file)
  img = as(img, "SpatialGridDataFrame")         #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  outfile = sub("archv", "reflectance", file)
  writeGDAL(img, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
}