#' Convert DN values to surface reflectance
#'
#' Convert DN values to surface reflectance using the COST model with dark object subtraction
#' @param file The full path name of the *archv file
#' @param drkobjv dark object values from mssdofinder
#' @import raster
#' @export


msscost = function(file, drkobjv){
  
  #link to the equations to convert DN to TOA and TOA to SR
  #http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
  
  refl = function(file, band, gain, bias, sunzen, d, esun, drkobjv){
    #     file = "E:/llr_test/30m/mss/wrs1/036032/images/1973/LM10360321973191_archv.tif"
    #     band =1
    #     gain = info$b1gain
    #     bias = info$b1bias
    #     esun = esun[1]
    
    orig = raster(file, band)
    img = as.matrix(orig)   
    img = ((gain*img)+bias)-drkobjv
    img[img < 0] = 0
    img = (pi * img * (d^2))/(esun * cos(sunzenith))
    img = round(img * 10000)
    img = setValues(orig,img)
    return(img)
    #     dataType(img) = "INT2S"
    #     projection(img) = set_projection(file)
    #     img = as(img, "SpatialGridDataFrame")         #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
    #     outfile = sub("archv", paste("refl_temp_",band,sep=""), file)
    #     writeGDAL(img, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  }
  
  
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
#   img = as.array(b) 
#   
#   drkobjv[1] = (info$b1gain*drkobjv[1])+info$b1bias
#   drkobjv[2] = (info$b2gain*drkobjv[2])+info$b2bias
#   drkobjv[3] = (info$b3gain*drkobjv[3])+info$b3bias
#   drkobjv[4] = (info$b4gain*drkobjv[4])+info$b4bias
  
  
  ############################################
  
 
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
  sunzenith = info$sunzen/57.2958 
  
  #apply the conversion to reflectance
  b1 = refl(file,1,info$b1gain, info$b1bias, sunzen, d, esun[1], drkobjv[1]) 
  b2 = refl(file,2,info$b2gain, info$b2bias, sunzen, d, esun[2], drkobjv[2])
  b3 = refl(file,3,info$b3gain, info$b3bias, sunzen, d, esun[3], drkobjv[3])
  b4 = refl(file,4,info$b4gain, info$b4bias, sunzen, d, esun[4], drkobjv[4])
  
  #stack the bands
  img = stack(b1,b2,b3,b4)
  
  
  
  ##############################
#   
#   
#   img[,,1] = ((info$b1gain*img[,,1])+info$b1bias)-drkobjv[1]
#   img[,,2] = ((info$b2gain*img[,,2])+info$b2bias)-drkobjv[2]
#   img[,,3] = ((info$b3gain*img[,,3])+info$b3bias)-drkobjv[3]
#   img[,,4] = ((info$b4gain*img[,,4])+info$b4bias)-drkobjv[4]
#   
#   img[img<0] = 0
#   
#   #convert solar zenith angle to radians for cos calculation (r expects radians)
#   sunzenith = info$sunzen/57.2958 
#   
#   #calulate SR from DOS TOA
#   img[,,1] = (pi * img[,,1] * (d^2))/(esun[1] * cos(sunzenith))
#   img[,,2] = (pi * img[,,2] * (d^2))/(esun[2] * cos(sunzenith))
#   img[,,3] = (pi * img[,,3] * (d^2))/(esun[3] * cos(sunzenith))
#   img[,,4] = (pi * img[,,4] * (d^2))/(esun[4] * cos(sunzenith))
#   
#   #scale
#   img = round(img * 10000)
#   
#   img = setValues(b,img)
#   dataType(img) = "INT2S"
  projection(img) = set_projection(file)
  img = as(img, "SpatialGridDataFrame")         #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  outfile = sub("archv", "dos_sr", file)
  writeGDAL(img, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
}