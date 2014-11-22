#' Convert MSS DN values to radiance
#'
#' Convert MSS DN values to radiance
#' @param file character. full file path of *archv.bsq file
#' @import raster
#' @export


mssdn2rad = function(file){
  
  info = get_metadata(file)
  b = brick(file)
  img = as.array(b)  
  
  img[,,1] = round(100*((info$b1gain*img[,,1])+info$b1bias))
  img[,,2] = round(100*((info$b2gain*img[,,2])+info$b2bias))
  img[,,3] = round(100*((info$b3gain*img[,,3])+info$b3bias))
  img[,,4] = round(100*((info$b4gain*img[,,4])+info$b4bias))
  
  img[img<0] = 0 #you can't have negative radiance - set negative values to 0
  
  img = setValues(b,img)
  dataType(img) = "INT2S"
  projection(img) = set_projection(file)
  img = as(img, "SpatialGridDataFrame")         #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  outfile = sub("archv", "radiance", file)
  writeGDAL(img, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  return(1)
}
  
