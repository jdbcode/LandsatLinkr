#' Determine MSS band Dark Object Values
#'
#' Determine the Dark Object Values in MSS imagery per band
#' @param file The full path name of the *archv file
#' @param pngout Logical. If TRUE then a 4-panel .png image file will be created that for each MSS band displays a truncated histogram and the estimated Dark Object Value
#' @import raster
#' @export


mssdos = function(file, drkobjv, outtype){
  
  info = get_metadata(file)
  b = brick(file)
  img = as.array(b)
  
  if(outtype == "dn"){
    img[,,1] = img[,,1]-drkobjv[1]
    img[,,2] = img[,,2]-drkobjv[2]
    img[,,3] = img[,,3]-drkobjv[3]
    img[,,4] = img[,,4]-drkobjv[4]
    img[img < 0] = 0  #you can't have negative radiance - set negative values to 0
    
    img = setValues(b, img)
    projection(img) = set_projection(file)
    img = as(img, "SpatialGridDataFrame")
    type = "Byte"
    mvFlag = 0
    outfile = sub("archv.tif", "dn_dos.tif", file)
  }
  
  if(outtype == "radiance"){
    img[,,1] = round(100*(((info$b1gain*img[,,1])+info$b1bias) - abs(((info$b1gain*drkobjv[1])+info$b1bias))))
    img[,,2] = round(100*(((info$b2gain*img[,,2])+info$b2bias) - abs(((info$b1gain*drkobjv[2])+info$b1bias))))
    img[,,3] = round(100*(((info$b3gain*img[,,3])+info$b3bias) - abs(((info$b1gain*drkobjv[3])+info$b1bias))))
    img[,,4] = round(100*(((info$b4gain*img[,,4])+info$b4bias) - abs(((info$b1gain*drkobjv[4])+info$b1bias))))
    img[img<0] = 0 #you can't have negative radiance - set negative values to 0
    
    img = setValues(b, img)
    projection(img) = set_projection(file)
    img = as(img, "SpatialGridDataFrame")
    type = "Int16"
    mvFlag = -32768
    outfile = sub("archv.tif", "radiance_dos.tif", file)
  }
  
  writeGDAL(img, outfile, drivername = "GTiff", type = type, mvFlag = mvFlag)  
}

