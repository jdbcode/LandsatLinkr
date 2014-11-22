#' write out calibrated mss tc files
#'
#' write out calibrated mss tc files
#' @param mss_file filename mss radiance file 
#' @param index character what tc transformation to be written out
#' @param coef array calibrartion coefficents for transformation to tc index
#' @param outfile filename full path of the outfile 


write_tc = function(mss_file, index, coef, outfile){
  if(index == "tcw"){#create wetness
    b1 = as.matrix(raster(mss_file, 1))
    b2 = as.matrix(raster(mss_file, 2))
    b3 = as.matrix(raster(mss_file, 3))
    b4 = as.matrix(raster(mss_file, 4))
    tc = coef$yint + (b1*coef$b1c + b2*coef$b2c + b3*coef$b3c + b4*coef$b4c)
  } else {
    if(index == "tca" | index == "tcb"){band=1}
    if(index == "tcg"){band=2}
    if(index == "tcw"){band=3}
    tc = as.matrixraster(mss_file, band)
    tc = coef$yint + (tc*coef$coef)
  }
  #tc[tc == coef$yint] = -32768
  tc = setValue(raster(mss_file), tc)
  projection(tc) = set_projection(mss_file)
  tc = as(tc, "SpatialGridDataFrame")
  writeGDAL(tc, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768)
}

