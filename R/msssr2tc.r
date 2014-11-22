#' Create MSS TC from surface reflectance and modeled coefficients
#'
#' Create MSS TC from surface reflectance and modeled coefficients
#' @param mss_file character. full path name to MSS surface reflectance file
#' @param bcoef numeric. numeric array containing the tcb coefficients
#' @param gcoef numeric. numeric array containing the tcg coefficients
#' @param wcoef numeric. numeric array containing the tcw coefficients
#' @param mode. character. how to deal with the outputs options: "calibrate" or "apply"
#' @import raster
#' @export


msssr2tc = function(mss_file,bcoef,gcoef,wcoef,mode){
  ref = raster(mss_file)
  b = brick(mss_file)
  img = as.array(b)
  
  bright = round((((img[,,1]*bcoef[2])+(img[,,2]*bcoef[3])+(img[,,3]*bcoef[4])+(img[,,4]*bcoef[5])) + bcoef[1]))
  green  = round((((img[,,1]*gcoef[2])+(img[,,2]*gcoef[3])+(img[,,3]*gcoef[4])+(img[,,4]*gcoef[5])) + gcoef[1]))
  wet    = round((((img[,,1]*wcoef[2])+(img[,,2]*wcoef[3])+(img[,,3]*wcoef[4])+(img[,,4]*wcoef[5])) + wcoef[1]))
  
  if(mode == "calibrate"){
    dir = substr(dirname(mss_file),1,nchar(dirname(mss_file))-12)
    tcname = sub("dos_sr.tif","tc.tif", basename(mss_file))
    tcaname = sub("dos_sr.tif","tca.tif", basename(mss_file))
    tcfiledir = file.path(dir,"calibration","composite_model_tc_imgs")
    dir.create(tcfiledir, recursive=T, showWarnings=F)
    tcfile = file.path(tcfiledir,sub("dos_sr.tif","tc.tif", basename(mss_file)))
    tcafile = file.path(tcfiledir,sub("dos_sr.tif","tca.tif", basename(mss_file)))
  }
  if(mode == "apply"){
    tcfile = sub("dos_sr.tif","tc.tif", mss_file)
    tcafile = sub("dos_sr.tif","tca.tif", mss_file)
  }
  #calc tc and convert to a raster
  brt = setValues(ref,bright)
  grn = setValues(ref,green)
  wet = setValues(ref,wet)
  tc = stack(brt,grn,wet)
  projection(tc) = set_projection(mss_file)
  tc = as(tc, "SpatialGridDataFrame")
  writeGDAL(tc, tcfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  tca = atan(green/bright) * (180/pi) * 100 #need to multiply by (180/pi) to get degrees because atan() returns radians, 100 is a scalar to preserve two decimal places
  tca = setValues(ref,tca)
  projection(tca) = set_projection(mss_file)
  tca = as(tca, "SpatialGridDataFrame")
  writeGDAL(tca, tcafile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  return(1)
}
