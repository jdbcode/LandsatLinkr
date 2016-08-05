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
  
  if(mode == "calibrate"){
    dir = substr(dirname(mss_file),1,nchar(dirname(mss_file))-12)
    tcfiledir = file.path(dir,"calibration","aggregate_model_tc_imgs")
    dir.create(tcfiledir, recursive=T, showWarnings=F)
    tcfile = file.path(tcfiledir,sub("dos_sr_30m.tif","tc_30m.tif", basename(mss_file)))
    tcafile = file.path(tcfiledir,sub("dos_sr_30m.tif","tca_30m.tif", basename(mss_file)))
  }
  if(mode == "apply"){
    tcfile = sub("dos_sr_30m.tif","tc_30m.tif", mss_file)
    tcafile = sub("dos_sr_30m.tif","tca_30m.tif", mss_file)
  }
  
  #read in the dos sr image
  ref = raster(mss_file)
  b1=as.matrix(raster(mss_file,1))
  b2=as.matrix(raster(mss_file,2))
  b3=as.matrix(raster(mss_file,3))
  b4=as.matrix(raster(mss_file,4))
  
  #transform to tc
  tcb = round((((b1*bcoef[2])+(b2*bcoef[3])+(b3*bcoef[4])+(b4*bcoef[5])) + bcoef[1]))
  tcg  = round((((b1*gcoef[2])+(b2*gcoef[3])+(b3*gcoef[4])+(b4*gcoef[5])) + gcoef[1]))
  tcw = round((((b1*wcoef[2])+(b2*wcoef[3])+(b3*wcoef[4])+(b4*wcoef[5])) + wcoef[1]))
  
  b1=b2=b3=b4=0 #save memory  

  tca = atan(tcg/tcb) * (180/pi) * 100 #need to multiply by (180/pi) to get degrees because atan() returns radians, 100 is a scalar to preserve two decimal places
  
  tca = matrix_to_raster(mss_file, tca)
  tcb = matrix_to_raster(mss_file, tcb)
  tcg = matrix_to_raster(mss_file, tcg)
  tcw = matrix_to_raster(mss_file, tcw)
  
  projection(tca) = set_projection(mss_file)
  tca = as(tca, "SpatialGridDataFrame")
  writeGDAL(tca, tcafile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  tca=0 #recover some memory

  #stack tc and write out
  tc = stack(tcb,tcg,tcw)
  projection(tc) = set_projection(mss_file)
  tc = as(tc, "SpatialGridDataFrame")
  writeGDAL(tc, tcfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  return(1)
}
