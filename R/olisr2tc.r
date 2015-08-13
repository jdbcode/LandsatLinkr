#' Create oli TC from surface reflectance and modeled coefficients
#'
#' Create oli TC from surface reflectance and modeled coefficients
#' @param oli_file character. full path name to oli surface reflectance file
#' @param bcoef numeric. numeric array containing the tcb coefficients
#' @param gcoef numeric. numeric array containing the tcg coefficients
#' @param wcoef numeric. numeric array containing the tcw coefficients
#' @param mode. character. how to deal with the outputs options: "calibrate" or "apply"
#' @import raster
#' @export


olisr2tc = function(oli_file,bcoef,gcoef,wcoef,mode,overwrite=F){
  
  if(mode == "calibrate"){
    dir = substr(dirname(oli_file),1,nchar(dirname(oli_file))-12)
    tcfiledir = file.path(dir,"calibration","aggregate_model_tc_imgs")
    dir.create(tcfiledir, recursive=T, showWarnings=F)
    tcfile = file.path(tcfiledir,sub("l8sr.tif","tc.tif", basename(oli_file)))
    tcafile = file.path(tcfiledir,sub("l8sr.tif","tca.tif", basename(oli_file)))
  }
  if(mode == "apply"){
    tcfile = sub("l8sr.tif","tc.tif", oli_file)
    tcafile = sub("l8sr.tif","tca.tif", oli_file)
  }
  
  check = file_check(oli_file,"l8sr_tc.tif",overwrite)
  if(check == 0){return(0)}
  
  #tasseled cap
  b2 = as.matrix(raster(oli_file, 2))
  b3 = as.matrix(raster(oli_file, 3))
  b4 = as.matrix(raster(oli_file, 4))
  b5 = as.matrix(raster(oli_file, 5))
  b6 = as.matrix(raster(oli_file, 6))
  b7 = as.matrix(raster(oli_file, 7))
  
  bright = (b2*bcoef[1])+(b3*bcoef[2])+(b4*bcoef[3])+(b5*bcoef[4])+(b6*bcoef[5])+(b7*bcoef[6])
  green = (b2*gcoef[1])+(b3*gcoef[2])+(b4*gcoef[3])+(b5*gcoef[4])+(b6*gcoef[5])+(b7*gcoef[6])
  wet = (b2*wcoef[1])+(b3*wcoef[2])+(b4*wcoef[3])+(b5*wcoef[4])+(b6*wcoef[5])+(b7*wcoef[6])
  
  b2=b3=b4=b5=b6=b7=0


  tcb = matrix_to_raster(oli_file,bright)
  tcg = matrix_to_raster(oli_file,green)
  tcw = matrix_to_raster(oli_file,wet)
  wet=0  
  tca = atan(green/bright) * (180/pi) * 100
  bright=green=0  
  
  tca = matrix_to_raster(oli_file,tca)
  
  outbase = substr(basename(oli_file),1,16)
  outdir = dirname(oli_file)
  temptcb = sub("l8sr.tif","_temptcb.tif", oli_file)
  temptcg = sub("l8sr.tif","_temptcg.tif", oli_file)
  temptcw = sub("l8sr.tif","_temptcw.tif", oli_file)
  projection(tcb) = set_projection(tcfile)
  projection(tcg) = set_projection(tcfile)
  projection(tcw) = set_projection(tcfile)
  tc = as(tcb, "SpatialGridDataFrame")
  tcb=0
  writeGDAL(tc, temptcb, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  tc = as(tcg, "SpatialGridDataFrame")
  tcg=0
  writeGDAL(tc, temptcg, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  tc = as(tcw, "SpatialGridDataFrame")
  tcw=0
  writeGDAL(tc, temptcw, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  bands = c(temptcb,temptcg,temptcw)
  temptcs = sub("l8sr.tif","_temptcstack.vrt", oli_file) 
  gdalbuildvrt(gdalfile=bands, output.vrt = temptcs, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=temptcs, dst_dataset=tcfile, of = "GTiff", co="INTERLEAVE=BAND")
  
  projection(tca) = set_projection(tcfile)
  tc = as(tca, "SpatialGridDataFrame")
  tca=0
  writeGDAL(tc, tcafile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  #delete temporary files
  unlink(c(temptcb,temptcg,temptcw,temptcs))
  
  return(1)
}
