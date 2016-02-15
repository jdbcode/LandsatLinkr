#' Make a usearea file in ENVI (bsq) format
#'
#' Make a usearea file in ENVI (bsq) format
#' @param infile character. full path name to a usearea mask file
#' @param projref charcter. full path to an image file from the scene relevent to the usearea file - can be any .tif file in the "images" directory
#' @import raster
#' @import rgdal
#' @import gdalUtils
#' @export

make_usearea_file_bsq = function(infile, projref){
  
  #infile = "K:/test/composite/useareafile.tif"
  #projref = "K:/test/mss/wrs1/041029/images/1976/LM10410291976181_archv.tif"
  print("Making a copy of use-area file as .bsq for optional use in LandTrendr")
  
  tempout = paste(infile,"_temp.tif",sep="")
  template = r = raster(infile)
  r = as.matrix(r)
  nas = which(is.na(r) == T)
  ones = which(r != 0)
  zero = which(r == 0)
  r[nas,zero] = 0
  r[ones] = 1
  r = setValues(template,r)
  projection(r) = set_projection(projref)
  r = as(r, "SpatialGridDataFrame")
  writeGDAL(r, tempout, drivername = "GTiff", type = "Byte", options="INTERLEAVE=BAND")
  
  outfile = paste(substr(infile,1,(nchar(infile)-3)),"bsq",sep="")
  
  s_srs = projection(raster(tempout))
  t_srs = set_projection(projref)
  gdalwarp(srcfile=tempout, dstfile=outfile, s_srs=s_srs, t_srs=t_srs, order=1, tr=c(30,30), r="near", of="ENVI", dstnodata="none")
  unlink(tempout)
}