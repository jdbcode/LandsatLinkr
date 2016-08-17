#' Resamples MSS dos_sr and cloudmask images to 30m
#'
#' Resamples MSS dos_sr and cloudmask images to 30m
#' @param file character. full path to either an MSS *reflectance.tif or *cloudmask.tif file
#' @import gdalUtils
#' @export

mss_resample = function(file, overwrite=F){
  type = c(length(grep("dos_sr", file)), length(grep("cloudmask", file)))
  
  if(type[1] == 1){
    check = file_check(file,"dos_sr_30m.tif",overwrite)
    if(check == 0){return(0)}
    
    newfile = sub("dos_sr", "dos_sr_30m", file)
    gdalwarp(srcfile=file, dstfile=newfile,tr=c(30,30),
             srcnodata=-32768, dstnodata=-32768, multi=T, r="cubic")
    
    return(1)
  }
  
  
  if(type[2] == 1){
    check = file_check(file,"cloudmask_30m.tif",overwrite)
    if(check == 0){return(0)}
    
    newfile = sub("cloudmask", "cloudmask_30m", file)
    gdalwarp(srcfile=file, dstfile=newfile,tr=c(30,30),
             srcnodata=255, dstnodata=255, multi=T)
    
    return(1)
  }
}

