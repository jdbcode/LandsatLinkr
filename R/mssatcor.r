#' Apply dark object subtraction atmospheric correction
#'
#' Apply dark object subtraction atmospheric correction with option of out putting the image files in units of: DN, radiance, or surface reflectance
#' @param file The full path name of the *archv file
#' @param outtype coded integer designating what units are desired in the output image: 1=DN, 2=radiance, 3=surface reflectance
#' @import raster
#' @export


mssatcor = function(file, outtype){
  
  drkobjv = mssdofindr(file,pngout=T)
    
  if(outtype == 1){
    mssdos(file, drkobjv, "dn")
  }
  if(outtype == 2){
    mssdos(file, drkobjv, "radiance")
  }
  if(outtype == 3){
    msscost(file, drkobjv)
  }  
}