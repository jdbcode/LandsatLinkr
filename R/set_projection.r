#' Set the projection of raster
#'
#' Set the projection of a raster from the corresponding *proj.txt file. This solves a problem with incorrect projection parameters for albers when read by 'raster'
#' @param file Filename of MSS image with DN values
#' @export

set_projection = function(file){
  projfile = file.path(dirname(file),paste(substr(basename(file),1,16),"_proj.txt",sep=""))
  return(readLines(projfile))
}
