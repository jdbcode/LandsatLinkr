#' Find the spatial intersection of a collection of images 
#'
#' Find the spatial intersection of a collection of images
#' @param files array of image files names 
#' @import raster
#' @export


get_intersection = function(files){
  int = intersect(extent(raster(files[1])),extent(raster(files[2])))
  if(length(files) >= 3){for(i in 3:length(files))int = intersect(extent(raster(files[i])), int)}
  return(int)
}