#' Find intersecting images
#'
#' From a list of images, find those that interest a given image
#' @param files character. Array of image filenames 
#' @param refimg character. filename of selected image that defines spatial extent of intersection search
#' @import raster
#' @export


find_intersecting_images = function(files, refimg){
  
  #get the extents of the files
  info = matrix(ncol = 4, nrow=length(files))
  print("Getting image extents")
  for(i in 1:length(files)){ 
    print(i)
    img = raster(files[i])
    ext = extent(img)
    info[i,1] = ext@xmin
    info[i,2] = ext@xmax
    info[i,3] = ext@ymin
    info[i,4] = ext@ymax
  }
  
  text = extent(raster(refimg))  
  these = which(info[,3] < text@ymax & info[,4] > text@ymin & info[,2] > text@xmin & info[,1] < text@xmax) 
  goods = files[these]
  return(goods)
}