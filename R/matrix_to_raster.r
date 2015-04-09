#' Converts a matrix to a raster file
#'
#' Converts a matrix to a raster file
#' @param rfile character. full path name of a reference raster file
#' @param rmatrix matrix. a numeric 2-d matrix to be converted to a raster
#' @import raster
#' @export

matrix_to_raster = function(rfile, rmatrix){
  r = raster(rfile)
  cres = 0.5*res(r)[1]
  xmin = xFromCol(r, col=1)-cres
  xmax = xFromCol(r, col=ncol(r))+cres
  ymin = yFromRow(r, row=nrow(r))-cres
  ymax = yFromRow(r, row=1)+cres
  r = raster(rmatrix,xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax,crs=projection(r))
  return(r)
}