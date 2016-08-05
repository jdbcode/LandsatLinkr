#' Trim NA rows and columns from an image
#'
#' Trim NA rows and columns from an image
#' @param imgfile character. input image filename
#' @param outimg character. output image filename
#' @param maskfile character. input mask filename (optional) 
#' @param outmask character. input mask filename (optional)
#' @import raster
#' @import gdalUtils
#' @export

trim_na_rowcol = function(imgfile, outimg, maskfile, outmask){
  #trim the stack
  x = raster(imgfile, band=1)
  cres = 0.5*res(x)
  y = x
  x = matrix(as.array(x),nrow=nrow(x),ncol=ncol(x))
  r.na = c.na <- c()
  for(i in 1:nrow(x)) r.na <- c(r.na, all(is.na(x[i,])))
  for(i in 1:ncol(x)) c.na <- c(c.na, all(is.na(x[,i])))
  r1 = 1 + which(diff(which(r.na))>1)[1] 
  r2 = nrow(x) -  which(diff(which(rev(r.na)))>1)[1]
  c1 = 1 + which(diff(which(c.na))>1)[1] 
  c2 = ncol(x) - which(diff(which(rev(c.na)))>1)[1]
  
  #if there are no NA rows and cols, then set to default row and col start and end
  if(is.na(r1)){
    if(r.na[1] == T){r1 = 1+length(which(r.na))} else {r1 = 1}
  }
  if(is.na(r2)){
    if(rev(r.na)[1] == T){r2 = nrow(x) - length(which(r.na))} else {r2 = nrow(x)}
  }
  if(is.na(c1)){
    if(c.na[1] == T){c1 = 1+length(which(c.na))} else {c1 = 1}
  }
  if(is.na(c2)){
    if(rev(c.na)[1] == T){c2 = ncol(x) - length(which(c.na))} else {c2 = ncol(x)}
  }
  
  xs = xFromCol(y,col=c(c1,c2)) + c(-1,1)*cres[1]
  ys = yFromRow(y,row=c(r2,r1)) + c(-1,1)*cres[2]
    
  #write out the trimmed file
  gdal_translate(src_dataset=imgfile, dst_dataset=outimg, of="GTiff", co="INTERLEAVE=BAND", projwin=c(xs[1],ys[2],xs[2],ys[1]))
  if(file.exists(maskfile) == T){gdal_translate(src_dataset=maskfile, dst_dataset=outmask, of="GTiff", co="INTERLEAVE=BAND", projwin=c(xs[1],ys[2],xs[2],ys[1]))}
}