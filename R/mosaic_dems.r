#' Create a DEM mosaic from a direcory of DEM's
#'
#' Create a DEM mosaic from a direcory of DEM's
#' @param dir character. The path to a directory containing DEM files to be mosaicked
#' @import raster
#' @export


mosaic_dems = function(dir, proj){
  
  align = function(img, refimg){
    img = raster(img)
    imgex = alignExtent(img, refimg, snap="near")
    extent(img) = imgex
    return(img)
  }
  
  reso = 60
  demfiles = list.files(dir, full.names=T)
  
  for(i in 1:length(demfiles)){
    demfile = demfiles[i]
    print(paste("reprojecting file:",demfile))
    bname = basename(demfile)
    extension = substr(bname,(nchar(bname)-3),nchar(bname))
    dstfile = sub(extension, "_reprojected.tif", demfile) 
  
    gdalwarp(srcfile=demfile, dstfile=dstfile, 
      t_srs=proj, of="GTiff",
      r="near", dstnodata=-32768, multi=T,
      tr=c(reso,reso), co="INTERLEAVE=BAND")
  }
  
  demfiles = list.files(dir, pattern="reprojected.tif$", full.names=T)
  len = length(demfiles)
  refimg = raster(demfiles[1])
  
  for(i in 1:len){
    if(i == 1){mergeit = "r1"} else {mergeit = paste(mergeit,",r",i, sep="")}
    #open the image as raster and aligns it
    dothis = paste("r",i,"=align(demfiles[",i,"], refimg)", sep="")
    eval(parse(text=dothis))
    if(i == len){mergeit = paste("big = mosaic(",mergeit,",fun=mean,na.rm=T,tolerance=0.5)", sep="")}
  }
  
  #run the merge function
  print("calculating DEM mosiac")
  if(len == 1){big = r1} else {eval(parse(text=mergeit))} #only run merge it if there are multiple files to merge
  
  #write out the DEM mosiac
  print("writing DEM mosiac")
  outfile = file.path(dir,"dem_mosaic.tif")
  writeRaster(big, outfile, format="GTiff", datatype = "INT2S",overwrite=T, bandorder = "BSQ",options=c("COMPRESS=NONE"))

  
}

