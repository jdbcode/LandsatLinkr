#' Prepare topographic layers: elevation, slope, aspect for use in MSS cloud masking 
#'
#' Prepare topographic layers: elevation, slope, aspect for use in MSS cloud masking 
#' @param imgdir character. full path to "images" directory for scene
#' @param demfile character. full path to scene-corresponding DEM file
#' @export


prepare_topo = function(imgdir, demfile){
  print("...Preparing DEM")
  files = list.files(imgdir, pattern="reflectance", full.names=T, recursive=T) #"radiance.tif"
  if(length(files) == 0){
    print(paste("Error - could not find any reflectance (*reflectance.tif) files in this directory:",imgdir))
    stop
  }
  examplefile = files[1]
  dname = dirname(examplefile)
  scenedir = substr(dname,1,nchar(dname)-12)
  topodir = file.path(scenedir,"topo")
  dir.create(topodir, showWarnings=F)
  info = get_metadata(examplefile)
  template = raster(examplefile)
  demname = paste(info$wrstype,"_",info$ppprrr,"_60m","_dem.tif",sep="")
  newdem = file.path(topodir,demname)
  tempdem = sub("dem.tif","temp_dem.tif",newdem)
  newslope = file.path(topodir,sub("dem","slope",demname))
  newasp = file.path(topodir,sub("dem","aspect",demname))
  newill = file.path(topodir,sub("dem","illumination",demname))
  s_srs = projection(raster(demfile)) #template
  t_srs = set_projection(examplefile)
  
  demfiles = list.files(topodir,"dem",full.names=T)
  unlink(demfiles)
  gdalwarp(srcfile=demfile,dstfile=tempdem,
           s_srs=s_srs,t_srs=t_srs, tr=c(60,60), dstnodata=-32768, ot="Int16") #should nodata be set here???
  
  extholder = matrix(ncol = 4, nrow=length(files))
  print("...Making sure DEM is large enough")
  print("......Getting MSS image extents")
  for(i in 1:length(files)){ 
    img = raster(files[i])
    ext = extent(img)
    extholder[i,1] = ext@xmin
    extholder[i,2] = ext@xmax
    extholder[i,3] = ext@ymin
    extholder[i,4] = ext@ymax
  }
  adj=1500
  xmin = min(extholder[,1]) - adj
  xmax = max(extholder[,2]) + adj
  ymin = min(extholder[,3]) - adj
  ymax = max(extholder[,4]) + adj
  
  dem = raster(tempdem)
  demext = extent(dem)
  
  xminokay = demext@xmin <= xmin
  xmaxokay = demext@xmax >= xmax
  yminokay = demext@ymin <= ymin
  ymaxokay = demext@ymax >= ymax
  
  print(paste(".........DEM x minimum is okay:",xminokay))
  print(paste(".........DEM x maximum is okay:",xmaxokay))
  print(paste(".........DEM y minimum is okay:",yminokay))
  print(paste(".........DEM y maximum is okay:",ymaxokay))
  
  if(sum(c(xminokay,xmaxokay,yminokay,ymaxokay)) != 4){
    print("Error - Please make sure DEM has minimum dimensions:")
    print(paste("x minimum:", xmin))
    print(paste("x maximum:", xmax))
    print(paste("y minimum:", ymin))
    print(paste("y maximum:", ymax))
    print(paste("For projection:",t_srs))
    return("Stopping LLR")
  }
  
  #crop the dem
  print("...Croppping the DEM to MSS image set union plus 25 pixel buffer")
  gdal_translate(src_dataset=tempdem, dst_dataset=newdem, projwin=c(xmin,ymax,xmax,ymin))
  tempfiles = list.files(topodir, "temp", full.names=T)
  unlink(tempfiles)
  dem = raster(newdem)
  
  #making slope
  slopefiles = list.files(topodir,"slope",full.names=T)
  unlink(slopefiles)
  print("...Preparing Slope")
  img = terrain(dem, opt="slope")
  projection(img) = set_projection(examplefile)
  img = as(img, "SpatialGridDataFrame")
  writeGDAL(img, newslope, drivername = "GTiff", type = "Float32", options="INTERLEAVE=BAND") #, mvFlag = -32768
  
  #making slope aspect
  aspfiles = list.files(topodir,"aspect",full.names=T)
  unlink(aspfiles)
  print("...Preparing Aspect")
  img = terrain(dem, opt="aspect")
  projection(img) = set_projection(examplefile)
  img = as(img, "SpatialGridDataFrame")
  writeGDAL(img, newasp, drivername = "GTiff", type = "Float32", options="INTERLEAVE=BAND") #, mvFlag = -32768  
  
  return(newdem)
}
