#' Prepare image MSS and TM images for calibration/compositing  
#'
#' Prepare image MSS and TM images for calibration/compositing  
#' @param scenedir character. scene file path
#' @param demfile character. full path to scene-corresponding DEM file
#' @param process numeric. integer or vector specifying which processes to run 1=mssunpackr, 2=msswarp, 3=mssdn2rad, 4=mssatcor, 5=msscvm, 6=tmunpackr
#' @param cores numeric. number of cores to use for parallel processing
#' @import foreach
#' @import doParallel
#' @import raster
#' @export


prepare_images = function(scenedir, demfile=NULL, proj="default", reso=60, process=seq(1:5), cores=2){
  
  targzdir = file.path(scenedir,"targz")
  imgdir = file.path(scenedir,"images")
  
  #mssunpackr
  if(all(is.na(match(process,1))) == F){
    print("Running mssunpackr")
    files = list.files(targzdir,"tar.gz",full.names=T)
    reso = 60
    if(reso == 30){cores = 1}
    t=proc.time()
    if(cores == 2){
      print("...in parallel")
      cl = makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% mssunpackr(files[i], proj=proj, reso=reso)
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssunpackr(files[i], proj=proj, reso=reso)}}
    print(proc.time()-t)
  }
  
  #geowarp
  if(all(is.na(match(process,2))) == F){
    print("Running msswarp")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    diagfiles = list.files(imgdir, pattern="cloud_rmse.csv", full.names=T, recursive=T)
    tbl = do.call(rbind, lapply(diagfiles, read.table, header = F,sep = ','))
    reffile = as.character(tbl[order(round(tbl[,3],digits=1), tbl[,2]),][1,1])
    if(xres(raster(files[1])) == 30){cores = 1}
    t = proc.time()
    if(cores == 2){
      print("...in parallel")
      cl=makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% msswarp(reffile=reffile, fixfile=files[i], sample=1000)
      stopCluster(cl)
    } else {for(i in 1:length(files)){msswarp(reffile=reffile, fixfile=files[i], sample=1000)}}
    print(proc.time()-t)
  }
  
  #convert to toa reflectance
  if(all(is.na(match(process,3))) == F){
    print("Running mssdn2refl")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    if(xres(raster(files[1])) == 30){cores = 1}
    t = proc.time()
    if(cores == 2){
      print("...in parallel")
      cl=makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% mssdn2refl(files[i]) #mssdn2rad
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssdn2refl(files[i])}} #mssdn2rad
    print(proc.time()-t)
  }
  
  #convert to surface reflectance
  if(all(is.na(match(process,4))) == F){
    print("Running mssatcor")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    if(xres(raster(files[1])) == 30){cores = 1}
    t = proc.time()
    if(cores == 2){
      print("...in parallel")
      cl=makeCluster(cores) #high memory with 2
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% mssatcor(files[i], outtype=3)
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssatcor(files[i], outtype=3)}}
    print(proc.time()-t)
  }
  
  #cloudmask
  if(all(is.na(match(process,5))) == F){
    print("Running msscvm")
    files = list.files(imgdir, pattern="reflectance", full.names=T, recursive=T) #"radiance.tif"
    
    print("...Preparing DEM")
    #file = "E:/llr_test/30m/mss/wrs1/036032/images/1973/LM10360321973191_reflectance.tif"
    #demfile = "K:/gis_data/dems/wrs_dem/wrs1_036032_60m_dem.tif"
    
    examplefile = files[1]
    dname = dirname(examplefile)
    scenedir = substr(dname,1,nchar(dname)-12)
    topodir = file.path(scenedir,"topo")
    dir.create(topodir)
    info = get_metadata(examplefile)
    template = raster(examplefile)
    reso = xres(template)
    demname = paste(info$wrstype,"_",info$ppprrr,"_",reso,"m","_dem.tif",sep="")
    newdem = file.path(topodir,demname)
    newslope = sub("dem","slope",newdem)
    newasp = sub("dem","aspect",newdem)
    newill = sub("dem","illumination",newdem)
    s_srs = projection(template)
    t_srs = set_projection(examplefile)
    
    gdalwarp(srcfile=demfile,dstfile=newdem,
             s_srs=s_srs,t_srs=t_srs, tr=c(reso,reso), dstnodata=-32768, ot="Int16")
    
    print("...Preparing Slope")
    dem = raster(newdem)
    img = terrain(dem, opt="slope")
    projection(img) = set_projection(examplefile)
    img = as(img, "SpatialGridDataFrame")         #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
    writeGDAL(img, newslope, drivername = "GTiff", type = "Float32", options="INTERLEAVE=BAND") #, mvFlag = -32768
    
    print("...Preparing Aspect")
    img = terrain(dem, opt="aspect")
    projection(img) = set_projection(examplefile)
    img = as(img, "SpatialGridDataFrame")         #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
    writeGDAL(img, newasp, drivername = "GTiff", type = "Float32", options="INTERLEAVE=BAND") #, mvFlag = -32768  
    
    img=0
    
    print("...Making masks")
    t = proc.time()
#     if(cores == 2){
#       cl=makeCluster(cores) #high memory with 2
#       registerDoParallel(cl)
#       o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% msscvm2(files[i], demfile)
#       stopCluster(cl)
#     } else {for(i in 1:length(files)){msscvm2(files[i], demfile)}}
    for(i in 1:length(files)){msscvm2(files[i], newdem, topoprep=T, test=F)} #demfile
    print(proc.time()-t)
  }
  
  #unpack tm
  if(all(is.na(match(process,6))) == F){
    print("Running tmunpackr")
    files = list.files(targzdir, pattern="tar.gz", full.names=T, recursive=T)
    reso = 30
    if(reso == 30){cores = 1}
    t = proc.time()
    if(cores == 2){
      cl=makeCluster(cores) #high memory with 2
      registerDoParallel(cl)
      cfun <- function(a, b) NULL
      t = proc.time()
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% tmunpackr(files[i], proj=proj, reso=reso)
      stopCluster(cl)
    } else {for(i in 1:length(files)){tmunpackr(files[i], proj=proj, reso=reso)}}
    print(proc.time()-t)
  }
}





