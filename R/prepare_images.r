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


prepare_images = function(scenedir, demfile=NULL, proj="default", process=seq(1:5), cores=2, overwrite=F){
  
  targzdir = file.path(scenedir,"targz")
  imgdir = file.path(scenedir,"images")
  
  cfun = function(a, b) NULL
  
  #mssunpackr
  if(all(is.na(match(process,1))) == F){
    print("Running mssunpackr")
    files = list.files(targzdir,"tar.gz",full.names=T)
    t=proc.time()
    if(cores == 2){
      print("...in parallel")
      cl = makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% mssunpackr(files[i], proj=proj, overwrite=overwrite)
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssunpackr(files[i], proj=proj, overwrite=overwrite)}}
    print(proc.time()-t)
  }
  
  #geowarp
  if(all(is.na(match(process,2))) == F){
    print("Running msswarp")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    diagfiles = list.files(imgdir, pattern="cloud_rmse.csv", full.names=T, recursive=T)
    tbl = do.call(rbind, lapply(diagfiles, read.table, header = F,sep = ','))
    reffile = as.character(tbl[order(round(tbl[,3],digits=1), tbl[,2]),][1,1])
    t = proc.time()
    if(cores == 2){
      print("...in parallel")
      cl=makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% msswarp(reffile=reffile, fixfile=files[i], sample=1000)
      stopCluster(cl)
    } else {for(i in 1:length(files)){msswarp(reffile=reffile, fixfile=files[i], sample=1000)}}
    print(proc.time()-t)
  }
  
  #convert to toa reflectance
  if(all(is.na(match(process,3))) == F){
    print("Running mssdn2refl")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    t = proc.time()
    if(cores == 2){
      print("...in parallel")
      cl=makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% mssdn2refl(files[i], overwrite)
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssdn2refl(files[i], overwrite)}}
    print(proc.time()-t)
  }
  
  #convert to surface reflectance
  if(all(is.na(match(process,4))) == F){
    print("Running msscost")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    t = proc.time()
    if(cores == 2){
      print("...in parallel")
      cl=makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% msscost(files[i], overwrite)
      stopCluster(cl)
    } else {for(i in 1:length(files)){msscost(files[i], overwrite)}}
    print(proc.time()-t)
  }
  
  #cloudmask
  if(all(is.na(match(process,5))) == F){
    print("Running msscvm")
    #prepare the topo layers required by msscvm
    newdem = prepare_topo(imgdir, demfile)
    print("...Making masks")
    files = list.files(imgdir, pattern="reflectance", full.names=T, recursive=T)
    t = proc.time()
    for(i in 1:length(files)){msscvm(files[i], newdem, topoprep=T, test=F, overwrite=overwrite)}
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
      cl=makeCluster(cores)
      registerDoParallel(cl)
      t = proc.time()
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% tmunpackr(files[i], proj=proj, overwrite=overwrite)
      stopCluster(cl)
    } else {for(i in 1:length(files)){tmunpackr(files[i], proj=proj, overwrite=overwrite)}}
    print(proc.time()-t)
  }
  
  #unpack oli
  if(all(is.na(match(process,7))) == F){
    print("Running oliunpackr")
    files = list.files(targzdir, pattern="tar.gz", full.names=T, recursive=T)
    reso = 30
    if(reso == 30){cores = 1}
    t = proc.time()
    if(cores == 2){
      cl=makeCluster(cores)
      registerDoParallel(cl)
      t = proc.time()
      o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% oliunpackr(files[i], proj=proj, overwrite=overwrite)
      stopCluster(cl)
    } else {for(i in 1:length(files)){oliunpackr(files[i], proj=proj, overwrite=overwrite)}}
    print(proc.time()-t)
  }
}





