#' Prepare image MSS and TM images for calibration/compositing  
#'
#' Prepare image MSS and TM images for calibration/compositing  
#' @param scenedir character. scene file path
#' @param demfile character. full path to scene-corresponding DEM file
#' @param process numeric. integer or vector specifying which processes to run 1=mssunpackr, 2=msswarp, 3=mssdn2rad, 4=mssatcor, 5=msscvm, 6=tmunpackr
#' @param cores numeric. number of cores to use for parallel processing
#' @import foreach
#' @import doParallel
#' @export


prepare_images = function(scenedir, demfile=NULL, proj="default", reso=60, process=seq(1:5), cores=2){
  
  targzdir = file.path(scenedir,"targz")
  imgdir = file.path(scenedir,"images")
  
  #mssunpackr
  if(all(is.na(match(process,1))) == F){
    print("Running mssunpackr")
    files = list.files(targzdir,"tar.gz",full.names=T)
    t=proc.time()
    if(cores == 2){
      cl = makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% mssunpackr(files[i], proj=proj, reso=reso)
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssunpackr(files[i], proj=proj, reso=reso)}}
    print(proc.time()-t)
  }
  
  #geowarp
  if(all(is.na(match(process,2))) == F){
    #dir = "E:/mss/wrs1/036032/images"
    print("Running msswarp")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    diagfiles = list.files(imgdir, pattern="cloud_rmse.csv", full.names=T, recursive=T)
    tbl = do.call(rbind, lapply(diagfiles, read.table, header = F,sep = ','))
    reffile = as.character(tbl[order(round(tbl[,3],digits=1), tbl[,2]),][1,1])
    t = proc.time()
    if(cores == 2){
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
    t = proc.time()
    if(cores == 2){
      cl=makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% mssdn2refl(files[i]) #mssdn2rad
      stopCluster(cl)
    } else {for(i in 1:length(files)){mssdn2refl(files[i])}} #mssdn2rad
    print(proc.time()-t)
  }
  
  #surface reflectance
  if(all(is.na(match(process,4))) == F){
    print("Running mssatcor")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    t = proc.time()
    if(cores == 2){
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
    #demfile = "K:/gis_data/dems/wrs_dem/wrs1_036032_60m_dem.tif"
    t = proc.time()
    if(cores == 2){
      cl=makeCluster(cores) #high memory with 2
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% msscvm(files[i], demfile)
      stopCluster(cl)
    } else {for(i in 1:length(files)){msscvm(files[i], demfile)}}
    print(proc.time()-t)
  }
  
  #unpack tm
  if(all(is.na(match(process,6))) == F){
    print("Running tmunpackr")
    files = list.files(targzdir, pattern="tar.gz", full.names=T, recursive=T)
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





