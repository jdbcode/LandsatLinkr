#' Run prepareMSS in parallel 
#'
#' Run prepareMSS in parallel using 
#' @param scenedir character. scene file path
#' @param demfile character. full path to scene-corresponding DEM file
#' @param process numeric. integer or vector specifying which processes to run 1=mssunpackr, 2=msswarp, 3=mssdn2rad, 4=mssatcor, 5=msscvm, 6=tmunpackr
#' @param cores numeric. number of cores to use for parallel processing
#' @import foreach
#' @import doParallel
#' @export


doinparallel = function(scenedir, demfile, process=seq(1:5),cores=2){
  
  #library(doParallel)
  
  targzdir = file.path(scenedir,"targz")
  imgdir = file.path(scenedir,"images")
  
  #mssunpackr
  if(all(is.na(match(process,1))) == F){
    print("Running mssunpackr")
    files = list.files(targzdir,"tar.gz",full.names=T)
    cl = makeCluster(cores)
    registerDoParallel(cl)
    t=proc.time()
    o = foreach(i=1:length(files), .combine="c",.packages="MSScvm") %dopar% mssunpackr(files[i], proj="albers") #
    print(proc.time()-t)
    stopCluster(cl)
  }
  
  #geowarp
  if(all(is.na(match(process,2))) == F){
    #dir = "E:/mss/wrs1/036032/images"
    print("Running msswarp")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    diagfiles = list.files(imgdir, pattern="cloud_rmse.csv", full.names=T, recursive=T)
    tbl = do.call(rbind, lapply(diagfiles, read.table, header = F,sep = ','))
    reffile = as.character(tbl[order(round(tbl[,3],digits=1), tbl[,2]),][1,1])
    cl=makeCluster(cores)
    registerDoParallel(cl)
    t = proc.time()
    o = foreach(i=1:length(files), .combine="c",.packages="MSScvm") %dopar% msswarp(reffile=reffile, fixfile=files[i], sample=1000)
    print(proc.time()-t)
    stopCluster(cl)
  }
  
  #convert to radiance
  if(all(is.na(match(process,3))) == F){
    print("Running mssdn2rad")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    cl=makeCluster(cores)
    registerDoParallel(cl)
    t = proc.time()
    o = foreach(i=1:length(files), .combine="c",.packages="MSScvm") %dopar% mssdn2rad(files[i])
    print(proc.time()-t)
    stopCluster(cl)
  }
  
  #surface reflectance
  if(all(is.na(match(process,4))) == F){
    print("Running mssatcor")
    files = list.files(imgdir, pattern="archv.tif", full.names=T, recursive=T)
    cl=makeCluster(cores) #high memory with 2
    registerDoParallel(cl)
    t = proc.time()
    o = foreach(i=1:length(files), .combine="c",.packages="MSScvm") %dopar% mssatcor(files[i], outtype=3)
    print(proc.time()-t)
    stopCluster(cl)
  }
  
  #cloudmask
  if(all(is.na(match(process,5))) == F){
    print("Running msscvm")
    files = list.files(imgdir, pattern="radiance.tif", full.names=T, recursive=T)
    #demfile = "K:/gis_data/dems/wrs_dem/wrs1_036032_60m_dem.tif"
    cl=makeCluster(cores) #high memory with 2
    registerDoParallel(cl)
    t = proc.time()
    o = foreach(i=1:length(files), .combine="c",.packages="MSScvm") %dopar% msscvm(files[i], demfile)
    print(proc.time()-t)
    stopCluster(cl)
  }
  
  #unpack tm
  if(all(is.na(match(process,6))) == F){
    print("Running msscvm")
    files = list.files(targzdir, pattern="tar.gz", full.names=T, recursive=T)
    cl=makeCluster(cores) #high memory with 2
    registerDoParallel(cl)
    cfun <- function(a, b) NULL
    t = proc.time()
    o = foreach(i=1:length(files), .combine="cfun",.packages="MSScvm") %dopar% tmunpackr(files[i], proj="albers", reso=60)
    print(proc.time()-t)
    stopCluster(cl)
  }
}





