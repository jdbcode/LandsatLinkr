#' Calibrate MSS imagery to TM and make cloud-free composites 
#'
#' Calibrate MSS imagery to TM and make cloud-free composites   
#' @param msswrs1dir character. mss wrs1 directory path
#' @param msswrs2dir character. mss wrs2 directory path
#' @param tmwrs2dir character. tm wrs2 directory path
#' @param index character. spectral index to make composites for. options: "tca", "tcb", "tcg", "tcw"
#' @param outdir character. path to output directory
#' @param runname character. unique name for the composite set
#' @param useareafile character. path to usearea file
#' @param doyears ??? what years to composite
#' @param order character. how to order the images options "sensor_and_doy", "doy", and "none"
#' @param overlap character. how to deal with overlapping images. options: "mean"
#' @param cores numeric. Number of cores to process with options: 1 or 2
#' @param process numeric. integer or vector specifying which processes to run: 1=msscal, 2=mixel
#' @import foreach
#' @import doParallel
#' @export


calibrate_and_composite = function(msswrs1dir,msswrs2dir,tmwrs2dir,oliwrs2dir,index,outdir,runname,useareafile,doyears="all",order="none",overlap="mean", cores=2, process, overwrite=F ,startday, endday, yearadj){
  
  #msscal
  if(1 %in% process ==T){
    #resample MSS
    print("Resampling MSS reflectance and cloudmask images")
    msswrs1srfiles = list.files(msswrs1dir, "dos_sr.tif", recursive=T, full.names=T)
    msswrs1cloudfiles = list.files(msswrs1dir, "cloudmask.tif", recursive=T, full.names=T)
    msswrs2srfiles = list.files(msswrs2dir, "dos_sr.tif", recursive=T, full.names=T)
    msswrs2cloudfiles = list.files(msswrs2dir, "cloudmask.tif", recursive=T, full.names=T)
    files = c(msswrs1srfiles,msswrs1cloudfiles,msswrs2srfiles,msswrs2cloudfiles)
    #cores=2
    if(cores == 2){
      print("...in parallel")
      cl = makeCluster(cores)
      registerDoParallel(cl)
      o = foreach(i=1:length(files), .combine="c",.packages="LandsatLinkr") %dopar% mss_resample(files[i], overwrite=F) #hardwired to not overwrite
      stopCluster(cl)
    } else {for(i in 1:length(files)){o = mss_resample(files[i], overwrite=F)}} #hardwired to not overwrite
    
    
    print("Running msscal")
    t=proc.time()
    msscal(msswrs1dir, msswrs2dir, tmwrs2dir, cores=cores)
    print(proc.time()-t)
  }
  
  #olical
  if(2 %in% process ==T){
    print("Running olical")
    t=proc.time()
    olical(oliwrs2dir, tmwrs2dir, cores=cores, overwrite=overwrite)
    print(proc.time()-t)
  }

  #mixel
  if(3 %in% process ==T){
    print("Running mixel")
    t=proc.time()
    if(index == "all"){
      index = c("tca", "tcb", "tcg", "tcw")
      outdir = c(file.path(outdir,"tca"),file.path(outdir,"tcb"),file.path(outdir,"tcg"),file.path(outdir,"tcw"))
      for(i in 1:length(index)){mixel(msswrs1dir,msswrs2dir,tmwrs2dir,oliwrs2dir,index[i],outdir[i],runname,useareafile,doyears="all",order="none",overlap=overlap, startday=startday, endday=endday, yearadj=yearadj)} #overlap="mean"
    } else {
      outdir = file.path(outdir,index)
      mixel(msswrs1dir,msswrs2dir,tmwrs2dir,oliwrs2dir,index,outdir,runname,useareafile,doyears="all",order="none",overlap=overlap, startday=startday, endday=endday, yearadj=yearadj) #overlap="mean"
    }
    print(proc.time()-t)
  }  
}