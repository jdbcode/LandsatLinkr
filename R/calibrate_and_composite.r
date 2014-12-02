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
#' @param order character. how to order the images options "sensor_and_doy" and "doy"
#' @param overlap character. how to deal with overlapping images. options: "mean"
#' @param cores numeric. Number of cores to process with options: 1 or 2
#' @param process numeric. integer or vector specifying which processes to run: 1=msscal, 2=mixel
#' @import foreach
#' @import doParallel
#' @export


calibrate_and_composite = function(msswrs1dir,msswrs2dir,tmwrs2dir,index,outdir,runname,useareafile,doyears="all",order="sensor_and_doy",overlap="mean", cores=2, process=c(1,2)){
  
  #msscal
  if(all(is.na(match(process,1))) == F){
    print("Running msscal")
    t=proc.time()
    msscal(msswrs1dir, msswrs2dir, tmwrs2dir, cores=cores)
    print(proc.time()-t)
  }
  
  #mssunpackr
  if(all(is.na(match(process,2))) == F){
    print("Running mixel")
    t=proc.time()
    mixel(msswrs1dir,msswrs2dir,tmwrs2dir,index,outdir,runname,useareafile,doyears="all",order="sensor_and_doy",overlap="mean")
    print(proc.time()-t)
  }
}