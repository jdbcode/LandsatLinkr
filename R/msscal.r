#' Calibrate MSS images to TM images
#'
#' Calibrate MSS images to TM images using linear regression
#' @param msswrs1dir character. MSS WRS-1 scene directory path
#' @param msswrs2dir character. MSS WRS-2 scene directory path
#' @param tmwrs2dir character. TM WRS-2 scene directory path
#' @param cores numeric. Number of cores to process with options: 1 or 2
#' @export


msscal = function(msswrs1dir, msswrs2dir, tmwrs2dir, cores=2){
  
  t=proc.time()
  
  mssfiles = list.files(msswrs2dir, "dos_sr.tif", recursive=T, full.names=T)
  tmfiles = list.files(tmwrs2dir, "tc.tif", recursive=T, full.names=T)
  
  mssimgid = substr(basename(mssfiles),3,16)
  tmimgid = substr(basename(tmfiles),3,16)
  
  thesemss = which(mssimgid %in% tmimgid)
  mssfiles = mssfiles[thesemss]
  
  thesetm = which(tmimgid %in% mssimgid)
  tmfiles = tmfiles[thesetm]
  
  if(length(mssfiles) == length(tmfiles)){
    mssimgid = substr(basename(mssfiles),3,16)
    tmimgid = substr(basename(tmfiles),3,16)
    mssf = mssfiles[order(mssimgid)]
    tmf = tmfiles[order(tmimgid)]
  }
  
  print("single image pair modeling")
  if(cores==2){
    cl = makeCluster(cores)
    registerDoParallel(cl)
    cfun <- function(a, b) NULL
    o = foreach(i=1:length(mssf), .combine="cfun",.packages="LandsatLinkr") %dopar% msscal_single(mssf[i], tmf[i]) #
    stopCluster(cl)
  } else {for(i in 1:length(mssf)){msscal_single(mssf[i], tmf[i])}}
  
  dir = file.path(msswrs2dir,"calibration")
  print("aggregate image pair modeling")
  cal_mss_tc_composite_model(dir)
  
  dir = file.path(dir,"composite_model")
  bcoef = as.numeric(read.csv(file.path(dir,"tcb_cal_composite_coef.csv"))[1,2:6])
  gcoef = as.numeric(read.csv(file.path(dir,"tcg_cal_composite_coef.csv"))[1,2:6])
  wcoef = as.numeric(read.csv(file.path(dir,"tcw_cal_composite_coef.csv"))[1,2:6])
  
  print("MSS image pair ")
  if(cores==2){
    cl = makeCluster(cores)
    registerDoParallel(cl)
    t=proc.time()
    cfun <- function(a, b) NULL
    o = foreach(i=1:length(mssf), .combine="cfun",.packages="LandsatLinkr") %dopar% msssr2tc(mssf[i],bcoef,gcoef,wcoef,"calibrate") #
    stopCluster(cl)
  } else {for(i in 1:length(mssf)){msssr2tc(mssf[i],bcoef,gcoef,wcoef,"calibrate")}}
  
  
  msswrs1dir = file.path(msswrs1dir,"images")
  msswrs2dir = file.path(msswrs2dir,"images")
  
  msswrs1files = list.files(msswrs1dir, "dos_sr.tif", recursive=T, full.names=T)
  msswrs2files = list.files(msswrs2dir, "dos_sr.tif", recursive=T, full.names=T)
  
  files = c(msswrs1files,msswrs2files)
  dir = file.path(msswrs2dir,"calibration","composite_model")
  bcoef = as.numeric(read.csv(file.path(dir,"tcb_cal_composite_coef.csv"))[1,2:6])
  gcoef = as.numeric(read.csv(file.path(dir,"tcg_cal_composite_coef.csv"))[1,2:6])
  wcoef = as.numeric(read.csv(file.path(dir,"tcw_cal_composite_coef.csv"))[1,2:6])
  
  if(cores==2){
    cl = makeCluster(cores)
    registerDoParallel(cl)
    t=proc.time()
    cfun <- function(a, b) NULL
    o = foreach(i=1:length(files), .combine="cfun",.packages="MSScvm") %dopar% msssr2tc(files[i],bcoef,gcoef,wcoef,"apply") #
    stopCluster(cl)
  } else {for(i in 1:length(files)){msssr2tc(files[i],bcoef,gcoef,wcoef,"apply")}}
  
  print(proc.time()-t)
}




