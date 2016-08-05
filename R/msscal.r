#' Calibrate MSS images to TM images
#'
#' Calibrate MSS images to TM images using linear regression
#' @param msswrs1dir character. MSS WRS-1 scene directory path
#' @param msswrs2dir character. MSS WRS-2 scene directory path
#' @param tmwrs2dir character. TM WRS-2 scene directory path
#' @param cores numeric. Number of cores to process with options: 1 or 2
#' @export


msscal = function(msswrs1dir, msswrs2dir, tmwrs2dir, cores=2){  
  
  mssfiles = list.files(msswrs2dir, "dos_sr_30m.tif", recursive=T, full.names=T)
  tmfiles = list.files(tmwrs2dir, "tc.tif", recursive=T, full.names=T)
  
  mssimgid = substr(basename(mssfiles),3,16)
  tmimgid = substr(basename(tmfiles),3,16)
  
  thesemss = which(mssimgid %in% tmimgid)
  mssfiles = mssfiles[thesemss]
  
  thesetm = which(tmimgid %in% mssimgid)
  tmfiles = tmfiles[thesetm]
  
  #check for matching mss/tm files to calibrate on
  if(length(mssfiles) == length(tmfiles)){
    mssimgid = substr(basename(mssfiles),3,16)
    tmimgid = substr(basename(tmfiles),3,16)
    mssf = mssfiles[order(mssimgid)]
    tmf = tmfiles[order(tmimgid)]
  } else {stop("There are no matching MSS and TM image dates to use for calibration")}
  
  #check to see if aggregated models have already been made
  tcacheck = list.files(msswrs2dir, "tca_cal_aggregate_coef.csv", recursive=T, full.names=T)
  tcbcheck = list.files(msswrs2dir, "tcb_cal_aggregate_coef.csv", recursive=T, full.names=T)
  tcgcheck = list.files(msswrs2dir, "tcg_cal_aggregate_coef.csv", recursive=T, full.names=T)
  tcwcheck = list.files(msswrs2dir, "tcw_cal_aggregate_coef.csv", recursive=T, full.names=T)
  checks = c(length(tcacheck),length(tcbcheck),length(tcgcheck),length(tcwcheck))
  check = sum(checks)
  if(check > 0){
    print("calibration models have already been create for:")
    if(checks[1]==1){print("...tca")} #don't actually need to model TCA because we make it from orig tc
    if(checks[2]==1){print("...tcb")}
    if(checks[3]==1){print("...tcg")}
    if(checks[4]==1){print("...tcw")}
  }

  if(check != 4){
    print("...single image pair modeling")
    if(cores==2){
      cl = makeCluster(cores)
      registerDoParallel(cl)
      cfun <- function(a, b) NULL
      o = foreach(i=1:length(mssf), .combine="cfun",.packages="LandsatLinkr") %dopar% msscal_single(mssf[i], tmf[i]) #
      stopCluster(cl)
    } else {for(i in 1:length(mssf)){msscal_single(mssf[i], tmf[i])}}
    
    dir = file.path(msswrs2dir,"calibration")
    print("...aggregate image pair modeling")
    cal_mss_tc_aggregate_model(dir)
  }

  msswrs1imgdir = file.path(msswrs1dir,"images")
  msswrs2imgdir = file.path(msswrs2dir,"images")
  
  msswrs1files = list.files(msswrs1imgdir, "dos_sr_30m.tif", recursive=T, full.names=T)
  msswrs2files = list.files(msswrs2imgdir, "dos_sr_30m.tif", recursive=T, full.names=T)
  
  files = c(msswrs1files,msswrs2files)
  dir = file.path(msswrs2dir,"calibration","aggregate_model")
  bcoef = as.numeric(read.csv(file.path(dir,"tcb_cal_aggregate_coef.csv"))[1,2:6])
  gcoef = as.numeric(read.csv(file.path(dir,"tcg_cal_aggregate_coef.csv"))[1,2:6])
  wcoef = as.numeric(read.csv(file.path(dir,"tcw_cal_aggregate_coef.csv"))[1,2:6])
  
  print("...applying model to all MSS images")
  cores = 1
  if(cores==2){
    cl = makeCluster(cores)
    registerDoParallel(cl)
    cfun <- function(a, b) NULL
    o = foreach(i=1:length(files), .combine="cfun",.packages="LandsatLinkr") %dopar% msssr2tc(files[i],bcoef,gcoef,wcoef,"apply") #
    stopCluster(cl)
  } else {for(i in 1:length(files)){msssr2tc(files[i],bcoef,gcoef,wcoef,"apply")}}
  
}




