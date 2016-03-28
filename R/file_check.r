#' Handles file existence checking and overwriting
#'
#' Handles file existence checking and overwriting
#' @param file Filename. of file being worked on 
#' @param output Filename. what output file is it checking for? ("archv.tif", "reflectance.tif", etc)
#' @param overwrite logical. if the output file exists should it be deleted
#' @export

file_check = function(file, output, overwrite){
  bname = basename(file)
  dname = dirname(file) 
  imgid = substr(bname, 1, 16)
  if(output == "archv.tif"){
    #print("archv")
    ppprrrdir = substr(dname,1,(nchar(dname)-6)) 
    search = paste(imgid,"_archv.tif",sep="")
    result = list.files(ppprrrdir, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) > 0 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "reflectance.tif"){ #if running mssdn2refl
    #print("reflectance")
    search = paste(imgid,"_reflectance.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "cloudmask.tif"){ #if running mssunpackr
    #print("cloudmask")
    search = paste(imgid,"_cloudmask.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "dos_sr.tif"){
    #print("dos_sr")
    search = paste(imgid,"_dos_sr.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "dos_sr_30m.tif"){
    #print("dos_sr_30m")
    search = paste(imgid,"_dos_sr_30m.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "cloudmask_30m.tif"){
    #print("cloudmask_30m")
    search = paste(imgid,"_cloudmask_30m.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "tca_30m.tif"){
    #print("tca_30m")
    search = paste(imgid,"_tca_30m.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "tc_30m.tif"){
    #print("tc_30m")
    search = paste(imgid,"_tc_30m.tif",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) == 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "ledaps.tif"){
    #print("ledaps")
    ppprrrdir = substr(dname,1,(nchar(dname)-6)) 
    search = paste(imgid,"_ledaps.tif",sep="")
    result = list.files(ppprrrdir, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) >= 1 & overwrite == T){
      result = list.files(file.path(ppprrrdir,"images"), imgid, recursive=T, full.names=T)
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "l8sr.tif"){
    #print("l8sr")
    ppprrrdir = substr(dname,1,(nchar(dname)-6)) 
    search = paste(imgid,"_l8sr.tif",sep="")
    result = list.files(ppprrrdir, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) >= 1 & overwrite == T){
      result = list.files(file.path(ppprrrdir,"images"), imgid, recursive=T, full.names=T)
      unlink(result)
      return(2)
    } else {return(0)}
  } else if(output == "l8sr_tc.tif"){
    #print("l8sr_tc")
    search = paste(imgid,"_tc",sep="")
    result = list.files(dname, search, recursive=T, full.names=T)
    if(length(result) == 0){return(1)} else if(length(result) >= 1 & overwrite == T){
      unlink(result)
      return(2)
    } else {return(0)}
  }
}


