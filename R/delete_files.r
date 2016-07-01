#' Delete files
#'
#' Delete all files associated with a particular image ID
#' @param file character. any file associated with a particular image ID
#' @param reason numeric or character. a reason for deleting the image.
#' @export


delete_files = function(file,reason){
  if(regexpr(".tar.gz", file)[1] == -1){
    dir = substr(dirname(file),1,nchar(dirname(file))-12)
  } else {dir = substr(dirname(file),1,nchar(dirname(file))-6)} 
      
  imgid = substr(basename(file),1,16)
  files = list.files(path=dir, pattern=imgid, full.names=T, recursive=T)
  if(reason == 1){reason = "L1G"}
  if(reason == 2){reason = "Could not find enough geowarp tie-point"}
  if(reason == 3){reason = "Poor TC regression"}
  if(reason == 4){reason = "Poor geowarping"}
  if(reason == 5){reason = "Empty band(s)"}
  if(reason == 6){reason = "MSS year >= 1995: missing band 4"}
  outfile = file.path(dir,"images_deleted",paste(imgid,"_delete_record.csv",sep=""))
  outdir = dirname(outfile)
  outdirtargz = file.path(outdir,"targz")
  dir.create(outdir, recursive=T, showWarnings=F)
  dir.create(outdirtargz, recursive=T, showWarnings=F)
  write(c(imgid,reason), file=outfile, ncolumns=2, sep=",")
  targz = grep(".tar.gz", files, value=T)
  if(length(targz) == 1){
    newtargz = file.path(outdirtargz, basename(targz))
    file.rename(targz,newtargz)
  }
  unlink(files, force=T)
}
