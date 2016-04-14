#' Changes a file's extension
#'
#' Changes a file's extension
#' @param old character. old extension
#' @param new character. new extension
#' @param file character. full path of the file to change extension
#' @export


change_extension = function(old, new, file){
  end = nchar(file)-nchar(old)
  newfile = paste(substr(file,1,end),new,sep="")
  return(newfile)
}


