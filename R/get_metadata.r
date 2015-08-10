#' Retrieve Landsat image metadata
#' 
#' Uses the image file name to find the corresponding *MTL.txt image metadata file provided with LPSG Landsat images and returns a dataframe with pertinent image information
#' @param file LPGS-processed Landsat image filename (full system path to file) 
#' @return Dataframe with pertinent image information
#' @export

get_metadata = function(file){
  mtlfile = file.path(dirname(file),paste(substr(basename(file),1,17),"MTL.txt", sep=""))
  tbl = unlist(read.delim(mtlfile, header=F, skipNul=T))
  
  #get the path row "ppprrr"
  ppprrr = substr(bname, 4,9)
  
  #get the day-of-year
  doy = as.numeric(substr(bname, 14,16))
  
  #get the year
  year = as.numeric(substr(bname, 10,13))
  
  #get the yearday
  yearday = as.numeric(substr(bname, 10,16))
  
  #get the image id
  imgid = substr(bname, 1, 21)
  
  #get the data type
  string = as.character(grep("DATA_TYPE = ", tbl, value=T))
  pieces = unlist(strsplit(string, " "))
  datatype = pieces[7]
  
  #get the sensor
  string = as.character(grep("SPACECRAFT_ID =", tbl, value=T))
  pieces = unlist(strsplit(string, " "))
  sensor = pieces[7]
  
  #get the sun elevation
  string = as.character(grep("SUN_ELEVATION = ", tbl, value=T))
  pieces = unlist(strsplit(string, " "))
  sunelev = as.numeric(pieces[7])
  sunzen = 90 - sunelev
  
  #get the sun azimuth
  string = as.character(grep("SUN_AZIMUTH = ", tbl, value=T))
  pieces = unlist(strsplit(string, " "))
  sunaz = as.numeric(pieces[7])
  
  # get min and max radiance; gain and bias
  if(sensor == "LANDSAT_1" | sensor == "LANDSAT_2" | sensor == "LANDSAT_3"){
    bands = c(4,5,6,7)} else{bands = c(1,2,3,4)}
  
  maxrad = array(0,4)
  minrad = array(0,4)
  gain = array(0,4)
  bias = array(0,4)
  
  for(i in 1:4){
    maxradsearch = paste("RADIANCE_MAXIMUM_BAND_", bands[i], " =", sep="")
    string = as.character(grep(maxradsearch, tbl, value=T))
    pieces = unlist(strsplit(string, " "))
    maxrad[i] = as.numeric(pieces[7])
    
    minradsearch = paste("RADIANCE_MINIMUM_BAND_", bands[i], " =", sep="")
    string = as.character(grep(minradsearch, tbl, value=T))
    pieces = unlist(strsplit(string, " "))
    minrad[i] = as.numeric(pieces[7])
    
    radmultsearch = paste("RADIANCE_MULT_BAND_", bands[i], " =", sep="")
    string = as.character(grep(radmultsearch, tbl, value=T))
    pieces = unlist(strsplit(string, " "))
    gain[i] = as.numeric(pieces[7])
    
    radaddsearch = paste("RADIANCE_ADD_BAND_", bands[i], " =", sep="")
    string = as.character(grep(radaddsearch, tbl, value=T))
    pieces = unlist(strsplit(string, " "))
    bias[i] = as.numeric(pieces[7])
  } 
  
  #prepare variables for inclusion in output table
  b1minrad = minrad[1]
  b2minrad = minrad[2]
  b3minrad = minrad[3]
  b4minrad = minrad[4]
  b1maxrad = maxrad[1]
  b2maxrad = maxrad[2]
  b3maxrad = maxrad[3]
  b4maxrad = maxrad[4]
  
  b1gain = gain[1]
  b2gain = gain[2]
  b3gain = gain[3]
  b4gain = gain[4]
  b1bias = bias[1]
  b2bias = bias[2]
  b3bias = bias[3]
  b4bias = bias[4]
  
  #get the wrs type
  if(sensor == "LANDSAT_1"){wrstype = "wrs1"}
  if(sensor == "LANDSAT_2"){wrstype = "wrs1"}
  if(sensor == "LANDSAT_3"){wrstype = "wrs1"}
  if(sensor == "LANDSAT_4"){wrstype = "wrs2"}
  if(sensor == "LANDSAT_5"){wrstype = "wrs2"}
  
  #fill in the output table
  df = data.frame(
    ppprrr,
    doy,
    year,
    yearday,
    imgid,
    sensor,
    datatype,
    wrstype,
    sunelev,
    sunzen,
    sunaz,
    b1minrad,
    b2minrad,
    b3minrad,
    b4minrad,
    b1maxrad,
    b2maxrad,
    b3maxrad,
    b4maxrad,
    b1gain,
    b2gain,
    b3gain,
    b4gain,
    b1bias,
    b2bias,
    b3bias,
    b4bias
  )
  
  return(df)
}