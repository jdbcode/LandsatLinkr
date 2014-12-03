#' GUI for prepare_images function 
#'
#' GUI for prepare_images function
#' @export


run_prepare_images = function(){
  
  choices = c("prepare mss", "prepare tm", "mssunpackr", "msswarp", "mssdn2rad", "mssatcor", "msscvm", "tmunpackr")
  selection = select.list(choices, title = "Select a function to run")
  if(selection == "prepare mss"){process = seq(1:5)}
  if(selection == "prepare tm"){process = 6}
  if(selection == "mssunpackr"){process = 1}
  if(selection == "msswarp"){process = 2}
  if(selection == "mssdn2rad"){process = 3}
  if(selection == "mssatcor"){process = 4}
  if(selection == "msscvm"){process = 5}
  if(selection == "tmunpackr"){process = 6}
  
  scenedir = choose.dir(caption = "Select a scene directory. ex. 'C:/mss/036032'")
  
  choices = c("No", "Yes")
  selection = select.list(choices, title = "Process in parallel using 2 cores?")
  if(selection == "No"){cores = 1}
  if(selection == "Yes"){cores = 2}
  
  reso = 60 #set default
  proj = "albers" #set default
  if(sum(process %in% c(1,6)) > 0){
    choices = c("30 meter", "60 meter")
    selection = select.list(choices, title = "Select a pixel resolution to use")
    if(selection == "30 meter"){reso = 30}
    if(selection == "60 meter"){reso = 60}
    
    choices = c("NAD83 UTM Default", "USGS North American Albers")
    selection = select.list(choices, title = "Select a map projection to use")
    if(selection == "NAD83 UTM Default"){proj = "default"}
    if(selection == "USGS North American Albers"){proj = "albers"}
  }
  
  demfile=NULL
  if(all(is.na(match(process,5))) == F){
    demfile = choose.files(caption = "Select a scene corresponding DEM file", multi=F)
  }
  
  prepare_images(scenedir, demfile, proj=proj, reso=reso, process=process,cores=cores)
  
}
