#' GUI for calibrate_and_composite function 
#'
#' GUI for calibrate_and_composite function
#' @export


run_calibrate_and_composite = function(){
  
  choices = c("calibrate and composite", "msscal", "mixel")
  selection = select.list(choices, title = "Select a function to run")
  if(selection == "calibrate and composite"){process = c(1,2)}
  if(selection == "msscal"){process = 1}
  if(selection == "mixel"){process = 2}
  
  msswrs1dir = choose.dir(caption = "Select a MSS WRS-1 scene directory. ex. 'C:/mss/wrs1/036032'")
  msswrs2dir = choose.dir(caption = "Select a MSS WRS-2 scene directory. ex. 'C:/mss/wrs2/034032'")
  tmwrs2dir = choose.dir(caption = "Select a TM WRS-2 scene directory. ex. 'C:/tm/wrs2/034032'")
  
  outdir=NULL
  index=NULL
  runname=NULL
  if(sum(process %in% c(2)) > 0){
    outdir = choose.dir(caption = "Select a directory to write the outputs to. ex. 'C:/composites/wrs2_034032'")
    
    choices = c("tca", "tcb", "tcg", "tcw")
    selection = select.list(choices, title = "Select an index to create composites for")
    index=selection
    
    runname = readline("provide a unique name for the composite series. ex. project1: ")
    
    choices = c("from file")
    selection = select.list(choices, title = "What area do you want to create composites for?")
    if(selection == "from file"){useareafile = choose.files(caption = "Select a 'usearea' file", multi=F)}
  }
  
  choices = c("No", "Yes")
  selection = select.list(choices, title = "Process in parallel using 2 cores?")
  if(selection == "No"){cores = 1}
  if(selection == "Yes"){cores = 2}
  
  calibrate_and_composite(msswrs1dir,msswrs2dir,tmwrs2dir,index,outdir,runname,useareafile,doyears="all",order="sensor_and_doy",overlap="mean", cores=cores, process=process)
  
}