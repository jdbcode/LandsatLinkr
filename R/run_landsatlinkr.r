#' GUI for running LandsatLinkr
#'
#' GUI for running LandsatLinkr
#' @export

run_landsatlinkr = function(){
  
  #################################################################################################################
  #get the process to run
  
  correct = "No"
  while(correct == "No"){
    choices = c("Prepare MSS",
                "Prepare TM/ETM+",
                "Prepare OLI",
                "Calibrate MSS to TM",
                "Calibrate OLI to ETM+",
                "Composite imagery")#,
                #"Fit neat lines")
    
    selection = select.list(choices, title = "Select a process to run")
    if(selection == "Prepare MSS"){process = seq(1:5)}
    if(selection == "Prepare TM/ETM+"){process = 6}
    if(selection == "Prepare OLI"){process = 7}
    if(selection == "Calibrate MSS to TM"){process = 8}
    if(selection == "Calibrate OLI to ETM+"){process = 9}
    if(selection == "Composite imagery"){process = 10}
    
    cat(paste("You have selected:",selection,"\n"))
    correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
    if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
  }
  
  #################################################################################################################
  #define the number of cores to use
  
  if(sum(process %in% seq(1,10)) > 0){  
    
    correct = "No"
    while(correct == "No"){
      choices = c("Yes", "No")
      selection = select.list(choices, title = "Process in parallel using 2 cores when possible?")
      if(selection == "Yes"){cores = 2}
      if(selection == "No"){cores = 1}
      cat(paste("You have selected:",selection,"\n"))
      correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
      if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
    }
  }
  
  
  if(sum(process %in% seq(1,10)) > 0){  
    
    correct = "No"
    while(correct == "No"){
      choices = c("Yes", "No")
      selection = select.list(choices, title = "If files from a process already exist, should they be overwritten?")
      if(selection == "Yes"){overwrite = T}
      if(selection == "No"){overwrite = F}
      cat(paste("You have selected:",selection,"\n"))
      correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
      if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
    }
  }
  
  #################################################################################################################
  #prepare mss and tm/etm+ data
  #################################################################################################################
  
  #################################################################################################################
  #select a directory to process
  
  if(sum(process %in% seq(1,7)) > 0){
    correct = "No"
    while(correct == "No"){
      if(sum(process %in% seq(1,5)) > 0){
        cat("please select an MSS directory to process\n\n")
        scenedir = choose.dir(caption = "Select an MSS scene directory. ex. 'C:/mss/wrs1/036032'")
      }
      if(sum(process %in% 6) > 0){
        cat("please select a TM/ETM+ directory to process\n\n")
        scenedir = choose.dir(caption = "Select an TM/ETM+ scene directory. ex. 'C:/tm/wrs2/036032'")
      }
      if(sum(process %in% 7) > 0){
        cat("please select a OLI directory to process\n\n")
        scenedir = choose.dir(caption = "Select an OLI scene directory. ex. 'C:/oli/wrs2/036032'")
      }
      
      cat(paste("You have selected:",scenedir,"\n"))
      correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
      if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      
      # try to figure out if this is a possible scenedir
      checkdir = basename(scenedir)
      
      # did a targz folder get selected
      if(checkdir == 'targz'){
        cat(paste0("It looks like you selected a scene's 'targz' folder.\n",
                   "Please select the scene head folder - something like: 'C:/tm/wrs2/036032'\n\n"))
        correct = "No"
      } else if(nchar(checkdir) != 6){
        cat(paste0("LLR is expecting a folder with 6 characters.\n",
                   "Please select the scene head folder - something like: 'C:/tm/wrs2/036032'\n\n"))
        correct = "No"
      } else if(is.na(as.numeric(checkdir))){
        cat(paste0("LLR is expecting a folder with all numeric characters.\n",
                   "Please select the scene head folder - something like: 'C:/tm/wrs2/036032'\n\n"))
        correct = "No"
      }
    }
    
    #################################################################################################################
    #define a projection
    
    correct = "No"
    while(correct == "No"){  
      cat("Please define a projection to use for this project\n")
      cat("If you need assistance, see the 'Selecting a geographic projection' section in the user manual\n")
      cat("If you want to exit, press the 'return' key and then select 'Exist'\n")

      proj = readline("  Provide PROJ.4 string: ") #http://spatialreference.org/

      cat(paste("You have selected:",proj,"\n"))
      correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
      if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
    }
    
    #################################################################################################################
    #provide a DEM
    
    if(sum(process %in% seq(1,5)) > 0){
      correct = "No"
      while(correct == "No"){
        cat("Please select a scene corresponding DEM file\n")
        demfile = choose.files(caption = "Select a scene corresponding DEM file", multi=F)
        cat(paste("You have selected:",demfile,"\n"))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
    }
    prepare_images(scenedir, demfile, proj=proj, process=process,cores=cores, overwrite=overwrite)
  }
  

  #################################################################################################################
  #msscal and mixel
  #################################################################################################################

  if(sum(process %in% 8:10) > 0){
    
    #################################################################################################################
    #select directories for calibration
    if(sum(process %in% 8) > 0){
      cat("Beginning scene selection for MSS to TM calibration\n")
      cat("Please see the 'Running LLR Step 3 - MSS to TM calibration' section of the user guide for assistance\n")
      cat("-----------------------------------------------------------------------\n")
      correct = "No"
      while(correct == "No"){
        cat("Select a MSS WRS-1 scene directory.\n")
        msswrs1dir = choose.dir(caption = "Select a MSS WRS-1 scene directory. ex. 'C:/mss/wrs1/036032'")
        cat("Select a MSS WRS-2 scene directory.\n")
        msswrs2dir = choose.dir(caption = "Select a MSS WRS-2 scene directory. ex. 'C:/mss/wrs2/034032'")
        cat("Select a TM/ETM+ WRS-2 scene directory.\n")
        tmwrs2dir = choose.dir(caption = "Select a TM/ETM+ WRS-2 scene directory. ex. 'C:/tm/wrs2/034032'")
        oliwrs2dir=NULL
        outdir=NULL
        index=NULL
        runname=NULL
        calcomprocess = 1
        
        cat("-----------------------------------------------------------------------\n")
        cat(paste("You have selected MSS WRS-1 scene:",msswrs1dir,'\n'))
        cat(paste("You have selected MSS WRS-2 scene:",msswrs2dir,'\n'))
        cat(paste("You have selected TM/ETM+ WRS-2 scene:",tmwrs2dir,'\n'))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
    }
    
    #################################################################################################################
    #select directories for calibration
    if(sum(process %in% 9) > 0){
      cat("Beginning scene selection for OLI to ETM+ calibration\n")
      cat("Please see the 'Running LLR Step 3 - OLI to ETM+ calibration' section of the user guide for assistance\n")
      cat("-----------------------------------------------------------------------\n")
      correct = "No"
      while(correct == "No"){
        cat("Select a OLI WRS-2 scene directory.\n")
        oliwrs2dir = choose.dir(caption = "Select a OLI WRS-2 scene directory. ex. 'C:/oli/wrs1/036032'")
        cat("Select a TM/ETM+ WRS-2 scene directory.\n")
        tmwrs2dir = choose.dir(caption = "Select a TM/ETM+ WRS-2 scene directory. ex. 'C:/tm/wrs2/034032'")
        msswrs1dir=NULL
        msswrs2dir=NULL
        outdir=NULL
        index=NULL
        runname=NULL
        calcomprocess = 2
        
        cat("-----------------------------------------------------------------------\n")
        cat(paste("You have selected OLI WRS-2 scene:",oliwrs2dir,'\n'))
        cat(paste("You have selected TM/ETM+ WRS-2 scene:",tmwrs2dir,'\n'))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
    }
    
    
    #################################################################################################################
    #select directories for mosaicking
    
    if(sum(process %in% 10) > 0){
    
      msswrs1dir=NA
      msswrs2dir=NA
      tmwrs2dir=NA
      oliwrs2dir=NA
    
      correct = "No"
      while(correct == "No"){
        choices = c("Yes", "No")
        

        doType = select.list(c("Yes","No","Exit"), title = "Are there MSS WRS-1 images you would like to include in composites?")
        if(doType == "Exit"){return(cat("Stopping LLR","\n\n"))} else 
        if(doType == "Yes"){
          msswrs1dir = choose.dir(caption = "Select a MSS WRS-1 scene directory. ex. 'C:/mss/wrs1/036032'")
          answer = "Yes"
          while(answer == "Yes"){
            cat("Here is what you have so far - MSS WRS-1 scene:\n")
            print(msswrs1dir)
            answer = select.list(choices, title = "Is there another MSS WRS-1 scene directory to add?")
            if(answer == "Yes"){msswrs1dir = c(msswrs1dir, choose.dir(caption = "Select a MSS WRS-1 scene directory. ex. 'C:/mss/wrs1/036032'"))}
          }
        }

        

        doType = select.list(c("Yes","No","Exit"), title = "Are there MSS WRS-2 images you would like to include in composites?")
        if(doType == "Exit"){return(cat("Stopping LLR","\n\n"))} else 
        if(doType == "Yes"){        
          msswrs2dir = choose.dir(caption = "Select a MSS WRS-2 scene directory. ex. 'C:/mss/wrs2/034032'")
          answer = "Yes"
          while(answer == "Yes"){
            cat("Here is what you have so far - MSS WRS-2 scene:\n")
            print(msswrs2dir)
            answer = select.list(choices, title = "Is there another MSS WRS-2 scene directory to add?")
            if(answer == "Yes"){msswrs2dir = c(msswrs2dir, choose.dir(caption = "Select a MSS WRS-2 scene directory. ex. 'C:/mss/wrs2/034032'"))}
          }
        }
        
        
        doType = select.list(c("Yes","No","Exit"), title = "Are there TM/ETM+ WRS-2 images you would like to include in composites?")
        if(doType == "Exit"){return(cat("Stopping LLR","\n\n"))} else 
        if(doType == "Yes"){            
          tmwrs2dir = choose.dir(caption = "Select a TM/ETM+ WRS-2 scene directory. ex. 'C:/tm/wrs2/034032'")
          answer = "Yes"
          while(answer == "Yes"){
            cat("Here is what you have so far - TM/ETM+ WRS-2 scene:\n")
            print(tmwrs2dir)
            answer = select.list(choices, title = "Is there another TM/ETM+ WRS-2 scene directory to add?")
            if(answer == "Yes"){tmwrs2dir = c(tmwrs2dir, choose.dir(caption = "Select a TM/ETM+ WRS-2 scene directory. ex. 'C:/tm/wrs2/034032'"))}
          }
        }
        
        
        doType = select.list(c("Yes","No","Exit"), title = "Are there OLI WRS-2 images you would like to include in composites?")
        if(doType == "Exit"){return(cat("Stopping LLR","\n\n"))} else 
        if(doType == "Yes"){  
          oliwrs2dir = choose.dir(caption = "Select a OLI WRS-2 scene directory. ex. 'C:/oli/wrs2/034032'")
          answer = "Yes"
          while(answer == "Yes"){
            cat("Here is what you have so far - OLI WRS-2 scene:\n")
            print(oliwrs2dir)
            answer = select.list(choices, title = "Is there another OLI WRS-2 scene directory to add?")
            if(answer == "Yes"){oliwrs2dir = c(oliwrs2dir, choose.dir(caption = "Select a OLI WRS-2 scene directory. ex. 'C:/oli/wrs2/034032'"))}
          } 
        }
        
        cat(paste("You have selected MSS WRS-1 scene:",msswrs1dir,'\n'))
        cat(paste("You have selected MSS WRS-2 scene:",msswrs2dir,'\n'))
        cat(paste("You have selected TM/ETM+ WRS-2 scene:",tmwrs2dir,'\n'))
        cat(paste("You have selected OLI WRS-2 scene:",oliwrs2dir,'\n'))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
      
      correct = "No"
      while(correct == "No"){
        cat("Please select a directory to write the outputs to. ex. 'C:/composites/study_area'\n")
        outdir = choose.dir(caption = "Select a directory to write the outputs to. ex. 'C:/composites/study_area'")
        
        cat(paste("You have selected:",outdir,"\n"))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
      

      correct = "No"
      while(correct == "No"){
        runname = readline("Provide a unique name for the composite series. ex. project1: ")
        cat(paste("You have selected:",runname,"\n"))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
      
      #################################################################################################################
      #use area file or coordinates
      

      
      
      choices = c("From file",
                  "Provide coordinates")
      selection = select.list(choices, title = "What area do you want to create composites for?")
      
      #find a files to get the project reference
      dirToCheck = c(msswrs1dir[1],msswrs2dir[1],tmwrs2dir[1],oliwrs2dir[1])
      projfiles = list.files(path = file.path(dirToCheck[1], "images"), pattern=".tif$", full.names=T, recursive=T)
      if(length(projfiles) == 0){projfiles = list.files(path = file.path(dirToCheck[2], "images"), pattern=".tif$", full.names=T, recursive=T)}
      if(length(projfiles) == 0){projfiles = list.files(path = file.path(dirToCheck[3], "images"), pattern=".tif$", full.names=T, recursive=T)} 
      if(length(projfiles) == 0){projfiles = list.files(path = file.path(dirToCheck[4], "images"), pattern=".tif$", full.names=T, recursive=T)}
      projfile = projfiles[1]
      
      #deal with the selection
      if(selection == "From file"){
        correct = "No"
        while(correct == "No"){
          useareafile = choose.files(caption = "Select a 'usearea' file", multi=F)
          cat(paste("You have selected:",useareafile,"\n"))
          correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
          if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
          
          make_usearea_file_bsq(useareafile, projfile)
        }
      } else if(selection == "Provide coordinates"){
        
        correct = "No"
        while(correct == "No"){
          check = 0
          while(check != 2){
            cat("Please provide min and max xy coordinates (in image set projection units) defining a study area that intersects the image set\n")
            xmx = as.numeric(readline("max x coordinate: "))
            xmn = as.numeric(readline("min x coordinate: "))
            ymx = as.numeric(readline("max y coordinate: "))
            ymn = as.numeric(readline("min y coordinate: "))
            check = (xmx > xmn) + (ymx > ymn)
            if(check != 2){cat("error - coordinates do not create a square - try again\n")}
          }
          cat("You have selected:\n")
          cat(paste("max x coordinate:",xmx,"\n"))
          cat(paste("min x coordinate:",xmn,"\n"))
          cat(paste("max y coordinate:",ymx,"\n"))
          cat(paste("min y coordinate:",ymn,"\n"))
          correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
          if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
        }
        useareafile = file.path(outdir, paste(runname,"_usearea.tif",sep=""))
        make_usearea_file(c(msswrs1dir[1],msswrs2dir[1],tmwrs2dir[1],oliwrs2dir[1]), useareafile, xmx, xmn, ymx, ymn)
      }
      
      #dates
      cat("\n\n")
      cat("The following two parameter inputs are the minimum and maximum day-of-year to include in your composites.\nExample: June 15th through August 31st would be 166 and 243.\nFor more information see the LLR guide page: http://landsatlinkr.jdbcode.com/guide.html#doy")
      cat("\n\n")
      correct = "No"
      while(correct == "No"){
        startday = as.numeric(readline("Define a minimum day-of-year to include in the composites: "))
        cat("You have selected:",startday,"\n")
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
      
      correct = "No"
      while(correct == "No"){
        endday = as.numeric(readline("Define a maximum day-of-year to include in the composites: "))
        cat("You have selected:",endday,"\n")
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
      }
      
      
      #figure out how to name the files
      yearadj = 0
      if((endday-startday) < 0){
        cat("\n")
        correct = "No"
        while(correct == "No"){
          yearadjanswer = select.list(c("Pre","Post"), title = "Your date range crosses the year divide, would you like to label the annual composite files with the pre- or post- divide year?")
          cat("You have selected:",yearadjanswer,"\n")
          correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
          if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
          if(yearadjanswer == "Post"){yearadj = 1}
        }
      }
      
      #overlap methods
      correct = "No"
      while(correct == "No"){
        choices = c("Mean",
                    "Maximum",
                    "Minimum",
                    "Median *long processing time*")#,
        #"Fit neat lines")
        
        overlap = select.list(choices, title = "Select a method to summarize the value of overlapping pixels")
        
        cat(paste("You have selected:",overlap,"\n"))
        correct = select.list(c("Yes","No","Exit"), title = "Is that correct?")
        if(correct == "Exit"){return(cat("Stopping LLR","\n\n"))}
        
        if(overlap == "Mean"){overlap = "mean"}
        if(overlap == "Maximum"){overlap = "max"}
        if(overlap == "Minimum"){overlap = "min"}
        if(overlap == "Median *long processing time*"){overlap = "median"}
        
      }
      calcomprocess = 3
    }
    
    #################################################################################################################
    #execute the calibrate_and_composite function
    
    calibrate_and_composite(msswrs1dir,msswrs2dir,tmwrs2dir,oliwrs2dir,index="all",outdir,runname,useareafile,doyears="all",order="none",overlap=overlap, cores=cores, process=calcomprocess, overwrite=overwrite, startday=startday, endday=endday, yearadj=yearadj)  #overlap="mean"
  }
}
