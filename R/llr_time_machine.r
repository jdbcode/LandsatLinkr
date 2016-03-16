#' Decompress, stack, and reproject LPSG MSS images
#'
#' Decompresses, stacks, and optionally reprojects LPGS MSS images recieved from USGS EROS as .tar.gz files
#' @param imgdir direcory path. full path to directory containing LandsatLinkr annual composites 
#' @param outdir direcory path. full path to the directory where you want LLR-TimeMachine data to be written
#' @param coordfile csv file path. full path to a comma delimited file containing the plot number, x, and y coordinates pixels you want to view in LLR-TimeMachine
#' @import raster
#' @export



llr_time_machine = function(imgdir,outdir,coordfile){
  
  ##########################################################################################
  writeplot = function(file,plot,gmBounds,date,tcb,tcg,tcw,tca, first=F,last=F, end=F, finalPlot=F){
    chipstrip = paste('"imgs/plot_',plot,'_chipstrip.png",',sep="")
    start=paste('{"plotID": ',plot,',','"chipStrip":',chipstrip,'"LatLon":[',gmBounds[1,2],',',gmBounds[1,1],'],','"bounds":[',gmBounds[2,2],',',gmBounds[3,2],',',gmBounds[3,1],',',gmBounds[2,1],'],','"Values": [')
    int=', '
    if(end == T & finalPlot == F){end=']},'} else if(end == T & finalPlot == T){end=']}'}
    imginfo = paste(
      paste('{"Year": ', date,',',sep=""),
      paste(' "TCB": ', tcb,',',sep=""),
      paste(' "TCG": ', tcg,',',sep=""),
      paste(' "TCW": ', tcw,',',sep=""),
      paste(' "TCA": ', tca,sep=""),
      '}',
      sep="")
    if(first == T){
      write(start, file, append=TRUE)
      write(paste(imginfo,int,sep=""), file, append=TRUE)
    }
    if(first != T & last != T){write(paste(imginfo,int,sep=""), file, append=TRUE)}
    if(last == T){
      write(imginfo, file, append=TRUE)
      write(end, file, append=TRUE)
    }
  }
  
  makeChipStrip = function(imgStack){ 
    singleShift = xres(imgStack)*255
    for(i in 1:nlayers(imgStack)){
      if(i == 1){
        fullSeries = subset(imgStack,i)
      } else{
        shifted = shift(subset(imgStack,i),x=0,y=singleShift*(i-1)*-1)
        fullSeries = merge(fullSeries,shifted)
      }
    }
    return(fullSeries)
  }
  
  stretchBGW = function(template, index){
    img = as.matrix(template)
    n_stdev = 2
    bmin = 3098-(1247*n_stdev)
    bmax = 3098+(1247*n_stdev)
    gmin = 1549-(799*n_stdev)
    gmax = 1549+(799*n_stdev)
    wmin = -701-(772*n_stdev)
    wmax = -701+(772*n_stdev)
    
    if(index == "tcb"){
      these = which(img < bmin)
      if(length(these != 0)){img[these] = bmin}
      these = which(img > bmax)
      if(length(these != 0)){img[these] = bmax}  
      img = ((img-bmin)/(bmax-bmin))*255
    }
    if(index == "tcg"){
      these = which(img < gmin)
      if(length(these != 0)){img[these] = gmin}
      these = which(img > gmax)
      if(length(these != 0)){img[these] = gmax}
      img = ((img-gmin)/(gmax-gmin))*255
      
    }
    if(index == "tcw"){
      these = which(img < wmin)
      if(length(these != 0)){img[these] = wmin}
      these = which(img > wmax)
      if(length(these != 0)){img[these] = wmax}
      img = ((img-wmin)/(wmax-wmin))*255
    }
    img = setValues(template, round(img))
    return(img)
  }
  
  
  chronoInfo = function(imgdir, outdir, plot, x, y, tcbb, tcgb, tcwb, tcab, years, jsfile, finalPlot=F){
    coords = data.frame(x,y)
    reso = xres(tcbb)
    half = reso*127.5
    ul = matrix(c(x-half, y+half), nrow=1)
    lr = matrix(c(x+half, y-half), nrow=1)
    proj = projection(tcbb)
    spCoordProj = SpatialPoints(matrix(c(x, y), nrow=1),proj4string=CRS(proj))
    SpCoordLonLat = spTransform(spCoordProj, CRS("+init=epsg:4269"))
    SpCoordLonLat = coordinates(SpCoordLonLat)
    spUlProj = SpatialPoints(ul,proj4string=CRS(proj))
    SpUlLonLat = spTransform(spUlProj, CRS("+init=epsg:4269"))
    SpUlLonLat = coordinates(SpUlLonLat)
    spLrProj = SpatialPoints(lr,proj4string=CRS(proj))
    SpLrLonLat = spTransform(spLrProj, CRS("+init=epsg:4269"))
    SpLrLonLat = coordinates(SpLrLonLat)
    gmBounds = matrix(c(SpCoordLonLat,SpUlLonLat,SpLrLonLat),ncol = 2,byrow = T)
    
    tcbv = extract(tcbb, coords)
    tcgv = extract(tcgb, coords)
    tcwv = extract(tcwb, coords)
    tcav = extract(tcab, coords)
    
    #for making year composites
    len = length(years)
    for(y in 1:len){
      thisyear=as.numeric(years[y])
      if(y == 1){writeplot(jsfile,plot,gmBounds,thisyear,tcbv[y],tcgv[y],tcwv[y],tcav[y], first=T, last=F, end=F)}
      if(y != 1 & y != len){writeplot(jsfile,plot,gmBounds,thisyear,tcbv[y],tcgv[y],tcwv[y],tcav[y], first=F, last=F, end=F)}
      if(y == len){writeplot(jsfile,plot,gmBounds,thisyear,tcbv[y],tcgv[y],tcwv[y],tcav[y], first=F, last=T, end=T, finalPlot=finalPlot)}
    }
    
    thisCell = cellFromXY(tcbb,coords)
    thisRowCol = rowColFromCell(tcbb,thisCell)
    e = extent(tcbb, thisRowCol[1,1]-127,thisRowCol[1,1]+127,thisRowCol[1,2]-127,thisRowCol[1,2]+127)
    
    tcbv = crop(tcbb, e)
    tcgv = crop(tcgb, e)
    tcwv = crop(tcwb, e)
    tcav = crop(tcab, e)
    
    tcbvs = stretchBGW(tcbv, "tcb")
    tcgvs = stretchBGW(tcgv, "tcg")
    tcwvs = stretchBGW(tcwv, "tcw")
    
    tcbseries = makeChipStrip(tcbvs)
    tcgseries = makeChipStrip(tcgvs)
    tcwseries = makeChipStrip(tcwvs)
    
    s = stack(tcbseries,tcgseries,tcwseries)
    
    outpng = file.path(outdir,"imgs",paste("plot_",plot,"_chipstrip.png",sep=""))
    png(outpng, width = 255, height=nrow(s))
    plotRGB(s,r=1, g=2, b=3, ext=NULL)
    dev.off()  
  }
  ################################################################################################################################################
  

  #read the csv plot file
  plots = read.csv(coordfile) #plot#,x,y
  
  colnames(plots) = tolower(colnames(plots))
  
  #create the javascript file path
  jsfile = file.path(outdir,"LLR-TimeMachine.js")
  
  #make sure output dirs exist/created
  dir.create(outdir, recursive=T, showWarnings = F)
  dir.create(file.path(outdir,"imgs"), recursive=T, showWarnings = F)
  
  
  tcbstack = list.files(imgdir, "tcb_composite_stack.bsq$", recursive = T, full.names = T)
  tcgstack = list.files(imgdir, "tcg_composite_stack.bsq$", recursive = T, full.names = T)
  tcwstack = list.files(imgdir, "tcw_composite_stack.bsq$", recursive = T, full.names = T)
  tcastack = list.files(imgdir, "tca_composite_stack.bsq$", recursive = T, full.names = T)
  
  n_tcbstack = length(tcbstack)
  n_tcgstack = length(tcgstack)
  n_tcwstack = length(tcwstack)
  n_tcastack = length(tcastack)
  if(sum(n_tcbstack,n_tcgstack,n_tcwstack,n_tcastack) != 4){
    print("Error! You are missing some required files, or the provided directory is wrong.")
    if(n_tcbstack == 0){print(paste("There was no *tcb_composite_stack.bsq file found at this directory:",imgdir,"(recursive search)"))}
    if(n_tcgstack == 0){print(paste("There was no *tcg_composite_stack.bsq file found at this directory:",imgdir,"(recursive search)"))}
    if(n_tcwstack == 0){print(paste("There was no *tcw_composite_stack.bsq file found at this directory:",imgdir,"(recursive search)"))}
    if(n_tcastack == 0){print(paste("There was no *tca_composite_stack.bsq file found at this directory:",imgdir,"(recursive search)"))}
    return
    
  }
  
  
  files = list.files(file.path(imgdir,"tcb"), "composite.bsq")
  years = sort(unique(substr(files,1,4)))
  
  ext = extent(raster(tcbstack))
  tcbb = extend(brick(tcbstack),130)
  tcgb = extend(brick(tcgstack),130)
  tcwb = extend(brick(tcwstack),130)
  tcab = extend(brick(tcastack),130)

  #for each plot create the data
  nplots = nrow(plots)
  for(i in 1:nplots){
    print(paste(i,"/",nplots,sep=""))
    if(i == 1){
      unlink(jsfile)
      write('var allData = [', jsfile, append=TRUE)
    }
    
    #check to make sure the point is not outside the image
    x = plots$x[i]
    y = plots$y[i]
    plot = plots$plotid[i]
    bad = 0
    if(x < ext[1]){
      print(paste("x coordinate for plot:",plot,"is outside the image to the west"))
      bad = 1
    }  
    if(x > ext[2]){
      print(paste("x coordinate for plot:",plot,"is outside the image to the east"))
      bad = 1
    }
    if(y < ext[3]){
      print(paste("y coordinate for plot:",plot,"is outside the image to the south"))
      bad = 1
    }
    if(y > ext[4]){
      print(paste("y coordinate for plot:",plot,"is outside the image to the north"))
      bad = 1
    }
    if(bad == 1){
      print(paste("plot:",plot,"is outside the image, skipping it..."))
      next
    }
    
    #if everything is okay, extract the data for the point
    if(i == nplots){
      chronoInfo(imgdir, outdir, plot, x, y, tcbb, tcgb, tcwb, tcab, years, jsfile, finalPlot=T)
      write(']', jsfile, append=TRUE)
    } else{
      chronoInfo(imgdir, outdir, plot, x, y, tcbb, tcgb, tcwb, tcab, years, jsfile)
    }
  }
}



