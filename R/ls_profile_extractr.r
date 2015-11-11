library(raster)

writeplot = function(file,date,tcb,tcg,tcw,tca,first=F,last=F){
  start='var data = {"PlotID": 1, "Values": ['
  int=', '
  end=']}'
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

chronoInfo = function(dir,file,outpng){
  tcbstack = list.files(dir, "tcb_composite_stack.bsq", recursive = T, full.names = T)
  tcgstack = list.files(dir, "tcg_composite_stack.bsq", recursive = T, full.names = T)
  tcwstack = list.files(dir, "tcw_composite_stack.bsq", recursive = T, full.names = T)
  tcastack = list.files(dir, "tca_composite_stack.bsq", recursive = T, full.names = T)
  
  files = list.files(file.path(dir,"tcb"), "composite.bsq")
  years = sort(unique(substr(files,1,4)))
  
  coords = data.frame(x,y)
  
  tcbb = extend(brick(tcbstack),130)
  tcgb = extend(brick(tcgstack),130)
  tcwb = extend(brick(tcwstack),130)
  tcab = extend(brick(tcastack),130)
  
  tcbv = extract(tcbb, coords)
  tcgv = extract(tcgb, coords)
  tcwv = extract(tcwb, coords)
  tcav = extract(tcab, coords)
  
  #for making year composites
  len = length(years)
  for(i in 1:len){
    print(paste(i,"/",len,sep=""))
    thisyear=as.numeric(years[i])
    if(i == 1){writeplot(file,thisyear,tcbv[i],tcgv[i],tcwv[i],tcav[i], first=T, last=F)}
    if(i != 1 & i != len){writeplot(file,thisyear,tcbv[i],tcgv[i],tcwv[i],tcav[i], first=F, last=F)}
    if(i == len){writeplot(file,thisyear,tcbv[i],tcgv[i],tcwv[i],tcav[i], first=F, last=T)}
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
  
  png(outpng, width = 255, height=nrow(s))
  plotRGB(s,r=1, g=2, b=3, ext=NULL)
  dev.off()  
}



################################################################

x=-1175187.0000
y= 2533375.0000
dir = "L:/composites/test10"
file = "L:/composites/test10/test10.js"
outpng = "L:/composites/test10/chipstrip.png"

chronoInfo(dir,file,outpng)

