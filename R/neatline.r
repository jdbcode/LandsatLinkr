#' Segment spectral profiles into lines of best fit 
#'
#' Segment spectral profiles into lines of best fit 
#' @param info data.frame. chunk information
#' @param tca raster. tasseled cap angle raster
#' @param tcb raster. tasseled cap brightness raster
#' @param tcg raster. tasseled cap greenness raster
#' @param tcw raster. tasseled cap wetness raster
#' @param outdir character. full path to desired output directory
#' @import ecp
#' @import zoo
#' @import raster

neatline = function(info, tca,tcb,tcg,tcw, outdir, bname){
  
  #define neatline functions
  #despike
  despike = function(value, len){
    for(i in 1:1){
      #value = tcgv
      #before = c(NA, value[-1] - value[-length(value)])*-1
      #after = c(value[-1] - value[-length(value)], NA)
      before = c(NA, value[-1] - value[-len])*-1
      after = c(value[-1] - value[-len], NA)
      absba = rowSums(cbind(abs(before),abs(after)), na.rm=T)
      difba = abs(rowSums(cbind(before,after*-1), na.rm=T))
      ratio1 = absba/median(absba)
      ratio2 = difba/median(difba)
      these = which(ratio1 > 4.5 & ratio2 < 4)
      #adj = rowMeans(cbind(before[these],after[these]))
      #value[these] = adj+value[these]
      value[these] = NA
    }
    return(value)
  }
  
  fitit = function(value,year,estimate, breaks, len){    
    df = data.frame(value,year,fitted=NA)  
    for(i in 1:(len-1)){
      these = seq(estimate[i], estimate[i+1])
      newdf = df[these,]
      line = lm(value[these] ~ year[these])
      v = fitted(line)
      breaks[i] = mean(c(breaks[i], v[1]),na.rm=T)
      breaks[i+1] = v[length(these)]
    }
    for(i in 1:(len-1)){
      value = breaks[c(i,i+1)]
      year = df$year[c(estimate[i],estimate[i+1])]
      these = seq(estimate[i], estimate[i+1])
      line = lm(value ~ year)
      v = predict(line, df[these,])
      df$fitted[these] = v
    }
    return(df)
  }
  
  writeout = function(values, template, index, outdir, info){
    values = setValues(template,values)
    values = as(values, "SpatialGridDataFrame")
    outfile = file.path(outdir,paste(index,"_chunk_",info$chunk,".tif",sep=""))
    writeGDAL(values, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  }
  
  merge_fitted_composite = function(outdir,index, bname){
    files = sort(list.files(outdir, index, full.names=T))
    
    if(index == "tca"){bname = bname[1]}
    if(index == "tcb"){bname = bname[2]}
    if(index == "tcg"){bname = bname[3]}
    if(index == "tcw"){bname = bname[4]}
    
    for(m in 1:length(files)){
      if(m == 1){mergeit = "r1"} else {mergeit = paste(mergeit,",r",m, sep="")}
      dothis = paste("r",m,"=brick(files[",m,"])", sep="")
      eval(parse(text=dothis))
      if(m == length(files)){mergeit = paste("newimg = merge(",mergeit,")", sep="")}
    }
    if(length(files) == 1){newimg = r1} else {eval(parse(text=mergeit))}
    
    newimg = as(newimg, "SpatialGridDataFrame")
    #outfile = file.path(outdir,paste(index,"_fitted_stack.bsq",sep=""))
    outfile = file.path(outdir,bname)
    writeGDAL(newimg, outfile, drivername = "ENVI", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  }
  
  #start the goods
  tca = template = crop(tca, extent(tca, info$from, info$to, 1, info$n_col))
  tca = getValues(tca)
  tcb = getValues(crop(tcb, extent(tcb, info$from, info$to, 1, info$n_col)))
  tcg = getValues(crop(tcg, extent(tcg, info$from, info$to, 1, info$n_col)))
  tcw = getValues(crop(tcw, extent(tcw, info$from, info$to, 1, info$n_col)))
  
  for(c in 1:nrow(tca)){
    tcav = as.vector(tca[c,]) 
    tcbv = as.vector(tcb[c,])
    tcgv = as.vector(tcg[c,])
    tcwv = as.vector(tcw[c,])
    
    tcav[tcav == 0] = NA
    tcbv[tcbv == 0] = NA
    tcgv[tcgv == 0] = NA
    tcwv[tcwv == 0] = NA
    
    check = length(which(is.na(tcav) == F))
    if(check > 6){
      
      len = length(tcav)
      tcad=despike(tcav, len)
      tcbd=despike(tcbv, len)
      tcgd=despike(tcgv, len)
      tcwd=despike(tcwv, len)
      
      tcad = na.approx(tcad, na.rm = F)
      tcbd = na.approx(tcbd, na.rm = F)
      tcgd = na.approx(tcgd, na.rm = F)
      tcwd = na.approx(tcwd, na.rm = F)
      
      nona = which(is.na(tcad) == F)
      first = nona[1]
      last = nona[length(nona)]
      
      mx = cbind(tcad,tcbd,tcgd,tcwd)
      
      #breakpoints
      y1 = e.divisive(X=mx,sig.lvl=0.05,R=75,k=NULL,min.size=2,alpha=1)
      estimate = y1$estimate
      len = length(estimate)
      #estimate[len] = estimate[len]-1
      estimate[1] = first
      estimate[len] = last
      
      
      #add fast bp
#       before = c(NA, tcad[-1] - tcad[-length(tcad)])*-1
#       for(i in 1:3){
#         befm = median(before[before > 0], na.rm=T)
#         #aftm = median(after[after < 0], na.rm=T)
#         these = which(before[estimate] > befm)
#         estimate = c(estimate, estimate[these]-1)
#       }
#       estimate = sort(unique(estimate))
#       
     before = abs(c(NA, tcad[-1] - tcad[-length(tcad)]))  
     befm = median(before,na.rm=T)
      
      for(i in 1:3){
        #i=1
        #befm = median(before[before > 0], na.rm=T)
        #aftm = median(after[after < 0], na.rm=T)
        these = which(before[estimate] > befm*2.6) #1.5
        estimate = c(estimate, estimate[these]-1)
      }
      estimate = sort(unique(estimate))
      

      len = length(estimate)
      breaks = array(NA, len)
      
      df = fitit(tcad,yearseries,estimate, breaks, len)
      tca[c,] = round(df$fitted)
      
      df = fitit(tcbd,yearseries,estimate, breaks, len)
      tcb[c,] = round(df$fitted)
      
      df = fitit(tcgd,yearseries,estimate, breaks, len)
      tcg[c,] = round(df$fitted)
      
      df = fitit(tcwd,yearseries,estimate, breaks, len)
      tcw[c,] = round(df$fitted)
      
    } else {
      tca[c,] = NA
      tcb[c,] = NA
      tcg[c,] = NA
      tcw[c,] = NA
    }
    
  }
  
  writeout(tca, template, "tca", outdir, info)
  writeout(tcb, template, "tcb", outdir, info)
  writeout(tcg, template, "tcg", outdir, info)
  writeout(tcw, template, "tcw", outdir, info)
  
  merge_fitted_composite(outdir, "tca")
  merge_fitted_composite(outdir, "tcb")
  merge_fitted_composite(outdir, "tcg")
  merge_fitted_composite(outdir, "tcw")

}

