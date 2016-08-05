#' Create an aggregate oli TC model
#'
#' Create an aggregate oli TC model with diagnostic figures
#' @param dirname character. the directory to the calibration folder
#' @export


cal_oli_tc_aggregate_model = function(dir, overwrite=F){

  
  #make new directory
  outdir = file.path(dir,"aggregate_model")
  dir.create(outdir, recursive=T, showWarnings=F)
  
  #check to see if single cal has already been run
  files = list.files(outdir)
  thesefiles = c("tca_aggregate_mean_dif.png","tca_aggregate_regression.png","tca_cal_aggregate_coef.csv","tca_cal_aggregate_sample.csv","tca_cal_combined_coef.csv","tca_single_mean_dif.png","tca_single_regression.png",
                 "tcb_aggregate_mean_dif.png","tcb_aggregate_regression.png","tcb_cal_aggregate_coef.csv","tcb_cal_aggregate_sample.csv","tcb_cal_combined_coef.csv","tcb_single_mean_dif.png","tcb_single_regression.png",
                 "tcg_aggregate_mean_dif.png","tcg_aggregate_regression.png","tcg_cal_aggregate_coef.csv","tcg_cal_aggregate_sample.csv","tcg_cal_combined_coef.csv","tcg_single_mean_dif.png","tcg_single_regression.png",
                 "tcw_aggregate_mean_dif.png","tcw_aggregate_regression.png","tcw_cal_aggregate_coef.csv","tcw_cal_aggregate_sample.csv","tcw_cal_combined_coef.csv","tcw_single_mean_dif.png","tcw_single_regression.png")
  
  results = rep(NA,length(thesefiles))
  for(i in 1:length(results)){
    test = grep(thesefiles[i], files)
    results[i] = length(test) > 0
  }
  if(all(results) == T & overwrite == F){return(0)}
  
  ###########define some functions############
  plot_multi_cal = function(df,index, model, outfile){
    
    if(index == "tcb"){limits = c(-500, 10000)}
    if(index == "tcg"){limits = c(-500, 5000)}
    if(index == "tcw"){limits = c(-6000, 1000)}
    if(index == "tca"){limits = c(-500, 5000)}
    if(model == "single"){
      
      
      tmp = with(df,by(df, oli_img, function(x) rlm(refsamp ~ singlepred, data = x)))
      tmp = sapply(tmp, coef)
      
      title = paste("single model scatterplot and robust linear regression lines for",index)
      png(outfile, width = 1100, height=800)
      cols = seq(from=2,to=ncol(tmp)+1)
      plot(x=df$singlepred,y=df$refsamp,
           main=title,
           xlab=paste("oli-predicted",index),
           ylab=paste("tm",index),
           xlim=limits,
           ylim=limits)
      for(i in 1:ncol(tmp)){abline(a = tmp[1,i], b=tmp[2,i], col=cols[i])}
      legend(x=limits[1]+25,y=limits[2]-25,
             colnames(tmp),
             lty=rep(1,ncol(tmp)),
             col=cols)
      
      dev.off()
    }
    if(model == "aggregate"){
      
      tmp = with(df,by(df, oli_img, function(x) rlm(refsamp ~ comppred, data = x)))
      tmp = sapply(tmp, coef)
      
      title = paste("aggregated model scatterplot and robust linear regression lines for",index)
      cols = seq(from=2,to=ncol(tmp)+1)
      
      png(outfile, width = 1100, height=800)
      plot(x=df$comppred,y=df$refsamp,
           main=title,
           xlab=paste("oli-predicted",index),
           ylab=paste("tm",index),
           xlim=limits,
           ylim=limits)
      
      m = rlm(refsamp ~ comppred, data=df)
      abline(coef = m$coefficients, lty=2, col="gray48", lwd=2.5)
      for(i in 1:ncol(tmp)){abline(a = tmp[1,i], b=tmp[2,i], col=cols[i])}
      legend(x=limits[1]+25,y=limits[2]-25,
             c("mean",colnames(tmp)),
             lty=c(2,rep(1,ncol(tmp))),
             col=c("gray48",cols))
      
      dev.off()
    }
  }
  
  plot_multi_cal_dif = function(df, index, model, outfile){
    if(index == "tcb"){limits = c(-500, 10000)}
    if(index == "tcg"){limits = c(-500, 5000)}
    if(index == "tcw"){limits = c(-6000, 1000)}
    if(index == "tca"){limits = c(-500, 5000)}
    meanrefsamp = mean(df$refsamp)
    if(model == "single"){
      newdf = aggregate(singlepreddif ~ oli_img , data = df, mean)
      newdf$diffm = meanrefsamp + newdf$singlepreddif
      title = paste("single model mean prediction differences from actual values for",index)
    }
    if(model == "aggregate"){
      newdf = aggregate(comppreddif ~ oli_img , data = df, mean)
      newdf$diffm = meanrefsamp + newdf$comppreddif
      title = paste("aggregate model mean prediction differences from actual values for",index)
    }
    
    png(outfile, width = 1100, height=800)
    d = density(df$refsamp)
    plot(d,
         xlim=limits,
         main=title,
         xlab=paste(index))
    abline(v=meanrefsamp,lty=2, col=1, lwd=2.5)
    abline(v=newdf$diffm, col=c(2:(length(newdf$diffm)+1)))
    legend(x=limits[2]-3000,y=max(d$y),
           c("mean",as.character(newdf$oli_img)),
           lty=c(2,rep(1,length(newdf$diffm))),
           col=c(1,2:(length(newdf$diffm)+1)))
    
    dev.off()
  }
  
  aggregate_cal_diag = function(sample_files, coef_files, index, outdir){
    
    if(class(sample_files) != "data.frame"){
      tbl = do.call("rbind", lapply(sample_files, read.csv, header = TRUE))
    } else {tbl = sample_files}
    #if(index == "tca"){model = rlm(refsamp ~ comppred, data=tbl)} else {
      model = rlm(refsamp ~ b2samp + b3samp + b4samp + b5samp + b6samp + b7samp, data=tbl)
      tbl$comppred = round(predict(model))
    #}
    tbl$singlepreddif = tbl$refsamp - tbl$singlepred
    tbl$comppreddif = tbl$refsamp - tbl$comppred
    tblcoef = model$coefficients
    r = cor(tbl$refsamp,tbl$comppred)
    coef = data.frame(index=index,yint=tblcoef[1],b2c=tblcoef[2],b3c=tblcoef[3],b4c=tblcoef[4],b5c=tblcoef[5],b6c=tblcoef[6],b7c=tblcoef[7],r=r)
    
    coeftbl = do.call("rbind", lapply(coef_files, read.csv, header = TRUE))
    outfile = file.path(outdir,paste(index,"_cal_combined_coef.csv",sep=""))
    write.csv(coeftbl, outfile, row.names=F)
    
    outfile = file.path(outdir,paste(index,"_cal_aggregate_sample.csv",sep=""))
    write.csv(tbl, outfile, row.names=F)
    outfile = file.path(outdir,paste(index,"_cal_aggregate_coef.csv",sep=""))
    write.csv(coef, outfile, row.names=F)
    
    outfile = file.path(outdir,paste(index,"_aggregate_mean_dif.png",sep=""))
    plot_multi_cal_dif(tbl,index,"aggregate",outfile)
    outfile = file.path(outdir,paste(index,"_single_mean_dif.png",sep=""))
    plot_multi_cal_dif(tbl,index,"single",outfile)
    
    outfile = file.path(outdir,paste(index,"_aggregate_regression.png",sep=""))
    plot_multi_cal(tbl,index,"aggregate",outfile)
    outfile = file.path(outdir,paste(index,"_single_regression.png",sep=""))
    plot_multi_cal(tbl,index,"single",outfile)
    
    return(tbl)
  }
  
  find_good_samples = function(coeffiles){
    len = length(coeffiles)
    id = substr(basename(coeffiles),1,16)
    r = array(NA,len)
    df = data.frame(id,r)
    for(i in 1:len){
      r = read.csv(coeffiles[i], header = TRUE)$r
      df$r[i] = r
    }
    n_goods = 0
    thresh = 0.8
    while(n_goods < 1){
      goods = which(df$r > thresh)
      n_goods = length(goods)
      thresh = thresh-0.05
    }
    return(as.character(df$id[goods]))
  }
  
  extract_good_samples = function(files, goodids){
    len = length(goodids)
    these = 0
    for(i in 1:len){
      match = grep(goodids[i], files)
      if(length(match) == 1){these = c(these,files[match])}
    }
    return(these[2:length(these)])
  }
  
  #######run the functions
  tcbsamps = list.files(dir,"tcb_cal_samp.csv",recursive=T,full.names=T)
  tcgsamps = list.files(dir,"tcg_cal_samp.csv",recursive=T,full.names=T)
  tcwsamps = list.files(dir,"tcw_cal_samp.csv",recursive=T,full.names=T)
  tcasamps = list.files(dir,"tca_cal_samp.csv",recursive=T,full.names=T)
  
  tcbcoef = list.files(dir,"tcb_cal_coef.csv",recursive=T,full.names=T)
  tcgcoef = list.files(dir,"tcg_cal_coef.csv",recursive=T,full.names=T)
  tcwcoef = list.files(dir,"tcw_cal_coef.csv",recursive=T,full.names=T)
  tcacoef = list.files(dir,"tca_cal_coef.csv",recursive=T,full.names=T)
  
  tcbgoods = find_good_samples(tcbcoef)
  tcggoods = find_good_samples(tcgcoef)
  tcwgoods = find_good_samples(tcwcoef)
  tcagoods = find_good_samples(tcacoef)
  
  tcbsamps = extract_good_samples(tcbsamps, tcbgoods)
  tcgsamps = extract_good_samples(tcgsamps, tcggoods)
  tcwsamps = extract_good_samples(tcwsamps, tcwgoods)
  tcasamps = extract_good_samples(tcasamps, tcagoods)
  
  tcbcoef = extract_good_samples(tcbcoef, tcbgoods)
  tcgcoef = extract_good_samples(tcgcoef, tcggoods)
  tcwcoef = extract_good_samples(tcwcoef, tcwgoods)
  tcacoef = extract_good_samples(tcacoef, tcagoods)
  
  btbl = aggregate_cal_diag(tcbsamps, tcbcoef, "tcb", outdir)
  gtbl = aggregate_cal_diag(tcgsamps, tcgcoef, "tcg", outdir)
  wtbl = aggregate_cal_diag(tcwsamps, tcwcoef, "tcw", outdir)
  atbl = aggregate_cal_diag(tcasamps, tcacoef, "tca", outdir)
  
  return(1)
}


