#' Create an aggregate MSS TC model
#'
#' Create an aggregate MSS TC model with diagnostic figures
#' @param dirname character. the directory to the calibration folder
#' @import ggplot2
#' @export


cal_mss_tc_aggregate_model = function(dir){
  plot_multi_cal = function(df,index, model, outfile){
    
    if(index == "tcb"){limits = c(-500, 10000)}
    if(index == "tcg"){limits = c(-500, 5000)}
    if(index == "tcw"){limits = c(-6000, 1000)}
    if(index == "tca"){limits = c(-500, 5000)}
    if(model == "single"){
      
      
      tmp = with(df,by(df, mss_img, function(x) rlm(refsamp ~ singlepred, data = x)))
      tmp = sapply(tmp, coef)
      
      title = paste("single model scatterplot and robust linear regression lines for",index)
      png(outfile, width = 1100, height=800)
      cols = seq(from=2,to=ncol(tmp)+1)
      plot(x=df$singlepred,y=df$refsamp,
           main=title,
           xlab=paste("mss-predicted",index),
           ylab=paste("tm",index),
           xlim=limits,
           ylim=limits)
      for(i in 1:ncol(tmp)){abline(a = tmp[1,i], b=tmp[2,i], col=cols[i])}
      legend(x=limits[1]+25,y=limits[2]-25,
             colnames(tmp),
             lty=rep(1,ncol(tmp)),
             col=cols)
      
      
      dev.off()
      
      
      #       g = ggplot(df, aes(x=singlepred, y=refsamp, colour=mss_img)) +
      #         geom_density2d(bins=3)+
      #         geom_smooth(method="rlm", se = FALSE)+
      #         xlim(limits)+
      #         ylim(limits)+
      #         xlab(paste("mss predicted",index))+
      #         ylab(paste("tm",index))+
      #         ggtitle(paste("single model 2d data density contours and robust linear regression lines for",index))+
      #         coord_fixed(ratio = 1)+
      #         theme_bw()
    }
    if(model == "aggregate"){
      
      tmp = with(df,by(df, mss_img, function(x) rlm(refsamp ~ comppred, data = x)))
      tmp = sapply(tmp, coef)
      
      title = paste("aggregated model scatterplot and robust linear regression lines for",index)
      cols = seq(from=2,to=ncol(tmp)+1)
      
      png(outfile, width = 1100, height=800)
      plot(x=df$comppred,y=df$refsamp,
           main=title,
           xlab=paste("mss-predicted",index),
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
      
      
      #print(g)
      dev.off()
      
      #       g=ggplot()+
      #         geom_density2d(bins=3,data=df, aes(x=comppred, y=refsamp, colour=mss_img))+
      #         geom_smooth(method="rlm", data=df, aes(x=comppred, y=refsamp, colour=mss_img), se = FALSE)+
      #         geom_smooth(method="rlm", data=df, aes(x=comppred, y=refsamp), colour="black", size=1, se = FALSE, linetype="dashed")+
      #         xlim(limits)+
      #         ylim(limits)+
      #         xlab(paste("mss predicted",index))+
      #         ylab(paste("tm",index))+
      #         ggtitle(paste("aggregate model 2d data density contours and robust linear regression lines for",index))+
      #         coord_fixed(ratio = 1)+
      #         theme_bw()
    }
    
    #     png(outfile, width = 1100, height=800)
    #     print(g)
    #     dev.off()
    #return(g)
  }
  
  plot_multi_cal_dif = function(df, index, model, outfile){
    if(index == "tcb"){limits = c(-500, 10000)}
    if(index == "tcg"){limits = c(-500, 5000)}
    if(index == "tcw"){limits = c(-6000, 1000)}
    if(index == "tca"){limits = c(-500, 5000)}
    meanrefsamp = mean(df$refsamp)
    if(model == "single"){
      newdf = aggregate(singlepreddif ~ mss_img , data = df, mean)
      newdf$diffm = meanrefsamp + newdf$singlepreddif
      title = paste("single model mean prediction differences from actual values for",index)
    }
    if(model == "aggregate"){
      newdf = aggregate(comppreddif ~ mss_img , data = df, mean)
      newdf$diffm = meanrefsamp + newdf$comppreddif
      title = paste("aggregate model mean prediction differences from actual values for",index)
    }
    
    #     g = ggplot(df, aes(x=refsamp))+
    #       geom_density(colour="black",linetype="dashed")+
    #       geom_vline(data = newdf,aes(xintercept = diffm, colour = mss_img), show_guide = TRUE, size=0.5) +
    #       geom_vline(xintercept = meanrefsamp , colour = "black", size=1, linetype="dashed") +
    #       xlim(limits)+
    #       xlab(paste(index))+
    #       ggtitle(title)+
    #       theme_bw()
    
    png(outfile, width = 1100, height=800)
    d = density(df$refsamp)
    plot(d,
         xlim=limits,
         main=title,
         xlab=paste(index))
    abline(v=meanrefsamp,lty=2, col=1,lwd=2.5)
    abline(v=newdf$diffm, col=c(2:(length(newdf$diffm)+1)))
    legend(x=limits[2]-3000,y=max(d$y),
           c("mean",as.character(newdf$mss_img)),
           lty=c(2,rep(1,length(newdf$diffm))),
           col=c(1,2:(length(newdf$diffm)+1)))
    
    #print(g)
    dev.off()
    #return(g)
  }
  
  aggregate_cal_diag = function(sample_files, coef_files, index, outdir){
    
    if(class(sample_files) != "data.frame"){
      tbl = do.call("rbind", lapply(sample_files, read.csv, header = TRUE))
    } else {tbl = sample_files}
    if(index == "tca"){model = rlm(refsamp ~ comppred, data=tbl)} else {
      model = rlm(refsamp ~ b1samp + b2samp + b3samp + b4samp, data=tbl)
      tbl$comppred = round(predict(model))
    }
    tbl$singlepreddif = tbl$refsamp - tbl$singlepred
    tbl$comppreddif = tbl$refsamp - tbl$comppred
    tblcoef = model$coefficients
    r = cor(tbl$refsamp,tbl$comppred)
    if(index == "tca"){coef = data.frame(index=index,yint=tblcoef[1],b1c=tblcoef[2],r=r)} else {
      coef = data.frame(index=index,yint=tblcoef[1],b1c=tblcoef[2],b2c=tblcoef[3],b3c=tblcoef[4],b4c=tblcoef[5],r=r)}
    
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
  
  outdir = file.path(dir,"aggregate_model")
  dir.create(outdir, recursive=T, showWarnings=F)
  
  tcbsamps = list.files(dir,"tcb_cal_samp.csv",recursive=T,full.names=T)
  tcgsamps = list.files(dir,"tcg_cal_samp.csv",recursive=T,full.names=T)
  tcwsamps = list.files(dir,"tcw_cal_samp.csv",recursive=T,full.names=T)
  tcasamps = list.files(dir,"tca_cal_samp.csv",recursive=T,full.names=T)
  
  tcbcoef = list.files(dir,"tcb_cal_coef.csv",recursive=T,full.names=T)
  tcgcoef = list.files(dir,"tcg_cal_coef.csv",recursive=T,full.names=T)
  tcwcoef = list.files(dir,"tcw_cal_coef.csv",recursive=T,full.names=T)
  tcacoef = list.files(dir,"tca_cal_coef.csv",recursive=T,full.names=T)
  
  btbl = aggregate_cal_diag(tcbsamps, tcbcoef, "tcb", outdir)
  gtbl = aggregate_cal_diag(tcgsamps, tcgcoef, "tcg", outdir)
  wtbl = aggregate_cal_diag(tcwsamps, tcwcoef, "tcw", outdir)
  
  tcasamps = do.call("rbind", lapply(tcasamps, read.csv, header = TRUE))
  comppred = atan(gtbl$comppred/btbl$comppred) * (180/pi) * 100
  tcasamps = data.frame(tcasamps,comppred)
  atbl = aggregate_cal_diag(tcasamps, tcacoef, "tca", outdir)
  
  #create plane plots
  tcb_samp_file = file.path(outdir,"tcb_cal_aggregate_sample.csv")
  tcg_samp_file = file.path(outdir,"tcg_cal_aggregate_sample.csv")
  tcw_samp_file = file.path(outdir,"tcw_cal_aggregate_sample.csv")
  outfile = file.path(outdir,"aggregate_cal_tc_planes_comparison.png") 
  #make_tc_planes_comparison(tcb_samp_file, tcg_samp_file, tcw_samp_file, outfile)
}


