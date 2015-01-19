#' Run the neatline program
#'
#' Run the neatline program   
#' @param dir character. full path to a composite directory
#' @param cores numeric. the number of cores to be used for parellel processing
#' @import foreach
#' @import doParallel


run_neatline = function(dir, cores){
 
  dir = "L:/composites/test9"
  files = list.files(dir, "composite_stack.bsq$", recursive=T, full.name=T)
  tcafile = files[grep("tca_composite_stack.bsq$", files)]
  tcbfile = files[grep("tcb_composite_stack.bsq$", files)]
  tcgfile = files[grep("tcg_composite_stack.bsq$", files)]
  tcwfile = files[grep("tcw_composite_stack.bsq$", files)]
  
  yearseries = sort(as.numeric(unique(substr(basename(list.files(dir, "composite.bsq$", recursive=T)),1,4))))
  tcafile = list.files(dir, "tca_composite_stack.bsq$", recursive=T, full.name=T)
  tca = brick(tcafile)
  tcb = brick(tcbfile)
  tcg = brick(tcgfile)
  tcw = brick(tcwfile)
bases = basename(c(tcafile,tcbfile,tcgfile,tcwfile))
bname = sub("composite", "fitted", bases)

#chunck up the data
n_rows = floor(nrow(tca)/25)
chunks1 = seq(from=1,to=nrow(tca),n_rows)
cdf = data.frame(chunk=seq(1:length(chunks1)),from=NA,to=NA,n_col=ncol(tca))
for(i in 1:nrow(cdf)){
  cdf$from[i] = chunks1[i]
  cdf$to[i] = chunks1[i+1]-1
}
cdf$to[nrow(cdf)] = nrow(tca)
outdir = file.path(dir,"fitted")
dir.create(outdir)

#pass chunks off to be processed by neatline
cores=3
cl=makeCluster(cores)
registerDoParallel(cl)
o = foreach(i=1:nrow(cdf), .combine="c",.packages=c("raster","ecp", "sp", "rgdal", "zoo")) %dopar% neatline(cdf[i,],tca,tcb,tcg,tcw,outdir,bname)
stopCluster(cl)

}