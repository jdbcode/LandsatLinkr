#' Create TC plane comparison plots
#'
#' Create TC plane comparison plots between modeled MSS TC and corresponding TM TC sample
#' @param tcb_samp_file character. full path name to TC brightness sample
#' @param tcg_samp_file character. full path name to TC greenness sample
#' @param tcw_samp_file character. full path name to TC wetness sample
#' @param model character. which model to use options: "single" or "composite"
#' @param outfile. character. define the full path of the .png file to be created
#' @import gridExtra
#' @import ggplot2
#' @export

make_tc_planes_comparison = function(tcbsamp, tcgsamp, tcwsamp, outfile){
  
  if(class(tcbsamp) == "character"){tcbsamp = read.csv(tcbsamp)}
  if(class(tcgsamp) == "character"){tcgsamp = read.csv(tcgsamp)}
  if(class(tcwsamp) == "character"){tcwsamp = read.csv(tcwsamp)}
  
  mssbgdf = data.frame(x=tcbsamp$singlepred,y=tcgsamp$singlepred)
  mssbwdf = data.frame(x=tcbsamp$singlepred,y=tcwsamp$singlepred)
  msswgdf = data.frame(x=tcwsamp$singlepred,y=tcgsamp$singlepred)
  
  tmbgdf = data.frame(x=tcbsamp$refsamp,y=tcgsamp$refsamp)
  tmbwdf = data.frame(x=tcbsamp$refsamp,y=tcwsamp$refsamp)
  tmwgdf = data.frame(x=tcwsamp$refsamp,y=tcgsamp$refsamp)
  
  mssbg = plot_tc_planes(mssbgdf,"bg","mss")
  mssbw = plot_tc_planes(mssbwdf,"bw","mss")
  msswg = plot_tc_planes(msswgdf,"wg","mss")
  
  tmbg = plot_tc_planes(tmbgdf,"bg","tm")
  tmbw = plot_tc_planes(tmbwdf,"bw","tm")
  tmwg = plot_tc_planes(tmwgdf,"wg","tm")
  
  plot_list = c(list(mssbg),list(mssbw),list(msswg),list(tmbg),list(tmbw),list(tmwg))
  png(outfile, width = 1100, height=700)
  do.call(grid.arrange, c(plot_list, list(ncol=3)))
  dev.off()
}