#' draw a random sample of image pixels
#'
#' draw a histogram equalized random sample of image pixels
#' @param img matrix image matrix to be sampled 
#' @param bins integer number of histogram bins from which to sample
#' @param n integer number of pixels to sample per histogram bin


sample_it = function(img, bins, n){
  
  mi = min(img, na.rm=T)
  ma = max(img, na.rm=T)
  
  step = (ma - mi)/bins
  breaks = seq(mi,ma,step)
  
  min_samp = array(n, bins)
  for(i in 1:(length(breaks)-1)){
    these = which(img > breaks[i] & img <= breaks[i+1])
    if(i == 1){samp = sample(these, size=min(min_samp[i],length(these)))} else {
      samp = c(samp, sample(these, size=min(min_samp[i],length(these))))
    } 
  }
  return(samp)
}