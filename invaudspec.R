######################################
## invaudspec
######################################
invaudspec <- function(aspectrum, sr, nfft, fbtype, minfreq, 
                       maxfreq, sumpower, bwidth){
  nargin <- length(as.list(match.call())) -1
  
  if(nargin < 2){
    sr = 16000
  }
  if(nargin < 3){
    nfft = 512
  }
  if(nargin < 4){ 
    fbtype = 'bark'
  }
  if(nargin < 5){ 
    minfreq = 0
  }
  if(nargin < 6){  

    maxfreq = sr/2
  }
  if(nargin < 7){
    sumpower = TRUE
  }
  if(nargin < 8){
    bwidth = 1.0
  }
  
  nfilts = ncol(aspectrum)
  nframes = nrow(aspectrum);
  
  if (fbtype == 'bark'){
    wts = tuneR:::fft2barkmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq)$wts;
  }else if(fbtype == 'mel'){
    wts = tuneR:::fft2melmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq)$wts;
  }else if (fbtype == 'htkmel'){
    wts = tuneR:::fft2melmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq, 1, 1)$wts
  }else if (fbtype == 'fcmel'){
    wts = tuneR:::fft2melmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq, 1)$wts
  }else{
    print(paste('fbtype', fbtype, 'not recognized', sep = " "))
  }
  
  ## Cut off 2nd half
  wts = wts[,1:(nfft/2+1)]
  
  ## Just transpose, fix up 
  ww = t(wts)%*%wts;

  max_vec = colSums(ww)
  max_vec[which(mean(diag(ww))/100 > max_vec)] = mean(diag(ww))/100

  iwts = t(wts)/matrix(rep(t(max_vec),nfilts), nrow = length(max_vec), ncol = nfilts, byrow = FALSE);
  ## Apply weights
  if (sumpower){
    spec = t(iwts%*%t(aspectrum))
  }else{
    spec = t((iwts%*%sqrt(t(aspectrum)))^2);
  }
  
  result = list()
  result$spec = spec
  result$wts = wts
  result$iwts = iwts
  return(result)
}
