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
  
  nfilts = nrow(aspectrum)
  nframes = ncol(aspectrum);
  
  if (fbtype == 'bark'){
    wts = tuneR:::fft2barkmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq)$wts;
  }else if(fbtype == 'mel'){
    wts = tuneR:::fft2melmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq)$wts;
  }else if (fbtype == 'htkmel'){
    wts = tuneR:::fft2melmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq, 1, 1)$wts
  }else if (fbtype == 'fcmel'){
    wts = tuneR:::fft2melmx(nfft, sr, nfilts, bwidth, minfreq, maxfreq, 1)$wts
  }else{
    print(paste('fbtype ', fbtype, 'not recognized', sep = " "))
  }
  
  ## Cut off 2nd half
  wts = wts[,1:(nfft/2+1)]
  
  ## Just transpose, fix up 
  ww = t(wts)%*%wts;
  iwts = t(wts)/max(mean(diag(ww))/100, sum(ww))
  #iwts = t(wts)/kronecker(matrix(1,1,nfilts),t(max(mean(diag(ww))/100, sum(ww))));
  ## Apply weights
  if (sumpower){
    spec = iwts%*%aspectrum
  }else{
    spec = (iwts%*%sqrt(aspectrum))^2;
  }
  
  result = list()
  result$spec = spec
  result$wts = wts
  result$iwts = iwts
  return(result)
}
