######################################
## invpowspec
######################################
invpowspec <- function(y, sr, wintime, steptime, excit){
  
    nrow = nrow(y) 
    ncol = ncol(y);
    
    nargin <- length(as.list(match.call())) -1
    
    if(nargin < 2){
      sr = 8000;
    }
    if(nargin < 3){
      wintime = 0.025;
    }
    if(nargin < 4){
      steptime = 0.010;
    }
    if(nargin < 5){
      r = c();
    }
    else{
      r = excit;
    }
    
    winpts = round(wintime*sr);
    steppts = round(steptime*sr);
    
    NFFT = 2^(ceiling(log(winpts)/log(2)));
    
    if(NFFT != 2*(nrow-1)){
      print('Inferred FFT size does not match specgram');
    }
    
    NOVERLAP = winpts - steppts;
    SAMPRATE = sr;
    
    ## Values coming out of rasta treat samples as integers, 
    ## not range -1..1, hence scale up here to match (approx)
    ##y = abs(specgram(x*32768,NFFT,SAMPRATE,WINDOW,NOVERLAP)).^2;
    
    xlen = winpts + steppts*(ncol - 1);
    
    if(length(r) == 0){
      r = rnorm(xlen);
    }
    r = r[1:xlen];
    ### ? S, F, T
    R = specgram(r/32768/12, NFFT, SAMPRATE, winpts, NOVERLAP)$S
    R = R%*%sqrt(y) ### dimension??
    x = ispecgram(R, NFFT, SAMPRATE, winpts, NOVERLAP)
    return(x)
}
