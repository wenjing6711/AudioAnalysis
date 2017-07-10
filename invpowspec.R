######################################
## invpowspec
######################################
invpowspec <- function(y, sr, wintime, steptime, excit){
      
    nrow = ncol(y) 
    ncol = nrow(y);
    
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
      r = c()
    }
    else{
      r = excit;
    }
    
    winpts = ceiling(wintime*sr);
    steppts = round(steptime*sr);
    
    NFFT = 2^(ceiling(log(winpts)/log(2)));
    
    if(NFFT != 2*(nrow-1)){
      print('Inferred FFT size does not match specgram');
    }
    
    NOVERLAP = winpts - steppts;
    SAMPRATE = sr;
    
       
    xlen = winpts + steppts*(ncol - 1);
    
    if(length(r) == 0){
      r = rnorm(xlen);
    }
    r = r[1:xlen];
  
    R = specgram(r, NFFT+1, SAMPRATE, winpts, NOVERLAP+0.1)$S
    R = R%*%sqrt(y) 
    x = ispecgram(R, NFFT, SAMPRATE, winpts, NOVERLAP)
    return(x)
}
