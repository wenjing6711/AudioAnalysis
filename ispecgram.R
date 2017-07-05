######################################
## ispecgram
######################################
ispecgram <- function(d, ftsize, sr, win, nov){
  
  nspec = nrow(d)
  ncol = ncol(d);
  
  nargin <- length(as.list(match.call())) -1
  
  if(nargin < 2){
    ftsize = 2*(nspec-1);
  }
  if(nargin < 3){}
  ## who cares?
  
  if(nargin < 4){
    win = ftsize;  ## doesn't matter either - assume it added up OK
  }
  if(nargin < 5){
    nov = ftsize/2;
  }
  
  hop = win - nov;
  
  if(nspec != (ftsize/2)+1){
    print('number of rows should be fftsize/2+1')
  }
  
  xlen = ftsize + (ncol-1) * hop;
  x = rep(0,xlen);
  
  halff = ftsize/2;   ## midpoint of win
  
  ## No reconstruction win (for now...)
  
  for(c in 1:ncol){
    ft = d[,c];
    ft = cbind(ft[1:(ftsize/2+1)], Conj(ft[c((ftsize/2):2)]));
  
    if(max(image(ifft(ft))) > 10^(-5)){
      print('imag oflow');
    }
  
    px = Re(ifft(ft));  ## no shift in specgram
  
    b = (c-1)*hop;
    x[b+c(1:ftsize)] = x[b+c(1:ftsize)] + px;
  }
  
  x = x*win/ftsize;  ## scale amplitude
  return(x)
} 
