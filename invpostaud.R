######################################
## invpostaud
######################################
invpostaud <- function(y,fmax,fbtype,broaden){

  nargin <- length(as.list(match.call())) -1
  
  if(nargin < 3){
    fbtype = 'bark';
  }
  if(nargin < 4){
  ## By default, don't add extra flanking bands
    broaden = 0;
  }
  
  nbands = ncol(y)
  nframes = nrow(y);
 
  if(fbtype == 'bark'){
    bandcfhz = tuneR:::bark2hz(seq(0, hz2bark(fmax), nbands));
  }
  else if(fbtype == 'mel'){
    bandcfhz = tuneR:::mel2hz(seq(0, hz2mel(fmax), nbands));
  }
  else if(fbtype == 'htkmel' || fbtype == 'fcmel'){
    bandcfhz = tuneR:::mel2hz(seq(0, hz2mel(fmax,1), nbands),1);
  }
  else{
    print(paste('unknown fbtype', fbtype, sep = " "))
  }
  
  ## Remove extremal bands (the ones that got duplicated)
  bandcfhz = bandcfhz[(1+broaden):(nbands-broaden)];
  
  ## Hynek's magic equal-loudness-curve formula
  fsq = bandcfhz^2;
  ftmp = fsq + 1.6*(10^5);
  eql = ((fsq/ftmp)^2)%*%((fsq + 1.44*(10^6))/(fsq + 9.61*(10^6)));


  ## cube expand
  x = y^(1/.33);

  #### squash the zero in the eql curve
  if(eql[1] == 0){  ## or maybe always
    eql[1] = eql[2];
    eql[length(eql)] = eql[length(eql)-1];
  }

  ## weight the critical bands
  x = x[(1+broaden):(nbands-broaden),]/matrix(rep(t(eql),nframes), nrow = length(eql), ncol = nframes, byrow = FALSE);
  result = list()
  result$x = x
  result$eql = eql
  return(result)                                            
}

