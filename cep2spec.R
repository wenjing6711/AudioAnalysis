cep2spec <- function (cep, nfreq, type){
  ## spec = cep2spec(cep, nfreq, type)
  ##   Reverse the cepstrum to recover a spectrum.
  ##   i.e. converse of spec2cep
  ##   nfreq is how many points to reconstruct in spec
  nargin <- length(as.list(match.call())) -1
  
  if(nargin < 2){
    nfreq = 21 
  }
  if(nargin < 3){   
    type = 2 ## type of DCT
  }
  
  ncep = ncol(cep)
  ncol = nrow(cep)
  
  ## Make the DCT matrix
  dctm = matrix(0, nrow = nfreq, ncol = ncep);
  idctm = matrix(0, nrow = ncep, ncol = nfreq);

  if( type == "t2" || type == "t3"){
  ## this is the orthogonal one, so inv matrix is same as fwd matrix
    for(i in 1:ncep){
      dctm[,i] = cos((i-1)*seq(1,2*nfreq-1,2)/(2*nfreq)*pi)*sqrt(2/nfreq);
    }
	
    if(type == "t2"){ 
  	## make it unitary! (but not for HTK type 3)
      dctm[,1] = dctm[,1]/sqrt(2);
    }
    else{
      dctm[,1] = dctm[,1]/2;    
    }
	
    idctm = t(dctm);
  }
  else if(type == "t4"){ 
  ## type 1 with implicit repetition of first, last bins
  ## so all we do is reconstruct the middle nfreq rows of an nfreq+2 row idctm
    
	for(i in 1:ncep){
    	## 2x to compensate for fact that only getting +ve freq half
      	idctm[i,] = 2*cos((i-1)*t(c(1:nfreq))/(nfreq+1)*pi);
    	}
	## fixup 'non-repeated' basis fns 
    	idctm[c(1, ncep),] = idctm[c(1, ncep),]/2;
  }
  else{ 
	## dpwe type 1 - idft of cosine terms
    	for(i in 1:ncep){
  	## 2x to compensate for fact that only getting +ve freq half
    		idctm[i,] = 2*cos((i-1)*t(c(0:(nfreq-1)))/(nfreq-1)*pi);
    	}
  	## fixup non-repeated basis fns 
    	idctm[c(1, ncep),] = 0.5* idctm[c(1, ncep),];
  }  

  spec = exp(cep%*%idctm);
  
  result = list()
  result$spec = spec
  result$idctm = idctm
  return(result)
}
