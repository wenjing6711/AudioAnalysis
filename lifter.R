lifter <- function(x, lift, invs){

	nargin <- length(as.list(match.call())) -1

	if(nargin < 2){
		lift = 0.6   % liftering exponent

	}
	if(nargin < 3){  
		invs = 0      % flag to undo liftering
	}
	
	ncep = ncol(x)
	nfrm = nrow(x)

	if (lift == 0){
  		y = x;
	}
	else{

  		if (lift > 0){
    			if (lift > 10){
      			print(paste('Unlikely lift exponent of', toString(lift),'(did you mean -ve?)', sep = " "));
    			}	
    			liftwts = c(1, (c(1:(ncep-1))^lift));
  		}
		else if (lift < 0){
    			% Hack to support HTK liftering
    			L = -lift;
    			if (L != round(L)){ 
      			print(paste('HTK liftering value', toString(L),'must be integer', sep = " "))
    			}
   	 		liftwts = c(1, c(1+L/2*sin(c(1:(ncep-1)*pi/L)));
  		}
	}

  	if (invs){
    		liftwts = 1/liftwts;
 	}

  	y = diag(liftwts)%*%x;
	return(y)
}

