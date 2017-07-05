######################################
## invmelfcc
######################################
invmelfcc <- function(cep, sr, wintime = 0.025, hoptime = 0.010,
                      numcep = 13, lifterexp = 0.6, sumpower = TRUE, 
                      preemph = 0.97, dither = FALSE, minfreq = 0, 
                      maxfreq = 4000, nbands = 40, bwidth = 1.0, dcttype = 2,
                      fbtype = 'mel', usecmp = FALSE, modelorder = NULL, 
                      broaden = 0, excitation = c()){


  winpts = ceiling(wintime*sr);
  nfft = 2^(ceiling(log(winpts)/log(2)));

  cep = lifter(cep, lifterexp, inv = TRUE);   ## 3rd arg nonzero means undo liftering

## Need to reconstruct the two extra flanking bands for invpostaud to delete
## (if we're doing usecmp)
   pspc = cep2spec(cep, nbands+2*broaden, dcttype)$spec;
   
   if (usecmp){
    aspc = invpostaud(pspc, maxfreq, fbtype, broaden)$x;
   }
   else{
    aspc = pspc;
   }
   
   ## Undo the auditory spectrum
   spec = invaudspec(aspc, sr, nfft, fbtype, minfreq, maxfreq, sumpower, bwidth)$spec;
   
   ## Back to waveform (modulate white noise, or specified excitation)
   x = invpowspec(spec, sr, wintime, hoptime, excitation);
   
   if(preemph != 0){
   ## Undo the original preemphasis
    x = rtf(c(1,-preemph),1, x) 
   }
   
   result = list()
   result$x = x
   result$aspc = aspc
   result$spec = spec
   return(result)
}


