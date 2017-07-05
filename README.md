# AudioAnalysis
R code for MFCC inversion function

- invmelfcc
Attempt to invert plp cepstra back to a full spectrum
and even a waveform.  Takes all the same options as melfcc.

- ispecgram 
Overlap-add the inverse of the output of specgram

- lifter 
Apply lifter to matrix of cepstra (one per column)

- cep2spec 
Reverse the cepstrum to recover a spectrum. i.e. converse of spec2cep

- invpostaud 
invert the effects of postaud (loudness equalization and cube root compression)

- invaudspec
Invert (as best we can) the effects of audspec()

- invpowspec
Attempt to go back from specgram-like powerspec to audio waveform by scaling specgram of white noise
