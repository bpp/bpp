 seed =  -1

 seqfile = yu2001.txt
 outfile = out.txt
 mcmcfile = mcmc.txt

 # fixed species delimitation and species tree
 speciesdelimitation = 0 
 speciestree = 0

 species&tree = 1  H
                 61  # max number of sequences
 # 0: no data (prior); 1:seq like
 usedata = 1    

 # number of data sets in seqfile
 nloci = 1    

 # remove sites with ambiguity data (1:yes, 0:no)?
 cleandata = 0    

 # gamma(a, b) for theta
 thetaprior = gamma 2 2000   
 
 # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
 finetune = 1: 2 0.00001 0.0001  0.0005 0.5 0.2 1.0  

 print = 1 0 0 0  * MCMC samples, locusrate, heredityscalars, Genetrees
 burnin = 4000
 sampfreq = 2
 nsample = 10000
