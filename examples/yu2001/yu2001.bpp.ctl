 seed =  -1

 seqfile = yu2001.txt
 jobname = out

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
 
 # auto-tune step lengths during burnin (append key:val pairs such as
 #   Gage:5 Gspr:0.001 mix:0.3   to override defaults)
 finetune = 1

 print = 1 0 0 0  * MCMC samples, locusrate, heredityscalars, Genetrees
 #burnin = 4000
 burnin = 20
 #sampfreq = 2
 #nsample = 10000
 nsample = 10
