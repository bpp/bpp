seed =  1
*seed =  -1

seqfile = simulate_frogs.txt
treefile=simulate_trees.txt
Imapfile = frogs.Imap.txt
datefile = dates_sim2phase.txt
*outfile = out.txt
*mcmcfile = mcmc.txt
seqDates = seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 2  KK C
                  1  1
		  (KK #0.04, C #0.03):.07 #.035;

# phased data for population
*phase =   0 0
phase =   1 1

# use sequence likelihood
*usedata = 1

loci&length = 100000 1000
*loci&length = 10 1000
clock = 1
locusrate =0
model = 0

# invgamma(a, b) for root tau & Dirichlet(a) for other tau's
*tauprior = invgamma 3 0.002

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
*finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  

# MCMC samples, locusrate, heredityscalars, Genetrees
*print = 1 0 0 0   * 
*burnin = 8000
*sampfreq = 2
*nsample = 100000
