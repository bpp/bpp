seed =  1
*seed =  -1

seqfile = simulate_4.txt
treefile=simulate_trees.txt
Imapfile = simple.Imap.txt
datefile = dates_sim4.txt
*outfile = out.txt
*mcmcfile = mcmc.txt
seqDates = seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 5  A B C AB ABC
                  2  1 1 1  0
		  ((A #0.04, B #0.03)AB:.07 #.035, C)ABC:.08 #.05;

# phased data for population
phase =   0 0 0 0 0

# use sequence likelihood
*usedata = 1

loci&length = 500000 1000
*loci&length = 1 1000
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
