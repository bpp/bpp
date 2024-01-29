seed =  1
*seed =  -1

seqfile = simulate_IM.txt
treefile=simulate_trees.txt
Imapfile = simple.Imap.txt
datefile = datesIM.txt
*outfile = out.txt
*mcmcfile = mcmc.txt
seqDates = seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 2  A B
                  2  1
		  (A #0.04, B #0.04):.9 #.035;

migration = 2
A B 0.3
B A 0.3

# phased data for population
phase =   0 0
*phase =   1 1

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
