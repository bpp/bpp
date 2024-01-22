seed =  1
*seed =  -1

seqfile = IM1/simulate_IM.txt
treefile=IM1/simulate_trees.txt
Imapfile = IM1/simple.Imap.txt
datefile = IM1/datesIM2.txt
seqDates = IM1/seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 2  A B
                  1  1
		  (A #0.004, B #0.004):.04 #.0035;

migration = 2
A B 0.2
B A 0.2

# phased data for population
phase =   0 0
*phase =   1 1

# use sequence likelihood
*usedata = 1

loci&length = 100000 100
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
