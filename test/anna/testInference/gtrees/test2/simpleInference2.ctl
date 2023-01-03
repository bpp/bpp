seed =  2
*seed =  -1

seqfile = simulate2.txt
*treefile=simulate_trees2.txt
Imapfile = simple2.Imap.txt
outfile = out2.txt
mcmcfile = mcmc2.txt
datefile = seqDates.txt

# fixed number of species/populations 
speciesdelimitation = 0

# fixed species tree

speciestree = 0
species&tree = 3  A B C
                  5  5 5
		  ((A, B), C);

*migration = 4
*A B 
*B A 
*A C 
*C B 
*
*migprior = .05 12

# phased data for population
phase =   0 0 0

# use sequence likelihood
usedata = 0

nloci = 10
*nloci = 1000
clock = 1
locusrate =0
*model = 0

# invgamma(a, b) for root tau & Dirichlet(a) for other tau's
tauprior = invgamma 3.0 0.08

thetaprior = gamma 8 2000 # gamma(a, b) for theta (estimate theta)

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  

# MCMC samples, locusrate, heredityscalars, Genetrees
print = 1 0 0 1   * 
burnin = 10000
*burnin = 8000
sampfreq = 4
nsample = 100000
*nsample = 1000
