    seed =  -1

seqfile  = frogs.txt
Imapfile = frogs.Imap.txt
jobname  = out

# fixed number of species/populations 
speciesdelimitation = 0

# fixed species tree
speciestree = 0

species&tree = 4  K  C  L  H
                  9  7 14  2
                 (((K, C), L), H);

# unphased data for all 4 populations
phase =   1  1  1  1

# use sequence likelihood
usedata = 1

nloci = 5

# do not remove sites with ambiguity data
cleandata = 0

thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  

# MCMC samples, locusrate, heredityscalars, Genetrees
print = 1 0 0 0   * 
burnin = 8000
sampfreq = 2
nsample = 100000
