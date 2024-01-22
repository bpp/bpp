seed = 1

seqfile = mammoth_nuclear.txt
Imapfile = map.txt
outfile = out.txt
mcmcfile = mcmc.txt
datefile = dates.txt

# fixed number of species/populations 
speciesdelimitation = 0

# fixed species tree
speciestree = 0

species&tree = 5 ASIAN2 MAM1 SAV1 FOR MAST
                  1 1 1 1 1 
		  (((ASIAN2,MAM1),(SAV1,FOR)),MAST);

# unphased data for all 4 populations
phase =   1  1  1  1 0

# use sequence likelihood
usedata = 1

# There are 347 loci in the alignment. Only 100 loci are 
# used as an example to decrease the run time. 
nloci = 100 #347

# do not remove sites with ambiguity data
cleandata = 0
locusrate = 3 5 10000000000 # The mutation rate prior is gamma(5, 10000000000)

thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
tauprior = gamma 35 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  

# MCMC samples, locusrate, heredityscalars, Genetrees
print = 1 0 0 0   * 
burnin = 40000
sampfreq = 2
nsample = 400000
