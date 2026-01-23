seed = 12345

seqfile = sim_veryhighmig_seqs.txt
Imapfile = sim.Imap.txt
jobname = delimit_mig

# Species delimitation with fixed guide tree (METHOD_10)
# speciesdelimitation = 1: enable delimitation
# speciestree = 0: fixed species tree topology
speciesdelimitation = 1 0 2
speciestree = 0

# Two populations A and B - guide tree
species&tree = 2  A B
                  10 10
                  (A, B);

# Enable migration estimation
migration = 2
A B
B A

# Prior on migration rate M ~ Gamma(alpha, beta)
# Mean = alpha/beta, so Gamma(2, 0.2) has mean 10
wprior = 2 0.2

# Phased data
phase = 0 0

# Use sequence data
usedata = 1
nloci = 50

cleandata = 0

# Priors
thetaprior = gamma 2 100    # Gamma(2, 100) -> mean = 0.02
tauprior = gamma 2 200      # Gamma(2, 200) -> mean = 0.01

# Finetune parameters (auto-adjust during burnin)
finetune = 1

# MCMC settings
print = 1 0 0 0
burnin = 4000
sampfreq = 2
nsample = 10000
