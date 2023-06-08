          seed =  -1

       seqfile = frogs.txt
      Imapfile = frogs.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

# estimation of parameters under the multispecies coalescent model (Yang, 2002; Rannala and Yang, 2003)
  speciesdelimitation = 0  # fixed number of species/populations 
          speciestree = 0  # fixed species tree

  species&tree = 4  K  C  L  H
                    9  7 14  2
                 (((K, C), L), H);
         phase =    1  1  1  1  # unphased data for all 4 populations

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 5  * number of data sets in seqfile

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
      tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
      finetune = 1: 5 0.001 0.001  0.001 0.3 0.33 1.0

         print = 1 0 0 0   # MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 100000
