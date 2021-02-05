          seed =  -1

       seqfile = yu2001.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

 speciesdelimitation = 0 * fixed species tree

  species&tree = 1  H
                 61  * max number of sequences

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 1    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 2000   # gamma(a, b) for theta
*      tauprior = 3 0.004   # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

       finetune = 1: 2 0.00001 0.0001  0.0005 0.5 0.2 1.0  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0  * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 4000
      sampfreq = 2
       nsample = 10000
