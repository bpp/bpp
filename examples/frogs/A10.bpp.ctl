          seed =  -1

       seqfile = frogs.txt
      Imapfile = frogs.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

# species delimitation using a fixed guide tree (Yang and Rannala, 2010; Rannala and Yang, 2013)
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
  speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
          speciestree = 0        * fixed species guide tree

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  K  C  L  H
                    9  7 14  2
                 ((K, C), (L, H));
         phase =    1  1  1  1

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 1  * number of data sets in seqfile

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
      tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

# finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr
      finetune = 1: 5 0.001 0.001  0.001 0.3 0.33 1.0

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 100000
