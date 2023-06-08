          seed =  -1

       seqfile = frogs.txt
      Imapfile = frogs.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

# inference of species tree when the assignment and delimitation are given (Rannala and Yang, 2017)
  speciesdelimitation = 0        * fixed species delimitation
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
          speciestree = 1 0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  K  C  L  H
                    9  7 14  2
                 ((K, C), (L, H));
         phase =    1  1  1  1

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 5  * number of data sets in seqfile

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
