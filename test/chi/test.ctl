          seed = 12345

     traitfile = test.morph.txt
       seqfile = test.seq.txt
      Imapfile = test.Imap.txt
       jobname = out

# estimation of parameters under the multispecies coalescent model (Yang, 2002; Rannala and Yang, 2003)
* speciesdelimitation = 0  # fixed number of species/populations 
*         speciestree = 0  # fixed species tree

# joint species delimitation and species-tree inference (Yang and Rannala, 2014)
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
  speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
*         speciestree = 1 0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
    speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 3  K  C  H
                    3  4  2
                  ((K, C), H);

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 3  * number of data sets in seqfile

     cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 2000 # gamma(a, b) for theta (estimate theta)
      tauprior = gamma 2 1000 # gamma(a, b) for root tau & Dirichlet(a) for other tau's

      finetune = 1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 1000
      sampfreq = 2
       nsample = 10000
