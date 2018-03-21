          seed = 1234567

       seqfile = testbed/ziheng/4/data/test4s.txt
      Imapfile = testbed/ziheng/4/data/Imap.4s.txt
       outfile = testbed/ziheng/4/out/out.txt
      mcmcfile = testbed/ziheng/4/out/mcmc.txt


*       diploid = 1        * 0: phased sequences, 1: unphased diploid sequences
*    breakpoint = 0        * 0: nothing;  1 : save;  2: read

* speciesdelimitation = 0                * fixed species tree
* speciesdelimitation = 1 0 2            * speciesdelimitation algorithm0 And finetune(e)
* speciesdelimitation = 1 1 2 1          * speciesdelimitation algorithm1 finetune (a m)
*         speciestree = 1  0.4 0.1 0.1   * speciestree pSlider ExpandRatio ShrinkRatio
          speciestree = 1 0               * species tree fixed

speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 4  A  B  C  D
                    2  2  2  2
                 (((A, B), C), D);

       usedata = 1 * 0: no data(prior); 1:seq Like
         nloci = 2 * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.002 e  # invgamma(a, b) for theta
      tauprior = 3 0.02    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's
*    thetaprior = 2 1000  # gamma(a, b) for theta
*      tauprior = 2 2000  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet
*      heredity = 1 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma
* sequenceerror = 0 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

       finetune = 1: 2 0.01 0.1 0.02 0.8 0 0 # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 10000
