          seed = 123

       seqfile = testbed/ziheng/1/data/test3s.txt
*       seqfile =my3s.het.txt
      Imapfile = testbed/ziheng/1/data/Imap.3s.txt
       outfile = testbed/ziheng/1/out/out.txt
      mcmcfile = testbed/ziheng/1/out/mcmc.txt

*        diploid = 1        * 0: phased sequences, 1: unphased diploid sequences
*     breakpoint = 2        * 0: nothing;  1 : save;  2: read

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2  * speciesdelimitation algorithm0 and finetune(e)
 speciesdelimitation = 1 1 2 1 * speciesdelimitation algorithm1 finetune (a m)
*         speciestree = 0 * species tree fixed
         speciestree = 1  0 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted
  species&tree = 3  A  B  C
                    2  1  1
                  ((A, B), C);
       diploid =    1  0  0         * 0: phased sequences, 1: unphased diploid sequences

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 1    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.02 e   # invgamma(a, b) for theta
      tauprior = 3 0.04    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

*     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
*      heredity = 1 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
* sequenceerror = 0 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

       finetune = 1: 0.2 0.002 0.01 0.001 0.1 0.1 0.1  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 100000
