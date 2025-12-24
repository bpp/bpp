          seed =  -1

       seqfile = Rokas2003-5species.phy
      Imapfile = Rokas2003-5species-Imap.txt
       jobname = out

*   speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * speciesdelimitation algorithm0 and finetune(e)
*  speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)

*  uniformrootedtrees = 1         * 0 means uniform labeled histories
 * speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees; 2:user probs

  species&tree = 5 Scer Spar Smik Skud Sbay
                   1    1    1    1    1
                   ((((Scer,Spar) a, Smik) b, Skud) c,Sbay) r;

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 106  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 100  # gamma(a, b) for theta
      tauprior = gamma 2 20   # gamma(a, b) for root tau & Dirichlet(a) for other tau's
        wprior = 20 2
     migration = 1
	             Skud Sbay

*      heredity = 0 heredity.txt   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
*     locusrate = 1 5.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)

# finetune: auto-adjust step lengths during burnin
       finetune = 1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 10000
      sampfreq = 2
       nsample = 100000
       threads = 2 1 1
