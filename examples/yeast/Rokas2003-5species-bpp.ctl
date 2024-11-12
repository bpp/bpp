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
                   ((((Scer,Spar),Smik),(Skud,(Sbay)H[&phi=0.5,&tau-parent=no])),H[&tau-parent=yes])R;


       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 106  * 29    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.04 e   # gamma(a, b) for theta
      tauprior = 3 0.2  # gamma(a, b) for root tau & Dirichlet(a) for other tau's
       phiprior = 1 1

*      heredity = 0 heredity.txt   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
*     locusrate = 1 5.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)

       finetune = 1:  2.60867 0.00833 0.00327 0.00084 0.01261 0.11280 .01 .01 # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars Genetrees
        burnin = 10000
      sampfreq = 2
       nsample = 100000
#       scaling = 1
#       threads = 4
