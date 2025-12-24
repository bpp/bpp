          seed = -1

       seqfile = loci_realign.txt
      Imapfile = Imap.txt
       jobname = out

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
*         speciestree = 1 0.1 0.1 0.2    * species tree SPR/SNL

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 6  G  C  R  L  A  Q
                    2  2  2  2  2  2
                 (R, ((C, G) b, ((A, Q) d, L) c ) a) o;

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 100 * number of data sets in seqfile

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

*    thetaprior = 3 0.04 e  # Inv-gamma(a, b) for theta (integrated out by default; add E to also sample theta)
*      tauprior = 3 0.2     # Inv-gamma(a, b) for root tau
    thetaprior = gamma 2 100  # gamma(a, b) for theta
      tauprior = gamma 2 10   # gamma(a, b) for root tau
 
       wprior = 20 1
      migration = 2 
	              A b
                      R Q

# finetune: auto-adjust step lengths during burnin
      finetune = 1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 32000
      sampfreq = 2
       nsample = 200000
       threads = 2 1 1
