          seed =  -1

       seqfile = loci_realign.txt
      Imapfile = Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

  speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2    * species delimitation rjMCMC algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1 * species delimitation rjMCMC algorithm1 finetune (a m)
*          speciestree = 1        * species tree NNI/SPR

*   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 6  G  C  R  L  A  Q
                    2  2  2  2  2  2
        ((R, (Q)h[&phi=0.3,&tau-parent=no]) g, (f[&tau-parent=yes,&phi=0.3], (((((G, C)b)f[&tau-parent=no], A)e, h[&tau-parent=yes])d, L)c)a)o;
#                 (R, ((C, G), ((A, Q), L)));

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 100 * number of data sets in seqfile

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.04 e  # Inv-gamma(a, b) for theta (integrated out by default; add E to also sample theta)
      tauprior = 3 0.2     # Inv-gamma(a, b) for root tau
    phiprior = 1 1  # Beta(a, b) for root tau & Dirichlet(a) for other tau's

      finetune =  1: 3 0.003 0.002 0.00002 0.005 0.9 0.001 0.001 # finetune for GBtj, GBspr, theta, tau, mix

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 16000
      sampfreq = 2
       nsample = 500000
