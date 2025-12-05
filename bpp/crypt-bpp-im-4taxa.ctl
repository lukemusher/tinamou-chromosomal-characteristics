          seed =  -1

       seqfile = 10kloci.4taxa.txt
      Imapfile = imap.4taxa.txt
       jobname = MSC-M_A00.rand10k.4taxa.output

  speciesdelimitation = 0 * fixed species tree
  speciestree = 0    * fixed species

  species&tree = 4  Cratr   Cstr Ceery  Ccgol
                    1  1  1  1
                 (Cratr,((Cstr,Ceery),Ccgol));

       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 9998 * number of data sets in seqfile
	model = gtr
	alphaprior = 1 1 4

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = gamma 2 100  # gamma(a, b) for theta
      tauprior = gamma 2 10   # gamma(a, b) for root tau
 
       wprior = 2 1
      migration = 2 
	              Cratr Cstr
	              Cstr Cratr

      finetune =  1

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 11000
      sampfreq = 2
       nsample = 100000
       threads = 10 42 1

