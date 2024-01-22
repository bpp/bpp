seed =  1
*seed =  -1

seqfile = simulate2.txt
treefile=simulate_trees2.txt
Imapfile = simple2.Imap.txt
datefile = dates.txt
*outfile = out.txt
*mcmcfile = mcmc.txt
seqDates = seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 4  A B C D
                  5  5 5 5
		  ((A #0.007, B #0.004):.1 #0.008, (C #0.005, D #.003):.1 #0.006)#.007:.2;
#((A #0.007, B #0.004):.008 #0.008, (C #0.005, D #.003):.005 #0.006)#.007:.016;

*migration = 4
*A B 0.03
*B A 0.02
*A C 0.05
*C B 0.07

# phased data for population
phase =   0 0 0 0
*phase =   1 1

# use sequence likelihood
*usedata = 1

*loci&length = 1000 1000
loci&length = 10 1000
clock = 1
locusrate =0
model = 0


