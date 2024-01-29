seed =  1
*seed =  -1

seqfile = simulate.txt
treefile=simulate_trees.txt
Imapfile = simple.Imap.txt
datefile = dates.txt
*outfile = out.txt
*mcmcfile = mcmc.txt
seqDates = seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 3  A B C
                  5  5 5
		  ((A #0.004, B #0.004):.1 #0.0035, C #0.005):.2 #0.006;

*migration = 4
*A B 0.03
*B A 0.02
*A C 0.05
*C B 0.07

# phased data for population
phase =   0 0 0
*phase =   1 1

# use sequence likelihood
*usedata = 1

*loci&length = 1000 1000
loci&length = 10 1000
clock = 1
locusrate =0
model = 0


