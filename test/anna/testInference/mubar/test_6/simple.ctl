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

species&tree = 4  A B C D
                  5  4 6 5
		  ((A #0.004, B #0.004):.1 #0.0035, (C #0.005, D #0.006):.07 #.0024)#.008:.2;


phase =   0 0 0 0


*loci&length = 1000 1000
loci&length = 10 1000
clock = 1
locusrate =0
model = 0


