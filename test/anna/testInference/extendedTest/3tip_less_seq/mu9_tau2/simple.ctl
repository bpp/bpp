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
                   10 10 10
		   ((A #0.002, B #0.003):.0035 #0.0035, C #0.005):.01 #.004 ;


phase =   0 0 0 


loci&length = 1000 1000
*loci&length = 10 1000
clock = 1
locusrate =0
model = 0


