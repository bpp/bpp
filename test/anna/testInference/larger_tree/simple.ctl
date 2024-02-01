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

species&tree = 10  A B C D E F G H I J
                   2 3 2 3 2 3 2 2 2 2   
		   ((((A #0.004, B #0.004):.02 #0.0035, (C #0.005, D #0.006):.04 #.0024):.05 #.0085, E #.002 )#.04:.08, (F #.0033, (G #.0037, (H #.005, (I #.006, J #.007)#.003:.01)#.0035:.03)#.005:.05)#.0044:.09)#.003:.1;


phase =   0 0 0 0 0 0 0 0 0 0


loci&length = 1000 1000
*loci&length = 10 1000
clock = 1
locusrate =0
model = 0


