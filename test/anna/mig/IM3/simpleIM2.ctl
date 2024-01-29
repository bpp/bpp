seed =  2
*seed =  -1

seqfile = IM3/simulate_IM.txt
treefile=IM3/simulate_trees.txt
Imapfile = IM3/simple.Imap.txt
datefile = IM3/datesIM2.txt
seqDates = IM3/seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 2  A B
                  1  1
		  (A #0.04, B #0.04):.1 #.03;

migration = 2
A B 0.2
B A 0.2

# phased data for population
phase =   0 0
*phase =   1 1

# use sequence likelihood
*usedata = 1

loci&length = 8000000 100
*loci&length = 10 1000
clock = 1
locusrate =0
model = 0

