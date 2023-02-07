seed =  2
*seed =  -1

seqfile = simulate_IM.txt
treefile=simulate_trees.txt
Imapfile = simple.Imap.txt
datefile = datesIM2.txt
seqDates = seqDates.txt

# fixed number of species/populations 
*speciesdelimitation = 0

# fixed species tree

species&tree = 4  A B C D
                  5 5 5 5
		  ((A #0.04, B #0.04):.05 #.03, (C #0.05, D #0.02 ):.07 #0.06):.1 #.03;

migration = 2
A B 0.2
B A 0.2

# phased data for population
phase =   0 0 0 0 
*phase =   1 1

# use sequence likelihood
*usedata = 1

loci&length = 900000 100
*loci&length = 10 1000
clock = 1
locusrate =0
model = 0

