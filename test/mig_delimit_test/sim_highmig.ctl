seed = 54321

seqfile = sim_highmig_seqs.txt
treefile = sim_highmig_trees.txt
Imapfile = sim.Imap.txt

# Two populations A and B
species&tree = 2  A B
                  10 10
                  (A #0.01, B #0.01): 0.005 #0.01;

# Higher bidirectional migration with M = 2.0 (strong migration)
# This should give ~2-4 migration events per locus on average
migration = 2
A B 2.0
B A 2.0

# Phased data (diploid, unphased)
phase = 0 0

# 50 loci, 500 bp each
loci&length = 50 500

# Strict molecular clock
clock = 1
locusrate = 0

# JC69 model
model = 0
