seed = 98765

seqfile = sim_veryhighmig_seqs.txt
treefile = sim_veryhighmig_trees.txt
Imapfile = sim.Imap.txt

# Two populations A and B with larger theta and tau
species&tree = 2  A B
                  10 10
                  (A #0.02, B #0.02): 0.01 #0.02;

# Very high bidirectional migration M = 10
migration = 2
A B 10.0
B A 10.0

# Phased data (diploid, unphased)
phase = 0 0

# 50 loci, 500 bp each
loci&length = 50 500

# Strict molecular clock
clock = 1
locusrate = 0

# JC69 model
model = 0
