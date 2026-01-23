seed = 12345

seqfile = sim_seqs.txt
treefile = sim_trees.txt
Imapfile = sim.Imap.txt

# Two populations A and B
species&tree = 2  A B
                  10 10
                  (A #0.004, B #0.004): 0.002 #0.004;

# Bidirectional migration with M = 0.5 (moderate migration)
# M = 4*N*m, so this represents about 0.5 migrants per generation
migration = 2
A B 0.5
B A 0.5

# Phased data (diploid, unphased)
phase = 0 0

# 50 loci, 500 bp each
loci&length = 50 500

# Strict molecular clock
clock = 1
locusrate = 0

# JC69 model
model = 0
