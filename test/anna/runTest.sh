#!/bin/bash

# In order to get the appropriate output to screen to run these tests, the "debug_opt_sim" flag needs to be turned on in the source code
# The output of the R script is the difference between the expected and realized means for coalescent times

#../../src/./bpp --simulate simple.ctl &> testSimple.txt 
../../src/./bpp --simulate simple.ctl | grep "Coalesce ("  | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testSimple.txt

../../src/./bpp --simulate simple2.ctl | grep "Coalesce ("  | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testSimple2.txt

../../src/./bpp --simulate simple2phase.ctl | grep "Coalesce ("  | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testSimple2phase.txt

../../src/./bpp --simulate simple3.ctl | grep "Coalesce ("  | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testSimple3.txt

../../src/./bpp --simulate simple4.ctl | grep "Coalesce ("  | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testSimple4.txt

Rscript testProgram.R
