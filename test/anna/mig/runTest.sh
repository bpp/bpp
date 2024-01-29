#!/bin/bash

# In order to get the appropriate output to screen to run these tests, the "debug_opt_sim" flag needs to be turned on in the source code
# The output of the R script is the difference between the expected and realized means for coalescent times

../../../src/bpp --simulate IM1/simpleIM2.ctl | grep "Coalesce (" | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testIM1.txt

../../../src/bpp --simulate IM2/simpleIM2.ctl | grep "Coalesce (" | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testIM2.txt

../../../src/bpp --simulate IM3/simpleIM2.ctl | grep "Coalesce (" | awk '{print $4 $5 $9}' | sed 's/)/,/g' |sed 's/.$//g' > testIM3.txt

cd IM4
~/bpp/src/bpp --simulate simpleIM2.ctl | grep 'A^a1' > test4
~/bpp/test/anna/parseGtree/./compare_tree test4 A^a1 B^b2 A^a2 B^b2 > coal &
wait 

grep a1 coal > coal_a1_b2
grep a2 coal > coal_a2_b2
