bpp=~/bpp/src/bpp
cd ex_1
$bpp --simulate simulate.ctl  | grep "Coalesce (" | awk -F : '{print $3}' | sed 's/)//g' > coal_time

cd ../ex_2
$bpp --simulate simulate.ctl  | grep "Coalesce (" | awk -F : '{print $3}' | sed 's/)//g' > coal_time

cd ../ex_3
$bpp --simulate simulate.ctl > all
grep "Coalesce (" all | awk -F : '{print $3}' | sed 's/)//g' > coal_time
#grep "age:"  all  > debug
grep "age:"  all  | grep X | awk '{print $4}'> debug_X
grep "age:"  all  | grep Y | awk '{print $4}'> debug_Y

cd ../ex_4
$bpp --simulate simulate.ctl > all
grep "Coalesce (" all | awk -F : '{print $3}' | sed 's/)//g' > coal_time
#grep "age:"  all  > debug
grep "age:"  all  | grep X | awk '{print $4}'> debug_X
grep "age:"  all  | grep Y | awk '{print $4}'> debug_Y

cd ../ex_5
$bpp --simulate simulate.ctl  > all
grep "Coalesce (" all | awk -F : '{print $3}' | sed 's/)//g' > coal_time
grep "age:"  all  | grep X | awk '{print $4}'> debug_X
grep "age:"  all  | grep Y | awk '{print $4}'> debug_Y

cd ../ex_6
$bpp --simulate simulate.ctl |  grep 'A^a1' > trees
~/bpp/test/anna/parseGtree/compare_tree trees A^a1 A^a2 B^b1 B^b2 B^b2 A^a2 > coal_time
grep A^a1_A^a2 coal_time > A_A
grep B^b1_B^b2 coal_time > B_B 
grep A^a2_B^b2 coal_time > A_B 

#
# rm ex_*/coal_time ex_*/A_A ex_*/B_B ex_*/A_B ex_*/debug_*  ex_*/trees  ex_*/all
