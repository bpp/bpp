seed = 1
treefile = mytree.tre
Imapfile = myimap.txt
seqfile = mydata.txt
datefile = dates.txt
seqDates = seqDates.txt
species&tree = 2 A B
		 2 0
		((A #0.01, (B #0.01) Y[&phi=0.3,&tau-parent=no] #0.01 :.035)X #0.01:.035, (X[&phi=0.1,&tau-parent=no] # 0.01)Y # 0.035 )R #0.01 :.045;
		 #(A # 0.01, B #0.01)#0.01 : .1;

#(A # 0.01, B #0.01)#0.01 : .1;
phase = 0 0 
loci&length= 1000000 500
model = 0
clock = 1
locusrate = 0 

# Failed attempts at Newick strings
#((A #0.01 , (B #0.01)[&phi=0.3])X #0.01 :.1, (X#0.01 [&phi=0.1]Y#0.01 )R #0.01  : .15;
#((A, Y)X, (B, X)Y) R;	
#((A #0.01, (B #0.01)Y[&phi=0.2,&tau-parent=no]:0.01 #0.01)X[&phi=0.1,&tau-parent=no]:0.01 #0.01, (X)Y:0.02 #0.01)R #0.01;
#		(A #0.01 , B #0.01)R #0.01: .1;
#((A, Y)X, (B, X)Y) R;	
#((A #0.01 , (B #0.01)[&phi=0.3])X #0.05 :.05, (X#0.01 [&phi=0.1]Y#0.01 )R #0.01  : .1;
