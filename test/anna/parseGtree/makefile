CC = gcc
CFLAGS=`pkg-config --cflags glib-2.0` -g -Wall -O3
LDADD=`pkg-config --libs glib-2.0`

compare_tree: compare_tree.o dna_sim_functions.o string_utilities.o
	$(CC) compare_tree.o dna_sim_functions.o string_utilities.o -o compare_tree $(LDADD) -lgsl -lgslcblas -lm -g

dna_sim_functions.o: dna_sim_functions.c dna_sim.h
	$(CC) $(CFLAGS) -c dna_sim_functions.c

string_utilities.o: string_utilities.c string_utilities.h
	$(CC) $(CFLAGS) -c string_utilities.c

sasprintf.o: sasprintf.c
	$(CC) $(CFLAGS) -c sasprintf.c -g -Wall -std=gnu11 -O3


compare_tree.o: compare_tree.c dna_sim.h string_utilities.h
	$(CC) -c compare_tree.c

clean:
	$(RM) sasprintf.o
	$(RM) sasprintf
	$(RM) string_utilities.o
	$(RM) compare_tree
	$(RM) dna_sim_functions.o
	$(RM) compare_tree.o
