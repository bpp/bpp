CC = gcc
CFLAGS=`pkg-config --cflags glib-2.0` -g -Wall -O3
LDADD=`pkg-config --libs glib-2.0`

testPrior: testPrior.c
	$(CC) $(CFLAGS) testPrior.c -lgsl -lgslcblas -lm -g -o testPrior

testPriorAsymmetric: testPriorAsymmetric.c
	$(CC) $(CFLAGS) testPriorAsymmetric.c -lgsl -lgslcblas -lm -g -o testPriorAsymmetric

clean:
	$(RM) testPrior
