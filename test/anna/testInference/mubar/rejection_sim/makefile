CC = gcc
CFLAGS=`pkg-config --cflags glib-2.0` -g -Wall -O3
LDADD=`pkg-config --libs glib-2.0`

testPrior: testPrior.c
	$(CC) $(CFLAGS) testPrior.c -lgsl -lgslcblas -lm -g -o testPrior

testPrior_4_tip_asymmetric: testPrior_4_tip_asymmetric.c
	$(CC) $(CFLAGS) testPrior_4_tip_asymmetric.c -lgsl -lgslcblas -lm -g -o testPrior_4_tip_asymmetric
testPrior_4_tip_symmetric: testPrior_4_tip_symmetric.c
	$(CC) $(CFLAGS) testPrior_4_tip_symmetric.c -lgsl -lgslcblas -lm -g -o testPrior_4_tip_symmetric
clean:
	$(RM) testPrior
	$(RM) testPrior_4_tip_asymmetric
	$(RM) testPrior_4_tip_symmetric
