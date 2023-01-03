#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <assert.h>

int main (int argc, char **argv) {
	/* Tree shape ((A,B)AB,(C,D)CD)ABCD */

	  /* Set up random number generator */
	gsl_rng * r;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
 
	gsl_rng_set(r, time(0));


	double alpha = 3, beta = 10; 
	double ABLower = .05 , CLower = .06, DLower = .05, ABCLower, rootLower = 0, time1, time2, timeRoot, tmp; 
	double ** ages; 
	long int samples = 100000000;
	long int rejected = 0;
	int i = 0;
	double mean[3] = {0,0,0};

	ages = malloc(sizeof(double *) * 3);
	for (i = 0; i < 3; i++) {
		ages[i] = malloc(sizeof(double) * samples);
	}

	i = 0;

	ABCLower = (ABLower > CLower ? ABLower : CLower);
	rootLower = (ABCLower > DLower ? ABCLower : DLower);

	while (i < samples) {
		timeRoot = gsl_ran_gamma(r, alpha, 1/beta);
		assert(timeRoot != 0);
		//printf("timeRoot %.15f %f\n", timeRoot, rootLower);
		if (timeRoot < rootLower) {
			rejected++;
			continue;
		}

                time1 = gsl_ran_flat(r, 0, timeRoot);
                time2 = gsl_ran_flat(r, 0, timeRoot);

                if (time1 > time2) {
                        tmp = time2;
                        time2 = time1;
                        time1 = tmp;
                }


		if (time1 < ABLower) {
			rejected++;
			continue;
		}
		if (time2 < ABCLower) {
			rejected++;
			continue;
		}
			
		ages[0][i] = timeRoot;
		ages[1][i] = time1; 
		ages[2][i] = time2; 

		//printf("%d %f %f %f\n",i, ages[0][i], ages[1][i], ages[2][i]); 
		i++;
	
	}

	for (i = 0; i < samples; i++) {
		for (int j = 0; j < 3; j++) 
			mean[j] += ages[j][i]/samples;

	}
	printf("Means: %f %f %f\n", mean[0], mean[1], mean[2]);
	printf("Rejected %ld of %ld \n", rejected, samples);

}
