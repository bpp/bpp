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


	double mu_alpha = 20, mu_beta = 1000000; 
	double tau_alpha = 3, tau_beta = .4; 
	double tau_root, tau_AB, tau_ABC, mu, tmp; 
	double AB_max = 2000, ABC_max = 2000, ABCD_max = 2200; 

	double ** samples; 
	long int reps = 10000000;
	long int rejected = 0;
	int i;
	double mean[4] = {0,0,0,0};

	samples = malloc(sizeof(double *) * 5);
	/* Mu, theta_a, theta_b, coal_a, coal_b*/
	for (i = 0; i < 4; i++) {
		samples[i] = malloc(sizeof(double) * reps);
	}

	i = 0;

	while (i < reps) {
		mu = gsl_ran_gamma(r, mu_alpha, 1/mu_beta);

		tau_root = 1/gsl_ran_gamma(r, tau_alpha, 1/tau_beta);
		tau_AB = gsl_ran_flat(r,0, tau_root);
		tau_ABC = gsl_ran_flat(r,0, tau_root);
		if (tau_AB > tau_ABC) {
			tmp = tau_AB;
			tau_AB = tau_ABC;
			tau_ABC = tmp;
		}
		

		if ((AB_max * mu > tau_AB) || (ABC_max * mu > tau_ABC) || 
				(ABCD_max * mu > tau_root)) {
			rejected++;
			continue;
		} 
			samples[0][i] = mu;
			samples[1][i] = tau_root;
			samples[2][i] = tau_ABC;
			samples[3][i] = tau_AB;

		i++;
	}

	for (i = 0; i < reps; i++) {
		for (int j = 0; j < 4; j++) 
			mean[j] += samples[j][i]/reps;

	}
	printf("Means: %.10f %.10f %.10f %.10f\n", mean[0], mean[1], mean[2], mean[3]);
	//printf("Rejected %ld of %ld \n", rejected, reps);

}
